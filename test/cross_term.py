"""Is the compression error orthogonal to the error the product already carries?

Section V-H expresses the storage verdict in units of the transfer error and
converts a compression target into a growth of the product's total error with

    D_total / D_0 = sqrt(1 + kappa^2),      kappa = D_c / D_0,

which is exact only if the two errors are orthogonal. They are not obviously
so: both are produced by operations applied to the same spatial structures, and
the general relation carries a cross term,

    D_total / D_0 = sqrt(1 + kappa^2 + 2 rho kappa),
    rho = <e_0, e_c> / (D_0 D_c).

A negative rho means the codec happens to cancel part of the transfer error and
the quadrature figure overstates the damage; a positive rho means it compounds
it and the quadrature figure understates it. At kappa = 0.5 the difference
between rho = -1 and rho = +1 is the difference between a product whose error
*halves* and one whose error grows by half, so the assumption is not a detail.

This module measures rho instead of assuming it. Both errors are taken on the
same support -- the level-20 cells -- and in the same norm, which is the other
thing the assumption needs and which a comparison made after projecting to UTM
would not give.

    e_0 = x - truth      the transfer error, x being the transferred cells
    e_c = x~ - x         the compression error, x~ the decoded cells
    e_0 + e_c = x~ - truth

The scene is the synthetic one of `resampling_error.py`, built through the same
`synthetic_case` so that the two measurements cannot drift apart; the codec is
run on the tile that case produces. What is measured is therefore rho for a
synthetic field with each region's own fitted spectrum and noise, which is an
estimate of rho for the real tiles rather than a measurement of it.

Usage
-----
    python cross_term.py --limit 4              # smoke test
    python cross_term.py                        # all 40 regions
    python cross_term.py --kernels clough_tocher psf
"""

from __future__ import annotations

import argparse
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt           # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402

import paper_config as cfg                # noqa: E402
import codec_lib                          # noqa: E402
import region_sites                       # noqa: E402
import build_tiles                        # noqa: E402
import resampling_error as R              # noqa: E402

# One image codec and one floating-point codec. The point is the geometry of
# the two error fields, not a codec ranking, and a full sweep on 40 synthetic
# tiles costs as much as the paper's main sweep for no extra information.
DEFAULT_CODECS = ["jpeg2000_lossy", "zfp"]
KAPPA_REPORT = (0.1, 0.2, 0.5, 1.0)


def rms(a):
    return float(np.sqrt(np.mean(np.asarray(a, dtype=np.float64) ** 2)))


def region_rows(scene, store, spectrum, kernels, codecs, lam, children,
                noise_scale, seed, reference_csv=None, verbose=True):
    from hcc.compression.zorder import expand_tile, flatten_tile

    case = R.synthetic_case(scene, store, spectrum, noise_scale=noise_scale,
                            seed=seed, children=children)
    rows = []
    for kernel in kernels:
        estimand = R.ESTIMAND.get(kernel, "observed")
        truth = (case.latent_cells if estimand == "latent"
                 else case.observed_cells)
        try:
            val, n_filled = R.transfer_cells(case, kernel, lam=lam)
        except Exception as exc:                          # noqa: BLE001
            print(f"    {kernel:14s} FAILED {type(exc).__name__}: {exc}")
            continue
        e0 = np.asarray(val, dtype=np.float64) - np.asarray(truth, np.float64)
        d0 = rms(e0)
        if not np.isfinite(d0) or d0 == 0:
            print(f"    {kernel:14s} skipped: zero transfer error")
            continue

        # The refactor that made this script possible moved the field
        # construction out of `measure`. If the published measurement is on
        # disk, check we reproduce its transfer error before trusting anything
        # built on top of it.
        drift = np.nan
        if reference_csv is not None:
            m = reference_csv[(reference_csv.scene == scene)
                              & (reference_csv.kernel == kernel)]
            if len(m):
                ref = float(m.e_resample.iloc[0])
                drift = abs(d0 - ref) / ref if ref else np.nan

        tile = expand_tile(np.asarray(val, dtype=np.float32)).astype(np.float32)
        for run in codec_lib.iter_runs(tile, codecs):
            if np.shape(run.reconstruction) != np.shape(tile):
                continue
            back = flatten_tile(np.asarray(run.reconstruction, dtype=np.float32))
            ec = np.asarray(back, dtype=np.float64) - np.asarray(val, np.float64)
            dc = rms(ec)
            if not np.isfinite(dc) or dc == 0:
                continue
            rho = float(np.mean(e0 * ec) / (d0 * dc))
            dt = rms(e0 + ec)
            rows.append(dict(
                scene=scene, cls=scene.split("__")[1], kernel=kernel,
                estimand=estimand, codec=run.codec,
                d0=d0, dc=dc, kappa=dc / d0, rho=rho,
                inflation_measured=dt / d0 - 1.0,
                inflation_quadrature=float(np.sqrt(1.0 + (dc / d0) ** 2)) - 1.0,
                total_bytes=int(run.total_bytes),
                cells_filled=int(n_filled), d0_drift_vs_published=drift))
        if verbose and rows:
            k = [r for r in rows if r["kernel"] == kernel]
            print(f"    {kernel:14s} D0 {d0:.4e}"
                  + (f" (published within {100*drift:.2g}%)"
                     if drift == drift else "")
                  + f" | {len(k)} operating point(s), "
                    f"rho median {np.median([r['rho'] for r in k]):+.3f}")
    return rows


def summarize(df):
    out = ["# Is the compression error orthogonal to the transfer error?", "",
           f"- regions: {df.scene.nunique()}, kernels: "
           f"{', '.join(sorted(df.kernel.unique()))}, codecs: "
           f"{', '.join(sorted(df.codec.unique()))}",
           f"- operating points: {len(df)}", ""]
    d = df.d0_drift_vs_published.dropna()
    if len(d):
        out += [f"- transfer error reproduces the published measurement to "
                f"{100*float(d.max()):.2g}% at worst over {len(d)} "
                f"(region, kernel) pairs" + (" **-- investigate**"
                                             if float(d.max()) > 1e-6 else ""),
                ""]
    out += ["## rho by kernel, over all operating points", "",
            "| kernel | n | median rho | 10th | 90th | share |rho| > 0.1 |",
            "|---|---|---|---|---|---|"]
    for k, g in df.groupby("kernel"):
        out.append(f"| `{k}` | {len(g)} | {g.rho.median():+.3f} | "
                   f"{g.rho.quantile(.1):+.3f} | {g.rho.quantile(.9):+.3f} | "
                   f"{100*float((g.rho.abs() > 0.1).mean()):.0f}% |")
    out += ["", "## what the quadrature assumption costs, near the reported "
            "kappa values", "",
            "| kernel | kappa | n | measured inflation | quadrature | "
            "difference |", "|---|---|---|---|---|---|"]
    for k, g in df.groupby("kernel"):
        for kap in KAPPA_REPORT:
            w = g[(g.kappa > kap / 1.5) & (g.kappa < kap * 1.5)]
            if len(w) < 5:
                continue
            m, q = w.inflation_measured.median(), w.inflation_quadrature.median()
            out.append(f"| `{k}` | ~{kap:g} | {len(w)} | {100*m:+.1f}% | "
                       f"{100*q:+.1f}% | {100*(m - q):+.1f} pt |")
    out.append("")
    rho_all = float(df.rho.median())
    small = float((df.rho.abs() <= 0.1).mean())
    if abs(rho_all) <= 0.1 and small >= 0.8:
        out.append(f"- **The cross term is small**: median $\\rho$ = "
                   f"{rho_all:+.3f}, and {100*small:.0f}% of operating points "
                   f"have $|\\rho| \\le 0.1$. Quadrature is then an accurate "
                   f"description rather than an assumption, and Section V-H "
                   f"can be stated as a measurement.")
    else:
        out.append(f"- **The cross term is not negligible**: median $\\rho$ = "
                   f"{rho_all:+.3f}, and only {100*small:.0f}% of operating "
                   f"points have $|\\rho| \\le 0.1$. Section V-H must either "
                   f"carry $\\rho$ explicitly or be declared a scenario under "
                   f"an orthogonality assumption.")
    return "\n".join(out)


def figure(df, publish=True, name="fig_cross_term"):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.1, 2.9))
    for k, g in df.groupby("kernel"):
        ax1.scatter(g.kappa, g.rho, s=6, alpha=0.35, label=k)
        gg = g.sort_values("kappa")
        ax2.scatter(gg.kappa, 100 * gg.inflation_measured, s=6, alpha=0.35,
                    label=k)
    kk = np.logspace(-2, 0.6, 200)
    ax2.plot(kk, 100 * (np.sqrt(1 + kk ** 2) - 1), "k-", lw=1.2,
             label=r"$\sqrt{1+\kappa^2}-1$")
    ax1.axhline(0.0, color="0.6", ls="--", lw=1.0)
    for ax, yl in ((ax1, r"$\rho$"), (ax2, "total-error inflation (%)")):
        ax.set_xscale("log")
        ax.set_xlabel(r"$\kappa = e_{\rm comp}/e_{\rm resample}$", fontsize=9)
        ax.set_ylabel(yl, fontsize=9)
        ax.tick_params(labelsize=8)
    ax2.set_yscale("symlog", linthresh=1.0)
    ax1.set_ylim(-1.05, 1.05)
    ax1.legend(fontsize=8, frameon=False)
    ax2.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    cfg.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(cfg.FIGURES_DIR / f"{name}.{ext}", bbox_inches="tight",
                    dpi=160)
    if publish:
        cfg.TEX_FIGURES.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(cfg.FIGURES_DIR / f"{name}.pdf",
                        cfg.TEX_FIGURES / f"{name}.pdf")
    plt.close(fig)
    print(f"  {name}")


def latex_table(df):
    body = []
    for k, g in sorted(df.groupby("kernel")):
        for kap in KAPPA_REPORT:
            w = g[(g.kappa > kap / 1.5) & (g.kappa < kap * 1.5)]
            if len(w) < 5:
                continue
            m, q = w.inflation_measured.median(), w.inflation_quadrature.median()
            body.append(
                f"\\texttt{{{k.replace('_', chr(92)+'_')}}} & ${kap:g}$ & "
                f"{len(w)} & {w.rho.median():+.3f} & "
                f"[{w.rho.quantile(.1):+.2f}, {w.rho.quantile(.9):+.2f}] & "
                f"{100*m:+.1f} & {100*q:+.1f} \\\\")
    return (
        "\\begin{table}[t]\n\\centering\n"
        "\\caption{The cross term between the transfer error and the "
        "compression error, measured rather than assumed. Both errors are "
        "taken on the level-20 cells and in the same norm; $\\rho$ is their "
        "normalised inner product, so the exact relation is "
        "$D_{\\mathrm{total}}/D_0=\\sqrt{1+\\kappa^2+2\\rho\\kappa}$ and the "
        "quadrature form used in Section~\\ref{sec:error_budget} is its "
        "$\\rho=0$ case. Operating points are grouped around each reported "
        "$\\kappa$ within a factor of $1.5$. The last two columns are the "
        "growth of the product's total error as measured, and as quadrature "
        "predicts it.}\n"
        "\\label{tab:cross_term}\n\\footnotesize\n"
        "\\begin{tabular}{@{}llrrrrr@{}}\n\\toprule\n"
        "Transfer & $\\kappa$ & $n$ & median $\\rho$ & $[10,90]\\%$ & "
        "Measured (\\%) & Quadrature (\\%) \\\\\n\\midrule\n"
        + "\n".join(body) + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kernels", nargs="*", default=["clough_tocher"])
    ap.add_argument("--codecs", nargs="*", default=DEFAULT_CODECS)
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--noise", type=float, default=1.0)
    ap.add_argument("--lam", type=float, default=1e-3)
    ap.add_argument("--children", type=int, default=16, choices=[1, 4, 16])
    ap.add_argument("--seed", type=int, default=23)
    ap.add_argument("--no-publish", action="store_true")
    ap.add_argument("--no-tex", action="store_true")
    args = ap.parse_args(argv)

    spec_path = cfg.RESULTS_DIR / "support_controls_spectrum.csv"
    spectrum = pd.read_csv(spec_path) if spec_path.exists() else None
    ref_path = cfg.RESULTS_DIR / "resampling_error.csv"
    reference = pd.read_csv(ref_path) if ref_path.exists() else None
    if reference is None:
        print("note: no resampling_error.csv, so the transfer error cannot be "
              "checked against the published measurement")

    scenes = args.regions or list(region_sites.load_sites()["scene"])
    scenes = [s for s in scenes
              if (cfg.REGION_DIR / f"{s}_data.zarr").exists()
              and build_tiles.tile_path(s).exists()]
    if args.limit:
        scenes = scenes[:args.limit]
    if not scenes:
        print("no region with both a window and a tile")
        return 1

    print(f"{len(scenes)} region(s), kernels {args.kernels}, "
          f"codecs {args.codecs}")
    rows = []
    for i, scene in enumerate(scenes, 1):
        print(f"[{i:>2}/{len(scenes)}] {scene}")
        try:
            rows += region_rows(scene, cfg.REGION_DIR / f"{scene}_data.zarr",
                                spectrum, args.kernels, args.codecs, args.lam,
                                args.children, args.noise, args.seed,
                                reference_csv=reference)
        except Exception as exc:                          # noqa: BLE001
            print(f"    FAILED {type(exc).__name__}: {exc}")
    if not rows:
        print("no operating point produced")
        return 1

    df = pd.DataFrame(rows)
    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    csv = cfg.RESULTS_DIR / "cross_term.csv"
    df.to_csv(csv, index=False)
    md = summarize(df)
    (cfg.RESULTS_DIR / "cross_term.md").write_text(md, encoding="utf-8")
    print("\n" + md)
    print(f"\nwrote {csv} ({len(df)} rows)")
    figure(df, publish=not args.no_publish)
    if not args.no_tex:
        cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
        out = cfg.TEX_TABLES / "table_cross_term.tex"
        out.write_text(latex_table(df), encoding="utf-8")
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
