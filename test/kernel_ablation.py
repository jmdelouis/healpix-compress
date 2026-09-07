"""Ablation of the transfer kernel that populates the carrier grid (V3).

Section~\\ref{sec:resampler} of the paper reports, qualitatively, that an
inverse-distance-weighted resampler left excess high-frequency power in the
carrier tile and that replacing it with Clough--Tocher removed most of it. That
observation decided the pipeline but was never measured: no spectral distance
was computed, and no rate--distortion consequence was quantified. Everything
the paper says about the kernel therefore rests on a visual comparison, which
is why the review asked for either the ablation or the removal of the claims.

This runs the ablation. For each kernel in

    nearest, linear, bicubic, clough_tocher, psf

the carrier tile is rebuilt from the same native window, and three quantities
are measured against the same observation:

    floor        the RMSE of the uncompressed tile read on the native pixels,
                 over the two footprints' intersection. This is the paper's
                 representation floor, now as a function of the kernel. It is
                 the transfer error and the reading error together, but the
                 read is identical across kernels, so differences between rows
                 are the transfer alone -- which is the decomposition the
                 operator characterisation of `operator_tests.py` makes
                 possible;

    log-PSD      the radially averaged power spectrum of the projected tile
                 against that of the observed crop, as an RMS distance in
                 log10 over the band, plus the ratio in the top octave alone.
                 This is the quantitative spectral metric the paper says it
                 never computed, and it is what tells excess high-frequency
                 content from a deficit;

    rate         with --rd, the compressed size at each target RMSE, so the
                 kernel's effect on the headline ratios is bounded rather than
                 assumed. This is the expensive half and is off by default.

Only the default kernel writes to the unsuffixed tile path, so nothing the
paper's results are built on is overwritten.

Usage
-----
    python kernel_ablation.py --limit 3                # floor + spectrum
    python kernel_ablation.py --limit 3 --rd           # + rate-distortion
    python kernel_ablation.py --kernels nearest psf --limit 2

Writes results/kernel_ablation.csv, results/kernel_ablation.md,
figures/fig_kernel_psd.pdf and ../tex/tables/table_kernel_ablation.tex.
"""

from __future__ import annotations

import argparse
import math
import shutil
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt           # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402

import paper_config as cfg                # noqa: E402
import paper_stats as st                  # noqa: E402
import build_tiles                        # noqa: E402
import common_reference                   # noqa: E402
import region_sites                       # noqa: E402
import support_controls as sup            # noqa: E402

REFERENCE_KERNEL = build_tiles.DEFAULT_KERNEL


def ensure_tile(scene: str, kernel: str, store: Path, build: bool, force=False):
    """Return the tile path, building it if it is not there yet."""
    path = build_tiles.tile_path(scene, kernel)
    if path.exists() and not force:
        return path, 0.0
    if not build:
        raise FileNotFoundError(f"{path} is missing and --no-build was given")
    if not store.exists():
        raise FileNotFoundError(f"missing window store {store}")
    t0 = time.perf_counter()
    build_tiles.build_tile(scene, store, kernel=kernel, verbose=False)
    return path, time.perf_counter() - t0


def spectral_distance(projected, observed, mask):
    """Log-PSD distance between a projected tile and the observation.

    The comparison is made on the largest rectangle inscribed in the two
    footprints' intersection: a spectrum taken over a masked octagon would
    measure the mask's edges rather than the field, and the whole point of this
    metric is to be about content.
    """
    rect = sup.largest_true_rectangle(mask)
    if rect is None:
        return dict(logpsd_rms=np.nan, logpsd_top_octave=np.nan,
                    psd_rect_px=0)
    r0, r1, c0, c1 = rect
    n = min(r1 - r0, c1 - c0)
    p = np.asarray(projected)[r0:r0 + n, c0:c0 + n]
    o = np.asarray(observed)[r0:r0 + n, c0:c0 + n]
    if not (np.all(np.isfinite(p)) and np.all(np.isfinite(o))):
        p = np.nan_to_num(p, nan=float(np.nanmean(p)))
        o = np.nan_to_num(o, nan=float(np.nanmean(o)))
    k, psd_p = sup.radial_psd(p, cfg.NATIVE_PIXEL_M)
    _k, psd_o = sup.radial_psd(o, cfg.NATIVE_PIXEL_M)
    k_nyq = math.pi / cfg.NATIVE_PIXEL_M
    band = (k > 0) & (k <= k_nyq) & (psd_p > 0) & (psd_o > 0)
    if band.sum() < 6:
        return dict(logpsd_rms=np.nan, logpsd_top_octave=np.nan, psd_rect_px=n)
    d = np.log10(psd_p[band] / psd_o[band])
    top = (k >= 0.5 * k_nyq) & band
    return dict(
        logpsd_rms=float(np.sqrt(np.mean(d ** 2))),
        # >1: the transfer leaves more high-frequency power than the
        # observation has, which is the failure mode the paper describes
        logpsd_top_octave=(float(psd_p[top].sum() / psd_o[top].sum())
                           if top.any() else np.nan),
        psd_rect_px=int(n))


def run_region(scene: str, kernel: str, store: Path, codecs, targets,
               build=True, force=False, do_rd=False, verbose=True):
    path, build_s = ensure_tile(scene, kernel, store, build, force)
    cell_ids, cell_data, utm_crop, meta, npz = build_tiles.load_tile(scene, kernel)
    cr = common_reference.CommonReference(cell_ids, np.asarray(npz["utm_lon"]),
                                          np.asarray(npz["utm_lat"]),
                                          mode="nearest")
    proj = cr.predict(cell_data)
    floor = cr.rmse_on_mask(proj, utm_crop)
    obs = np.asarray(utm_crop, dtype=np.float64)
    scene_std = float(np.std(obs[cr.mask]))

    row = dict(scene=scene, kernel=kernel, floor=floor,
               floor_frac_of_std=floor / scene_std if scene_std else np.nan,
               scene_std=scene_std, build_seconds=build_s,
               constant_check=meta.get("resampler_constant_check", np.nan),
               hpx_std=float(np.std(cell_data)))
    row.update(spectral_distance(proj, obs, cr.mask))

    if do_rd:
        tile = sup.expanded(cell_data)
        sweep = sup.sweep_array(tile, scene, codecs, label=f"hpx_{kernel}")
        for t in targets:
            row[f"size_at_{t:g}"] = st.best_of_codecs_size(sweep, scene, codecs, t)
            row[f"bracketed_at_{t:g}"] = st.any_bracketed(sweep, scene, codecs,
                                                          t, "rmse")
        sweep.to_csv(cfg.RESULTS_DIR / "kernel_ablation_sweeps.csv", mode="a",
                     header=not (cfg.RESULTS_DIR /
                                 "kernel_ablation_sweeps.csv").exists(),
                     index=False)
    if verbose:
        extra = ""
        if do_rd and np.isfinite(row.get(f"size_at_{targets[0]:g}", np.nan)):
            extra = (f" | size@{targets[0]:g} "
                     f"{row[f'size_at_{targets[0]:g}']/cfg.SHARED_RAW:.3f}")
        print(f"    {kernel:14s} floor {floor:.3e} "
              f"({100*row['floor_frac_of_std']:.1f}% of std) | "
              f"log-PSD {row['logpsd_rms']:.3f} | "
              f"top octave {row['logpsd_top_octave']:.2f}x{extra}")
    return row


# --- reporting ---------------------------------------------------------------
def paired(df, col, reference=REFERENCE_KERNEL):
    """Per-region ratio of each kernel's column to the reference kernel's."""
    ref = df[df.kernel == reference].set_index("scene")[col]
    out = []
    for kernel, g in df.groupby("kernel"):
        v = g.set_index("scene")[col]
        common = v.index.intersection(ref.index)
        if not len(common):
            continue
        r = (v.loc[common] / ref.loc[common]).to_numpy(dtype=float)
        r = r[np.isfinite(r)]
        out.append(dict(kernel=kernel, n=len(r),
                        median=float(np.median(r)) if r.size else np.nan,
                        lo=float(np.min(r)) if r.size else np.nan,
                        hi=float(np.max(r)) if r.size else np.nan))
    return pd.DataFrame(out)


def summarize(df, targets, do_rd):
    out = ["# Transfer-kernel ablation (V3)", "",
           f"- regions: {df.scene.nunique()}, kernels: "
           f"{', '.join(sorted(df.kernel.unique()))}",
           f"- reference kernel: `{REFERENCE_KERNEL}` (the one the paper uses)",
           ""]
    out += ["## Absolute measurements (median over regions)", "",
            "| kernel | floor | floor / scene std | log-PSD RMS | top-octave "
            "power ratio | build (s) |", "|---|---|---|---|---|---|"]
    for kernel, g in df.groupby("kernel"):
        out.append(f"| `{kernel}` | {g.floor.median():.3e} | "
                   f"{100*g.floor_frac_of_std.median():.1f}% | "
                   f"{g.logpsd_rms.median():.3f} | "
                   f"{g.logpsd_top_octave.median():.2f} | "
                   f"{g.build_seconds.replace(0, np.nan).median():.0f} |")
    out.append("")
    out.append("A top-octave ratio above 1 means the transfer leaves more "
               "high-frequency power than the observation carries, which is "
               "the failure the paper describes qualitatively; below 1 means "
               "it smooths.")
    out.append("")

    out += ["## Paired against the reference kernel", "",
            "| kernel | n | floor ratio | range |", "|---|---|---|---|"]
    for _, r in paired(df, "floor").iterrows():
        out.append(f"| `{r.kernel}` | {int(r.n)} | {r['median']:.3f} | "
                   f"[{r.lo:.3f}, {r.hi:.3f}] |")
    out.append("")

    if do_rd:
        out += ["## Compressed size at matched RMSE, paired against the "
                "reference kernel", "",
                "| kernel | " + " | ".join(f"{t:g}" for t in targets) + " |",
                "|---" * (len(targets) + 1) + "|"]
        for kernel in sorted(df.kernel.unique()):
            cells = []
            for t in targets:
                col = f"size_at_{t:g}"
                if col not in df:
                    cells.append("--")
                    continue
                p = paired(df, col)
                v = p[p.kernel == kernel]
                cells.append(f"{v['median'].iloc[0]:.3f}" if len(v) else "--")
            out.append(f"| `{kernel}` | " + " | ".join(cells) + " |")
        out.append("")
        out.append("Values are per-region ratios of compressed size to the "
                   "reference kernel's, so 1.00 means the kernel changes "
                   "nothing at that target.")
        out.append("")

    ref = df[df.kernel == REFERENCE_KERNEL]
    others = df[df.kernel != REFERENCE_KERNEL]
    if len(ref) and len(others):
        best = others.groupby("kernel").floor.median().idxmin()
        bestv = others.groupby("kernel").floor.median().min()
        refv = ref.floor.median()
        out.append(f"- lowest floor: `{best}` at {bestv:.3e}, against "
                   f"{refv:.3e} for `{REFERENCE_KERNEL}` "
                   f"({100*(bestv/refv - 1):+.1f}%)")
    return "\n".join(out)


def figure(df, publish=True, name="fig_kernel_psd"):
    if not len(df):
        return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.6, 2.9))
    kernels = sorted(df.kernel.unique())
    cmap = plt.get_cmap("viridis")
    colors = {k: cmap(i / max(len(kernels) - 1, 1)) for i, k in enumerate(kernels)}
    for i, k in enumerate(kernels):
        g = df[df.kernel == k]
        ax1.scatter(np.full(len(g), i), 100 * g.floor_frac_of_std,
                    color=colors[k], s=18, alpha=0.75, linewidths=0)
        ax1.plot([i], [100 * g.floor_frac_of_std.median()], marker="_", ms=18,
                 mew=2.2, color=colors[k])
        ax2.scatter(np.full(len(g), i), g.logpsd_top_octave, color=colors[k],
                    s=18, alpha=0.75, linewidths=0)
        ax2.plot([i], [g.logpsd_top_octave.median()], marker="_", ms=18,
                 mew=2.2, color=colors[k])
    for ax, ylab in ((ax1, "floor (% of scene std)"),
                     (ax2, "top-octave power vs observation")):
        ax.set_xticks(range(len(kernels)))
        ax.set_xticklabels([k.replace("clough_tocher", "C--T") for k in kernels],
                           rotation=30, ha="right", fontsize=8)
        ax.set_ylabel(ylab, fontsize=9)
        ax.tick_params(labelsize=8)
    ax2.axhline(1.0, color="0.6", ls="--", lw=1.0)
    ax2.set_yscale("log")
    fig.suptitle("Transfer kernel: representation floor and spectral fidelity",
                 fontsize=10)
    fig.tight_layout()
    cfg.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(cfg.FIGURES_DIR / f"{name}.{ext}", bbox_inches="tight", dpi=160)
    if publish:
        cfg.TEX_FIGURES.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(cfg.FIGURES_DIR / f"{name}.pdf",
                        cfg.TEX_FIGURES / f"{name}.pdf")
    plt.close(fig)
    print(f"  {name}")


def latex_table(df, targets, do_rd):
    rows = []
    fp = paired(df, "floor").set_index("kernel")
    for kernel, g in df.groupby("kernel"):
        ratio = (f"{fp.loc[kernel,'median']:.2f}" if kernel in fp.index else "--")
        cells = ""
        if do_rd:
            for t in targets[:3]:
                col = f"size_at_{t:g}"
                p = paired(df, col).set_index("kernel") if col in df else None
                cells += (f" & {p.loc[kernel,'median']:.2f}"
                          if p is not None and kernel in p.index else " & --")
        rows.append(
            f"\\texttt{{{kernel.replace('_', chr(92)+'_')}}} & "
            f"{g.floor.median():.2e} & {100*g.floor_frac_of_std.median():.1f} & "
            f"{ratio} & {g.logpsd_rms.median():.2f} & "
            f"{g.logpsd_top_octave.median():.2f}" + cells + " \\\\")
    extra_h = "".join(f" & ${t:g}$" for t in targets[:3]) if do_rd else ""
    spec = "lrrrrr" + ("r" * (3 if do_rd else 0))
    return (
        "\\begin{table}[t]\n\\centering\n"
        "\\caption{Transfer-kernel ablation. Each row rebuilds the carrier tile "
        "from the same native window with a different kernel and scores it "
        "against the same observation with the same reading operator, so "
        "differences between rows are the transfer alone. \\emph{Floor} is the "
        "RMSE of the uncompressed tile read on the native pixels over the two "
        "footprints' intersection, absolute and as a fraction of the scene's "
        "standard deviation; \\emph{rel.} is the per-region ratio to "
        "\\texttt{clough\\_tocher}, the kernel used throughout this paper. "
        "\\emph{log-PSD} is the RMS distance in $\\log_{10}$ between the "
        "projected tile's radial power spectrum and the observation's over the "
        "band up to the native Nyquist, and \\emph{top oct.} the ratio of "
        "their power in the highest octave alone: above 1 the transfer leaves "
        "more high-frequency power than the observation carries, which is the "
        "failure this paper reports qualitatively for its first resampler."
        + (" The last columns are per-region compressed-size ratios to the "
           "same reference at three RMSE targets." if do_rd else "") + "}\n"
        "\\label{tab:kernel_ablation}\n\\footnotesize\n"
        f"\\begin{{tabular}}{{@{{}}{spec}@{{}}}}\n\\toprule\n"
        "Kernel & Floor & \\% std & rel. & log-PSD & top oct." + extra_h
        + " \\\\\n\\midrule\n" + "\n".join(rows)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kernels", nargs="*", default=list(build_tiles.KERNELS))
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=3)
    ap.add_argument("--codecs", nargs="*", default=cfg.COMPRESSOR_CODECS)
    ap.add_argument("--targets", type=float, nargs="*", default=cfg.TARGET_RMSE)
    ap.add_argument("--rd", action="store_true",
                    help="also run the codec sweep per kernel (the expensive "
                         "half; without it only the floor and the spectrum are "
                         "measured)")
    ap.add_argument("--no-build", action="store_true",
                    help="fail instead of building a missing tile")
    ap.add_argument("--force", action="store_true",
                    help="rebuild tiles even when they exist")
    ap.add_argument("--no-publish", action="store_true")
    ap.add_argument("--no-tex", action="store_true")
    args = ap.parse_args(argv)

    sites = region_sites.load_sites()
    scenes = args.regions or list(sites["scene"])
    scenes = [s for s in scenes
              if (cfg.REGION_DIR / f"{s}_data.zarr").exists()
              or build_tiles.tile_path(s).exists()]
    if not scenes:
        print(f"no region windows under {cfg.REGION_DIR} and no tiles under "
              f"{cfg.TILE_DIR}")
        return 1
    if args.limit:
        scenes = scenes[:args.limit]

    print(f"{len(scenes)} region(s), kernels {args.kernels}, "
          f"rate-distortion {'on' if args.rd else 'off'}")
    rows, failures = [], []
    for scene in scenes:
        print(f"  {scene}")
        store = cfg.REGION_DIR / f"{scene}_data.zarr"
        for kernel in args.kernels:
            try:
                rows.append(run_region(scene, kernel, store, args.codecs,
                                       args.targets, build=not args.no_build,
                                       force=args.force, do_rd=args.rd))
            except Exception as exc:                       # noqa: BLE001
                print(f"    {kernel:14s} FAILED {type(exc).__name__}: {exc}")
                failures.append((scene, kernel, f"{type(exc).__name__}: {exc}"))
    if not rows:
        print("no (region, kernel) pair produced a row")
        return 1

    df = pd.DataFrame(rows)
    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    csv = cfg.RESULTS_DIR / "kernel_ablation.csv"
    df.to_csv(csv, index=False)
    md = summarize(df, args.targets, args.rd)
    (cfg.RESULTS_DIR / "kernel_ablation.md").write_text(md, encoding="utf-8")
    print("\n" + md)
    print(f"wrote {csv} ({len(df)} rows)")

    figure(df, publish=not args.no_publish)
    if not args.no_tex:
        out = cfg.TEX_TABLES / "table_kernel_ablation.tex"
        cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
        out.write_text(latex_table(df, args.targets, args.rd), encoding="utf-8")
        print(f"wrote {out}")

    if failures:
        print(f"\n{len(failures)} (region, kernel) pair(s) failed:")
        for s, k, why in failures:
            print(f"  {s} {k}: {why}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
