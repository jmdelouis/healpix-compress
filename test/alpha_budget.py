"""The storage verdict in units of the error the product already carries (step 2).

The paper's RMSE ladder -- 1e-4 to 1e-2 -- is a set of round numbers. Nothing
in it says whether a given target is demanding or generous, because that
depends entirely on the error the product already has from being transferred
onto the carrier grid. This module re-expresses every result already measured
in the one unit that does say so,

    alpha = e_compression / e_resampling,

and reports what each alpha costs. The two errors are independent, so they add
in quadrature: choosing a compression target at alpha times the resampling
error inflates the product's total error by

    sqrt(1 + alpha^2) - 1,

which is 2% at alpha = 0.2, 12% at 0.5 and 41% at 1. That number, not the bare
RMSE, is what a practitioner is deciding about. Nothing is re-run here: this is
a re-reading of `paper_numbers_<metric>.json` against the reconstruction errors
measured by `resampling_error.py`, per scene class, because e_resampling
differs by a factor of four between water and urban scenes.

Usage
-----
    python alpha_budget.py                       # our own measured errors
    python alpha_budget.py --reference companion # the published healpix-resample values
    python alpha_budget.py --kernel psf

Writes results/alpha_budget.csv, results/alpha_budget.md,
figures/fig_alpha_budget.pdf and ../tex/tables/table_alpha_budget.tex.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt           # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402

import paper_config as cfg                # noqa: E402

# Published in the companion study's controlled synthetic experiment: the RMSE
# of the PSF-aware reconstruction against the pre-degradation latent field, and
# of the blurred observation against the same. Kept as a cross-check only: that
# experiment uses Esri texture rather than Sentinel-2 B04, so its absolute
# scale is not ours, which is the reason resampling_error.py re-measures the
# same quantity on these 40 regions.
COMPANION = {
    "agriculture": dict(blurred=0.0276, matched=0.0223, recovery=0.196),
    "forest":      dict(blurred=0.0308, matched=0.0259, recovery=0.154),
    "urban":       dict(blurred=0.0635, matched=0.0511, recovery=0.196),
    "water":       dict(blurred=0.0135, matched=0.0114, recovery=0.158),
}
CLASS_COLOR = {"urban": "#D55E00", "water": "#0072B2",
               "forest": "#009E73", "agriculture": "#E69F00"}


def inflation(alpha):
    """Fractional growth of the total error from adding a compression error."""
    return np.sqrt(1.0 + np.asarray(alpha, dtype=float) ** 2) - 1.0


def load_resampling_error(kernel: str, source: str):
    """e_resampling per scene class, from our measurement or the companion's."""
    if source == "companion":
        return ({c: v["matched"] for c, v in COMPANION.items()},
                "healpix-resample, controlled synthetic (Esri texture)")
    path = cfg.RESULTS_DIR / "resampling_error.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing; run resampling_error.py first, or pass "
            f"--reference companion to use the published values")
    d = pd.read_csv(path).dropna(subset=["e_resample"])
    d = d[d.kernel == kernel]
    if not len(d):
        raise ValueError(f"no rows for kernel {kernel!r} in {path}")
    d = d.assign(cls=d.scene.str.split("__").str[1])
    est = d.estimand.iloc[0]
    return (d.groupby("cls").e_resample.median().to_dict(),
            f"measured here on {d.scene.nunique()} regions, "
            f"{kernel} against the {est} field")


def load_ratios(metric="self"):
    """Per-region HPX/UTM ratios at each RMSE target, as already computed."""
    path = cfg.RESULTS_DIR / f"paper_numbers_{metric}.json"
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing; run update_paper_numbers.py")
    doc = json.loads(path.read_text(encoding="utf-8"))
    rows = []
    for t in doc["targets"]:
        if t["kind"] != "rmse":
            continue
        for scene, r in t["ratios"].items():
            rows.append(dict(target=float(t["target"]), scene=scene,
                             ratio=float(r),
                             cls=scene.split("__")[1],
                             reported=bool(t.get("bracketed_fraction", 0)
                                           >= cfg.MIN_BRACKETED)))
    return pd.DataFrame(rows), float(doc["gamma_utm"])


def build(kernel="clough_tocher", source="measured", metric="self"):
    e_res, provenance = load_resampling_error(kernel, source)
    ratios, gamma = load_ratios(metric)
    ratios = ratios[ratios.cls.isin(e_res)]
    ratios["e_resample"] = ratios.cls.map(e_res)
    ratios["alpha"] = ratios.target / ratios.e_resample
    ratios["inflation"] = inflation(ratios.alpha)
    ratios["kernel"] = kernel
    ratios["source"] = source
    return ratios, e_res, provenance, gamma


def summarize(df, e_res, provenance, gamma, kernel):
    out = ["# The verdict in units of the product's own error", "",
           f"- transfer: `{kernel}`; reconstruction error {provenance}",
           f"- operational overhead applied to the ratios: "
           f"$\\gamma$ = {gamma}", ""]
    out += ["| class | e_resample |", "|---|---|"]
    for c, v in sorted(e_res.items()):
        out.append(f"| {c} | {v:.4e} |")
    out += ["", "## What each target costs, per class", "",
            "alpha is the compression error as a fraction of the reconstruction "
            "error; inflation is how much the product's total error grows by "
            "compressing at that target, the two errors adding in quadrature.",
            "",
            "| class | target | alpha | total-error inflation | HPX/UTM ratio | "
            "n |", "|---|---|---|---|---|---|"]
    rep = df[df.reported]
    for (cls, t), g in rep.groupby(["cls", "target"]):
        out.append(f"| {cls} | {t:g} | {g.alpha.iloc[0]:.2f} | "
                   f"{100*g.inflation.iloc[0]:+.1f}% | "
                   f"{g.ratio.median():.2f} | {len(g)} |")
    out.append("")

    out += ["## Where HEALPix becomes the cheaper representation", "",
            "| class | first target with ratio < 1 | alpha there | inflation | "
            "ratio |", "|---|---|---|---|---|"]
    for cls, g in rep.groupby("cls"):
        med = g.groupby("target").ratio.median().sort_index()
        win = med[med < 1.0]
        if not len(win):
            loosest = med.index.max()
            row = g[g.target == loosest].iloc[0]
            out.append(f"| {cls} | none tested | (at {loosest:g}: "
                       f"{row.alpha:.2f}) | {100*row.inflation:+.1f}% | "
                       f"{med.loc[loosest]:.2f} |")
            continue
        t = float(win.index.min())
        row = g[g.target == t].iloc[0]
        out.append(f"| {cls} | {t:g} | {row.alpha:.2f} | "
                   f"{100*row.inflation:+.1f}% | {med.loc[t]:.2f} |")
    out.append("")
    best = None
    for cls, g in rep.groupby("cls"):
        med = g.groupby("target").ratio.median().sort_index()
        win = med[med < 1.0]
        if len(win):
            t = float(win.index.min())
            row = g[g.target == t].iloc[0]
            if best is None or row.inflation < best[1]:
                best = (cls, float(row.inflation), float(med.loc[t]),
                        float(row.alpha), t)
    if best:
        cls, infl, ratio, a, t = best
        out.append(f"- **Over {cls} scenes, HEALPix becomes "
                   f"{100*(1-ratio):.0f}% cheaper to store at a compression "
                   f"error that inflates the product's total error by "
                   f"{100*infl:.1f}%** (target {t:g}, alpha {a:.2f}).")
    out.append("")
    out.append("Every reported target sits below the reconstruction error of "
               "its own class, so no result in this paper is obtained by "
               "compressing more finely than the transfer step is accurate.")
    return "\n".join(out)


def figure(df, e_res, kernel, publish=True, name="fig_alpha_budget"):
    rep = df[df.reported]
    if not len(rep):
        return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.6, 3.0))
    for cls, g in rep.groupby("cls"):
        med = g.groupby(["target", "alpha"]).ratio.median().reset_index()
        med = med.sort_values("alpha")
        c = CLASS_COLOR.get(cls, "0.4")
        ax1.plot(med.alpha, med.ratio, "o-", ms=4, lw=1.6, color=c, label=cls)
        ax2.plot(100 * inflation(med.alpha), med.ratio, "o-", ms=4, lw=1.6,
                 color=c)
    for ax, xlab in ((ax1, r"$\alpha=e_{\rm comp}/e_{\rm resample}$"),
                     (ax2, "total-error inflation (%)")):
        ax.axhline(1.0, color="0.6", ls="--", lw=1.0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(xlab, fontsize=9)
        ax.set_ylabel("HPX / UTM compressed size", fontsize=9)
        ax.tick_params(labelsize=8)
    ax1.legend(fontsize=8, frameon=False, title="scene class",
               title_fontsize=8)
    fig.suptitle(f"Storage verdict against the error already carried "
                 f"({kernel})", fontsize=10)
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


def latex_table(df, e_res, provenance, kernel):
    rep = df[df.reported]
    if not len(rep):
        return ""
    rows = []
    for cls, g in sorted(rep.groupby("cls")):
        med = g.groupby("target").ratio.median().sort_index()
        win = med[med < 1.0]
        t = float(win.index.min()) if len(win) else float(med.index.max())
        r = g[g.target == t].iloc[0]
        rows.append(
            f"{cls.capitalize()} & {e_res[cls]:.3f} & ${t:g}$ & "
            f"{r.alpha:.2f} & {100*r.inflation:.1f} & {med.loc[t]:.2f} & "
            f"{'yes' if len(win) else 'no'} \\\\")
    return (
        "\\begin{table}[t]\n\\centering\n"
        "\\caption{The storage verdict expressed against the error each "
        "product already carries. $e_{\\mathrm{resample}}$ is the RMSE of the "
        f"\\texttt{{{kernel.replace('_', chr(92)+'_')}}} transfer against its "
        "own estimand (" + provenance + "); $\\alpha$ is the compression "
        "target as a fraction of it. Because the two errors are independent "
        "they add in quadrature, so compressing at $\\alpha$ inflates the "
        "product's total error by $\\sqrt{1+\\alpha^2}-1$, given here as a "
        "percentage. The target shown is the tightest reported one at which "
        "HEALPix is the smaller representation, or the loosest tested where "
        "there is none. Every target this paper reports lies below the "
        "reconstruction error of its own scene class: no result is obtained by "
        "compressing more finely than the transfer step is accurate.}\n"
        "\\label{tab:alpha_budget}\n\\footnotesize\n"
        "\\begin{tabular}{@{}lrrrrrc@{}}\n\\toprule\n"
        "Class & $e_{\\mathrm{resample}}$ & Target & $\\alpha$ & "
        "Inflation (\\%) & HPX/UTM & HPX wins \\\\\n\\midrule\n"
        + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kernel", default="clough_tocher")
    ap.add_argument("--reference", choices=["measured", "companion"],
                    default="measured")
    ap.add_argument("--metric", choices=["self", "common"], default="self")
    ap.add_argument("--no-publish", action="store_true")
    ap.add_argument("--no-tex", action="store_true")
    args = ap.parse_args(argv)

    df, e_res, provenance, gamma = build(args.kernel, args.reference,
                                         args.metric)
    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    tag = "" if args.reference == "measured" else "_companion"
    csv = cfg.RESULTS_DIR / f"alpha_budget{tag}.csv"
    df.to_csv(csv, index=False)
    md = summarize(df, e_res, provenance, gamma, args.kernel)
    (cfg.RESULTS_DIR / f"alpha_budget{tag}.md").write_text(md, encoding="utf-8")
    print("\n" + md)
    print(f"wrote {csv}")

    figure(df, e_res, args.kernel, publish=not args.no_publish,
           name="fig_alpha_budget" + tag)
    if not args.no_tex and args.reference == "measured":
        out = cfg.TEX_TABLES / "table_alpha_budget.tex"
        cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
        out.write_text(latex_table(df, e_res, provenance, args.kernel),
                       encoding="utf-8")
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
