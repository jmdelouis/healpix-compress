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
    """Growth of the total error if the two errors were orthogonal.

    This is the rho = 0 case of sqrt(1 + kappa^2 + 2 rho kappa). cross_term.py
    measures rho and finds it neither zero nor of one sign: against the observed
    field it is negative (the codec discards part of the sensor noise the
    transfer faithfully carried, so the errors partly cancel), and against the
    latent scene it is positive (there the transfer error is essentially a
    deficit of high-frequency content -- the blur it did not remove -- and a
    wavelet coder's error is another, so the two align). The published figure is
    therefore the measured inflation, with this one beside it.
    """
    return np.sqrt(1.0 + np.asarray(alpha, dtype=float) ** 2) - 1.0


# The oracle of the main sweep picks a wavelet coder at essentially every
# operating point (JPEG 2000 at 165 of 400, JPEG XL at 150), so the inflation a
# reader would actually incur is the wavelet one. ZFP stays in cross_term.py as
# the control that makes the mechanism a mechanism: it does not sort by
# frequency, its rho is ~0, and it reproduces the orthogonal formula.
INFLATION_CODEC = "jpeg2000_lossy"


def measured_inflation(kappas, reference="latent", codec=INFLATION_CODEC,
                       n_bins=14):
    """Interpolate the measured total-error inflation at these kappa values.

    NaN outside the kappa range the cross-term run actually sampled, rather
    than an extrapolation: a crossing beyond the sweep's reach has no
    measurement, and saying so is the point.
    """
    path = cfg.RESULTS_DIR / "cross_term.csv"
    nan = np.full(np.shape(kappas), np.nan, dtype=float)
    if not path.exists():
        return nan, None
    d = pd.read_csv(path)
    if "reference" in d:
        d = d[d.reference == reference]
    d = d[(d.codec == codec) & (d.kappa > 0)]
    if len(d) < 20:
        return nan, None
    g = (d.assign(_b=pd.cut(np.log10(d.kappa.to_numpy()), n_bins))
         .groupby("_b", observed=True)
         .agg(k=("kappa", "median"), infl=("inflation_measured", "median"),
              n=("kappa", "size")))
    g = g[g.n >= 5].sort_values("k")
    if len(g) < 3:
        return nan, None
    kx, ky = g.k.to_numpy(), np.clip(g.infl.to_numpy(), 1e-9, None)
    q = np.asarray(kappas, dtype=float)
    out = np.exp(np.interp(np.log(q), np.log(kx), np.log(ky)))
    out = np.where((q >= kx.min()) & (q <= kx.max()), out, np.nan)
    return out, (f"{codec} against the {reference} field, "
                 f"{int(d.scene.nunique())} regions")


def load_resampling_error(kernel: str, source: str, denominator: str):
    """e_resampling per scene class, and what it is the error *of*.

    Two denominators are defensible for a non-deconvolving transfer and they
    answer different questions, so neither is chosen silently:

        own    the transfer's error against its own estimand -- for
               Clough--Tocher, the observed field it claims to represent. This
               is the strict reading, and it is small, so a compression target
               is large relative to it.

        blur   the error against the latent scene, which for a transfer that
               does not deconvolve is the blur itself. This is the yardstick if
               the carrier grid is used as if it were the scene, which is what
               a reader treating the HEALPix cube as "the data" does.

    The choice moves alpha by more than an order of magnitude, so the paper has
    to report both rather than pick one.
    """
    if source == "companion":
        return ({c: v["matched"] for c, v in COMPANION.items()},
                "healpix-resample, controlled synthetic (Esri texture), "
                "PSF-aware reconstruction against the latent field")
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
    if denominator == "blur":
        col, what = "e_against_other_estimand", "the latent scene (i.e. the blur)"
    else:
        col, what = "e_resample", f"its own estimand, the {est} field"
    if col not in d:
        raise ValueError(f"{path} has no column {col!r}; re-run "
                         f"resampling_error.py")
    return (d.groupby("cls")[col].median().to_dict(),
            f"measured here on {d.scene.nunique()} regions in "
            f"{d.cls.nunique()} class(es), {kernel} against {what}")


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


def build(kernel="clough_tocher", source="measured", metric="self",
          denominator="own"):
    e_res, provenance = load_resampling_error(kernel, source, denominator)
    ratios, gamma = load_ratios(metric)
    missing = sorted(set(ratios.cls.unique()) - set(e_res))
    if missing:
        print(f"WARNING: no reconstruction error for {', '.join(missing)}; "
              f"those classes are dropped. Re-run resampling_error.py over all "
              f"40 regions before reading anything below as a result.")
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
    above = rep[rep.alpha > 1.0]
    if len(above):
        pairs = sorted({(r.cls, float(r.target)) for _, r in above.iterrows()})
        out.append(f"- **{len(pairs)} of {len(set((r.cls, float(r.target)) for _, r in rep.iterrows()))} "
                   f"(class, target) pairs have alpha > 1**, i.e. the "
                   f"compression error there exceeds the reconstruction error "
                   f"itself: "
                   + ", ".join(f"{c} at {t:g}" for c, t in pairs[:8])
                   + ("…" if len(pairs) > 8 else "")
                   + ". Under this denominator those targets are not "
                     "conservative, and the inflation figures beside them say "
                     "by how much.")
    else:
        out.append("- every reported target sits below the reconstruction "
                   "error of its own class, so no result here is obtained by "
                   "compressing more finely than the transfer step is "
                   "accurate.")
    if len(e_res) < 4:
        out.append(f"- **only {len(e_res)} scene class(es) available**; this is "
                   f"a partial reading, not a result.")
    return "\n".join(out)


def figure(df, e_res, kernel, publish=True, name="fig_alpha_budget"):
    rep = df[df.reported]
    if not len(rep):
        return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.1, 2.9))
    for cls, g in rep.groupby("cls"):
        med = g.groupby(["target", "alpha"]).ratio.median().reset_index()
        med = med.sort_values("alpha")
        c = CLASS_COLOR.get(cls, "0.4")
        ax1.plot(med.alpha, med.ratio, "o-", ms=4, lw=1.6, color=c, label=cls)
        ax2.plot(100 * inflation(med.alpha), med.ratio, "o-", ms=4, lw=1.6,
                 color=c)
    for ax, xlab in ((ax1, r"$\kappa=e_{\rm comp}/e_{\rm resample}$"),
                     (ax2, "total-error inflation (%)")):
        ax.axhline(1.0, color="0.6", ls="--", lw=1.0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(xlab, fontsize=9)
        ax.set_ylabel("HPX / UTM compressed size", fontsize=9)
        ax.tick_params(labelsize=8)
    # The panels are dense at the crossover, which is where the legend used to
    # sit; putting it outside keeps both readable when the float is set at
    # full text width.
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8, frameon=False, ncol=4,
               loc="lower center", bbox_to_anchor=(0.5, -0.06),
               title="scene class", title_fontsize=8)
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


def _tex_escape(t: str) -> str:
    """Underscores in a provenance string are text-mode LaTeX errors."""
    return t.replace("\\", "").replace("_", r"\_").replace("%", r"\%")


def _infl_cell(meas, kap):
    """Measured inflation with the orthogonal value in parentheses."""
    q = 100.0 * float(inflation(kap))
    meas = float(meas)
    return f"--- ({q:.1f})" if meas != meas else f"{100*meas:.1f} ({q:.1f})"


def latex_table(df, e_res, e_other, denominator, provenance, other_prov,
                kernel):
    """One table, both denominators.

    The crossover target depends only on the HPX/UTM ratios, so it is the same
    under either denominator; what changes is what compressing there costs the
    product. Reporting one denominator alone would let the reader mistake a
    choice of yardstick for a result, which is why both columns are here.
    """
    rep = df[df.reported]
    if not len(rep):
        return ""
    n_reg = int(rep.scene.nunique())
    own, blur = ((e_res, e_other) if denominator == "own"
                 else (e_other, e_res))
    rows, worst_own, worst_blur, any_win = [], 0.0, 0.0, False
    over_blur = []
    for cls, g in sorted(rep.groupby("cls")):
        med = g.groupby("target").ratio.median().sort_index()
        win = med[med < 1.0]
        t = float(win.index.min()) if len(win) else float(med.index.max())
        any_win = any_win or bool(len(win))
        k_own = t / own[cls] if cls in own else float("nan")
        k_blur = t / blur[cls] if cls in blur else float("nan")
        worst_own = max(worst_own, 0.0 if k_own != k_own else k_own)
        worst_blur = max(worst_blur, 0.0 if k_blur != k_blur else k_blur)
        if k_blur == k_blur and k_blur > 1.0:
            over_blur.append(cls)
        rows.append(
            f"{cls.capitalize()} & ${t:g}$ & {med.loc[t]:.2f} & "
            f"{'yes' if len(win) else 'no'} & "
            f"{own.get(cls, float('nan')):.4f} & {k_own:.2f} & "
            f"{_infl_cell(measured_inflation([k_own], 'observed')[0][0], k_own)} & "
            f"{blur.get(cls, float('nan')):.4f} & {k_blur:.2f} & "
            f"{_infl_cell(measured_inflation([k_blur], 'latent')[0][0], k_blur)}"
            f" \\\\")
    # The closing sentence is computed, never asserted: under one denominator
    # some crossings do cost more than the transfer error itself.
    lo, hi = min(worst_own, worst_blur), max(worst_own, worst_blur)
    verdict = (f"The choice of denominator moves $\\kappa$ at the crossing by "
               f"a factor of {hi/lo:.0f} ({lo:.2f} against {hi:.2f} at the "
               f"most expensive class), which is why both are reported rather "
               f"than one elected. ")
    if not over_blur:
        verdict += ("Against the latent scene no crossing is obtained by "
                    "compressing more finely than the transfer step is "
                    "accurate.")
    else:
        verdict += ("Against the latent scene the crossing costs less than "
                    "the transfer error itself in "
                    f"{len(rows) - len(over_blur)} of {len(rows)} classes; "
                    f"the exception{'s are' if len(over_blur) > 1 else ' is'} "
                    + ", ".join(sorted(over_blur))
                    + ", whose inflation is given rather than glossed.")
    return (
        "\\begin{table*}[t]\n\\centering\n"
        "\\caption{The storage verdict expressed against the error each "
        "product already carries, under both defensible denominators. The "
        f"transfer is \\texttt{{{kernel.replace('_', chr(92)+'_')}}}, which does "
        f"not deconvolve, measured on the {n_reg} regions of "
        "Table~\\ref{tab:benchmark_scenes}. \\emph{Own estimand}: its RMSE "
        "against the observed field it claims to represent. "
        "\\emph{Latent scene}: its RMSE against the unblurred scene, which for "
        "a non-deconvolving transfer is the instrument blur it does not remove "
        "-- the yardstick that applies when the carrier grid is read as if it "
        "were the scene. $\\kappa$ is the compression target as a fraction of "
        "the denominator; the two errors being independent, they add in "
        "quadrature only if their cross term vanishes, and it does not: "
        "Appendix~\\ref{app:cross_term} measures it negative against the "
        "observed field and positive against the latent scene. The inflation "
        "columns are therefore the \\emph{measured} growth of the product's "
        "total error, with the orthogonal value $\\sqrt{1+\\kappa^2}-1$ in "
        "parentheses; a dash marks a crossing beyond the range the cross-term "
        "run sampled. The target shown is "
        "the tightest reported one at which HEALPix is the smaller "
        "representation, or the loosest tested where there is none; it is set "
        "by the ratios alone and so is common to both denominators. "
        + verdict + "}\n"
        "\\label{tab:alpha_budget}\n\\footnotesize\n"
        "\\begin{tabular}{@{}lrrc rrr rrr@{}}\n\\toprule\n"
        " & \\multicolumn{3}{c}{Crossover} & "
        "\\multicolumn{3}{c}{Own estimand} & "
        "\\multicolumn{3}{c}{Latent scene} \\\\\n"
        "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}\n"
        "Class & Target & HPX/UTM & HPX wins & "
        "$e_{\\mathrm{resample}}$ & $\\kappa$ & Infl.\\ \\% m.\\ (o.) & "
        "$e_{\\mathrm{resample}}$ & $\\kappa$ & Infl.\\ \\% m.\\ (o.) \\\\\n"
        "\\midrule\n"
        + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n\\end{table*}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kernel", default="clough_tocher")
    ap.add_argument("--reference", choices=["measured", "companion"],
                    default="measured")
    ap.add_argument("--metric", choices=["self", "common"], default="self")
    ap.add_argument("--denominator", choices=["own", "blur"], default="own",
                    help="what the compression error is measured against: the "
                         "transfer's own estimand ('own'), or the latent scene "
                         "('blur'), which for a non-deconvolving transfer is "
                         "the blur it does not remove. They differ by more "
                         "than an order of magnitude; report both")
    ap.add_argument("--no-publish", action="store_true")
    ap.add_argument("--no-tex", action="store_true")
    args = ap.parse_args(argv)

    df, e_res, provenance, gamma = build(args.kernel, args.reference,
                                         args.metric, args.denominator)
    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    tag = ("" if args.reference == "measured" else "_companion")
    tag += "" if args.denominator == "own" else f"_{args.denominator}"
    csv = cfg.RESULTS_DIR / f"alpha_budget{tag}.csv"
    df.to_csv(csv, index=False)
    md = summarize(df, e_res, provenance, gamma, args.kernel)
    (cfg.RESULTS_DIR / f"alpha_budget{tag}.md").write_text(md, encoding="utf-8")
    print("\n" + md)
    print(f"wrote {csv}")

    figure(df, e_res, args.kernel, publish=not args.no_publish,
           name="fig_alpha_budget" + tag)
    if not args.no_tex and args.reference == "measured":
        # The published table carries both denominators, so it is identical
        # whichever one was asked for on the command line and neither run can
        # silently overwrite the other with half the story.
        other = "blur" if args.denominator == "own" else "own"
        e_other, other_prov = load_resampling_error(args.kernel,
                                                    args.reference, other)
        out = cfg.TEX_TABLES / "table_alpha_budget.tex"
        cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
        out.write_text(latex_table(df, e_res, e_other, args.denominator,
                                   provenance, other_prov, args.kernel),
                       encoding="utf-8")
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
