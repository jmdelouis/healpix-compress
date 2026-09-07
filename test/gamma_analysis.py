"""Report the verdict as a function of the operational overhead (V4).

The paper's headline comparison multiplies the native UTM size by an assumed
granule overhead of 30% and then reports a verdict. The review's objection is
that picking a value makes the conclusion an artefact of that choice, and that
the overhead is applied asymmetrically: the UTM side is charged for padding,
overlap and per-granule metadata while the HEALPix side is charged for nothing.

Both are answerable without re-running a single sweep, because everything the
comparison says about an overhead is already contained in one per-region
quantity. Write the intrinsic (overhead-free) ratio of a region as r, and its
critical overhead as

    gamma_star = r,

the overhead at which the two representations reach equal compressed size. At
any operational overhead gamma, that region's operational ratio is
gamma_star/gamma, HEALPix is the smaller of the two exactly when
gamma_star <= gamma, and the paired test across regions at overhead gamma is
the signed-rank test on log(gamma_star/gamma). So the whole verdict, its
significance and its spread are functions of the empirical distribution of
gamma_star, and reporting them over a continuous range of gamma costs nothing
and commits to nothing.

The asymmetry is answered by reading gamma as a *ratio* of overheads rather
than as a property of UTM alone: with an overhead gamma_hpx on the carrier side
and gamma_utm on the native side, the operational ratio is
(S_hpx gamma_hpx)/(S_utm gamma_utm), so the single parameter that decides the
verdict is gamma = gamma_utm/gamma_hpx. Every number below is a statement about
that ratio. `--containers` additionally measures the per-chunk container cost
actually paid on each side, which bounds the part of gamma that comes from
container format rather than from granule geometry.

Usage
-----
    python gamma_analysis.py                    # self metric
    python gamma_analysis.py --metric common
    python gamma_analysis.py --containers       # + measured container bytes

Writes results/gamma_analysis_<metric>.{csv,md},
figures/fig_verdict_vs_gamma.pdf and ../tex/tables/table_gamma_continuous.tex.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt           # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402

import paper_config as cfg                # noqa: E402
import paper_stats as st                  # noqa: E402

RMSE_LABEL = {1e-4: r"$10^{-4}$", 3e-4: r"$3\times10^{-4}$", 1e-3: r"$10^{-3}$",
              3e-3: r"$3\times10^{-3}$", 1e-2: r"$10^{-2}$"}
PLAIN_LABEL = {1e-4: "1e-4", 3e-4: "3e-4", 1e-3: "1e-3",
               3e-3: "3e-3", 1e-2: "1e-2"}
# Report gamma over a range that brackets every plausible operational layout.
# It starts below 1: gamma is a *ratio* of the two sides' overheads, so
# gamma < 1 is the case where the carrier grid carries the heavier container,
# and a fair share of regions have a critical overhead there -- starting the
# axis at 1 would hide them and make the win fractions unreadable.
GAMMA_GRID = np.linspace(0.8, 2.6, 361)
REFERENCE_GAMMAS = [1.0, 1.15, 1.30, 1.50, 2.00]


def load_targets(metric: str):
    """Per-region critical overheads, for every reported RMSE target."""
    path = cfg.RESULTS_DIR / f"paper_numbers_{metric}.json"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing; run update_paper_numbers.py"
            + ("" if metric == "self" else f" --metric {metric}"))
    doc = json.loads(path.read_text(encoding="utf-8"))
    gamma_utm = float(doc["gamma_utm"])
    out = []
    for t in doc["targets"]:
        if t["kind"] != "rmse":
            continue
        scenes = sorted(t["ratios"])
        # The stored ratios are operational, i.e. already divided by the
        # assumed overhead; multiplying it back gives each region's critical
        # overhead, which is the quantity everything below is a function of.
        gs = np.array([t["ratios"][s] for s in scenes], dtype=float) * gamma_utm
        ok = np.isfinite(gs) & (gs > 0)
        out.append(dict(target=float(t["target"]), scenes=list(np.array(scenes)[ok]),
                        gamma_star=gs[ok],
                        bracketed_fraction=float(t.get("bracketed_fraction", np.nan)),
                        reported=bool(t.get("bracketed_fraction", 0)
                                      >= cfg.MIN_BRACKETED)))
    return out, gamma_utm, doc


def curves(gamma_star, grid=GAMMA_GRID):
    """Win fraction, median operational ratio and paired p-value, over gamma."""
    gs = np.asarray(gamma_star, dtype=float)
    win = np.array([float(np.mean(gs <= g)) for g in grid])
    med = np.array([float(np.median(gs)) / g for g in grid])
    p = np.array([st.paired_ratio_test(gs / g)["p"] for g in grid])
    return win, med, p


def inconclusive_interval(grid, p, alpha=0.05):
    """The gamma range over which the paired test does not reject parity."""
    m = np.isfinite(p) & (p > alpha)
    if not m.any():
        return (np.nan, np.nan)
    return float(grid[m].min()), float(grid[m].max())


def analyse(metric="self"):
    targets, gamma_utm, doc = load_targets(metric)
    rows = []
    for t in targets:
        gs = t["gamma_star"]
        if gs.size < cfg.MIN_SCENES:
            continue
        win, med, p = curves(gs)
        lo, hi = inconclusive_interval(GAMMA_GRID, p)
        stats = st.paired_ratio_test(gs)          # test of gamma_star vs 1
        row = dict(metric=metric, target=t["target"], n=int(gs.size),
                   reported=t["reported"],
                   bracketed_fraction=t["bracketed_fraction"],
                   gamma_star_median=float(np.median(gs)),
                   gamma_star_q25=float(np.percentile(gs, 25)),
                   gamma_star_q75=float(np.percentile(gs, 75)),
                   gamma_star_ci_lo=stats.get("ci_lo", np.nan),
                   gamma_star_ci_hi=stats.get("ci_hi", np.nan),
                   inconclusive_lo=lo, inconclusive_hi=hi)
        for g in REFERENCE_GAMMAS:
            row[f"win_at_{g:g}"] = float(np.mean(gs <= g))
            row[f"median_ratio_at_{g:g}"] = float(np.median(gs)) / g
        rows.append(dict(row, _gamma_star=gs, _win=win, _med=med, _p=p))
    return pd.DataFrame(rows), gamma_utm, doc


# --- figure ------------------------------------------------------------------
def figure(df, gamma_utm, publish=True, name="fig_verdict_vs_gamma"):
    d = df[df.reported]
    if not len(d):
        print("  no reported target; figure skipped")
        return
    cmap = plt.get_cmap("viridis")
    colors = {t: cmap(i / max(len(d) - 1, 1))
              for i, t in enumerate(sorted(d.target))}

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(6.4, 6.2), sharex=True,
        gridspec_kw=dict(height_ratios=[2.0, 1.35], hspace=0.12))

    for _, r in d.sort_values("target").iterrows():
        c = colors[r.target]
        ax1.plot(GAMMA_GRID, 100 * r._win, color=c, lw=1.8,
                 label=f"RMSE {PLAIN_LABEL.get(r.target, f'{r.target:g}')}")
        # where the paired test cannot separate the two grids
        m = np.isfinite(r._p) & (r._p > 0.05)
        if m.any():
            ax1.plot(GAMMA_GRID[m], 100 * np.asarray(r._win)[m], color=c,
                     lw=7.0, alpha=0.30, solid_capstyle="butt", zorder=0)

    ax1.axvline(gamma_utm, color="0.35", ls="--", lw=1.1)
    ax1.annotate(f"illustrative\n$\\gamma={gamma_utm:g}$", xy=(gamma_utm, 96),
                 xytext=(gamma_utm + 0.04, 96), fontsize=8, color="0.3",
                 va="top")
    ax1.axhline(50, color="0.75", lw=0.8, zorder=-1)
    ax1.axvline(1.0, color="0.75", lw=0.8, zorder=-1)
    ax1.set_ylabel("regions where HEALPix is smaller (%)")
    ax1.set_ylim(-3, 103)
    ax1.legend(fontsize=8, loc="upper left", frameon=False, ncol=2,
               title="target", title_fontsize=8, bbox_to_anchor=(0.02, 0.98))
    ax1.set_title("Verdict as a function of the overhead ratio "
                  r"$\gamma=\gamma_{\mathrm{UTM}}/\gamma_{\mathrm{HPX}}$",
                  fontsize=10)

    order = list(sorted(d.target))
    for i, t in enumerate(order):
        gs = d[d.target == t].iloc[0]._gamma_star
        c = colors[t]
        ax2.scatter(gs, np.full(gs.size, i) + np.random.default_rng(0)
                    .uniform(-0.14, 0.14, gs.size), s=7, color=c, alpha=0.55,
                    linewidths=0)
        ax2.plot([np.median(gs)], [i], marker="|", ms=16, mew=2.0, color=c)
    ax2.axvline(gamma_utm, color="0.35", ls="--", lw=1.1)
    ax2.set_yticks(range(len(order)))
    ax2.set_yticklabels([PLAIN_LABEL.get(t, f"{t:g}") for t in order], fontsize=8)
    ax2.set_ylabel("target RMSE")
    ax2.set_xlabel(r"overhead ratio $\gamma$")
    ax2.set_ylim(-0.6, len(order) - 0.4)
    ax2.set_xlim(GAMMA_GRID.min(), GAMMA_GRID.max())
    ax2.set_title(r"per-region critical overhead $\gamma_\star$ "
                  "(bar: median)", fontsize=9)

    cfg.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(cfg.FIGURES_DIR / f"{name}.{ext}", bbox_inches="tight", dpi=160)
    if publish:
        cfg.TEX_FIGURES.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(cfg.FIGURES_DIR / f"{name}.pdf",
                        cfg.TEX_FIGURES / f"{name}.pdf")
    plt.close(fig)
    print(f"  {name}")


# --- measured container cost on each side ------------------------------------
def container_costs(scene=None):
    """Bytes each side pays for its container, measured rather than assumed.

    This is the part of an overhead that comes from the storage format, and it
    is the part that can be settled here. It is *not* the 30% of the paper's
    illustrative scenario, which is granule geometry -- swath padding,
    inter-granule overlap, per-granule metadata -- and is a property of a
    particular product layout rather than of a format. Separating the two is
    the point: what is measurable is measured, what is assumed stays a
    parameter.
    """
    import io
    out = {}

    # HEALPix side: one level-20 chunk in a zarr v3 store.
    try:
        import zarr
        store = zarr.storage.MemoryStore()
        root = zarr.api.synchronous.create_group(store, zarr_format=3)
        data = np.zeros(cfg.N_CHILDREN, dtype=np.float32)
        root.create_array("cells", data=data, chunks=(cfg.N_CHILDREN,),
                          fill_value=float("nan"), dimension_names=["cells"])
        meta = sum(len(v) for k, v in store._store_dict.items()
                   if str(k).endswith("zarr.json"))
        out["hpx_container_bytes"] = int(meta)
        out["hpx_container_note"] = (
            "zarr v3 group and array metadata for one level-20 chunk; the "
            "cell identifiers are implicit for a contiguous NESTED block, so "
            "no index is stored")
    except Exception as exc:                                # noqa: BLE001
        out["hpx_container_bytes"] = np.nan
        out["hpx_container_note"] = f"not measured ({type(exc).__name__})"

    # Native side: the header of a single-tile GeoTIFF of the same crop size.
    try:
        import rasterio
        from rasterio.io import MemoryFile
        n = cfg.UTM_SIDE_PX
        with MemoryFile() as mem:
            with mem.open(driver="GTiff", width=n, height=n, count=1,
                          dtype="float32", crs="EPSG:32631",
                          transform=rasterio.transform.from_origin(
                              0, 0, cfg.NATIVE_PIXEL_M, cfg.NATIVE_PIXEL_M),
                          compress="deflate") as dst:
                dst.write(np.zeros((n, n), dtype=np.float32), 1)
            total = len(mem.read())
        # an all-zero float32 plane compresses to a few hundred bytes, so what
        # is left is essentially the header and the tag directory
        out["utm_container_bytes"] = int(total)
        out["utm_container_note"] = (
            "complete single-tile deflate GeoTIFF of an all-zero "
            f"{cfg.UTM_SIDE_PX}x{cfg.UTM_SIDE_PX} float32 plane, i.e. header, "
            "tag directory and a near-empty payload")
    except Exception as exc:                                # noqa: BLE001
        out["utm_container_bytes"] = np.nan
        out["utm_container_note"] = f"not measured ({type(exc).__name__})"

    for side, raw in (("hpx", cfg.SHARED_RAW), ("utm", cfg.UTM_RAW)):
        b = out.get(f"{side}_container_bytes", np.nan)
        out[f"{side}_container_fraction_of_raw"] = (b / raw if np.isfinite(b)
                                                    else np.nan)
    return out


# --- reporting ---------------------------------------------------------------
def latex_table(df, gamma_utm):
    d = df[df.reported].sort_values("target")
    if not len(d):
        return ""
    body = []
    for _, r in d.iterrows():
        label = RMSE_LABEL.get(r.target, f"${r.target:g}$")
        inc = ("--" if not np.isfinite(r.inconclusive_lo)
               else f"${r.inconclusive_lo:.2f}$--${r.inconclusive_hi:.2f}$")
        body.append(
            f"{label} & {int(r.n)} & {r.gamma_star_median:.2f} & "
            f"[{r.gamma_star_q25:.2f}, {r.gamma_star_q75:.2f}] & "
            f"{100*r['win_at_1']:.0f} & {100*r['win_at_1.3']:.0f} & "
            f"{100*r['win_at_2']:.0f} & {inc} \\\\")
    return (
        "\\begin{table*}[t]\n\\centering\n"
        "\\caption{The verdict over a continuous range of the overhead ratio "
        "$\\gamma=\\gamma_{\\mathrm{UTM}}/\\gamma_{\\mathrm{HPX}}$, rather than "
        "at one assumed value. Each region's critical overhead "
        "$\\gamma_\\star$ is the value at which the two representations reach "
        "equal compressed size for that region, so HEALPix is the smaller of "
        "the two exactly where $\\gamma\\ge\\gamma_\\star$; the win columns are "
        "the fraction of the 40 regions satisfying that at three illustrative "
        "values. \\emph{Inconclusive} is the range of $\\gamma$ over which the "
        "paired signed-rank test does not separate the two grids at the $5\\%$ "
        "level. No value of $\\gamma$ is privileged: the table is a map from an "
        "assumption a reader may hold to the conclusion it implies.}\n"
        "\\label{tab:gamma_continuous}\n\\footnotesize\n"
        "\\begin{tabular}{lrrlrrrl}\n\\toprule\n"
        "Target RMSE & $n$ & $\\gamma_\\star$ & IQR & \\multicolumn{3}{c}{HPX "
        "smaller (\\% of regions)} & Inconclusive \\\\\n"
        "\\cmidrule(lr){5-7}\n"
        " & & median & & $\\gamma{=}1$ & $\\gamma{=}1.3$ & $\\gamma{=}2$ & "
        "$\\gamma$ range \\\\\n\\midrule\n"
        + "\n".join(body) + "\n\\bottomrule\n\\end{tabular}\n\\end{table*}\n")


def summarize(df, gamma_utm, containers=None):
    d = df[df.reported].sort_values("target")
    out = ["# Verdict as a function of the overhead (V4)", "",
           "Every column below is a function of the per-region critical "
           "overhead, so nothing here required re-running a sweep.", "",
           "| target | n | gamma* median | IQR | HPX wins at g=1 | at 1.3 | "
           "at 2.0 | inconclusive range |",
           "|---|---|---|---|---|---|---|---|"]
    for _, r in d.iterrows():
        inc = ("--" if not np.isfinite(r.inconclusive_lo)
               else f"{r.inconclusive_lo:.2f}-{r.inconclusive_hi:.2f}")
        out.append(f"| {r.target:g} | {int(r.n)} | {r.gamma_star_median:.2f} | "
                   f"[{r.gamma_star_q25:.2f}, {r.gamma_star_q75:.2f}] | "
                   f"{100*r['win_at_1']:.0f}% | {100*r['win_at_1.3']:.0f}% | "
                   f"{100*r['win_at_2']:.0f}% | {inc} |")
    out.append("")
    lo = d.gamma_star_median.min()
    hi = d.gamma_star_median.max()
    out.append(f"- median critical overhead spans **{lo:.2f} to {hi:.2f}** "
               f"across the reported targets")
    never = d[d["win_at_2"] < 0.5]
    if len(never):
        out.append("- targets where HEALPix is still behind for a majority of "
                   "regions even at an overhead ratio of 2.0: "
                   + ", ".join(f"{t:g}" for t in never.target))
    out.append("")
    if containers:
        out += ["## Measured container cost on each side", ""]
        for side in ("hpx", "utm"):
            b = containers.get(f"{side}_container_bytes", np.nan)
            f = containers.get(f"{side}_container_fraction_of_raw", np.nan)
            out.append(f"- **{side.upper()}**: "
                       + (f"{int(b)} B, {100*f:.4f}% of the raw tile"
                          if np.isfinite(b) else "not measured")
                       + f" — {containers.get(f'{side}_container_note','')}")
        both = all(np.isfinite(containers.get(f"{s_}_container_bytes", np.nan))
                    for s_ in ("hpx", "utm"))
        worst = max((containers.get(f"{s_}_container_fraction_of_raw", 0.0) or 0.0)
                    for s_ in ("hpx", "utm")) if both else np.nan
        # Only claim what was actually measured: an unmeasured side cannot
        # support a conclusion about container cost being negligible.
        out.append(
            (f"Container format therefore accounts for at most "
             f"{100*worst:.4f}% on either side, a negligible part of any "
             f"plausible gamma."
             if both else
             "One or both sides could not be measured in this environment, so "
             "no conclusion about container cost is drawn here; install zarr "
             "and rasterio and re-run to obtain it.")
            + " The illustrative 30% is granule geometry -- swath padding, "
              "inter-granule overlap, per-granule metadata -- which is a "
              "property of a product layout and not of a format, and is why it "
              "is reported as a parameter rather than measured here.")
        out.append("")
    return "\n".join(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metric", choices=["self", "common"], default="self")
    ap.add_argument("--containers", action="store_true",
                    help="also measure the per-chunk container cost on each side")
    ap.add_argument("--no-publish", action="store_true")
    ap.add_argument("--no-tex", action="store_true")
    args = ap.parse_args(argv)

    df, gamma_utm, _doc = analyse(args.metric)
    if not len(df):
        print("no target with enough regions")
        return 1

    cont = container_costs() if args.containers else None

    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    flat = df.drop(columns=[c for c in df.columns if c.startswith("_")])
    csv = cfg.RESULTS_DIR / f"gamma_analysis_{args.metric}.csv"
    flat.to_csv(csv, index=False)
    md = summarize(df, gamma_utm, cont)
    (cfg.RESULTS_DIR / f"gamma_analysis_{args.metric}.md").write_text(
        md, encoding="utf-8")
    print("\n" + md)
    print(f"wrote {csv}")

    figure(df, gamma_utm, publish=not args.no_publish,
           name="fig_verdict_vs_gamma" + ("" if args.metric == "self"
                                          else f"_{args.metric}"))

    if not args.no_tex:
        tex = latex_table(df, gamma_utm)
        if tex:
            tag = "" if args.metric == "self" else f"_{args.metric}"
            out = cfg.TEX_TABLES / f"table_gamma_continuous{tag}.tex"
            cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
            out.write_text(tex, encoding="utf-8")
            print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
