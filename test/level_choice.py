"""Is level 20 the right carrier level, or is level 19 enough? (idea 1)

The paper reaches for level 20 on a nominal argument: level 19's \\SI{12.4}{m}
cell is coarser than the \\SI{10}{m} product, so a level-20 carrier is used to
avoid undersampling it, at a cost of $2.58\\times$ in cells. That is a statement
about resolution, not about information, and the paper says as much -- it notes
that claiming no detail is lost would need a transfer-error analysis it does not
provide. This is that analysis.

A spectral argument makes the question look live and gets it wrong, which is
why the measurement here is direct. The observed field carries a median
$0.21\\%$ of its variance above half the native Nyquist and, propagating the
fitted slope and response width, about $0.03\\%$ above the level-19 Nyquist --
from which one would conclude that level 19 discards three parts in ten
thousand. It discards $3.07\\%$, a hundred times more. Aggregation is not a
brick wall at the coarser Nyquist: the residual of a $2\\times2$ mean has
response $1-\\cos^2(k_xd/2)\\cos^2(k_yd/2)$, which grows like $(kd)^4$ and so
feeds mainly on the shoulder *below* the native Nyquist. Variance above a
Nyquist is not the error of sampling below it.

Variance above a Nyquist is not the error of sampling below it, though, so the
measurement here is direct. Two quantities are separated, because a coarser
level is charged twice and only one of the charges is about information:

    discarded   the RMS of what aggregating the stored level-20 tile to a
                coarser level throws away, i.e. the intra-cell variance. This is
                the information question, and it is exact: HEALPix cells being
                equal-area, a parent's value *is* the mean of its four children,
                so the aggregate is the coarser carrier rather than a model of
                it.

    floor       the RMSE of that coarser carrier, read at its nearest containing
                cell, against the same Sentinel-2 observation over the same
                mask. This is the operational number, and it also pays for a
                coarser read -- at level 19 the read jitters over a
                \\SI{12.4}{m} cell instead of a \\SI{6.2}{m} one.

The two charges separate cleanly in the result: the coarser floor is the level-20
floor and the discarded energy in quadrature (median measured/quadrature $1.07$
at level 19), so the excess is information the aggregation removes, not a worse
read. The decision does not follow from the raw ratio, which never appears in a
compressed pipeline -- at the level-20 floor the coded file is $2.24\\%$ of the
raw tile. What the per-level energies are good for is the multiresolution coder
of idea 2: they are the subband energies of the exact HEALPix Haar pyramid.

Usage
-----
    python level_choice.py --limit 4
    python level_choice.py
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
import build_tiles                        # noqa: E402
import common_reference                   # noqa: E402
import region_sites                       # noqa: E402

DEPTHS = (0, 1, 2, 3)          # levels 20, 19, 18, 17
CLASS_COLOR = {"urban": "#D55E00", "water": "#0072B2",
               "forest": "#009E73", "agriculture": "#E69F00"}


def rms(a):
    a = np.asarray(a, dtype=np.float64)
    a = a[np.isfinite(a)]
    return float(np.sqrt(np.mean(a ** 2))) if a.size else np.nan


def region_rows(scene, mode="nearest", verbose=True):
    cell_ids, cell_data, utm_crop, meta, npz = build_tiles.load_tile(scene)
    cr = common_reference.build_for_tile(npz, mode=mode)
    obs = np.asarray(utm_crop, dtype=np.float64)
    flat = np.asarray(cell_data, dtype=np.float64).ravel()
    sd = float(np.std(flat)) or np.nan
    rows = []
    for d in DEPTHS:
        level = cfg.LEVEL - d
        agg, prolonged = common_reference.CommonReference.aggregate(flat, d)
        read = cr.read_at_level(flat, d).reshape(cr.shape)
        floor = cr.rmse_on_mask(read, obs)
        discarded = rms(flat - prolonged)
        # Raw storage, the geometric constant this experiment is about.
        hpx_bytes = int(agg.size * 4)
        utm_bytes = int(np.asarray(utm_crop).nbytes)
        rows.append(dict(
            scene=scene, cls=scene.split("__")[1], level=level, depth=d,
            n_cells=int(agg.size), cell_m=cfg.DELTA_EQ * (2 ** d),
            floor=floor, floor_frac_of_std=floor / sd if sd == sd else np.nan,
            discarded=discarded,
            discarded_frac_of_std=discarded / sd if sd == sd else np.nan,
            discarded_var_frac=(discarded ** 2 / sd ** 2) if sd == sd else np.nan,
            hpx_raw_bytes=hpx_bytes, utm_raw_bytes=utm_bytes,
            raw_ratio=hpx_bytes / utm_bytes,
            raw_ratio_gamma=hpx_bytes / (utm_bytes * cfg.GAMMA_UTM),
            tile_std=sd, mask_fraction=cr.mask_fraction))
        if verbose:
            r = rows[-1]
            print(f"    level {level}  cell {r['cell_m']:5.2f} m  "
                  f"floor {floor:.4e} ({100*r['floor_frac_of_std']:5.1f}% of std)  "
                  f"discarded {100*r['discarded_frac_of_std']:6.3f}% of std "
                  f"({100*r['discarded_var_frac']:.4f}% of variance)  "
                  f"raw HPX/UTM {r['raw_ratio']:.3f}")
    return rows


def summarize(df):
    out = ["# Which carrier level? Level 20 against coarser ones", "",
           f"- regions: {df.scene.nunique()}; the read is the nearest "
           f"containing cell at every level", ""]
    out += ["| level | cell (m) | cells | raw HPX/UTM | floor | floor % of std "
            "| discarded % of variance |", "|---|---|---|---|---|---|---|"]
    for lvl, g in df.groupby("level", sort=False):
        out.append(f"| {lvl} | {g.cell_m.iloc[0]:.2f} | {g.n_cells.iloc[0]} | "
                   f"{g.raw_ratio.iloc[0]:.3f} | {g.floor.median():.3e} | "
                   f"{100*g.floor_frac_of_std.median():.1f}% | "
                   f"{100*g.discarded_var_frac.median():.4f}% |")
    out.append("")
    top = df[df.depth == 0]
    one = df[df.depth == 1]
    if len(top) and len(one):
        f20, f19 = top.floor.median(), one.floor.median()
        out += [f"- level 20 floor {f20:.3e}, level 19 floor {f19:.3e}: "
                f"a factor of {f19/f20:.2f}",
                f"- aggregating to level 19 discards "
                f"{100*one.discarded_var_frac.median():.4f}% of the tile's "
                f"variance. This is *not* the power above the level-19 "
                f"Nyquist: a $2\\times2$ mean is a high-pass whose response "
                f"$1-\\cos^2(k_xd/2)\\cos^2(k_yd/2)$ grows like $(kd)^4$, so "
                f"it feeds mainly on the shoulder below the native Nyquist",
                f"- raw storage: level 20 is {top.raw_ratio.iloc[0]:.3f} of the "
                f"equal-area UTM crop, level 19 is {one.raw_ratio.iloc[0]:.3f}",
                ""]
        # The decision this measurement is for.
        if f19 / f20 < 1.2:
            out.append("- **Level 19 is the better carrier on this evidence**: "
                       "it costs a quarter of the cells and its floor is within "
                       "20% of level 20's, so the $2.58\\times$ the paper "
                       "charges for oversampling buys almost nothing against "
                       "the observation.")
        else:
            # Is the extra error the information the aggregation removed, or
            # a coarser read on top of it?  Quadrature settles it per scene.
            a = top.set_index("scene")
            b = one.set_index("scene")
            j = b.join(a[["floor"]].rename(columns={"floor": "floor20"}),
                       how="inner").dropna(subset=["floor", "floor20",
                                                   "discarded"])
            quad = (j.floor20 ** 2 + j.discarded ** 2) ** 0.5
            over = float((j.floor / quad).median())
            excess = j.floor ** 2 - j.floor20 ** 2
            expl = float((j.discarded ** 2 / excess).median())
            out.append(f"- **The level-19 floor is {f19/f20:.2f} times level "
                       f"20's**, so a coarser carrier is not free against the "
                       f"observation. That excess is the information the "
                       f"aggregation removes, not a worse read: the level-19 "
                       f"floor is the level-20 floor and the discarded "
                       f"intra-cell energy in quadrature to within "
                       f"{100*abs(over-1):.0f}% (median measured/quadrature "
                       f"{over:.3f}), and the discarded energy alone accounts "
                       f"for {100*expl:.0f}% of the excess variance.")
    out += ["", "## By scene class, level 19 against level 20", "",
            "| class | floor 20 | floor 19 | ratio | discarded % of variance |",
            "|---|---|---|---|---|"]
    for c in sorted(df.cls.unique()):
        a = df[(df.cls == c) & (df.depth == 0)].floor.median()
        b = df[(df.cls == c) & (df.depth == 1)]
        out.append(f"| {c} | {a:.3e} | {b.floor.median():.3e} | "
                   f"{b.floor.median()/a:.2f} | "
                   f"{100*b.discarded_var_frac.median():.4f}% |")
    return "\n".join(out)


def figure(df, publish=True, name="fig_level_choice"):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.1, 2.9))
    for c, g in df.groupby("cls"):
        m = g.groupby("level")[["floor", "discarded"]].median().sort_index()
        col = CLASS_COLOR.get(c, "0.4")
        ax1.plot(m.index, m.floor, "o-", ms=4, lw=1.6, color=col, label=c)
        ax2.plot(m.index, m.discarded, "o-", ms=4, lw=1.6, color=col, label=c)
    for ax, yl in ((ax1, "floor: RMSE vs the observation"),
                   (ax2, "discarded by aggregation (RMS)")):
        ax.set_yscale("log")
        ax.set_xlabel("carrier level", fontsize=9)
        ax.set_ylabel(yl, fontsize=9)
        ax.invert_xaxis()
        ax.tick_params(labelsize=8)
    ax1.legend(fontsize=8, frameon=False)
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
    for lvl, g in df.groupby("level", sort=False):
        body.append(
            f"{lvl} & {g.cell_m.iloc[0]:.2f} & {g.n_cells.iloc[0]:,} & "
            f"{g.raw_ratio.iloc[0]:.3f} & {g.floor.median():.2e} & "
            f"{100*g.floor_frac_of_std.median():.1f} & "
            f"{100*g.discarded_var_frac.median():.4f} \\\\".replace(",", "\\,"))
    return (
        "\\begin{table}[t]\n\\centering\n"
        "\\caption{What a coarser carrier level costs. The level-$L$ carrier is "
        "obtained from the stored level-20 tile by exact aggregation --- "
        "HEALPix cells being equal-area, a parent's value is the mean of its "
        "four children --- so each row is the carrier that transfer would have "
        "produced, not a model of it. \\emph{Floor} is the RMSE of that carrier, "
        "read at its nearest containing cell, against the same Sentinel-2 "
        "observation over the two footprints' intersection; it charges a coarser "
        "level twice, once for what the aggregation discards and once for a "
        "coarser read. \\emph{Discarded} is the intra-cell variance the "
        "aggregation removes, as a percentage of the tile's variance, and "
        "isolates the information question from the read. \\emph{Raw HPX/UTM} "
        "is the uncompressed byte ratio against the equal-area native crop, a "
        "geometric constant.}\n"
        "\\label{tab:level_choice}\n\\footnotesize\n"
        "\\begin{tabular}{@{}rrrrrrr@{}}\n\\toprule\n"
        "Level & Cell (m) & Cells & Raw HPX/UTM & Floor & Floor \\% of std & "
        "Discarded \\% of var. \\\\\n\\midrule\n"
        + "\n".join(body) + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--mode", choices=["nearest", "knn", "psf"],
                    default="nearest",
                    help="the read. Only 'nearest' is defined at every level "
                         "without rebuilding an operator, and it is the read "
                         "the paper's floor uses")
    ap.add_argument("--no-publish", action="store_true")
    ap.add_argument("--no-tex", action="store_true")
    args = ap.parse_args(argv)

    scenes = args.regions or list(region_sites.load_sites()["scene"])
    scenes = [s for s in scenes if build_tiles.tile_path(s).exists()]
    if args.limit:
        scenes = scenes[:args.limit]
    if not scenes:
        print(f"no tile under {cfg.TILE_DIR}; run build_tiles.py first")
        return 1

    print(f"{len(scenes)} region(s), levels "
          f"{[cfg.LEVEL - d for d in DEPTHS]}, read {args.mode!r}")
    rows = []
    for i, scene in enumerate(scenes, 1):
        print(f"[{i:>2}/{len(scenes)}] {scene}")
        try:
            rows += region_rows(scene, mode=args.mode)
        except Exception as exc:                          # noqa: BLE001
            print(f"    FAILED {type(exc).__name__}: {exc}")
    if not rows:
        print("no region produced a row")
        return 1

    df = pd.DataFrame(rows)
    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    csv = cfg.RESULTS_DIR / "level_choice.csv"
    df.to_csv(csv, index=False)
    md = summarize(df)
    (cfg.RESULTS_DIR / "level_choice.md").write_text(md, encoding="utf-8")
    print("\n" + md)
    print(f"\nwrote {csv} ({len(df)} rows)")
    figure(df, publish=not args.no_publish)
    if not args.no_tex:
        cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
        out = cfg.TEX_TABLES / "table_level_choice.tex"
        out.write_text(latex_table(df), encoding="utf-8")
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
