"""Fig. 1: the study workflow, from native UTM to HEALPix carrier and back.

Drawn with matplotlib at column width so that it is typeset at the size it
is drawn (revision 5, review point R1). No data enters this figure; it is a
diagram of the chain the paper measures, and the caption in the text carries
the cross-references.

Usage
-----
    python fig_workflow.py            # writes figures/fig_workflow.pdf and copies it to tex/figures
    python fig_workflow.py --no-publish
"""

from __future__ import annotations

import argparse
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch   # noqa: E402

import paper_config as cfg                            # noqa: E402

COLUMN_IN = 3.5
INK = "#222222"
GRID = "#0072B2"      # native UTM side
CARRIER = "#D55E00"   # HEALPix side
CTRL = "#666666"      # controls
FILL = "#F4F4F4"


def box(ax, xy, w, h, text, color=INK, lw=0.8, fs=6.5, fill=FILL, style="round,pad=0.02,rounding_size=0.08"):
    x, y = xy
    p = FancyBboxPatch((x, y), w, h, boxstyle=style, linewidth=lw,
                       edgecolor=color, facecolor=fill)
    ax.add_patch(p)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color=INK, linespacing=1.15)
    return (x, y, w, h)


def arrow(ax, a, b, color=INK, lw=0.7, style="-|>", text=None, fs=5.5, dx=0, dy=0, conn="arc3,rad=0"):
    p = FancyArrowPatch(a, b, arrowstyle=style, mutation_scale=7, linewidth=lw,
                        color=color, connectionstyle=conn, shrinkA=1, shrinkB=1)
    ax.add_patch(p)
    if text:
        ax.text((a[0] + b[0]) / 2 + dx, (a[1] + b[1]) / 2 + dy, text, fontsize=fs,
                ha="center", va="center", color=color,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none"))


def draw():
    fig, ax = plt.subplots(figsize=(COLUMN_IN, 3.9))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 11)
    ax.axis("off")

    # --- source ------------------------------------------------------------
    src = box(ax, (2.6, 9.7), 4.8, 1.0,
              "Sentinel-2 L2A B04, 10 m\n40 regions, 1512$^2$-px native window",
              color=INK, fs=6.5)

    # --- two branches ------------------------------------------------------
    utm = box(ax, (0.3, 7.4), 4.1, 1.35,
              "Native UTM square\n637$^2$ px, equal nominal area\n(never resampled)",
              color=GRID, fs=6.3)
    hpx = box(ax, (5.6, 7.4), 4.1, 1.35,
              "HEALPix level-20 tile\n1024$^2$ cells, 2.58$\\times$ the samples\nMorton de-interleave to a raster",
              color=CARRIER, fs=6.3)
    arrow(ax, (4.0, 9.7), (2.35, 8.75), color=GRID, text="crop", dx=-0.55, dy=0.05)
    arrow(ax, (6.0, 9.7), (7.65, 8.75), color=CARRIER,
          text="transfer\nC\u2013T / PSF / IDW", dx=0.95, dy=0.1)

    # --- codec sweep ----------------------------------------------------------
    cod = box(ax, (1.5, 5.45), 7.0, 1.2,
              "Five codecs, full quality ladders\n"
              "ZFP, SZ3, JPEG 2000, JPEG XL, DWT+zstd\n"
              "size and RMSE at every setting",
              color=INK, fs=6.3)
    arrow(ax, (2.35, 7.4), (3.2, 6.65), color=GRID)
    arrow(ax, (7.65, 7.4), (6.8, 6.65), color=CARRIER)

    # --- references -----------------------------------------------------------
    ref1 = box(ax, (0.3, 3.55), 4.4, 1.25,
               "Within-representation reference\neach raster vs its own\nuncompressed form",
               color=INK, fs=6.0)
    ref2 = box(ax, (5.3, 3.55), 4.4, 1.25,
               "Common reference\nread back on native pixels\nvs the observation",
               color=INK, fs=6.0)
    arrow(ax, (3.6, 5.45), (2.5, 4.8), color=INK)
    arrow(ax, (6.4, 5.45), (7.5, 4.8), color=INK)

    # --- outcome --------------------------------------------------------------
    out = box(ax, (1.0, 1.55), 8.0, 1.35,
              "Sizes at equal RMSE, oracle-best codec, log\u2013log interpolation\n"
              "HPX/UTM ratio, critical overhead $\\gamma_\\star(e)$, recovered fraction $f(e)$\n"
              "paired tests over regions, Holm-corrected",
              color=INK, fs=6.0, fill="#E8EEF7")
    arrow(ax, (2.5, 3.55), (3.6, 2.9), color=INK)
    arrow(ax, (7.5, 3.55), (6.4, 2.9), color=INK)

    # --- controls, side column ------------------------------------------------
    ctl = box(ax, (0.3, 0.05), 9.4, 1.1,
              "Controls: analytic field on both lattices (direct, via the transfer),\n"
              "smoothed native, level-19 carrier compressed, same codec on both grids,\n"
              "tiling overhead from the product geometry",
              color=CTRL, fs=5.6, fill="white", style="round,pad=0.02,rounding_size=0.05")
    arrow(ax, (5.0, 1.55), (5.0, 1.1), color=CTRL, style="<|-|>")

    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    return fig


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-publish", action="store_true")
    args = ap.parse_args(argv)
    fig = draw()
    cfg.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(cfg.FIGURES_DIR / f"fig_workflow.{ext}", dpi=200)
    if not args.no_publish:
        cfg.TEX_FIGURES.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(cfg.FIGURES_DIR / "fig_workflow.pdf",
                        cfg.TEX_FIGURES / "fig_workflow.pdf")
    print(f"  fig_workflow -> {cfg.FIGURES_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
