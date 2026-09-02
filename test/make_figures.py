"""Regenerate the paper figures for the 40-region protocol.

    fig_ratio_boxplots        per-region HPX/UTM ratio distributions (Fig. 3)
    fig_utm_crop_footprints   HEALPix tile vs equal-area UTM square (Fig. 1)
    fig_region_map            where the 40 regions sit, and in which HEALPix regime
    fig_latitude_effect       intrinsic ratio against |latitude|

Figures are written to test/figures/ and, unless --no-publish, copied into
tex/figures/ as PDF.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402

import paper_config as cfg                # noqa: E402
import paper_stats as st                  # noqa: E402

CLASS_COLOR = {"urban": "#D55E00", "water": "#0072B2",
               "forest": "#009E73", "agriculture": "#E69F00"}
RMSE_TICK = {1e-4: "1e-4", 3e-4: "3e-4", 1e-3: "1e-3", 3e-3: "3e-3", 1e-2: "1e-2"}


def _save(fig, name, publish=True):
    cfg.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(cfg.FIGURES_DIR / f"{name}.{ext}", bbox_inches="tight", dpi=160)
    if publish:
        cfg.TEX_FIGURES.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(cfg.FIGURES_DIR / f"{name}.pdf",
                        cfg.TEX_FIGURES / f"{name}.pdf")
    plt.close(fig)
    print(f"  {name}")


def _load(metric="self"):
    hpx = pd.read_csv(cfg.RESULTS_DIR / "hpx_codec_bench.csv")
    utm = pd.read_csv(cfg.RESULTS_DIR / "utm_codec_bench.csv")
    rmse_col = "rmse_common" if metric == "common" else "rmse"
    scenes = sorted(set(hpx.scene) & set(utm.scene))
    idx_path = cfg.TILE_DIR / "tiles_index.csv"
    index = pd.read_csv(idx_path) if idx_path.exists() else pd.DataFrame()
    codecs = [c for c in cfg.COMPRESSOR_CODECS if c in set(hpx.codec)]
    return hpx, utm, scenes, codecs, rmse_col, index


def ratios_per_target(hpx, utm, scenes, codecs, rmse_col, gamma=cfg.GAMMA_UTM):
    out = {}
    for t in cfg.TARGET_RMSE:
        d = {}
        for s in scenes:
            h = st.best_of_codecs_size(hpx, s, codecs, t, rmse_col=rmse_col)
            u = st.best_of_codecs_size(utm, s, codecs, t, gamma=gamma)
            if np.isfinite(h) and np.isfinite(u) and u > 0:
                d[s] = h / u
        out[("rmse", t)] = d
    for f in cfg.BUDGETS:
        d = {}
        for s in scenes:
            h = st.best_of_codecs_rmse(hpx, s, codecs, f * cfg.SHARED_RAW, rmse_col=rmse_col)
            u = st.best_of_codecs_rmse(utm, s, codecs, f * cfg.SHARED_RAW / gamma)
            if np.isfinite(h) and np.isfinite(u) and u > 0:
                d[s] = h / u
        out[("budget", f)] = d
    return out


def fig_ratio_boxplots(ratios, index, publish=True):
    cls = index.set_index("scene")["scene_class"].to_dict() if "scene_class" in index else {}
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, kind, targets, xlabel in (
            (axL, "rmse", cfg.TARGET_RMSE, "Target RMSE"),
            (axR, "budget", cfg.BUDGETS, "Storage budget (of shared raw)")):
        data, labels, kept = [], [], []
        for t in targets:
            d = ratios.get((kind, t), {})
            if len(d) < cfg.MIN_SCENES:
                continue
            kept.append(t)
            data.append(list(d.values()))
            labels.append(RMSE_TICK.get(t, f"{100*t:g}%") if kind == "rmse"
                          else f"{100*t:g}%")
        if not data:
            continue
        try:
            ax.boxplot(data, tick_labels=labels, showfliers=False,
                       medianprops=dict(color="black"))
        except TypeError:                       # matplotlib < 3.9
            ax.boxplot(data, labels=labels, showfliers=False,
                       medianprops=dict(color="black"))
        for i, (t, vals) in enumerate(zip(kept, data), 1):
            scenes_i = list(ratios[(kind, t)])
            x = np.random.default_rng(i).normal(i, 0.05, len(vals))
            colors = [CLASS_COLOR.get(cls.get(s, ""), "0.4") for s in scenes_i]
            ax.scatter(x, vals, s=14, c=colors, alpha=0.8, zorder=3, linewidths=0)
        ax.axhline(1.0, color="grey", linestyle="--", linewidth=1)
        ax.set_yscale("log")
        ax.set_xlabel(xlabel)
    axL.set_ylabel("HPX / UTM compressed-size ratio")
    axR.set_ylabel("HPX / UTM RMSE ratio")
    axL.set_title("Storage ratio at constant RMSE")
    axR.set_title("Error ratio at constant budget")
    handles = [plt.Line2D([], [], marker="o", linestyle="", color=c, label=k)
               for k, c in CLASS_COLOR.items()]
    if cls:
        fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
                   bbox_to_anchor=(0.5, -0.06))
    _save(fig, "fig_ratio_boxplots", publish)


def fig_footprints(scenes, publish=True, max_panels=8):
    import build_tiles
    shown = scenes[:max_panels]
    n = len(shown)
    if n == 0:
        print("  fig_utm_crop_footprints: skipped (no tiles under data/tiles)")
        return
    ncol = min(4, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 3.6 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, scene in zip(axes, shown):
        z = np.load(build_tiles.tile_path(scene))
        ax.fill(z["hull_lon"], z["hull_lat"], color="#0072B2", alpha=0.8,
                edgecolor="black", linewidth=0.5)
        lon, lat = z["utm_lon"], z["utm_lat"]
        cx = np.concatenate([lon[0, :], lon[:, -1], lon[-1, ::-1], lon[::-1, 0]])
        cy = np.concatenate([lat[0, :], lat[:, -1], lat[-1, ::-1], lat[::-1, 0]])
        ax.plot(cx, cy, color="#D55E00", linewidth=1.8)
        ax.set_title(scene.replace("region__", "").replace("__", " / "), fontsize=9)
        ax.set_aspect("equal")
        ax.tick_params(labelsize=7)
    for ax in axes[n:]:
        ax.axis("off")
    handles = [plt.Line2D([], [], color="#0072B2", linewidth=6,
                          label="HEALPix tile footprint"),
               plt.Line2D([], [], color="#D55E00", linewidth=2,
                          label=f"equal-area UTM square ({cfg.UTM_SIDE_PX}$^2$ px)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    _save(fig, "fig_utm_crop_footprints", publish)


def fig_region_map(index, publish=True):
    if "lat" not in index:
        return
    fig, ax = plt.subplots(figsize=(9, 4.4))
    for cls, g in index.groupby("scene_class"):
        ax.scatter(g.lon, g.lat, s=42, label=cls,
                   color=CLASS_COLOR.get(cls, "0.4"),
                   edgecolor="black", linewidth=0.4, zorder=3)
    for lat in (cfg.EQUATORIAL_LAT_DEG, -cfg.EQUATORIAL_LAT_DEG):
        ax.axhline(lat, color="grey", linestyle=":", linewidth=1)
    ax.text(-178, cfg.EQUATORIAL_LAT_DEG + 2,
            r"HEALPix equatorial belt boundary  $|\mathrm{lat}|=41.81^\circ$",
            fontsize=8, color="grey")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-70, 75)
    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")
    ax.grid(alpha=0.25, linewidth=0.5)
    ax.legend(ncol=4, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.32))
    ax.set_title(f"{len(index)} evaluation regions, shared with healpix-resample")
    _save(fig, "fig_region_map", publish)


def fig_latitude(hpx, utm, scenes, codecs, rmse_col, index, publish=True):
    if "abs_lat" not in index:
        return
    lat = index.set_index("scene")["abs_lat"].to_dict()
    intrinsic = ratios_per_target(hpx, utm, scenes, codecs, rmse_col, gamma=1.0)
    targets = [t for t in cfg.TARGET_RMSE
               if len(intrinsic.get(("rmse", t), {})) >= cfg.MIN_SCENES]
    if not targets:
        return
    fig, axes = plt.subplots(1, len(targets), figsize=(3.2 * len(targets), 3.4),
                             sharey=True)
    axes = np.atleast_1d(axes)
    for ax, t in zip(axes, targets):
        d = intrinsic[("rmse", t)]
        xs = [lat.get(s, np.nan) for s in d]
        ys = list(d.values())
        ax.scatter(xs, ys, s=26, color="#0072B2", edgecolor="black", linewidth=0.3)
        ax.axvline(cfg.EQUATORIAL_LAT_DEG, color="grey", linestyle=":")
        ax.axhline(1.0, color="grey", linestyle="--", linewidth=1)
        sp = st.spearman(xs, ys)
        ax.set_title(f"RMSE {t:g}\n" + rf"$\rho$={sp['rho']:.2f}, p={sp['p']:.3f}",
                     fontsize=9)
        ax.set_xlabel("|latitude| (deg)")
        ax.set_yscale("log")
    axes[0].set_ylabel(r"intrinsic HPX/UTM size ratio ($\gamma=1$)")
    _save(fig, "fig_latitude_effect", publish)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metric", choices=["self", "common"], default="self")
    ap.add_argument("--no-publish", action="store_true")
    ap.add_argument("--footprint-regions", nargs="*", default=None)
    args = ap.parse_args(argv)

    hpx, utm, scenes, codecs, rmse_col, index = _load(args.metric)
    publish = not args.no_publish
    print(f"{len(scenes)} regions; figures -> {cfg.FIGURES_DIR}"
          + (f" and {cfg.TEX_FIGURES}" if publish else ""))

    ratios = ratios_per_target(hpx, utm, scenes, codecs, rmse_col)
    fig_ratio_boxplots(ratios, index, publish)
    fp = args.footprint_regions or scenes
    fig_footprints([s for s in fp if (cfg.TILE_DIR / f"{s}.npz").exists()], publish)
    fig_region_map(index, publish)
    fig_latitude(hpx, utm, scenes, codecs, rmse_col, index, publish)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
