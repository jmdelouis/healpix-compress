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
    index = _region_index(hpx, scenes)
    codecs = [c for c in cfg.COMPRESSOR_CODECS if c in set(hpx.codec)]
    return hpx, utm, scenes, codecs, rmse_col, index


def _region_index(hpx, scenes):
    """One row per benchmarked region: class, lon, lat, |lat|, regime.

    Built from the frozen manifest (region_sites) and the benchmark itself,
    never from data/tiles/tiles_index.csv alone: a partial build_tiles run
    once left that file with a single row, and the region map, the box-plot
    class colours and the latitude figure all silently inherited it (a map
    titled "1 evaluation regions"). tiles_index.csv is only merged in for
    tile-specific columns it may carry.
    """
    cols = ["scene", "scene_class", "lat", "abs_lat", "hpx_regime"]
    base = hpx.drop_duplicates("scene")[[c for c in cols if c in hpx.columns]]
    try:
        import region_sites
        sites = region_sites.load_sites()
        idx = sites[sites.scene.isin(scenes)][
            ["scene", "scene_class", "lon", "lat", "abs_lat", "hpx_regime"]]
    except Exception as exc:                              # noqa: BLE001
        print(f"  region index: manifest unavailable ({type(exc).__name__}: "
              f"{exc}); falling back to the benchmark columns (no longitude)")
        idx = base[base.scene.isin(scenes)]
    idx_path = cfg.TILE_DIR / "tiles_index.csv"
    if idx_path.exists():
        tiles = pd.read_csv(idx_path)
        extra = [c for c in tiles.columns if c not in idx.columns and c != "scene"]
        if extra:
            idx = idx.merge(tiles[["scene"] + extra], on="scene", how="left")
    missing = sorted(set(scenes) - set(idx.scene))
    if missing:
        print(f"  region index: {len(missing)} benchmarked region(s) absent from "
              f"the manifest: {', '.join(missing[:5])}{'...' if len(missing) > 5 else ''}")
    if len(idx) != len(scenes):
        print(f"  region index: {len(idx)} rows for {len(scenes)} benchmarked regions")
    return idx.sort_values("scene").reset_index(drop=True)


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


def _local_km(lon, lat, lon0, lat0):
    """Offsets in km from (lon0, lat0), local flat-Earth conversion.

    Plotting footprints in degrees is what made the previous version of this
    figure unreadable: one degree of longitude is not one degree of latitude, so
    `set_aspect("equal")` squashed every panel by cos(lat), differently at each
    site, and the axis ticks mixed a longitude offset with an absolute
    coordinate. In a local metric frame a square looks square everywhere, all
    panels share one scale, and the axes read in metres of ground.

    Over a 9 km tile the flat-Earth conversion is accurate to far better than a
    10 m pixel, which is all a figure needs; the quantitative overlap below is
    computed from HEALPix indices, not from this projection.
    """
    m_per_deg_lat = 111132.95
    m_per_deg_lon = 111320.0 * np.cos(np.radians(lat0))
    # healpix_geo returns cell-centre longitudes in [0, 360) while the UTM
    # crop carries [-180, 180): west of Greenwich the two footprints of one
    # region were 360 degrees (40 000 km) apart and every panel shared that
    # scale, which is what reduced the footprints to dots. Wrap the offset.
    dlon = _wrap_deg(np.asarray(lon, dtype=np.float64) - lon0)
    x = dlon * m_per_deg_lon
    y = (np.asarray(lat, dtype=np.float64) - lat0) * m_per_deg_lat
    return x / 1000.0, y / 1000.0


def _wrap_deg(d):
    """Longitude difference wrapped to (-180, 180]."""
    return (np.asarray(d, dtype=np.float64) + 180.0) % 360.0 - 180.0


def _overlap_fractions(z):
    """Intersection of the two footprints, as a fraction of the crop, and IoU.

    A UTM pixel is inside the carrier tile exactly when the level-20 cell
    containing it belongs to the tile's child block -- the same definition the
    common-reference mask uses, so the figure and Section V-E cannot disagree.
    Both footprints have the same nominal area, so IoU = f / (2 - f).
    """
    try:
        import healpix_geo
    except Exception:                                     # noqa: BLE001
        return np.nan, np.nan
    parent = int(z["parent_id"])
    ipix = healpix_geo.nested.lonlat_to_healpix(
        np.asarray(z["utm_lon"], dtype=np.float64).ravel(),
        np.asarray(z["utm_lat"], dtype=np.float64).ravel(),
        cfg.LEVEL, ellipsoid="WGS84").astype(np.int64)
    f = float((ipix // cfg.N_CHILDREN == parent).mean())
    return f, f / (2.0 - f)


def _spread_over_classes(scenes, max_panels):
    """Pick regions round-robin over the land-cover classes, in scene order,
    so that eight panels show two of each class rather than eight agriculture
    regions (the first eight alphabetically)."""
    by_cls = {}
    for s in scenes:
        by_cls.setdefault(s.split("__")[1] if "__" in s else "", []).append(s)
    out, queues = [], [list(v) for _, v in sorted(by_cls.items())]
    while len(out) < max_panels and any(queues):
        for q in queues:
            if q and len(out) < max_panels:
                out.append(q.pop(0))
    return out


def fig_footprints(scenes, publish=True, max_panels=8, ncol=4):
    import build_tiles
    avail = [s for s in scenes if build_tiles.tile_path(s).exists()]
    shown = _spread_over_classes(avail, max_panels)
    n = len(shown)
    if n == 0:
        print("  fig_utm_crop_footprints: skipped (no tiles under data/tiles)")
        return
    ncol = min(ncol, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.5 * ncol, 2.8 * nrow),
                             squeeze=False)
    axes = axes.ravel()

    panels, lim = [], 0.0
    for scene in shown:
        z = np.load(build_tiles.tile_path(scene))
        lon0 = float(_wrap_deg(np.mean(_wrap_deg(z["hull_lon"]))))
        lat0 = float(np.mean(z["hull_lat"]))
        hx, hy = _local_km(z["hull_lon"], z["hull_lat"], lon0, lat0)
        # close the hull ring
        hx, hy = np.append(hx, hx[0]), np.append(hy, hy[0])
        lon, lat = z["utm_lon"], z["utm_lat"]
        bx = np.concatenate([lon[0, :], lon[:, -1], lon[-1, ::-1], lon[::-1, 0]])
        by = np.concatenate([lat[0, :], lat[:, -1], lat[-1, ::-1], lat[::-1, 0]])
        sx, sy = _local_km(bx, by, lon0, lat0)
        frac, iou = _overlap_fractions(z)
        panels.append((scene, hx, hy, sx, sy, frac, iou))
        lim = max(lim, np.abs(hx).max(), np.abs(hy).max(),
                  np.abs(sx).max(), np.abs(sy).max())
    lim = float(np.ceil(lim * 1.08))

    for ax, (scene, hx, hy, sx, sy, frac, iou) in zip(axes, panels):
        ax.fill(hx, hy, color="#0072B2", alpha=0.75, edgecolor="#023d63",
                linewidth=0.7, zorder=2)
        ax.plot(sx, sy, color="#D55E00", linewidth=1.6, zorder=3)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_title(scene.replace("region__", "").replace("__", " / "),
                     fontsize=7.5, pad=3)
        if np.isfinite(frac):
            ax.text(0.03, 0.03, f"∩ {100*frac:.0f}%  IoU {iou:.2f}",
                    transform=ax.transAxes, fontsize=6.5, va="bottom",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none",
                              alpha=0.8))
        ax.tick_params(labelsize=6, length=2)
        ax.grid(alpha=0.2, linewidth=0.4)
    for ax in axes[n:]:
        ax.axis("off")
    for i, ax in enumerate(axes[:n]):
        if i % ncol:
            ax.set_yticklabels([])
        else:
            ax.set_ylabel("km", fontsize=7)
        if i < n - ncol:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("km", fontsize=7)

    handles = [plt.Line2D([], [], color="#0072B2", linewidth=6,
                          label="HEALPix carrier tile (level-20 cell centres)"),
               plt.Line2D([], [], color="#D55E00", linewidth=1.6,
                          label=f"equal-area UTM square ({cfg.UTM_SIDE_PX}"
                                r"$^2$ px)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=7.5, bbox_to_anchor=(0.5, -0.015))
    fig.suptitle("Footprints in a local metric frame, common scale for every "
                 "panel", fontsize=8.5, y=1.0)
    _save(fig, "fig_utm_crop_footprints", publish)
    fr = [p[5] for p in panels if np.isfinite(p[5])]
    if fr:
        print(f"    intersection over the crop: {100*np.median(fr):.1f}% median "
              f"[{100*min(fr):.1f}, {100*max(fr):.1f}]")


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
    ap.add_argument("--footprint-panels", type=int, default=8)
    args = ap.parse_args(argv)

    hpx, utm, scenes, codecs, rmse_col, index = _load(args.metric)
    publish = not args.no_publish
    print(f"{len(scenes)} regions; figures -> {cfg.FIGURES_DIR}"
          + (f" and {cfg.TEX_FIGURES}" if publish else ""))

    ratios = ratios_per_target(hpx, utm, scenes, codecs, rmse_col)
    fig_ratio_boxplots(ratios, index, publish)
    fp = args.footprint_regions or scenes
    fig_footprints([s for s in fp if (cfg.TILE_DIR / f"{s}.npz").exists()],
                   publish, max_panels=args.footprint_panels)
    fig_region_map(index, publish)
    fig_latitude(hpx, utm, scenes, codecs, rmse_col, index, publish)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
