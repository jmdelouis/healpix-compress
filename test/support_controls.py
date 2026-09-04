"""Three controls for the support mismatch between the two footprints (V1).

The HEALPix carrier tile is a square rotated by about 45 degrees in UTM; the
reference crop is axis-aligned and of equal nominal area. They therefore
enclose the same amount of ground but not the same ground: the intersection is
81% of the crop in median, an intersection-over-union of about 0.68. Every
ratio in the paper is computed between two rasters that see partly different
scenes, and the review is right that this is not obviously negligible --
compressed size depends on texture, edges and heterogeneity, so two supports
that differ over a fifth of their extent can differ in entropy by an amount
comparable to the grid effect under study.

Nothing in the main pipeline can settle that, because it never compresses the
same ground twice. This module does, three ways, from the cheapest bound to the
cleanest isolation:

    shift       Re-crop the UTM square at eight offsets (+/- 8 and +/- 16
                native pixels in x and y) and measure how far the ratio moves.
                Each offset changes the crop's content by roughly as much as
                the footprint mismatch does, while the HEALPix side stays put.
                If the ratio barely moves, the mismatch is bounded and we can
                say by how much. This bounds the effect; it does not remove it.

    common      Compress, on each grid's own lattice, the largest region of
                that lattice whose ground lies inside the intersection -- the
                largest axis-aligned pixel rectangle inscribed in the mask on
                the UTM side, the largest Morton-aligned cell block on the
                HEALPix side -- and normalise by the ground area each actually
                covers, in bytes per square kilometre. Same ground, native
                lattices, smaller area.

    synthetic   Evaluate one analytic random field, band-limited and with a
                controlled power-law spectrum, at the sample positions of both
                grids. No resampling on either side and no content difference
                whatsoever, so whatever ratio survives is the grid effect
                alone. This is the control that isolates rather than bounds.

Usage
-----
    python support_controls.py --control synthetic            # cheapest, ~2 min/region
    python support_controls.py --control common --limit 10
    python support_controls.py --control shift --limit 5      # the expensive one
    python support_controls.py --control all --limit 5

Writes results/support_controls_<control>.csv, results/support_controls.md and
../tex/tables/table_support_controls.tex.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

import paper_config as cfg
import paper_stats as st
import codec_lib
import build_tiles
import region_sites


# Equal-area cell: the whole point of the pixelization, so the area a block of
# cells covers is exact and needs no geometry beyond a count.
CELL_AREA_M2 = 4.0 * math.pi * cfg.R_EARTH ** 2 / (12.0 * 4 ** cfg.LEVEL)
NATIVE_PIXEL_AREA_M2 = cfg.NATIVE_PIXEL_M ** 2


# --- sweeping ----------------------------------------------------------------
def sweep_array(arr, scene: str, codecs=None, label: str = "") -> pd.DataFrame:
    """Rate-distortion sweep of one raster, in the column shape paper_stats wants."""
    a = np.asarray(arr, dtype=np.float32)
    rows = []
    for run in codec_lib.iter_runs(a, codecs=codecs, on_error="warn"):
        rows.append(dict(scene=scene, label=label, codec=run.codec,
                         config=json.dumps(run.config, sort_keys=True),
                         total_bytes=int(run.total_bytes),
                         rmse=codec_lib.rmse(run.reconstruction, a),
                         n_samples=int(a.size)))
    return pd.DataFrame(rows)


def oracle_size(df, scene, codecs, target):
    return st.best_of_codecs_size(df, scene, codecs, target)


def ratio_at_targets(hpx_df, utm_df, scene, codecs, targets, gamma=1.0):
    """HPX/UTM oracle-best size ratio at each target RMSE."""
    out = {}
    for t in targets:
        h = st.best_of_codecs_size(hpx_df, scene, codecs, t)
        u = st.best_of_codecs_size(utm_df, scene, codecs, t)
        out[t] = float(h / (u * gamma)) if (np.isfinite(h) and np.isfinite(u)
                                            and u > 0) else np.nan
    return out


# --- geometry helpers --------------------------------------------------------
def tile_cell_ids(parent: int) -> np.ndarray:
    return parent * cfg.N_CHILDREN + np.arange(cfg.N_CHILDREN, dtype=np.int64)


def cell_lonlat(cell_ids):
    import healpix_geo
    lon, lat = healpix_geo.nested.healpix_to_lonlat(
        np.asarray(cell_ids, dtype=np.int64), cfg.LEVEL, ellipsoid="WGS84")
    return np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64)


def expanded(cell_data):
    """Morton de-interleave a flat NESTED block into its 2-D raster."""
    from hcc.compression.zorder import expand_tile
    return np.asarray(expand_tile(np.asarray(cell_data, dtype=np.float32)),
                      dtype=np.float32)


def expand_bool(flat_bool):
    """Same de-interleaving, for a per-cell boolean."""
    from hcc.compression.zorder import expand_tile
    out = expand_tile(np.asarray(flat_bool, dtype=np.float32))
    return np.asarray(out) > 0.5


def largest_true_rectangle(mask):
    """Largest all-True axis-aligned rectangle in a boolean image.

    Standard maximal-rectangle-in-histogram sweep, O(H*W). Returns
    (r0, r1, c0, c1) as half-open slices, or None if the mask is empty.
    """
    m = np.asarray(mask, dtype=bool)
    H, W = m.shape
    heights = np.zeros(W + 1, dtype=np.int64)   # sentinel column of height 0
    best = (0, None)
    for r in range(H):
        heights[:W] = np.where(m[r], heights[:W] + 1, 0)
        stack = []
        for c in range(W + 1):
            start = c
            while stack and stack[-1][1] > heights[c]:
                s, h = stack.pop()
                area = h * (c - s)
                if area > best[0]:
                    best = (area, (r - h + 1, r + 1, s, c))
                start = s
            stack.append((start, heights[c]))
    return best[1]


def largest_aligned_block(inside, min_side=32):
    """Largest Morton-aligned square block of cells entirely inside `inside`.

    Alignment matters: a block of side 2^k aligned on a multiple of 2^k in the
    de-interleaved raster is exactly one HEALPix cell at level 20-k, i.e. a real
    chunk of the hierarchy rather than an arbitrary window. Returns
    (r0, c0, side) or None.
    """
    m = np.asarray(inside, dtype=bool)
    n = m.shape[0]
    side = n
    while side >= min_side:
        k = n // side
        blocks = m.reshape(k, side, k, side).all(axis=(1, 3))
        if blocks.any():
            r, c = np.argwhere(blocks)[0]
            return int(r) * side, int(c) * side, side
        side //= 2
    return None


def crop_pixel_index(lon, lat, crs_wkt, utm_x, utm_y):
    """Fractional (col, row) of lon/lat on the crop's own pixel grid."""
    import pyproj
    crs = pyproj.CRS.from_wkt(crs_wkt)
    to_crs = pyproj.Transformer.from_crs(pyproj.CRS.from_epsg(4326), crs,
                                         always_xy=True)
    cx, cy = to_crs.transform(lon, lat)
    return build_tiles._to_index(np.asarray(cx), np.asarray(cy),
                                 np.asarray(utm_x), np.asarray(utm_y))


def load_region(scene: str):
    """Everything the controls need for one region, from its tile npz."""
    cell_ids, cell_data, utm_crop, meta, npz = build_tiles.load_tile(scene)
    return dict(scene=scene, cell_ids=cell_ids, cell_data=cell_data,
                utm_crop=np.asarray(utm_crop, dtype=np.float32), meta=meta,
                utm_x=np.asarray(npz["utm_x"], dtype=np.float64),
                utm_y=np.asarray(npz["utm_y"], dtype=np.float64),
                utm_lon=np.asarray(npz["utm_lon"], dtype=np.float64),
                utm_lat=np.asarray(npz["utm_lat"], dtype=np.float64),
                parent=int(npz["parent_id"]),
                crs_wkt=str(meta.get("crs_wkt", "")))


# --- control A: shifted-crop ensemble ---------------------------------------
SHIFTS = [(0, 0), (8, 0), (-8, 0), (0, 8), (0, -8),
          (16, 0), (-16, 0), (0, 16), (0, -16)]


def control_shift(scene, codecs, targets, shifts=SHIFTS, verbose=True):
    """Re-crop the UTM square at several offsets; the HEALPix side never moves."""
    reg = load_region(scene)
    store = cfg.REGION_DIR / f"{scene}_data.zarr"
    if not store.exists():
        raise FileNotFoundError(
            f"{scene}: {store} is missing. The shift control re-crops from the "
            f"native window, so build_region_windows.py must have been run and "
            f"its output kept (the tile .npz alone only holds the 637 px crop).")
    _ds, img, _lon, _lat, x, y, _crs = build_tiles._open_window(store, cfg.BAND)
    H, W = img.shape

    half = cfg.UTM_SIDE_PX // 2
    ccol, crow = build_tiles._to_index(reg["utm_x"][half], reg["utm_y"][half], x, y)
    ccol, crow = int(round(ccol)), int(round(crow))

    hpx = sweep_array(expanded(reg["cell_data"]), scene, codecs, label="hpx")

    rows = []
    for dx, dy in shifts:
        x1, x2 = ccol + dx - half, ccol + dx + half + 1
        y1, y2 = crow + dy - half, crow + dy + half + 1
        if x1 < 0 or y1 < 0 or x2 > W or y2 > H:
            if verbose:
                print(f"    shift ({dx:+d},{dy:+d}): outside the window, skipped")
            continue
        crop = np.asarray(img[y1:y2, x1:x2], dtype=np.float32)
        if not np.all(np.isfinite(crop)):
            filled, _ = build_tiles._fill_small_holes(crop)
            crop = np.asarray(filled, dtype=np.float32)
            if not np.all(np.isfinite(crop)):
                if verbose:
                    print(f"    shift ({dx:+d},{dy:+d}): no-data, skipped")
                continue
        utm = sweep_array(crop, scene, codecs, label=f"utm_{dx:+d}_{dy:+d}")
        ratios = ratio_at_targets(hpx, utm, scene, codecs, targets)
        for t, r in ratios.items():
            rows.append(dict(scene=scene, control="shift", dx=dx, dy=dy,
                             target=t, ratio=r,
                             crop_std=float(np.std(crop)),
                             crop_mean=float(np.mean(crop))))
        if verbose:
            shown = ", ".join(f"{t:g}:{ratios[t]:.3f}" for t in targets
                              if np.isfinite(ratios[t]))
            print(f"    shift ({dx:+d},{dy:+d}): {shown}")
    return pd.DataFrame(rows)


# --- control B: common central support --------------------------------------
def control_common(scene, codecs, targets, verbose=True):
    """Compress the largest common-ground region on each native lattice."""
    import common_reference
    reg = load_region(scene)
    cr = common_reference.CommonReference(reg["cell_ids"], reg["utm_lon"],
                                          reg["utm_lat"], mode="nearest")

    rect = largest_true_rectangle(cr.mask)
    if rect is None:
        raise RuntimeError(f"{scene}: the intersection mask is empty")
    r0, r1, c0, c1 = rect
    utm_sub = reg["utm_crop"][r0:r1, c0:c1]
    utm_km2 = (r1 - r0) * (c1 - c0) * NATIVE_PIXEL_AREA_M2 / 1e6

    # which cells fall inside the UTM crop at all?
    lon, lat = cell_lonlat(reg["cell_ids"])
    col, row = crop_pixel_index(lon, lat, reg["crs_wkt"], reg["utm_x"], reg["utm_y"])
    n = cfg.UTM_SIDE_PX
    inside_flat = (col >= 0) & (col <= n - 1) & (row >= 0) & (row <= n - 1)
    block = largest_aligned_block(expand_bool(inside_flat))
    if block is None:
        raise RuntimeError(f"{scene}: no aligned cell block lies inside the crop")
    br, bc, side = block
    hpx_sub = expanded(reg["cell_data"])[br:br + side, bc:bc + side]
    hpx_km2 = side * side * CELL_AREA_M2 / 1e6

    if verbose:
        print(f"    UTM rectangle {r1-r0}x{c1-c0} px ({utm_km2:.2f} km2, "
              f"{100*(r1-r0)*(c1-c0)/reg['utm_crop'].size:.0f}% of the crop); "
              f"HEALPix block {side}x{side} cells ({hpx_km2:.2f} km2, level "
              f"{cfg.LEVEL - int(round(math.log2(side)))} chunk)")

    hpx = sweep_array(hpx_sub, scene, codecs, label="hpx_common")
    utm = sweep_array(utm_sub, scene, codecs, label="utm_common")

    rows = []
    for t in targets:
        h = st.best_of_codecs_size(hpx, scene, codecs, t)
        u = st.best_of_codecs_size(utm, scene, codecs, t)
        rows.append(dict(
            scene=scene, control="common", target=t,
            hpx_bytes=h, utm_bytes=u, hpx_km2=hpx_km2, utm_km2=utm_km2,
            hpx_bytes_per_km2=h / hpx_km2 if np.isfinite(h) else np.nan,
            utm_bytes_per_km2=u / utm_km2 if np.isfinite(u) else np.nan,
            ratio_per_km2=((h / hpx_km2) / (u / utm_km2)
                           if np.isfinite(h) and np.isfinite(u) and u > 0
                           else np.nan),
            hpx_side=side, utm_rows=r1 - r0, utm_cols=c1 - c0,
            mask_fraction=cr.mask_fraction))
    return pd.DataFrame(rows)


# --- control C: one analytic field on both geometries ------------------------
def synthetic_field(px, py, *, beta=2.0, k_min, k_max, n_modes=512, seed=0):
    """A band-limited random field evaluable at arbitrary positions.

    Built as a sum of plane waves with amplitudes following |k|^(-beta/2), so
    the power spectrum is a controlled power law and the field has no grid of
    its own: it is defined on the plane and merely *sampled* by each geometry.
    That is what makes the comparison content-free -- both grids see the same
    function, so any difference in compressed size is the sampling geometry and
    nothing else.
    """
    rng = np.random.default_rng(seed)
    k = np.exp(rng.uniform(math.log(k_min), math.log(k_max), n_modes))
    theta = rng.uniform(0, 2 * math.pi, n_modes)
    phase = rng.uniform(0, 2 * math.pi, n_modes)
    amp = k ** (-0.5 * beta)
    kx, ky = k * np.cos(theta), k * np.sin(theta)

    x = np.asarray(px, dtype=np.float64).ravel()
    y = np.asarray(py, dtype=np.float64).ravel()
    out = np.zeros(x.size, dtype=np.float64)
    # One mode at a time: a 1024^2 sample set times even a few dozen modes
    # would be a several-hundred-megabyte temporary, and the accumulation is
    # cheap enough that there is nothing to gain by batching it.
    buf = np.empty_like(out)
    for a, kxi, kyi, ph in zip(amp, kx, ky, phase):
        np.multiply(x, kxi, out=buf)
        buf += y * kyi
        buf += ph
        np.cos(buf, out=buf)
        out += a * buf
    return out.reshape(np.shape(px))


def control_synthetic(scene, codecs, targets, betas=(1.5, 2.0, 2.5),
                      band="utm", seed=0, verbose=True):
    """Sample one analytic field on both lattices and compress both."""
    import pyproj
    reg = load_region(scene)

    crs = pyproj.CRS.from_wkt(reg["crs_wkt"])
    to_crs = pyproj.Transformer.from_crs(pyproj.CRS.from_epsg(4326), crs,
                                         always_xy=True)
    lon, lat = cell_lonlat(reg["cell_ids"])
    cx, cy = to_crs.transform(lon, lat)
    x0 = float(np.mean(cx))
    y0 = float(np.mean(cy))
    # Both geometries are expressed in the same projected frame, so the field
    # is evaluated on one plane and the two samplings are directly comparable.
    hpx_x = np.asarray(cx, dtype=np.float64) - x0
    hpx_y = np.asarray(cy, dtype=np.float64) - y0
    gx, gy = np.meshgrid(reg["utm_x"] - x0, reg["utm_y"] - y0)

    k_min = 2.0 * math.pi / cfg.TILE_SIDE_M
    # Band-limit at the coarser grid's Nyquist by default: if the field carried
    # power above what 10 m sampling can represent, the UTM side would be
    # penalised for aliasing rather than for its geometry, and the control
    # would measure the band limit instead of the grids.
    nyq = {"utm": cfg.NATIVE_PIXEL_M, "hpx": cfg.DELTA_EQ}[band]
    k_max = math.pi / nyq

    target_mean = float(np.mean(reg["cell_data"]))
    target_std = float(np.std(reg["cell_data"]))

    rows = []
    for beta in betas:
        f_h = synthetic_field(hpx_x, hpx_y, beta=beta, k_min=k_min, k_max=k_max,
                              seed=seed)
        f_u = synthetic_field(gx, gy, beta=beta, k_min=k_min, k_max=k_max,
                              seed=seed)
        # One common normalisation for both, derived from the field itself, so
        # the two rasters stay the same function rather than two rescalings.
        mu = float(np.mean(f_h))
        sd = float(np.std(f_h)) or 1.0
        f_h = target_mean + target_std * (f_h - mu) / sd
        f_u = target_mean + target_std * (f_u - mu) / sd

        hpx = sweep_array(expanded(f_h.astype(np.float32)), scene, codecs,
                          label=f"hpx_synth_b{beta}")
        utm = sweep_array(f_u.astype(np.float32), scene, codecs,
                          label=f"utm_synth_b{beta}")
        ratios = ratio_at_targets(hpx, utm, scene, codecs, targets)
        for t, r in ratios.items():
            rows.append(dict(scene=scene, control="synthetic", beta=beta,
                             band=band, target=t, ratio=r,
                             hpx_std=float(np.std(f_h)),
                             utm_std=float(np.std(f_u))))
        if verbose:
            shown = ", ".join(f"{t:g}:{ratios[t]:.3f}" for t in targets
                              if np.isfinite(ratios[t]))
            print(f"    beta={beta}: {shown}")
    return pd.DataFrame(rows)


# --- reporting ---------------------------------------------------------------
def summarize(frames: dict) -> str:
    out = ["# Support-mismatch controls (V1)", ""]

    if "shift" in frames and len(frames["shift"]):
        d = frames["shift"]
        out += ["## Shifted-crop ensemble", "",
                "How far the ratio moves when the UTM crop is re-cut at "
                "offsets of up to 16 native pixels, with the HEALPix side "
                "fixed. Spread is across offsets, within a region.", "",
                "| target | regions | offsets | median ratio | median spread "
                "(max-min) | spread / (ratio-1) |", "|---|---|---|---|---|---|"]
        for t, g in d.groupby("target"):
            per = g.groupby("scene").ratio
            spread = (per.max() - per.min()).median()
            med = g.ratio.median()
            effect = abs(med - 1.0)
            rel = spread / effect if effect > 1e-6 else np.nan
            out.append(f"| {t:g} | {g.scene.nunique()} | "
                       f"{g.groupby('scene').size().median():.0f} | {med:.3f} | "
                       f"{spread:.3f} | "
                       f"{'--' if not np.isfinite(rel) else f'{rel:.2f}'} |")
        out.append("")

    if "common" in frames and len(frames["common"]):
        d = frames["common"]
        out += ["## Common central support", "",
                "Largest common-ground region on each native lattice, "
                "normalised by the ground area each covers.", "",
                f"- UTM rectangle: {d.utm_rows.median():.0f}x"
                f"{d.utm_cols.median():.0f} px median "
                f"({d.utm_km2.median():.2f} km2)",
                f"- HEALPix block: {d.hpx_side.median():.0f} cells per side "
                f"({d.hpx_km2.median():.2f} km2)", "",
                "| target | regions | HPX B/km2 | UTM B/km2 | ratio |",
                "|---|---|---|---|---|"]
        for t, g in d.groupby("target"):
            ok = g.dropna(subset=["ratio_per_km2"])
            if not len(ok):
                continue
            out.append(f"| {t:g} | {len(ok)} | "
                       f"{ok.hpx_bytes_per_km2.median():.3g} | "
                       f"{ok.utm_bytes_per_km2.median():.3g} | "
                       f"{ok.ratio_per_km2.median():.3f} |")
        out.append("")

    if "synthetic" in frames and len(frames["synthetic"]):
        d = frames["synthetic"]
        out += ["## One analytic field on both geometries", "",
                "Same function sampled by both lattices: no content difference "
                "at all, so the ratio is the grid effect alone.", "",
                "| beta | target | regions | median ratio | range |",
                "|---|---|---|---|---|"]
        for (beta, t), g in d.groupby(["beta", "target"]):
            ok = g.dropna(subset=["ratio"])
            if not len(ok):
                continue
            out.append(f"| {beta:g} | {t:g} | {len(ok)} | "
                       f"{ok.ratio.median():.3f} | "
                       f"[{ok.ratio.min():.3f}, {ok.ratio.max():.3f}] |")
        out.append("")
    return "\n".join(out)


def latex_table(frames: dict) -> str:
    rows = []
    if "synthetic" in frames and len(frames["synthetic"]):
        d = frames["synthetic"]
        for (beta, t), g in d.groupby(["beta", "target"]):
            ok = g.dropna(subset=["ratio"])
            if len(ok):
                rows.append(f"Synthetic, $\\beta={beta:g}$ & ${t:g}$ & {len(ok)} & "
                            f"{ok.ratio.median():.2f} & "
                            f"[{ok.ratio.min():.2f}, {ok.ratio.max():.2f}] \\\\")
    if "common" in frames and len(frames["common"]):
        d = frames["common"].dropna(subset=["ratio_per_km2"])
        for t, g in d.groupby("target"):
            rows.append(f"Common support, B/km$^2$ & ${t:g}$ & {len(g)} & "
                        f"{g.ratio_per_km2.median():.2f} & "
                        f"[{g.ratio_per_km2.min():.2f}, "
                        f"{g.ratio_per_km2.max():.2f}] \\\\")
    if "shift" in frames and len(frames["shift"]):
        d = frames["shift"]
        for t, g in d.groupby("target"):
            per = g.groupby("scene").ratio
            rows.append(f"Shifted crops, spread & ${t:g}$ & {g.scene.nunique()} & "
                        f"{(per.max()-per.min()).median():.2f} & "
                        f"median ratio {g.ratio.median():.2f} \\\\")
    if not rows:
        return ""
    return (
        "\\begin{table}[t]\n\\centering\n"
        "\\caption{Controls for the support mismatch between the rotated "
        "HEALPix tile and the axis-aligned equal-area crop. The synthetic rows "
        "sample one analytic band-limited field at both lattices' positions, so "
        "no content difference remains and the ratio is the grid effect alone; "
        "the common-support rows compress the largest common-ground region on "
        "each native lattice and normalise by the ground area covered; the "
        "shifted-crop rows report how far the ratio moves when the UTM crop is "
        "re-cut at offsets of up to 16 native pixels, which bounds the "
        "mismatch's contribution without removing it.}\n"
        "\\label{tab:support_controls}\n"
        "\\begin{tabular}{llrrl}\n\\toprule\n"
        "Control & Target RMSE & $n$ & Median & Range \\\\\n\\midrule\n"
        + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n")


# --- driver ------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--control", nargs="*", default=["synthetic"],
                    choices=["shift", "common", "synthetic", "all"])
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--codecs", nargs="*", default=cfg.COMPRESSOR_CODECS)
    ap.add_argument("--targets", type=float, nargs="*", default=cfg.TARGET_RMSE)
    ap.add_argument("--betas", type=float, nargs="*", default=[1.5, 2.0, 2.5])
    ap.add_argument("--band", choices=["utm", "hpx"], default="utm",
                    help="Nyquist to band-limit the synthetic field at "
                         "(default: the coarser grid, so UTM is not penalised "
                         "for aliasing rather than for its geometry)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--no-tex", action="store_true")
    args = ap.parse_args(argv)

    controls = (["shift", "common", "synthetic"] if "all" in args.control
                else list(dict.fromkeys(args.control)))

    scenes = args.regions or list(region_sites.load_sites()["scene"])
    scenes = [s for s in scenes if (cfg.TILE_DIR / f"{s}.npz").exists()]
    if not scenes:
        print(f"no tiles under {cfg.TILE_DIR}; run build_tiles.py first")
        return 1
    if args.limit:
        scenes = scenes[:args.limit]

    print(f"{len(scenes)} region(s), controls {controls}, "
          f"codecs {args.codecs}")

    frames, failures = {}, []
    for control in controls:
        fn = {"shift": control_shift, "common": control_common,
              "synthetic": control_synthetic}[control]
        parts = []
        for scene in scenes:
            print(f"  [{control}] {scene}")
            try:
                if control == "synthetic":
                    df = fn(scene, args.codecs, args.targets,
                            betas=args.betas, band=args.band, seed=args.seed)
                else:
                    df = fn(scene, args.codecs, args.targets)
            except Exception as exc:                       # noqa: BLE001
                print(f"    FAILED {type(exc).__name__}: {exc}")
                failures.append((control, scene, f"{type(exc).__name__}: {exc}"))
                continue
            parts.append(df)
        if parts:
            frames[control] = pd.concat(parts, ignore_index=True)
            path = cfg.RESULTS_DIR / f"support_controls_{control}.csv"
            cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
            frames[control].to_csv(path, index=False)
            print(f"  wrote {path} ({len(frames[control])} rows)")

    if not frames:
        print("no control produced any row")
        return 1

    md = summarize(frames)
    (cfg.RESULTS_DIR / "support_controls.md").write_text(md, encoding="utf-8")
    print("\n" + md)

    if not args.no_tex:
        tex = latex_table(frames)
        if tex:
            out = cfg.TEX_TABLES / "table_support_controls.tex"
            cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
            out.write_text(tex, encoding="utf-8")
            print(f"wrote {out}")

    if failures:
        print(f"\n{len(failures)} region/control pair(s) failed:")
        for c, s, why in failures:
            print(f"  {c} {s}: {why}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
