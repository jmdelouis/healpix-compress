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

    smoothed    Read the native crop through the same Clough-Tocher operator
                at points shifted by half a pixel, at the native density, and
                compress that. The native side then carries the carrier's
                smoothing and nothing else: smoothed/native is what the
                transfer alone buys, HEALPix/smoothed what remains of the
                carrier's cost once both sides are equally smooth (revision 5).

Usage
-----
    # the paper's controls (revision 5): direct and transferred synthetic
    # rows, plus the smoothed-native rows, in one run
    python -B support_controls.py --control synthetic smoothed --fitted --transfer clough_tocher
    python support_controls.py --control common --limit 10
    python support_controls.py --control shift --limit 5      # the expensive one
    python support_controls.py --control all --limit 5

The paper's synthetic rows (Table XI) use --fitted: one field with the slope
and response width read off the imagery (the medians of
results/support_controls_spectrum.csv, beta = 1.9 and 11.9 m FWHM). Without
it the defaults below are the exploratory grid beta = 3, 4, 5 at 12.5 m,
which is a different experiment; on 16 September that run silently replaced
the paper's rows, which is why an existing CSV is now archived before it is
overwritten, and why controls that are not re-run are reloaded from their
CSVs so that the table keeps every control it had.

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


def _save_sweeps(control: str, scene: str, frames):
    """Append the raw rate-distortion points behind every reported ratio.

    They are small, and without them a surprising ratio cannot be told from a
    sweep that simply failed to reach the target.
    """
    frames = [f for f in frames if f is not None and len(f)]
    if not frames:
        return
    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    path = cfg.RESULTS_DIR / f"support_controls_sweeps_{control}.csv"
    df = pd.concat(frames, ignore_index=True)
    # Appending across runs left the file holding several sweeps of the same
    # region under one label, which cannot be re-analysed: 8000 rows for 40
    # regions where one run writes 5680. One run, one file.
    first = control not in _save_sweeps.started
    df.to_csv(path, mode="w" if first else "a", header=first, index=False)
    _save_sweeps.started.add(control)


_save_sweeps.started = set()


def oracle_size(df, scene, codecs, target):
    return st.best_of_codecs_size(df, scene, codecs, target)


def ratio_at_targets(hpx_df, utm_df, scene, codecs, targets, gamma=1.0):
    """HPX/UTM oracle-best size ratio at each target RMSE, with bracketing.

    A ratio built from two censored bounds is not an estimate of anything, and
    on a synthetic field the sweeps can miss a target on one side only -- which
    would show up as a spurious jump in the ratio and be indistinguishable from
    a grid effect. Both flags travel with every ratio for that reason.
    """
    out = {}
    for t in targets:
        h = st.best_of_codecs_size(hpx_df, scene, codecs, t)
        u = st.best_of_codecs_size(utm_df, scene, codecs, t)
        out[t] = dict(
            ratio=(float(h / (u * gamma)) if (np.isfinite(h) and np.isfinite(u)
                                              and u > 0) else np.nan),
            hpx_bytes=float(h), utm_bytes=float(u),
            bracketed_hpx=st.oracle_bracketed(hpx_df, scene, codecs, t, "rmse"),
            bracketed_utm=st.oracle_bracketed(utm_df, scene, codecs, t, "rmse"))
    return out


def _both_bracketed(d):
    """Both grids under the paper's rule: the codec the oracle selects must
    bracket the target with two measured points (paper_stats.oracle_bracketed).
    Until revision 5 this control asked only whether *some* codec did, the
    weaker rule Section IV-B rejects; the synthetic rows of Table XI therefore
    have to be regenerated before their loose-target entries are quoted."""
    return bool(d["bracketed_hpx"] and d["bracketed_utm"])


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

    rows, sweeps = [], []
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
        sweeps.append(utm)
        for t, r in ratios.items():
            rows.append(dict(scene=scene, control="shift", dx=dx, dy=dy,
                             target=t, **r, both_bracketed=_both_bracketed(r),
                             crop_std=float(np.std(crop)),
                             crop_mean=float(np.mean(crop))))
        if verbose:
            shown = ", ".join(
                f"{t:g}:{ratios[t]['ratio']:.3f}"
                + ("" if _both_bracketed(ratios[t]) else "*")
                for t in targets if np.isfinite(ratios[t]["ratio"]))
            print(f"    shift ({dx:+d},{dy:+d}): {shown}")
    _save_sweeps("shift", scene, [hpx] + sweeps)
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

    _save_sweeps("common", scene, [hpx, utm])
    rows = []
    for t in targets:
        h = st.best_of_codecs_size(hpx, scene, codecs, t)
        u = st.best_of_codecs_size(utm, scene, codecs, t)
        rows.append(dict(
            scene=scene, control="common", target=t,
            bracketed_hpx=st.oracle_bracketed(hpx, scene, codecs, t, "rmse"),
            bracketed_utm=st.oracle_bracketed(utm, scene, codecs, t, "rmse"),
            both_bracketed=bool(st.oracle_bracketed(hpx, scene, codecs, t, "rmse")
                                and st.oracle_bracketed(utm, scene, codecs, t,
                                                     "rmse")),
            hpx_bytes=h, utm_bytes=u, hpx_km2=hpx_km2, utm_km2=utm_km2,
            hpx_bytes_per_km2=h / hpx_km2 if np.isfinite(h) else np.nan,
            utm_bytes_per_km2=u / utm_km2 if np.isfinite(u) else np.nan,
            ratio_per_km2=((h / hpx_km2) / (u / utm_km2)
                           if np.isfinite(h) and np.isfinite(u) and u > 0
                           else np.nan),
            hpx_side=side, utm_rows=r1 - r0, utm_cols=c1 - c0,
            mask_fraction=cr.mask_fraction))
    return pd.DataFrame(rows)


# --- the transfer, evaluated where the paper needs it -------------------------
def _to_lonlat(crs_wkt, x, y):
    """Projected (x, y) of the region's CRS back to geographic lon/lat."""
    import pyproj
    crs = pyproj.CRS.from_wkt(crs_wkt)
    tr = pyproj.Transformer.from_crs(crs, pyproj.CRS.from_epsg(4326),
                                     always_xy=True)
    lon, lat = tr.transform(np.asarray(x, dtype=np.float64),
                            np.asarray(y, dtype=np.float64))
    return np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64)


def ct_transfer_to_points(lon, lat, val, q_lon, q_lat, xy=None, q_xy=None,
                          verbose=False):
    """Clough--Tocher transfer from scattered samples to arbitrary points.

    The library's CloughTocherResampler only targets HEALPix cell centres. The
    smoothed-native control (revision 5, review point "free smoothing") needs
    the same operator evaluated on a UTM lattice shifted by a fraction of a
    pixel, so the operator is assembled here from the library's own building
    blocks -- gnomonic frame, Delaunay triangulation, 1-ring gradient
    operators, sparse M -- exactly as its constructor does, with the query
    points in place of the cell centres and nothing else changed. The package
    is pinned at tag jstars-submission, so these internals are stable for the
    paper.

    If they are unavailable, scipy's CloughTocher2DInterpolator on the
    projected plane is used instead and the implementation is recorded in the
    CSV: the two estimate gradients differently, and a row must say which
    operator produced it.

    Returns (values at the query points, NaN where outside the hull, impl).
    """
    lon = np.asarray(lon, dtype=np.float64).ravel()
    lat = np.asarray(lat, dtype=np.float64).ravel()
    val = np.asarray(val, dtype=np.float64).ravel()
    q_lon = np.asarray(q_lon, dtype=np.float64).ravel()
    q_lat = np.asarray(q_lat, dtype=np.float64).ravel()
    try:
        import torch
        from scipy.spatial import Delaunay
        from healpix_resample import clough_tocher as ct
        radius = 6371000.0                     # the library's default
        t, east, north = ct._gnomonic_tangent_frame(lon, lat)
        x, y, _den = ct._gnomonic_project(lon, lat, t, east, north, radius=radius)
        pts = np.stack([x, y], axis=-1)
        tri = Delaunay(pts)
        dtype, device = torch.float64, torch.device("cpu")
        pts_t = torch.as_tensor(pts, dtype=dtype, device=device)
        (_gx, _gy, edges_j, cgx, cgy, dgx, dgy, indptr) = \
            ct._build_gradient_operators(pts_t, tri, device, dtype,
                                         det_floor_rel=1e-9)
        qx, qy, qden = ct._gnomonic_project(q_lon, q_lat, t, east, north,
                                            radius=radius)
        simplex = tri.find_simplex(np.stack([qx, qy], axis=-1))
        simplex = np.where(qden <= 0.05, -1, simplex)
        keep = simplex >= 0
        M = ct._assemble_M(points2d=pts_t, tri=tri, simplex_idx=simplex[keep],
                           qx=qx[keep], qy=qy[keep], N=int(lon.size),
                           K=int(keep.sum()), edges_j=edges_j,
                           coef_gx_edge=cgx, coef_gy_edge=cgy, diag_gx=dgx,
                           diag_gy=dgy, nbr_indptr=indptr, device=device,
                           dtype=dtype)
        yv = torch.as_tensor(val, dtype=dtype, device=device)[None, :]
        out = np.full(q_lon.size, np.nan, dtype=np.float64)
        out[keep] = (yv @ M)[0].cpu().numpy()
        return out, "healpix_resample"
    except Exception as exc:                                   # noqa: BLE001
        if verbose:
            print(f"    library operator unavailable ({type(exc).__name__}: "
                  f"{exc}); falling back to scipy's Clough-Tocher")
    from scipy.interpolate import CloughTocher2DInterpolator
    if xy is None or q_xy is None:
        raise RuntimeError("scipy fallback needs projected sample and query "
                           "coordinates (xy, q_xy)")
    p = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
    q = np.asarray(q_xy, dtype=np.float64).reshape(-1, 2)
    f = CloughTocher2DInterpolator(p, val, fill_value=np.nan)
    return np.asarray(f(q), dtype=np.float64), "scipy"


# --- control E: the transfer alone, on the native lattice ---------------------
def control_smoothed(scene, codecs, targets, shift=(0.5, 0.5), verbose=True):
    """Native crop -> Clough--Tocher -> UTM lattice shifted by a fraction of a
    pixel, at the native density, compressed like the crop.

    The carrier is scored against its own smooth tile, and the review objects
    that the transfer therefore hands the carrier a free smoothing the native
    side never receives. This control gives the native side the same
    smoothing and nothing else: the crop is read through the same operator at
    points a fraction of a pixel away, so the result is a raster of the same
    density as the crop, filled by the same interpolant as the tile. Two
    ratios follow. Smoothed/native says what the transfer alone buys at each
    target; HEALPix/smoothed is what remains of the carrier's cost once both
    sides are equally smooth, i.e. the density and the lattice.

    Half a pixel is the phase at which the interpolant departs most from the
    samples; cell centres of the carrier sit at every phase, so this is the
    strongest smoothing the transfer can give, and the ratio is a bound in
    the carrier's favour. `shift` selects another phase.
    """
    reg = load_region(scene)
    crop = np.asarray(reg["utm_crop"], dtype=np.float64)
    n = crop.shape[0]
    x, y = reg["utm_x"], reg["utm_y"]
    fx, fy = float(shift[0]), float(shift[1])
    # Interior lattice, so that every query point lies inside the sample hull:
    # Clough--Tocher does not extrapolate.
    qx1 = x[:-1] + fx * (x[1:] - x[:-1])
    qy1 = y[:-1] + fy * (y[1:] - y[:-1])
    qx, qy = np.meshgrid(qx1, qy1)
    gx, gy = np.meshgrid(x, y)
    q_lon, q_lat = _to_lonlat(reg["crs_wkt"], qx, qy)
    sm, impl = ct_transfer_to_points(
        reg["utm_lon"], reg["utm_lat"], crop, q_lon, q_lat,
        xy=np.stack([gx.ravel(), gy.ravel()], -1),
        q_xy=np.stack([qx.ravel(), qy.ravel()], -1), verbose=verbose)
    sm = sm.reshape(qx.shape)
    n_missing = int((~np.isfinite(sm)).sum())
    if n_missing:
        if n_missing / sm.size > 1e-4:
            raise RuntimeError(f"{scene}: {n_missing} shifted points fell outside "
                               f"the hull ({100*n_missing/sm.size:.3f}%)")
        filled, _ = build_tiles._fill_small_holes(sm.astype(np.float32))
        sm = np.asarray(filled, dtype=np.float64)
    if verbose:
        print(f"    smoothed native {sm.shape[0]}x{sm.shape[1]} px at phase "
              f"({fx:g},{fy:g}) via {impl}; std ratio "
              f"{np.std(sm)/np.std(crop):.4f}")

    utm = sweep_array(crop.astype(np.float32), scene, codecs, label="utm")
    smo = sweep_array(sm.astype(np.float32), scene, codecs, label="utm_smoothed")
    hpx = sweep_array(expanded(reg["cell_data"]), scene, codecs, label="hpx")
    _save_sweeps("smoothed", scene, [utm, smo, hpx])

    # Sizes per sample: the shifted lattice has one row and one column fewer.
    n_utm, n_smo = float(crop.size), float(sm.size)
    rows = []
    for t in targets:
        u = st.best_of_codecs_size(utm, scene, codecs, t)
        s = st.best_of_codecs_size(smo, scene, codecs, t)
        h = st.best_of_codecs_size(hpx, scene, codecs, t)
        bu = st.oracle_bracketed(utm, scene, codecs, t, "rmse")
        bs = st.oracle_bracketed(smo, scene, codecs, t, "rmse")
        bh = st.oracle_bracketed(hpx, scene, codecs, t, "rmse")
        ok_su = np.isfinite(u) and np.isfinite(s) and u > 0
        ok_hs = np.isfinite(h) and np.isfinite(s) and s > 0
        rows.append(dict(
            scene=scene, control="smoothed", shift_x=fx, shift_y=fy,
            transfer_impl=impl, target=t,
            utm_bytes=float(u), smoothed_bytes=float(s), hpx_bytes=float(h),
            bracketed_utm=bool(bu), bracketed_smoothed=bool(bs),
            bracketed_hpx=bool(bh),
            both_bracketed=bool(bu and bs),
            hpx_both_bracketed=bool(bh and bs),
            ratio_smoothed_native=((s / n_smo) / (u / n_utm)) if ok_su else np.nan,
            ratio_hpx_smoothed=(h / (s * n_utm / n_smo)) if ok_hs else np.nan,
            ratio_hpx_native=(h / u) if (np.isfinite(h) and np.isfinite(u)
                                         and u > 0) else np.nan,
            std_ratio=float(np.std(sm) / np.std(crop))))
    if verbose:
        shown = ", ".join(
            f"{r['target']:g}:{r['ratio_smoothed_native']:.3f}/"
            f"{r['ratio_hpx_smoothed']:.3f}"
            + ("" if r["both_bracketed"] and r["hpx_both_bracketed"] else "*")
            for r in rows if np.isfinite(r["ratio_smoothed_native"]))
        print(f"    smoothed/native, HPX/smoothed: {shown}")
    return pd.DataFrame(rows)


# --- control C: one analytic field on both geometries ------------------------
def psf_sigma_from_fwhm(fwhm_m: float) -> float:
    """Gaussian standard deviation of the instrument response, in metres.

    Defined here, not asked from the library. `healpix_resample.fwhm_to_scale`
    returns the operator *scale* s (FWHM/sqrt(2 ln 2)), not the standard
    deviation: the companion study writes its Gaussian as exp(-2 d^2 / s^2),
    whose standard deviation is s/2 = FWHM/(2 sqrt(2 ln 2)). The library's
    convention also flipped mid-project (it once returned s/2 and now returns
    s, exactly twice), which is what silently doubled this control's PSF on
    16 September: the field was band-limited at sigma = 10.1 m, an MTF of
    0.006 at the native Nyquist against the ESA B04 spec of 0.15-0.30, and
    the synthetic ratios fell to 1.8-1.1 instead of 2.2-1.2. So sigma is
    computed here and common_reference.py's reader uses the same formula; the
    library value is kept only as a recorded cross-check.
    """
    if not fwhm_m:
        return 0.0
    sigma = float(fwhm_m) / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    try:
        from healpix_resample.psf_geometry import fwhm_to_scale
        lib = float(fwhm_to_scale(float(fwhm_m)))
        # s/2 within 5%: the library returns a standard deviation, so trust it.
        # ~s (about twice): the new scale convention, ignore it.
        if abs(lib - sigma) / sigma < 0.05:
            return lib
    except Exception:                                      # noqa: BLE001
        pass
    return sigma


def synthetic_field(px, py, *, beta=2.0, k_min, k_max, n_modes=512, seed=0,
                    psf_sigma_m=0.0):
    """A band-limited random field evaluable at arbitrary positions.

    Built as a sum of plane waves with amplitudes following |k|^(-beta/2), so
    the power spectrum is a controlled power law and the field has no grid of
    its own: it is defined on the plane and merely *sampled* by each geometry.
    That is what makes the comparison content-free -- both grids see the same
    function, so any difference in compressed size is the sampling geometry and
    nothing else.

    `psf_sigma_m` convolves the field with the instrument's Gaussian response
    before sampling, and this matters more than anything else in this module.
    Sentinel-2's B04 response is about 12.5 m FWHM, so the real imagery carries
    very little power between the native 10 m sampling and the 6.22 m HEALPix
    cell: the carrier grid is oversampling a field that is already smooth at
    its own scale. A synthetic field with power all the way to the native
    Nyquist is a far harder object than anything the instrument delivers, and
    comparing the two grids on it measures the chosen band limit rather than
    the geometry. Because the field is a sum of plane waves the convolution is
    exact and costs nothing -- each mode's amplitude is scaled by
    exp(-k^2 sigma^2 / 2) -- and the field stays analytic and grid-free.
    """
    rng = np.random.default_rng(seed)
    k = np.exp(rng.uniform(math.log(k_min), math.log(k_max), n_modes))
    theta = rng.uniform(0, 2 * math.pi, n_modes)
    phase = rng.uniform(0, 2 * math.pi, n_modes)
    # `beta` is the slope of the 2-D power spectrum, so that it can be set
    # directly from what `--control spectrum` fits on the real imagery. With
    # modes drawn log-uniformly in k, the power in an annulus is amp^2 dk / k
    # and the 2-D PSD is that over 2 pi k dk, i.e. P(k) ~ amp^2 / k^2; taking
    # amp ~ k^((2-beta)/2) therefore gives P(k) ~ k^(-beta). The naive
    # amp ~ k^(-beta/2) instead yields P ~ k^(-beta-2), which would make every
    # fitted slope unusable as an input.
    amp = k ** (0.5 * (2.0 - beta))
    if psf_sigma_m:
        amp = amp * np.exp(-0.5 * (k * float(psf_sigma_m)) ** 2)
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


def _synthetic_transfer(reg, x0, y0, field_kwargs, mu, sd, target_mean,
                        target_std, kernel="clough_tocher", halo_px=32,
                        cell_xy=None):
    """The field sampled on the native lattice, then transferred to the cells.

    Reproduces the paper's own chain on the analytic field: native samples on
    a UTM lattice large enough to enclose the rotated tile (the crop alone
    leaves the tile's corners outside the sample hull), then the library's
    resampler to the level-20 cell centres, then the same hole rule as
    build_tiles. Returns the transferred cell values in the crop's units,
    the number of cells the resampler dropped, and the lattice shape.

    The lattice is sized on the projected cell centres themselves
    (`cell_xy`, absolute UTM coordinates) plus a halo: a square set from the
    tile's equivalent side is too small where the cells are anisotropic and
    the projected tile is a rhombus, and the resampler then drops the
    corners (0.9% of the cells on the first region tried, 17 September).
    """
    step = float(cfg.NATIVE_PIXEL_M)
    x_axis = np.asarray(reg["utm_x"], dtype=np.float64)
    y_axis = np.asarray(reg["utm_y"], dtype=np.float64)
    if cell_xy is None:
        xc, yc = float(np.mean(x_axis)), float(np.mean(y_axis))
        half = cfg.TILE_SIDE_M / math.sqrt(2.0)
        xmin, xmax, ymin, ymax = xc - half, xc + half, yc - half, yc + half
    else:
        cx, cy = cell_xy
        xmin, xmax = float(np.min(cx)), float(np.max(cx))
        ymin, ymax = float(np.min(cy)), float(np.max(cy))
    margin = halo_px * step
    # Lattice nodes on the crop's own pixel grid, so the native samples of
    # the control sit exactly where the product's pixels sit.
    sx = float(x_axis[1] - x_axis[0])
    sy = float(y_axis[1] - y_axis[0])
    i0 = int(math.floor((xmin - margin - x_axis[0]) / sx))
    i1 = int(math.ceil((xmax + margin - x_axis[0]) / sx))
    j0 = int(math.floor((ymin - margin - y_axis[0]) / sy))
    j1 = int(math.ceil((ymax + margin - y_axis[0]) / sy))
    if i1 < i0:
        i0, i1 = i1, i0
    if j1 < j0:
        j0, j1 = j1, j0
    xe = x_axis[0] + sx * np.arange(i0, i1 + 1, dtype=np.float64)
    ye = y_axis[0] + sy * np.arange(j0, j1 + 1, dtype=np.float64)
    gxe, gye = np.meshgrid(xe, ye)
    f_ue = synthetic_field(gxe - x0, gye - y0, **field_kwargs)
    f_ue = target_mean + target_std * (f_ue - mu) / sd
    lon_e, lat_e = _to_lonlat(reg["crs_wkt"], gxe, gye)
    rs = build_tiles.make_resampler(kernel, lon_e.ravel(), lat_e.ravel(),
                                    reg["cell_ids"])
    res = rs.resample(f_ue.ravel().astype(np.float64))
    full, n_missing = build_tiles.complete_block(res.cell_ids, res.cell_data,
                                                 reg["cell_ids"])
    return full, n_missing, (int(ye.size), int(xe.size))


def control_synthetic(scene, codecs, targets, betas=(3.0, 4.0, 5.0),
                      band="hpx", seed=0, psf_fwhm_m=12.5, verbose=True,
                      transfer=None):
    """Sample one analytic field on both lattices and compress both.

    With `transfer` set (e.g. "clough_tocher"), the field is also sampled on
    the native lattice and carried to the cells by that resampler, and the
    transferred carrier is swept as well. The direct rows isolate the sampling
    geometry; the transferred rows add the one step of the real chain the
    direct rows omit, and the difference between the two is the transfer's
    own contribution to the ratio, on a field whose truth is known.
    """
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
    # The instrument response, not an arbitrary cutoff, is what band-limits the
    # real imagery, so it is what band-limits this field. Modes are generated
    # up to the finer grid's Nyquist -- so the HEALPix side never aliases -- and
    # the PSF then decides how much survives at each scale, leaving the UTM
    # side to alias exactly the small residual it aliases in reality.
    nyq = {"utm": cfg.NATIVE_PIXEL_M, "hpx": cfg.DELTA_EQ}[band]
    k_max = math.pi / nyq
    sigma = psf_sigma_from_fwhm(psf_fwhm_m)
    k_nyq_utm = math.pi / cfg.NATIVE_PIXEL_M
    k_nyq_hpx = math.pi / cfg.DELTA_EQ
    mtf_utm = math.exp(-0.5 * (k_nyq_utm * sigma) ** 2) if sigma else 1.0
    mtf_hpx = math.exp(-0.5 * (k_nyq_hpx * sigma) ** 2) if sigma else 1.0
    if verbose:
        print(f"    PSF {psf_fwhm_m:g} m FWHM (sigma {sigma:.2f} m): "
              f"MTF {mtf_utm:.3f} at the UTM Nyquist, {mtf_hpx:.3f} at the "
              f"HEALPix Nyquist"
              if sigma else "    no PSF: the field is hard to the band limit")
    # Refuse to run with an out-of-band PSF rather than write silently wrong
    # ratios. A doubled sigma (the library scale used as a standard deviation,
    # see psf_sigma_from_fwhm) puts the MTF at 0.006 and drops the synthetic
    # ratios from 2.2-1.2 to 1.8-1.1; the Sentinel-2 B04 MTF at the native
    # Nyquist is 0.15-0.30. On 16 September a stale module produced exactly
    # that, and the numbers looked plausible enough to nearly reach the paper.
    if sigma and not (0.10 <= mtf_utm <= 0.40):
        raise ValueError(
            f"synthetic PSF out of band: FWHM {psf_fwhm_m:g} m gives sigma "
            f"{sigma:.2f} m and MTF {mtf_utm:.3f} at the native Nyquist, "
            f"outside the Sentinel-2 B04 window 0.15-0.30. The library's "
            f"fwhm_to_scale most likely returned a scale, not a standard "
            f"deviation, and it was used as sigma. psf_sigma_from_fwhm should "
            f"have caught this; check that the running module is the current "
            f"one (no stale .pyc, no shadowing copy on PYTHONPATH).")

    target_mean = float(np.mean(reg["cell_data"]))
    target_std = float(np.std(reg["cell_data"]))

    rows, sweeps = [], []
    for beta in betas:
        f_h = synthetic_field(hpx_x, hpx_y, beta=beta, k_min=k_min, k_max=k_max,
                              seed=seed, psf_sigma_m=sigma)
        f_u = synthetic_field(gx, gy, beta=beta, k_min=k_min, k_max=k_max,
                              seed=seed, psf_sigma_m=sigma)
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
        sweeps += [hpx, utm]
        common = dict(scene=scene, control="synthetic", beta=beta, band=band,
                      psf_fwhm_m=psf_fwhm_m, psf_sigma_m=sigma,
                      mtf_utm_nyquist=mtf_utm, mtf_hpx_nyquist=mtf_hpx,
                      hpx_std=float(np.std(f_h)), utm_std=float(np.std(f_u)))
        for t, r in ratios.items():
            rows.append(dict(**common, transfer="none", target=t, **r,
                             both_bracketed=_both_bracketed(r),
                             transfer_rmse_rel=np.nan, cells_filled=0))
        if verbose:
            shown = ", ".join(
                f"{t:g}:{ratios[t]['ratio']:.3f}"
                + ("" if _both_bracketed(ratios[t]) else "*")
                for t in targets if np.isfinite(ratios[t]["ratio"]))
            print(f"    beta={beta}: {shown}   (* = not bracketed on both sides)")

        if transfer:
            f_t, n_filled, shape = _synthetic_transfer(
                reg, x0, y0, dict(beta=beta, k_min=k_min, k_max=k_max,
                                  seed=seed, psf_sigma_m=sigma),
                mu, sd, target_mean, target_std, kernel=transfer,
                cell_xy=(cx, cy))
            # Transfer error against the field's exact value at the cells, in
            # units of the scene's standard deviation: the same quantity the
            # kernel ablation measures, here on a field with no noise.
            t_rmse = codec_lib.rmse(f_t, f_h) / (target_std or 1.0)
            hpx_t = sweep_array(expanded(f_t.astype(np.float32)), scene, codecs,
                                label=f"hpx_synth_{transfer}_b{beta}")
            ratios_t = ratio_at_targets(hpx_t, utm, scene, codecs, targets)
            sweeps.append(hpx_t)
            for t, r in ratios_t.items():
                rows.append(dict(**common, transfer=transfer, target=t, **r,
                                 both_bracketed=_both_bracketed(r),
                                 transfer_rmse_rel=float(t_rmse),
                                 cells_filled=int(n_filled)))
            if verbose:
                shown = ", ".join(
                    f"{t:g}:{ratios_t[t]['ratio']:.3f}"
                    + ("" if _both_bracketed(ratios_t[t]) else "*")
                    for t in targets if np.isfinite(ratios_t[t]["ratio"]))
                print(f"    beta={beta} via {transfer} (lattice {shape[0]}x{shape[1]}, "
                      f"transfer RMSE {100*t_rmse:.2f}% of std, "
                      f"{n_filled} cells filled): {shown}")
    _save_sweeps("synthetic", scene, sweeps)
    return pd.DataFrame(rows)


# --- control D: calibrate the synthetic field against the real imagery -------
def radial_psd(img, pixel_m):
    """Radially averaged power spectrum of a square image, in rad/m."""
    a = np.asarray(img, dtype=np.float64)
    a = a - a.mean()
    n = a.shape[0]
    # Hann window in both axes: without it the crop's edges inject a k^-2 tail
    # that would be read as scene power and bias every slope fitted below.
    w = np.hanning(n)
    a = a * np.outer(w, w)
    f = np.fft.fftshift(np.abs(np.fft.fft2(a)) ** 2)
    kx = np.fft.fftshift(np.fft.fftfreq(n, d=pixel_m)) * 2.0 * math.pi
    kk = np.hypot(*np.meshgrid(kx, kx))
    kmax = math.pi / pixel_m
    edges = np.linspace(0, kmax, 64)
    idx = np.digitize(kk.ravel(), edges) - 1
    ok = (idx >= 0) & (idx < len(edges) - 1)
    power = np.bincount(idx[ok], weights=f.ravel()[ok], minlength=len(edges) - 1)
    count = np.bincount(idx[ok], minlength=len(edges) - 1)
    centres = 0.5 * (edges[1:] + edges[:-1])
    with np.errstate(invalid="ignore", divide="ignore"):
        psd = power / np.maximum(count, 1)
    good = count > 0
    return centres[good], psd[good]


def fit_slope_and_psf(k, psd, k_lo, k_hi):
    """Fit log P(k) = c - beta log k - k^2 sigma^2 over [k_lo, k_hi].

    Two parameters, linear in (log k, k^2), so an ordinary least squares. The
    Gaussian term is what the instrument response contributes; recovering its
    width from the imagery is the point -- the nominal 12.5 m FWHM is a
    specification, and the control should use what the data actually show.
    """
    m = (k >= k_lo) & (k <= k_hi) & np.isfinite(psd) & (psd > 0)
    if m.sum() < 6:
        return np.nan, np.nan
    A = np.column_stack([np.ones(m.sum()), -np.log(k[m]), -(k[m] ** 2)])
    coef, *_ = np.linalg.lstsq(A, np.log(psd[m]), rcond=None)
    beta = float(coef[1])
    sig2 = float(coef[2])
    sigma = math.sqrt(sig2) if sig2 > 0 else np.nan
    return beta, sigma


def control_spectrum(scene, codecs, targets, verbose=True, **_ignored):
    """What slope and what instrument width does the real crop actually show?"""
    reg = load_region(scene)
    k, psd = radial_psd(reg["utm_crop"], cfg.NATIVE_PIXEL_M)
    k_nyq_utm = math.pi / cfg.NATIVE_PIXEL_M
    k_lo = 2.0 * math.pi / (cfg.UTM_SIDE_PX * cfg.NATIVE_PIXEL_M) * 4
    beta, sigma = fit_slope_and_psf(k, psd, k_lo, k_nyq_utm)
    fwhm = sigma * 2.0 * math.sqrt(2.0 * math.log(2.0)) if np.isfinite(sigma) else np.nan

    # How much of the crop's variance sits in the octave the carrier grid adds,
    # i.e. between the UTM Nyquist and the HEALPix Nyquist? On band-limited
    # imagery this is what the extra HEALPix samples can possibly be carrying.
    tot = float(np.sum(psd))
    hi = float(np.sum(psd[k >= 0.5 * k_nyq_utm]))
    if verbose:
        print(f"    slope beta={beta:.2f}, fitted PSF sigma={sigma:.2f} m "
              f"(FWHM {fwhm:.1f} m), variance above half-Nyquist "
              f"{100*hi/tot:.2f}%")
    return pd.DataFrame([dict(scene=scene, control="spectrum", beta_fit=beta,
                              psf_sigma_m=sigma, psf_fwhm_m=fwhm,
                              frac_var_above_half_nyquist=hi / tot,
                              k_lo=k_lo, k_hi=k_nyq_utm)])


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

    if "spectrum" in frames and len(frames["spectrum"]):
        d = frames["spectrum"].dropna(subset=["psf_fwhm_m"])
        if len(d):
            out += ["## What the real imagery's spectrum says", "",
                    "Radial power spectrum of each UTM crop, fitted as a power "
                    "law times a Gaussian instrument response. The fitted FWHM "
                    "is what the synthetic control should use: the nominal "
                    "12.5 m is a specification, this is what the data show.", "",
                    f"- fitted slope beta: median {d.beta_fit.median():.2f}, "
                    f"range [{d.beta_fit.min():.2f}, {d.beta_fit.max():.2f}]",
                    f"- fitted PSF FWHM: median {d.psf_fwhm_m.median():.1f} m, "
                    f"range [{d.psf_fwhm_m.min():.1f}, "
                    f"{d.psf_fwhm_m.max():.1f}] m",
                    f"- variance above half the UTM Nyquist: median "
                    f"{100*d.frac_var_above_half_nyquist.median():.2f}%", "",
                    "The same estimator is applied to the synthetic field, and "
                    "it carries a small known bias from the Hann window and the "
                    "radial binning (about +0.2 in slope and -1 m in FWHM on a "
                    "field of known parameters). Matching fitted-to-fitted "
                    "cancels it to first order, which is what feeding these "
                    "medians straight back in does:", "",
                    f"    python support_controls.py --control synthetic "
                    f"--betas {d.beta_fit.median():.1f} "
                    f"--psf-fwhm {d.psf_fwhm_m.median():.1f}", ""]

    if "synthetic" in frames and len(frames["synthetic"]):
        d = _with_transfer(frames["synthetic"])
        out += ["## One analytic field on both geometries", "",
                "Same function sampled by both lattices: no content difference "
                "at all, so the ratio is the grid effect alone. Rows where the "
                "target is bracketed on both sides for fewer than "
                f"{100*cfg.MIN_BRACKETED:.0f}% of regions are marked withheld: "
                "there the sweeps did not reach the target and the ratio is a "
                "quotient of bounds, not an estimate. `transfer` = none is the "
                "field evaluated at the cell centres; otherwise the field was "
                "sampled on the native lattice and carried to the cells by "
                "that resampler, the one step of the real chain the direct "
                "rows omit.", "",
                "| transfer | beta | target | n | bracketed | median ratio | "
                "range | reportable |", "|---|---|---|---|---|---|---|---|"]
        for (tr, beta, t), g in d.groupby(["transfer", "beta", "target"]):
            ok = g.dropna(subset=["ratio"])
            if not len(ok):
                continue
            br = ok[ok.both_bracketed]
            frac = len(br) / len(ok)
            use = br if len(br) else ok
            keep = frac >= cfg.MIN_BRACKETED
            out.append(f"| {tr} | {beta:g} | {t:g} | {len(ok)} | {frac:.0%} | "
                       f"{use.ratio.median():.3f} | "
                       f"[{use.ratio.min():.3f}, {use.ratio.max():.3f}] | "
                       f"{'yes' if keep else '**withheld**'} |")
        tr_rows = d[(d.transfer != "none") & d.transfer_rmse_rel.notna()]
        if len(tr_rows):
            per = tr_rows.groupby("scene").transfer_rmse_rel.first()
            out.append("")
            out.append(f"Transfer error on the analytic field (RMSE against the "
                       f"exact cell values, in units of the scene std): median "
                       f"{100*per.median():.2f}%, range [{100*per.min():.2f}%, "
                       f"{100*per.max():.2f}%] over {len(per)} regions.")
        out.append("")

    if "smoothed" in frames and len(frames["smoothed"]):
        d = frames["smoothed"]
        impl = ", ".join(sorted(d.transfer_impl.unique()))
        out += ["## The transfer alone, on the native lattice", "",
                "The native crop read through the same Clough-Tocher operator "
                "at points shifted by a fraction of a pixel, at the native "
                f"density (operator: {impl}). Smoothed/native is what the "
                "transfer alone buys at each target; HPX/smoothed is what "
                "remains of the carrier's cost once both sides carry the same "
                "smoothing. Sizes are per sample (the shifted lattice has one "
                "row and one column fewer).", "",
                f"- std ratio smoothed/native: median {d.std_ratio.median():.4f}",
                "", "| phase, target | n | smoothed/native | range | n | HPX/smoothed | "
                "range | HPX/native |", "|---|---|---|---|---|---|---|---|"]
        for (fx, fy, t), g in d.groupby([d.shift_x.round(3), d.shift_y.round(3),
                                         "target"]):
            a = g[g.both_bracketed].dropna(subset=["ratio_smoothed_native"])
            b = g[g.hpx_both_bracketed].dropna(subset=["ratio_hpx_smoothed"])
            c = g[g.bracketed_hpx & g.bracketed_utm].dropna(subset=["ratio_hpx_native"])
            if not len(a) and not len(b):
                continue
            out.append(
                f"| ({fx:g},{fy:g}) {t:g} | {len(a)} | "
                + (f"{a.ratio_smoothed_native.median():.3f} | "
                   f"[{a.ratio_smoothed_native.min():.3f}, "
                   f"{a.ratio_smoothed_native.max():.3f}]" if len(a) else "-- | --")
                + f" | {len(b)} | "
                + (f"{b.ratio_hpx_smoothed.median():.3f} | "
                   f"[{b.ratio_hpx_smoothed.min():.3f}, "
                   f"{b.ratio_hpx_smoothed.max():.3f}]" if len(b) else "-- | --")
                + " | " + (f"{c.ratio_hpx_native.median():.3f}" if len(c) else "--")
                + " |")
        out.append("")
    return "\n".join(out)


def _with_transfer(d: pd.DataFrame) -> pd.DataFrame:
    """Older synthetic CSVs predate the `transfer` column: they are direct rows."""
    d = d.copy()
    if "transfer" not in d.columns:
        d["transfer"] = "none"
    d["transfer"] = d["transfer"].fillna("none").astype(str)
    if "transfer_rmse_rel" not in d.columns:
        d["transfer_rmse_rel"] = np.nan
    return d


def latex_table(frames: dict) -> str:
    rows = []
    if "synthetic" in frames and len(frames["synthetic"]):
        d = _with_transfer(frames["synthetic"])
        names = {"none": "direct", "clough_tocher": "via C--T", "psf": "via PSF",
                 "linear": "via bilinear", "nearest": "via nearest",
                 "bicubic": "via bicubic"}
        for (tr, beta, t), g in d.groupby(["transfer", "beta", "target"]):
            ok = g.dropna(subset=["ratio"])
            br = ok[ok.both_bracketed] if len(ok) else ok
            # Same reporting rule as the paper's tables: a target that is not
            # bracketed for most regions is withheld rather than printed.
            if len(ok) and len(br) / len(ok) >= cfg.MIN_BRACKETED:
                rows.append(f"Synthetic, {names.get(tr, tr)}, $\\beta={beta:g}$ & "
                            f"${t:g}$ & {len(br)} & {br.ratio.median():.2f} & "
                            f"[{br.ratio.min():.2f}, {br.ratio.max():.2f}] \\\\")
    if "smoothed" in frames and len(frames["smoothed"]):
        d = frames["smoothed"]
        phases = sorted(set(zip(d.shift_x.round(3), d.shift_y.round(3))))
        for fx, fy in phases:
            dp = d[np.isclose(d.shift_x, fx) & np.isclose(d.shift_y, fy)]
            tag = "" if len(phases) == 1 else f", phase $({fx:g},{fy:g})$"
            for t, g in dp.groupby("target"):
                a = g[g.both_bracketed].dropna(subset=["ratio_smoothed_native"])
                if len(g) and len(a) / len(g) >= cfg.MIN_BRACKETED:
                    rows.append(f"Smoothed native / native{tag} & ${t:g}$ & {len(a)} & "
                                f"{a.ratio_smoothed_native.median():.2f} & "
                                f"[{a.ratio_smoothed_native.min():.2f}, "
                                f"{a.ratio_smoothed_native.max():.2f}] \\\\")
            for t, g in dp.groupby("target"):
                b = g[g.hpx_both_bracketed].dropna(subset=["ratio_hpx_smoothed"])
                if len(g) and len(b) / len(g) >= cfg.MIN_BRACKETED:
                    rows.append(f"HEALPix / smoothed native{tag} & ${t:g}$ & {len(b)} & "
                                f"{b.ratio_hpx_smoothed.median():.2f} & "
                                f"[{b.ratio_hpx_smoothed.min():.2f}, "
                                f"{b.ratio_hpx_smoothed.max():.2f}] \\\\")
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
            spread = float((per.max() - per.min()).median())
            med = float(g.ratio.median())
            eff = abs(med - 1.0)
            frac = spread / eff if eff > 1e-6 else float("nan")
            rows.append(
                f"Shifted crops, spread & ${t:g}$ & {g.scene.nunique()} & "
                f"{spread:.3f} & "
                + ("--" if frac != frac else f"{100*frac:.0f}\\% of effect")
                + " \\\\")
    if not rows:
        return ""
    # The caption describes the controls that actually ran: a caption listing
    # controls absent from the table is how a reader is told the wrong thing
    # about what was done.
    parts = ["The rotated HEALPix tile and the axis-aligned equal-area crop "
             "enclose the same nominal area but not the same ground; these "
             "controls bound what that mismatch contributes."]
    if "synthetic" in frames and len(frames["synthetic"]):
        d = _with_transfer(frames["synthetic"])
        fw = d["psf_fwhm_m"].median() if "psf_fwhm_m" in d else float("nan")
        parts.append(
            "The synthetic rows evaluate one analytic random field --- a power "
            "law of slope $\\beta$ times a Gaussian instrument response of "
            f"{fw:.1f}\\,m FWHM, both fitted to the imagery itself --- at both "
            "lattices' sample positions. In the direct rows neither side is "
            "resampled and both see the same function, so no content "
            "difference remains and the ratio is the sampling geometry alone.")
        if (d.transfer != "none").any():
            parts.append(
                "In the rows marked via a transfer, the field is sampled on the "
                "native lattice and carried to the cell centres by that "
                "resampler, as in the real chain; the difference from the "
                "direct rows is the transfer's own contribution.")
        parts.append(
            "Ratios are intrinsic ($\\gamma_{\\mathrm{tile}}=1$); $n$ counts "
            "the regions where the target is bracketed on both grids.")
    if "smoothed" in frames and len(frames["smoothed"]):
        d = frames["smoothed"]
        phases = sorted(set(zip(d.shift_x.round(3), d.shift_y.round(3))))
        ph = " and ".join(f"({fx:g}, {fy:g})" for fx, fy in phases)
        parts.append(
            "The smoothed-native rows read the native crop through the same "
            "Clough--Tocher operator at points shifted by "
            f"{ph} pixel, at the native density, so the native "
            "side receives the carrier's smoothing and nothing else; sizes "
            "are per sample. Smoothed native / native is what the transfer "
            "alone buys; HEALPix / smoothed native is the carrier's cost once "
            "both sides are equally smooth.")
    if "common" in frames and len(frames["common"]):
        parts.append("The common-support rows compress the largest "
                     "common-ground region on each native lattice and "
                     "normalise by the ground area each covers.")
    if "shift" in frames and len(frames["shift"]):
        parts.append("The shifted-crop rows report how far the ratio moves "
                     "when the UTM crop is re-cut at eight offsets of up to 16 "
                     "native pixels with the HEALPix side held fixed, which "
                     "bounds the mismatch's contribution without removing it; "
                     "their \\emph{Median} column is that spread, not a ratio.")
    return (
        "\\begin{table}[t]\n\\centering\n"
        + "\n".join(cfg.tex_caption(
            "Controls for the support mismatch", "tab:support_controls"))
        + "\n\\footnotesize\n\\begin{tabular}{@{}llrrl@{}}\n\\toprule\n"
        "Control & Target & $n$ & Median & Range \\\\\n\\midrule\n"
        + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n"
        + "\n".join(cfg.tex_note(" ".join(parts))) + "\n\\end{table}\n")


# --- driver ------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--control", nargs="*", default=["synthetic"],
                    choices=["spectrum", "shift", "common", "synthetic",
                             "smoothed", "all"])
    ap.add_argument("--transfer", choices=["none"] + list(build_tiles.KERNELS),
                    default="none",
                    help="for the synthetic control: also sample the field on "
                         "the native lattice and carry it to the cells with "
                         "this resampler, and sweep that carrier too (the "
                         "paper's rows use clough_tocher)")
    ap.add_argument("--smooth-shift", type=float, nargs=2, default=(0.5, 0.5),
                    metavar=("FX", "FY"),
                    help="for the smoothed control: phase of the shifted "
                         "native lattice, in pixels (default half a pixel in "
                         "both axes, the strongest smoothing the interpolant "
                         "gives)")
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--codecs", nargs="*", default=cfg.COMPRESSOR_CODECS)
    ap.add_argument("--targets", type=float, nargs="*", default=cfg.TARGET_RMSE)
    ap.add_argument("--betas", type=float, nargs="*", default=[3.0, 4.0, 5.0],
                    help="slopes of the 2-D power spectrum to synthesise. Set "
                         "these from what --control spectrum fits on the real "
                         "crops rather than leaving the defaults")
    ap.add_argument("--band", choices=["utm", "hpx"], default="hpx",
                    help="Nyquist up to which synthetic modes are generated "
                         "(default: the finer grid, so HEALPix never aliases; "
                         "the PSF, not this cutoff, is what band-limits the "
                         "field)")
    ap.add_argument("--psf-fwhm", type=float, default=12.5,
                    help="instrument response FWHM in metres applied to the "
                         "synthetic field (Sentinel-2 B04 is about 12.5 m). "
                         "Pass 0 for an unconvolved field, which is a much "
                         "harder object than the instrument ever delivers")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--fitted", action="store_true",
                    help="use the slope and response width fitted to the "
                         "imagery (medians of support_controls_spectrum.csv, "
                         "rounded to 0.1) instead of --betas/--psf-fwhm; this "
                         "is the paper's synthetic control")
    ap.add_argument("--no-tex", action="store_true")
    args = ap.parse_args(argv)

    controls = (["spectrum", "shift", "common", "synthetic", "smoothed"]
                if "all" in args.control else list(dict.fromkeys(args.control)))
    if args.fitted:
        spec = cfg.RESULTS_DIR / "support_controls_spectrum.csv"
        if not spec.exists():
            print(f"--fitted needs {spec.name}; run --control spectrum first")
            return 1
        sp = pd.read_csv(spec)
        args.betas = [round(float(sp.beta_fit.median()), 1)]
        args.psf_fwhm = round(float(sp.psf_fwhm_m.median()), 1)
        print(f"--fitted: beta = {args.betas[0]:g}, PSF {args.psf_fwhm:g} m FWHM "
              f"(medians over {len(sp)} regions)")
    # Any CSV about to be replaced is archived first, with a timestamp, so a
    # run with the wrong parameters can never again erase the paper's rows.
    import shutil, datetime
    stamp = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
    for control in controls:
        for name in (f"support_controls_{control}.csv",
                     f"support_controls_sweeps_{control}.csv"):
            src = cfg.RESULTS_DIR / name
            if src.exists():
                dst = cfg.RESULTS_DIR / "archive" / f"{src.stem}.{stamp}.csv"
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)
                print(f"  archived {name} -> archive/{dst.name}")

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
              "synthetic": control_synthetic, "smoothed": control_smoothed,
              "spectrum": control_spectrum}[control]
        parts = []
        for scene in scenes:
            print(f"  [{control}] {scene}")
            try:
                if control == "synthetic":
                    df = fn(scene, args.codecs, args.targets,
                            betas=args.betas, band=args.band, seed=args.seed,
                            psf_fwhm_m=args.psf_fwhm,
                            transfer=(None if args.transfer == "none"
                                      else args.transfer))
                elif control == "smoothed":
                    df = fn(scene, args.codecs, args.targets,
                            shift=tuple(args.smooth_shift))
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
            if control == "smoothed" and path.exists():
                # One CSV holds every phase swept: rows of other phases are
                # kept, rows of this phase are replaced.
                prev = pd.read_csv(path)
                if {"shift_x", "shift_y"} <= set(prev.columns):
                    fx, fy = tuple(args.smooth_shift)
                    keep = prev[~(np.isclose(prev.shift_x, fx)
                                  & np.isclose(prev.shift_y, fy))]
                    if len(keep):
                        frames[control] = pd.concat([keep, frames[control]],
                                                    ignore_index=True)
                        print(f"  kept {len(keep)} rows of other phases")
            frames[control].to_csv(path, index=False)
            print(f"  wrote {path} ({len(frames[control])} rows)")

    if not frames:
        print("no control produced any row")
        return 1
    # The table reports every control that has results, not only the ones
    # this invocation ran: a synthetic-only re-run must not drop the
    # common-support and shifted-crop rows from Table XI.
    for control in ("synthetic", "smoothed", "common", "shift"):
        path = cfg.RESULTS_DIR / f"support_controls_{control}.csv"
        if control not in frames and path.exists():
            frames[control] = pd.read_csv(path)
            print(f"  reloaded {path.name} ({len(frames[control])} rows) for the table")

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
