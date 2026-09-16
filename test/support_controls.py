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
    python support_controls.py --control synthetic --fitted   # the paper's control
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


# --- control C: one analytic field on both geometries ------------------------
def psf_sigma_from_fwhm(fwhm_m: float) -> float:
    """Gaussian sigma of the instrument response, in metres.

    Uses healpix-resample's own conversion when it is importable, so the
    control and the common-reference operator of common_reference.py describe
    the same instrument rather than two slightly different Gaussians.
    """
    if not fwhm_m:
        return 0.0
    try:
        from healpix_resample.psf_geometry import fwhm_to_scale
        return float(fwhm_to_scale(float(fwhm_m)))
    except Exception:                                      # noqa: BLE001
        return float(fwhm_m) / (2.0 * math.sqrt(2.0 * math.log(2.0)))


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


def control_synthetic(scene, codecs, targets, betas=(3.0, 4.0, 5.0),
                      band="hpx", seed=0, psf_fwhm_m=12.5, verbose=True):
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
        for t, r in ratios.items():
            rows.append(dict(scene=scene, control="synthetic", beta=beta,
                             band=band, psf_fwhm_m=psf_fwhm_m,
                             psf_sigma_m=sigma, mtf_utm_nyquist=mtf_utm,
                             mtf_hpx_nyquist=mtf_hpx, target=t, **r,
                             both_bracketed=_both_bracketed(r),
                             hpx_std=float(np.std(f_h)),
                             utm_std=float(np.std(f_u))))
        if verbose:
            shown = ", ".join(
                f"{t:g}:{ratios[t]['ratio']:.3f}"
                + ("" if _both_bracketed(ratios[t]) else "*")
                for t in targets if np.isfinite(ratios[t]["ratio"]))
            print(f"    beta={beta}: {shown}   (* = not bracketed on both sides)")
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
        d = frames["synthetic"]
        out += ["## One analytic field on both geometries", "",
                "Same function sampled by both lattices: no content difference "
                "at all, so the ratio is the grid effect alone. Rows where the "
                "target is bracketed on both sides for fewer than "
                f"{100*cfg.MIN_BRACKETED:.0f}% of regions are marked withheld: "
                "there the sweeps did not reach the target and the ratio is a "
                "quotient of bounds, not an estimate.", "",
                "| beta | target | n | bracketed | median ratio | range | "
                "reportable |", "|---|---|---|---|---|---|---|"]
        for (beta, t), g in d.groupby(["beta", "target"]):
            ok = g.dropna(subset=["ratio"])
            if not len(ok):
                continue
            br = ok[ok.both_bracketed]
            frac = len(br) / len(ok)
            use = br if len(br) else ok
            keep = frac >= cfg.MIN_BRACKETED
            out.append(f"| {beta:g} | {t:g} | {len(ok)} | {frac:.0%} | "
                       f"{use.ratio.median():.3f} | "
                       f"[{use.ratio.min():.3f}, {use.ratio.max():.3f}] | "
                       f"{'yes' if keep else '**withheld**'} |")
        out.append("")
    return "\n".join(out)


def latex_table(frames: dict) -> str:
    rows = []
    if "synthetic" in frames and len(frames["synthetic"]):
        d = frames["synthetic"]
        for (beta, t), g in d.groupby(["beta", "target"]):
            ok = g.dropna(subset=["ratio"])
            br = ok[ok.both_bracketed] if len(ok) else ok
            # Same reporting rule as the paper's tables: a target that is not
            # bracketed for most regions is withheld rather than printed.
            if len(ok) and len(br) / len(ok) >= cfg.MIN_BRACKETED:
                rows.append(f"Synthetic, $\\beta={beta:g}$ & ${t:g}$ & {len(br)} & "
                            f"{br.ratio.median():.2f} & "
                            f"[{br.ratio.min():.2f}, {br.ratio.max():.2f}] \\\\")
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
        d = frames["synthetic"]
        fw = d["psf_fwhm_m"].median() if "psf_fwhm_m" in d else float("nan")
        parts.append(
            "The synthetic rows evaluate one analytic random field --- a power "
            "law of slope $\\beta$ times a Gaussian instrument response of "
            f"{fw:.1f}\\,m FWHM, both fitted to the imagery itself --- at both "
            "lattices' sample positions. Neither side is resampled and both see "
            "the same function, so no content difference remains and the ratio "
            "is the sampling geometry alone. Ratios are intrinsic "
            "($\\gamma_{\\mathrm{tile}}=1$); $n$ counts the regions where the "
            "target is bracketed on both grids.")
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
                    choices=["spectrum", "shift", "common", "synthetic", "all"])
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

    controls = (["spectrum", "shift", "common", "synthetic"]
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
              "synthetic": control_synthetic,
              "spectrum": control_spectrum}[control]
        parts = []
        for scene in scenes:
            print(f"  [{control}] {scene}")
            try:
                if control == "synthetic":
                    df = fn(scene, args.codecs, args.targets,
                            betas=args.betas, band=args.band, seed=args.seed,
                            psf_fwhm_m=args.psf_fwhm)
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
    # The table reports every control that has results, not only the ones
    # this invocation ran: a synthetic-only re-run must not drop the
    # common-support and shifted-crop rows from Table XI.
    for control in ("synthetic", "common", "shift"):
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
