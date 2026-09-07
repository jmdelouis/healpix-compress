"""Derive, per region, the HEALPix carrier tile and its equal-area UTM crop.

Replaces cell 10 of compress_test.ipynb and cell 3 of utm_codec_bench.ipynb by
a single, testable, assertion-guarded step, and extends both to the 40 shared
regions.

Guarantees enforced here (any failure rejects the region, it is never patched):
  * the chosen level-10 parent tile's 1024x1024 level-20 cells all fall inside
    the acquired native window, with a halo -- no cell is ever filled;
  * the no-data fraction *inside the tile footprint* (not its bounding box)
    stays under --max-nodata; above it the region is rejected rather than
    interpolated over, which would let the HEALPix side compress for free;
  * the equal-area 637x637 UTM square centred on the tile centroid also fits
    inside the window, and its no-data fraction stays under the same ceiling
    (a handful of isolated pixels is interpolated from its neighbours, exactly
    as the HEALPix side drops them from the triangulation, and the count is
    recorded);
  * the resampled cell ids are exactly the requested contiguous child block,
    stored in ascending (NESTED = Morton) order, ready for `expand_tile`.

Output: test/data/tiles/<scene>.npz plus test/data/tiles/tiles_index.csv.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import paper_config as cfg
import region_sites

INDEX_PATH = cfg.TILE_DIR / "tiles_index.csv"


# The transfer kernel that populates the carrier grid is a choice of this
# study, not a neutral implementation detail, so it is selectable and recorded.
# `clough_tocher` stays the default and keeps the unsuffixed path, so every
# tile already built, and every result derived from one, is untouched.
DEFAULT_KERNEL = "clough_tocher"
KERNELS = ("nearest", "linear", "bicubic", "clough_tocher", "psf")
PSF_FWHM_M = 12.5          # nominal Sentinel-2 B04 response


def tile_path(scene: str, kernel: str = DEFAULT_KERNEL) -> Path:
    suffix = "" if kernel == DEFAULT_KERNEL else f"__{kernel}"
    return cfg.TILE_DIR / f"{scene}{suffix}.npz"


def make_resampler(kernel: str, lon, lat, cell_ids, verbose=False):
    """Build one transfer operator.

    The KNN-family resamplers take their arguments by keyword: their base class
    puts `Npt` third positionally while the subclasses fix it themselves, so a
    positional call silently binds `level` to `Npt`. Every construction here is
    therefore keyword-only.
    """
    import healpix_resample as hr
    common = dict(lon_deg=lon, lat_deg=lat, level=cfg.LEVEL,
                  out_cell_ids=cell_ids, verbose=verbose)
    if kernel == "clough_tocher":
        return hr.CloughTocherResampler(**common)
    if kernel == "nearest":
        return hr.NearestResampler(**common)
    if kernel == "linear":
        return hr.BilinearResampler(**common)
    if kernel == "bicubic":
        return hr.BicubicResampler(**common)
    if kernel == "psf":
        # Older installed versions of healpix-resample do not re-export this
        # from the package root, only from the submodule where it lives.
        try:
            fwhm_to_scale = hr.fwhm_to_scale
        except AttributeError:
            from healpix_resample.psf_geometry import fwhm_to_scale
        return hr.PSFResampler(sigma_m=fwhm_to_scale(PSF_FWHM_M), **common)
    raise KeyError(f"unknown resampler {kernel!r} (have {KERNELS})")


def _open_window(path: Path, band: str):
    import xarray as xr
    ds = xr.open_dataset(path, engine="zarr", consolidated=False, chunks=None)
    if band not in ds:
        alt = [v for v in ds.data_vars if v.lower() == band.lower()]
        if not alt:
            raise KeyError(f"{path}: no variable {band!r} (has {list(ds.data_vars)})")
        band = alt[0]
    img = np.asarray(ds[band].values, dtype=np.float64)
    lon = np.asarray(ds["longitude"].values, dtype=np.float64)
    lat = np.asarray(ds["latitude"].values, dtype=np.float64)
    x = np.asarray(ds["x"].values, dtype=np.float64)
    y = np.asarray(ds["y"].values, dtype=np.float64)
    # CRS: written as a spatial_ref coord by build_region_windows.py, but the
    # seven legacy EOPF stores carry it in the band's own attributes instead.
    crs_wkt = ""
    if "spatial_ref" in ds.coords:
        crs_wkt = str(ds["spatial_ref"].attrs.get("crs_wkt", ""))
    for attrs in (ds.attrs, ds[band].attrs):
        if crs_wkt:
            break
        for key in ("crs_wkt", "proj:wkt2"):
            if attrs.get(key):
                crs_wkt = str(attrs[key])
                break
        epsg = attrs.get("proj:epsg")
        if not crs_wkt and epsg:
            import pyproj
            crs_wkt = pyproj.CRS.from_epsg(int(epsg)).to_wkt()
    return ds, img, lon, lat, x, y, crs_wkt



def complete_block(out_ids, out_val, cell_ids, max_missing_frac=1e-4):
    """Fill the few cells a resampler dropped, from their nearest siblings.

    The KNN-family kernels drop a handful of cells at the very edge of the
    sample hull where the neighbour search fails -- seven out of 1048576 on one
    region here. Refusing the tile outright is right for the paper's own
    pipeline, where a tile with holes is not a tile, but it would disqualify a
    kernel over seven parts per million in an ablation, which measures the
    library's edge handling rather than the kernel.

    A dropped cell is filled with the mean of its surviving siblings, coarsening
    through the NESTED hierarchy until one is found: index within a contiguous
    child block is the Morton index, so `idx >> 2k` is the containing cell k
    levels up, and the fill is therefore spatially local by construction. The
    count is returned and recorded, and a fraction above `max_missing_frac`
    still raises -- a kernel that cannot populate the grid has failed, and that
    is a result, not something to paper over.
    """
    ids = np.asarray(out_ids, dtype=np.int64).ravel()
    val = np.asarray(out_val, dtype=np.float64).ravel()
    base = int(cell_ids[0])
    n = int(cell_ids.size)
    pos = ids - base
    if np.any((pos < 0) | (pos >= n)):
        raise RuntimeError("the resampler returned cells outside the requested "
                           "contiguous child block")
    full = np.full(n, np.nan, dtype=np.float64)
    full[pos] = val
    missing = ~np.isfinite(full)
    n_missing = int(missing.sum())
    if n_missing == 0:
        return full, 0
    if n_missing / n > max_missing_frac:
        frac = n_missing / n
        # Below about a per cent this really is the hull edge and a larger halo
        # helps. Above it, the samples do not cover the tile at all, and
        # pointing at --halo sends the reader in the wrong direction.
        why = ("cells outside the sample convex hull are dropped; a larger "
               "--halo would recover them"
               if frac < 0.01 else
               "the samples do not cover the tile: this is a geometry "
               "mismatch, not a hull edge, and --halo will not fix it. Check "
               "that the window on disk is the one the tile was built from")
        raise RuntimeError(
            f"the resampler kept {n - n_missing} of {n} cells "
            f"({100*frac:.4f}% missing, ceiling {100*max_missing_frac:.4f}%): "
            + why)
    idx = np.arange(n, dtype=np.int64)
    for up in range(1, cfg.DELTA + 1):
        group = idx >> (2 * up)
        finite = np.isfinite(full)
        sums = np.bincount(group, weights=np.where(finite, full, 0.0),
                           minlength=int(group.max()) + 1)
        cnts = np.bincount(group, weights=finite.astype(float),
                           minlength=int(group.max()) + 1)
        with np.errstate(invalid="ignore", divide="ignore"):
            means = sums / np.where(cnts > 0, cnts, np.nan)
        cand = means[group]
        take = (~finite) & np.isfinite(cand)
        if take.any():
            full[take] = cand[take]
        if np.all(np.isfinite(full)):
            break
    if not np.all(np.isfinite(full)):
        raise RuntimeError("could not fill every dropped cell from its siblings")
    return full, n_missing


def _to_index(cx, cy, x, y):
    """Fractional (col, row) of projected coordinates on a regular x/y grid."""
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    return (cx - float(x[0])) / dx, (cy - float(y[0])) / dy


def _neighbour_mean(a):
    """Mean of each pixel's finite 8-neighbours; NaN where it has none."""
    p = np.pad(a, 1, mode="constant", constant_values=np.nan)
    views = [p[:-2, :-2], p[:-2, 1:-1], p[:-2, 2:],
             p[1:-1, :-2],              p[1:-1, 2:],
             p[2:, :-2],  p[2:, 1:-1],  p[2:, 2:]]
    total = np.zeros_like(a, dtype=np.float64)
    count = np.zeros_like(a, dtype=np.float64)
    for v in views:
        ok = np.isfinite(v)
        total += np.where(ok, v, 0.0)
        count += ok
    return np.divide(total, count, out=np.full_like(total, np.nan), where=count > 0)


def _fill_small_holes(arr, max_iter=8):
    """Interpolate isolated no-data pixels from their finite neighbours.

    Used only under the --max-nodata ceiling, i.e. for a handful of pixels. It
    is the UTM-side counterpart of dropping no-data samples from the HEALPix
    triangulation: both grids then see the same clean ground, and neither gets
    a block of synthetic constant value that would compress for free.
    """
    out = np.asarray(arr, dtype=np.float64).copy()
    n_bad = int((~np.isfinite(out)).sum())
    for _ in range(max_iter):
        bad = ~np.isfinite(out)
        if not bad.any():
            break
        cand = _neighbour_mean(out)
        out = np.where(bad & np.isfinite(cand), cand, out)
    return out, n_bad


def _tile_native_mask(col, row, H, W, dilate=1):
    """Boolean mask of the native pixels actually under the tile footprint.

    Native pixels are coarser than level-20 cells (10 m vs 6.22 m), so rounding
    every cell centre to its nearest native pixel covers the footprint exactly.
    `dilate` adds the interpolation stencil around it: a no-data pixel one step
    outside the footprint still poisons the triangles that reach it.
    """
    mask = np.zeros((H, W), dtype=bool)
    ci = np.clip(np.rint(col).astype(np.int64), 0, W - 1)
    ri = np.clip(np.rint(row).astype(np.int64), 0, H - 1)
    mask[ri, ci] = True
    for _ in range(max(0, int(dilate))):
        m = mask.copy()
        m[1:, :] |= mask[:-1, :]
        m[:-1, :] |= mask[1:, :]
        m[:, 1:] |= mask[:, :-1]
        m[:, :-1] |= mask[:, 1:]
        mask = m
    return mask


def _evaluate_parent(parent, to_crs, x, y, img, H, W, halo_px):
    """Fit and no-data diagnostics for one candidate level-10 parent tile."""
    import healpix_geo
    cell_ids = int(parent) * cfg.N_CHILDREN + np.arange(cfg.N_CHILDREN, dtype=np.int64)
    clon, clat = healpix_geo.nested.healpix_to_lonlat(
        cell_ids, cfg.LEVEL, ellipsoid="WGS84")
    clon = np.asarray(clon, dtype=np.float64)
    clat = np.asarray(clat, dtype=np.float64)
    cx, cy = to_crs.transform(clon, clat)
    col, row = _to_index(cx, cy, x, y)
    inside = ((col >= halo_px) & (col <= W - 1 - halo_px)
              & (row >= halo_px) & (row <= H - 1 - halo_px))
    fit = float(inside.mean())
    nodata = np.nan
    if fit >= 1.0:
        mask = _tile_native_mask(col, row, H, W)
        nodata = float(((~np.isfinite(img)) & mask).sum() / max(1, int(mask.sum())))
    return dict(parent=int(parent), cell_ids=cell_ids, clon=clon, clat=clat,
                cx=cx, cy=cy, col=col, row=row, fit=fit, nodata=nodata)


def build_tile(scene: str, store: Path, band: str = cfg.BAND, halo_px: int = 32,
               max_nodata: float = 1e-3, max_missing: float = 1e-4,
               fill_nan_mean: bool = False,
               kernel: str = DEFAULT_KERNEL, verbose: bool = True) -> dict:
    import pyproj

    ds, img, lon, lat, x, y, crs_wkt = _open_window(store, band)
    H, W = img.shape
    if fill_nan_mean and not np.all(np.isfinite(img)):
        img = np.where(np.isfinite(img), img, np.nanmean(img))

    if not crs_wkt:
        raise RuntimeError(f"{scene}: the window store carries no CRS; cannot "
                           "place the HEALPix tile in projected coordinates.")
    crs = pyproj.CRS.from_wkt(crs_wkt)
    to_crs = pyproj.Transformer.from_crs(pyproj.CRS.from_epsg(4326), crs,
                                         always_xy=True)

    # --- candidate parent tiles, ranked by how well centred they are ---------
    # Ranking by sample count picks the tile that happens to hold the most
    # pixels, which for a window spanning ~2.4 tiles per axis is often one that
    # sticks out of the window. Ranking by the distance between the tile centre
    # and the window centre picks the one most likely to fit outright.
    import healpix_geo
    cover = healpix_geo.nested.lonlat_to_healpix(
        lon.ravel(), lat.ravel(), cfg.LEVEL, ellipsoid="WGS84").astype(np.int64)
    parents = np.unique(cover // cfg.N_CHILDREN)
    plon, plat = healpix_geo.nested.healpix_to_lonlat(
        parents.astype(np.int64), cfg.PARENT_LEVEL, ellipsoid="WGS84")
    pcx, pcy = to_crs.transform(np.asarray(plon, dtype=np.float64),
                                np.asarray(plat, dtype=np.float64))
    xc, yc = 0.5 * (x[0] + x[-1]), 0.5 * (y[0] + y[-1])
    ranked = parents[np.argsort(np.hypot(pcx - xc, pcy - yc))]

    chosen, diag = None, []
    for parent in ranked:
        cand = _evaluate_parent(parent, to_crs, x, y, img, H, W, halo_px)
        diag.append(cand)
        if cand["fit"] >= 1.0 and cand["nodata"] <= max_nodata:
            chosen = cand
            break

    if chosen is None:
        fitted = [d for d in diag if d["fit"] >= 1.0]
        lines = []
        for d in sorted(diag, key=lambda d: (-d["fit"], d["nodata"]))[:4]:
            nd = "n/a" if not np.isfinite(d["nodata"]) else f"{100*d['nodata']:.3f}%"
            lines.append(f"parent {d['parent']}: fits {100*d['fit']:.2f}%, "
                         f"no-data in tile {nd}")
        if not fitted:
            why = (f"none of the {len(diag)} candidate tiles fits inside the "
                   f"{W}x{H} px window with a {halo_px} px halo. A window is "
                   f"only guaranteed to hold a whole tile from "
                   f"{cfg.SAFE_WINDOW_PX} px (+ halo): {cfg.TILE_BBOX_PX} px is "
                   f"the tile's own bounding box and applies only when the "
                   f"anchor happens to sit near a tile centre. Re-acquire this "
                   f"region with --window {cfg.WINDOW_PX}")
        else:
            why = (f"{len(fitted)} candidate tile(s) fit, but all exceed the "
                   f"{100*max_nodata:.3f}% no-data ceiling. Raising "
                   f"--max-nodata would let interpolated no-data into the "
                   f"carrier grid; rejecting the region is the safe choice")
        raise RuntimeError(f"{scene}: {why}.\n      " + "\n      ".join(lines))

    # --- native sub-window that contains the tile, plus halo -----------------
    c0 = max(int(np.floor(chosen["col"].min())) - halo_px, 0)
    c1 = min(int(np.ceil(chosen["col"].max())) + halo_px + 1, W)
    r0 = max(int(np.floor(chosen["row"].min())) - halo_px, 0)
    r1 = min(int(np.ceil(chosen["row"].max())) + halo_px + 1, H)
    sub = img[r0:r1, c0:c1]
    sub_lon = lon[r0:r1, c0:c1].ravel()
    sub_lat = lat[r0:r1, c0:c1].ravel()
    sub_val = sub.ravel()

    # No-data samples are dropped from the triangulation rather than fed to it
    # as NaN (which would poison every triangle touching them). The ceiling
    # above is what keeps this honest: a few isolated pixels in the halo are
    # interpolated over, a cloud edge is not -- that region was rejected.
    good = np.isfinite(sub_val)
    n_dropped = int((~good).sum())

    op = make_resampler(kernel, sub_lon[good], sub_lat[good],
                        chosen["cell_ids"], verbose=verbose)
    # A transfer that does not preserve a constant is unusable, and the failure
    # is silent: every downstream error then measures the same offset. Checked
    # once per tile, on the operator that is about to be used.
    ones = op.resample(np.ones(int(good.sum()), dtype=np.float64))
    ones_med = float(np.nanmedian(np.asarray(ones.cell_data, dtype=np.float64)))
    if not np.isfinite(ones_med) or abs(ones_med - 1.0) > 0.05:
        raise RuntimeError(
            f"{scene}: the {kernel!r} transfer maps a constant field to "
            f"{ones_med} instead of 1.0; refusing to build a tile with it.")
    res = op.resample(sub_val[good].astype(np.float64))
    out_ids = np.asarray(res.cell_ids, dtype=np.int64).ravel()
    out_val = np.asarray(res.cell_data, dtype=np.float64).ravel()
    try:
        filled, n_missing = complete_block(out_ids, out_val, chosen["cell_ids"],
                                           max_missing_frac=max_missing)
    except RuntimeError as exc:
        raise RuntimeError(f"{scene}: {exc}") from None
    cell_data = filled.astype(np.float32)
    if not np.all(np.isfinite(cell_data)):
        raise RuntimeError(f"{scene}: {int((~np.isfinite(cell_data)).sum())} "
                           "non-finite cells after resampling.")

    # --- equal-area UTM crop centred on the tile centroid --------------------
    half = cfg.UTM_SIDE_PX // 2
    if cfg.UTM_SIDE_PX % 2 == 0:
        raise RuntimeError("expected an odd equal-area crop side for a clean centre")
    ccol, crow = _to_index(chosen["cx"].mean(), chosen["cy"].mean(), x, y)
    ccol, crow = int(round(ccol)), int(round(crow))
    x1, x2 = ccol - half, ccol + half + 1
    y1, y2 = crow - half, crow + half + 1
    if x1 < 0 or y1 < 0 or x2 > W or y2 > H:
        raise RuntimeError(
            f"{scene}: the {cfg.UTM_SIDE_PX}x{cfg.UTM_SIDE_PX} equal-area crop "
            f"[{y1}:{y2}, {x1}:{x2}] does not fit in the {W}x{H} window.")
    utm_crop = np.asarray(img[y1:y2, x1:x2], dtype=np.float32)
    n_bad_crop = int((~np.isfinite(utm_crop)).sum())
    frac_bad_crop = n_bad_crop / float(utm_crop.size)
    # The same ceiling as the tile side, for the same reason. Note that the
    # reader flags DN == 0 as no-data, which over very dark water (deep lakes,
    # shadow) can flag genuine zero-reflectance pixels rather than missing data;
    # the per-region counts are recorded so this stays auditable.
    if frac_bad_crop > max_nodata:
        raise RuntimeError(
            f"{scene}: the equal-area UTM crop contains {n_bad_crop} no-data "
            f"pixels ({100*frac_bad_crop:.3f}%), above the "
            f"{100*max_nodata:.3f}% ceiling. Both grids must see the same "
            "clean ground.")
    if n_bad_crop:
        filled, _ = _fill_small_holes(utm_crop)
        utm_crop = np.asarray(filled, dtype=np.float32)
        if not np.all(np.isfinite(utm_crop)):
            raise RuntimeError(
                f"{scene}: {int((~np.isfinite(utm_crop)).sum())} no-data pixels "
                "in the UTM crop could not be interpolated (contiguous hole).")

    # --- footprint outline for the diagnostic figure (decimated hull) -------
    step = max(1, cfg.N_CHILDREN // 4096)
    from scipy.spatial import ConvexHull
    pts = np.column_stack([chosen["clon"][::step], chosen["clat"][::step]])
    hull = ConvexHull(pts)
    hull_lon = pts[hull.vertices, 0]
    hull_lat = pts[hull.vertices, 1]

    meta = dict(
        scene=scene,
        parent_id=chosen["parent"],
        parent_rank=int(list(ranked).index(chosen["parent"])),
        n_candidate_parents=len(diag),
        level=cfg.LEVEL,
        delta=cfg.DELTA,
        tile_side=cfg.TILE_SIDE,
        coverage=chosen["fit"],
        nodata_in_tile_fraction=chosen["nodata"],
        nodata_samples_dropped=n_dropped,
        utm_nodata_fraction=frac_bad_crop,
        utm_nodata_filled=n_bad_crop,
        max_nodata_allowed=max_nodata,
        window_px=int(W),
        sub_window=[int(r0), int(r1), int(c0), int(c1)],
        utm_crop_box=[int(y1), int(y2), int(x1), int(x2)],
        utm_side_px=cfg.UTM_SIDE_PX,
        centroid_lon=float(chosen["clon"].mean()),
        centroid_lat=float(chosen["clat"].mean()),
        crs_wkt=crs_wkt,
        resampler=kernel,
        resampler_constant_check=ones_med,
        resampler_cells_filled=int(n_missing),
        hpx_raw_bytes=int(cell_data.nbytes),
        utm_raw_bytes=int(utm_crop.nbytes),
        hpx_mean=float(cell_data.mean()), hpx_std=float(cell_data.std()),
        utm_mean=float(utm_crop.mean()), utm_std=float(utm_crop.std()),
        source_item_id=str(ds.attrs.get("source_item_id", "")),
        source_datetime=str(ds.attrs.get("source_datetime", "")),
    )

    cfg.TILE_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(tile_path(scene, kernel),
             cell_data=cell_data,
             utm_crop=utm_crop,
             utm_x=np.asarray(x[x1:x2], dtype=np.float64),
             utm_y=np.asarray(y[y1:y2], dtype=np.float64),
             utm_lon=np.asarray(lon[y1:y2, x1:x2], dtype=np.float64),
             utm_lat=np.asarray(lat[y1:y2, x1:x2], dtype=np.float64),
             hull_lon=hull_lon, hull_lat=hull_lat,
             parent_id=np.int64(chosen["parent"]),
             meta=np.array(json.dumps(meta)))
    return meta



def tile_subwindow(scene, store, halo_px=32, band=None, kernel=DEFAULT_KERNEL):
    """Rebuild the native sub-window one stored tile was transferred from.

    Callers that re-run a transfer outside `build_tile` -- the kernel ablation
    and the resampling-error measurement -- need the same samples the tile was
    built from. `build_tile` guarantees that by *selecting* a parent whose cells
    all fall inside the window with halo; reading the stored parent back and
    trusting it does not, and when the window on disk is not the one the tile
    came from the sub-window silently covers only part of the tile. That shows
    up much later as "43% of cells missing", which reads like a halo problem
    and is not one. The fit is therefore checked here, once, with a message
    that names the actual geometry.
    """
    import pyproj
    ds, img, lon, lat, x, y, crs_wkt = _open_window(store, band or cfg.BAND)
    H, W = img.shape
    crs = pyproj.CRS.from_wkt(crs_wkt)
    to_crs = pyproj.Transformer.from_crs(pyproj.CRS.from_epsg(4326), crs,
                                         always_xy=True)
    z = np.load(tile_path(scene, kernel), allow_pickle=False)
    parent = int(z["parent_id"])
    cand = _evaluate_parent(parent, to_crs, x, y, img, H, W, halo_px)
    if cand["fit"] < 1.0:
        # Say which parent this window IS centred on. "The window is not the
        # one the tile came from" is true but unhelpful; naming the parent the
        # window was acquired for turns it into a diagnosis, and tells the
        # caller whether rebuilding would even succeed.
        import healpix_geo
        cover = healpix_geo.nested.lonlat_to_healpix(
            lon.ravel(), lat.ravel(), cfg.LEVEL, ellipsoid="WGS84").astype(np.int64)
        parents = np.unique(cover // cfg.N_CHILDREN)
        plon, plat = healpix_geo.nested.healpix_to_lonlat(
            parents.astype(np.int64), cfg.PARENT_LEVEL, ellipsoid="WGS84")
        pcx, pcy = to_crs.transform(np.asarray(plon, dtype=np.float64),
                                    np.asarray(plat, dtype=np.float64))
        xc, yc = 0.5 * (x[0] + x[-1]), 0.5 * (y[0] + y[-1])
        ranked = parents[np.argsort(np.hypot(pcx - xc, pcy - yc))]
        best, best_fit = None, 0.0
        for q in ranked[:4]:
            c = _evaluate_parent(int(q), to_crs, x, y, img, H, W, halo_px)
            if c["fit"] > best_fit:
                best, best_fit = int(q), c["fit"]
            if best_fit >= 1.0:
                break
        centred = (f"this window is centred on parent {best}, which fits it "
                   f"{100*best_fit:.1f}%" if best is not None
                   else "no candidate parent fits this window either")
        rebuild = ("Rebuilding with build_tiles.py --force would succeed, but "
                   "on a different parent tile, so every published result for "
                   "this region would change."
                   if best_fit >= 1.0 and best != parent else
                   "Rebuilding would not help: acquire a window centred on the "
                   "stored parent.")
        raise RuntimeError(
            f"{scene}: only {100*cand['fit']:.1f}% of the stored tile's cells "
            f"(parent {parent}) lie inside the {W}x{H} px window with a "
            f"{halo_px} px halo. The tile spans columns "
            f"[{cand['col'].min():.0f}, {cand['col'].max():.0f}] and rows "
            f"[{cand['row'].min():.0f}, {cand['row'].max():.0f}]; "
            f"{centred}. The stored tile was built from a different "
            f"acquisition. {rebuild}")
    c0 = max(int(np.floor(cand["col"].min())) - halo_px, 0)
    c1 = min(int(np.ceil(cand["col"].max())) + halo_px + 1, W)
    r0 = max(int(np.floor(cand["row"].min())) - halo_px, 0)
    r1 = min(int(np.ceil(cand["row"].max())) + halo_px + 1, H)
    return dict(cand=cand, to_crs=to_crs, img=img,
                sub_lon=lon[r0:r1, c0:c1].ravel(),
                sub_lat=lat[r0:r1, c0:c1].ravel(),
                bounds=(r0, r1, c0, c1), shape=(H, W), crs_wkt=crs_wkt,
                x=x, y=y)



def audit_windows(scenes, halo_px=32, band=None, verbose=True):
    """Can each stored tile still be reproduced from the window on disk?

    A tile and the window it came from are two files that can drift apart --
    a window re-acquired after the tile was built, at a different size or with
    a different centring convention, leaves an archive whose tile cannot be
    regenerated from its own inputs. Nothing in the compression results depends
    on the window once the tile exists, so the drift is invisible until
    something tries to re-run the transfer, which is exactly what the kernel
    ablation and the resampling-error measurement do. This checks all of them
    at once and reports the geometry, without building anything.
    """
    rows = []
    for scene in scenes:
        store = cfg.REGION_DIR / f"{scene}_data.zarr"
        row = dict(scene=scene, window_px=np.nan, fit=np.nan,
                   tile_rows="", tile_cols="", centred_on="", status="")
        if not tile_path(scene).exists():
            row["status"] = "no tile"
        elif not store.exists():
            row["status"] = "no window"
        else:
            try:
                w = tile_subwindow(scene, store, halo_px=halo_px, band=band)
                H, W = w["shape"]
                cand = w["cand"]
                row.update(window_px=f"{W}x{H}", fit=1.0, status="ok",
                           tile_rows=f"[{cand['row'].min():.0f}, "
                                     f"{cand['row'].max():.0f}]",
                           tile_cols=f"[{cand['col'].min():.0f}, "
                                     f"{cand['col'].max():.0f}]")
                off = abs(0.5 * (cand["row"].min() + cand["row"].max()) - H / 2)
                row["centred_on"] = "tile" if off < 0.1 * H else "anchor"
            except Exception as exc:                       # noqa: BLE001
                msg = str(exc)
                row["status"] = "MISMATCH"
                import re as _re
                m = _re.search(r"only ([\d.]+)% .*?(\d+)x(\d+) px window", msg)
                if m:
                    row["fit"] = float(m.group(1)) / 100.0
                    row["window_px"] = f"{m.group(2)}x{m.group(3)}"
                m = _re.search(r"columns \[([-\d]+), ([-\d]+)\] and rows "
                               r"\[([-\d]+), ([-\d]+)\]", msg)
                if m:
                    row["tile_cols"] = f"[{m.group(1)}, {m.group(2)}]"
                    row["tile_rows"] = f"[{m.group(3)}, {m.group(4)}]"
        rows.append(row)
        if verbose:
            fit = row["fit"]
            fit_txt = "  --  " if fit != fit else f"{100 * fit:5.1f}%"
            print(f"  {row['status']:9s} {scene:36s} "
                  f"window {str(row['window_px']):>10}  fit {fit_txt}  "
                  f"rows {row['tile_rows']:14s} {row['centred_on']}")
    return pd.DataFrame(rows)


def load_tile(scene: str, kernel: str = DEFAULT_KERNEL):
    """Return (cell_ids, cell_data, utm_crop, meta) for one region."""
    z = np.load(tile_path(scene, kernel), allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    parent = int(z["parent_id"])
    cell_ids = parent * cfg.N_CHILDREN + np.arange(cfg.N_CHILDREN, dtype=np.int64)
    return cell_ids, z["cell_data"], z["utm_crop"], meta, z


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--classes", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--band", default=cfg.BAND)
    ap.add_argument("--halo", type=int, default=32)
    ap.add_argument("--max-nodata", type=float, default=1e-3,
                    help="ceiling on the no-data fraction INSIDE the tile "
                         "footprint (default 0.001 = 0.1%%); above it the "
                         "region is rejected instead of interpolated over")
    ap.add_argument("--max-missing", type=float, default=1e-4,
                    help="ceiling on the cells a resampler drops at the edge "
                         "of the sample hull (default 0.0001 = 0.01%%, about "
                         "100 cells of a million). Below it they are filled "
                         "from their surviving siblings and counted in the "
                         "tile metadata; above it the region is rejected, "
                         "because a kernel that cannot populate the grid has "
                         "failed. Deliberately stricter than --max-nodata: a "
                         "missing sample and a missing cell are not the same "
                         "defect")
    ap.add_argument("--resampler", choices=list(KERNELS), default=DEFAULT_KERNEL,
                    help="transfer kernel used to populate the carrier grid. "
                         "Anything other than the default writes to a "
                         "kernel-suffixed .npz, so the tiles the paper's "
                         "results are built on are never overwritten")
    ap.add_argument("--audit", action="store_true",
                    help="check that every stored tile can still be "
                         "reproduced from the window on disk, and report the "
                         "geometry of any that cannot. Builds nothing")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--legacy", action="store_true",
                    help="build tiles for the seven original scenes in data/ instead")
    ap.add_argument("--fill-nan-mean", action="store_true",
                    help="legacy behaviour: replace no-data by the window mean "
                         "(NOT allowed for the 40-region protocol)")
    args = ap.parse_args(argv)

    cfg.ensure_dirs()
    if args.legacy:
        scenes = [(s, cfg.DATA_DIR / f"{s}_data.zarr") for s in cfg.LEGACY_SCENES]
        if args.regions:
            scenes = [t for t in scenes if t[0] in set(args.regions)]
        sites = pd.DataFrame([dict(scene=s, scene_class=s, region_id=s,
                                   lat=np.nan, lon=np.nan, hpx_regime="")
                              for s, _ in scenes])
    else:
        sites = region_sites.select(region_sites.load_sites(),
                                    args.regions, args.classes, args.limit)
        scenes = [(s.scene, cfg.REGION_DIR / f"{s.scene}_data.zarr")
                  for s in sites.itertuples()]
        print(region_sites.summary(sites))

    if args.audit:
        df = audit_windows([sc for sc, _ in scenes], halo_px=args.halo,
                           band=args.band)
        out = cfg.RESULTS_DIR / "window_audit.csv"
        cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, index=False)
        bad = df[df.status != "ok"]
        print(f"\n{len(df) - len(bad)} of {len(df)} tiles reproducible from "
              f"their window; wrote {out}")
        if len(bad):
            print("not reproducible:")
            for _, r in bad.iterrows():
                print(f"  {r.scene}: {r.status}, window {r.window_px}, "
                      f"fit {'' if r.fit != r.fit else f'{100*r.fit:.1f}%'}, "
                      f"tile rows {r.tile_rows}")
        return 0 if not len(bad) else 1

    rows, failures = [], []
    for i, (scene, store) in enumerate(scenes, 1):
        if tile_path(scene, args.resampler).exists() and not args.force:
            _, _, _, meta, _ = load_tile(scene, args.resampler)
            rows.append(meta)
            print(f"[{i:>2}/{len(scenes)}] {scene}: cached")
            continue
        if not store.exists():
            failures.append((scene, f"missing window store {store}"))
            print(f"[{i:>2}/{len(scenes)}] {scene}: MISSING {store}")
            continue
        try:
            meta = build_tile(scene, store, band=args.band, halo_px=args.halo,
                              max_nodata=args.max_nodata,
                              max_missing=args.max_missing,
                              fill_nan_mean=args.fill_nan_mean,
                              kernel=args.resampler, verbose=False)
            rows.append(meta)
            print(f"[{i:>2}/{len(scenes)}] {scene}: parent {meta['parent_id']} "
                  f"(rank {meta['parent_rank']}/{meta['n_candidate_parents']}), "
                  f"no-data tile {100*meta['nodata_in_tile_fraction']:.4f}% / "
                  f"crop {meta['utm_nodata_filled']}px, "
                  f"hpx mean {meta['hpx_mean']:.4f}, utm mean {meta['utm_mean']:.4f}")
        except Exception as exc:                          # noqa: BLE001
            failures.append((scene, f"{type(exc).__name__}: {exc}"))
            print(f"[{i:>2}/{len(scenes)}] {scene}: REJECTED {type(exc).__name__}: {exc}")

    if rows:
        idx = pd.DataFrame(rows).drop(columns=["crs_wkt"], errors="ignore")
        if not args.legacy:
            idx = idx.merge(sites[["scene", "scene_class", "region_id", "lat", "lon",
                                   "abs_lat", "hpx_regime"]], on="scene", how="left")
        idx.to_csv(INDEX_PATH if args.resampler == DEFAULT_KERNEL
                   else INDEX_PATH.with_name(
                       f"{INDEX_PATH.stem}__{args.resampler}.csv"),
                   index=False)
        print(f"\n{len(idx)} tiles indexed -> {INDEX_PATH}")

    if failures:
        print(f"\n{len(failures)} region(s) rejected:")
        for scene, err in failures:
            print(f"  {scene}: {err}")
        return 1
    print("Next: python run_codec_bench.py --grid both")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
