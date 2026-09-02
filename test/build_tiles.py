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
    inside the window and is fully finite;
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


def tile_path(scene: str) -> Path:
    return cfg.TILE_DIR / f"{scene}.npz"


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


def _to_index(cx, cy, x, y):
    """Fractional (col, row) of projected coordinates on a regular x/y grid."""
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    return (cx - float(x[0])) / dx, (cy - float(y[0])) / dy


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
               max_nodata: float = 1e-3, fill_nan_mean: bool = False,
               verbose: bool = True) -> dict:
    import pyproj
    from healpix_resample import CloughTocherResampler

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
                   f"{W}x{H} px window with a {halo_px} px halo. A larger "
                   f"window (>= {cfg.MIN_WINDOW_PX} px + halo) or a smaller "
                   f"halo would be needed")
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

    op = CloughTocherResampler(sub_lon[good], sub_lat[good], cfg.LEVEL,
                               out_cell_ids=chosen["cell_ids"], verbose=verbose)
    res = op.resample(sub_val[good].astype(np.float64))
    out_ids = np.asarray(res.cell_ids, dtype=np.int64).ravel()
    out_val = np.asarray(res.cell_data, dtype=np.float64).ravel()
    if out_ids.size != cfg.N_CHILDREN:
        raise RuntimeError(
            f"{scene}: the resampler kept {out_ids.size} of {cfg.N_CHILDREN} "
            "cells (cells outside the sample convex hull are dropped). "
            "Increase --halo or acquire a larger window.")
    order = np.argsort(out_ids)
    if not np.array_equal(out_ids[order], chosen["cell_ids"]):
        raise RuntimeError(f"{scene}: resampled cell ids differ from the requested "
                           "contiguous child block.")
    cell_data = out_val[order].astype(np.float32)
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
    if n_bad_crop:
        raise RuntimeError(
            f"{scene}: the equal-area UTM crop contains {n_bad_crop} no-data "
            f"pixels ({100*n_bad_crop/utm_crop.size:.3f}%). Both grids must see "
            "the same clean ground.")

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
        max_nodata_allowed=max_nodata,
        window_px=int(W),
        sub_window=[int(r0), int(r1), int(c0), int(c1)],
        utm_crop_box=[int(y1), int(y2), int(x1), int(x2)],
        utm_side_px=cfg.UTM_SIDE_PX,
        centroid_lon=float(chosen["clon"].mean()),
        centroid_lat=float(chosen["clat"].mean()),
        crs_wkt=crs_wkt,
        resampler="CloughTocherResampler",
        hpx_raw_bytes=int(cell_data.nbytes),
        utm_raw_bytes=int(utm_crop.nbytes),
        hpx_mean=float(cell_data.mean()), hpx_std=float(cell_data.std()),
        utm_mean=float(utm_crop.mean()), utm_std=float(utm_crop.std()),
        source_item_id=str(ds.attrs.get("source_item_id", "")),
        source_datetime=str(ds.attrs.get("source_datetime", "")),
    )

    cfg.TILE_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(tile_path(scene),
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


def load_tile(scene: str):
    """Return (cell_ids, cell_data, utm_crop, meta) for one region."""
    z = np.load(tile_path(scene), allow_pickle=False)
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

    rows, failures = [], []
    for i, (scene, store) in enumerate(scenes, 1):
        if tile_path(scene).exists() and not args.force:
            _, _, _, meta, _ = load_tile(scene)
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
                              fill_nan_mean=args.fill_nan_mean, verbose=False)
            rows.append(meta)
            print(f"[{i:>2}/{len(scenes)}] {scene}: parent {meta['parent_id']} "
                  f"(rank {meta['parent_rank']}/{meta['n_candidate_parents']}), "
                  f"no-data in tile {100*meta['nodata_in_tile_fraction']:.4f}%, "
                  f"hpx mean {meta['hpx_mean']:.4f}, utm mean {meta['utm_mean']:.4f}")
        except Exception as exc:                          # noqa: BLE001
            failures.append((scene, f"{type(exc).__name__}: {exc}"))
            print(f"[{i:>2}/{len(scenes)}] {scene}: REJECTED {type(exc).__name__}: {exc}")

    if rows:
        idx = pd.DataFrame(rows).drop(columns=["crs_wkt"], errors="ignore")
        if not args.legacy:
            idx = idx.merge(sites[["scene", "scene_class", "region_id", "lat", "lon",
                                   "abs_lat", "hpx_regime"]], on="scene", how="left")
        idx.to_csv(INDEX_PATH, index=False)
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
