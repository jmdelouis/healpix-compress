"""Derive, per region, the HEALPix carrier tile and its equal-area UTM crop.

Replaces cell 10 of compress_test.ipynb and cell 3 of utm_codec_bench.ipynb by
a single, testable, assertion-guarded step, and extends both to the 40 shared
regions.

Guarantees enforced here (any failure rejects the region, it is never patched):
  * the chosen level-10 parent tile's 1024x1024 level-20 cells all fall inside
    the acquired native window, with a halo -- no cell is ever filled;
  * every native pixel in the tile's bounding box is finite -- no no-data is
    silently interpolated into the carrier grid;
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


def build_tile(scene: str, store: Path, band: str = cfg.BAND, halo_px: int = 32,
               n_candidates: int = 6, fill_nan_mean: bool = False,
               verbose: bool = True) -> dict:
    import healpix_geo
    import pyproj
    from healpix_resample import CloughTocherResampler

    ds, img, lon, lat, x, y, crs_wkt = _open_window(store, band)
    H, W = img.shape
    if fill_nan_mean and not np.all(np.isfinite(img)):
        img = np.where(np.isfinite(img), img, np.nanmean(img))

    crs = pyproj.CRS.from_wkt(crs_wkt) if crs_wkt else None
    if crs is None:
        raise RuntimeError(f"{scene}: the window store carries no CRS; cannot "
                           "place the HEALPix tile in projected coordinates.")
    to_crs = pyproj.Transformer.from_crs(pyproj.CRS.from_epsg(4326), crs,
                                         always_xy=True)

    # --- which level-10 parent tiles does this window cover? -----------------
    cover = healpix_geo.nested.lonlat_to_healpix(
        lon.ravel(), lat.ravel(), cfg.LEVEL, ellipsoid="WGS84").astype(np.int64)
    parents, counts = np.unique(cover // cfg.N_CHILDREN, return_counts=True)
    ranked = parents[np.argsort(-counts)][:n_candidates]

    chosen, diag = None, []
    for parent in ranked:
        cell_ids = int(parent) * cfg.N_CHILDREN + np.arange(cfg.N_CHILDREN, dtype=np.int64)
        clon, clat = healpix_geo.nested.healpix_to_lonlat(
            cell_ids, cfg.LEVEL, ellipsoid="WGS84")
        clon = np.asarray(clon, dtype=np.float64)
        clat = np.asarray(clat, dtype=np.float64)
        cx, cy = to_crs.transform(clon, clat)
        col, row = _to_index(cx, cy, x, y)
        inside = ((col >= halo_px) & (col <= W - 1 - halo_px)
                  & (row >= halo_px) & (row <= H - 1 - halo_px))
        coverage = float(inside.mean())
        diag.append((int(parent), coverage))
        if coverage >= 1.0:
            chosen = dict(parent=int(parent), cell_ids=cell_ids, clon=clon, clat=clat,
                          cx=cx, cy=cy, col=col, row=row, coverage=coverage)
            break

    if chosen is None:
        best = max(diag, key=lambda t: t[1]) if diag else (None, 0.0)
        raise RuntimeError(
            f"{scene}: no level-{cfg.PARENT_LEVEL} tile is fully inside the "
            f"{W}x{H} px window with a {halo_px} px halo "
            f"(best candidate {best[0]} covers {100*best[1]:.2f}%). "
            f"Acquire a larger window (need >= {cfg.MIN_WINDOW_PX} px + halo) "
            "rather than filling the missing cells.")

    # --- native sub-window that contains the tile, plus halo -----------------
    c0 = int(np.floor(chosen["col"].min())) - halo_px
    c1 = int(np.ceil(chosen["col"].max())) + halo_px + 1
    r0 = int(np.floor(chosen["row"].min())) - halo_px
    r1 = int(np.ceil(chosen["row"].max())) + halo_px + 1
    c0, r0 = max(c0, 0), max(r0, 0)
    c1, r1 = min(c1, W), min(r1, H)
    sub = img[r0:r1, c0:c1]
    finite_sub = float(np.isfinite(sub).mean())
    if finite_sub < 1.0:
        raise RuntimeError(
            f"{scene}: {100*(1-finite_sub):.3f}% of the native pixels under the "
            "tile are no-data. Rejecting the region: filling them would make the "
            "HEALPix side compress for free.")

    # --- resample onto the tile (Clough-Tocher, the paper's reference) -------
    op = CloughTocherResampler(lon[r0:r1, c0:c1].ravel(),
                               lat[r0:r1, c0:c1].ravel(),
                               cfg.LEVEL, out_cell_ids=chosen["cell_ids"],
                               verbose=verbose)
    res = op.resample(sub.ravel().astype(np.float64))
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
    if not np.all(np.isfinite(utm_crop)):
        raise RuntimeError(f"{scene}: the equal-area UTM crop contains no-data.")

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
        level=cfg.LEVEL,
        delta=cfg.DELTA,
        tile_side=cfg.TILE_SIDE,
        coverage=chosen["coverage"],
        finite_native_fraction=finite_sub,
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
                              fill_nan_mean=args.fill_nan_mean, verbose=False)
            rows.append(meta)
            print(f"[{i:>2}/{len(scenes)}] {scene}: parent {meta['parent_id']}, "
                  f"coverage {100*meta['coverage']:.2f}%, "
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
