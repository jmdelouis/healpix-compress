"""Acquire the 1512 x 1512 native windows for the 40 shared regions.

Why this script exists
----------------------
healpix-resample publishes 256 x 256 px (2.56 km) patches. A level-20 HEALPix
carrier tile of 1024 x 1024 cells is 6.37 km on a side and is a *rotated*
square in UTM, so it needs a native window of at least
``paper_config.MIN_WINDOW_PX`` (901) pixels whatever the rotation. This script
re-extracts a larger window at exactly the same region centres, from exactly
the same Sentinel-2 products, so the two papers share sites, products and
radiometry.

The reading itself is delegated to healpix-resample's own acquisition helper
(``tests/real_groundtruth_common_tools.py``): same COG access, same no-data
handling, same processing-baseline BOA offset, same reflectance scaling. This
script never writes into the healpix-resample bundle -- output goes to
``test/data/regions/`` in this repository.

Usage
-----
    python build_region_windows.py --list
    python build_region_windows.py --check
    python build_region_windows.py --regions paris amazon      # try a few first
    python build_region_windows.py                             # all 40
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import paper_config as cfg
import region_sites

LOG_PATH = cfg.REGION_DIR / "acquisition_log.csv"


def load_gt():
    """Import healpix-resample's acquisition helper by file path."""
    root = cfg.resample_root()
    module_path = root / "tests" / "real_groundtruth_common_tools.py"
    if not module_path.exists():
        raise FileNotFoundError(
            f"{module_path} not found -- this script reuses healpix-resample's "
            "own Sentinel-2 reader so the two papers share radiometry."
        )
    spec = importlib.util.spec_from_file_location("real_groundtruth_common_tools",
                                                  module_path)
    gt = importlib.util.module_from_spec(spec)
    sys.modules["real_groundtruth_common_tools"] = gt
    spec.loader.exec_module(gt)
    gt.OFFLINE = False          # this module only ever reads, never writes there
    return gt


def store_path(scene: str) -> Path:
    return cfg.REGION_DIR / f"{scene}_data.zarr"


def _item_by_id(gt, product_id: str):
    """Resolve one STAC item by its exact product id (no cloud/date filter)."""
    import pystac_client
    catalog = pystac_client.Client.open(gt.EARTH_SEARCH_URL)
    items = list(catalog.search(collections=[gt.EARTH_SEARCH_COLLECTION],
                                ids=[product_id], limit=1).items())
    return items[0] if items else None


def fetch_window(gt, site, window_px: int, band: str):
    """Read one native window centred on the region anchor. Returns (ds, info)."""
    import pyproj
    import xarray as xr

    lon0, lat0 = float(site.lon), float(site.lat)
    pinned = str(site.product_id or "")
    item, how = None, ""
    if pinned:
        item = _item_by_id(gt, pinned)
        how = "pinned" if item is not None else ""
    if item is None:
        coords = {"wgs84": {"lat": lat0, "lon": lon0},
                  "recommended_date": gt.REGION_DATE_WINDOW,
                  "cloud": gt.REGION_CLOUD_RELAXED}
        items = gt._search_earth_search(lon0, lat0, coords)
        if not items:
            raise RuntimeError(f"{site.scene}: no earth-search item at "
                               f"({lat0:.4f}, {lon0:.4f})")
        item = items[0]
        how = "resolved" if not pinned else "pinned-missing->resolved"

    props = item.properties
    if gt.EARTH_SEARCH_ASSET not in item.assets:
        raise KeyError(f"{site.scene}: asset {gt.EARTH_SEARCH_ASSET!r} missing "
                       f"from {item.id}")

    refl, x, y, crs = gt._read_cog_patch(
        item.assets[gt.EARTH_SEARCH_ASSET].href, lon0, lat0, window_px, props)

    to_wgs = pyproj.Transformer.from_crs(crs, pyproj.CRS.from_epsg(4326),
                                         always_xy=True)
    xx, yy = np.meshgrid(x, y)
    lon, lat = to_wgs.transform(xx, yy)

    ds = xr.Dataset(
        {band: (("y", "x"), refl)},
        coords={"x": x, "y": y,
                "longitude": (("y", "x"), lon),
                "latitude": (("y", "x"), lat),
                "spatial_ref": ((), 0)},
    )
    ds.spatial_ref.attrs["crs_wkt"] = crs.to_wkt()
    finite = int(np.isfinite(refl).sum())
    ds.attrs.update({
        "scene": site.scene,
        "scene_class": site.scene_class,
        "region_id": site.region_id,
        "anchor_lat": lat0,
        "anchor_lon": lon0,
        "window_px": int(window_px),
        "native_pixel_m": cfg.NATIVE_PIXEL_M,
        "band": band,
        "source_item_id": item.id,
        "source_collection": gt.EARTH_SEARCH_COLLECTION,
        "source_stac_endpoint": gt.EARTH_SEARCH_URL,
        "source_datetime": str(props.get("datetime", "")),
        "source_cloud_cover": str(props.get("eo:cloud_cover", "")),
        "source_processing_baseline": str(props.get("s2:processing_baseline", "")),
        "boa_offset_dn_applied": str(gt._boa_offset_dn(props)),
        "product_pinned_from_resample": pinned,
        "product_resolution": how,
        "resample_archive_doi": str(site.archive_doi or ""),
        "resample_relative_path": str(site.resample_relative_path or ""),
        "acquired_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "note": ("Window centred on the same anchor as the healpix-resample "
                 "256x256 patch, so that patch is the central crop of this one."),
    })
    info = dict(
        scene=site.scene, scene_class=site.scene_class, lat=lat0, lon=lon0,
        product_id=item.id, product_pinned=pinned, product_resolution=how,
        datetime=str(props.get("datetime", "")),
        cloud_cover=props.get("eo:cloud_cover", np.nan),
        window_px=int(window_px),
        finite_fraction=finite / float(refl.size),
        vmin=float(np.nanmin(refl)) if finite else np.nan,
        vmax=float(np.nanmax(refl)) if finite else np.nan,
        vmean=float(np.nanmean(refl)) if finite else np.nan,
    )
    return ds, info


def check(sites) -> pd.DataFrame:
    """Offline inventory of what is already on disk."""
    import zarr
    rows = []
    for site in sites.itertuples():
        p = store_path(site.scene)
        row = dict(scene=site.scene, path=str(p), present=p.exists(),
                   window_px=np.nan, product_id="", finite_fraction=np.nan)
        if p.exists():
            try:
                attrs = zarr.open_group(str(p), mode="r").attrs.asdict()
                row["window_px"] = attrs.get("window_px", np.nan)
                row["product_id"] = attrs.get("source_item_id", "")
            except Exception as exc:                     # noqa: BLE001
                row["product_id"] = f"<unreadable: {type(exc).__name__}>"
        rows.append(row)
    df = pd.DataFrame(rows)
    ok = df.present.sum()
    print(f"{ok}/{len(df)} region windows present under {cfg.REGION_DIR}")
    bad = df[df.present & (df.window_px != cfg.WINDOW_PX)]
    if len(bad):
        print(f"WARNING: {len(bad)} store(s) have a window_px != {cfg.WINDOW_PX}:")
        print(bad[["scene", "window_px"]].to_string(index=False))
    return df


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", nargs="*", default=None,
                    help="region ids or full scene ids; default: all 40")
    ap.add_argument("--classes", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--window", type=int, default=cfg.WINDOW_PX)
    ap.add_argument("--band", default=cfg.BAND)
    ap.add_argument("--force", action="store_true", help="re-download existing stores")
    ap.add_argument("--check", action="store_true", help="offline inventory only")
    ap.add_argument("--list", action="store_true", help="print the site table and exit")
    args = ap.parse_args(argv)

    if args.window < cfg.MIN_WINDOW_PX:
        raise SystemExit(
            f"--window {args.window} < {cfg.MIN_WINDOW_PX}: a {cfg.TILE_SIDE}-cell "
            f"level-{cfg.LEVEL} tile cannot fit whatever its rotation. Refusing.")

    cfg.ensure_dirs()
    sites = region_sites.select(region_sites.load_sites(),
                                args.regions, args.classes, args.limit)
    print(region_sites.summary(sites))

    if args.list:
        print(sites.to_string(index=False))
        return 0
    if args.check:
        check(sites)
        return 0

    gt = load_gt()
    gt.warn_about_aws_credentials()
    print(f"[source] {gt.EARTH_SEARCH_COLLECTION} @ {gt.EARTH_SEARCH_URL}, "
          f"asset {gt.EARTH_SEARCH_ASSET!r}, window {args.window} px")

    logs, failures = [], []
    for i, site in enumerate(sites.itertuples(), 1):
        path = store_path(site.scene)
        if path.exists() and not args.force:
            print(f"[{i:>2}/{len(sites)}] {site.scene}: cached")
            continue
        try:
            ds, info = fetch_window(gt, site, args.window, args.band)
            ds.to_zarr(path, mode="w", **gt._zarr_v2_kwargs())
            info["status"] = "ok"
            logs.append(info)
            print(f"[{i:>2}/{len(sites)}] {site.scene}: {info['product_id']} "
                  f"({info['product_resolution']}, cloud {info['cloud_cover']}%, "
                  f"finite {100*info['finite_fraction']:.2f}%)")
        except Exception as exc:                          # noqa: BLE001
            failures.append((site.scene, f"{type(exc).__name__}: {exc}"))
            logs.append(dict(scene=site.scene, status="failed",
                             error=f"{type(exc).__name__}: {exc}"))
            print(f"[{i:>2}/{len(sites)}] {site.scene}: FAILED {type(exc).__name__}: {exc}")

    if logs:
        new = pd.DataFrame(logs)
        if LOG_PATH.exists():
            new = pd.concat([pd.read_csv(LOG_PATH), new], ignore_index=True)
        new.to_csv(LOG_PATH, index=False)
        print(f"acquisition log -> {LOG_PATH}")

    if failures:
        print(f"\n{len(failures)} region(s) failed:")
        for scene, err in failures:
            print(f"  {scene}: {err}")
        return 1
    print("\nAll requested regions acquired. Next: python build_tiles.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
