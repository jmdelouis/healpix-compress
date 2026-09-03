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
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import paper_config as cfg
import region_sites

LOG_PATH = cfg.REGION_DIR / "acquisition_log.csv"


class NoCoveringProduct(RuntimeError):
    """No Sentinel-2 product covers the requested window at this anchor."""


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


COG_ASSET_KEYS = ("red", "B04", "b04")
# The requester-pays JP2 archive. earth-search sometimes hands out `red` hrefs
# pointing there instead of the public COG bucket; those reads fail with
# RasterioIOError even though the item itself is perfectly valid. This is what
# defeated the same four regions in the companion healpix-resample run
# (notebooks/tables/real_groundtruth_multiregion_products.csv: four `failed`
# rows, 9 lattice positions each, same hrefs), so it is a property of the
# catalogue for these sites, not a transient error worth retrying.
REQUESTER_PAYS_TOKENS = ("sentinel-s2-l2a.s3", "//sentinel-s2-l2a/")
DEFAULT_COLLECTIONS = ("sentinel-2-l2a", "sentinel-2-c1-l2a", "sentinel-2-pre-c1-l2a")
_EOPF_ID = re.compile(r"S2(\w)_MSIL2A_(\d{8})T\d{6}_N\d+_R\d+_T(\w{5})")


def _catalog(gt):
    import pystac_client
    return pystac_client.Client.open(gt.EARTH_SEARCH_URL)


def _usable_asset(gt, item):
    """Return (asset_key, href, reason_unusable) for one item."""
    for key in COG_ASSET_KEYS:
        asset = item.assets.get(key)
        if asset is None:
            continue
        href = gt._public_https_href(asset.href)
        low = href.lower()
        if low.endswith(".jp2"):
            return key, href, "JP2 asset, not a COG"
        if any(tok in low for tok in REQUESTER_PAYS_TOKENS):
            return key, href, "asset in the requester-pays sentinel-s2-l2a archive"
        return key, href, ""
    return None, None, f"no {'/'.join(COG_ASSET_KEYS)} asset"


def _window_ring_lonlat(lon0, lat0, window_px, n_per_side=8, margin_px=8):
    """Approximate lon/lat samples along the window's boundary, plus its centre.

    A local flat-Earth conversion is enough here: the test is containment in a
    granule footprint, with a margin of a few pixels, over a window under 20 km.
    """
    half_m = (window_px / 2.0 + margin_px) * cfg.NATIVE_PIXEL_M
    dlat = half_m / 111320.0
    dlon = half_m / (111320.0 * max(0.15, math.cos(math.radians(lat0))))
    pts = [(lon0, lat0)]
    for i in range(n_per_side + 1):
        t = -1.0 + 2.0 * i / n_per_side
        pts += [(lon0 + t * dlon, lat0 - dlat), (lon0 + t * dlon, lat0 + dlat),
                (lon0 - dlon, lat0 + t * dlat), (lon0 + dlon, lat0 + t * dlat)]
    return pts


def _point_in_ring(px, py, ring) -> bool:
    inside = False
    n = len(ring)
    for i in range(n):
        x1, y1 = ring[i]
        x2, y2 = ring[(i + 1) % n]
        if (y1 > py) != (y2 > py):
            xin = (x2 - x1) * (py - y1) / (y2 - y1) + x1
            if px < xin:
                inside = not inside
    return inside


def _footprint_covers(geom, pts):
    """Does a GeoJSON footprint contain every test point? None if unknown.

    Sentinel-2 granules at the edge of a swath are only partially filled: the
    STAC item's geometry is the *valid data* polygon, not the 110 km MGRS
    square. Testing it before downloading turns 'read 1800x1800 px, discover a
    no-data wedge, discard' into a metadata comparison -- and it is the only
    thing that distinguishes a fixed granule edge (identical valid fraction for
    every acquisition of that tile) from cloud (which varies).
    """
    if not geom:
        return None
    kind = geom.get("type")
    if kind == "Polygon":
        polys = [geom.get("coordinates") or []]
    elif kind == "MultiPolygon":
        polys = geom.get("coordinates") or []
    else:
        return None
    for px, py in pts:
        ok = False
        for poly in polys:
            if not poly:
                continue
            outer = [(c[0], c[1]) for c in poly[0]]
            if _point_in_ring(px, py, outer):
                in_hole = any(_point_in_ring(px, py, [(c[0], c[1]) for c in ring])
                              for ring in poly[1:])
                if not in_hole:
                    ok = True
                    break
        if not ok:
            return False
    return True


def _search_collections(gt, lon0, lat0, collections, cloud_max, date_window, limit=200):
    """Cloud-sorted items over several collections, deterministic ordering."""
    cat = _catalog(gt)
    found = []
    for coll in collections:
        try:
            res = cat.search(
                collections=[coll],
                intersects={"type": "Point", "coordinates": [float(lon0), float(lat0)]},
                datetime=date_window,
                query={"eo:cloud_cover": {"lt": float(cloud_max)}},
                limit=limit)
            found += [(coll, it) for it in res.items()]
        except Exception as exc:                          # noqa: BLE001
            print(f"    [{coll}] search failed: {type(exc).__name__}: {exc}")
    found.sort(key=lambda t: (t[1].properties.get("eo:cloud_cover", 100.0),
                              str(t[1].properties.get("datetime", "")), t[1].id))
    return found


def _item_by_id(gt, product_id, collections=DEFAULT_COLLECTIONS):
    """Resolve one STAC item by its exact id (no cloud or date filter)."""
    cat = _catalog(gt)
    for coll in collections:
        try:
            items = list(cat.search(collections=[coll], ids=[product_id],
                                    limit=1).items())
        except Exception:                                 # noqa: BLE001
            continue
        if items:
            return coll, items[0]
    return None, None


def resample_store_item_id(site) -> str:
    """The product id recorded inside healpix-resample's own frozen patch.

    Only available when the Zenodo bundle has been extracted in that checkout;
    it is the one place the four EOPF-sourced regions carry a product id at all.
    """
    rel = str(getattr(site, "resample_relative_path", "") or "")
    if not rel:
        return ""
    try:
        import zarr
        path = cfg.resample_root() / rel
        if not path.exists():
            return ""
        return str(zarr.open_group(str(path), mode="r").attrs.asdict()
                   .get("source_item_id", "") or "")
    except Exception:                                     # noqa: BLE001
        return ""


def _same_acquisition(gt, item_id, lon0, lat0, collections):
    """Candidates for the same acquisition as an EOPF product id.

    EOPF ids (S2A_MSIL2A_<date>T<time>_N....._T31TGM_...) and earth-search ids
    (S2A_31TGM_<date>_0_L2A) name the same granule differently, so the match is
    made on (MGRS tile, acquisition day) rather than on the id itself.
    """
    m = _EOPF_ID.search(item_id or "")
    if not m:
        return []
    _sat, date, tile = m.groups()
    day = f"{date[:4]}-{date[4:6]}-{date[6:]}"
    cands = _search_collections(gt, lon0, lat0, collections, 100.0, f"{day}/{day}")
    tile_l = tile.lower()
    same = [(c, it) for c, it in cands
            if tile_l in it.id.lower()
            or str(it.properties.get("grid:code", "")).lower().endswith(tile_l)]
    return same or cands


def fetch_window(gt, site, window_px: int, band: str, *,
                 collections=DEFAULT_COLLECTIONS, max_candidates: int = 12,
                 min_finite: float = 0.999, cloud_ladder=None):
    """Read one native window centred on the region anchor. Returns (ds, info).

    Candidates are tried in order -- the product pinned by healpix-resample, then
    the same acquisition as its EOPF patch when that bundle is available, then
    the cloud-sorted catalogue search -- and each is skipped if its asset is not
    a readable public COG. Taking `items[0]` unconditionally is what made four
    regions fail in the companion study.
    """
    import pyproj
    import xarray as xr

    lon0, lat0 = float(site.lon), float(site.lat)
    pinned = str(site.product_id or "")
    candidates, seen = [], set()

    def _add(coll, item, how):
        if item is None or item.id in seen:
            return
        seen.add(item.id)
        candidates.append((coll, item, how))

    if pinned:
        _add(*_item_by_id(gt, pinned, collections), "pinned")
    eopf_id = resample_store_item_id(site)
    if eopf_id and not pinned:
        for coll, item in _same_acquisition(gt, eopf_id, lon0, lat0, collections):
            _add(coll, item, f"eopf-acquisition-match({eopf_id})")
    for cloud_max in (cloud_ladder or gt.REGION_CLOUD_LADDER):
        for coll, item in _search_collections(gt, lon0, lat0, collections,
                                              cloud_max, gt.REGION_DATE_WINDOW):
            _add(coll, item, f"resolved(cloud<{cloud_max:g})")
        if len(candidates) >= max_candidates:
            break
    if not candidates:
        raise RuntimeError(f"{site.scene}: no catalogue item at "
                           f"({lat0:.4f}, {lon0:.4f}) in {collections}")

    ring = _window_ring_lonlat(lon0, lat0, window_px)
    attempts, chosen, n_uncovered = [], None, 0
    for coll, item, how in candidates:
        if len(attempts) - n_uncovered >= max_candidates:
            break
        covers = _footprint_covers(getattr(item, "geometry", None), ring)
        if covers is False:
            n_uncovered += 1
            attempts.append(f"{item.id}[{coll}]: granule footprint does not "
                            f"cover the {window_px} px window")
            continue
        key, href, problem = _usable_asset(gt, item)
        if problem:
            attempts.append(f"{item.id}[{coll}]: {problem}")
            continue
        try:
            refl, x, y, crs = gt._read_cog_patch(href, lon0, lat0, window_px,
                                                item.properties)
        except Exception as exc:                          # noqa: BLE001
            attempts.append(f"{item.id}[{coll}]: {type(exc).__name__}: "
                            f"{str(exc)[:90]}")
            continue
        # A window straddling a granule edge is mostly valid but carries a
        # no-data wedge that would later force build_tiles.py to reject the
        # region. Rejecting the *product* here instead, and moving to the next
        # acquisition, is both cheaper and better: the criterion is data
        # availability, evaluated before anything is measured.
        finite_frac = float(np.isfinite(refl).mean())
        if finite_frac < min_finite:
            attempts.append(f"{item.id}[{coll}]: only {100*finite_frac:.2f}% "
                            f"valid pixels (granule edge or cloud mask)")
            continue
        chosen = (coll, item, how, key, refl, x, y, crs)
        break

    if chosen is None:
        summary = (f"{site.scene}: none of the {len(candidates)} candidate "
                   f"products covers a clean {window_px} px window")
        if n_uncovered:
            summary += (f" ({n_uncovered} rejected on footprint alone -- the "
                        f"anchor sits near a granule edge)")
        raise NoCoveringProduct(summary + ".\n      "
                                + "\n      ".join(attempts[:6]))

    coll, item, how, asset_key, refl, x, y, crs = chosen
    props = item.properties
    if attempts:
        print(f"    skipped {len(attempts)} unusable candidate(s) before "
              f"{item.id}: {attempts[0]}")

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
        "source_collection": coll,
        "source_asset": asset_key,
        "source_stac_endpoint": gt.EARTH_SEARCH_URL,
        "source_datetime": str(props.get("datetime", "")),
        "source_cloud_cover": str(props.get("eo:cloud_cover", "")),
        "source_processing_baseline": str(props.get("s2:processing_baseline", "")),
        "boa_offset_dn_applied": str(gt._boa_offset_dn(props)),
        "product_pinned_from_resample": pinned,
        "product_resolution": how,
        "rejected_candidates": " | ".join(attempts),
        "resample_archive_doi": str(site.archive_doi or ""),
        "resample_relative_path": str(site.resample_relative_path or ""),
        "acquired_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "note": ("Window centred on the same anchor as the healpix-resample "
                 "256x256 patch, so that patch is the central crop of this one."),
    })
    info = dict(
        scene=site.scene, scene_class=site.scene_class, lat=lat0, lon=lon0,
        product_id=item.id, collection=coll, asset=asset_key,
        product_pinned=pinned, product_resolution=how,
        n_rejected_candidates=len(attempts),
        rejected_candidates=" | ".join(attempts),
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
    ap.add_argument("--collections", nargs="*", default=list(DEFAULT_COLLECTIONS),
                    help="STAC collections to try, in order")
    ap.add_argument("--max-candidates", type=int, default=12,
                    help="how many catalogue items to try before giving up")
    ap.add_argument("--window-fallback", type=int, nargs="*", default=[],
                    help="smaller windows to fall back to when no product "
                         f"covers the requested one. Anything below "
                         f"{cfg.SAFE_WINDOW_PX} px + halo is refused: a smaller "
                         "window can only hold a tile when the anchor happens "
                         "to sit near a tile centre, which is luck, not a "
                         "protocol")
    ap.add_argument("--min-finite", type=float, default=0.999,
                    help="reject a product whose window has fewer valid pixels "
                         "than this fraction (default 0.999) and try the next")
    ap.add_argument("--cloud-max", type=float, nargs="*", default=None,
                    help="override the cloud-cover ladder, e.g. --cloud-max 5 10")
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

    floor = cfg.MIN_WINDOW_PX + 64
    windows = [args.window] + [w for w in sorted(args.window_fallback, reverse=True)
                               if w < args.window]
    bad = [w for w in windows if w < floor]
    if bad:
        raise SystemExit(
            f"window size(s) {bad} below the {floor} px floor. A window is only "
            f"guaranteed to contain a whole level-{cfg.PARENT_LEVEL} tile from "
            f"{cfg.SAFE_WINDOW_PX} px (+ halo): the anchor is an arbitrary point "
            f"of the HEALPix grid and can land on a four-tile junction, where "
            f"the nearest tile centre is a half-diagonal "
            f"({cfg.TILE_HALF_DIAGONAL_M:.0f} m) away.")

    gt = load_gt()
    gt.warn_about_aws_credentials()
    print(f"[source] {gt.EARTH_SEARCH_URL} | collections {args.collections} | "
          f"assets {COG_ASSET_KEYS} | window {args.window} px | "
          f"up to {args.max_candidates} candidates per region")
    if len(windows) > 1:
        print(f"[window] {windows[0]} px, falling back to {windows[1:]} px "
              f"where no product covers it")

    logs, failures = [], []
    for i, site in enumerate(sites.itertuples(), 1):
        path = store_path(site.scene)
        if path.exists() and not args.force:
            print(f"[{i:>2}/{len(sites)}] {site.scene}: cached")
            continue
        try:
            ds = info = None
            last = None
            for w in windows:
                try:
                    ds, info = fetch_window(gt, site, w, args.band,
                                            collections=tuple(args.collections),
                                            max_candidates=args.max_candidates,
                                            min_finite=args.min_finite,
                                            cloud_ladder=args.cloud_max)
                    break
                except NoCoveringProduct as exc:
                    last = exc
                    if w != windows[-1]:
                        print(f"    no product covers a {w} px window here; "
                              f"falling back to {windows[windows.index(w)+1]} px")
            if ds is None:
                raise last
            ds.to_zarr(path, mode="w", **gt._zarr_v2_kwargs())
            info["status"] = "ok"
            logs.append(info)
            print(f"[{i:>2}/{len(sites)}] {site.scene}: {info['product_id']} "
                  f"({info['product_resolution']}, cloud {info['cloud_cover']}%, "
                  f"finite {100*info['finite_fraction']:.2f}%, "
                  f"window {info['window_px']} px)")
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
