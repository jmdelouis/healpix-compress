"""Does a stored tile's UTM crop still hold the ground the window holds?

Why this exists
---------------
Each tile's `.npz` also carries the equal-area UTM crop the compression results
compare it against, written in the same call. The two are therefore consistent
with each other by construction; what this asks is whether the imagery they
hold is the imagery the window on disk holds now, which matters for every
measurement that sets one tile against another or rebuilds a transfer.

The obvious test is wrong. `meta["utm_crop_box"]` records pixel indices into
the window the tile was built from, and several of these windows have since
been re-acquired at a different size, hence with a different pixel origin. The
same indices then point at different ground, and comparing there reports a
mismatch on 35 of 40 regions that is an artefact of the indexing and says
nothing about the imagery. The crop box is therefore recomputed from the
tile's own centroid --- a position on the Earth, which survives any re-cutting
of the window --- exactly as `build_tile` derives it.

Read the result alongside `window_radiometry.py`, which asks the same question
of the carrier rather than of the crop and is immune to the same trap, being
geolocated throughout. A region that fails both has imagery that genuinely
differs; a region that fails only the recorded-box form of this test has merely
been re-windowed.

Usage
-----
    python tile_crop_consistency.py
    python tile_crop_consistency.py --regions region__water__lake_geneva
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

import build_tiles
import paper_config as cfg
import region_sites


def region_row(scene, kernel, band=None):
    z = np.load(build_tiles.tile_path(scene, kernel), allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    crop = np.asarray(z["utm_crop"], dtype=np.float64)
    y1, y2, x1, x2 = meta["utm_crop_box"]

    store = cfg.REGION_DIR / f"{scene}_data.zarr"
    _ds, img, _lon, _lat, x, y, crs_wkt = build_tiles._open_window(
        store, band or cfg.BAND)
    H, W = img.shape

    # The recorded box is a pair of pixel indices into the window the tile was
    # built from. A window re-acquired at a different size has a different
    # pixel origin, so the same indices then point at different ground and the
    # comparison fails for a reason that has nothing to do with the imagery.
    # The box is therefore recomputed here from the tile's own centroid, which
    # is a position on the Earth and survives any re-cutting of the window --
    # exactly as `build_tile` derives it.
    import pyproj
    import healpix_geo
    parent = int(meta["parent_id"])
    ids = parent * cfg.N_CHILDREN + np.arange(cfg.N_CHILDREN, dtype=np.int64)
    clon, clat = healpix_geo.nested.healpix_to_lonlat(ids, cfg.LEVEL,
                                                      ellipsoid="WGS84")
    to_crs = pyproj.Transformer.from_crs(
        pyproj.CRS.from_epsg(4326),
        pyproj.CRS.from_wkt(crs_wkt or meta.get("crs_wkt", "")), always_xy=True)
    cx, cy = to_crs.transform(np.asarray(clon, dtype=np.float64),
                              np.asarray(clat, dtype=np.float64))
    ccol, crow = build_tiles._to_index(float(np.mean(cx)), float(np.mean(cy)),
                                       x, y)
    half = cfg.UTM_SIDE_PX // 2
    ccol, crow = int(round(ccol)), int(round(crow))
    gy1, gy2 = crow - half, crow + half + 1
    gx1, gx2 = ccol - half, ccol + half + 1
    box_shifted = (gy1, gy2, gx1, gx2) != (y1, y2, x1, x2)
    y1, y2, x1, x2 = gy1, gy2, gx1, gx2
    win = np.asarray(img[y1:y2, x1:x2], dtype=np.float64)

    row = dict(scene=scene, kernel=kernel,
               crop_box=f"{y1}:{y2},{x1}:{x2}",
               recorded_box=",".join(str(v) for v in meta["utm_crop_box"]),
               box_shifted=bool(box_shifted),
               crop_shape=f"{crop.shape[0]}x{crop.shape[1]}",
               slice_shape=f"{win.shape[0]}x{win.shape[1]}",
               window_shape=f"{H}x{W}",
               window_px_meta=int(meta.get("window_px", -1)),
               utm_side_meta=int(meta.get("utm_side_px", -1)),
               utm_side_now=int(cfg.UTM_SIDE_PX),
               box_inside_window=bool(0 <= y1 and y2 <= H and 0 <= x1 and x2 <= W),
               nodata_filled=int(meta.get("utm_nodata_filled", 0)),
               shape_match=bool(win.shape == crop.shape))
    if not row["shape_match"]:
        # Two different failures look alike here and must not be confused.
        # Either the crop was cut at a side this pipeline no longer uses, or
        # the box is valid but falls outside the window now on disk, which
        # numpy silently clips -- and that means the window was re-acquired
        # smaller than the one the tile was built from.
        if not row["box_inside_window"]:
            why = (f"the tile's own centroid falls too near the edge of the "
                   f"{H}x{W} window now on disk for a {cfg.UTM_SIDE_PX} px "
                   f"crop (the tile was built from a {meta.get('window_px', '?')} "
                   f"px window)")
        elif crop.shape[0] != cfg.UTM_SIDE_PX:
            why = (f"the stored crop is {crop.shape[0]} px a side, not the "
                   f"{cfg.UTM_SIDE_PX} px this pipeline cuts")
        else:
            why = "shapes disagree for a reason neither test explains"
        row.update(identical=False, n_differ=-1, max_abs_diff=np.nan,
                   verdict=f"SHAPE MISMATCH: {why}")
        return row
    good = np.isfinite(win) & np.isfinite(crop)
    diff = np.abs(win[good] - crop[good])
    n_differ = int((diff > 0).sum())
    # Anything the hole filler touched is allowed to differ; nothing else is.
    tolerated = int(meta.get("utm_nodata_filled", 0))
    identical = n_differ <= tolerated
    row.update(identical=bool(identical), n_differ=n_differ,
               max_abs_diff=float(diff.max()) if diff.size else 0.0,
               frac_differ=float(n_differ / max(1, good.sum())),
               verdict=("crop comes from THIS window" if identical
                        else "crop does NOT come from this window"))
    return row


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--kernels", nargs="*", default=None,
                    help="default: every kernel whose tile exists")
    ap.add_argument("--band", default=None)
    args = ap.parse_args(argv)

    sites = region_sites.select(region_sites.load_sites(), args.regions,
                                None, None)
    kernels = args.kernels or list(build_tiles.KERNELS)
    rows, failures = [], []
    for scene in sites.scene:
        for kernel in kernels:
            if not build_tiles.tile_path(scene, kernel).exists():
                continue
            try:
                r = region_row(scene, kernel, band=args.band)
            except Exception as exc:                          # noqa: BLE001
                failures.append((f"{scene} [{kernel}]",
                                 f"{type(exc).__name__}: {exc}"))
                continue
            rows.append(r)
            flag = "ok " if r["identical"] else "NO "
            if r["shape_match"]:
                print(f"  {flag} {scene} [{kernel}]: {r['n_differ']} pixel(s) "
                      f"differ, max {r['max_abs_diff']:.2e}  -> {r['verdict']}")
            else:
                print(f"  {flag} {scene} [{kernel}]: stored crop "
                      f"{r['crop_shape']}, window slice {r['slice_shape']} at "
                      f"box {r['crop_box']} of a {r['window_shape']} window\n"
                      f"      -> {r['verdict']}")
    if not rows:
        for s, why in failures:
            print(f"  {s}: {why}")
        raise SystemExit("nothing could be checked")
    df = pd.DataFrame(rows)
    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    tag = "__subset" if args.regions else ""
    out = cfg.RESULTS_DIR / f"tile_crop_consistency{tag}.csv"
    df.to_csv(out, index=False)
    bad = df[~df.identical]
    print(f"\n{len(df)} tile(s) checked, {len(bad)} whose crop is not from "
          f"the window on disk")
    print(f"wrote {out}")
    if failures:
        print(f"\n{len(failures)} could not be checked:")
        for s, why in failures:
            print(f"  {s}: {why}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
