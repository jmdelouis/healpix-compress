"""Is a stored tile internally consistent with its own UTM crop?

Why this matters before rebuilding anything
-------------------------------------------
`window_radiometry.py` found that the Lake Geneva Clough--Tocher tile does not
match the window on disk, while that region's nearest, linear and PSF tiles do.
So the window is current and that one tile is not. What that costs depends on a
question the window audit cannot answer: each tile's `.npz` also holds the
equal-area UTM crop the compression results compare it against, written in the
same call, and the two are only comparable if they came from the same
acquisition.

    crop matches the window   the crop was re-extracted later than the tile,
                              so every published number for this region sets
                              an old carrier against a new native square. The
                              comparison is broken and the tile must be
                              rebuilt.

    crop does not match       tile and crop are both from the earlier
                              acquisition and agree with each other. Nothing
                              published is wrong: the region is simply a
                              different day's scene, which is no more a
                              problem than the forty regions being forty
                              different days. Only the measurements that
                              re-run a transfer from the current window are
                              affected, and those already exclude it.

The test is exact, not statistical. `meta["utm_crop_box"]` records the pixel
box the crop was cut from, so the stored crop and the window's pixels at that
box must agree to the bit, up to the few no-data pixels `build_tile`
interpolates (`meta["utm_nodata_filled"]`, normally zero).

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
    _ds, img, _lon, _lat, _x, _y, _crs = build_tiles._open_window(
        store, band or cfg.BAND)
    win = np.asarray(img[y1:y2, x1:x2], dtype=np.float64)

    row = dict(scene=scene, kernel=kernel,
               crop_box=f"{y1}:{y2},{x1}:{x2}",
               nodata_filled=int(meta.get("utm_nodata_filled", 0)),
               shape_match=bool(win.shape == crop.shape))
    if not row["shape_match"]:
        row.update(identical=False, n_differ=-1, max_abs_diff=np.nan,
                   verdict="SHAPE MISMATCH")
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
            print(f"  {flag} {scene} [{kernel}]: {r['n_differ']} pixel(s) "
                  f"differ, max {r['max_abs_diff']:.2e}  -> {r['verdict']}")
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
