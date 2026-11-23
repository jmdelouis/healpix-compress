"""What the loosest target means for an application: NDVI from B04 and B08.

The crossing of the main comparison is at an RMSE of 1e-2 in reflectance,
and the external review of 2026-11-23 asked what that error means in use.
This script answers with the simplest index a red band enters, NDVI, and
with the stability of an NDVI classification, on both grids:

1.  `--fetch`: read the B08 (10 m NIR) window of every region from the same
    Sentinel-2 product, at the same window centre and size as its B04 window
    (window_geometry.plan_acquisition, which reproduces the B04 window exactly),
    through healpix-resample's own reader, and check that the native crop
    lands on the same pixels as the B04 tile's crop. Then build the B08
    carrier tile with the paper's Clough--Tocher transfer and check that it
    selects the same level-10 parent. B08 windows are deleted once their tile
    exists unless --keep-windows is given (the hub has little disk).

2.  `--run`: on each grid, sweep B04 and B08 with the paper's codec ladders.
    At each target, each band takes the smallest file of any codec whose
    RMSE over the raster does not exceed the target (the oracle, at measured
    points: NDVI needs decoded images, not an interpolated size). NDVI is
    computed from the two decoded bands and compared with the NDVI of the
    uncompressed bands of the same grid, over the whole raster and over the
    intersection D of the two footprints (same_domain.domain_masks). Class
    stability is the fraction of samples whose NDVI class changes, with the
    classes < 0.2, 0.2-0.5 and > 0.5. Both grids are equal-area samplings, so
    sample fractions are area fractions. The size ratios of B04, of B08 and
    of the pair are also computed with the paper's rule (interpolated oracle
    size, bracketing on both grids).

Usage (hub)
-----
    python ndvi_validation.py --fetch --workers 6
    python ndvi_validation.py --run --workers 12
    python ndvi_validation.py --summary-only

Writes data/regions_b08/, data/tiles_b08/, results/ndvi_validation_parts/,
results/ndvi_validation_sweeps.csv, results/ndvi_validation.csv and
results/ndvi_validation.json.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import time
import traceback

import numpy as np
import pandas as pd

import paper_config as cfg

REGION_B08 = cfg.DATA_DIR / "regions_b08"
TILE_B08 = cfg.DATA_DIR / "tiles_b08"
PARTS = cfg.RESULTS_DIR / "ndvi_validation_parts"
NIR_ASSET_KEYS = ("nir", "B08", "b08")
CLASS_EDGES = (0.2, 0.5)
VALID_SUM = 0.02           # NDVI is computed where B04 + B08 exceeds this


# --- 1. B08 windows and tiles -----------------------------------------------------
def b08_tile_path(scene):
    return TILE_B08 / f"{scene}.npz"


def fetch_region(scene, keep_window=False, force=False):
    import xarray as xr
    import pyproj
    import build_region_windows as brw
    import build_tiles
    import region_sites
    import window_geometry

    if b08_tile_path(scene).exists() and not force:
        return f"{scene}: cached"
    t0 = time.time()
    _ids, _data, _crop, meta4, z4 = build_tiles.load_tile(scene)
    sites = region_sites.load_sites()
    site = sites[sites.scene == scene].iloc[0]
    plan = window_geometry.plan_acquisition(site.lon, site.lat, halo_px=32)
    if int(plan["window_px"]) != int(meta4["window_px"]):
        raise RuntimeError(f"{scene}: planned window {plan['window_px']} px, "
                           f"B04 window {meta4['window_px']} px")
    gt = brw.load_gt()
    coll, item = brw._item_by_id(gt, meta4["source_item_id"])
    if item is None:
        raise RuntimeError(f"{scene}: product {meta4['source_item_id']} not found")
    key = next((k for k in NIR_ASSET_KEYS if k in item.assets), None)
    if key is None:
        raise RuntimeError(f"{scene}: no NIR asset in {list(item.assets)}")
    href = gt._public_https_href(item.assets[key].href)
    low = href.lower()
    if low.endswith(".jp2") or any(t in low for t in brw.REQUESTER_PAYS_TOKENS):
        raise RuntimeError(f"{scene}: NIR asset is not a public COG ({href})")
    refl, x, y, crs = gt._read_cog_patch(href, plan["centre_lon"], plan["centre_lat"],
                                         int(plan["window_px"]), item.properties)
    if abs(float(x[1] - x[0]) - cfg.NATIVE_PIXEL_M) > 1e-6:
        raise RuntimeError(f"{scene}: NIR asset pixel is {x[1] - x[0]} m")
    y1, y2, x1, x2 = meta4["utm_crop_box"]
    if not (np.allclose(x[x1:x2], z4["utm_x"]) and np.allclose(y[y1:y2], z4["utm_y"])):
        raise RuntimeError(f"{scene}: B08 window does not land on the B04 crop")
    finite = float(np.isfinite(refl).mean())
    to_wgs = pyproj.Transformer.from_crs(crs, pyproj.CRS.from_epsg(4326),
                                         always_xy=True)
    xx, yy = np.meshgrid(x, y)
    lon, lat = to_wgs.transform(xx, yy)
    ds = xr.Dataset({"b08": (("y", "x"), refl)},
                    coords={"x": x, "y": y,
                            "longitude": (("y", "x"), lon),
                            "latitude": (("y", "x"), lat),
                            "spatial_ref": ((), 0)})
    ds.spatial_ref.attrs["crs_wkt"] = crs.to_wkt()
    ds.attrs.update(dict(scene=scene, band="b08", source_item_id=item.id,
                         source_collection=coll, source_asset=key,
                         source_datetime=str(item.properties.get("datetime", "")),
                         boa_offset_dn_applied=str(gt._boa_offset_dn(item.properties)),
                         window_centre_lon=float(plan["centre_lon"]),
                         window_centre_lat=float(plan["centre_lat"]),
                         window_px=int(plan["window_px"]),
                         finite_fraction=finite))
    REGION_B08.mkdir(parents=True, exist_ok=True)
    store = REGION_B08 / f"{scene}_data.zarr"
    ds.to_zarr(store, mode="w", **gt._zarr_v2_kwargs())

    TILE_B08.mkdir(parents=True, exist_ok=True)
    old = cfg.TILE_DIR
    cfg.TILE_DIR = TILE_B08
    try:
        meta8 = build_tiles.build_tile(scene, store, band="b08", verbose=False)
    finally:
        cfg.TILE_DIR = old
    z8 = np.load(b08_tile_path(scene), allow_pickle=False)
    same = (int(meta8["parent_id"]) == int(meta4["parent_id"])
            and list(meta8["utm_crop_box"]) == list(meta4["utm_crop_box"])
            and np.allclose(z8["utm_x"], z4["utm_x"])
            and np.allclose(z8["utm_y"], z4["utm_y"]))
    if not same:
        b08_tile_path(scene).unlink()
        raise RuntimeError(f"{scene}: the B08 tile does not match the B04 tile's "
                           f"parent or crop")
    if not keep_window:
        shutil.rmtree(store, ignore_errors=True)
    return (f"{scene}: {item.id} [{coll}/{key}], finite {100 * finite:.3f}%, "
            f"parent {meta8['parent_id']}, {time.time() - t0:.0f} s")


def check_rebuild(scene):
    """Rebuild the B04 tile here from its window and compare with the deposit's."""
    import build_tiles
    import tempfile
    store = cfg.REGION_DIR / f"{scene}_data.zarr"
    _ids, data4, _crop, meta4, _z = build_tiles.load_tile(scene)
    old = cfg.TILE_DIR
    with tempfile.TemporaryDirectory() as tmp:
        cfg.TILE_DIR = type(old)(tmp)
        try:
            build_tiles.build_tile(scene, store, verbose=False)
            z = np.load(build_tiles.tile_path(scene), allow_pickle=False)
            d = np.asarray(z["cell_data"], np.float64) - np.asarray(data4, np.float64)
        finally:
            cfg.TILE_DIR = old
    return dict(scene=scene, max_abs_diff=float(np.abs(d).max()),
                rms_diff=float(np.sqrt(np.mean(d * d))))


# --- 2. NDVI at the targets ----------------------------------------------------------
def ndvi(red, nir):
    red = np.asarray(red, np.float64)
    nir = np.asarray(nir, np.float64)
    s = red + nir
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(np.abs(s) > 1e-6, (nir - red) / s, np.nan), s


def ndvi_class(v):
    return np.digitize(v, CLASS_EDGES)


def band_sweep(arr, codecs):
    """All sweep points of one raster, reconstructions kept."""
    import codec_lib
    a = np.asarray(arr, dtype=np.float32)
    pts = []
    for run in codec_lib.iter_runs(a, codecs=codecs, on_error="warn"):
        rec = np.asarray(run.reconstruction, dtype=np.float32)
        d = rec.astype(np.float64) - a
        pts.append(dict(codec=run.codec, config=json.dumps(run.config, sort_keys=True),
                        total_bytes=int(run.total_bytes),
                        rmse=float(np.sqrt(np.mean(d * d))), rec=rec))
    return pts


def pick(pts, target):
    """Smallest file whose RMSE does not exceed the target; and whether the
    sweep reaches beyond it (otherwise the target is met at the coarsest
    point and is a bound)."""
    ok = [p for p in pts if p["rmse"] <= target]
    if not ok:
        return None, False
    best = min(ok, key=lambda p: p["total_bytes"])
    return best, any(p["rmse"] > target for p in pts)


def ndvi_metrics(ref_r, ref_n, rec_r, rec_n, mask):
    v0, s0 = ndvi(ref_r, ref_n)
    v1, s1 = ndvi(rec_r, rec_n)
    valid = (s0 > VALID_SUM) & (s1 > VALID_SUM / 2) & np.isfinite(v0) & np.isfinite(v1)
    out = {}
    for name, m in (("full", valid), ("domain", valid & mask)):
        if not m.any():
            # deep, dark water: B04 + B08 never reaches VALID_SUM, no NDVI
            out.update({f"{k}_{name}": np.nan for k in (
                "ndvi_rmse", "ndvi_mae", "ndvi_p99", "class_change",
                "dense_frac", "sparse_frac")})
            out[f"valid_frac_{name}"] = 0.0
            continue
        e = v1[m] - v0[m]
        c0, c1 = ndvi_class(v0[m]), ndvi_class(v1[m])
        out.update({
            f"ndvi_rmse_{name}": float(np.sqrt(np.mean(e * e))),
            f"ndvi_mae_{name}": float(np.mean(np.abs(e))),
            f"ndvi_p99_{name}": float(np.percentile(np.abs(e), 99)),
            f"class_change_{name}": float(np.mean(c0 != c1)),
            f"valid_frac_{name}": float(m.sum() / max(1, (mask if name == "domain" else np.ones_like(mask)).sum())),
            f"dense_frac_{name}": float(np.mean(c0 == 2)),
            f"sparse_frac_{name}": float(np.mean(c0 == 1)),
        })
    return out


def run_region(scene, codecs, targets, force=False):
    PARTS.mkdir(parents=True, exist_ok=True)
    out = PARTS / f"{scene}.csv"
    if out.exists() and not force:
        return f"{scene}: cached"
    import support_controls as sc
    import same_domain
    import paper_stats as st
    t0 = time.time()
    reg = sc.load_region(scene)
    z8 = np.load(b08_tile_path(scene), allow_pickle=False)
    if int(z8["parent_id"]) != reg["parent"] or not np.allclose(z8["utm_x"], reg["utm_x"]):
        raise RuntimeError(f"{scene}: B08 tile does not match the B04 tile")
    utm_mask, hpx_mask = same_domain.domain_masks(reg)
    rasters = {
        "hpx": (sc.expanded(reg["cell_data"]), sc.expanded(z8["cell_data"]), hpx_mask),
        "utm": (np.asarray(reg["utm_crop"], np.float32),
                np.asarray(z8["utm_crop"], np.float32), utm_mask),
    }
    rows, sweep_rows, sizes = [], [], {}
    for grid, (red, nir, mask) in rasters.items():
        sw = {"b04": band_sweep(red, codecs), "b08": band_sweep(nir, codecs)}
        for band, pts in sw.items():
            df = pd.DataFrame([{k: v for k, v in p.items() if k != "rec"} for p in pts])
            df["scene"], df["grid"], df["band"] = scene, grid, band
            sweep_rows.append(df)
            for t in targets:
                sizes[(grid, band, t)] = (
                    st.best_of_codecs_size(df, scene, codecs, t),
                    st.oracle_bracketed(df, scene, codecs, t, "rmse"))
        for t in targets:
            pr, br = pick(sw["b04"], t)
            pn, bn = pick(sw["b08"], t)
            row = dict(scene=scene, grid=grid, target=t,
                       reached_b04=br, reached_b08=bn,
                       std_b04=float(np.std(red)), std_b08=float(np.std(nir)))
            if pr is not None and pn is not None:
                row.update(codec_b04=pr["codec"], codec_b08=pn["codec"],
                           bytes_b04=pr["total_bytes"], bytes_b08=pn["total_bytes"],
                           rmse_b04=pr["rmse"], rmse_b08=pn["rmse"],
                           **ndvi_metrics(red, nir, pr["rec"], pn["rec"], mask))
            rows.append(row)
        del sw
    res = pd.DataFrame(rows)
    for band in ("b04", "b08"):
        for t in targets:
            (hb, bh), (ub, bu) = sizes[("hpx", band, t)], sizes[("utm", band, t)]
            sel = res.target == t
            res.loc[sel, f"ratio_{band}"] = hb / ub if np.isfinite(hb) and np.isfinite(ub) else np.nan
            res.loc[sel, f"ratio_{band}_bracketed"] = bool(bh and bu)
            res.loc[sel, f"hpx_size_{band}"] = hb
            res.loc[sel, f"utm_size_{band}"] = ub
    pd.concat(sweep_rows, ignore_index=True).to_csv(PARTS / f"{scene}__sweeps.csv", index=False)
    tmp = out.with_suffix(".tmp")
    res.to_csv(tmp, index=False)
    tmp.replace(out)
    return f"{scene}: {time.time() - t0:.0f} s"


def _fetch_worker(a):
    scene, keep, force = a
    try:
        return fetch_region(scene, keep, force)
    except Exception:                                         # noqa: BLE001
        return f"{scene}: FAILED\n{traceback.format_exc()[-1500:]}"


def _run_worker(a):
    scene, codecs, targets, force = a
    try:
        return run_region(scene, codecs, targets, force)
    except Exception:                                         # noqa: BLE001
        return f"{scene}: FAILED\n{traceback.format_exc()[-1500:]}"


def _pool(fn, jobs, workers):
    if workers > 1:
        from multiprocessing import get_context
        with get_context("spawn").Pool(workers) as pool:
            for msg in pool.imap_unordered(fn, jobs):
                print(msg, flush=True)
    else:
        for j in jobs:
            print(fn(j), flush=True)


# --- summary -------------------------------------------------------------------------
def summarize(r, targets):
    import paper_stats as st
    out = {"targets": [float(t) for t in targets], "class_edges": list(CLASS_EDGES),
           "valid_sum": VALID_SUM, "n_regions": int(r.scene.nunique()), "rows": []}
    for t in targets:
        d = r[r.target == t]
        # the paper's rule: on both grids, the oracle codec brackets the
        # target for B04 and for B08 (same regions as the size ratios)
        good = {s for s, g in d.groupby("scene")
                if len(g) == 2 and g.reached_b04.astype(bool).all()
                and g.reached_b08.astype(bool).all()
                and g.ratio_b04_bracketed.astype(bool).all()
                and g.ratio_b08_bracketed.astype(bool).all()
                and g.ndvi_rmse_full.notna().all()}
        dd = d[d.scene.isin(good)]
        row = dict(target=float(t), n=len(good))
        for grid in ("utm", "hpx"):
            g = dd[dd.grid == grid]
            for col in ("ndvi_rmse_full", "ndvi_rmse_domain", "ndvi_mae_full",
                        "ndvi_p99_full", "class_change_full", "class_change_domain",
                        "rmse_b04", "rmse_b08"):
                v = g[col].astype(float)
                row[f"{grid}_{col}"] = [float(v.quantile(.25)), float(v.median()),
                                        float(v.quantile(.75))] if len(v) else None
        h = dd[dd.grid == "hpx"].set_index("scene")
        u = dd[dd.grid == "utm"].set_index("scene")
        for col in ("ndvi_rmse_full", "ndvi_rmse_domain", "class_change_full",
                    "class_change_domain"):
            q = (h[col] / u[col]).replace([np.inf, -np.inf], np.nan).dropna().values
            test = st.paired_ratio_test(q)
            row[f"hpx_over_utm_{col}"] = dict(n=test["n"], median=test["median"],
                                              q25=test["q25"], q75=test["q75"],
                                              p=test["p"])
        # the paper's size ratios on each band and on the pair
        u_all = d[d.grid == "utm"].drop_duplicates("scene").set_index("scene")
        for band in ("b04", "b08"):
            sel = u_all[u_all[f"ratio_{band}_bracketed"].astype(bool)]
            test = st.paired_ratio_test(sel[f"ratio_{band}"].values)
            row[f"size_ratio_{band}"] = dict(n=test["n"], median=test["median"],
                                             p=test["p"], p_sign=test["p_sign"])
        both = u_all[u_all.ratio_b04_bracketed.astype(bool)
                     & u_all.ratio_b08_bracketed.astype(bool)]
        pair = ((both.hpx_size_b04 + both.hpx_size_b08)
                / (both.utm_size_b04 + both.utm_size_b08)).values
        test = st.paired_ratio_test(pair)
        row["size_ratio_pair"] = dict(n=test["n"], median=test["median"],
                                      p=test["p"], p_sign=test["p_sign"])
        out["rows"].append(row)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--summary-only", action="store_true")
    ap.add_argument("--check-rebuild", nargs="*", default=None,
                    help="rebuild these regions' B04 tiles from their windows and "
                         "compare with the deposited tiles")
    ap.add_argument("--keep-windows", action="store_true")
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--codecs", nargs="*", default=cfg.COMPRESSOR_CODECS)
    ap.add_argument("--targets", type=float, nargs="*", default=cfg.TARGET_RMSE)
    args = ap.parse_args(argv)

    import region_sites
    scenes = list(region_sites.select(region_sites.load_sites(), args.regions,
                                      None, None).scene)
    if args.check_rebuild is not None:
        for s in args.check_rebuild:
            print(check_rebuild(s), flush=True)
    if args.fetch:
        _pool(_fetch_worker, [(s, args.keep_windows, args.force) for s in scenes],
              args.workers)
    if args.run:
        ready = [s for s in scenes if b08_tile_path(s).exists()]
        _pool(_run_worker, [(s, args.codecs, args.targets, args.force) for s in ready],
              args.workers)
    if not (args.run or args.summary_only):
        return 0
    parts = [PARTS / f"{s}.csv" for s in scenes if (PARTS / f"{s}.csv").exists()]
    if parts:
        r = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
        r.to_csv(cfg.RESULTS_DIR / "ndvi_validation.csv", index=False)
        sw = [PARTS / f"{s}__sweeps.csv" for s in scenes
              if (PARTS / f"{s}__sweeps.csv").exists()]
        pd.concat([pd.read_csv(p) for p in sw], ignore_index=True).to_csv(
            cfg.RESULTS_DIR / "ndvi_validation_sweeps.csv", index=False)
    elif (cfg.RESULTS_DIR / "ndvi_validation.csv").exists():
        # a checkout with the assembled results only (the parts stay on the
        # machine that ran the sweeps)
        r = pd.read_csv(cfg.RESULTS_DIR / "ndvi_validation.csv")
    else:
        print("no results")
        return 1
    summary = summarize(r, args.targets)
    summary["codecs"] = list(args.codecs)
    (cfg.RESULTS_DIR / "ndvi_validation.json").write_text(json.dumps(summary, indent=1))
    for row in summary["rows"]:
        print(f"target {row['target']:g}: n={row['n']}")
        for grid in ("utm", "hpx"):
            print(f"   {grid}: NDVI RMSE {row[f'{grid}_ndvi_rmse_full']}, "
                  f"class change {row[f'{grid}_class_change_full']}, "
                  f"band RMSE b04 {row[f'{grid}_rmse_b04']}, b08 {row[f'{grid}_rmse_b08']}")
        print(f"   hpx/utm NDVI RMSE {row['hpx_over_utm_ndvi_rmse_full']}")
        print(f"   size ratios b04 {row['size_ratio_b04']}, b08 {row['size_ratio_b08']}, "
              f"pair {row['size_ratio_pair']}")
    return 0


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    raise SystemExit(main())
