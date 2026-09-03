"""Run the rate-distortion sweep on both grids, over the 40 shared regions.

Replaces the sweep cells of hpx_codec_bench.ipynb and utm_codec_bench.ipynb by
one resumable command. Every (region, grid) pair is written to its own part
file, so a run that dies at region 27 resumes without redoing the first 26.

Outputs (schema-compatible with the previous notebooks, plus new columns):
    results/bench_parts/<grid>__<scene>.csv
    results/hpx_codec_bench.csv
    results/utm_codec_bench.csv
    results/common_reference_floor.csv        (with --common-reference)

Usage
-----
    python run_codec_bench.py --limit 2                 # smoke test
    python run_codec_bench.py --grid both               # the full sweep
    python run_codec_bench.py --grid both --common-reference
    python run_codec_bench.py --assemble                # rebuild the two CSVs
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

import build_tiles
import codec_lib
import paper_config as cfg
import region_sites

META_COLS = ["scene", "scene_class", "lat", "abs_lat", "hpx_regime"]


def environment() -> dict:
    """Record exactly which codec libraries produced these numbers.

    Codec behaviour is version-dependent in ways that change results silently:
    PyWavelets, for one, reconstructs an odd-sized periodized decomposition to a
    different shape depending on the release, which is why the 637x637 UTM crop
    can round-trip in one environment and fail in another while the 1024x1024
    HEALPix tile never does. The paper lists 'codec library and version not
    reported' as an open item; this closes it.
    """
    import platform
    import sys as _sys
    import hashlib
    from pathlib import Path as _P
    here = _P(__file__).resolve().parent
    digest = hashlib.sha256()
    for name in sorted(("common_reference.py", "run_codec_bench.py",
                        "codec_lib.py", "build_tiles.py", "paper_config.py")):
        f = here / name
        if f.exists():
            digest.update(f.read_bytes())
    env = {"python": _sys.version.split()[0], "platform": platform.platform(),
           # Fingerprint of the code that produced the results, so a results
           # folder can always be matched to the version that wrote it. Three
           # debugging rounds were spent on a stale copy; this is cheaper.
           "code_sha": digest.hexdigest()[:12],
           "written_utc": __import__("datetime").datetime.now(
               __import__("datetime").timezone.utc).isoformat(timespec="seconds")}
    for mod in ("numpy", "scipy", "pandas", "zarr", "pywt", "zfpy", "pysz",
                "imagecodecs", "PIL", "zstandard", "healpix_resample",
                "healpix_geo", "torch", "hcc"):
        try:
            m = __import__(mod)
            env[mod] = str(getattr(m, "__version__", "unknown"))
        except Exception:                                 # noqa: BLE001
            env[mod] = "absent"
    return env


def part_path(grid: str, scene: str) -> Path:
    return cfg.PARTS_DIR / f"{grid}__{scene}.csv"


def _tile_arrays(scene: str):
    from hcc.compression.zorder import expand_tile
    cell_ids, cell_data, utm_crop, meta, npz = build_tiles.load_tile(scene)
    hpx_tile = expand_tile(np.asarray(cell_data, dtype=np.float32)).astype(np.float32)
    return cell_ids, cell_data, hpx_tile, np.asarray(utm_crop, dtype=np.float32), meta, npz


def run_region(scene: str, site: dict, grids, codecs, common_ref=False,
               cr_mode="nearest", force=False):
    """Return (rows, floor_row). Rows are also cached per (grid, scene)."""
    from hcc.compression.zorder import flatten_tile

    todo = [g for g in grids if force or not part_path(g, scene).exists()]
    if not todo and not common_ref:
        return [], None
    if common_ref and "utm" not in todo and "utm" in grids:
        todo.append("utm")   # the UTM side must be re-scored on the same mask

    cell_ids, cell_data, hpx_tile, utm_tile, meta, npz = _tile_arrays(scene)
    base = {k: site.get(k, "") for k in META_COLS}
    base["scene"] = scene

    cr, floor_row = None, None
    if common_ref:
        import common_reference
        cr = common_reference.build_for_tile(npz, mode=cr_mode)
        floor = cr.rmse_on_mask(cr.predict(cell_data), utm_tile)
        std_on_mask = float(np.asarray(utm_tile, dtype=np.float64)[cr.mask].std())
        floor_row = dict(base, code_sha=environment()["code_sha"],
                         mode=cr.mode, mode_requested=cr_mode,
                         npt=cr.npt, scale_m=cr.scale_m,
                         raw_inside_fraction=cr.raw_inside_fraction,
                         support_median=cr.support_median,
                         mask_fraction=cr.mask_fraction,
                         representation_floor_rmse=floor,
                         utm_std=float(utm_tile.std()),
                         utm_std_on_mask=std_on_mask)
        print(f"    common reference: {cr.describe()}  floor RMSE {floor:.3e} "
              f"({100*floor/std_on_mask:.1f}% of the comparable area's std)")

    all_rows = []
    for grid in todo:
        tile = hpx_tile if grid == "hpx" else utm_tile
        raw_bytes = int(tile.size * 4)
        rows, stats = [], {}
        for run in codec_lib.iter_runs(tile, codecs, stats=stats):
            if np.shape(run.reconstruction) != np.shape(tile):
                print(f"    [skip] {run.codec} {run.config}: decoded shape "
                      f"{np.shape(run.reconstruction)} != tile "
                      f"{np.shape(tile)}")
                continue
            r = dict(base, grid=grid, codec=run.codec,
                     config_json=json.dumps(run.config),
                     total_bytes=int(run.total_bytes),
                     rmse=codec_lib.rmse(tile, run.reconstruction),
                     encode_seconds=run.encode_seconds,
                     decode_seconds=run.decode_seconds,
                     raw_bytes=raw_bytes)
            # Both grids are scored on the same pixels -- the intersection of
            # the two footprints -- so the comparison never depends on how the
            # two pixel lattices are oriented with respect to each other.
            if cr is None:
                r["rmse_common"] = r["rmse"] if grid == "utm" else np.nan
            elif grid == "utm":
                r["rmse_common"] = cr.rmse_on_mask(run.reconstruction, utm_tile)
            else:
                flat = flatten_tile(np.asarray(run.reconstruction, dtype=np.float32))
                r["rmse_common"] = cr.rmse_on_mask(cr.predict(flat), utm_tile)
            r["common_mask_fraction"] = (cr.mask_fraction if cr is not None
                                          else np.nan)
            rows.append(r)
        df = pd.DataFrame(rows)
        df.to_csv(part_path(grid, scene), index=False)
        attempted = sum(v["attempted"] for v in stats.values())
        (part_path(grid, scene).with_suffix(".skips.json")
         ).write_text(json.dumps(stats, indent=2), encoding="utf-8")
        print(f"    {grid}: {len(df)}/{attempted} sweep points -> "
              f"{part_path(grid, scene).name}")
        # A codec that produced nothing on one grid but works on the other is
        # not a detail: the oracle-best comparison would then be decided by
        # which codec was evaluable rather than by which one was better.
        empty = [c for c, v in stats.items() if v["ok"] == 0]
        if empty:
            print(f"    WARNING [{grid}] no usable point for: {', '.join(empty)}"
                  f" -- see {part_path(grid, scene).with_suffix('.skips.json').name}")
        all_rows.extend(rows)
    return all_rows, floor_row


def assemble():
    """Concatenate the part files into the two CSVs the analysis reads."""
    out = {}
    for grid, name in (("hpx", "hpx_codec_bench.csv"), ("utm", "utm_codec_bench.csv")):
        parts = sorted(cfg.PARTS_DIR.glob(f"{grid}__*.csv"))
        if not parts:
            print(f"no {grid} parts found in {cfg.PARTS_DIR}")
            continue
        df = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
        path = cfg.RESULTS_DIR / name
        # Never overwrite a previous sweep in place: these two CSVs are what
        # every published table was computed from, and a re-run with a
        # different region set would otherwise erase the only raw record of it.
        if path.exists():
            old = pd.read_csv(path)
            if set(old.scene.unique()) != set(df.scene.unique()):
                stamp = datetime.fromtimestamp(path.stat().st_mtime).strftime(
                    "%Y%m%dT%H%M%S")
                archive = cfg.RESULTS_DIR / "archive"
                archive.mkdir(exist_ok=True)
                dest = archive / f"{path.stem}.{stamp}.csv"
                if not dest.exists():
                    dest.write_bytes(path.read_bytes())
                    print(f"  previous sweep ({old.scene.nunique()} regions) "
                          f"archived to {dest.relative_to(cfg.RESULTS_DIR)}")
        df.to_csv(path, index=False)
        n = df.scene.nunique()
        print(f"{grid}: {len(df)} rows, {n} regions -> {path}")
        out[grid] = df
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grid", choices=["hpx", "utm", "both"], default="both")
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--classes", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--codecs", nargs="*", default=None,
                    help=f"default: {sorted(codec_lib.CODEC_RUNS)}")
    ap.add_argument("--common-reference", action="store_true",
                    help="also score the HEALPix side against the UTM observation")
    ap.add_argument("--cr-mode", choices=["nearest", "knn", "psf"],
                    default="nearest",
                    help="how the decoded HEALPix tile is read on the UTM "
                         "pixels: 'nearest' takes the cell containing each "
                         "pixel (exact, no smoothing of its own); 'knn'/'psf' "
                         "apply healpix-resample's Gaussian observation "
                         "operator, as a sensitivity check")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--assemble", action="store_true", help="only rebuild the CSVs")
    ap.add_argument("--legacy", action="store_true",
                    help="run on the seven original scenes instead of the 40 regions")
    args = ap.parse_args(argv)

    cfg.ensure_dirs()
    if args.assemble:
        assemble()
        return 0

    env = environment()
    (cfg.RESULTS_DIR / "bench_environment.json").write_text(
        json.dumps(env, indent=2), encoding="utf-8")
    print("[env] " + ", ".join(f"{k} {env[k]}" for k in
                                ("python", "numpy", "pywt", "zfpy", "pysz",
                                 "imagecodecs", "zstandard") if k in env))

    grids = ["hpx", "utm"] if args.grid == "both" else [args.grid]
    if args.legacy:
        sites = pd.DataFrame([dict(scene=s, scene_class=s, lat=np.nan,
                                   abs_lat=np.nan, hpx_regime="")
                              for s in cfg.LEGACY_SCENES])
        if args.regions:
            sites = sites[sites.scene.isin(set(args.regions))]
    else:
        sites = region_sites.select(region_sites.load_sites(),
                                    args.regions, args.classes, args.limit)
        print(region_sites.summary(sites))

    floors, failures = [], []
    for i, site in enumerate(sites.to_dict("records"), 1):
        scene = site["scene"]
        if not build_tiles.tile_path(scene).exists():
            failures.append((scene, "no tile -- run build_tiles.py first"))
            print(f"[{i:>2}/{len(sites)}] {scene}: NO TILE")
            continue
        print(f"[{i:>2}/{len(sites)}] {scene}")
        try:
            _, floor_row = run_region(scene, site, grids, args.codecs,
                                      common_ref=args.common_reference,
                                      cr_mode=args.cr_mode, force=args.force)
            if floor_row:
                floors.append(floor_row)
        except Exception as exc:                          # noqa: BLE001
            failures.append((scene, f"{type(exc).__name__}: {exc}"))
            print(f"    FAILED {type(exc).__name__}: {exc}")

    if floors:
        path = cfg.RESULTS_DIR / "common_reference_floor.csv"
        df = pd.DataFrame(floors)
        if path.exists():
            old = pd.read_csv(path)
            df = pd.concat([old[~old.scene.isin(df.scene)], df], ignore_index=True)
        df.to_csv(path, index=False)
        print(f"\nrepresentation floors -> {path}")

    assemble()
    if failures:
        print(f"\n{len(failures)} region(s) failed:")
        for scene, err in failures:
            print(f"  {scene}: {err}")
        return 1
    print("\nNext: python update_paper_numbers.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
