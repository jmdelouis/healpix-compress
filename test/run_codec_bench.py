"""Run the rate-distortion sweep on both grids, over the 40 shared regions.

Replaces the sweep cells of hpx_codec_bench.ipynb and utm_codec_bench.ipynb by
one resumable command. Every (region, grid) pair is written to its own part
file, so a run that dies at region 27 resumes without redoing the first 26.

Outputs (schema-compatible with the previous notebooks, plus new columns):
    results/bench_parts/<grid>[__<kernel>]__<scene>.csv
    results/hpx_codec_bench[__<kernel>].csv
    results/utm_codec_bench[__<kernel>].csv
    results/common_reference_floor[__<kernel>].csv   (with --common-reference)

The kernel suffix is empty for the default transfer, so the published sweep
keeps its filenames and a sweep on another transfer's tiles writes beside it
rather than over it.

Usage
-----
    python run_codec_bench.py --limit 2                 # smoke test
    python run_codec_bench.py --grid both               # the full sweep
    python run_codec_bench.py --grid both --common-reference
    python run_codec_bench.py --grid hpx --resampler psf   # deconvolved tiles
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
    # healpix-resample has no __version__ that survives an editable install,
    # so the module line above reads "unknown" and says nothing. Its checkout
    # is what the paper cites, so record the checkout's own description and
    # whether it sits on the tag the companion paper was submitted from.
    rv = cfg.resample_version()
    env["healpix_resample_git"] = rv["git"]
    env["healpix_resample_expected_tag"] = rv["expected_tag"]
    env["healpix_resample_on_expected_tag"] = rv["on_expected_tag"]
    return env


def kernel_tag(kernel: str) -> str:
    """Suffix for every artefact a non-default transfer produces.

    The default kernel keeps the unsuffixed names, so nothing already published
    moves; any other kernel writes beside it. Without this a psf sweep would
    overwrite the Clough--Tocher part files region by region and the two would
    be indistinguishable afterwards.
    """
    return "" if kernel == build_tiles.DEFAULT_KERNEL else f"__{kernel}"


def reader_tag(cr_mode: str) -> str:
    """Suffix for every HEALPix-side artefact a non-default reader produces.

    `rmse_common` is computed by reading the decoded tile on the native pixels,
    so it depends on the reader, and a sweep scored through the interpolating
    read must not overwrite the parts of the sweep scored through the nearest
    read: they are two measurements of the same tiles. The nearest read keeps
    the unsuffixed names, so nothing already published moves.
    """
    return "" if cr_mode in (None, "nearest") else f"__{cr_mode}"


def run_tag(kernel: str, cr_mode: str = "nearest") -> str:
    return kernel_tag(kernel) + reader_tag(cr_mode)


def part_path(grid: str, scene: str,
              kernel: str = build_tiles.DEFAULT_KERNEL,
              cr_mode: str = "nearest") -> Path:
    # The UTM side is scored on the intersection mask alone, which does not
    # depend on how the HEALPix tile is read (the mask is the footprint
    # overlap, eroded), so its part is the same whatever the reader and is
    # never re-swept for a change of reader.
    tag = kernel_tag(kernel) if grid == "utm" else run_tag(kernel, cr_mode)
    return cfg.PARTS_DIR / f"{grid}{tag}__{scene}.csv"


def _tile_arrays(scene: str, kernel: str = build_tiles.DEFAULT_KERNEL):
    from hcc.compression.zorder import expand_tile
    cell_ids, cell_data, utm_crop, meta, npz = build_tiles.load_tile(scene,
                                                                     kernel)
    hpx_tile = expand_tile(np.asarray(cell_data, dtype=np.float32)).astype(np.float32)
    return cell_ids, cell_data, hpx_tile, np.asarray(utm_crop, dtype=np.float32), meta, npz


def run_region(scene: str, site: dict, grids, codecs, common_ref=False,
               cr_mode="nearest", force=False,
               kernel=build_tiles.DEFAULT_KERNEL):
    """Return (rows, floor_row). Rows are also cached per (grid, kernel, scene)."""
    from hcc.compression.zorder import flatten_tile

    todo = [g for g in grids
            if force or not part_path(g, scene, kernel, cr_mode).exists()]
    if not todo and not common_ref:
        return [], None
    if common_ref and "utm" not in todo and "utm" in grids:
        # The UTM side must carry rmse_common on the same mask. It is
        # reader-independent, so an existing part that already has the column
        # is reused rather than re-encoded; one produced without
        # --common-reference is not.
        try:
            cols = pd.read_csv(part_path("utm", scene, kernel), nrows=1).columns
            have = "rmse_common" in cols
        except Exception:                                  # noqa: BLE001
            have = False
        if not have:
            todo.append("utm")

    cell_ids, cell_data, hpx_tile, utm_tile, meta, npz = _tile_arrays(scene, kernel)
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
        pp = part_path(grid, scene, kernel, cr_mode)
        df.to_csv(pp, index=False)
        attempted = sum(v["attempted"] for v in stats.values())
        pp.with_suffix(".skips.json").write_text(json.dumps(stats, indent=2),
                                                 encoding="utf-8")
        print(f"    {grid}: {len(df)}/{attempted} sweep points -> {pp.name}")
        # A codec that produced nothing on one grid but works on the other is
        # not a detail: the oracle-best comparison would then be decided by
        # which codec was evaluable rather than by which one was better.
        empty = [c for c, v in stats.items() if v["ok"] == 0]
        if empty:
            print(f"    WARNING [{grid}] no usable point for: {', '.join(empty)}"
                  f" -- see {pp.with_suffix('.skips.json').name}")
        all_rows.extend(rows)
    return all_rows, floor_row


def family_parts(grid: str, tag: str):
    """The part files of exactly one sweep family.

    A glob such as ``hpx__*.csv`` also matches ``hpx__psf__*.csv`` and
    ``hpx__bilinear__*.csv``, so assembling the default family used to sweep
    the other transfers and the other reader into the published CSV whenever
    their parts existed (revision 4). The remainder of the name after the
    family prefix must therefore be a scene name and nothing else: it may not
    begin with a kernel or reader token, and backups are not parts.
    """
    prefix = f"{grid}{tag}__"
    foreign = set(build_tiles.KERNELS) | {"bilinear", "knn", "psf"}
    out = []
    for p in cfg.PARTS_DIR.glob(f"{prefix}*.csv"):
        rest = p.name[len(prefix):-4]
        if ".pre_extension" in p.name or rest.split("__")[0] in foreign:
            continue
        out.append(p)
    return out


def assemble(kernel: str = build_tiles.DEFAULT_KERNEL, cr_mode: str = "nearest"):
    """Concatenate the part files into the two CSVs the analysis reads.

    A non-default reader assembles its own HEALPix-side CSV and reuses the
    UTM-side one, which does not depend on the reader (see `part_path`).
    """
    out = {}
    for grid, name in (("hpx", f"hpx_codec_bench{run_tag(kernel, cr_mode)}.csv"),
                       ("utm", f"utm_codec_bench{kernel_tag(kernel)}.csv")):
        tag = kernel_tag(kernel) if grid == "utm" else run_tag(kernel, cr_mode)
        parts = sorted(family_parts(grid, tag))
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
    ap.add_argument("--resampler", choices=list(build_tiles.KERNELS),
                    default=build_tiles.DEFAULT_KERNEL,
                    help="which transfer's HEALPix tiles to sweep. Anything "
                         "other than the default reads kernel-suffixed tiles "
                         "and writes kernel-suffixed parts and CSVs, so the "
                         "published sweep is never overwritten. The UTM side "
                         "is unresampled and identical whichever is chosen; it "
                         "is re-swept only so the pair stays self-contained")
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--classes", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--codecs", nargs="*", default=None,
                    help=f"default: {sorted(codec_lib.CODEC_RUNS)}")
    ap.add_argument("--common-reference", action="store_true",
                    help="also score the HEALPix side against the UTM observation")
    ap.add_argument("--cr-mode", choices=["nearest", "bilinear", "knn", "psf"],
                    default="nearest",
                    help="how the decoded HEALPix tile is read on the UTM "
                         "pixels: 'nearest' takes the cell containing each "
                         "pixel (exact, no smoothing of its own); 'bilinear' "
                         "interpolates within the carrier lattice at sub-cell "
                         "position (Section V-F: removes two thirds of the "
                         "nearest read's floor); 'knn'/'psf' apply "
                         "healpix-resample's observation operator. Any reader "
                         "but 'nearest' writes reader-suffixed parts and CSVs "
                         "on the HEALPix side, so the published sweep is never "
                         "overwritten, and reuses the UTM side, which does not "
                         "depend on the reader")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--env-only", action="store_true",
                    help="only rewrite results/bench_environment.json, so the "
                         "record can be refreshed after a library is pinned "
                         "without rerunning a sweep")
    ap.add_argument("--assemble", action="store_true", help="only rebuild the CSVs")
    ap.add_argument("--legacy", action="store_true",
                    help="run on the seven original scenes instead of the 40 regions")
    args = ap.parse_args(argv)

    cfg.ensure_dirs()
    if args.assemble:
        assemble(args.resampler, args.cr_mode)
        return 0

    env = environment()
    (cfg.RESULTS_DIR / "bench_environment.json").write_text(
        json.dumps(env, indent=2), encoding="utf-8")
    print("[env] " + ", ".join(f"{k} {env[k]}" for k in
                                ("python", "numpy", "pywt", "zfpy", "pysz",
                                 "imagecodecs", "zstandard") if k in env))
    print(f"[env] healpix-resample checkout: {env['healpix_resample_git']} "
          f"(expected tag {env['healpix_resample_expected_tag']})")
    if env["healpix_resample_on_expected_tag"] is False:
        print("[env] WARNING: the healpix-resample checkout is not on the tag "
              "the paper cites. Its PSF operator width has changed between "
              "versions before; check out the tag or update "
              "paper_config.HEALPIX_RESAMPLE_TAG on purpose.")
    if args.env_only:
        print(f"wrote {cfg.RESULTS_DIR / 'bench_environment.json'}")
        return 0

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
        if not build_tiles.tile_path(scene, args.resampler).exists():
            failures.append((scene, f"no {args.resampler} tile -- run "
                                    f"build_tiles.py --resampler "
                                    f"{args.resampler} first"))
            print(f"[{i:>2}/{len(sites)}] {scene}: NO TILE")
            continue
        print(f"[{i:>2}/{len(sites)}] {scene}")
        try:
            _, floor_row = run_region(scene, site, grids, args.codecs,
                                      common_ref=args.common_reference,
                                      cr_mode=args.cr_mode, force=args.force,
                                      kernel=args.resampler)
            if floor_row:
                floors.append(floor_row)
        except Exception as exc:                          # noqa: BLE001
            failures.append((scene, f"{type(exc).__name__}: {exc}"))
            print(f"    FAILED {type(exc).__name__}: {exc}")

    if floors:
        path = (cfg.RESULTS_DIR /
                f"common_reference_floor{run_tag(args.resampler, args.cr_mode)}.csv")
        df = pd.DataFrame(floors)
        if path.exists():
            old = pd.read_csv(path)
            df = pd.concat([old[~old.scene.isin(df.scene)], df], ignore_index=True)
        df.to_csv(path, index=False)
        print(f"\nrepresentation floors -> {path}")

    assemble(args.resampler, args.cr_mode)
    if failures:
        print(f"\n{len(failures)} region(s) failed:")
        for scene, err in failures:
            print(f"  {scene}: {err}")
        # A region the transfer is documented not to cover (Section V-L,
        # paper_config.KNOWN_REJECTIONS) has no tile by design; that is the
        # published state, not a failure. Anything else still fails the run.
        known = set(cfg.KNOWN_REJECTIONS.get(args.resampler, []))
        if any(scene not in known or not err.startswith("no ")
               for scene, err in failures):
            return 1
        print(f"  expected: documented coverage of the {args.resampler} "
              f"transfer (paper_config.KNOWN_REJECTIONS); not an error")
    if args.common_reference and args.cr_mode != "nearest":
        print(f"\nNext: python update_paper_numbers.py --metric common "
              f"--cr-mode {args.cr_mode}")
    else:
        print("\nNext: python update_paper_numbers.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
