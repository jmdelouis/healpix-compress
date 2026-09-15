"""Complete every cached sweep with the points codec_lib now defines but the
part files lack, then reassemble the CSVs the analysis reads.

Why
---
A rate--distortion target that a codec's coarsest setting already satisfies
is returned by the protocol as a bound, not an estimate
(paper_stats.FEASIBLE_CENSORING). At the loosest target of the paper, RMSE
1e-2, the original ladders left the oracle codec censored on 16 of the 40
regions (censoring_sensitivity.py), so the sweeps were extended at their
coarse end (codec_lib.*_SWEEP). Re-running run_codec_bench.py --force would
recompute every point of every region for every transfer and every reader,
hours of work to add a dozen cheap points; this script adds only the missing
ones, region by region, to the same part files, with the same row layout
and the same common-reference scoring, and leaves everything already measured
untouched.

It covers every sweep family in results/bench_parts: the published
Clough--Tocher sweep and its interpolating-read re-score, the PSF-aware and
inverse-distance transfers, and the UTM side. Extending one family and not
the others would move the censoring of one column of a comparison and not
the other, which is worse than the censoring itself.

Usage
-----
    python extend_sweep.py --dry-run          # count what is missing, write nothing
    python extend_sweep.py                    # add the points, reassemble every CSV
    python extend_sweep.py --regions region__water__lake_tahoe ...   # a subset

Before assembling, the previous CSVs are copied to
results/archive/pre_extension_<timestamp>/ so that the tables of the
submitted revision can always be regenerated.
"""

from __future__ import annotations

import argparse
import json
import shutil
from datetime import datetime

import numpy as np
import pandas as pd

import build_tiles
import codec_lib
import paper_config as cfg
import run_codec_bench as rcb

READERS = ("bilinear", "knn", "psf")


def parse_part(path):
    """(grid, kernel, cr_mode, scene) from a part file name, or None."""
    tokens = path.stem.split("__")
    if "region" not in tokens:
        return None                       # legacy scenes: not part of the paper
    i = tokens.index("region")
    grid, tags, scene = tokens[0], tokens[1:i], "__".join(tokens[i:])
    kernel, cr_mode = build_tiles.DEFAULT_KERNEL, "nearest"
    for t in tags:
        if t in build_tiles.KERNELS and kernel == build_tiles.DEFAULT_KERNEL:
            kernel = t            # run_tag puts the kernel first
        elif t in READERS:
            cr_mode = t
        else:
            return None
    if grid == "utm":
        cr_mode = "nearest"
    return grid, kernel, cr_mode, scene


def _same(conf, stored):
    """Does a stored config_json describe the sweep point `conf`?"""
    try:
        d = json.loads(stored)
    except Exception:                                      # noqa: BLE001
        return False
    for k, v in conf.items():
        if k not in d:
            return False
        a, b = d[k], v
        if isinstance(v, (int, float)) and isinstance(a, (int, float)):
            if not np.isclose(float(a), float(b), rtol=1e-9, atol=0):
                return False
        elif a != b:
            return False
    return True


def missing_points(df, codecs):
    out = []
    for codec in codecs:
        fn, sweep = codec_lib.CODEC_RUNS[codec]
        stored = df[df.codec == codec].config_json.tolist()
        for conf in sweep:
            if not any(_same(conf, s) for s in stored):
                out.append((codec, fn, conf))
    return out


def extend_part(path, codecs, dry=False):
    parsed = parse_part(path)
    if parsed is None:
        return 0, None
    grid, kernel, cr_mode, scene = parsed
    df = pd.read_csv(path)
    todo = missing_points(df, codecs)
    if not todo:
        return 0, (kernel, cr_mode)
    if dry:
        return len(todo), (kernel, cr_mode)

    from hcc.compression.zorder import flatten_tile
    cell_ids, cell_data, hpx_tile, utm_tile, meta, npz = rcb._tile_arrays(scene, kernel)
    tile = hpx_tile if grid == "hpx" else utm_tile
    raw_bytes = int(tile.size * 4)
    base = {k: (df[k].iloc[0] if k in df.columns else "") for k in rcb.META_COLS}
    base["scene"] = scene

    # Common-reference scoring exactly as run_region does it, on the same mask.
    cr = None
    if "rmse_common" in df.columns and df.rmse_common.notna().any():
        import common_reference
        cr = common_reference.build_for_tile(npz, mode=cr_mode)

    rows = []
    for codec, fn, conf in todo:
        try:
            run = codec_lib.match_shape(tile, fn(tile, **conf))
        except Exception as exc:                          # noqa: BLE001
            print(f"    [skip] {codec} {conf}: {type(exc).__name__}: {exc}")
            continue
        if np.shape(run.reconstruction) != np.shape(tile):
            print(f"    [skip] {codec} {conf}: decoded shape mismatch")
            continue
        r = dict(base, grid=grid, codec=run.codec,
                 config_json=json.dumps(run.config),
                 total_bytes=int(run.total_bytes),
                 rmse=codec_lib.rmse(tile, run.reconstruction),
                 encode_seconds=run.encode_seconds,
                 decode_seconds=run.decode_seconds,
                 raw_bytes=raw_bytes)
        if cr is None:
            r["rmse_common"] = r["rmse"] if grid == "utm" else np.nan
            r["common_mask_fraction"] = np.nan
        elif grid == "utm":
            r["rmse_common"] = cr.rmse_on_mask(run.reconstruction, utm_tile)
            r["common_mask_fraction"] = cr.mask_fraction
        else:
            flat = flatten_tile(np.asarray(run.reconstruction, dtype=np.float32))
            r["rmse_common"] = cr.rmse_on_mask(cr.predict(flat), utm_tile)
            r["common_mask_fraction"] = cr.mask_fraction
        rows.append(r)
    if not rows:
        return 0, (kernel, cr_mode)
    backup = path.with_suffix(".pre_extension.csv")
    if not backup.exists():
        shutil.copyfile(path, backup)
    new = pd.concat([df, pd.DataFrame(rows)], ignore_index=True)
    new = new.reindex(columns=list(df.columns) + [c for c in new.columns if c not in df.columns])
    new.to_csv(path, index=False)
    return len(rows), (kernel, cr_mode)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--codecs", nargs="*", default=None,
                    help=f"default: {cfg.COMPRESSOR_CODECS}")
    ap.add_argument("--no-assemble", action="store_true")
    args = ap.parse_args(argv)
    codecs = args.codecs or cfg.COMPRESSOR_CODECS

    parts = sorted(cfg.PARTS_DIR.glob("*.csv"))
    parts = [p for p in parts if not p.name.endswith(".pre_extension.csv")]
    if args.regions:
        keep = set(args.regions)
        parts = [p for p in parts
                 if parse_part(p) and parse_part(p)[3] in keep]
    if not parts:
        print(f"no part files under {cfg.PARTS_DIR}")
        return 1

    families, added, touched = set(), 0, 0
    for i, p in enumerate(parts, 1):
        n, fam = extend_part(p, codecs, dry=args.dry_run)
        if fam is not None:
            families.add(fam)
        if n:
            touched += 1
            added += n
            print(f"[{i:>3}/{len(parts)}] {p.name}: "
                  f"{'would add' if args.dry_run else 'added'} {n} point(s)")
    print(f"\n{added} point(s) {'missing' if args.dry_run else 'added'} in "
          f"{touched} of {len(parts)} part file(s); sweep families: "
          + ", ".join(f"{k}/{m}" for k, m in sorted(families)))
    if args.dry_run or args.no_assemble or not added:
        return 0

    stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    archive = cfg.RESULTS_DIR / "archive" / f"pre_extension_{stamp}"
    archive.mkdir(parents=True, exist_ok=True)
    for f in cfg.RESULTS_DIR.glob("*_codec_bench*.csv"):
        shutil.copyfile(f, archive / f.name)
    print(f"previous CSVs copied to {archive}")
    for kernel, cr_mode in sorted(families):
        print(f"\nassembling kernel={kernel} reader={cr_mode}")
        rcb.assemble(kernel, cr_mode)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
