"""Build the level-19 carrier the same transfer would produce, and price it.

Why this exists
---------------
Section V-I answered "is level 19 enough?" from the stored level-20 tiles, by
aggregating them in groups of four. In NESTED ordering those four are the
children of one level-19 cell and HEALPix cells are equal-area, so the
aggregate is the exact cell-average of the carrier over each level-19 cell.
It is not, however, the carrier that the same transfer would have built at
level 19: the carrier stores point evaluations of the Clough--Tocher
interpolant at cell centres, and for a general interpolant

    f(x_parent) != (1/4) sum_i f(x_child,i)

The two objects differ by the curvature of the interpolant over a cell, which
is not obviously small compared with what the level choice is being asked to
decide. This script builds the other object -- the same interpolant, over the
same native samples, evaluated at the level-19 cell centres -- so that the
paper can say which quantity it is reporting and how far apart the two are.

What is measured, per region, on the same intersection mask as everywhere
else in the paper:

    floor20        RMSE of the level-20 carrier against the observation
    floor19_direct RMSE of the directly evaluated level-19 carrier
    floor19_agg    RMSE of the four-child (quadtree Haar) aggregate
    carrier_rms    RMS difference between the two level-19 carriers, and the
                   same in units of each one's own floor

Both reads are the nearest containing cell at the carrier's own level, which
is what `CommonReference.read_at_level` does for the aggregate; for the direct
carrier the same index is used, divided by four, so the two level-19 numbers
differ only in how the cell values were produced.

The transfer is rebuilt from the same native sub-window the stored tile came
from (`build_tiles.tile_subwindow`), on the same parent, with the same
Clough--Tocher operator and the same missing-cell completion, so nothing but
the output level changes.

Usage
-----
    python level19_direct.py                    # 40 regions, nearest read
    python level19_direct.py --reader bilinear  # and through the better read
    python level19_direct.py --regions region__urban__paris
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

import build_tiles
import common_reference
import paper_config as cfg
import region_sites


COARSE_DELTA = 1                      # one level up from the carrier
LEVEL19 = cfg.LEVEL - COARSE_DELTA
N19 = cfg.N_CHILDREN // (4 ** COARSE_DELTA)


def rms(a):
    a = np.asarray(a, dtype=np.float64).ravel()
    return float(np.sqrt(np.mean(a * a)))


def build_level19(scene, store, verbose=False):
    """The Clough--Tocher interpolant of the same samples, at level-19 centres."""
    import healpix_resample as hr

    sw = build_tiles.tile_subwindow(scene, store)
    r0, r1, c0, c1 = sw["bounds"]
    sub_val = np.asarray(sw["img"][r0:r1, c0:c1], dtype=np.float64).ravel()
    good = np.isfinite(sub_val)
    parent = int(sw["cand"]["parent"])

    ids19 = parent * N19 + np.arange(N19, dtype=np.int64)
    op = hr.CloughTocherResampler(lon_deg=sw["sub_lon"][good],
                                  lat_deg=sw["sub_lat"][good],
                                  level=LEVEL19, out_cell_ids=ids19,
                                  verbose=verbose)
    # The constant check build_tiles runs on every tile, for the same reason:
    # a transfer that does not preserve a constant makes every error below an
    # offset measurement.
    ones = op.resample(np.ones(int(good.sum()), dtype=np.float64))
    ones_med = float(np.nanmedian(np.asarray(ones.cell_data, dtype=np.float64)))
    if not np.isfinite(ones_med) or abs(ones_med - 1.0) > 0.05:
        raise RuntimeError(f"{scene}: the level-19 transfer maps a constant to "
                           f"{ones_med} instead of 1.0")

    res = op.resample(sub_val[good])
    out_ids = np.asarray(res.cell_ids, dtype=np.int64).ravel()
    out_val = np.asarray(res.cell_data, dtype=np.float64).ravel()
    filled, n_missing = build_tiles.complete_block(out_ids, out_val, ids19)
    return np.asarray(filled, dtype=np.float64), int(n_missing), ones_med


def region_row(scene, reader="nearest", verbose=False):
    store = cfg.REGION_DIR / f"{scene}_data.zarr"
    cell_ids, cell_data, utm_crop, meta, npz = build_tiles.load_tile(scene)
    cr = common_reference.build_for_tile(npz, mode=reader)
    obs = np.asarray(utm_crop, dtype=np.float64)
    flat20 = np.asarray(cell_data, dtype=np.float64).ravel()
    sd = float(flat20.std()) or np.nan

    v19, n_missing, ones_med = build_level19(scene, store, verbose=verbose)
    agg19 = flat20.reshape(-1, 4 ** COARSE_DELTA).mean(axis=1)

    idx19 = np.clip(cr._safe_index // (4 ** COARSE_DELTA), 0, N19 - 1)
    read19_direct = v19[idx19].reshape(cr.shape)
    read19_agg = cr.read_at_level(flat20, COARSE_DELTA).reshape(cr.shape)
    read20 = cr.read_at_level(flat20, 0).reshape(cr.shape)

    floor20 = cr.rmse_on_mask(read20, obs)
    floor19_direct = cr.rmse_on_mask(read19_direct, obs)
    floor19_agg = cr.rmse_on_mask(read19_agg, obs)
    return dict(
        scene=scene, cls=scene.split("__")[1], reader=reader,
        n_missing=n_missing, constant_check=ones_med, tile_std=sd,
        floor20=floor20, floor19_direct=floor19_direct, floor19_agg=floor19_agg,
        ratio_direct=floor19_direct / floor20 if floor20 else np.nan,
        ratio_agg=floor19_agg / floor20 if floor20 else np.nan,
        ratio_direct_over_agg=(floor19_direct / floor19_agg
                               if floor19_agg else np.nan),
        carrier_rms=rms(v19 - agg19),
        carrier_rms_frac_of_std=rms(v19 - agg19) / sd if sd == sd else np.nan,
        carrier_rms_over_floor19=(rms(v19 - agg19) / floor19_agg
                                  if floor19_agg else np.nan),
        mask_fraction=cr.mask_fraction)


def summarize(df):
    out = {}
    for reader, g in df.groupby("reader"):
        d = dict(n=int(len(g)))
        for c in ("floor20", "floor19_direct", "floor19_agg", "ratio_direct",
                  "ratio_agg", "ratio_direct_over_agg", "carrier_rms",
                  "carrier_rms_frac_of_std", "carrier_rms_over_floor19"):
            d[c] = float(np.nanmedian(g[c]))
        d["direct_better_on"] = int((g.floor19_direct < g.floor19_agg).sum())
        out[reader] = d
    return out


def latex_table(df, digest):
    """One row per carrier, per reader: what the level choice actually costs."""
    readers = [r for r in ("nearest", "bilinear") if r in digest]
    body = []
    for r in readers:
        d = digest[r]
        label = "nearest cell" if r == "nearest" else "interpolating"
        body.append(
            f"{label} & {d['n']:d} & {d['floor20']:.2e} & "
            f"{d['floor19_direct']:.2e} & {d['ratio_direct']:.2f} & "
            f"{d['floor19_agg']:.2e} & {d['ratio_agg']:.2f} & "
            f"{d['carrier_rms_over_floor19']:.2f} \\\\")
    note = (
        "Median over the regions of the RMSE of each carrier against the same "
        "Sentinel-2 observation, on the intersection mask of "
        "Section~\\ref{sec:common_reference}, read at the nearest containing "
        "cell of that carrier's own level. \\emph{Direct} evaluates the "
        "Clough--Tocher interpolant of the same native samples at the level-19 "
        "cell centres; \\emph{aggregate} averages the four level-20 children, "
        "which is the exact cell-average of the stored carrier and the object "
        "Section~\\ref{sec:level_choice} measures. The last column is the RMS "
        "difference between the two level-19 carriers in units of the "
        "aggregate's own floor: it says how much of the level-19 error budget "
        "is decided by which of the two objects is called the level-19 "
        "carrier. Ratios are to the level-20 floor of the same read.")
    return "\n".join([
        r"\begin{table}[t]", r"\centering",
        *cfg.tex_caption("Two level-19 carriers, against the observation",
                         "tab:level19_direct"),
        r"\footnotesize", r"\setlength{\tabcolsep}{4pt}",
        r"\begin{tabular}{@{}lrrrrrrr@{}}", r"\toprule",
        r" & & Level 20 & \multicolumn{2}{c}{Direct} & \multicolumn{2}{c}{Aggregate} & Carrier \\",
        r"\cmidrule(lr){4-5}\cmidrule(lr){6-7}",
        r"Read & $n$ & floor & floor & $/$L20 & floor & $/$L20 & diff. \\",
        r"\midrule", *body, r"\bottomrule", r"\end{tabular}",
        *cfg.tex_note(note), r"\end{table}"]) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--reader", nargs="*", default=["nearest", "bilinear"],
                    help="which read(s) to score through; default both")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--no-tex", action="store_true")
    args = ap.parse_args(argv)

    sites = region_sites.select(region_sites.load_sites(), args.regions,
                                None, args.limit)
    scenes = list(sites.scene)
    rows, failures = [], []
    for i, scene in enumerate(scenes, 1):
        if not build_tiles.tile_path(scene).exists():
            failures.append((scene, "no level-20 tile; run build_tiles.py"))
            continue
        for reader in args.reader:
            try:
                r = region_row(scene, reader=reader, verbose=args.verbose)
            except Exception as exc:                      # noqa: BLE001
                failures.append((scene, f"{type(exc).__name__}: {exc}"))
                break
            rows.append(r)
            print(f"[{i:>2}/{len(scenes)}] {scene} ({reader}): "
                  f"L20 {r['floor20']:.3e}  L19 direct {r['floor19_direct']:.3e} "
                  f"({r['ratio_direct']:.2f}x)  L19 agg {r['floor19_agg']:.3e} "
                  f"({r['ratio_agg']:.2f}x)  carriers differ by "
                  f"{100*r['carrier_rms_frac_of_std']:.2f}% of std")
    if not rows:
        for s, why in failures:
            print(f"  {s}: {why}")
        raise SystemExit("no region could be evaluated")

    df = pd.DataFrame(rows)
    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(cfg.RESULTS_DIR / "level19_direct.csv", index=False)
    digest = summarize(df)
    (cfg.RESULTS_DIR / "level19_direct.json").write_text(
        json.dumps(digest, indent=1))
    if not args.no_tex:
        cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
        (cfg.TEX_TABLES / "table_level19_direct.tex").write_text(
            latex_table(df, digest), encoding="utf-8")
        print(f"wrote {cfg.TEX_TABLES / 'table_level19_direct.tex'}")

    print("\n" + json.dumps(digest, indent=1))
    if failures:
        print(f"\n{len(failures)} region(s) could not be evaluated:")
        for s, why in failures:
            print(f"  {s}: {why}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
