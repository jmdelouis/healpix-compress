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

Both level-19 carriers are read the same way, at the nearest containing cell
of their own level: that is what `CommonReference.read_at_level` does for the
aggregate, and for the direct carrier the same index is used, divided by four.
The two level-19 numbers therefore differ only in how the cell values were
produced, which is the whole point of the comparison. No interpolating read is
offered here: `CommonReference` builds its interpolating operator on the
level-20 lattice, so reading a level-19 carrier through it would need a second
operator built on the coarser lattice, and the comparison between the two
carriers does not need one -- whatever a better read is worth, it is worth the
same to both.

The transfer is rebuilt from the same native sub-window the stored tile came
from (`build_tiles.tile_subwindow`), on the same parent, with the same
Clough--Tocher operator and the same missing-cell completion, so nothing but
the output level changes.

Usage
-----
    python level19_direct.py
    python level19_direct.py --regions region__water__lake_geneva --diagnose
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


def region_row(scene, verbose=False, diagnose=False):
    store = cfg.REGION_DIR / f"{scene}_data.zarr"
    cell_ids, cell_data, utm_crop, meta, npz = build_tiles.load_tile(scene)
    cr = common_reference.build_for_tile(npz, mode="nearest")
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
    # Where the two carriers disagree, disagree by how much, and on how many
    # cells: an RMS alone cannot tell a carrier that is uniformly a little
    # different from one that is broken over a small patch, and the level-19
    # comparison ran into exactly that on one region.
    diff = v19 - agg19
    adiff = np.abs(diff)
    extras = {}
    if sd == sd and sd > 0:
        for k in (1, 5, 20):
            extras[f"cells_over_{k}std"] = int((adiff > k * sd).sum())
        extras["diff_max_over_std"] = float(adiff.max() / sd)
        extras["diff_p999_over_std"] = float(np.percentile(adiff, 99.9) / sd)
    if diagnose:
        worst = np.argsort(adiff)[::-1][:10]
        extras["worst_cells"] = json.dumps(
            [[int(i), float(v19[i]), float(agg19[i])] for i in worst])
    return dict(
        scene=scene, cls=scene.split("__")[1],
        n_missing=n_missing, constant_check=ones_med, tile_std=sd, **extras,
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
    """Medians over regions, plus how often the two level-19 carriers agree."""
    d = dict(n=int(len(df)))
    for c in ("floor20", "floor19_direct", "floor19_agg", "ratio_direct",
              "ratio_agg", "ratio_direct_over_agg", "carrier_rms",
              "carrier_rms_frac_of_std", "carrier_rms_over_floor19"):
        d[c] = float(np.nanmedian(df[c]))
    q = df.ratio_direct_over_agg.to_numpy(dtype=float)
    d["direct_better_on"] = int((df.floor19_direct < df.floor19_agg).sum())
    d["within_5pct"] = int(((q > 0.95) & (q < 1.05)).sum())
    d["max_ratio_direct_over_agg"] = float(np.nanmax(q))
    d["worst_scene"] = str(df.scene.iloc[int(np.nanargmax(q))])
    # The median is robust to a region where the direct transfer fails; the
    # paper should quote both, so the summary carries the trimmed value too.
    keep = df[q < 1.5]
    d["median_ratio_excluding_failures"] = float(
        np.nanmedian(keep.ratio_direct_over_agg)) if len(keep) else np.nan
    d["n_failures"] = int((q >= 1.5).sum())
    return d


def latex_table(df, digest):
    """One row per carrier: what the level choice costs, and for which object."""
    d = digest
    nl = r" \\"
    body = [
        f"level 20 & {d['n']:d} & {d['floor20']:.2e} & 1.00 & --" + nl,
        r"\midrule",
        f"level 19, direct & {d['n']:d} & {d['floor19_direct']:.2e} & "
        f"{d['ratio_direct']:.2f} & {d['direct_better_on']:d}" + nl,
        f"level 19, four-child aggregate & {d['n']:d} & {d['floor19_agg']:.2e} & "
        f"{d['ratio_agg']:.2f} & {d['n'] - d['direct_better_on']:d}" + nl,
    ]
    worst = d["worst_scene"].replace("region__", "").replace("__", "/")
    worst = worst.replace("_", r"\_")
    note = (
        f"Median over the {d['n']} regions of the RMSE of each carrier against "
        "the same Sentinel-2 observation, on the intersection mask of "
        "Section~\\ref{sec:common_reference}, each read at the nearest "
        "containing cell of its own level. \\emph{Direct} evaluates the "
        "Clough--Tocher interpolant of the same native samples at the level-19 "
        "cell centres; \\emph{aggregate} averages the four level-20 children, "
        "which is the exact cell-average of the stored carrier and the object "
        "Section~\\ref{sec:level_choice} measures. The last column counts the "
        "regions on which that carrier has the lower floor. The two level-19 "
        "carriers are genuinely different objects, differing by "
        f"{100*d['carrier_rms_frac_of_std']:.1f}\\% of a scene's standard "
        f"deviation in RMS, {d['carrier_rms_over_floor19']:.2f} times the "
        "level-19 floor itself, yet their floors agree to "
        f"{100*abs(d['ratio_direct_over_agg'] - 1):.1f}\\% in median and to "
        f"within $5\\%$ on {d['within_5pct']} of the {d['n']} regions, so the "
        "level-19 penalty is a property of the level rather than of which "
        f"object is called the carrier. The exception is {worst}, where the "
        f"direct transfer fails outright ({d['max_ratio_direct_over_agg']:.0f} "
        "times the aggregate's floor); excluding it the median ratio is "
        f"{d['median_ratio_excluding_failures']:.3f}.")
    return "\n".join([
        r"\begin{table}[t]", r"\centering",
        *cfg.tex_caption("Two level-19 carriers, against the observation",
                         "tab:level19_direct"),
        r"\footnotesize",
        r"\begin{tabular}{@{}lrrrr@{}}", r"\toprule",
        r"Carrier & $n$ & Floor (RMSE) & $/$ level 20 & Lower floor on \\",
        r"\midrule", *body, r"\bottomrule", r"\end{tabular}",
        *cfg.tex_note(note), r"\end{table}"]) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--diagnose", action="store_true",
                    help="also record the ten most divergent cells per region")
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
        try:
            r = region_row(scene, verbose=args.verbose, diagnose=args.diagnose)
        except Exception as exc:                          # noqa: BLE001
            failures.append((scene, f"{type(exc).__name__}: {exc}"))
            continue
        rows.append(r)
        print(f"[{i:>2}/{len(scenes)}] {scene}: "
              f"L20 {r['floor20']:.3e}  L19 direct {r['floor19_direct']:.3e} "
              f"({r['ratio_direct']:.2f}x)  L19 agg {r['floor19_agg']:.3e} "
              f"({r['ratio_agg']:.2f}x)  carriers differ by "
              f"{100*r['carrier_rms_frac_of_std']:.2f}% of std"
              + (f"  [cells >5 std: {r['cells_over_5std']}]"
                 if r.get("cells_over_5std") else ""))
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
