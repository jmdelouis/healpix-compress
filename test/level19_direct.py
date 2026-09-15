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

Diagnosing a region where the direct carrier fails
--------------------------------------------------
`--diagnose` adds the control that separates the two things that can go
wrong. The direct level-19 carrier is built by re-running the transfer from
the native sub-window; if that rebuild is faithful, the stored level-20 tile
must come back out of it when the same operator is asked for level-20 cells.
So `--diagnose` rebuilds at level 20 as well and reports the RMS difference
against the tile on disk. Two outcomes, two different causes:

  rebuild20 matches the stored tile   the sub-window and the samples are the
                                      ones the tile came from, and whatever
                                      fails is in evaluating the interpolant
                                      at the coarser centres;
  rebuild20 does not match            the rebuild itself is not reproducing
                                      the tile for this region, and the
                                      level-19 number measures that, not the
                                      level.

`--dump` additionally writes results/level19_diag__<scene>.npz with both
level-19 carriers and the level-19 cell centres, so the difference can be
mapped without re-running anything.
"""

from __future__ import annotations

import argparse
import json
import shutil

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


def build_at_level(scene, sw, delta, verbose=False):
    """The Clough--Tocher interpolant of the same samples, at one level's centres.

    `delta` counts levels up from the carrier: 0 rebuilds the stored tile,
    1 builds the level-19 object this section is about. Everything else --
    sub-window, no-data handling, operator, missing-cell completion -- is what
    `build_tiles.build_tile` does, so that the only difference between the two
    calls is the output level.
    """
    import healpix_resample as hr

    level = cfg.LEVEL - delta
    n_cells = cfg.N_CHILDREN // (4 ** delta)
    r0, r1, c0, c1 = sw["bounds"]
    sub_val = np.asarray(sw["img"][r0:r1, c0:c1], dtype=np.float64).ravel()
    good = np.isfinite(sub_val)
    parent = int(sw["cand"]["parent"])

    ids = parent * n_cells + np.arange(n_cells, dtype=np.int64)
    op = hr.CloughTocherResampler(lon_deg=sw["sub_lon"][good],
                                  lat_deg=sw["sub_lat"][good],
                                  level=level, out_cell_ids=ids,
                                  verbose=verbose)
    # The constant check build_tiles runs on every tile, for the same reason:
    # a transfer that does not preserve a constant makes every error below an
    # offset measurement.
    ones = op.resample(np.ones(int(good.sum()), dtype=np.float64))
    ones_med = float(np.nanmedian(np.asarray(ones.cell_data, dtype=np.float64)))
    if not np.isfinite(ones_med) or abs(ones_med - 1.0) > 0.05:
        raise RuntimeError(f"{scene}: the level-{level} transfer maps a constant "
                           f"to {ones_med} instead of 1.0")

    res = op.resample(sub_val[good])
    out_ids = np.asarray(res.cell_ids, dtype=np.int64).ravel()
    out_val = np.asarray(res.cell_data, dtype=np.float64).ravel()
    n_out = len(np.unique(out_ids))
    filled, n_missing = build_tiles.complete_block(out_ids, out_val, ids)
    return dict(values=np.asarray(filled, dtype=np.float64), ids=ids,
                n_missing=int(n_missing), ones_med=ones_med,
                n_samples=int(good.sum()), n_dropped=int((~good).sum()),
                n_returned=int(n_out))


def region_row(scene, verbose=False, diagnose=False, dump=False, control=True):
    store = cfg.REGION_DIR / f"{scene}_data.zarr"
    cell_ids, cell_data, utm_crop, meta, npz = build_tiles.load_tile(scene)
    cr = common_reference.build_for_tile(npz, mode="nearest")
    obs = np.asarray(utm_crop, dtype=np.float64)
    flat20 = np.asarray(cell_data, dtype=np.float64).ravel()
    sd = float(flat20.std()) or np.nan

    sw = build_tiles.tile_subwindow(scene, store)
    b = build_at_level(scene, sw, COARSE_DELTA, verbose=verbose)
    v19, n_missing, ones_med = b["values"], b["n_missing"], b["ones_med"]
    n_used, n_dropped = b["n_samples"], b["n_dropped"]
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
    extras = dict(bias=float(diff.mean()),
                  bias_over_rms=float(abs(diff.mean()) / rms(diff))
                  if rms(diff) else np.nan)
    if sd == sd and sd > 0:
        for k in (1, 5, 20):
            extras[f"cells_over_{k}std"] = int((adiff > k * sd).sum())
        extras["diff_max_over_std"] = float(adiff.max() / sd)
        extras["diff_p999_over_std"] = float(np.percentile(adiff, 99.9) / sd)
    if diagnose:
        worst = np.argsort(adiff)[::-1][:10]
        extras["worst_cells"] = json.dumps(
            [[int(i), float(v19[i]), float(agg19[i])] for i in worst])
        # What kind of object is the direct carrier? A permutation of the right
        # values, a shifted copy, an overshoot, or a different field entirely:
        # the first two moments and the correlation separate those cases, and
        # none of them can be read off an RMS.
        extras.update(
            direct_mean=float(v19.mean()), direct_std=float(v19.std()),
            direct_min=float(v19.min()), direct_max=float(v19.max()),
            agg_mean=float(agg19.mean()), agg_std=float(agg19.std()),
            corr_direct_agg=float(np.corrcoef(v19, agg19)[0, 1]),
            n_returned=b["n_returned"], n_cells=int(v19.size))
        # The control: the same rebuild, asked for the level the tile is
        # stored at. If it does not return the stored tile, the level-19
        # number is measuring the rebuild and not the level. It costs a second
        # triangulation of the same ~700k native samples, which doubles the
        # run; --no-control skips it when only the moments are wanted.
        r20 = None
        if control:
            c = build_at_level(scene, sw, 0, verbose=verbose)
            r20 = c["values"]
            extras.update(
                rebuild20_rms=rms(r20 - flat20),
                rebuild20_rms_frac_of_std=(rms(r20 - flat20) / sd
                                           if sd == sd else np.nan),
                rebuild20_corr=float(np.corrcoef(r20, flat20)[0, 1]),
                rebuild20_n_missing=c["n_missing"])
        if dump:
            lon = lat = None
            try:
                import healpix_geo
                lon, lat = healpix_geo.nested.healpix_to_lonlat(
                    b["ids"], LEVEL19, ellipsoid="WGS84")
            except Exception:                              # noqa: BLE001
                pass
            out = dict(direct=v19.astype(np.float32),
                       aggregate=agg19.astype(np.float32),
                       stored20=flat20.astype(np.float32))
            if r20 is not None:
                out["rebuild20"] = r20.astype(np.float32)
            if lon is not None:
                out["lon"] = np.asarray(lon, dtype=np.float32)
                out["lat"] = np.asarray(lat, dtype=np.float32)
            p = cfg.RESULTS_DIR / f"level19_diag__{scene}.npz"
            cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(p, **out)
            print(f"    dumped {p}")
    return dict(
        scene=scene, cls=scene.split("__")[1],
        n_missing=n_missing, constant_check=ones_med, tile_std=sd,
        n_samples=n_used, n_dropped=n_dropped,
        dropped_fraction=n_dropped / max(1, n_used + n_dropped), **extras,
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


def excluded_scenes():
    """Regions whose window is not the product their stored tile came from.

    This comparison rebuilds the transfer from the native window, so it is
    only meaningful where the window still holds the acquisition the tile was
    built from. `window_radiometry.py` tests that for all 40 regions, and one
    fails: over Lake Geneva the window is related to the stored tile by a gain
    and an offset of its own, so a rebuild from it reproduces neither the
    level-19 carrier nor the level-20 one. Excluding it here removes an input
    mismatch from a measurement about levels; it is not trimming an
    inconvenient result, and the region is reported in its own right.
    """
    return (cfg.window_mismatched_scenes(),
            cfg.RESULTS_DIR / "window_radiometry.csv")


def summarize(df, dropped=None):
    """Medians over the regions the window audit keeps, plus what it dropped.

    `dropped` is the rows for the excluded regions. They are reported, not
    discarded: a region whose window is a different product is a fact about
    the archive that the paper states, and the reader can see what it would
    have contributed had it been kept.
    """
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
    d["n_excluded"] = 0 if dropped is None else int(len(dropped))
    d["excluded_scenes"] = [] if dropped is None else list(dropped.scene)
    if dropped is not None and len(dropped):
        d["excluded"] = [
            dict(scene=str(r["scene"]),
                 ratio_direct_over_agg=float(r["ratio_direct_over_agg"]),
                 ratio_direct=float(r["ratio_direct"]))
            for _, r in dropped.iterrows()]
        d["median_ratio_including_excluded"] = float(
            np.nanmedian(pd.concat([df, dropped]).ratio_direct_over_agg))
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
        f"object is called the carrier. The widest departure is {worst}, at "
        f"{d['max_ratio_direct_over_agg']:.3f}.")
    if d.get("n_excluded"):
        names = ", ".join(
            s.replace("region__", "").replace("__", "/").replace("_", r"\_")
            for s in d["excluded_scenes"])
        n_ex = d["n_excluded"]
        one = n_ex == 1
        note += (
            f" {'One region is' if one else f'{n_ex} regions are'} excluded and "
            f"reported separately ({names}): this comparison rebuilds the "
            "transfer from the native window, and the audit of "
            "Appendix~\\ref{app:window_radiometry} finds "
            f"{'its' if one else 'their'} window is not the product "
            f"{'its' if one else 'their'} stored tile was built from, so a "
            f"rebuild reproduces neither level. Kept, the median ratio would "
            f"be {d['median_ratio_including_excluded']:.3f}.")
    return "\n".join([
        r"\begin{table}[t]", r"\centering",
        *cfg.tex_caption("Two level-19 carriers, against the observation",
                         "tab:level19_direct"),
        r"\footnotesize",
        r"\begin{tabular}{@{}lrrrr@{}}", r"\toprule",
        r"Carrier & $n$ & Floor (RMSE) & $/$ level 20 & Lower floor on \\",
        r"\midrule", *body, r"\bottomrule", r"\end{tabular}",
        *cfg.tex_note(note), r"\end{table}"]) + "\n"


def write_outputs(df, subset=False, no_tex=False, archive=True) -> int:
    """Split the eligible regions from the excluded ones, then write."""
    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    tag = "__subset" if subset else ""
    out_csv = cfg.RESULTS_DIR / f"level19_direct{tag}.csv"
    # A full run replaces the file the paper's table is built from. Keep the
    # previous one: a run that dies on region 30 of 40, or a script change that
    # turns out to be wrong, should not cost the results already in hand.
    if archive and not subset and out_csv.exists():
        from datetime import datetime
        arch = cfg.RESULTS_DIR / "archive"
        arch.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
        for suffix in (".csv", ".json"):
            src = cfg.RESULTS_DIR / f"level19_direct{suffix}"
            if src.exists():
                shutil.copyfile(src, arch / f"level19_direct_{stamp}{suffix}")
        print(f"previous full run archived to {arch} as level19_direct_{stamp}.*")
        df.to_csv(out_csv, index=False)
    elif archive:
        df.to_csv(out_csv, index=False)
    # The CSV keeps every region that was measured; only the digest and the
    # table apply the window audit, so nothing is lost and the exclusion is
    # visible in one place.
    drop_names, audit_path = excluded_scenes()
    mask = df.scene.isin(drop_names)
    kept, dropped = df[~mask], df[mask]
    if not audit_path.exists():
        print(f"note: {audit_path.name} not found; no region excluded. Run "
              "window_radiometry.py -- this comparison rebuilds the transfer "
              "from the native window and is only valid where that window is "
              "the product the tile came from.")
    elif len(dropped):
        print(f"window audit excludes {len(dropped)} region(s): "
              + ", ".join(dropped.scene))
    digest = summarize(kept, dropped)
    (cfg.RESULTS_DIR / f"level19_direct{tag}.json").write_text(
        json.dumps(digest, indent=1))
    if subset:
        print(f"subset run: wrote level19_direct{tag}.csv; the full-run files "
              f"are untouched")
        no_tex = True                     # a table from a subset would mislead
    if not no_tex:
        cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
        (cfg.TEX_TABLES / "table_level19_direct.tex").write_text(
            latex_table(kept, digest), encoding="utf-8")
        print(f"wrote {cfg.TEX_TABLES / 'table_level19_direct.tex'}")
    print("\n" + json.dumps(digest, indent=1))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--diagnose", action="store_true",
                    help="record the ten most divergent cells, the moments of "
                         "both carriers, and the level-20 rebuild control")
    ap.add_argument("--no-control", action="store_true",
                    help="with --diagnose, skip the level-20 rebuild control "
                         "(it costs a second triangulation per region)")
    ap.add_argument("--dump", action="store_true",
                    help="with --diagnose, also write both carriers and the "
                         "cell centres to results/level19_diag__<scene>.npz")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--no-tex", action="store_true")
    ap.add_argument("--tables-only", action="store_true",
                    help="rebuild the digest and the table from the existing "
                         "level19_direct.csv, without recomputing a carrier; "
                         "use after the window audit changes which regions "
                         "are eligible")
    args = ap.parse_args(argv)

    if args.tables_only:
        src = cfg.RESULTS_DIR / "level19_direct.csv"
        if not src.exists():
            raise SystemExit(f"{src} does not exist; run without --tables-only")
        return write_outputs(pd.read_csv(src), subset=False,
                             no_tex=args.no_tex, archive=False)

    sites = region_sites.select(region_sites.load_sites(), args.regions,
                                None, args.limit)
    scenes = list(sites.scene)
    rows, failures = [], []
    for i, scene in enumerate(scenes, 1):
        if not build_tiles.tile_path(scene).exists():
            failures.append((scene, "no level-20 tile; run build_tiles.py"))
            continue
        try:
            r = region_row(scene, verbose=args.verbose, diagnose=args.diagnose,
                           dump=args.dump, control=not args.no_control)
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
        if args.diagnose:
            print(f"      direct  mean {r['direct_mean']:+.5f} std "
                  f"{r['direct_std']:.5f} range [{r['direct_min']:+.4f}, "
                  f"{r['direct_max']:+.4f}]  cells returned "
                  f"{r['n_returned']}/{r['n_cells']}")
            print(f"      aggreg. mean {r['agg_mean']:+.5f} std "
                  f"{r['agg_std']:.5f}  corr(direct, aggregate) = "
                  f"{r['corr_direct_agg']:+.4f}")
            if "rebuild20_rms" in r:
                print(f"      CONTROL: level-20 rebuild vs stored tile: RMS "
                      f"{r['rebuild20_rms']:.3e} = "
                      f"{100*r['rebuild20_rms_frac_of_std']:.2f}% of std, "
                      f"corr {r['rebuild20_corr']:+.5f}, "
                      f"{r['rebuild20_n_missing']} cell(s) completed")
    if not rows:
        for s, why in failures:
            print(f"  {s}: {why}")
        raise SystemExit("no region could be evaluated")

    write_outputs(pd.DataFrame(rows), subset=bool(args.regions or args.limit),
                  no_tex=args.no_tex)
    if failures:
        print(f"\n{len(failures)} region(s) could not be evaluated:")
        for s, why in failures:
            print(f"  {s}: {why}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
