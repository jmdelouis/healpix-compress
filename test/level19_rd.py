"""The level-19 carrier in the rate--distortion sweep (revision 5, run 2.2).

Every compressed-size ratio in the paper compares the level-20 carrier, which
oversamples the native product by 2.58, with the native crop. Level 19 is the
coarser neighbour, 0.65 times the native density, and it has only been
compared on reconstruction floors (Section V-K, level19_direct.py), never
compressed. This script closes that gap: it builds the level-19 carrier and
runs it through the same codec sweep as the level-20 tile, then reports its
size against the native crop and against the level-20 tile at every RMSE
target, under the paper's rules (oracle-best codec, log-log interpolation,
the oracle's codec bracketing the target on both sides).

What the result decides. If the coarser carrier is smaller than the native
crop at every target, the cost the paper measures for level 20 is the cost of
its density and nothing about the HEALPix lattice or the transfer stands in
the way of a carrier as such. If it stays larger at the tight targets, the
transfer costs bytes even when it discards samples, and the level-20 story
has to say so.

Two level-19 carriers are swept. `direct` is the Clough--Tocher interpolant
of the native samples evaluated at the level-19 centres, the honest object
(one triangulation per region, as in level19_direct.py). `aggregate` is the
mean of the four level-20 children of the stored tile, the cheap one;
Section V-K found the two within 0.2% in reconstruction floor, and this run
records how far apart they are in compressed size.

The within-representation RMSE (each carrier against its own raster, the
paper's main protocol) is the one the ratios use. The common-reference RMSE
(the decoded carrier read on the native pixels and scored against the
observation) is recorded too, with the level-19 read of common_reference.py,
for completeness; its floor is 1.8 times the level-20 floor and it is not
what this comparison is about.

Regions whose window is no longer the product their tile was built from
(window_radiometry.py, one region) are excluded from the direct carrier for
the reason level19_direct.py gives; the aggregate carrier needs no window and
is built for all 40.

Usage
-----
    python -B level19_rd.py                    # both carriers, all regions
    python -B level19_rd.py --no-direct        # aggregate only, no triangulation
    python -B level19_rd.py --tables-only      # recompute from the saved sweeps

Writes results/level19_rd_sweeps.csv (every sweep point), results/level19_rd.csv
(per region and target), results/level19_rd.json (the digest the checker
reads), results/level19_rd.md and ../tex/tables/table_level19_rd.tex.
"""

from __future__ import annotations

import argparse
import datetime
import json
import math
import shutil

import numpy as np
import pandas as pd

import paper_config as cfg
import paper_stats as st
import codec_lib
import build_tiles
import region_sites

COARSE_DELTA = 1
LEVEL19 = cfg.LEVEL - COARSE_DELTA
N19 = cfg.N_CHILDREN // (4 ** COARSE_DELTA)
SIDE19 = cfg.TILE_SIDE // (2 ** COARSE_DELTA)
CARRIERS = ("direct", "aggregate")


# --- carriers ----------------------------------------------------------------
def expanded(flat):
    from hcc.compression.zorder import expand_tile
    out = np.asarray(expand_tile(np.asarray(flat, dtype=np.float32)),
                     dtype=np.float32)
    if out.shape != (SIDE19, SIDE19):
        raise RuntimeError(f"level-{LEVEL19} raster is {out.shape}, expected "
                           f"{(SIDE19, SIDE19)}")
    return out


def flattened(raster):
    from hcc.compression.zorder import flatten_tile
    return np.asarray(flatten_tile(np.asarray(raster, dtype=np.float32)),
                      dtype=np.float32).ravel()


def carriers_for(scene, direct=True, verbose=False):
    """Level-19 carriers of one region, keyed by name, plus what the common
    reference needs."""
    import common_reference
    cell_ids, cell_data, utm_crop, meta, npz = build_tiles.load_tile(scene)
    flat20 = np.asarray(cell_data, dtype=np.float64).ravel()
    out = {"aggregate": flat20.reshape(-1, 4 ** COARSE_DELTA).mean(axis=1)}
    n_missing = 0
    if direct:
        import level19_direct
        store = cfg.REGION_DIR / f"{scene}_data.zarr"
        if not store.exists():
            raise FileNotFoundError(
                f"{scene}: {store} is missing; the direct carrier is rebuilt "
                f"from the native window (build_region_windows.py output)")
        sw = build_tiles.tile_subwindow(scene, store)
        b = level19_direct.build_at_level(scene, sw, COARSE_DELTA,
                                          verbose=verbose)
        out["direct"] = np.asarray(b["values"], dtype=np.float64)
        n_missing = int(b["n_missing"])
    cr = common_reference.build_for_tile(npz, mode="nearest")
    idx19 = np.clip(cr._safe_index // (4 ** COARSE_DELTA), 0, N19 - 1)
    return out, dict(cr=cr, idx19=idx19,
                     obs=np.asarray(utm_crop, dtype=np.float64),
                     n_missing=n_missing, tile_std=float(flat20.std()))


def sweep_carrier(scene, name, values, ref, codecs, verbose=False):
    raster = expanded(values)
    rows = []
    stats = {}
    for run in codec_lib.iter_runs(raster, codecs=codecs, on_error="warn",
                                   stats=stats):
        rec = np.asarray(run.reconstruction, dtype=np.float32)
        read19 = flattened(rec)[ref["idx19"]].reshape(ref["cr"].shape)
        rows.append(dict(
            scene=scene, scene_class=scene.split("__")[1], carrier=name,
            level=LEVEL19, codec=run.codec,
            config_json=json.dumps(run.config, sort_keys=True),
            total_bytes=int(run.total_bytes),
            rmse=codec_lib.rmse(raster, rec),
            rmse_common=ref["cr"].rmse_on_mask(read19, ref["obs"]),
            encode_seconds=run.encode_seconds, decode_seconds=run.decode_seconds,
            raw_bytes=int(raster.nbytes), n_cells=int(raster.size),
            n_missing=ref["n_missing"] if name == "direct" else 0))
    if not any(r["codec"] == "uncompressed" for r in rows):
        rows.append(dict(scene=scene, scene_class=scene.split("__")[1],
                         carrier=name, level=LEVEL19, codec="uncompressed",
                         config_json="{}", total_bytes=int(raster.nbytes),
                         rmse=0.0, rmse_common=np.nan, encode_seconds=0.0,
                         decode_seconds=0.0, raw_bytes=int(raster.nbytes),
                         n_cells=int(raster.size), n_missing=0))
    empty = [c for c, v in stats.items() if v.get("ok", 0) == 0]
    if empty and verbose:
        print(f"    WARNING [{name}] no usable point for: {', '.join(empty)}")
    if verbose:
        print(f"    {name}: {len(rows)} sweep points, raw {raster.nbytes/1e6:.3f} MB")
    return pd.DataFrame(rows)


# --- analysis ----------------------------------------------------------------
def load_reference_sweeps():
    utm = pd.read_csv(cfg.RESULTS_DIR / "utm_codec_bench.csv")
    hpx = pd.read_csv(cfg.RESULTS_DIR / "hpx_codec_bench.csv")
    return utm, hpx


def per_region_ratios(sweeps, utm, hpx, codecs, targets):
    rows = []
    for (scene, name), l19 in sweeps.groupby(["scene", "carrier"]):
        if scene not in set(utm.scene) or scene not in set(hpx.scene):
            continue
        raw19 = float(l19[l19.codec == "uncompressed"].total_bytes.iloc[0]) \
            if (l19.codec == "uncompressed").any() else np.nan
        raw_u = utm[(utm.scene == scene) & (utm.codec == "uncompressed")]
        raw_u = float(raw_u.total_bytes.iloc[0]) if len(raw_u) else np.nan
        for t in targets:
            s19 = st.best_of_codecs_size(l19, scene, codecs, t)
            su = st.best_of_codecs_size(utm, scene, codecs, t)
            s20 = st.best_of_codecs_size(hpx, scene, codecs, t)
            b19 = st.oracle_bracketed(l19, scene, codecs, t, "rmse")
            bu = st.oracle_bracketed(utm, scene, codecs, t, "rmse")
            b20 = st.oracle_bracketed(hpx, scene, codecs, t, "rmse")
            rows.append(dict(
                scene=scene, scene_class=scene.split("__")[1], carrier=name,
                target=t, l19_bytes=s19, utm_bytes=su, l20_bytes=s20,
                raw_l19_bytes=raw19, raw_utm_bytes=raw_u,
                bracketed_l19=bool(b19), bracketed_utm=bool(bu),
                bracketed_l20=bool(b20),
                eligible_utm=bool(b19 and bu), eligible_l20=bool(b19 and b20),
                ratio_l19_utm=(s19 / su if np.isfinite(s19) and np.isfinite(su)
                               and su > 0 else np.nan),
                ratio_l19_l20=(s19 / s20 if np.isfinite(s19) and np.isfinite(s20)
                               and s20 > 0 else np.nan),
                raw_ratio_l19_utm=(raw19 / raw_u if np.isfinite(raw19)
                                   and np.isfinite(raw_u) and raw_u > 0
                                   else np.nan)))
    return pd.DataFrame(rows)


def digest_of(df, targets):
    out = {"level": LEVEL19, "targets": [float(t) for t in targets],
           "raw_ratio_l19_utm": float(df.raw_ratio_l19_utm.median())
           if df.raw_ratio_l19_utm.notna().any() else None,
           "carriers": {}}
    for name, d in df.groupby("carrier"):
        entry = {"n_regions": int(d.scene.nunique()), "targets": {}}
        p_utm, p_l20 = [], []
        for t in targets:
            g = d[d.target == t]
            a = g[g.eligible_utm].ratio_l19_utm
            b = g[g.eligible_l20].ratio_l19_l20
            ta = st.paired_ratio_test(a)
            tb = st.paired_ratio_test(b)
            p_utm.append(ta["p"])
            p_l20.append(tb["p"])
            entry["targets"][f"{t:g}"] = {
                "l19_over_utm": {k: ta[k] for k in ("n", "median", "lo", "hi",
                                                    "ci_lo", "ci_hi", "p")},
                "l19_over_utm_win_fraction": (float((a < 1).mean())
                                              if len(a) else None),
                "l19_over_l20": {k: tb[k] for k in ("n", "median", "lo", "hi",
                                                    "ci_lo", "ci_hi", "p")},
                "bracketed_fraction_utm": float(g.eligible_utm.mean())
                if len(g) else None,
            }
        for key, ps in (("l19_over_utm", p_utm), ("l19_over_l20", p_l20)):
            hs = st.holm(ps)
            for t, h in zip(targets, hs):
                entry["targets"][f"{t:g}"][key]["p_holm"] = (
                    float(h) if np.isfinite(h) else None)
        out["carriers"][name] = entry
    if {"direct", "aggregate"} <= set(df.carrier.unique()):
        piv = df.pivot_table(index=["scene", "target"], columns="carrier",
                             values="l19_bytes")
        r = (piv["direct"] / piv["aggregate"]).dropna()
        out["direct_over_aggregate_size"] = {
            "median": float(r.median()), "lo": float(r.min()),
            "hi": float(r.max()), "n": int(r.size)}
        # per target, on the regions eligible for both carriers (revision 6,
        # run 2.4): the aggregate's advantage grows with the target, and its
        # compressed / raw size tracks the native crop's, which is why its
        # L19/UTM column is flat at the raw ratio.
        raw_u = df.pivot_table(index=["scene", "target"], values="raw_utm_bytes")
        ub = df.pivot_table(index=["scene", "target"], values="utm_bytes")
        raw19 = df.pivot_table(index=["scene", "target"], values="raw_l19_bytes")
        elig = df.pivot_table(index=["scene", "target"], columns="carrier",
                              values="eligible_utm")
        per = {}
        for t in targets:
            e = elig.xs(t, level="target")
            keep = e.index[(e.get("direct", False) == True) & (e.get("aggregate", False) == True)]
            if not len(keep):
                continue
            b = piv.xs(t, level="target").loc[keep]
            ad = (b["aggregate"] / b["direct"]).dropna()
            cu = (ub.xs(t, level="target").loc[keep, "utm_bytes"]
                  / raw_u.xs(t, level="target").loc[keep, "raw_utm_bytes"])
            r19 = raw19.xs(t, level="target").loc[keep, "raw_l19_bytes"]
            per[f"{t:g}"] = {
                "n": int(ad.size),
                "aggregate_over_direct_median": float(ad.median()),
                "aggregate_cheaper_fraction": float((ad < 1).mean()),
                "compact_utm": float(cu.median()),
                "compact_direct": float((b["direct"] / r19).median()),
                "compact_aggregate": float((b["aggregate"] / r19).median()),
            }
        out["aggregate_vs_direct"] = per
    return out


def _fmt_p(p):
    if p is None or not np.isfinite(p):
        return "--"
    return "$<10^{-4}$" if p < 1e-4 else f"{p:.3f}"


def latex_table(digest, targets):
    main = "direct" if "direct" in digest["carriers"] else "aggregate"
    other = "aggregate" if main == "direct" else None
    if other is not None and other not in digest["carriers"]:
        other = None
    c = digest["carriers"][main]
    cols = "lrrlrr" + ("r" if other else "")
    head = ("Target & $n$ & L19/UTM & Range & $p_{\\mathrm{Holm}}$ & L19/L20"
            + (" & L19/UTM (agg.)" if other else "") + " \\\\")
    rows = []
    for t in targets:
        e = c["targets"][f"{t:g}"]
        a, b = e["l19_over_utm"], e["l19_over_l20"]
        if a["n"] == 0 or (e["bracketed_fraction_utm"] or 0) < cfg.MIN_BRACKETED:
            continue
        row = (f"${t:g}$ & {a['n']} & {a['median']:.2f} & "
               f"[{a['lo']:.2f}, {a['hi']:.2f}] & {_fmt_p(a.get('p_holm'))} & "
               f"{b['median']:.2f}")
        if other:
            o = digest["carriers"][other]["targets"][f"{t:g}"]["l19_over_utm"]
            row += f" & {o['median']:.2f}"
        rows.append(row + " \\\\")
    if not rows:
        return ""
    raw = digest.get("raw_ratio_l19_utm")
    raw_txt = f", {raw:.2f} times the native crop" if raw else ""
    note = (
        f"The level-{LEVEL19} carrier (${SIDE19}^2$ cells, "
        f"{(4*N19)/1e6:.3f}\\,MB raw{raw_txt}) "
        "compressed under the same protocol as the level-20 tile: "
        "oracle-best codec, log--log interpolation, within-representation "
        "RMSE, $n$ the regions where the oracle's codec brackets the target "
        "on both sides. L19/UTM is the compressed-size ratio to the native "
        "crop at $\\gamma_{\\mathrm{tile}}=1$; L19/L20 the ratio to the "
        "level-20 tile. The main columns use the carrier evaluated directly "
        "at the level-19 centres by the Clough--Tocher transfer"
        + ("; the last column uses the mean of the four level-20 children "
           "of the stored tile." if other else "."))
    return (
        "\\begin{table}[t]\n\\centering\n"
        + "\n".join(cfg.tex_caption(
            f"The level-{LEVEL19} carrier under compression",
            "tab:level19_rd"))
        + "\n\\footnotesize\n\\setlength{\\tabcolsep}{3pt}\n"
        + f"\\begin{{tabular}}{{@{{}}{cols}@{{}}}}\n\\toprule\n{head}\n\\midrule\n"
        + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n"
        + "\n".join(cfg.tex_note(note)) + "\n\\end{table}\n")


def summarize(digest, targets):
    out = [f"# Level-{LEVEL19} carrier under compression", ""]
    raw = digest.get("raw_ratio_l19_utm")
    if raw:
        out.append(f"Raw size ratio L19/UTM: {raw:.3f} (level 20: 2.584).")
        out.append("")
    for name, c in digest["carriers"].items():
        out += [f"## carrier: {name} ({c['n_regions']} regions)", "",
                "| target | n | L19/UTM median | 95% CI | p | p_Holm | "
                "win frac | L19/L20 median | n |",
                "|---|---|---|---|---|---|---|---|---|"]
        for t in targets:
            e = c["targets"][f"{t:g}"]
            a, b = e["l19_over_utm"], e["l19_over_l20"]
            wf = e["l19_over_utm_win_fraction"]
            out.append(
                f"| {t:g} | {a['n']} | {a['median']:.3f} | "
                f"[{a['ci_lo']:.3f}, {a['ci_hi']:.3f}] | {_fmt_p(a['p'])} | "
                f"{_fmt_p(a.get('p_holm'))} | "
                f"{'--' if wf is None else f'{wf:.0%}'} | "
                f"{b['median']:.3f} | {b['n']} |")
        out.append("")
    d = digest.get("direct_over_aggregate_size")
    if d:
        out.append(f"Direct/aggregate compressed size at equal RMSE: median "
                   f"{d['median']:.3f}, range [{d['lo']:.3f}, {d['hi']:.3f}] "
                   f"over {d['n']} (region, target) pairs.")
        out.append("")
    return "\n".join(out)


# --- driver ------------------------------------------------------------------
def _archive(names):
    stamp = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")
    for name in names:
        src = cfg.RESULTS_DIR / name
        if src.exists():
            dst = cfg.RESULTS_DIR / "archive" / f"{src.stem}.{stamp}{src.suffix}"
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            print(f"  archived {name} -> archive/{dst.name}")


def write_outputs(sweeps, codecs, targets, no_tex=False):
    utm, hpx = load_reference_sweeps()
    df = per_region_ratios(sweeps, utm, hpx, codecs, targets)
    if not len(df):
        print("no region present in the level-19 sweeps and both reference sweeps")
        return 1
    digest = digest_of(df, targets)
    digest["codecs"] = list(codecs)
    df.to_csv(cfg.RESULTS_DIR / "level19_rd.csv", index=False)
    (cfg.RESULTS_DIR / "level19_rd.json").write_text(json.dumps(digest, indent=1))
    md = summarize(digest, targets)
    (cfg.RESULTS_DIR / "level19_rd.md").write_text(md, encoding="utf-8")
    print("\n" + md)
    if not no_tex:
        tex = latex_table(digest, targets)
        if tex:
            cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
            (cfg.TEX_TABLES / "table_level19_rd.tex").write_text(tex, encoding="utf-8")
            print(f"wrote {cfg.TEX_TABLES / 'table_level19_rd.tex'}")
    return 0


def diagnose(sweeps, codecs, targets):
    """Why is L19/UTM flat for the aggregate carrier? (revision 6, run 2.4)

    Review R1 of v6 noticed that the aggregate column of the level-19 table
    stays at 0.63--0.64 over the five targets, i.e. at the raw ratio 0.65,
    while the direct carrier rises from 0.67 to 0.78. This prints, per
    target, what could produce that: the populations (a different n per
    carrier), the pairing (agg/direct on the same regions), the
    compactability of each raster (compressed / own lossless size, so the
    aggregate can be compared with the crop on a ratio that divides the raw
    size out), and which codec the oracle picks. It computes nothing new
    from the tiles: everything comes from level19_rd_sweeps.csv and the two
    published sweeps.
    """
    utm, hpx = load_reference_sweeps()
    df = per_region_ratios(sweeps, utm, hpx, codecs, targets)
    lines = [f"# Level-{LEVEL19} diagnostic: aggregate against direct", ""]
    piv_b = df.pivot_table(index=["scene", "target"], columns="carrier",
                           values="l19_bytes")
    piv_e = df.pivot_table(index=["scene", "target"], columns="carrier",
                           values="eligible_utm")
    lines += ["| target | n direct | n agg | direct/UTM | agg/UTM | "
              "agg/direct (paired, same regions) | agg < direct | "
              "compact. UTM | compact. direct | compact. agg |",
              "|---|---|---|---|---|---|---|---|---|---|"]

    def compact(sw, scene, t):
        """compressed size at target / raw size of the same raster (the
        level-19 sweeps carry no lossless point, so raw is the common base)"""
        s = st.best_of_codecs_size(sw, scene, codecs, t)
        g = sw[sw.scene == scene]
        base = g.raw_bytes.iloc[0] if len(g) and "raw_bytes" in g else np.nan
        return s / base if np.isfinite(s) and np.isfinite(base) and base > 0 else np.nan

    for t in targets:
        d = df[df.target == t]
        nd = int(d[(d.carrier == "direct") & d.eligible_utm].scene.nunique())
        na = int(d[(d.carrier == "aggregate") & d.eligible_utm].scene.nunique())
        md_ = d[(d.carrier == "direct") & d.eligible_utm].ratio_l19_utm.median()
        ma_ = d[(d.carrier == "aggregate") & d.eligible_utm].ratio_l19_utm.median()
        both = piv_e.xs(t, level="target")
        both = both[(both.get("direct", False) == True) & (both.get("aggregate", False) == True)]
        bb = piv_b.xs(t, level="target").loc[both.index] if len(both) else piv_b.iloc[0:0]
        r = (bb["aggregate"] / bb["direct"]).dropna() if len(bb) else pd.Series(dtype=float)
        cu, cd, ca = [], [], []
        for scene in both.index:
            cu.append(compact(utm, scene, t))
            l19 = sweeps[sweeps.scene == scene]
            cd.append(compact(l19[l19.carrier == "direct"], scene, t))
            ca.append(compact(l19[l19.carrier == "aggregate"], scene, t))
        fmt = lambda v: f"{np.nanmedian(v):.3f}" if len(v) and np.isfinite(np.nanmedian(v)) else "--"
        lines.append(f"| {t:g} | {nd} | {na} | {md_:.3f} | {ma_:.3f} | "
                     f"{fmt(r.values) if len(r) else '--'} (n={len(r)}) | "
                     f"{(r < 1).mean():.2f} | {fmt(cu)} | {fmt(cd)} | {fmt(ca)} |"
                     if len(r) else
                     f"| {t:g} | {nd} | {na} | {md_:.3f} | {ma_:.3f} | -- | -- | "
                     f"{fmt(cu)} | {fmt(cd)} | {fmt(ca)} |")
    # which codec the oracle selects, per carrier and target
    lines += ["", "## Oracle codec per carrier (share of regions)", ""]
    for name, sw in sweeps.groupby("carrier"):
        for t in targets:
            picks = []
            for scene in sw.scene.unique():
                best, bsize = None, np.inf
                for c in codecs:
                    v = st.best_of_codecs_size(sw, scene, [c], t)
                    if np.isfinite(v) and v < bsize:
                        best, bsize = c, v
                if best:
                    picks.append(best)
            vc = pd.Series(picks).value_counts(normalize=True)
            lines.append(f"- {name}, target {t:g}: "
                         + ", ".join(f"{k} {100*v:.0f}%" for k, v in vc.items()))
    # sanity: the aggregate's own RMSE reference is the aggregate raster
    agg = sweeps[(sweeps.carrier == "aggregate") & (sweeps.codec != "uncompressed")]
    lines += ["", "## Reference check",
              f"- aggregate sweep rows: {len(agg)}, raw_bytes per row "
              f"{sorted(agg.raw_bytes.unique())} (expected {4 * N19} for a "
              f"{SIDE19}^2 float32 raster: the RMSE column is against that raster)",
              "- compactability columns are compressed / raw bytes of each raster; "
              "the aggregate is a 2x2 mean of the Clough--Tocher tile, hence "
              "smoother than the directly evaluated carrier, and a flat agg/UTM "
              "column with a falling agg/direct column means the aggregate gains "
              "compressibility with the target at the same rate as the native crop",
              ]
    md = "\n".join(lines)
    out = cfg.RESULTS_DIR / "level19_rd_diagnose.md"
    out.write_text(md, encoding="utf-8")
    print(md)
    print(f"\nwrote {out}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--codecs", nargs="*", default=cfg.COMPRESSOR_CODECS)
    ap.add_argument("--targets", type=float, nargs="*", default=cfg.TARGET_RMSE)
    ap.add_argument("--no-direct", action="store_true",
                    help="sweep the aggregate carrier only (no triangulation)")
    ap.add_argument("--tables-only", action="store_true",
                    help="recompute ratios, digest and table from the saved "
                         "level19_rd_sweeps.csv without sweeping")
    ap.add_argument("--no-tex", action="store_true")
    ap.add_argument("--diagnose", action="store_true",
                    help="explain the flat aggregate column from the saved "
                         "sweeps (revision 6, run 2.4); writes "
                         "results/level19_rd_diagnose.md, sweeps nothing")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args(argv)

    sweeps_path = cfg.RESULTS_DIR / "level19_rd_sweeps.csv"
    if args.diagnose:
        if not sweeps_path.exists():
            print(f"{sweeps_path} is missing; run without --diagnose first")
            return 1
        return diagnose(pd.read_csv(sweeps_path), args.codecs, args.targets)
    if args.tables_only:
        if not sweeps_path.exists():
            print(f"{sweeps_path} is missing; run without --tables-only first")
            return 1
        return write_outputs(pd.read_csv(sweeps_path), args.codecs,
                             args.targets, no_tex=args.no_tex)

    scenes = args.regions or list(region_sites.load_sites()["scene"])
    scenes = [s for s in scenes if (cfg.TILE_DIR / f"{s}.npz").exists()]
    if args.limit:
        scenes = scenes[:args.limit]
    if not scenes:
        print(f"no tiles under {cfg.TILE_DIR}; run build_tiles.py first")
        return 1
    excluded = set(cfg.window_mismatched_scenes()) if not args.no_direct else set()
    if excluded:
        print(f"direct carrier skipped for {sorted(excluded)} (window is not "
              f"the product the tile was built from; window_radiometry.py)")

    _archive(["level19_rd_sweeps.csv", "level19_rd.csv", "level19_rd.json",
              "level19_rd.md"])
    print(f"{len(scenes)} region(s), codecs {args.codecs}, "
          f"carriers {'aggregate only' if args.no_direct else 'direct + aggregate'}")
    frames, failures = [], []
    for scene in scenes:
        print(f"  {scene}")
        try:
            direct = not args.no_direct and scene not in excluded
            cars, ref = carriers_for(scene, direct=direct, verbose=args.verbose)
            for name, values in cars.items():
                frames.append(sweep_carrier(scene, name, values, ref,
                                            args.codecs, verbose=True))
        except Exception as exc:                                # noqa: BLE001
            print(f"    FAILED {type(exc).__name__}: {exc}")
            failures.append((scene, f"{type(exc).__name__}: {exc}"))
            continue
        # written after every region so an interrupted run keeps its work
        cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        pd.concat(frames, ignore_index=True).to_csv(sweeps_path, index=False)
    if not frames:
        print("no sweep produced any row")
        return 1
    rc = write_outputs(pd.concat(frames, ignore_index=True), args.codecs,
                       args.targets, no_tex=args.no_tex)
    if failures:
        print(f"\n{len(failures)} region(s) failed:")
        for s, why in failures:
            print(f"  {s}: {why}")
        return 2
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
