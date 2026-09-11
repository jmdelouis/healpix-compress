"""Recompute every number the paper quotes, and rewrite the LaTeX tables.

Reads results/hpx_codec_bench.csv + results/utm_codec_bench.csv (whatever
number of regions they contain) and writes:

    ../tex/tables/table_storage_accounting.tex        (Table III)
    ../tex/tables/table_storage_at_constant_rmse.tex  (Table IV)
    ../tex/tables/table_error_at_constant_size.tex    (Table V)
    ../tex/tables/table_gamma_star.tex                (critical overhead)
    ../tex/tables/table_codec_ranking.tex             (Table VI, audited counts)
    ../tex/tables/table_benchmark_regions.tex         (dataset summary)
    ../tex/tables/table_latitude_effect.tex           (new: HEALPix regime)
    ../tex/tables/table_fixed_codec.tex               (new: same-codec control)
    results/paper_numbers.json                        (machine-readable)
    results/paper_numbers.md                          (digest for the redaction)

    python update_paper_numbers.py                    # self RMSE (current paper)
    python update_paper_numbers.py --metric common    # common UTM reference
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import paper_config as cfg
import paper_stats as st

RMSE_LABEL = {1e-4: r"$10^{-4}$", 3e-4: r"$3\times10^{-4}$", 1e-3: r"$10^{-3}$",
              3e-3: r"$3\times10^{-3}$", 1e-2: r"$10^{-2}$"}


def _p(v):
    """p-values: never print 0.000, which reads as exactly zero."""
    if not np.isfinite(v):
        return "--"
    if v < 1e-4:
        return r"$<10^{-4}$"
    if v < 1e-3:
        return r"$<10^{-3}$"
    return f"{v:.3f}"


def _fmt(v, nd=2):
    return "--" if not np.isfinite(v) else f"{v:.{nd}f}"


def _table(caption, label, colspec, header, body_lines, note=None, wide=False):
    """Assemble one LaTeX table.

    `wide` emits a two-column float. A nine-column table does not fit an IEEE
    text column, and under `nonstopmode` pdflatex does not stop for that: it
    sets an overfull box and, if the declared column count is also short, drops
    a column without complaint. So the width is decided here rather than
    noticed in the proof.
    """
    env = "table*" if wide else "table"
    out = [f"\\begin{{{env}}}[t]", r"\centering", f"\\caption{{{caption}}}",
           f"\\label{{{label}}}"]
    if wide:
        out.append(r"\footnotesize")
    out += [f"\\begin{{tabular}}{{{colspec}}}", r"\toprule",
            header, r"\midrule", *body_lines, r"\bottomrule", r"\end{tabular}"]
    if note:
        out.append(f"\\\\[2pt]\\footnotesize {note}")
    out.append(f"\\end{{{env}}}")
    return "\n".join(out) + "\n"


def check_codec_symmetry(hpx, utm, scenes, codecs):
    """Refuse to compare grids that were not offered the same codecs.

    The oracle-best protocol picks the smallest result over the available
    codecs on each side independently. If a codec produced points on one grid
    and none on the other, the reported ratio partly measures which codec
    happened to run, not which grid compresses better -- and it does so in a
    direction nobody controls. This is cheap to check and expensive to miss.
    """
    problems = []
    for scene in scenes:
        h = set(hpx[(hpx.scene == scene)].codec.unique()) & set(codecs)
        u = set(utm[(utm.scene == scene)].codec.unique()) & set(codecs)
        if h != u:
            problems.append((scene, sorted(u - h), sorted(h - u)))
    return problems


def quantization_floor(df, scenes):
    """RMSE reached by the nominally lossless codecs, per grid.

    The image codecs are fed a 16-bit requantization of the float32 tile, so
    their `lossless` points are lossless *after* that step, not bit-exact
    against the original. The residual is the floor every codec on that side
    inherits, and it belongs in the paper next to the tightest RMSE target.
    """
    sub = df[df.codec.isin(cfg.LOSSLESS_CODECS) & df.scene.isin(scenes)]
    if sub.empty:
        return {}
    return dict(median=float(sub.rmse.median()), min=float(sub.rmse.min()),
                max=float(sub.rmse.max()), n=int(len(sub)))


def load(metric="self"):
    hpx = pd.read_csv(cfg.RESULTS_DIR / "hpx_codec_bench.csv")
    utm = pd.read_csv(cfg.RESULTS_DIR / "utm_codec_bench.csv")
    rmse_col = "rmse"
    if metric == "common":
        if "rmse_common" not in hpx.columns or hpx["rmse_common"].isna().all():
            raise SystemExit(
                "--metric common needs the rmse_common column: rerun\n"
                "    python run_codec_bench.py --grid hpx --common-reference --force")
        rmse_col = "rmse_common"
    scenes = sorted(set(hpx.scene) & set(utm.scene))
    if not scenes:
        raise SystemExit("no scene present on both grids")
    return hpx, utm, scenes, rmse_col


def table_storage_accounting(utm, scenes):
    raw = utm[utm.codec == "uncompressed"].set_index("scene")["total_bytes"].to_dict()
    ratios_g = np.array([cfg.SHARED_RAW / (raw[s] * cfg.GAMMA_UTM)
                         for s in scenes if s in raw])
    ratios_1 = np.array([cfg.SHARED_RAW / raw[s] for s in scenes if s in raw])
    t_g, t_1 = st.paired_ratio_test(ratios_g), st.paired_ratio_test(ratios_1)
    body = [
        f"HEALPix tile (level {cfg.LEVEL}, ${cfg.TILE_SIDE}^2$ cells) & "
        f"{cfg.SHARED_RAW/1e6:.3f} & 1.000 \\\\",
        f"Equal-area UTM crop (${cfg.UTM_SIDE_PX}^2$ px) & "
        f"{np.mean(list(raw.values()))/1e6:.3f} & {t_1['median']:.3f} \\\\",
        f"UTM crop $\\times\\gamma_{{\\mathrm{{tile}}}}={cfg.GAMMA_UTM}$ & "
        f"{np.mean(list(raw.values()))*cfg.GAMMA_UTM/1e6:.3f} & {t_g['median']:.3f} \\\\",
    ]
    tex = _table(
        f"Raw, uncompressed storage for the same nominal ground area, across the "
        f"{len(ratios_g)} regions. Both representations are fixed-size, so the ratio is a "
        f"geometric constant (measured range {t_1['lo']:.3f}--{t_1['hi']:.3f} at "
        f"$\\gamma_{{\\mathrm{{tile}}}}=1$).",
        "tab:storage_accounting", "lrr",
        r"Representation & Raw size (MB) & HEALPix / this \\", body)
    return tex, dict(n=int(t_g["n"]), ratio_gamma1=t_1, ratio_gamma=t_g,
                     utm_raw_bytes_mean=float(np.mean(list(raw.values()))))


def _ratio_rows(hpx, utm, scenes, rmse_col, codecs):
    """Per-target, per-scene HPX/UTM ratios for both protocols."""
    rows = []
    for target in cfg.TARGET_RMSE:
        per_scene = {}
        for s in scenes:
            h = st.best_of_codecs_size(hpx, s, codecs, target, rmse_col=rmse_col)
            u = st.best_of_codecs_size(utm, s, codecs, target, gamma=cfg.GAMMA_UTM)
            u1 = st.best_of_codecs_size(utm, s, codecs, target, gamma=1.0)
            per_scene[s] = (h, u, u1)
        rows.append(dict(kind="rmse", target=target, per_scene=per_scene))
    for frac in cfg.BUDGETS:
        budget = frac * cfg.SHARED_RAW
        per_scene = {}
        for s in scenes:
            h = st.best_of_codecs_rmse(hpx, s, codecs, budget, rmse_col=rmse_col)
            u = st.best_of_codecs_rmse(utm, s, codecs, budget / cfg.GAMMA_UTM)
            u1 = st.best_of_codecs_rmse(utm, s, codecs, budget)
            per_scene[s] = (h, u, u1)
        rows.append(dict(kind="budget", target=frac, per_scene=per_scene))
    for r in rows:
        kind = r["kind"]
        tgt = r["target"] if kind == "rmse" else r["target"] * cfg.SHARED_RAW
        br = [st.any_bracketed(hpx, s, codecs, tgt, kind, rmse_col)
              and st.any_bracketed(utm, s, codecs, tgt, kind)
              for s in scenes]
        r["bracketed_fraction"] = float(np.mean(br)) if br else np.nan
        r["ratios"] = {s: (h / u) for s, (h, u, _) in r["per_scene"].items()
                       if np.isfinite(h) and np.isfinite(u) and u > 0}
        r["ratios_gamma1"] = {s: (h / u1) for s, (h, _, u1) in r["per_scene"].items()
                              if np.isfinite(h) and np.isfinite(u1) and u1 > 0}
        r["test"] = st.paired_ratio_test(list(r["ratios"].values()))
        r["test_gamma1"] = st.paired_ratio_test(list(r["ratios_gamma1"].values()))
    dropped = [r for r in rows
               if r["test"]["n"] >= cfg.MIN_SCENES
               and not (r["bracketed_fraction"] >= cfg.MIN_BRACKETED)]
    for r in dropped:
        print(f"  target {r['kind']} {r['target']:g} not reported: only "
              f"{100*r['bracketed_fraction']:.0f}% of regions bracket it "
              f"(threshold {100*cfg.MIN_BRACKETED:.0f}%); it would be a set of "
              f"bounds, not an estimate")
    reported = [r for r in rows
                if r["test"]["n"] >= cfg.MIN_SCENES
                and r["bracketed_fraction"] >= cfg.MIN_BRACKETED]
    for r in rows:
        r["reported"] = r in reported
    adj = st.holm([r["test"]["p"] for r in reported])
    for r, a in zip(reported, adj):
        r["p_holm"] = float(a) if np.isfinite(a) else np.nan
    for r in rows:
        r.setdefault("p_holm", np.nan)
    return rows, len(reported)


def table_iv_v(rows, family_size, n_regions):
    def lines(kind, fmt_target):
        out = []
        for r in rows:
            if r["kind"] != kind or not r.get("reported"):
                continue
            t = r["test"]
            favors = "HPX" if t["median"] < 1 else "UTM"
            out.append(f"{fmt_target(r['target'])} & {t['n']:d} & {_fmt(t['median'])} & "
                       f"[{_fmt(t['lo'])}, {_fmt(t['hi'])}] & "
                       f"{_p(t['p'])} & {_p(r['p_holm'])} & {favors} \\\\")
        return out

    note = (f"Exact/asymptotic two-sided Wilcoxon on $\\log$(ratio); "
            f"$p_{{\\mathrm{{Holm}}}}$ corrects across the whole family of "
            f"{family_size} targets reported in Tables~\\ref{{tab:storage_at_rmse}} "
            f"and~\\ref{{tab:error_at_size}}.")
    header = (r"Target & $n$ & Median ratio & Range & $p$ & $p_{\mathrm{Holm}}$ "
              r"& Favors \\")
    tex_iv = _table(
        f"Per-region HPX/UTM compressed-size ratio at constant reconstruction RMSE "
        f"(oracle-best codec chosen independently per region per grid; the UTM side "
        f"already includes $\\gamma_{{\\mathrm{{tile}}}}={cfg.GAMMA_UTM}$). Ratio $<1$ "
        f"favors HEALPix. Median and [min, max] across the $n$ evaluable regions "
        f"(of {n_regions}).",
        "tab:storage_at_rmse", "lrrrrrl", header,
        lines("rmse", lambda t: RMSE_LABEL.get(t, f"${t:g}$")), note)
    tex_v = _table(
        f"Per-region HPX/UTM reconstruction-RMSE ratio at a constant shared storage "
        f"budget, expressed as a fraction of the shared raw HEALPix tile size "
        f"({cfg.SHARED_RAW} B); the UTM side receives that budget divided by "
        f"$\\gamma_{{\\mathrm{{tile}}}}={cfg.GAMMA_UTM}$. Ratio $<1$ favors HEALPix.",
        "tab:error_at_size", "lrrrrrl", header,
        lines("budget", lambda t: f"${100*t:g}\\%$"), note)
    return tex_iv, tex_v


def table_gamma_star(rows):
    body = []
    gammas = {}
    for r in rows:
        if r["kind"] != "rmse" or not r.get("reported"):
            continue
        g = r["test_gamma1"]["median"]
        gammas[r["target"]] = g
        cmp_ = (r"$\gamma_\star > \gamma_{\mathrm{tile}}$" if g > cfg.GAMMA_UTM
                else r"$\gamma_\star < \gamma_{\mathrm{tile}}$")
        label = RMSE_LABEL.get(r["target"], f"${r['target']:g}$")
        body.append(f"{label} & {_fmt(g)} & {cmp_} & "
                    f"{'UTM' if g > cfg.GAMMA_UTM else 'HPX'} \\\\")
    tex = _table(
        r"Critical operational overhead $\gamma_\star$ at which HEALPix and native UTM "
        r"reach equal compressed size, by RMSE target (median across regions, "
        r"oracle-best codec per side, $\gamma_{\mathrm{tile}}=1$). HEALPix is smaller "
        r"than the operational UTM comparison only where the real overhead of the "
        r"native product exceeds $\gamma_\star$.",
        "tab:gamma_star", "lrrr",
        r"Target RMSE & $\gamma_\star$ & vs.\ $\gamma_{\mathrm{tile}}$ & Favors \\",
        body)
    return tex, gammas


def table_codec_ranking(hpx, scenes, rmse_col, codecs):
    points, wins, evaluable, ranks = 0, {}, {}, {}
    for target in cfg.TARGET_RMSE:
        for s in scenes:
            vals = {c: st.interp_size_at_rmse(hpx[(hpx.scene == s) & (hpx.codec == c)],
                                              target, rmse_col) for c in codecs}
            vals = {c: v for c, v in vals.items() if np.isfinite(v)}
            if not vals:
                continue
            points += 1
            order = sorted(vals, key=vals.get)
            wins[order[0]] = wins.get(order[0], 0) + 1
            for i, c in enumerate(order, 1):
                evaluable[c] = evaluable.get(c, 0) + 1
                ranks.setdefault(c, []).append(i)
    for frac in cfg.BUDGETS:
        budget = frac * cfg.SHARED_RAW
        for s in scenes:
            vals = {c: st.interp_rmse_at_size(hpx[(hpx.scene == s) & (hpx.codec == c)],
                                              budget, rmse_col) for c in codecs}
            vals = {c: v for c, v in vals.items() if np.isfinite(v)}
            if not vals:
                continue
            points += 1
            order = sorted(vals, key=vals.get)
            wins[order[0]] = wins.get(order[0], 0) + 1
            for i, c in enumerate(order, 1):
                evaluable[c] = evaluable.get(c, 0) + 1
                ranks.setdefault(c, []).append(i)
    total_wins = sum(wins.values())
    assert total_wins == points, (
        f"codec-ranking bookkeeping broken: {total_wins} wins for {points} points")
    body = []
    for c in sorted(codecs, key=lambda c: -wins.get(c, 0)):
        body.append(f"\\texttt{{{c.replace('_', chr(92)+'_')}}} & {wins.get(c,0):d} & "
                    f"{evaluable.get(c,0):d} & "
                    f"{np.mean(ranks.get(c, [np.nan])):.2f} \\\\")
    tex = _table(
        f"Codec ranking on the HEALPix side. An operating point is one (region, "
        f"target) pair over the {len(scenes)} regions and the "
        f"{len(cfg.TARGET_RMSE)} RMSE plus {len(cfg.BUDGETS)} budget targets "
        f"({len(scenes)*(len(cfg.TARGET_RMSE)+len(cfg.BUDGETS))} in total), "
        f"retained when at least one codec is evaluable there, which leaves "
        f"{points}. Wins sum to the number of points by construction: "
        f"each point has exactly one best codec. `Evaluable' counts the points where a "
        f"codec produced a value at all.",
        "tab:codec_ranking", "lrrr",
        r"Codec & Wins & Evaluable points & Mean rank \\", body)
    return tex, dict(points=points, wins=wins, evaluable=evaluable,
                     mean_rank={c: float(np.mean(v)) for c, v in ranks.items()})


def table_latitude(rows, index):
    """Does the HPX/UTM gap depend on where the tile sits on the sphere?"""
    lat = index.set_index("scene")["abs_lat"].to_dict() if "abs_lat" in index else {}
    regime = index.set_index("scene")["hpx_regime"].to_dict() if "hpx_regime" in index else {}
    body, stats_out, order = [], {}, []
    # Half the sample at least: the gamma=1 variant of the loosest budget is
    # out of interpolation range for most regions, and a correlation on the
    # seven that remain would be noise dressed as a finding.
    min_n = max(cfg.MIN_SCENES, 20)
    for r in rows:
        if (not r.get("reported")) or r["test_gamma1"]["n"] < min_n or not lat:
            continue
        scenes = list(r["ratios_gamma1"])
        vals = [r["ratios_gamma1"][s] for s in scenes]
        lats = [lat.get(s, np.nan) for s in scenes]
        sp = st.spearman(lats, vals)
        eq = [v for s, v in zip(scenes, vals) if regime.get(s) == "equatorial"]
        cap = [v for s, v in zip(scenes, vals) if regime.get(s) == "cap"]
        mw = st.mannwhitney(eq, cap)
        key = (f"{r['kind']}:{r['target']:g}")
        stats_out[key] = dict(spearman=sp, equatorial_vs_cap=mw)
        order.append(key)
        label = (RMSE_LABEL.get(r["target"], f"${r['target']:g}$") if r["kind"] == "rmse"
                 else f"${100*r['target']:g}\\%$")
        body.append(f"{'RMSE ' if r['kind']=='rmse' else 'Budget '}{label} & "
                    f"{sp['n']:d} & {_fmt(sp['rho'])} & {_p(sp['p'])} & "
                    f"PHOLM_RHO & {_fmt(mw['median_a'])} & "
                    f"{_fmt(mw['median_b'])} & {_p(mw['p'])} & PHOLM_MW \\\\")
    # This is a second family of tests, and it needs its own correction: the
    # Holm adjustment applied to the storage targets says nothing about it.
    adj_rho = st.holm([stats_out[k]["spearman"]["p"] for k in order])
    adj_mw = st.holm([stats_out[k]["equatorial_vs_cap"]["p"] for k in order])
    for k, a, b in zip(order, adj_rho, adj_mw):
        stats_out[k]["spearman"]["p_holm"] = float(a) if np.isfinite(a) else np.nan
        stats_out[k]["equatorial_vs_cap"]["p_holm"] = (float(b) if np.isfinite(b)
                                                       else np.nan)
    body = [ln.replace("PHOLM_RHO", _p(a)).replace("PHOLM_MW", _p(b))
            for ln, a, b in zip(body, adj_rho, adj_mw)]

    tex = _table(
        r"Dependence of the \emph{intrinsic} (i.e.\ $\gamma_{\mathrm{tile}}=1$) "
        r"HPX/UTM ratio on the tile's position on the sphere. HEALPix cells are "
        r"equal-area but not equal-shape, and their local elongation differs between "
        r"the equatorial belt ($|\mathrm{lat}| < 41.81^\circ$) and the polar caps. "
        r"Spearman $\rho$ is against $|\mathrm{latitude}|$; the last four columns "
        r"compare the two regimes (two-sided Mann--Whitney).",
        "tab:latitude_effect", "lrrrrrrrr",
        r"Target & $n$ & $\rho$ & $p_\rho$ & $p_\rho^{\mathrm{Holm}}$ & "
        r"Median (belt) & Median (cap) & $p$ & $p^{\mathrm{Holm}}$ \\",
        body, wide=True)
    if not body:
        print("  latitude table skipped: no per-region latitudes "
              "(run build_tiles.py so data/tiles/tiles_index.csv exists)")
        return "", stats_out
    return tex, stats_out


def table_fixed_codec(hpx, utm, scenes, rmse_col, codecs):
    """Same codec on both grids -- isolates the grid effect from codec selection."""
    body, out = [], {}
    for c in codecs:
        cells = []
        for target in cfg.TARGET_RMSE:
            ratios = []
            for s in scenes:
                h = st.interp_size_at_rmse(hpx[(hpx.scene == s) & (hpx.codec == c)],
                                           target, rmse_col)
                u = st.interp_size_at_rmse(utm[(utm.scene == s) & (utm.codec == c)],
                                           target)
                if np.isfinite(h) and np.isfinite(u) and u > 0:
                    ratios.append(h / (u * cfg.GAMMA_UTM))
            t = st.paired_ratio_test(ratios)
            out.setdefault(c, {})[f"{target:g}"] = t
            cells.append(f"{_fmt(t['median'])}" if t["n"] >= cfg.MIN_SCENES else "--")
        body.append(f"\\texttt{{{c.replace('_', chr(92)+'_')}}} & " + " & ".join(cells)
                    + r" \\")
    header = "Codec & " + " & ".join(RMSE_LABEL.get(t, f"${t:g}$")
                                     for t in cfg.TARGET_RMSE) + r" \\"
    tex = _table(
        r"Same-codec control: median per-region HPX/UTM compressed-size ratio at each "
        r"target RMSE when the \emph{same} codec is used on both grids (UTM includes "
        r"$\gamma_{\mathrm{tile}}$). Complements the oracle-best comparison of "
        r"Table~\ref{tab:storage_at_rmse}, which lets codec selection differ between "
        r"grids.",
        "tab:fixed_codec", "l" + "r" * len(cfg.TARGET_RMSE), header, body)
    return tex, out


def table_representation_floor(fl):
    """Table XIII: the transfer-and-read floor by scene class.

    This table had no producing script -- it was written once by hand from
    common_reference_floor.csv and then drifted from the wording the text
    settled on. It is now emitted from the same CSV the floors are read from,
    and it names the quantity the way Section V-F established it: the floor of
    the transfer-and-read pair, an upper bound on what the carrier costs under
    the present reading operator, not the grid's own representation error.
    """
    if fl is None or not len(fl):
        return ""
    fl = fl.copy()
    fl["pct"] = 100.0 * fl.representation_floor_rmse / fl.utm_std_on_mask
    mode = sorted(set(fl.get("mode", pd.Series(["nearest"])).astype(str)))
    body = []
    for cls, g in fl.groupby("scene_class"):
        body.append(f"{cls} & {len(g)} & {g.representation_floor_rmse.median():.1e} & "
                    f"[{g.representation_floor_rmse.min():.1e}, "
                    f"{g.representation_floor_rmse.max():.1e}] & "
                    f"{g.pct.median():.1f} \\\\")
    body.append("\\midrule")
    body.append(f"all & {len(fl)} & {fl.representation_floor_rmse.median():.1e} & "
                f"[{fl.representation_floor_rmse.min():.1e}, "
                f"{fl.representation_floor_rmse.max():.1e}] & "
                f"{fl.pct.median():.1f} \\\\")
    caption = (
        "Transfer-and-read floor: the RMSE of the \\emph{uncompressed} level-20 "
        "carrier tile, read on the native UTM pixels of the two footprints' "
        "intersection with the " + "/".join(mode) + " read, against the "
        "Sentinel-2 observation itself. No codec setting on the HEALPix side "
        "can go below it, so any target tighter than it is unreachable by "
        "construction under the common reference. It is the error of the "
        "transfer and the read together, and Section~\\ref{sec:operator} "
        "shows that the read alone accounts for an error of this size, so it "
        "is an upper bound on what the carrier representation costs under the "
        "present reading operator, not the grid's own representation error. "
        "It scales with the scene's high-frequency content.")
    return _table(caption, "tab:representation_floor", "lrrrr",
                  "Class & $n$ & Median floor & Range & \\% of scene std \\\\",
                  body)


def table_regions(index):
    body = []
    if "scene_class" not in index:
        return "", {}
    for cls, g in index.groupby("scene_class"):
        regimes = g.hpx_regime.value_counts().to_dict() if "hpx_regime" in g else {}
        body.append(f"{cls} & {len(g):d} & {g.lat.min():.1f} & {g.lat.max():.1f} & "
                    f"{regimes.get('equatorial', 0)} / {regimes.get('cap', 0)} \\\\")
    tex = _table(
        f"Evaluation regions, shared with the companion \\texttt{{healpix-resample}} "
        f"study: one Sentinel-2 L2A B04 tile per region, {cfg.WINDOW_PX}$\\times$"
        f"{cfg.WINDOW_PX} native pixels re-extracted at the same anchors and from the "
        f"same products as its published 256$\\times$256 patches.",
        "tab:benchmark_scenes", "lrrrl",
        r"Class & Regions & Min lat & Max lat & Belt / cap \\", body)
    return tex, {c: int(n) for c, n in index.scene_class.value_counts().items()}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metric", choices=["self", "common"], default="self")
    ap.add_argument("--no-write-tex", action="store_true")
    ap.add_argument("--pareto", action="store_true",
                    help="reduce each sweep to its Pareto frontier before "
                         "interpolating (closes the 'dominated point' open item)")
    ap.add_argument("--allow-asymmetric-codecs", action="store_true",
                    help="produce numbers even if a codec ran on one grid only "
                         "(diagnostic use; never for the paper)")
    args = ap.parse_args(argv)

    st.USE_PARETO = bool(args.pareto)
    hpx, utm, scenes, rmse_col = load(args.metric)
    codecs = [c for c in cfg.COMPRESSOR_CODECS
              if c in set(hpx.codec) or c in set(utm.codec)]
    print(f"{len(scenes)} regions | codecs {codecs} | HPX metric column {rmse_col!r}")
    envf = cfg.RESULTS_DIR / "bench_environment.json"
    if envf.exists():
        e = json.loads(envf.read_text())
        print(f"  results written by code {e.get('code_sha', '?')} at "
              f"{e.get('written_utc', '?')}")
    floorf = cfg.RESULTS_DIR / "common_reference_floor.csv"
    if args.metric == "common" and floorf.exists():
        fl = pd.read_csv(floorf)
        modes = sorted(set(fl.get("mode", pd.Series(dtype=str)).astype(str)))
        print(f"  common-reference floors: {len(fl)} regions, mode(s) {modes}, "
              f"median {fl.representation_floor_rmse.median():.3e}")
        if "support_median" not in fl.columns:
            raise SystemExit(
                "common_reference_floor.csv predates the projection fix (no "
                "support_median column): these floors were produced by the "
                "all-zero projection and every rmse_common is the RMS of the "
                "observations, not a codec error. Rerun:\n"
                "    python run_codec_bench.py --grid both --common-reference --force")

    asym = check_codec_symmetry(hpx, utm, scenes, codecs)
    if asym:
        print("\nCODEC COVERAGE IS NOT SYMMETRIC -- the oracle comparison would "
              "be decided partly by which codec ran, not by which grid is "
              "better:")
        for scene, only_utm, only_hpx in asym[:10]:
            bits = []
            if only_utm:
                bits.append(f"missing on HPX: {', '.join(only_utm)}")
            if only_hpx:
                bits.append(f"missing on UTM: {', '.join(only_hpx)}")
            print(f"  {scene}: " + "; ".join(bits))
        if len(asym) > 10:
            print(f"  ... and {len(asym)-10} more region(s)")
        print("  Look at results/bench_parts/<grid>__<region>.skips.json for the "
              "reason, fix it, and rerun that grid with --force.")
        if not args.allow_asymmetric_codecs:
            raise SystemExit("refusing to produce paper numbers from an "
                             "asymmetric comparison (--allow-asymmetric-codecs "
                             "to override, for a diagnostic run only)")

    floors = {"hpx": quantization_floor(hpx, scenes),
              "utm": quantization_floor(utm, scenes)}
    for grid, fl in floors.items():
        if fl:
            print(f"16-bit requantization floor, {grid}: median "
                  f"{fl['median']:.2e} (range {fl['min']:.2e}-{fl['max']:.2e}). "
                  f"Targets at or below it are not achievable by any codec on "
                  f"that grid.")

    # Per-region latitude and class: the tile index carries them, but they are a
    # property of the sites, so fall back to the site table when a checkout has
    # the results without the intermediate tiles.
    index_path = cfg.TILE_DIR / "tiles_index.csv"
    index = pd.read_csv(index_path) if index_path.exists() else pd.DataFrame()
    if "abs_lat" not in index.columns:
        try:
            import region_sites
            sites = region_sites.load_sites()[
                ["scene", "scene_class", "region_id", "lat", "lon", "abs_lat",
                 "hpx_regime"]]
            index = (sites if index.empty
                     else index.merge(sites, on="scene", how="left"))
            print(f"  latitudes taken from the site table "
                  f"({len(index)} regions); tiles_index.csv not needed")
        except Exception as exc:                          # noqa: BLE001
            print(f"  no per-region latitude available: {type(exc).__name__}: {exc}")

    tex_iii, acc = table_storage_accounting(utm, scenes)
    rows, family = _ratio_rows(hpx, utm, scenes, rmse_col, codecs)
    tex_iv, tex_v = table_iv_v(rows, family, len(scenes))
    if args.metric != "self":
        suffix_lbl = "_" + args.metric
        for _i, _t in enumerate([tex_iv, tex_v]):
            pass
        tex_iv = tex_iv.replace("tab:storage_at_rmse", "tab:storage_at_rmse" + suffix_lbl)
        tex_v = tex_v.replace("tab:error_at_size", "tab:error_at_size" + suffix_lbl)
    tex_gs, gammas = table_gamma_star(rows)
    tex_vi, ranking = table_codec_ranking(hpx, scenes, rmse_col, codecs)
    # Two tables are built from the tile index rather than from the sweep, so a
    # partial build_tiles run used to make them disagree with every other table
    # in the same invocation -- n=1 beside n=40, published. Refuse rather than
    # emit a table that contradicts the run it belongs to.
    n_index = int(len(index))
    if n_index and n_index != len(scenes):
        raise SystemExit(
            f"the tile index describes {n_index} region(s) but the sweep covers "
            f"{len(scenes)}: the region and latitude tables would contradict "
            f"every other table in this run. Rebuild the index over all regions "
            f"(python build_tiles.py) or delete "
            f"{cfg.TILE_DIR / 'tiles_index.csv'} to fall back on the site table, "
            f"then re-run.")
    tex_lat, lat_stats = table_latitude(rows, index)
    tex_fix, fixed = table_fixed_codec(hpx, utm, scenes, rmse_col, codecs)
    tex_reg, per_class = table_regions(index)
    tex_floor = ""
    if args.metric == "common" and floorf.exists():
        tex_floor = table_representation_floor(pd.read_csv(floorf))

    if not args.no_write_tex:
        cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
        # The two metrics answer different questions and both are reported, so
        # they must not overwrite each other's tables.
        tag = "" if args.metric == "self" else f"_{args.metric}"
        for name, tex in [("table_storage_accounting", tex_iii),
                          ("table_storage_at_constant_rmse", tex_iv),
                          ("table_error_at_constant_size", tex_v),
                          ("table_gamma_star", tex_gs),
                          ("table_codec_ranking", tex_vi),
                          ("table_latitude_effect", tex_lat),
                          ("table_fixed_codec", tex_fix),
                          ("table_benchmark_regions", tex_reg)]:
            if tex:
                (cfg.TEX_TABLES / f"{name}{tag}.tex").write_text(tex, encoding="utf-8")
        if tex_floor:
            # one floor, one file: the common-reference run owns it
            (cfg.TEX_TABLES / "table_representation_floor.tex").write_text(
                tex_floor, encoding="utf-8")
        print(f"LaTeX tables -> {cfg.TEX_TABLES}")

    payload = dict(
        metric=args.metric, n_regions=len(scenes), regions=scenes, codecs=codecs,
        family_size=family, gamma_utm=cfg.GAMMA_UTM, shared_raw=cfg.SHARED_RAW,
        storage_accounting=acc, gamma_star=gammas, codec_ranking=ranking,
        quantization_floor=floors, codec_asymmetry=[
            dict(scene=a[0], missing_on_hpx=a[1], missing_on_utm=a[2])
            for a in asym],
        latitude=lat_stats, fixed_codec=fixed,
        targets=[dict(kind=r["kind"], target=r["target"], test=r["test"],
                      test_gamma1=r["test_gamma1"], p_holm=r["p_holm"],
                      bracketed_fraction=r["bracketed_fraction"],
                      ratios={k: float(v) for k, v in r["ratios"].items()})
                 for r in rows],
    )
    payload["pareto"] = bool(args.pareto)
    suffix = args.metric + ("_pareto" if args.pareto else "")
    out_json = cfg.RESULTS_DIR / f"paper_numbers_{suffix}.json"
    out_json.write_text(json.dumps(payload, indent=2, default=float), encoding="utf-8")

    # --- digest ready to paste into the abstract / results text --------------
    md = [f"# Paper numbers ({args.metric} RMSE reference"
          + (", Pareto frontier" if args.pareto else "") + ")", "",
          f"- regions: **{len(scenes)}**, codecs: {', '.join(codecs)}",
          f"- raw ratio HPX / (UTM x gamma={cfg.GAMMA_UTM}): "
          f"**{acc['ratio_gamma']['median']:.3f}** "
          f"(gamma=1: {acc['ratio_gamma1']['median']:.3f})",
          f"- multiple-testing family: {family} targets",
          (f"- 16-bit requantization floor: HPX "
           f"{floors['hpx'].get('median', float('nan')):.2e}, UTM "
           f"{floors['utm'].get('median', float('nan')):.2e} "
           f"(the image codecs are fed a 16-bit tile, so their `lossless` "
           f"points are lossless after that step, not bit-exact)"), "",
          "## Ratios (oracle codec, per region)", "",
          "| protocol | target | n | median | range | p | p_Holm |",
          "|---|---|---|---|---|---|---|"]
    for r in rows:
        t = r["test"]
        if not r.get("reported"):
            continue
        md.append(f"| {r['kind']} | {r['target']:g} | {t['n']} | {t['median']:.3f} | "
                  f"[{t['lo']:.3f}, {t['hi']:.3f}] | {t['p']:.4f} | {r['p_holm']:.4f} |")
    md += ["", "## Censoring (fraction of regions where the target is bracketed",
           "on both grids; the rest are conservative bounds, not estimates)", ""]
    for r in rows:
        if r["test"]["n"] < cfg.MIN_SCENES:
            continue
        md.append(f"- {r['kind']} {r['target']:g}: "
                  f"{100*r['bracketed_fraction']:.0f}% bracketed")
    md += ["", "## Intrinsic gap (gamma=1) and critical overhead", ""]
    for r in rows:
        if r["kind"] != "rmse" or not r.get("reported"):
            continue
        g = r["test_gamma1"]["median"]
        md.append(f"- RMSE {r['target']:g}: gap {100*(g-1):+.0f}%, gamma* = {g:.2f}")
    md += ["", "## HEALPix regime (latitude)", ""]
    for k, v in lat_stats.items():
        md.append(f"- {k}: Spearman rho={v['spearman']['rho']:.3f} "
                  f"(p={v['spearman']['p']:.3f}, "
                  f"Holm={v['spearman'].get('p_holm', float('nan')):.3f}), belt "
                  f"{v['equatorial_vs_cap']['median_a']:.3f} vs cap "
                  f"{v['equatorial_vs_cap']['median_b']:.3f} "
                  f"(p={v['equatorial_vs_cap']['p']:.3f}, "
                  f"Holm={v['equatorial_vs_cap'].get('p_holm', float('nan')):.3f})")
    md += ["", "## Codec ranking (HEALPix side)", "",
           f"- {ranking['points']} operating points, wins sum to the same number"]
    for c, w in sorted(ranking["wins"].items(), key=lambda kv: -kv[1]):
        md.append(f"- {c}: {w} wins / {ranking['evaluable'].get(c,0)} evaluable, "
                  f"mean rank {ranking['mean_rank'].get(c, float('nan')):.2f}")
    floor_path = cfg.RESULTS_DIR / "common_reference_floor.csv"
    if floor_path.exists():
        fl = pd.read_csv(floor_path)
        md += ["", "## HEALPix representation floor (common reference)", "",
               f"- median {fl.representation_floor_rmse.median():.3e}, "
               f"range [{fl.representation_floor_rmse.min():.3e}, "
               f"{fl.representation_floor_rmse.max():.3e}] over {len(fl)} regions",
               "- every target RMSE below this floor is unreachable on the HEALPix "
               "side by construction and must be reported as such."]
    out_md = cfg.RESULTS_DIR / f"paper_numbers_{suffix}.md"
    out_md.write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"digest -> {out_md}\njson   -> {out_json}")
    print("\n".join(md[:14]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
