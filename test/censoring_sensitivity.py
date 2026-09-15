"""How much of the matched-error verdict rests on censored sweep points?

The rate--distortion protocol (paper_stats.FEASIBLE_CENSORING) returns, for a
codec whose coarsest setting already meets a target, the size of that coarsest
point. That value is admissible but it is a *bound*: a coarser setting, had it
been swept, could only have been smaller. The tables of the paper therefore mix
interpolated estimates with bounds, and a signed-rank test on such a mixture
is only as good as the share of bounds in it. This script measures that share
and reports what survives when the bounds are treated as what they are.

Three readings are produced, per target, from the same sweeps that produced
the paper's tables (nothing is re-measured):

  as reported     the ratios exactly as update_paper_numbers.py computes them,
                  for reference;
  bracketed only  the same test restricted to the regions whose selected
                  (oracle-best) codec brackets the target on BOTH grids, so that
                  every ratio in the test is an interpolated estimate;
  bound interval  for every region the reported ratio is classified as exact,
                  an upper bound (the HEALPix side is censored, the true ratio
                  can only be lower), a lower bound (the UTM side is censored)
                  or undetermined (both). The median is then bracketed by
                  replacing lower-bound ratios by +inf (upper bound on the
                  median) or upper-bound ratios by 0 (lower bound on it), and a
                  sign test is run with every undetermined region assigned to
                  whichever side is least favourable to the reported verdict.

A verdict is called robust to censoring when the bracketed-only test agrees
in sign and survives Holm, and the median interval lies entirely on one side
of parity.

The same-codec table is treated the same way, codec by codec, because there a
single censored codec is not rescued by the others.

Sources: hpx_codec_bench.csv, utm_codec_bench.csv (within-representation
metric, as Tables VII--VIII and the same-codec table).

Usage
-----
    python censoring_sensitivity.py            # tables + JSON digest
    python censoring_sensitivity.py --gamma 1  # ratios with no overhead
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

import paper_config as cfg
import paper_stats as st

RMSE_LABEL = {1e-4: r"$10^{-4}$", 3e-4: r"$3\times10^{-4}$", 1e-3: r"$10^{-3}$",
              3e-3: r"$3\times10^{-3}$", 1e-2: r"$10^{-2}$"}


# --- per-region classification ----------------------------------------------
def _best_with_status(df, scene, codecs, target, kind, rmse_col="rmse"):
    """(value, selected codec, selected codec brackets?, all feasible bracket?)."""
    sub = df[df.scene == scene]
    if kind == "rmse":
        val_fn, br_fn = st.interp_size_at_rmse, st.bracketed_size_at_rmse
    else:
        val_fn, br_fn = st.interp_rmse_at_size, st.bracketed_rmse_at_size
    best, name, best_br, all_br = np.inf, None, False, True
    for c in codecs:
        d = sub[sub.codec == c]
        v = val_fn(d, target, rmse_col)
        if not np.isfinite(v):
            continue
        br = br_fn(d, target, rmse_col)
        all_br &= br
        if v < best:
            best, name, best_br = v, c, br
    if name is None:
        return np.nan, None, False, False
    return float(best), name, bool(best_br), bool(all_br)


def classify(hpx, utm, scenes, codecs, rmse_col="rmse", gamma=cfg.GAMMA_UTM):
    """One row per (target, region) with the ratio and its censoring status."""
    rows = []
    targets = [("rmse", t, t) for t in cfg.TARGET_RMSE] + \
              [("budget", f, f * cfg.SHARED_RAW) for f in cfg.BUDGETS]
    for kind, label, tgt in targets:
        for s in scenes:
            if kind == "rmse":
                h, hc, hb, hall = _best_with_status(hpx, s, codecs, tgt, kind, rmse_col)
                u, uc, ub, uall = _best_with_status(utm, s, codecs, tgt, kind)
                u_scaled = u * gamma
            else:
                h, hc, hb, hall = _best_with_status(hpx, s, codecs, tgt, kind, rmse_col)
                u, uc, ub, uall = _best_with_status(utm, s, codecs, tgt / gamma, kind)
                u_scaled = u
            if not (np.isfinite(h) and np.isfinite(u_scaled) and u_scaled > 0):
                continue
            r = h / u_scaled
            # censored HPX value is an upper bound on HPX -> ratio is an upper bound
            if hb and ub:
                status = "exact"
            elif (not hb) and ub:
                status = "upper"
            elif hb and (not ub):
                status = "lower"
            else:
                status = "undetermined"
            rows.append(dict(kind=kind, target=label, scene=s,
                             cls=s.split("__")[1], ratio=r, status=status,
                             hpx_codec=hc, utm_codec=uc,
                             strict=bool(hall and uall)))
    return pd.DataFrame(rows)


# --- per-target summaries ----------------------------------------------------
def _median_interval(df):
    r = df.ratio.to_numpy(dtype=float)
    hi = r.copy(); hi[df.status.isin(["lower", "undetermined"]).to_numpy()] = np.inf
    lo = r.copy(); lo[df.status.isin(["upper", "undetermined"]).to_numpy()] = 0.0
    return float(np.median(lo)), float(np.median(hi))


def _robust_signs(df):
    """Regions whose side of parity survives any interpretation of the bound."""
    r, s = df.ratio.to_numpy(dtype=float), df.status.to_numpy()
    sure_hpx = int(((r < 1) & np.isin(s, ["exact", "upper"])).sum())
    sure_utm = int(((r > 1) & np.isin(s, ["exact", "lower"])).sum())
    undet = int(len(r) - sure_hpx - sure_utm)
    return sure_hpx, sure_utm, undet


def _sign_p(k, n):
    from scipy import stats
    return float(stats.binomtest(k, n, 0.5).pvalue) if n else np.nan


def summarize(cl):
    out = []
    for (kind, target), g in cl.groupby(["kind", "target"], sort=False):
        rep = st.paired_ratio_test(g.ratio)
        ex = g[g.status == "exact"]
        ex_t = st.paired_ratio_test(ex.ratio)
        strict = g[g.strict]
        st_t = st.paired_ratio_test(strict.ratio)
        lo, hi = _median_interval(g)
        sure_hpx, sure_utm, undet = _robust_signs(g)
        n = len(g)
        # least favourable to the reported verdict
        verdict_hpx = rep["median"] < 1
        if verdict_hpx:
            p_worst = _sign_p(sure_hpx, n)            # undetermined -> UTM
        else:
            p_worst = _sign_p(sure_utm, n)            # undetermined -> HPX
        out.append(dict(
            kind=kind, target=target, n=n,
            n_exact=int(len(ex)), n_strict=int(len(strict)),
            n_upper=int((g.status == "upper").sum()),
            n_lower=int((g.status == "lower").sum()),
            n_undetermined=int((g.status == "undetermined").sum()),
            median_reported=rep["median"], p_reported=rep["p"],
            median_exact=ex_t["median"], p_exact=ex_t["p"],
            median_strict=st_t["median"], p_strict=st_t["p"],
            median_lo=lo, median_hi=hi,
            sure_hpx=sure_hpx, sure_utm=sure_utm, undetermined=undet,
            p_sign_worst=p_worst,
            robust=bool(np.isfinite(ex_t["median"])
                        and (ex_t["median"] < 1) == verdict_hpx
                        and ((hi < 1) if verdict_hpx else (lo > 1)))))
    df = pd.DataFrame(out)
    # Holm within the same family the paper corrects: the targets it reports
    # (update_paper_numbers withholds a target that fewer than MIN_BRACKETED of
    # the regions bracket; the flag is read back from its digest so that the
    # family is the same here).
    df["reported"] = True
    digest = cfg.RESULTS_DIR / "paper_numbers_self.json"
    if digest.exists():
        doc = json.loads(digest.read_text())
        rep = {(t["kind"], float(t["target"])):
               bool(t.get("bracketed_fraction", 1.0) >= cfg.MIN_BRACKETED
                    and t["test"]["n"] >= cfg.MIN_SCENES)
               for t in doc["targets"]}
        df["reported"] = [rep.get((k, float(t)), True)
                          for k, t in zip(df.kind, df.target)]
    for col in ("p_reported", "p_exact", "p_strict"):
        adj = np.full(len(df), np.nan)
        m = df.reported.to_numpy()
        adj[m] = st.holm(df.loc[m, col].to_numpy())
        df[col + "_holm"] = adj
    return df


# --- same-codec table --------------------------------------------------------
def fixed_codec(hpx, utm, scenes, codecs, rmse_col="rmse", gamma=cfg.GAMMA_UTM):
    rows = []
    for c in codecs:
        for target in cfg.TARGET_RMSE:
            recs = []
            for s in scenes:
                dh = hpx[(hpx.scene == s) & (hpx.codec == c)]
                du = utm[(utm.scene == s) & (utm.codec == c)]
                h = st.interp_size_at_rmse(dh, target, rmse_col)
                u = st.interp_size_at_rmse(du, target)
                if not (np.isfinite(h) and np.isfinite(u) and u > 0):
                    continue
                hb, ub = (st.bracketed_size_at_rmse(dh, target, rmse_col),
                          st.bracketed_size_at_rmse(du, target))
                recs.append((h / (u * gamma), hb, ub))
            if not recs:
                continue
            r = np.array([x[0] for x in recs])
            both = np.array([x[1] and x[2] for x in recs])
            hpx_c = np.array([not x[1] for x in recs])
            utm_c = np.array([not x[2] for x in recs])
            t_all, t_br = st.paired_ratio_test(r), st.paired_ratio_test(r[both])
            rows.append(dict(codec=c, target=target, n=len(r),
                             n_bracketed=int(both.sum()),
                             frac_hpx_censored=float(hpx_c.mean()),
                             frac_utm_censored=float(utm_c.mean()),
                             median_reported=t_all["median"],
                             median_bracketed=t_br["median"], p_bracketed=t_br["p"]))
    return pd.DataFrame(rows)


# --- LaTeX -------------------------------------------------------------------
def _f(x, nd=2):
    if not np.isfinite(x):
        return "--"
    if np.isinf(x):
        return r"$\infty$"
    return f"{x:.{nd}f}"


def _p(x):
    if not np.isfinite(x):
        return "--"
    if x < 1e-4:
        return "$<10^{-4}$"
    if x >= 1e-3:
        return f"{x:.3f}"
    mant, exp = f"{x:.0e}".split("e")
    return f"${mant}\\times10^{{{int(exp)}}}$"


def _tlabel(kind, t):
    return RMSE_LABEL.get(t, f"${t:g}$") if kind == "rmse" else f"{100*t:g}\\%"


def latex_table(summary, gamma):
    body = []
    for _, r in summary.iterrows():
        if r.n < cfg.MIN_SCENES or not r.reported:
            continue
        body.append(
            f"{'RMSE' if r.kind == 'rmse' else 'Budget'} {_tlabel(r.kind, r.target)} & "
            f"{r.n_exact:d} / {r.n:d} & {_f(r.median_reported)} & "
            f"{_f(r.median_exact)} & {_p(r.p_exact_holm)} & "
            f"[{_f(r.median_lo)}, {_f(r.median_hi)}] & "
            f"{r.sure_hpx:d} / {r.sure_utm:d} / {r.undetermined:d} & "
            f"{_p(r.p_sign_worst)} \\\\")
    note = (
        "Per target: number of regions whose oracle-best codec brackets the "
        "target with two measured points on both grids (against all regions); "
        "the median HPX/UTM ratio as reported in "
        "Tables~\\ref{tab:storage_at_rmse}--\\ref{tab:error_at_size}; the "
        "same median and its Holm-adjusted signed-rank $p$ on the bracketed "
        "regions only; the interval that bounds the median of all regions "
        "once every censored value is taken as the bound it is (a censored "
        "HEALPix value can only be smaller, a censored UTM value only larger); "
        "the number of regions whose side of parity holds under any "
        "interpretation of the bounds, for HEALPix, for UTM, and undetermined; "
        "and the sign-test $p$ with every undetermined region assigned to the "
        f"side least favourable to the reported verdict. The UTM side carries "
        f"$\\gamma_{{\\mathrm{{tile}}}}={gamma:g}$ as in those tables. RMSE "
        "targets compare sizes, budget targets compare errors.")
    return "\n".join([
        r"\begin{table*}[t]", r"\centering",
        *cfg.tex_caption("Sensitivity of the matched comparisons to censored sweep points",
                         "tab:censoring"),
        r"\footnotesize",
        r"\begin{tabular}{@{}lrrrrlrr@{}}", r"\toprule",
        r"Target & Bracketed & \multicolumn{3}{c}{Median ratio} & Median & Sure & Sign $p$ \\",
        r"\cmidrule(lr){3-5}",
        r" & regions & reported & bracketed & $p^{\mathrm{Holm}}$ & interval & HPX / UTM / ? & (worst case) \\",
        r"\midrule", *body, r"\bottomrule", r"\end{tabular}",
        *cfg.tex_note(note, wide=True), r"\end{table*}"]) + "\n"


def latex_fixed(fc, gamma):
    body = []
    for c in cfg.COMPRESSOR_CODECS:
        g = fc[fc.codec == c]
        if g.empty:
            continue
        cells = []
        for t in cfg.TARGET_RMSE:
            r = g[np.isclose(g.target, t)]
            if r.empty:
                cells.append("--")
                continue
            r = r.iloc[0]
            cells.append(f"{_f(r.median_bracketed)} ({r.n_bracketed:d})")
        body.append(f"\\texttt{{{c.replace('_', chr(92) + '_')}}} & " + " & ".join(cells) + r" \\")
    note = (
        "The same-codec ratios of Table~\\ref{tab:fixed_codec} recomputed on "
        "the regions where that codec brackets the target on both grids, with "
        "the number of such regions in parentheses (out of 40). A cell whose "
        "count is small is a statement about a few scenes; a cell of "
        "Table~\\ref{tab:fixed_codec} whose count here is small is a set of "
        f"bounds rather than an estimate. UTM carries $\\gamma_{{\\mathrm{{tile}}}}={gamma:g}$.")
    header = "Codec & " + " & ".join(RMSE_LABEL[t] for t in cfg.TARGET_RMSE) + r" \\"
    return "\n".join([
        r"\begin{table}[t]", r"\centering",
        *cfg.tex_caption("Same-codec ratios on bracketed regions only",
                         "tab:fixed_codec_bracketed"),
        r"\footnotesize",
        r"\begin{tabular}{@{}l" + "r" * len(cfg.TARGET_RMSE) + "@{}}", r"\toprule",
        header, r"\midrule", *body, r"\bottomrule", r"\end{tabular}",
        *cfg.tex_note(note), r"\end{table}"]) + "\n"


# --- driver ------------------------------------------------------------------
def load():
    hpx = pd.read_csv(cfg.RESULTS_DIR / "hpx_codec_bench.csv")
    utm = pd.read_csv(cfg.RESULTS_DIR / "utm_codec_bench.csv")
    scenes = sorted(set(hpx.scene) & set(utm.scene))
    codecs = [c for c in cfg.COMPRESSOR_CODECS
              if c in set(hpx.codec) or c in set(utm.codec)]
    return hpx, utm, scenes, codecs


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gamma", type=float, default=cfg.GAMMA_UTM)
    args = ap.parse_args(argv)
    hpx, utm, scenes, codecs = load()
    print(f"{len(scenes)} regions, codecs {codecs}, gamma={args.gamma:g}")

    cl = classify(hpx, utm, scenes, codecs, gamma=args.gamma)
    summary = summarize(cl)
    fc = fixed_codec(hpx, utm, scenes, codecs, gamma=args.gamma)

    tag = "" if np.isclose(args.gamma, cfg.GAMMA_UTM) else f"_gamma{args.gamma:g}"
    cl.to_csv(cfg.RESULTS_DIR / f"censoring_regions{tag}.csv", index=False)
    summary.to_csv(cfg.RESULTS_DIR / f"censoring_summary{tag}.csv", index=False)
    fc.to_csv(cfg.RESULTS_DIR / f"censoring_fixed_codec{tag}.csv", index=False)
    if not tag:
        cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
        (cfg.TEX_TABLES / "table_censoring.tex").write_text(
            latex_table(summary, args.gamma), encoding="utf-8")
        (cfg.TEX_TABLES / "table_fixed_codec_bracketed.tex").write_text(
            latex_fixed(fc, args.gamma), encoding="utf-8")

    digest = {}
    for _, r in summary.iterrows():
        digest[f"{r.kind}:{r.target:g}"] = {
            k: (None if (isinstance(v, float) and not np.isfinite(v)) else
                (float(v) if isinstance(v, (float, np.floating)) else
                 int(v) if isinstance(v, (np.integer,)) else v))
            for k, v in r.items() if k not in ("kind", "target")}
    (cfg.RESULTS_DIR / f"censoring_digest{tag}.json").write_text(
        json.dumps(digest, indent=1, default=str))

    pd.set_option("display.width", 200)
    cols = ["kind", "target", "n", "n_exact", "n_strict", "n_upper", "n_lower",
            "n_undetermined", "median_reported", "p_reported_holm",
            "median_exact", "p_exact_holm", "median_lo", "median_hi",
            "sure_hpx", "sure_utm", "undetermined", "p_sign_worst", "robust"]
    print(summary[cols].to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print()
    print(fc.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
