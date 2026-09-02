"""Statistics helpers shared by the paper-number scripts.

Extracted from stats_significance.ipynb so the interpolation, the oracle-codec
selection and the multiple-testing correction have exactly one implementation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


# --- rate-distortion interpolation ------------------------------------------
def _xy(df, xcol, ycol):
    """Sorted, duplicate-averaged (x, y) from a raw sweep. NOT a Pareto front."""
    d = df[[xcol, ycol]].replace([np.inf, -np.inf], np.nan).dropna()
    d = d[(d[xcol] > 0) & (d[ycol] > 0)]
    if len(d) < 2:
        return None
    d = d.groupby(xcol, as_index=False)[ycol].mean().sort_values(xcol)
    return d[xcol].to_numpy(), d[ycol].to_numpy()


def interp_size_at_rmse(df, target, rmse_col="rmse"):
    xy = _xy(df, rmse_col, "total_bytes")
    if xy is None:
        return np.nan
    x, y = xy
    if not (x.min() <= target <= x.max()):
        return np.nan
    return float(np.exp(np.interp(np.log(target), np.log(x), np.log(y))))


def interp_rmse_at_size(df, target_bytes, rmse_col="rmse"):
    xy = _xy(df, "total_bytes", rmse_col)
    if xy is None:
        return np.nan
    x, y = xy
    if not (x.min() <= target_bytes <= x.max()):
        return np.nan
    return float(np.exp(np.interp(np.log(target_bytes), np.log(x), np.log(y))))


def best_of_codecs_size(df, scene, codecs, target, gamma=1.0, rmse_col="rmse"):
    sub = df[df.scene == scene]
    vals = [interp_size_at_rmse(sub[sub.codec == c], target, rmse_col) for c in codecs]
    vals = [v * gamma for v in vals if np.isfinite(v)]
    return min(vals) if vals else np.nan


def best_of_codecs_rmse(df, scene, codecs, target_bytes, rmse_col="rmse"):
    sub = df[df.scene == scene]
    vals = [interp_rmse_at_size(sub[sub.codec == c], target_bytes, rmse_col)
            for c in codecs]
    vals = [v for v in vals if np.isfinite(v)]
    return min(vals) if vals else np.nan


def argbest_codec_size(df, scene, codecs, target, rmse_col="rmse"):
    sub = df[df.scene == scene]
    best, name = np.inf, None
    for c in codecs:
        v = interp_size_at_rmse(sub[sub.codec == c], target, rmse_col)
        if np.isfinite(v) and v < best:
            best, name = v, c
    return name, (best if np.isfinite(best) else np.nan)


def argbest_codec_rmse(df, scene, codecs, target_bytes, rmse_col="rmse"):
    sub = df[df.scene == scene]
    best, name = np.inf, None
    for c in codecs:
        v = interp_rmse_at_size(sub[sub.codec == c], target_bytes, rmse_col)
        if np.isfinite(v) and v < best:
            best, name = v, c
    return name, (best if np.isfinite(best) else np.nan)


# --- inference ---------------------------------------------------------------
def paired_ratio_test(ratios):
    """Two-sided Wilcoxon signed-rank that a set of ratios differs from 1.

    Ranks are taken on log(ratio): a ratio and its inverse then carry the same
    magnitude, which a test on (ratio - 1) does not give. `p_linear` keeps the
    previous convention for continuity.
    """
    r = np.asarray([x for x in ratios if np.isfinite(x) and x > 0], dtype=float)
    out = dict(n=int(r.size), median=np.nan, lo=np.nan, hi=np.nan,
               p=np.nan, p_linear=np.nan, ci_lo=np.nan, ci_hi=np.nan)
    if r.size == 0:
        return out
    out.update(median=float(np.median(r)), lo=float(r.min()), hi=float(r.max()))
    if r.size >= 5:
        method = "exact" if r.size <= 25 else "auto"
        try:
            out["p"] = float(stats.wilcoxon(np.log(r), alternative="two-sided",
                                            method=method).pvalue)
            out["p_linear"] = float(stats.wilcoxon(r - 1.0, alternative="two-sided",
                                                   method=method).pvalue)
        except ValueError:
            pass
    rng = np.random.default_rng(0)
    boot = np.median(rng.choice(r, size=(10000, r.size), replace=True), axis=1)
    out["ci_lo"], out["ci_hi"] = (float(np.percentile(boot, 2.5)),
                                  float(np.percentile(boot, 97.5)))
    return out


def holm(pvalues):
    """Holm-Bonferroni adjusted p-values, NaN-safe and order-preserving."""
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan)
    idx = np.where(np.isfinite(p))[0]
    if idx.size == 0:
        return out
    m = idx.size
    order = idx[np.argsort(p[idx])]
    running = 0.0
    for rank, i in enumerate(order):
        adj = min(1.0, (m - rank) * p[i])
        running = max(running, adj)
        out[i] = running
    return out


def spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5:
        return dict(n=int(m.sum()), rho=np.nan, p=np.nan)
    r = stats.spearmanr(x[m], y[m])
    return dict(n=int(m.sum()), rho=float(r.statistic), p=float(r.pvalue))


def mannwhitney(a, b):
    a = np.asarray([v for v in a if np.isfinite(v)], dtype=float)
    b = np.asarray([v for v in b if np.isfinite(v)], dtype=float)
    if a.size < 3 or b.size < 3:
        return dict(n_a=int(a.size), n_b=int(b.size), p=np.nan,
                    median_a=np.nan, median_b=np.nan)
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    return dict(n_a=int(a.size), n_b=int(b.size), p=float(u.pvalue),
                median_a=float(np.median(a)), median_b=float(np.median(b)))
