"""Statistics helpers shared by the paper-number scripts.

Extracted from stats_significance.ipynb so the interpolation, the oracle-codec
selection and the multiple-testing correction have exactly one implementation.
"""

from __future__ import annotations

import math
import numpy as np
import pandas as pd
from scipy import stats


# --- rate-distortion interpolation ------------------------------------------
# When True, each (region, codec) sweep is reduced to its Pareto frontier
# before interpolation: sorted by increasing RMSE, only points that also
# decrease in size are kept. A raw sweep can be non-monotonic (a coarser codec
# setting occasionally costs *more* bytes for *more* error), and interpolating
# through such a dominated point produces a target estimate no real
# configuration achieves. The paper lists this as an open item; setting this
# flag is how the reported numbers are checked against it.
USE_PARETO = False


def _pareto(x, y):
    """Keep points that are not dominated: increasing x must decrease y."""
    keep_x, keep_y, best = [], [], np.inf
    for xi, yi in zip(x, y):
        if yi < best:
            keep_x.append(xi)
            keep_y.append(yi)
            best = yi
    return np.asarray(keep_x), np.asarray(keep_y)


def _xy(df, xcol, ycol):
    """Sorted, duplicate-averaged (x, y) from a sweep; Pareto if USE_PARETO."""
    d = df[[xcol, ycol]].replace([np.inf, -np.inf], np.nan).dropna()
    d = d[(d[xcol] > 0) & (d[ycol] > 0)]
    if len(d) < 2:
        return None
    d = d.groupby(xcol, as_index=False)[ycol].mean().sort_values(xcol)
    x, y = d[xcol].to_numpy(), d[ycol].to_numpy()
    if USE_PARETO:
        # The frontier is always defined on (rmse, bytes): more error must cost
        # fewer bytes. Callers that interpolate the other way pass (bytes,
        # rmse), so the same rule applies after swapping.
        x, y = _pareto(x, y)
        if len(x) < 2:
            return None
    return x, y


# A rate-distortion curve is R(D) = min{ S(q) : d(q) <= D }. Requiring the
# target to be *bracketed* by two measured points, as the first version of this
# code did, silently drops a scene whenever a codec's coarsest setting is still
# more accurate than the target -- even though that setting is a perfectly
# admissible solution. Those dropped scenes are the smoothest ones, and they
# were dropped precisely at the loose targets where the verdict changes sides,
# so the censoring was not random. With FEASIBLE_CENSORING the constrained
# formulation is used instead: bracketed targets are interpolated as before,
# and an unbracketed but feasible target returns the best achieved point, which
# is a conservative bound rather than a missing value.
FEASIBLE_CENSORING = True


def interp_size_at_rmse(df, target, rmse_col="rmse"):
    """Smallest size achieving RMSE <= target."""
    xy = _xy(df, rmse_col, "total_bytes")
    if xy is None:
        return np.nan
    x, y = xy
    if x.min() <= target <= x.max():
        return float(np.exp(np.interp(np.log(target), np.log(x), np.log(y))))
    if not FEASIBLE_CENSORING:
        return np.nan
    if target > x.max():
        # The codec never gets that coarse: its coarsest point already meets
        # the target, so it is feasible at that point's size.
        return float(y[np.argmax(x)])
    return np.nan          # even the finest setting misses the target


def interp_rmse_at_size(df, target_bytes, rmse_col="rmse"):
    """Lowest RMSE achievable within a size budget."""
    xy = _xy(df, "total_bytes", rmse_col)
    if xy is None:
        return np.nan
    x, y = xy
    if x.min() <= target_bytes <= x.max():
        return float(np.exp(np.interp(np.log(target_bytes), np.log(x), np.log(y))))
    if not FEASIBLE_CENSORING:
        return np.nan
    if target_bytes > x.max():
        # Budget larger than anything the codec produces: its most accurate
        # point fits, so that is the achievable error.
        return float(y[np.argmax(x)])
    return np.nan          # even the smallest file exceeds the budget


def bracketed_size_at_rmse(df, target, rmse_col="rmse") -> bool:
    """Was the target actually bracketed by two measured points?

    Reported per target so a row built largely from censored bounds is not read
    as an estimate. See FEASIBLE_CENSORING.
    """
    xy = _xy(df, rmse_col, "total_bytes")
    return bool(xy is not None and xy[0].min() <= target <= xy[0].max())


def bracketed_rmse_at_size(df, target_bytes, rmse_col="rmse") -> bool:
    xy = _xy(df, "total_bytes", rmse_col)
    return bool(xy is not None and xy[0].min() <= target_bytes <= xy[0].max())


def any_bracketed(df, scene, codecs, target, kind, rmse_col="rmse") -> bool:
    """True if at least one codec brackets the target for this scene."""
    sub = df[df.scene == scene]
    fn = bracketed_size_at_rmse if kind == "rmse" else bracketed_rmse_at_size
    return any(fn(sub[sub.codec == c], target, rmse_col) for c in codecs)


# When True, a (region, target) pair enters the paired comparison only if the
# codec the oracle selects brackets the target with two measured points on
# BOTH grids, so that every ratio tested is an interpolated estimate rather
# than the bound FEASIBLE_CENSORING returns. Pairs that fail are excluded from
# that target and counted (see update_paper_numbers), not bounded. After the
# coarse end of every sweep was extended (extend_sweep.py), the pairs that
# still fail are those whose target exceeds the scene's own standard
# deviation, for which a constant fill already meets the target and the
# comparison is void.
REQUIRE_ORACLE_BRACKETED = True


def oracle_bracketed(df, scene, codecs, target, kind, rmse_col="rmse"):
    """Does the codec the oracle picks bracket the target on this grid?"""
    sub = df[df.scene == scene]
    if kind == "rmse":
        name, _ = argbest_codec_size(df, scene, codecs, target, rmse_col)
        fn = bracketed_size_at_rmse
    else:
        name, _ = argbest_codec_rmse(df, scene, codecs, target, rmse_col)
        fn = bracketed_rmse_at_size
    return bool(name is not None and fn(sub[sub.codec == name], target, rmse_col))


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
               q25=np.nan, q75=np.nan, p=np.nan, p_linear=np.nan,
               p_sign=np.nan, n_above_one=np.nan, ci_lo=np.nan, ci_hi=np.nan)
    if r.size == 0:
        return out
    out.update(median=float(np.median(r)), lo=float(r.min()), hi=float(r.max()),
               q25=float(np.percentile(r, 25)), q75=float(np.percentile(r, 75)))
    # Paired sign test: makes no symmetry assumption at all, which matters here
    # because the ratio distributions are strongly skewed over water scenes.
    n_above = int((r > 1).sum())
    n_tied = int((r == 1).sum())
    if r.size - n_tied > 0:
        out["p_sign"] = float(stats.binomtest(n_above, r.size - n_tied,
                                              0.5).pvalue)
        out["n_above_one"] = n_above
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


def partial_spearman(x, y, z):
    """Rank correlation of x and y with z held fixed (first-order partial).

    Computed on the ranks: rho_xy.z = (r_xy - r_xz r_yz) / sqrt((1-r_xz^2)(1-r_yz^2)),
    with the usual t approximation on n-3 degrees of freedom. Used to ask
    whether an association with latitude survives once the scene's spectral
    slope, which is measured per region, is controlled for.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(m.sum())
    if n < 6:
        return dict(n=n, rho=np.nan, p=np.nan, rho_xz=np.nan, rho_yz=np.nan)
    rx, ry, rz = (stats.rankdata(v[m]) for v in (x, y, z))
    r_xy = float(np.corrcoef(rx, ry)[0, 1])
    r_xz = float(np.corrcoef(rx, rz)[0, 1])
    r_yz = float(np.corrcoef(ry, rz)[0, 1])
    den = math.sqrt(max((1 - r_xz ** 2) * (1 - r_yz ** 2), 1e-300))
    rho = (r_xy - r_xz * r_yz) / den
    rho = max(-1.0, min(1.0, rho))
    df = n - 3
    t = rho * math.sqrt(df / max(1 - rho ** 2, 1e-300))
    p = float(2 * stats.t.sf(abs(t), df))
    return dict(n=n, rho=float(rho), p=p, rho_xz=r_xz, rho_yz=r_yz)


def mannwhitney(a, b):
    a = np.asarray([v for v in a if np.isfinite(v)], dtype=float)
    b = np.asarray([v for v in b if np.isfinite(v)], dtype=float)
    if a.size < 3 or b.size < 3:
        return dict(n_a=int(a.size), n_b=int(b.size), p=np.nan,
                    median_a=np.nan, median_b=np.nan)
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    return dict(n_a=int(a.size), n_b=int(b.size), p=float(u.pvalue),
                median_a=float(np.median(a)), median_b=float(np.median(b)))
