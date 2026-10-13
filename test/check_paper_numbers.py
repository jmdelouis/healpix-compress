"""Refuse to let the manuscript drift away from the run that produced it.

Every headline figure in the text is derived from `results/paper_numbers_*.json`
and from the generated tables. Nothing keeps the two in step: on 7 September the
sweep was re-run, the tables were regenerated, and the prose was not, so the
manuscript carried ratios from the previous sweep for two days without anything
complaining. This script complains.

It works by computing, from the current results, the exact string each claim
should contain, and asserting that string appears in the manuscript. A claim
whose value moved therefore fails loudly, and the report names the replacement.
It is deliberately literal: a check that tried to parse the prose would be
fragile, while a check that only asks "does the current number appear here"
cannot pass on a stale one.

Usage
-----
    python check_paper_numbers.py                 # report
    python check_paper_numbers.py --quiet         # exit status only
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys

import numpy as np
import pandas as pd

import paper_config as cfg

TEX = cfg.REPO_ROOT / "tex" / "tgrs_healpix_compress.tex"


def read_manuscript(path):
    """The manuscript with every \\input{sections/...} expanded in place.

    Since the split of revision 3 the text lives in one file per section; the
    figures the checks look for may sit in any of them, so the check runs on
    the concatenation, in document order. Table inputs are not expanded: the
    tables are generated, and checked at their source.
    """
    root = path.parent
    def expand(p):
        out = []
        for line in p.read_text(encoding="utf-8").splitlines(keepends=True):
            m = re.match(r"\\input\{(sections/[^}]+)\}", line.strip())
            if m:
                sub = root / (m.group(1) + ("" if m.group(1).endswith(".tex") else ".tex"))
                out.append(expand(sub) if sub.exists() else line)
            else:
                out.append(line)
        return "".join(out)
    return expand(path)


def _fmt(x, nd=2):
    return f"{x:.{nd}f}"


def load(metric):
    p = cfg.RESULTS_DIR / f"paper_numbers_{metric}.json"
    if not p.exists():
        raise SystemExit(f"{p} is missing; run update_paper_numbers.py --metric {metric}")
    return json.loads(p.read_text(encoding="utf-8"))


def target_median(doc, kind, target):
    for t in doc["targets"]:
        if t["kind"] == kind and abs(float(t["target"]) - target) < 1e-12:
            return float(np.median(list(t["ratios"].values()))), len(t["ratios"])
    return float("nan"), 0


def checks(tex: str):
    """(label, must_appear, note) triples computed from the current results."""
    self_, common = load("self"), load("common")
    out = []

    # --- within-representation ratios at each RMSE target -------------------
    # The text quotes the tight targets as a range rather than one by one, so
    # only its endpoints have to appear verbatim; an interior value is reported
    # for inspection but never asserted, or the check fires on prose that is
    # perfectly correct.
    meds = [(t,) + target_median(self_, "rmse", t) for t in cfg.TARGET_RMSE]
    tight = [m for t, m, n in meds if n and t <= 1e-3]
    for t, m, n in meds:
        if not n:
            continue
        interior = (t <= 1e-3 and len(tight) > 2
                    and min(tight) < m < max(tight))
        out.append((f"self ratio @ RMSE {t:g}"
                    + ("  [range interior, not asserted]" if interior else ""),
                    _fmt(m), f"median over {n} regions"))

    # --- the gamma = 1 verdict at the loosest target rests on one paired
    # test at p = 0.03 uncorrected; the prose gives its confidence interval
    # so that a reader sees how close to parity the upper bound sits.
    for t in self_["targets"]:
        if t["kind"] == "rmse" and abs(float(t["target"]) - cfg.TARGET_RMSE[-1]) < 1e-12:
            g1 = t.get("test_gamma1") or {}
            if "ci_lo" in g1 and "ci_hi" in g1:
                out.append(("gamma=1 CI at the loosest target [in sections/05b_rd_main.tex]",
                            f"$0.89$--$0.99$".replace("0.89", f"{g1['ci_lo']:.2f}").replace("0.99", f"{g1['ci_hi']:.2f}"),
                            "95% CI of the median ratio, from paper_numbers_self"))
                out.append(("gamma=1 p at the loosest target [in sections/05b_rd_main.tex]",
                            f"$p={g1['p']:.2f}$", "paired signed-rank, uncorrected"))
    # --- the two the abstract and conclusion quote by name ------------------
    m_loose, _ = target_median(self_, "rmse", cfg.TARGET_RMSE[-1])
    out.append(("loosest-target ratio quoted in abstract/conclusion",
                _fmt(m_loose), "appears as $0.xx$ in several places"))

    # --- common-reference budgets -------------------------------------------
    for b in (0.02, 0.05, 0.10):
        m, n = target_median(common, "budget", b)
        if n:
            out.append((f"common ratio @ budget {100*b:g}%",
                        _fmt(m) if m < 10 else _fmt(m, 1),
                        f"median over {n} regions"))

    # --- critical overhead ---------------------------------------------------
    for k, v in self_.get("gamma_star", {}).items():
        out.append((f"gamma_star @ {k}", _fmt(float(v)),
                    "must match table_gamma_continuous; re-run gamma_analysis.py"))

    # --- raw accounting ------------------------------------------------------
    acc = self_.get("storage_accounting", {})
    for key, nd in (("ratio_gamma1", 3), ("ratio_gamma", 2)):
        v = acc.get(key)
        if isinstance(v, dict):
            v = v.get("median")
        if v is not None:
            out.append((f"raw {key}", _fmt(float(v), nd), "geometric constant"))
    if acc.get("ratio_uint16") is not None:
        out.append(("raw ratio, native as distributed (uint16) [in sections/05a_raw.tex]",
                    f"${acc['ratio_uint16']:.2f}\\times$",
                    "2 x the float32 ratio; Table IV last row"))

    # --- how many latitude associations survive Holm -------------------------
    lat = self_.get("latitude", {})
    survivors = sorted(k for k, v in lat.items()
                       if min(v.get("spearman", {}).get("p_holm", 1.0),
                              v.get("equatorial_vs_cap", {}).get("p_holm", 1.0)) < 0.05)
    out.append(("number of latitude associations surviving Holm",
                str(len(survivors)),
                "survivors: " + (", ".join(survivors) or "none")))
    # Revision 5: the association at the tightest budget, with the spectral
    # slope held fixed (partial Spearman), and the slope's own relation to
    # latitude; both quoted in Section V-D.
    tb = lat.get("budget:0.02", {})
    pb, bl = tb.get("partial_beta"), tb.get("beta_vs_lat")
    if pb and bl:
        where = "[in sections/05d_latitude.tex]"
        out.append((f"latitude | beta, partial rho at the 2% budget {where}",
                    f"$\\rho={pb['rho']:.2f}$", "paper_numbers_self latitude.partial_beta"))
        out.append((f"latitude | beta, partial p at the 2% budget {where}",
                    f"$p={pb['p']:.3f}$", ""))
        out.append((f"latitude | beta, n at the 2% budget {where}",
                    f"$n={pb['n']}$", ""))
        out.append((f"beta vs latitude, rho {where}",
                    f"$\\rho={bl['rho']:.2f}$", ""))
        out.append((f"beta vs latitude, p {where}",
                    f"$p={bl['p']:.2f}$", ""))
    else:
        out.append(("latitude | beta  [skipped: rerun update_paper_numbers.py "
                    "with support_controls_spectrum.csv present]", "", ""))
    out += level_and_reader_checks()
    out += gamma_continuous_checks()
    out += headline_checks()
    out += censoring_checks()
    out += same_codec_checks()
    out += transfer_mechanism_checks()
    out += transfer_cost_prose_checks()
    out += scoped_checks()
    out += level19_checks()
    out += window_radiometry_checks()
    out += resampling_error_checks()
    out += target_scale_checks()
    out += bicubic_diag_checks()
    out += tile_crop_checks()
    out += halo_layout_checks()
    out += coauthor_review_checks()
    out += appendix_table_checks()
    out += consistency_checks(tex)
    out += wording_checks(tex)
    # Revision 7: the PSF-aware carrier is no longer scored against the
    # observation (its estimand is the latent field), so the PSF-chain and
    # matched-reader blocks of revision 6 are gone with the prose they guarded.
    return out, survivors



def gamma_continuous_checks():
    """The win-fraction ranges and the inconclusive band of Section V-B.

    They are read off the same per-region critical overheads that
    gamma_analysis.py turns into table_gamma_continuous, so the prose cannot
    drift from the table without this firing (revision 3: it had, by 6 points
    at gamma = 1.3).
    """
    try:
        import gamma_analysis as ga
        df, _, _ = ga.analyse("self")
    except Exception as exc:                      # noqa: BLE001
        return [("gamma-continuous ranges  [skipped: "
                 f"{type(exc).__name__}]", "", str(exc)[:60])]
    rep = df[df.reported] if "reported" in df else df
    out = []
    for g, label in ((1.0, "gamma=1"), (1.3, "gamma=1.3"), (2.0, "gamma=2")):
        w = 100 * rep[f"win_at_{g:g}"]
        out.append((f"HPX-smaller range at {label}",
                    f"${w.min():.0f}\\%$ to ${w.max():.0f}\\%$",
                    "over the reported RMSE targets"))
    r = rep[np.isclose(rep.target, 1e-3)]
    if len(r) and np.isfinite(r.inconclusive_lo.iloc[0]):
        out.append(("inconclusive gamma band at RMSE 1e-3",
                    f"${r.inconclusive_lo.iloc[0]:.2f}$--${r.inconclusive_hi.iloc[0]:.2f}$",
                    "paired test does not reject parity"))
    if "recovered_fraction" in rep and len(rep):
        s_ = rep.sort_values("target")
        where = "[in sections/05b_rd_main.tex]"
        out.append((f"recovered fraction f(e) at the tightest target {where}",
                    f"${s_.recovered_fraction.iloc[0]:.2f}$", "(alpha - gamma_star)/(alpha - 1)"))
        out.append((f"recovered fraction f(e) at the loosest target {where}",
                    f"${s_.recovered_fraction.iloc[-1]:.2f}$", ""))
    return out


def headline_checks():
    """The plain-language sizes of Section V-B, the abstract and the conclusion.

    Every figure comes from headline.py, which assembles them from the results
    the other sections are checked against; the abstract quotes them by value,
    so each one is asserted verbatim.
    """
    try:
        import headline as hl
        per_target, raw, transfers, std, gamma = hl.load()
    except Exception as exc:                      # noqa: BLE001
        return [("headline sizes  [skipped: "
                 f"{type(exc).__name__}]", "", str(exc)[:60])]
    out = []
    ts = sorted(per_target)
    for name, t in (("tightest", ts[0]), ("loosest", ts[-1])):
        bc = per_target[t]["by_class"].values()
        out.append((f"headline class range, {name} target",
                    f"${min(bc):.2f}$--${max(bc):.2f}$",
                    "range of the four class medians, gamma=1"))
    pct = [100 * ts[0] / std, 100 * ts[-1] / std]
    words = "quarter of a percent" if abs(pct[0] - 0.25) < 0.03 else f"{pct[0]:.2g}\\%"
    out.append(("tightest target as share of contrast", words,
                f"{pct[0]:.3f}% of the median tile std {std:.4f}"))
    words = "quarter of the contrast" if abs(pct[1] - 25) < 3 else f"{pct[1]:.2g}\\%"
    out.append(("loosest target as share of contrast", words,
                f"{pct[1]:.2f}% of the median tile std {std:.4f}"))
    out.append(("median tile contrast (glossary)", f"${std:.2f}$",
                "median std of the level-20 tiles"))
    if transfers:
        vals = [transfers[t] for t in sorted(transfers)]
        for key, label in (("psf", "PSF-aware"), ("idw", "inverse-distance")):
            lo, hi = min(v[key] for v in vals), max(v[key] for v in vals)
            out.append((f"{label} tile / UTM at gamma=1, range",
                        f"${lo:.1f}$--${hi:.1f}$", "abstract and conclusion"))
            out.append((f"{label} tile / UTM at gamma=1, tightest",
                        f"${vals[0][key]:.2f}$", "Section V-J, final sizes"))
            out.append((f"{label} tile / UTM at gamma=1, loosest",
                        f"${vals[-1][key]:.2f}$", "Section V-J, final sizes"))
        out.append(("PSF-aware tile / UTM at gamma=1.3, range",
                    f"${min(v['psf'] for v in vals)/gamma:.1f}$ to "
                    f"${max(v['psf'] for v in vals)/gamma:.1f}$",
                    "Section V-B, illustrative overhead"))
    return out


def scoped_checks():
    """Figures asserted in ONE section rather than anywhere in the manuscript.

    The default rule -- does this value appear somewhere in the text -- cannot
    see a stale copy of a figure as long as the fresh value occurs somewhere
    else. That is exactly how the headline subsection kept the pre-extension
    ratios for a day while the checker passed: Section V-C quoted the new
    loosest-target ratio, so the assertion was satisfied, and Section V-B
    still said the old one. A claim that a given section must carry a given
    number closes that hole, and the sections listed here are the ones a
    reader meets first.

    Labels carry ``[in <path>]``; main() restricts the search to that file.
    """
    out = []
    try:
        import headline as hl
        per_target, raw, transfers, std, gamma = hl.load()
    except Exception as exc:                              # noqa: BLE001
        return [(f"scoped headline figures  [skipped: {type(exc).__name__}]",
                 "", str(exc)[:60])]
    ts = sorted(per_target)
    loose, tight = per_target[ts[-1]], per_target[ts[0]]
    hd = "[in sections/05a2_headline.tex]"
    out.append((f"headline section: loosest ratio {hd}",
                _fmt(loose["median"]), "gamma=1"))
    out.append((f"headline section: tightest ratio {hd}",
                _fmt(tight["median"]), "gamma=1"))
    bc = loose["by_class"].values()
    out.append((f"headline section: class range at the loosest target {hd}",
                f"${min(bc):.2f}$--${max(bc):.2f}$", "four class medians"))
    out.append((f"headline section: loosest ratio at gamma_tile {hd}",
                _fmt(loose["median"] / gamma), f"gamma={gamma:g}"))
    if transfers:
        vals = [transfers[t] for t in sorted(transfers)]
        for key, label in (("psf", "PSF-aware"), ("idw", "inverse-distance")):
            lo, hi = min(v[key] for v in vals), max(v[key] for v in vals)
            out.append((f"headline section: {label} range {hd}",
                        f"${lo:.1f}$--${hi:.1f}$", "gamma=1"))
        # Revision 5: the C--T median on the three-way population is a column
        # of the merged transfer table (merged_tables.py) and is no longer
        # quoted in the headline prose.
    return out


def level19_checks():
    """Section V-J: the two level-19 carriers, from level19_direct.py's digest.

    Both floors and the ratio between them are quoted in the prose, and the
    whole point of the subsection is that the two agree, so a drift in either
    would change the argument rather than only a figure.
    """
    p = cfg.RESULTS_DIR / "level19_direct.json"
    if not p.exists():
        return [("level-19 carriers  [skipped: no digest]", "",
                 "run level19_direct.py")]
    d = json.loads(p.read_text())
    if "n" not in d:                       # older per-reader digest
        d = d.get("nearest", {})
    if not d:
        return [("level-19 carriers  [skipped: unreadable digest]", "", str(p))]
    # Revision 6: the section was folded into Section V-K and Appendix B-F.
    where = "[in sections/08_appendices.tex]"
    out = [
        (f"level-19 direct floor {where}", f"${d['floor19_direct']:.2e}$"
         .replace("e-0", r"\times10^{-").replace("$", "$", 1) + "",
         "median over the regions"),
        (f"level-19 direct / aggregate {where}",
         f"${d['ratio_direct_over_agg']:.3f}$", "median ratio"),
        (f"level-19 direct / level 20 {where}", f"${d['ratio_direct']:.2f}$", ""),
        (f"level-19 aggregate / level 20 {where}", f"${d['ratio_agg']:.2f}$", ""),
        (f"level-19 carriers within 5% {where}", f"{d['within_5pct']:d}",
         f"of {d['n']} regions"),
        (f"level-19 eligible regions {where}", f"{d['n']:d}",
         "regions the window audit keeps"),
        (f"level-19 widest departure {where}",
         f"${d['max_ratio_direct_over_agg']:.3f}$", "over the eligible regions"),
        (f"level-19 direct better on {where}", f"{d['direct_better_on']:d}",
         "regions"),
        (f"level-19 carrier difference, % of std {where}",
         f"${100*d['carrier_rms_frac_of_std']:.1f}\\%$", "RMS"),
        (f"level-19 carrier difference / floor {where}",
         f"${d['carrier_rms_over_floor19']:.2f}$", "RMS over the L19 floor"),
    ]
    # the floor itself is quoted in scientific notation in the prose
    mant, exp = f"{d['floor19_direct']:.2e}".split("e")
    out[0] = (f"level-19 direct floor {where}",
              f"${mant}\\times10^{{{int(exp)}}}$", "median over the regions")
    mant, exp = f"{d['floor19_agg']:.2e}".split("e")
    out.append((f"level-19 aggregate floor {where}",
                f"${mant}\\times10^{{{int(exp)}}}$", "median over the regions"))
    return out


def _table_rows(name):
    """Data rows of a generated table, split on & with the \\ stripped."""
    q = cfg.TEX_TABLES / f"{name}.tex"
    if not q.exists():
        return []
    rows = []
    body = q.read_text(encoding="utf-8")
    body = body.split("\\midrule")[-1].split("\\bottomrule")[0]
    for line in body.splitlines():
        line = line.strip()
        if not line or line.startswith("\\") or "&" not in line:
            continue
        rows.append([c.strip() for c in line.rstrip("\\").split("&")])
    return rows


def appendix_table_checks():
    """Appendix C and D quote their own tables; make the prose answer to them.

    This is the pass that was missing when the water class's crossover moved
    from the loosest target to 1e-3 and the text kept the old reading for two
    revisions. The source of truth here is the generated table itself rather
    than a results file, because that is what the reader sees beside the
    paragraph, and a paragraph that disagrees with the table on the same page
    is the failure this is for.
    """
    out = []
    where = "[in sections/08_appendices.tex]"

    # Table XXVI: class, crossover target, HPX/UTM, wins, own e/kappa/infl,
    # latent e/kappa/infl.
    for r in _table_rows("table_alpha_budget"):
        if len(r) < 10:
            continue
        cls = r[0].lower()
        own_k, latent_k = r[5], r[8]
        own_i = r[6].split("(")[0].strip()
        latent_i = r[9].split("(")[0].strip()
        # The prose writes these as $\\kappa=0.19$ and, in a list, as a bare
        # $2.28$, so the assertion is on the decimal itself; the section scope
        # is what keeps it from matching an unrelated number elsewhere.
        out += [
            (f"alpha budget, {cls}: latent kappa {where}", latent_k,
             "crossover kappa against the latent scene"),
            (f"alpha budget, {cls}: own-estimand kappa {where}", own_k,
             "crossover kappa against the own estimand"),
        ]
        # The measured inflation needs the cross term of cross_term.py; a
        # class the table leaves as "---" has none and the prose quotes only
        # the orthogonal value there.
        if latent_i not in ("---", ""):
            out.append((f"alpha budget, {cls}: latent inflation {where}",
                        f"{latent_i}\\%", "measured, not orthogonal"))
        del own_i

    # Table XXV: the psf recovery columns, whose class range the prose quotes.
    rec = [(r[0].lower(), r[-2], r[-1]) for r in _table_rows("table_resampling_error")
           if len(r) >= 9]
    if rec:
        meas = sorted(float(m.replace("+", "")) for _c, m, _n in rec)
        out += [
            (f"psf recovery, worst class {where}",
             f"${meas[0]:+.1f}\\%$".replace("+-", "-"), "measured noise"),
            (f"psf recovery, best class {where}",
             f"${meas[-1]:+.1f}\\%$", "measured noise"),
        ]
        ctrl = [float(n.replace("+", "")) for _c, _m, n in rec if n not in ("---", "")]
        if ctrl:
            ctrl.sort()
            out += [
                (f"psf recovery, noiseless, worst class {where}",
                 f"${ctrl[0]:+.1f}\\%$".replace("+-", "-"), "noiseless control"),
                (f"psf recovery, noiseless, best class {where}",
                 f"${ctrl[-1]:+.1f}\\%$", "noiseless control"),
            ]
    out += psf_recovery_checks()
    return out


def psf_recovery_checks():
    """Appendix D: the recovery per region against the fitted response width.

    Revision 8. With the harness corrected (field blurred at the right sigma,
    inverse handed the operator scale s = 2 sigma), the noiseless recovery
    splits cleanly on the region's fitted FWHM: at or above 11 m every region
    recovers, below it the blur is a few 1e-4 and the inversion's own error
    exceeds it. The prose quotes that split, the counts of regions recovering
    with and without noise, and the rank correlation.
    """
    import pandas as pd
    q = cfg.RESULTS_DIR / "resampling_error_noise0.csv"
    m = cfg.RESULTS_DIR / "resampling_error.csv"
    if not q.exists() or not m.exists():
        return [("psf recovery per region  [skipped: resampling_error CSVs absent]", "", "")]
    where = "[in sections/08_appendices.tex]"
    def rec(path):
        d = pd.read_csv(path)
        d = d[(d.kernel == "psf")].dropna(subset=["e_resample", "blur_error"]).copy()
        d["rec"] = 1.0 - d.e_resample / d.blur_error
        return d
    n0, nm = rec(q), rec(m)
    wide, narrow = n0[n0.psf_fwhm_m >= 11.0], n0[n0.psf_fwhm_m < 11.0]
    out = [
        (f"noiseless recovery, regions with FWHM >= 11 m, count {where}", f"{len(wide)}", ""),
        (f"noiseless recovery, regions with FWHM >= 11 m, low {where}",
         f"${100*wide.rec.min():.0f}\\%$", ""),
        (f"noiseless recovery, regions with FWHM >= 11 m, high {where}",
         f"${100*wide.rec.max():.0f}\\%$", ""),
        (f"noiseless recovery, regions with FWHM < 11 m, count {where}", f"{len(narrow)}", ""),
        (f"noiseless recovery, regions with FWHM < 11 m, high {where}",
         f"${100*narrow.rec.max():+.0f}\\%$", "the best of them"),
        (f"noiseless recovery, positive regions {where}",
         f"{int((n0.rec > 0).sum())} of the {len(n0)}", ""),
        (f"measured-noise recovery, positive regions {where}",
         f"{int((nm.rec > 0).sum())} of the {len(nm)}", ""),
        (f"noiseless recovery vs FWHM, Spearman rho {where}",
         f"$\\rho={n0[['psf_fwhm_m', 'rec']].corr(method='spearman').iloc[0, 1]:.2f}$", ""),
        (f"narrowest fitted FWHM {where}", f"\\SI{{{n0.psf_fwhm_m.min():.1f}}}{{m}}", ""),
    ]
    return out


def tile_crop_checks():
    """The exact crop test the window appendix quotes as its confirmation.

    Two numbers carry the argument -- how many stored tiles agree with their
    window to the bit, and how many pixels differ on the one that does not --
    plus the count of stale recorded boxes, which is the trap the passage
    warns about. All three come from the CSV rather than from the prose.
    """
    p = cfg.RESULTS_DIR / "tile_crop_consistency.csv"
    if not p.exists():
        return [("tile/crop audit  [skipped: no CSV]", "",
                 "run tile_crop_consistency.py")]
    import pandas as pd
    d = pd.read_csv(p)
    where = "[in sections/08_appendices.tex]"
    out = [(f"tile/crop audit, tiles checked {where}", f"{len(d):d}",
            "one row per (region, transfer)"),
           (f"tile/crop audit, agreeing {where}", f"{int(d.identical.sum()):d}",
            "agree with the window to the bit")]
    bad = d[~d.identical]
    for _, r in bad.iterrows():
        n = int(r["n_differ"])
        out.append((f"tile/crop audit, pixels differing on {r['scene']} {where}",
                    f"{n // 1000}\\,{n % 1000:03d}", "of the crop"))
    if "box_shifted" in d.columns:
        # Per region, not per tile: the recorded box is stale because the
        # window was re-acquired, which is a property of the region.
        n_shift = int(d[d.box_shifted].scene.nunique())
        n_rows = int(d.box_shifted.sum())
        # A bare "36" is satisfied by any 36 in the appendix (a percentage,
        # an n), which is how the prose said "39 of the 40 windows" -- the
        # row count -- for a while while this check passed. Lock the phrase.
        out.append((f"tile/crop audit, stale recorded boxes {where}",
                    f"{n_shift:d} of the {d.scene.nunique():d} windows",
                    "regions re-acquired at another size"))
        out.append((f"tile/crop audit, mismatches at recorded indices {where}",
                    f"reports {n_rows:d} mismatched tiles",
                    "rows, i.e. (region, transfer) pairs"))
    return out


def halo_layout_checks():
    """Section VI-B (revision 8): the native pixels that build one tile.

    A store that keeps the native product and runs the transfer at read time
    needs the tile footprint plus a halo for the kernel support. The halo
    overhead is pure geometry: (side + 2h)^2 / side^2 - 1, quoted for a halo
    of four cells on the 1024-cell tile and four native pixels on the
    637-pixel equal-area crop. The bounding box of the ~45-degree footprint
    (paper_config.TILE_BBOX_PX) is the rectangle a 2-D codec would need.
    """
    where = "[in sections/06_discussion.tex]"
    h = 4
    cells = (cfg.TILE_SIDE + 2 * h) ** 2 / cfg.TILE_SIDE ** 2 - 1.0
    px = (cfg.UTM_SIDE_PX + 2 * h) ** 2 / cfg.UTM_SIDE_PX ** 2 - 1.0
    bbox = cfg.TILE_BBOX_PX ** 2 / cfg.UTM_SIDE_PX ** 2
    return [
        (f"halo of {h} cells, tile enlargement {where}",
         f"$({cfg.TILE_SIDE + 2 * h}/{cfg.TILE_SIDE})^2-1\\approx{100 * cells:.1f}\\%$",
         "pure geometry"),
        (f"halo of {h} native pixels, crop enlargement {where}",
         f"${100 * px:.1f}\\%$", f"on the {cfg.UTM_SIDE_PX}-pixel crop"),
        (f"footprint bounding box / footprint area {where}",
         f"${bbox:.1f}$ times its area", f"{cfg.TILE_BBOX_PX}^2 / {cfg.UTM_SIDE_PX}^2"),
        ("halo overhead, conclusion [in sections/07_conclusion.tex]",
         f"about ${round(100 * (cells + px) / 2):d}\\%$", "rounded midpoint of the two halos"),
    ]


def coauthor_review_checks():
    """Figures rewritten after the coauthor review on Overleaf (revision 9).

    Three sentences quoted a figure loosely or from an earlier run: the share
    of regions where native UTM stays ahead at gamma = 2 and the tightest
    target ("three quarters" -> the table's own 72 %), the area of the native
    rectangle of the common-support control (a 405 x 401 rectangle from an
    earlier run; the stored control has one rectangle per region), and the
    reader-floor halving, where the paired per-region ratio (0.582) and the
    ratio of the two medians are different statistics and are now both named.
    """
    out = []
    g = _csv("gamma_analysis_self.csv")
    if g is not None and "win_at_2" in g:
        r = g[(g.metric == "self") & (g.target == g.target.min())]
        if len(r):
            out.append(("UTM ahead at gamma = 2, tightest target, % of regions "
                        "[in sections/05b_rd_main.tex]",
                        f"${100 * (1 - float(r.win_at_2.iloc[0])):.0f}\\%$",
                        "1 - win_at_2 of the synthesis table"))
    c = _csv("support_controls_common.csv")
    if c is not None and "utm_km2" in c:
        a = c.groupby("scene").utm_km2.first()
        h = c.groupby("scene").hpx_km2.first()
        where = "[in sections/08_appendices.tex]"
        out += [(f"common-support control, native rectangle area range {where}",
                 f"\\SIrange{{{a.min():.1f}}}{{{a.max():.1f}}}{{\\kilo\\metre\\squared}}",
                 f"over {len(a)} regions"),
                (f"common-support control, native rectangle median area {where}",
                 f"\\SI{{{a.median():.1f}}}{{\\kilo\\metre\\squared}}", ""),
                (f"common-support control, HEALPix block area {where}",
                 f"\\SI{{{h.iloc[0]:.1f}}}{{\\kilo\\metre\\squared}}", "256 x 256 cells")]
    rf = _csv("reader_floor.csv")
    if rf is not None:
        ok = rf[rf.error.fillna("").astype(str).str.len() == 0]
        ct = ok[ok.kernel == "clough_tocher"].pivot_table(
            index="scene", columns="reader", values="floor")
        if "nearest" in ct and "bilinear" in ct:
            out.append(("bilinear / nearest floor, ratio of the two medians "
                        "[in sections/05i_common_reference.tex]",
                        f"${ct.bilinear.median() / ct.nearest.median():.2f}$",
                        "not the paired median, which is 0.582"))
    return out


def target_scale_checks():
    """Section IV-B: the RMSE targets read against each transfer's own error.

    Revision 7. The Clough--Tocher error against the observed field and the
    blur the PSF-aware transfer is charged against the latent field both come
    from resampling_error.csv (resampling_error.py). The compression error
    that would cancel the companion's recovery is derived here: with recovery
    r = 1 - e_psf/e_blur and an orthogonal compression error e_c, the total
    returns to e_blur when e_c = e_blur * sqrt(2r - r^2); r is the companion
    study's lower figure, 15.4%, which the prose quotes next to it.
    """
    p = cfg.RESULTS_DIR / "resampling_error.csv"
    if not p.exists():
        return [("target scale  [skipped: resampling_error.csv absent]", "",
                 "run resampling_error.py")]
    import pandas as pd
    d = pd.read_csv(p)
    d["cls"] = d.scene.str.split("__").str[1]
    where = "[in sections/04_design.tex]"
    e_ct = d[d.kernel == "clough_tocher"].e_resample.median()
    targets = (1e-4, 3e-3, 1e-2)
    out = [
        (f"C-T error vs observed, median {where}",
         f"${e_ct*1e3:.1f}\\times10^{{-3}}$", "same value as Appendix D"),
        (f"tightest target / C-T error, percent {where}",
         f"${100*targets[0]/e_ct:.0f}\\%$", ""),
        (f"3e-3 target / C-T error {where}", f"${targets[1]/e_ct:.1f}$ times", ""),
        (f"loosest target / C-T error {where}", f"${targets[2]/e_ct:.1f}$ times", ""),
    ]
    psf = d[d.kernel == "psf"]
    if len(psf):
        b = psf.groupby("cls").blur_error.median()
        def sci(x):
            e = int(np.floor(np.log10(x)))
            return f"${x/10**e:.1f}\\times10^{{{e}}}$"
        out += [(f"PSF blur vs latent, smallest class median {where}", sci(b.min()),
                 f"class {b.idxmin()}"),
                (f"PSF blur vs latent, largest class median {where}", sci(b.max()),
                 f"class {b.idxmax()}")]
        r = 0.154
        cancel = b * np.sqrt(2 * r - r * r)
        out += [(f"compression error cancelling the 15.4% recovery, smallest {where}",
                 sci(cancel.min()), f"class {cancel.idxmin()}"),
                (f"compression error cancelling the 15.4% recovery, largest {where}",
                 sci(cancel.max()), f"class {cancel.idxmax()}")]
    return out


def bicubic_diag_checks():
    """Section V-L: why the bicubic row of Table XIII(c) is marked invalid.

    Revision 8. diag_bicubic.py runs four probe fields through the radial
    bicubic operator on one region per class and a separable bicubic
    (scipy, order 3) on the same analytic field. The prose quotes the share
    of cells on which a constant is not reproduced and the separable
    kernel's transfer error; both are medians over the diagnosed regions.
    """
    p = cfg.RESULTS_DIR / "diag_bicubic.csv"
    if not p.exists():
        return [("bicubic diagnostic  [skipped: diag_bicubic.csv absent]", "",
                 "run diag_bicubic.py")]
    import pandas as pd
    d = pd.read_csv(p)
    where = "[in sections/05l_transfer_cost.tex]"
    c = d[(d.probe == "constant") & (d.variant == "radial_default")]
    s = d[(d.probe == "analytic") & (d.variant == "separable_scipy_order3")]
    out = []
    if len(c):
        out.append((f"bicubic: cells failing the constant test, percent (several places)",
                    f"${100*c.bad_frac_5std.median():.1f}\\%$",
                    "V-L and the table footnote; |value - 0.3| above float noise"))
    if len(s):
        out.append((f"separable bicubic transfer error, percent of std {where}",
                    f"${100*s.rmse_over_std.median():.1f}\\%$", "median over regions"))
    return out


def resampling_error_checks():
    """The appendix figures on each transfer's error against its own estimand.

    These were quoted in the prose and produced by resampling_error.py, but
    nothing tied the two together, so the region count could change under them
    without anything complaining -- which is what happened when the window
    audit removed a region. Every figure the passage states and this script
    produces is asserted here.
    """
    p = cfg.RESULTS_DIR / "resampling_error.csv"
    if not p.exists():
        return [("resampling error  [skipped: no CSV]", "",
                 "run resampling_error.py")]
    import pandas as pd
    d = pd.read_csv(p)
    d["cls"] = d.scene.str.split("__").str[1]
    ct = d[d.kernel == "clough_tocher"]
    where = "[in sections/08_appendices.tex]"
    by_cls = 100 * ct.groupby("cls").e_resample_frac_of_std.median()
    noise = 100 * ct.groupby("cls").noise_frac_of_std.median()
    out = [
        (f"C-T error vs observed, median {where}",
         f"${ct.e_resample.median()*1e3:.1f}\\times10^{{-3}}$", ""),
        (f"C-T error vs latent, median {where}",
         f"${ct.e_against_other_estimand.median()*1e2:.1f}\\times10^{{-2}}$", ""),
        (f"C-T error, class range low {where}", f"${by_cls.min():.1f}\\%$", ""),
        (f"C-T error, class range high {where}", f"${by_cls.max():.1f}\\%$", ""),
        # Written as one range in the prose, so asserted as one range.
        (f"noise level, class range {where}",
         f"${noise.min():.1f}$--${noise.max():.1f}\\%$",
         "Immerkaer estimate, by class"),
    ]
    # The blur the non-deconvolving transfer leaves, as a share of the latent
    # field's contrast. Normalised on the PSF rows, whose estimand IS the
    # latent field, so that truth_std is the latent standard deviation: taken
    # from the Clough--Tocher rows it would divide by the observed field's
    # contrast instead and read 23.4% and 141.1%.
    psf = d[d.kernel == "psf"].copy()
    if len(psf):
        psf["blur_frac"] = 100 * psf.blur_error / psf.truth_std
        b = psf.groupby("cls").blur_frac.median()
        for c in ("forest", "urban"):
            if c in b.index:
                out.append((f"blur / latent std, {c} {where}",
                            f"${b[c]:.1f}\\%$", "PSF rows: truth_std is latent"))
    q = cfg.RESULTS_DIR / "resampling_error_noise0.csv"
    if q.exists():
        n = pd.read_csv(q)
        n["cls"] = n.scene.str.split("__").str[1]
        nc = 100 * n[n.kernel == "clough_tocher"].groupby("cls") \
            .e_resample_frac_of_std.median()
        out += [
            (f"noiseless C-T error, class range low {where}",
             f"${nc.min():.1f}\\%$", ""),
            (f"noiseless C-T error, class range high {where}",
             f"${nc.max():.1f}\\%$", ""),
        ]
    return out


def window_radiometry_checks():
    """Appendix: the audit that decides which regions may be rebuilt.

    The appendix argues from three numbers -- the median of the identity, the
    interval the eligible regions occupy, and where the excluded one sits --
    and the exclusion it licenses changes the level-19 table's n. If the audit
    is re-run and any of them moves, the appendix is wrong and so is the
    table's note, so they are checked against the CSV rather than the prose.
    """
    p = cfg.RESULTS_DIR / "window_radiometry.csv"
    if not p.exists():
        return [("window radiometry  [skipped: no audit]", "",
                 "run window_radiometry.py")]
    import pandas as pd
    a = pd.read_csv(p)
    if "slope_over_corr2" not in a.columns:
        return [("window radiometry  [skipped: old audit format]", "",
                 "re-run window_radiometry.py")]
    # The audit now covers one tile per transfer. The appendix reports the
    # carrier's, which is what every published number rests on; the other
    # transfers' tiles are a separate question (the identity assumes a
    # transfer that does not systematically alter the field, which a
    # deconvolving one does by design).
    if "kernel" in a.columns:
        a = a[a.kernel == getattr(cfg, "CARRIER_KERNEL", "clough_tocher")]
    ok, bad = a[a.status == "ok"], a[a.status == "MISMATCH"]
    where = "[in sections/08_appendices.tex]"
    out = [
        (f"window audit, median identity {where}",
         f"${a.slope_over_corr2.median():.4f}$", "slope / corr^2 over 40"),
        # The prose writes the interval as one math group, so assert the
        # interval and not its endpoints: two separate assertions would pass
        # on a text that quoted the two numbers anywhere at all.
        (f"window audit, eligible interval {where}",
         f"[{ok.slope_over_corr2.min():.3f}, {ok.slope_over_corr2.max():.3f}]",
         "slope / corr^2 over the eligible regions"),
        (f"window audit, eligible correlation low {where}",
         f"${ok['corr'].min():.3f}$", ""),
        (f"window audit, eligible correlation high {where}",
         f"${ok['corr'].max():.3f}$", ""),
    ]
    for _, r in bad.iterrows():
        out.append((f"window audit, excluded {r['scene']} {where}",
                    f"${r['slope_over_corr2']:.3f}$", "breaks the identity"))
    return out


def transfer_mechanism_checks():
    """The lossless-size correlation Section V-K quotes as evidence.

    It is a rank correlation over regions between how much a raster carries
    and the surcharge it pays, and it had drifted once already because no
    script produced it; superres_cost.py now writes it into the digest.
    """
    p = cfg.RESULTS_DIR / "transfer_cost_digest.json"
    if not p.exists():
        return [("lossless-size correlation  [skipped: no digest]", "", str(p))]
    d = json.loads(p.read_text()).get("lossless_size_vs_ratio")
    if not d:
        return [("lossless-size correlation  [skipped: not in digest]", "",
                 "re-run superres_cost.py")]
    out = []
    for key, label in (("psf", "PSF"), ("linear", "inverse-distance")):
        v = d.get(key)
        if v:
            out.append((f"lossless size vs surcharge, {label}",
                        f"${v['rho']:.2f}$", f"n={v['n']}, p={v['p']:.3f}"))
            # The prose leans on the p-values as much as on the sign -- one of
            # the two does not reach 0.05 and the sentence says so -- and a
            # rho that moved without its p would leave that reading stale.
            out.append((f"lossless size vs surcharge, {label}, p",
                        f"$p={v['p']:.3f}$", f"n={v['n']}"))
            out.append((f"lossless size vs surcharge, {label}, n",
                        f"$n={v['n']}$", "regions"))
    return out


def transfer_cost_prose_checks():
    """Every figure Section V-M quotes from Table XVIII and Table XIX.

    Until revision 5 none of these was asserted, and the two tables were
    built under a weaker eligibility rule than the rest of the paper; when
    the rule was harmonised, eight prose figures moved and only the two the
    headline check happened to cover were caught. Everything here is read
    from the digest superres_cost.py writes beside the tables.
    """
    p = cfg.RESULTS_DIR / "transfer_cost_digest.json"
    if not p.exists():
        return [("transfer-cost prose  [skipped: no digest]", "", str(p))]
    d = json.loads(p.read_text())
    rows = {float(k): v for k, v in d.items()
            if isinstance(v, dict) and "n" in v}
    if not rows:
        return [("transfer-cost prose  [skipped: digest has no targets]", "", "")]
    ts = sorted(rows)
    where = "[in sections/05l_transfer_cost.tex]"
    out = []
    for key, label in (("psf_over_base", "PSF / C-T"),
                       ("linear_over_base", "IDW / C-T"),
                       ("deconv_over_control", "PSF / IDW")):
        v = [rows[t][key] for t in ts]
        out.append((f"transfer cost, {label} range {where}",
                    f"${min(v):.2f}$ to ${max(v):.2f}$", "paired medians over targets"))
    both = [rows[t][k] for t in ts for k in ("psf_over_base", "linear_over_base")]
    # Round outwards, not to nearest: the two-decimal table shows 2.45, and a
    # one-decimal "2.4" beside it reads as a contradiction (review R1 of v6).
    lo = math.floor(min(both) * 10) / 10
    hi = math.ceil(max(both) * 10) / 10
    out.append(("transfer cost, either transfer / C-T, one-decimal range (abstract, intro, conclusion)",
                f"${lo:.1f}$ to ${hi:.1f}$",
                "min and max over both transfers and all targets, rounded outwards"))
    v = [rows[t]["deconv_over_control"] for t in ts]
    pct = f"${100*(min(v)-1):.0f}$--${100*(max(v)-1):.0f}\\%$"
    for scope in (where, "[in sections/07_conclusion.tex]",
                  "[in sections/08_appendices.tex]"):
        out.append((f"transfer cost, PSF/IDW differential as percent {scope}",
                    pct, "(ratio - 1), min and max over targets"))
    v = [rows[t]["psf_over_base"] for t in ts]
    out.append(("transfer cost, PSF / C-T range [in sections/08_appendices.tex]",
                f"${min(v):.2f}$ to ${max(v):.2f}$", "Appendix E item"))
    for key, label in (("clough_tocher_over_utm", "C-T"), ("psf_over_utm", "PSF"),
                       ("linear_over_utm", "IDW")):
        a, b = rows[ts[-2]][key], rows[ts[-1]][key]
        out.append((f"transfer cost, {label} / UTM at the two loosest targets {where}",
                    f"${a:.2f}$ and ${b:.2f}$", "gamma_tile on the UTM side"))
    out.append((f"transfer cost, n at the loosest target {where}",
                f"falls to {rows[ts[-1]]['n']}", "oracle-bracketed on all four sweeps"))
    comp = _csv("compactability.csv")
    if comp is not None:
        col = f"t{ts[-1]:g}"
        names = {"clough_tocher": "C-T", "psf": "PSF", "linear": "IDW",
                 "native UTM": "native"}
        for _, r in comp.iterrows():
            if col in r and np.isfinite(r[col]):
                out.append((f"compactability at loosest, {names.get(r.transfer, r.transfer)} {where}",
                            f"${100*r[col]:.1f}\\%$", "size / own lossless size"))
    return out


def censoring_checks():
    """Everything Section IV-B, Section V-C and Appendix B-D say about censoring.

    The exclusion counts and the bracketed-only medians come from the same
    two sweeps as the tables, through censoring_sensitivity.py, so a ladder
    extended (or not) on one side cannot drift away from the prose.
    """
    out = []
    try:
        doc = json.loads((cfg.RESULTS_DIR / "paper_numbers_self.json").read_text())
    except Exception as exc:                              # noqa: BLE001
        return [(f"censoring counts  [skipped: {type(exc).__name__}]", "", str(exc)[:60])]
    per = {}
    for t in doc["targets"]:
        if t["kind"] != "rmse":
            continue
        per[float(t["target"])] = len(t.get("excluded", []))
    ts = sorted(per)
    if ts:
        out.append(("regions excluded at the loosest RMSE target",
                    str(per[ts[-1]]), "oracle codec censored on a grid"))
        out.append(("regions remaining at the loosest RMSE target",
                    str(40 - per[ts[-1]]), "n of table_storage_at_rmse"))
        if len(ts) >= 3:
            out.append(("exclusion counts, three loosest, 'one/three/eight' form",
                        f"one at $10^{{-3}}$, three at $3\\times10^{{-3}}$ and "
                        f"{'eight' if per[ts[-1]] == 8 else per[ts[-1]]} at $10^{{-2}}$",
                        f"{per[ts[-3]]}, {per[ts[-2]]}, {per[ts[-1]]} excluded"))
    d = cfg.RESULTS_DIR / "censoring_digest.json"
    if d.exists():
        cz = json.loads(d.read_text())
        loose = cz.get("rmse:0.01")
        if loose:
            lo, hi = loose.get("median_lo"), loose.get("median_hi")
            if lo and hi:
                out.append(("censoring interval at the loosest target",
                            f"$[{lo:.2f}, {hi:.2f}]$",
                            "median under any reading of the bounds"))
    return out


def same_codec_checks():
    """The per-codec medians Section V-E quotes, from the same digest as its table."""
    try:
        doc = json.loads((cfg.RESULTS_DIR / "paper_numbers_self.json").read_text())
        fc = doc["fixed_codec"]
    except Exception as exc:                              # noqa: BLE001
        return [(f"same-codec medians  [skipped: {type(exc).__name__}]", "", str(exc)[:60])]
    loose = f"{max(cfg.TARGET_RMSE):g}"
    out = []
    for codec in ("jpegxl_lossy", "jpeg2000_lossy", "dwt_zstd", "sz3"):
        t = fc.get(codec, {}).get(loose)
        if t and np.isfinite(t.get("median", np.nan)):
            out.append((f"same-codec median @ loosest, {codec}",
                        _fmt(t["median"]), f"n={t['n']}"))
    z = fc.get("zfp", {})
    meds = [z[k]["median"] for k in z if np.isfinite(z[k].get("median", np.nan))]
    if meds:
        out.append(("same-codec range, zfp (no crossover)",
                    f"${min(meds):.2f}$--${max(meds):.2f}$",
                    "over the five RMSE targets"))
    return out


# Phrases the revision removed. They are checked rather than trusted because
# each was a claim the results do not support: an exact level-19 carrier, a
# variance decomposition of the floor, a bound on the support mismatch, and
# a grid comparison that is really a carrier comparison.
BANNED = [
    ("drop-in carrier", "kept: say what the measurement shows"),
    ("that the same transfer would have produced", "four-child averaging is not a level-19 transfer"),
    ("two thirds of the variance", "use: removes 66% of the transfer-and-read MSE"),
    ("grid-intrinsic", "use: overhead-free carrier/native"),
    ("at equal accuracy", "use: at the same coding RMSE against each transfer's own carrier"),
    ("at the same accuracy", "use: at the same coding RMSE against each transfer's own carrier"),
    ("less than a fifth", "the control does not decompose the residual"),
    ("operational UTM size", "use: the native product under the illustrative overhead"),
]


def wording_checks(tex):
    """Report every retired phrase that has come back into the manuscript."""
    return [(f"retired phrase: {phrase!r}", "ABSENT", why)
            for phrase, why in BANNED if phrase in tex]


def _smoothed_rows(sm, where, tag, _st):
    """One phase of the smoothed-native control."""
    out = []
    ts = sorted(sm.target.unique())
    a = {t: sm[(sm.target == t) & sm.both_bracketed] for t in ts}
    b = {t: sm[(sm.target == t) & sm.hpx_both_bracketed] for t in ts}
    if ts and len(a[ts[0]]) and len(a[ts[-1]]):
        r0 = a[ts[0]].ratio_smoothed_native.median()
        r1 = a[ts[-1]].ratio_smoothed_native.median()
        out.append((f"smoothed/native at the tightest target{tag} {where}",
                    f"${r0:.2f}$", "median, both bracketed"))
        out.append((f"smoothed/native at the loosest target{tag} {where}",
                    f"${r1:.2f}$", ""))
        out.append((f"smoothed/native n at tightest and loosest{tag} {where}",
                    f"${len(a[ts[0]])}$ and ${len(a[ts[-1]])}$", ""))
        if not tag:
            out.append((f"fixed-phase control, native size reduction range {where}",
                        f"${100*(1-r0):.0f}\\%$ to ${100*(1-r1):.0f}\\%$", "1 - ratio"))
    if ts and len(b[ts[0]]) and len(b[ts[-1]]):
        h0 = b[ts[0]].ratio_hpx_smoothed.median()
        h1 = b[ts[-1]].ratio_hpx_smoothed.median()
        out.append((f"HPX/smoothed at the tightest target{tag} {where}",
                    f"${h0:.2f}$", ""))
        out.append((f"HPX/smoothed at the loosest target{tag}"
                    + ("" if tag else " (several places)"),
                    f"{h1:.2f}", "abstract, V-B, V-G, discussion, conclusion"))
        out.append((f"HPX/smoothed n at the loosest target{tag} {where}",
                    f"${len(b[ts[-1]])}$ regions", ""))
        pt = _st.paired_ratio_test(b[ts[-1]].ratio_hpx_smoothed)
        pv = pt["p"]
        if np.isfinite(pv):
            if pv < 1e-4:
                ptxt = "$p<10^{-4}$"
            else:
                e = int(np.floor(np.log10(pv)))
                ptxt = f"$p={pv/10**e:.0f}\\times10^{{{e}}}$"
            out.append((f"HPX/smoothed paired p at the loosest target{tag} {where}",
                        ptxt, "signed-rank on log ratio"))
    return out


def consistency_checks(tex):
    """The figures the second review found stale, recomputed at the source.

    Each of these had a generated table on one side and a hand-typed number on
    the other. The table is authoritative, so the check reads the same data the
    table is built from and asks the prose to agree.
    """
    out = []
    # --- codec wins and mean ranks (Table IX) ------------------------------
    tab = cfg.REPO_ROOT / "tex" / "tables" / "table_codec_ranking.tex"
    if tab.exists():
        for line in tab.read_text(encoding="utf-8").splitlines():
            m = re.match(r"\\texttt\{(\S+?)\} & (\d+) & (\d+) & ([\d.]+)", line)
            if m:
                name = m.group(1).replace("\\_", "_")
                out.append((f"codec wins, {name}", m.group(2), "Table IX"))
                out.append((f"codec mean rank, {name}", f"${m.group(4)}$", "Table IX"))
    else:
        out.append(("codec ranking  [skipped: table_codec_ranking.tex absent]", "", ""))

    # --- synthetic control medians and the control/measured agreement --------
    syn = _csv("support_controls_synthetic.csv")
    if syn is None:
        out.append(("synthetic control  [skipped: support_controls_synthetic.csv "
                    "absent]", "", ""))
    else:
        # The paper's rows are the field with the slope fitted to the imagery
        # (support_controls.py --fitted, beta = 1.9); a run on the exploratory
        # grid beta = 3, 4, 5 writes the same CSV and must not be mistaken for
        # them. Keep only the beta nearest the fitted median, and say so when
        # that beta is not in the file at all.
        # The instrument response must sit in the ESA MSI B04 MTF-at-Nyquist
        # window (0.15-0.30). A run that band-limited the field with the
        # library scale s (about twice the standard deviation) records an MTF
        # near 0.006 and its ratios are meaningless; reject it by the recorded
        # MTF rather than by trusting a number that survived the doubling.
        if "mtf_utm_nyquist" in syn.columns and len(syn):
            mtf = float(syn.mtf_utm_nyquist.median())
            if not (0.12 <= mtf <= 0.33):
                out.append((f"synthetic control  [skipped: MTF at the native "
                            f"Nyquist is {mtf:.3f}, outside the ESA B04 window "
                            "0.15-0.30 -- the PSF width was doubled (library "
                            "scale used as sigma); rerun support_controls.py "
                            "--control synthetic --fitted with the hardened "
                            "psf_sigma_from_fwhm]", "", ""))
                syn = syn.iloc[0:0]
        spec = _csv("support_controls_spectrum.csv")
        if len(syn) and spec is not None and "beta" in syn.columns:
            b_fit = round(float(spec.beta_fit.median()), 1)
            near = syn.beta.sub(b_fit).abs()
            if near.min() > 0.15:
                out.append((f"synthetic control  [skipped: no rows at the fitted "
                            f"beta {b_fit:g}; file holds beta = "
                            f"{sorted(syn.beta.unique())} -- run "
                            "support_controls.py --control synthetic --fitted]", "", ""))
                syn = syn.iloc[0:0]
            else:
                syn = syn[near <= 0.15]
        # Revision 5 adds rows where the field is carried to the cells by a
        # transfer (column `transfer`); the geometry-only statements of the
        # text are about the direct rows.
        if "transfer" in syn.columns and len(syn):
            syn = syn[syn.transfer.fillna("none").astype(str) == "none"]
        syn_all = syn
        syn = syn[syn.both_bracketed] if len(syn) else syn
        med = syn.groupby("target").ratio.median() if len(syn) else pd.Series(dtype=float)
        # n at the two loosest targets, quoted as "partly censored (n = a and b)"
        if len(syn):
            ns = syn.groupby("target").size()
            ts_ = sorted(ns.index)
            if len(ts_) >= 2:
                out.append(("synthetic control n at the two loosest targets [in sections/05e_support.tex]",
                            f"$n={ns[ts_[-2]]}$ and ${ns[ts_[-1]]}$", "both bracketed"))
            out.append(("synthetic control at the loosest target under gamma_tile [in sections/05e_support.tex]",
                        f"${med[ts_[-1]]/cfg.GAMMA_UTM:.2f}$", "median / 1.30"))
        gs = {float(k): float(v) for k, v in load("self").get("gamma_star", {}).items()}
        agree = []
        ts = sorted(med.index)
        for t, m in med.items():
            # the prose quotes the endpoints of the decline; interior values
            # are reported, as for the self ratios, but never asserted
            interior = t not in (ts[0], ts[-1])
            out.append((f"synthetic control ratio @ {t:g}"
                        + ("  [range interior, not asserted]" if interior else ""),
                        _fmt(m), "Table VIII"))
            if t in gs:
                agree.append(100 * abs(m - gs[t]) / gs[t])
        if agree:
            rng = f"${min(agree):.0f}$--${max(agree):.0f}\\%"
            n_found = tex.count(rng)
            out.append((f"control/measured agreement range (x{n_found} in text; "
                        f"needs every occurrence)", rng,
                        "min/max over targets of |control - gamma_star|/gamma_star"))
            # the same paragraph also states the excess in "X% to Y%" form
            # (the control sits above the measurement); a stale 7--17 % had
            # survived priority 0 in that form, so it is asserted separately
            rng2 = f"${min(agree):.0f}\\%$ to ${max(agree):.0f}\\%$"
            out.append(("control excess over measured, 'X% to Y%' form", rng2,
                        "same quantity as the agreement range"))

    # --- the synthetic field carried through the real transfer (revision 5) --
    syn2 = _csv("support_controls_synthetic.csv")
    if syn2 is not None and "transfer" in syn2.columns and len(syn2):
        where = "[in sections/05e_support.tex]"
        ok_mtf = ("mtf_utm_nyquist" not in syn2.columns
                  or 0.12 <= float(syn2.mtf_utm_nyquist.median()) <= 0.33)
        spec = _csv("support_controls_spectrum.csv")
        if ok_mtf and spec is not None:
            b_fit = round(float(spec.beta_fit.median()), 1)
            syn2 = syn2[syn2.beta.sub(b_fit).abs() <= 0.15]
            via = syn2[(syn2.transfer.fillna("none").astype(str) != "none") & syn2.both_bracketed]
            direct = syn2[(syn2.transfer.fillna("none").astype(str) == "none") & syn2.both_bracketed]
            if len(via) and len(direct):
                mv = via.groupby("target").ratio.median()
                md_ = direct.groupby("target").ratio.median()
                ts2 = sorted(mv.index)
                out.append((f"synthetic via transfer at the tightest target {where}",
                            f"${mv[ts2[0]]:.2f}$", "both bracketed"))
                out.append((f"synthetic via transfer at the loosest target {where}",
                            f"${mv[ts2[-1]]:.2f}$", ""))
                out.append((f"synthetic via transfer n at the loosest target {where}",
                            f"${len(via[via.target == ts2[-1]])}$", ""))
                # what the transfer removes at the loosest target, as percent
                out.append((f"transfer's share of the synthetic ratio at the loosest target {where}",
                            f"${100*(1-mv[ts2[-1]]/md_[ts2[-1]]):.0f}\\%$", "1 - via/direct"))
                tr = via.groupby("scene").transfer_rmse_rel.first().dropna()
                if len(tr):
                    out.append((f"transfer RMSE on the analytic field, % of std {where}",
                                f"${100*tr.median():.1f}\\%$", "median over regions"))
                # agreement of the via-transfer control with gamma_star, min/max
                gs2 = {float(k): float(v) for k, v in load("self").get("gamma_star", {}).items()}
                ag = [100*abs(mv[t]-gs2[t])/gs2[t] for t in ts2 if t in gs2]
                if ag:
                    out.append((f"via-transfer control / measured agreement range {where}",
                                f"${min(ag):.0f}$--${max(ag):.0f}\\%$", "|via - gamma_star|/gamma_star"))
                # agreement of HPX/smoothed-native with the direct control
                smx = _csv("support_controls_smoothed.csv")
                if smx is not None and len(smx):
                    if {"shift_x", "shift_y"} <= set(smx.columns):
                        smx = smx[np.isclose(smx.shift_x, 0.5) & np.isclose(smx.shift_y, 0.5)]
                    hs = smx[smx.hpx_both_bracketed].groupby("target").ratio_hpx_smoothed.median()
                    ag2 = [100*abs(hs[t]-md_[t])/md_[t] for t in ts2 if t in hs.index]
                    if ag2:
                        out.append(("HPX/smoothed vs direct synthetic, largest gap, percent [in sections/05e_support.tex]",
                                    f"${max(ag2):.0f}\\%$", "max over targets of |hpx/smoothed - direct|/direct"))

    # --- the smoothed-native control (revision 5): the transfer's smoothing
    # given to the native side, Section V-G ---------------------------------
    sm = _csv("support_controls_smoothed.csv")
    if sm is None:
        out.append(("smoothed-native control  [skipped: support_controls_smoothed.csv "
                    "absent]", "", ""))
    else:
        import paper_stats as _st
        where = "[in sections/05e_support.tex]"
        # Revision 6 (run 2.3) added a second phase: the native side is read
        # at quarter-pixel as well as at half-pixel offset. The half-pixel
        # block keeps the labels the published text was checked against; the
        # other phases get their own, so both are asserted.
        phases = [(0.5, "")]
        if {"shift_x", "shift_y"} <= set(sm.columns):
            other = sorted({(round(float(x), 3), round(float(y), 3))
                            for x, y in zip(sm.shift_x, sm.shift_y)}
                           - {(0.5, 0.5)})
            phases += [(x, f", phase {x:g}") for x, y in other]
        sm_all = sm
        for _phase, _tag in phases:
            if {"shift_x", "shift_y"} <= set(sm_all.columns):
                sm = sm_all[np.isclose(sm_all.shift_x, _phase)
                            & np.isclose(sm_all.shift_y, _phase)]
            else:
                sm = sm_all
            if not len(sm):
                continue
            out += _smoothed_rows(sm, where, _tag, _st)


    # --- the level-19 carrier under compression (revision 5, Section V-K) --
    p19 = cfg.RESULTS_DIR / "level19_rd.json"
    if not p19.exists():
        out.append(("level-19 compressed  [skipped: level19_rd.json absent]", "", ""))
    else:
        d19 = json.loads(p19.read_text())
        where = "[in sections/05h_level_choice.tex]"
        c = d19.get("carriers", {}).get("direct")
        if c:
            ts = sorted(float(k) for k in c["targets"])
            e0, e1 = c["targets"][f"{ts[0]:g}"], c["targets"][f"{ts[-1]:g}"]
            if d19.get("raw_ratio_l19_utm"):
                out.append((f"level-19 raw ratio {where}", f"${d19['raw_ratio_l19_utm']:.2f}$", ""))
            out.append((f"L19/UTM at the tightest target {where}",
                        f"${e0['l19_over_utm']['median']:.2f}$", "direct carrier"))
            out.append((f"L19/UTM at the loosest target {where}",
                        f"${e1['l19_over_utm']['median']:.2f}$", ""))
            out.append((f"L19/L20 at the tightest target {where}",
                        f"${e0['l19_over_l20']['median']:.2f}$", ""))
            out.append((f"L19/L20 at the loosest target {where}",
                        f"${e1['l19_over_l20']['median']:.2f}$", ""))
            out.append((f"L20/L19 at the loosest target {where}",
                        f"${1/e1['l19_over_l20']['median']:.2f}$", "inverse"))
            out.append((f"L20/L19 at the tightest target {where}",
                        f"${1/e0['l19_over_l20']['median']:.1f}$", "inverse, one decimal"))
            wf = [c["targets"][f"{t:g}"]["l19_over_utm_win_fraction"] for t in ts]
            wf = [w for w in wf if w is not None]
            if wf:
                out.append((f"L19 win-fraction range {where}",
                            f"${100*min(wf):.0f}\\%$ to ${100*max(wf):.0f}\\%$", ""))
        da = d19.get("direct_over_aggregate_size")
        if da:
            out.append((f"aggregate cheaper than direct, percent {where}",
                        f"${100*(1-1/da['median']):.0f}\\%$", "1 - 1/median"))
            out.append((f"direct/aggregate n {where}", f"${da['n']}$", ""))
        # Revision 6 (run 2.4): why the aggregate's L19/UTM column is flat.
        av = d19.get("aggregate_vs_direct") or {}
        if av:
            ts_av = sorted(float(k) for k in av)
            lo, hi = av[f"{ts_av[0]:g}"], av[f"{ts_av[-1]:g}"]
            out.append((f"aggregate/direct at the tightest target {where}",
                        f"${lo['aggregate_over_direct_median']:.2f}$",
                        "paired, regions eligible for both carriers"))
            out.append((f"aggregate/direct at the loosest target {where}",
                        f"${hi['aggregate_over_direct_median']:.2f}$", ""))
            out.append((f"aggregate compressed / raw at the loosest target {where}",
                        f"${100*hi['compact_aggregate']:.1f}\\%$",
                        "compare with the native crop"))
            out.append((f"native crop compressed / raw at the loosest target {where}",
                        f"${100*hi['compact_utm']:.1f}\\%$", ""))
            out.append((f"direct carrier compressed / raw at the loosest target {where}",
                        f"${100*hi['compact_direct']:.1f}\\%$", ""))

    # --- tiling overhead from the product geometry (revision 5, IV-B) -------
    pg = cfg.RESULTS_DIR / "granule_overhead.json"
    if not pg.exists():
        out.append(("granule overhead  [skipped: granule_overhead.json absent]", "", ""))
    else:
        g = json.loads(pg.read_text())
        where = "[in sections/04_design.tex]"
        out.append((f"tile overlap overhead {where}", f"$(109.8/100)^2={g['overlap']:.2f}$",
                    "granule_overhead.py"))
        sw = g.get("swath") or {}
        if sw:
            # Revision 6 (R1 3.6): a reader will redo the swath arithmetic,
            # so every step of it is asserted, not only the product.
            out.append((f"swath, mean tiles across track {where}",
                        f"${sw['n_tiles_mean']:.3f}$", "averaged over the swath edge"))
            out.append((f"swath, across-track factor {where}",
                        f"${sw['across_track']:.2f}$",
                        "stored raster / swath width"))
            out.append((f"swath, along-track factor {where}",
                        f"${sw['along_track']:.3f}$", "overlap only"))
            out.append((f"swath, no-data share {where}",
                        f"${100*sw['nodata_fraction_mean']:.0f}\\%$", ""))
            out.append((f"swath, raw-byte overhead {where}",
                        f"${sw['total']:.2f}$", "across x along"))
        out.append((f"tiled-swath raw overhead {where}", f"${g['swath']['total']:.2f}$", ""))
        out.append((f"no-data share of the tiled swath {where}",
                    f"${100*g['swath']['nodata_fraction_mean']:.0f}\\%$", ""))
        out.append(("overlap overhead as percent (several places)",
                    f"${100*(g['overlap']-1):.0f}\\%$", "V-B, limitations"))

    # --- block-size and odd-side bias (revision 5, IV-B) -------------------
    pb_ = cfg.RESULTS_DIR / "block_size_bias.json"
    if not pb_.exists():
        out.append(("block-size bias  [skipped: block_size_bias.json absent]", "", ""))
    else:
        b = json.loads(pb_.read_text())
        where = "[in sections/04_design.tex]"
        x = b.get("odd_side_excess_lossless")
        if x is not None:
            out.append((f"odd-side excess, lossless, percent {where}",
                        f"${100*abs(x-1):.1f}\\%$", "|excess - 1|"))
        r = b["ratios_bytes_per_px"]["lossless_bytes_per_px"].get("utm_640_over_utm_odd")
        if r:
            out.append((f"padding to 640, lossless change {where}",
                        f"${100*abs(r['median']-1):.1f}\\%$", ""))

    # --- the estimand's depth, as the table caption states it ----------------
    re_csv = _csv("resampling_error.csv")
    if re_csv is not None and "children" in re_csv:
        depths = sorted(set(int(c) for c in re_csv.children.dropna().unique()))
        if len(depths) == 1:
            words = {1: "one", 4: "four", 16: "sixteen", 64: "sixty-four"}
            n_ch = depths[0]
            lvl = 20 + int(round(np.log(n_ch) / np.log(4))) if n_ch > 1 else 20
            out.append(("estimand depth in the text",
                        f"{words.get(n_ch, n_ch)} level-{lvl} children",
                        "must match table_resampling_error caption"))

    # --- the interpolating reader (Table reader_floor) -------------------------
    rf = _csv("reader_floor.csv")
    if rf is None:
        out.append(("reader floor  [skipped: results/reader_floor.csv absent]", "",
                    "run reader_floor.py"))
    else:
        ok = rf[rf.error.fillna("").astype(str).str.len() == 0]
        ct = ok[ok.kernel == "clough_tocher"].pivot_table(
            index="scene", columns="reader", values="floor")
        if "nearest" in ct and "bilinear" in ct:
            r = (ct.bilinear / ct.nearest).dropna()
            share = 1.0 - r.median() ** 2
            b = ok[(ok.kernel == "clough_tocher") & (ok.reader == "bilinear")]
            out += [("bilinear / nearest floor, C-T tiles", _fmt(r.median(), 3),
                     f"paired median over {len(r)} regions"),
                    ("read share of floor variance, %", f"{100 * share:.0f}\\%",
                     "1 - ratio^2; the text says 'two thirds' when this is 60-70"),
                    ("bilinear floor, % of std", _fmt(100 * b.floor_frac_of_std.median(), 1) + "\\%", ""),
                    ("bilinear floor, absolute",
                     f"{b.floor.median():.1e}".replace("e-0", "\\times10^{-").replace("e-", "\\times10^{-") + "}",
                     "as $x.x\\times10^{-3}$")]

    # --- the transfer ablation (Table kernel_ablation) ---------------------------
    ka = _csv("kernel_ablation.csv")
    if ka is None:
        out.append(("kernel ablation  [skipped: results/kernel_ablation.csv absent]",
                    "", "run kernel_ablation.py"))
    else:
        # Same admission rule as merged_tables.table_transfer block (c): the
        # window audit's MISMATCH regions are left out (Lake Geneva).
        wr = _csv("window_radiometry.csv")
        if wr is not None and "status" in wr:
            bad = set(wr[wr.status.astype(str).str.upper() == "MISMATCH"].scene)
            ka = ka[~ka.scene.isin(bad)]
        med = ka.groupby("kernel").agg(tr=("transfer_frac_of_std", "median"),
                                       oc=("logpsd_top_octave", "median"))
        for k in ("clough_tocher", "nearest", "linear", "psf"):
            if k in med.index:
                out.append((f"transfer error vs analytic field, {k}, % of std",
                            _fmt(100 * med.loc[k, "tr"]) + "\\%", "Table kernel_ablation"))
        for k in ("psf", "linear"):
            if k in med.index:
                out.append((f"top-octave power ratio, {k}", _fmt(med.loc[k, "oc"]), ""))
        if {"clough_tocher", "linear"} <= set(med.index):
            out.append(("IDW / C-T transfer error", _fmt(med.loc["linear", "tr"] / med.loc["clough_tocher", "tr"], 1),
                        "the text says 'x times less accurate'"))

    # --- the common reference under both readers (Section V-H) ------------------
    for reader, fname in (("nearest", "paper_numbers_common.json"),
                          ("bilinear", "paper_numbers_common_bilinear.json")):
        pth = cfg.RESULTS_DIR / fname
        if not pth.exists():
            out.append((f"common reference, {reader}  [skipped: {fname} absent]", "", ""))
            continue
        doc = json.loads(pth.read_text(encoding="utf-8"))
        for tgt in doc["targets"]:
            r = list(tgt["ratios"].values())
            br = float(tgt["bracketed_fraction"])
            kind, tv = tgt["kind"], float(tgt["target"])
            if kind == "budget" and tv in (0.02, 0.05, 0.1) and r:
                m = float(np.median(r))
                out.append((f"common ratio @ budget {100*tv:g}%, {reader}",
                            _fmt(m) if m < 10 else _fmt(m, 1), f"n={len(r)}"))
            if kind == "rmse" and tv in (3e-4, 1e-3, 3e-3, 1e-2):
                out.append((f"bracketed @ RMSE {tv:g}, {reader}", f"{100*br:.0f}\\%",
                            "censoring figure quoted in V-H"))
            if kind == "rmse" and tv == 1e-2 and br >= 0.75 and r:
                g1 = (tgt.get("test_gamma1") or {}).get("median")
                out += [(f"common ratio @ RMSE 0.01, {reader}", _fmt(float(np.median(r))), "operational"),
                        (f"common p_Holm @ RMSE 0.01, {reader}", f"{tgt['p_holm']:.2f}", ""),
                        (f"intrinsic ratio @ RMSE 0.01, {reader}", _fmt(float(g1)), "gamma_star under this reference")]
        cr = doc.get("codec_ranking", {})
        if isinstance(cr, dict):
            for c in ("jpegxl_lossy", "sz3"):
                v = cr.get(c)
                if isinstance(v, dict) and "wins" in v:
                    out.append((f"common codec wins, {c}, {reader}", str(int(v["wins"])), ""))
                    if "mean_rank" in v:
                        out.append((f"common codec mean rank, {c}, {reader}", _fmt(float(v["mean_rank"])), ""))
        lat = doc.get("latitude", {}).get("budget:0.02", {})
        evc = lat.get("equatorial_vs_cap", {})
        if evc:
            for k, lab in (("belt", "belt"), ("cap", "cap")):
                if k in evc:
                    out.append((f"belt/cap ratio @ budget 2%, {lab}, {reader}", _fmt(float(evc[k])), ""))
            if "p_holm" in evc:
                out.append((f"belt vs cap Holm p @ budget 2%, {reader}", f"{evc['p_holm']:.2f}" if evc["p_holm"] >= 0.01 else f"{evc['p_holm']:.3f}", ""))
    fl = _csv("common_reference_floor__bilinear.csv")
    if fl is not None:
        out.append(("bilinear common floor, median",
                    f"{fl.representation_floor_rmse.median():.1e}".replace("e-0", "\\times10^{-") + "}", ""))

    # --- the floor is never called the grid's own -----------------------------
    prose = re.sub(r"tab:representation_floor", "", tex)
    n_bad = len(re.findall(r"representation floor", prose))
    out.append((f"'representation floor' in prose (found {n_bad}; must be 0)  "
                "[count, not asserted]" if n_bad == 0 else
                f"'representation floor' in prose (found {n_bad}; must be 0)",
                "transfer-and-read floor", "one term throughout"))
    return out


def _csv(name):
    p = cfg.RESULTS_DIR / name
    return pd.read_csv(p) if p.exists() else None


def level_and_reader_checks():
    """Figures quoted in Sections V-F (reader's kernel), V-G (level choice)
    and III-E (quadtree test), recomputed from the runs that produced them.

    Each block is skipped with a note when its result file is absent, so the
    check still runs on a checkout without those runs; it never passes a
    stale figure when the file is present.
    """
    out = []

    # --- level choice ---------------------------------------------------------
    lc = _csv("level_choice.csv")
    if lc is None:
        out.append(("level choice  [skipped: results/level_choice.csv absent]",
                    "", "run level_choice.py"))
    else:
        f20 = lc[lc.depth == 0].set_index("scene")
        f19 = lc[lc.depth == 1].set_index("scene")
        n = f19.index.size
        quad = np.sqrt(f20.floor ** 2 + f19.discarded ** 2)
        excess = f19.floor ** 2 - f20.floor ** 2
        out += [
            ("level 19 discards, % of variance",
             _fmt(100 * f19.discarded_var_frac.median()), f"median over {n}"),
            ("level 19 discards, % of std",
             _fmt(100 * f19.discarded_frac_of_std.median(), 1) + "\\%", ""),
            ("discarded / level-20 floor",
             _fmt((f19.discarded / f20.floor).median()), "paired median"),
            ("level-19 floor / level-20 floor",
             _fmt((f19.floor / f20.floor).median()), "paired median"),
            ("level-19 floor / quadrature",
             _fmt((f19.floor / quad).median()), "paired median"),
            ("discarded^2 / excess variance, %",
             f"{100 * (f19.discarded ** 2 / excess).median():.0f}\\%", ""),
        ]
        b = _csv("hpx_codec_bench.csv")
        if b is not None:
            b = b.dropna(subset=["rmse"])
            b = b[b.rmse > 0].copy()
            b["f20"] = b.scene.map(f20.floor)
            b["f19"] = b.scene.map(f19.floor)
            b = b.dropna(subset=["f20", "f19"])
            out += [
                ("sweep points below level-20 floor, %",
                 _fmt(100 * (b.rmse < b.f20).mean(), 1) + "\\%",
                 f"{len(b)} points"),
                ("sweep points below level-19 floor, %",
                 _fmt(100 * (b.rmse < b.f19).mean(), 1) + "\\%", ""),
            ]

    # --- the reader's kernel, nearest read on the Clough-Tocher tiles ---------
    rd = _csv("interscale_reader.csv")
    if rd is None:
        out.append(("reader kernel  [skipped: results/interscale_reader.csv "
                    "absent]", "", "run interscale.py"))
    else:
        if "error" in rd:
            rd = rd[rd["error"].fillna("").astype(str).str.len() == 0]
        g = rd[rd.reader == "nearest"]
        out += [
            ("nearest read: cells reached, %",
             _fmt(100 * g.support_fraction.median(), 1) + "\\%", ""),
            ("nearest read: error energy in kernel, %",
             _fmt(100 * g.kernel_energy_frac.median(), 1) + "\\%", ""),
            ("nearest read: transmission",
             _fmt(g.transmission.median()), "must read 1.00"),
            ("nearest read: rate saved by product target, %",
             _fmt(100 * g.gain_product_target.median(), 1) + "\\%", ""),
        ]

    # --- the quadtree test of Section III-E ------------------------------------
    for name, lab in (("interscale.csv", "Haar"), ("interscale_cdf53.csv", "CDF 5/3")):
        d = _csv(name)
        if d is None:
            out.append((f"quadtree bound, {lab}  [skipped: results/{name} absent]",
                        "", ""))
            continue
        h = d.dropna(subset=["ref_bytes"])
        out.append((f"quadtree bound / codec bytes, {lab}",
                    _fmt((h.bytes_both / h.ref_bytes).median()), "median"))
        if lab == "Haar":
            out += [("parent context saves, %",
                     _fmt(100 * h.gain_parent.median(), 1) + "\\%", ""),
                    ("neighbour context saves, %",
                     _fmt(100 * h.gain_neighbour.median(), 1) + "\\%", "")]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    if not TEX.exists():
        raise SystemExit(f"{TEX} not found")
    tex = read_manuscript(TEX)
    rows, survivors = checks(tex)

    bad = []
    for label, value, note in rows:
        if ("not asserted" in label or "[skipped" in label
                or label.startswith("number of latitude")):
            # a count, not a figure in the text: reported, never asserted
            continue
        if "needs every occurrence" in label:
            if tex.count(value) < 3:
                bad.append((label, value, note))
            continue
        if "[count, not asserted]" in label:
            continue
        m = re.search(r"\[in ([^\]]+)\]", label)
        if m:
            sub = TEX.parent / m.group(1)
            hay = sub.read_text(encoding="utf-8") if sub.exists() else ""
            if value not in hay:
                bad.append((label, value, note))
            continue
        if value not in tex:
            bad.append((label, value, note))

    if not args.quiet:
        print(f"manuscript: {TEX}")
        print(f"{len(rows)} figure(s) checked against the current results\n")
        for label, value, note in rows:
            mark = "MISSING" if (label, value, note) in bad else "ok     "
            print(f"  {mark}  {label:<46s} {value:>8s}   ({note})")
        print()
        if bad:
            print(f"{len(bad)} figure(s) from the current results do not appear in "
                  f"the manuscript. Either the prose is stale, or it rounds them "
                  f"differently -- check each before dismissing it:")
            for label, value, note in bad:
                print(f"  - {label}: expected to find {value}")
        else:
            print("every checked figure from the current results appears in the "
                  "manuscript.")
        print(f"\nlatitude associations surviving Holm: {len(survivors)} "
              f"({', '.join(survivors) or 'none'}) -- the text must say the same "
              f"number.")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
