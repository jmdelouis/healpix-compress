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
import re
import sys

import numpy as np
import pandas as pd

import paper_config as cfg

TEX = cfg.REPO_ROOT / "tex" / "tgrs_healpix_compress.tex"


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

    # --- how many latitude associations survive Holm -------------------------
    lat = self_.get("latitude", {})
    survivors = sorted(k for k, v in lat.items()
                       if min(v.get("spearman", {}).get("p_holm", 1.0),
                              v.get("equatorial_vs_cap", {}).get("p_holm", 1.0)) < 0.05)
    out.append(("number of latitude associations surviving Holm",
                str(len(survivors)),
                "survivors: " + (", ".join(survivors) or "none")))
    out += level_and_reader_checks()
    out += consistency_checks(tex)
    return out, survivors


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
        syn = syn[syn.both_bracketed]
        med = syn.groupby("target").ratio.median()
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
        ps = ok[ok.kernel == "psf"].pivot_table(index="scene", columns="reader", values="floor")
        if "psf" in ps and "bilinear" in ps:
            out.append(("psf / bilinear floor, PSF tiles",
                        _fmt((ps.psf / ps.bilinear).dropna().median(), 3), "paired median"))

    # --- the transfer ablation (Table kernel_ablation) ---------------------------
    ka = _csv("kernel_ablation.csv")
    if ka is None:
        out.append(("kernel ablation  [skipped: results/kernel_ablation.csv absent]",
                    "", "run kernel_ablation.py"))
    else:
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

    # --- the matched reader: PSF operator on the PSF tiles --------------------
    rp = _csv("interscale_reader_psf.csv")
    if rp is None:
        out.append(("matched PSF reader  [skipped: "
                    "results/interscale_reader_psf.csv absent]", "",
                    "run interscale.py --kernel psf"))
    else:
        if "error" in rp:
            rp = rp[rp["error"].fillna("").astype(str).str.len() == 0]
        g = rp[rp.reader == "psf"]
        piv = rp.pivot_table(index="scene", columns="reader",
                             values="reader_floor", aggfunc="median")
        q = pd.qcut(np.log10(g.rmse_cell), 4, labels=list("abcd"))
        gq = g.groupby(q, observed=True).gain_product_target.median()
        out += [
            ("PSF reader floor / nearest floor, PSF tiles",
             _fmt((piv.psf / piv.nearest).median(), 3),
             f"lower on {(piv.psf < piv.nearest).sum()}/{len(piv)} regions"),
            ("PSF reader: transmission", _fmt(g.transmission.median()), ""),
            ("PSF reader: rate saved by product target, %",
             f"{100 * g.gain_product_target.median():.0f}\\%", "median"),
            ("PSF reader: saved at tightest quartile, %",
             _fmt(100 * gq.iloc[0], 1) + "\\%", ""),
            ("PSF reader: saved at loosest quartile, %",
             f"{100 * gq.iloc[-1]:.0f}\\%", ""),
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
    tex = TEX.read_text(encoding="utf-8")
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
