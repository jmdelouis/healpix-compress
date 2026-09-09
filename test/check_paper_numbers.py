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
    return out, survivors


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
        if "not asserted" in label or label.startswith("number of latitude"):
            # a count, not a figure in the text: reported, never asserted
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
