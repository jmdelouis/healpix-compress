"""Rebuild every LaTeX table in tex/tables from the results already on disk.

Why this exists
---------------
The tables are generated, never hand-edited, so a change to how a table is
*presented* -- a caption reworded, a column renamed, the note moved under the
table -- would otherwise mean re-running the measurement that produced it.
Some of those runs take hours (the per-kernel codec sweeps) and none of them
would change a single number.

This script therefore calls each generator's own table function on the CSV it
already wrote, so the table text comes from the same source of truth as a full
run, with no measurement repeated. It rebuilds what the results directory can
support and names, per table, the command to run for the rest.

It is not a substitute for the generators: a table whose results file is
missing is reported, not invented.

Usage
-----
    python rebuild_tables.py              # every table it can rebuild
    python rebuild_tables.py --list       # say what is available, write nothing
"""

from __future__ import annotations

import argparse
import traceback

import pandas as pd

import paper_config as cfg

# table file -> (results file(s) it needs, the command that produces them)
NEEDS = {
    "table_transfer_cost / table_compactability":
        ("superres_cost_*.csv", "python run_codec_bench.py --grid hpx "
         "--resampler psf  (and --resampler linear), then "
         "python superres_cost.py --control linear --tex"),
    "table_alpha_budget": ("resampling_error.csv", "python alpha_budget.py"),
    "table_cross_term": ("cross_term.csv", "python cross_term.py"),
    "table_support_controls": ("support_controls_*.csv",
                               "python support_controls.py --control synthetic"),
    "table_resampling_error": ("resampling_error.csv",
                               "python resampling_error.py"),
    "table_operator_tests": ("operator_tests.csv", "python operator_tests.py"),
    "table_byte_audit": ("codec_audit.csv", "python codec_audit.py"),
}


def _csv(name):
    p = cfg.RESULTS_DIR / name
    return pd.read_csv(p) if p.exists() else None


def _write(name, text):
    cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
    out = cfg.TEX_TABLES / f"{name}.tex"
    out.write_text(text, encoding="utf-8")
    print(f"  wrote {out.name}")


def rebuild(dry=False):
    done, missing = [], []

    def step(label, fn, needs=()):
        absent = [n for n in needs if not (cfg.RESULTS_DIR / n).exists()]
        if absent:
            missing.append((label, ", ".join(absent)))
            return
        if dry:
            done.append(label)
            return
        try:
            fn()
            done.append(label)
        except RuntimeError as exc:      # a missing input, already explained
            missing.append((label, str(exc)))
        except Exception as exc:                          # noqa: BLE001
            missing.append((label, f"{type(exc).__name__}: {exc}"))
            traceback.print_exc()

    # --- the six tables update_paper_numbers owns, per reference ------------
    def _upn(argv):
        import update_paper_numbers as up
        return lambda: up.main(argv)

    step("update_paper_numbers (self)", _upn(["--metric", "self"]),
         ["hpx_codec_bench.csv", "utm_codec_bench.csv"])
    step("update_paper_numbers (common, nearest)", _upn(["--metric", "common"]),
         ["hpx_codec_bench.csv", "common_reference_floor.csv"])
    step("update_paper_numbers (common, bilinear)",
         _upn(["--metric", "common", "--cr-mode", "bilinear"]),
         ["hpx_codec_bench__bilinear.csv"])

    def _gamma():
        import gamma_analysis as ga
        ga.main(["--metric", "self"])
    step("gamma_analysis -> table_gamma_continuous", _gamma,
         ["paper_numbers_self.json"])

    def _kernel():
        import kernel_ablation as ka
        _write("table_kernel_ablation",
               ka.latex_table(_csv("kernel_ablation.csv"), [], False))
    step("kernel_ablation -> table_kernel_ablation", _kernel,
         ["kernel_ablation.csv"])

    def _reader_floor():
        import reader_floor as rf
        rf.latex_table(_csv("reader_floor.csv"),
                       cfg.TEX_TABLES / "table_reader_floor.tex")
    step("reader_floor -> table_reader_floor", _reader_floor,
         ["reader_floor.csv"])

    def _rescore():
        import reader_rescore as rr
        rr.main()
    step("reader_rescore -> table_common_two_readers", _rescore,
         ["paper_numbers_common.json", "paper_numbers_common_bilinear.json"])

    def _resampling():
        import resampling_error as re_
        df = _csv("resampling_error.csv")
        control = _csv("resampling_error_noise0.csv")
        if control is not None and set(control.children) != set(df.children):
            print("  note: the noiseless control sits at a different "
                  "cell-average depth; dropping it, as the generator does")
            control = None
        _write("table_resampling_error", re_.latex_table(df, control))
    step("resampling_error -> table_resampling_error", _resampling,
         ["resampling_error.csv"])

    def _support():
        import support_controls as sc
        frames = {}
        for key in ("synthetic", "common", "shift"):
            d = _csv(f"support_controls_{key}.csv")
            if d is not None and len(d):
                frames[key] = d
        tex = sc.latex_table(frames)
        if tex:
            _write("table_support_controls", tex)
    step("support_controls -> table_support_controls", _support,
         ["support_controls_synthetic.csv"])

    def _cross():
        import cross_term as ct
        _write("table_cross_term", ct.latex_table(_csv("cross_term.csv")))
    step("cross_term -> table_cross_term", _cross, ["cross_term.csv"])

    def _operator():
        import operator_tests as ot
        _write("table_operator_tests", ot.latex_table(_csv("operator_tests.csv")))
    step("operator_tests -> table_operator_tests", _operator,
         ["operator_tests.csv"])

    def _audit():
        import codec_audit as ca
        _write("table_byte_audit", ca.latex_table(_csv("codec_audit.csv")))
    step("codec_audit -> table_byte_audit", _audit, ["codec_audit.csv"])

    def _alpha():
        import alpha_budget as ab
        ab.main([])
    step("alpha_budget -> table_alpha_budget", _alpha,
         ["resampling_error.csv"])

    def _superres():
        import superres_cost as sr
        if sr.latex_tables(baseline="clough_tocher", tested=("psf", "linear"),
                           publish=True) is None:
            raise RuntimeError("the per-kernel sweeps are not all on disk")
    step("superres_cost -> table_transfer_cost, table_compactability",
         _superres, ["hpx_codec_bench.csv"])

    def _headline():
        import headline as hl
        hl.main([])
    step("headline -> table_headline, fig_headline", _headline,
         ["paper_numbers_self.json", "level_choice.csv"])

    return done, missing


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true",
                    help="report what can be rebuilt, write nothing")
    args = ap.parse_args(argv)

    print(f"results: {cfg.RESULTS_DIR}\ntables : {cfg.TEX_TABLES}\n")
    done, missing = rebuild(dry=args.list)

    print(f"\n{len(done)} generator step(s) "
          + ("available" if args.list else "rebuilt") + ":")
    for d in done:
        print(f"  {d}")
    if missing:
        print(f"\n{len(missing)} step(s) could not run:")
        for label, why in missing:
            print(f"  {label}: {why}")
        print("\nThe measurement behind each of those is produced by:")
        for table, (needs, cmd) in NEEDS.items():
            print(f"  {table}\n      needs {needs}\n      {cmd}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
