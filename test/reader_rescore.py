"""The common-reference verdict under two readers (revision 3, item 3.2).

Section V-H scores both grids against the same Sentinel-2 observation, which
requires reading the decoded HEALPix tile on the native pixels. The paper's
sweep did that with the nearest-cell read; `reader_floor.py` then measured
that an interpolating read removes two thirds of that read's floor, so the
sweep was re-scored through it (`run_codec_bench.py --common-reference
--cr-mode bilinear`, which writes reader-suffixed files beside the published
ones). This puts the two summaries side by side, target by target, from the
two `paper_numbers_common*.json` files `update_paper_numbers.py` wrote, and
emits the one table the paper needs.

Nothing is recomputed: the medians, censoring and Holm-adjusted p-values are
the ones each summary already carries, so the table cannot disagree with the
per-reader tables it condenses.

Usage
-----
    python reader_rescore.py
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

import paper_config as cfg

READERS = (("nearest", "paper_numbers_common.json"),
           ("bilinear", "paper_numbers_common_bilinear.json"))
MIN_BRACKETED = getattr(cfg, "MIN_BRACKETED", 0.75)


def load(name):
    p = cfg.RESULTS_DIR / name
    if not p.exists():
        raise SystemExit(f"{p} is missing; run update_paper_numbers.py "
                         f"--metric common{' --cr-mode bilinear' if 'bilinear' in name else ''}")
    return json.loads(p.read_text(encoding="utf-8"))


def rows(doc, reader):
    out = []
    for t in doc["targets"]:
        r = list(t["ratios"].values())
        test, g1 = t.get("test") or {}, t.get("test_gamma1") or {}
        out.append(dict(reader=reader, kind=t["kind"], target=float(t["target"]),
                        n=len(r), median=float(np.median(r)) if r else np.nan,
                        p_holm=t.get("p_holm"), bracketed=float(t["bracketed_fraction"]),
                        reported=float(t["bracketed_fraction"]) >= MIN_BRACKETED,
                        intrinsic=g1.get("median", np.nan), ci_lo=test.get("ci_lo"),
                        ci_hi=test.get("ci_hi")))
    return out


def _p(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "---"
    return "$<10^{-4}$" if x < 1e-4 else f"${x:.2g}$"


def _fmt_target(kind, t):
    return (f"RMSE ${t:g}$" if kind == "rmse" else f"budget ${100*t:g}\\%$")


def latex_table(df, floors):
    keep = df[df.reported].copy()
    targets = sorted(set(zip(keep.kind, keep.target)),
                     key=lambda kt: (kt[0] != "rmse", -kt[1] if kt[0] == "rmse" else kt[1]))
    note = ("Median per-region ratio of the HEALPix side to the operational "
            "UTM side ($\\gamma_{\\mathrm{tile}}=1.30$) at each target the sweep "
            f"brackets for at least {100*MIN_BRACKETED:.0f}\\% of the regions "
            "under that reader, with the Holm-adjusted paired signed-rank $p$; "
            "the last column is the same ratio at $\\gamma=1$ under the "
            "interpolating read, which for an RMSE target is the critical "
            "overhead $\\gamma_\\star$. \\emph{nearest} is the read the paper's "
            "sweep used; \\emph{bilinear} interpolates within the carrier "
            "lattice (Section~\\ref{sec:operator}). The UTM side is identical "
            "under both, since it depends on the intersection mask alone. The "
            "floor row is the uncompressed carrier's RMSE against the "
            "observation.")
    lines = ["\\begin{table}[t]", "\\centering",
             # The two median columns are at the illustrative overhead and the
             # last column at gamma=1; a reader who meets the table before the
             # prose has no other way to know which.
             *cfg.tex_caption("The common-reference verdict under the two "
                              f"readers, at $\\gamma_{{\\mathrm{{tile}}}}"
                              f"={cfg.GAMMA_UTM}$",
                              "tab:common_two_readers"),
             # Eight columns, two of them $p$-values, do not fit a column at
             # the body font; one step below \\footnotesize and two thirds of
             # the default padding do.
             "\\scriptsize", "\\setlength{\\tabcolsep}{4pt}",
             "\\begin{tabular}{@{}lrrlrrlr@{}}", "\\toprule",
             " & \\multicolumn{3}{c}{nearest} & \\multicolumn{3}{c}{bilinear} & \\\\",
             "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
             "Target & $n$ & median & $p_{\\mathrm{Holm}}$ & $n$ & median & "
             "$p_{\\mathrm{Holm}}$ & $\\gamma_\\star$ \\\\",
             " & & \\multicolumn{2}{c}{} & & \\multicolumn{2}{c}{} & "
             "($\\gamma{=}1$) \\\\", "\\midrule"]
    for kind, t in targets:
        cells = [_fmt_target(kind, t)]
        for reader in ("nearest", "bilinear"):
            g = df[(df.reader == reader) & (df.kind == kind) & (df.target == t)]
            if len(g) and bool(g.reported.iloc[0]):
                g = g.iloc[0]
                med = f"{g['median']:.2f}" if g["median"] < 10 else f"{g['median']:.1f}"
                cells += [str(int(g.n)), med, _p(g.p_holm)]
            else:
                cells += ["---", "censored", "---"]
        gb = df[(df.reader == "bilinear") & (df.kind == kind) & (df.target == t)]
        if len(gb) and bool(gb.reported.iloc[0]):
            v = float(gb.intrinsic.iloc[0])
            gs = f"{v:.2f}" if v < 10 else f"{v:.1f}"
        else:
            gs = "---"
        cells.append(gs)
        lines.append(" & ".join(cells) + " \\\\")
    lines += ["\\midrule",
              f"floor (median) & \\multicolumn{{3}}{{l}}{{{floors['nearest']:.1e}}} & "
              f"\\multicolumn{{3}}{{l}}{{{floors['bilinear']:.1e}}} & \\\\",
              "\\bottomrule", "\\end{tabular}", *cfg.tex_note(note),
              "\\end{table}"]
    return "\n".join(lines) + "\n"


def summarize(df, floors):
    out = ["# The common-reference verdict under two readers", "",
           f"- floors: nearest {floors['nearest']:.3e}, bilinear {floors['bilinear']:.3e}",
           f"- a target is reported when bracketed on both grids for >= "
           f"{100*MIN_BRACKETED:.0f}% of regions", "",
           "| target | nearest n | nearest median | p_Holm | bilinear n | bilinear median | p_Holm | gamma* (bilinear) |",
           "|---|---|---|---|---|---|---|---|"]
    piv = {(r.reader, r.kind, r.target): r for r in df.itertuples()}
    for kind, t in sorted(set(zip(df.kind, df.target))):
        a, b = piv.get(("nearest", kind, t)), piv.get(("bilinear", kind, t))
        def cell(r):
            if r is None:
                return "--- | --- | ---"
            tag = "" if r.reported else " (censored)"
            return f"{r.n} | {r.median:.3f}{tag} | {r.p_holm if r.p_holm is not None else float('nan'):.2g}"
        gs = f"{b.intrinsic:.2f}" if b is not None and b.reported else "---"
        out.append(f"| {kind} {t:g} | {cell(a)} | {cell(b)} | {gs} |")
    out.append("")
    b2 = piv.get(("bilinear", "budget", 0.02)); a2 = piv.get(("nearest", "budget", 0.02))
    r1 = piv.get(("bilinear", "rmse", 0.01))
    if a2 and b2:
        out.append(f"- at the 2% budget the interpolating read moves the ratio from "
                   f"{a2.median:.2f} to {b2.median:.2f}; native UTM still leads "
                   f"(p_Holm {b2.p_holm:.1e})")
    if r1 and r1.reported:
        out.append(f"- at RMSE 0.01 the two grids reach parity under the "
                   f"interpolating read: {r1.median:.3f}, p_Holm {r1.p_holm:.2f}, "
                   f"n={r1.n}; intrinsic ratio {r1.intrinsic:.2f}, so gamma* = "
                   f"{r1.intrinsic:.2f} sits at the illustrative 1.30")
    return "\n".join(out)


def main() -> int:
    docs = {r: load(f) for r, f in READERS}
    df = pd.DataFrame([x for r, d in docs.items() for x in rows(d, r)])
    floors = {}
    for r, tag in (("nearest", ""), ("bilinear", "__bilinear")):
        p = cfg.RESULTS_DIR / f"common_reference_floor{tag}.csv"
        floors[r] = float(pd.read_csv(p).representation_floor_rmse.median()) if p.exists() else np.nan
    df.to_csv(cfg.RESULTS_DIR / "reader_rescore.csv", index=False)
    md = summarize(df, floors)
    (cfg.RESULTS_DIR / "reader_rescore.md").write_text(md, encoding="utf-8")
    cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
    out = cfg.TEX_TABLES / "table_common_two_readers.tex"
    out.write_text(latex_table(df, floors), encoding="utf-8")
    print(md)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
