"""Three merged tables for the body of the paper (revision 5, phase 5).

The reviews asked for fewer, denser tables. Three groups of three said the
same thing from adjacent angles and are merged here, each from the digests
its source scripts already write, so that a re-run of any source script and
of this one keeps the merged table in step:

    table_synthesis.tex   the verdict per RMSE target in one place:
                          overhead-free ratio (= gamma_star), IQR, the
                          recovered fraction f(e), the win fractions at three
                          overheads, the inconclusive band, and the ratio,
                          corrected p-value and verdict at the illustrative
                          overhead. Replaces the headline, gamma-continuous
                          and constant-RMSE tables.
                          Sources: gamma_analysis.analyse("self"),
                          results/paper_numbers_self.json.

    table_floors.tex      every transfer-and-read floor against the native
                          observation: by class for the main chain, by reader
                          for both transfers, and the two level-19 carriers.
                          Replaces the representation-floor, reader-floor and
                          level-19-direct tables.
                          Sources: results/common_reference_floor.csv,
                          results/reader_floor.csv, results/level19_direct.json.

    table_transfer.tex    what the transfer costs: the three carriers against
                          the native crop at gamma = 1, their cost relative to
                          Clough--Tocher, the lossless size and compressed
                          fraction by transfer, and the accuracy and spectral
                          columns of the kernel ablation.
                          Replaces the transfer-cost, compactability and
                          kernel-ablation tables.
                          Sources: results/transfer_cost_digest.json,
                          results/compactability.csv, results/kernel_ablation.csv.

Usage
-----
    python merged_tables.py            # writes the three tables
    python merged_tables.py --only synthesis
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

import paper_config as cfg

RMSE_LABEL = {1e-4: r"$10^{-4}$", 3e-4: r"$3\times10^{-4}$", 1e-3: r"$10^{-3}$",
              3e-3: r"$3\times10^{-3}$", 1e-2: r"$10^{-2}$"}


def _label(t):
    return RMSE_LABEL.get(float(t), f"${float(t):g}$")


def _p(v):
    if v is None or not np.isfinite(v):
        return "--"
    if v < 1e-4:
        return r"$<10^{-4}$"
    return f"{v:.3f}"


def _caption(title, label):
    return "\n".join(cfg.tex_caption(title, label))


# --- 1. synthesis ------------------------------------------------------------
def table_synthesis():
    import gamma_analysis as ga
    df, gamma_utm, _ = ga.analyse("self")
    df = df[df.reported].sort_values("target")
    doc = json.loads((cfg.RESULTS_DIR / "paper_numbers_self.json").read_text())
    at_g = {}
    for t in doc["targets"]:
        # the JSON keeps every target; gamma_analysis says which are reported
        if t["kind"] == "rmse" and t.get("test", {}).get("n"):
            at_g[round(float(t["target"]), 12)] = t
    rows = []
    for _, r in df.iterrows():
        t = at_g.get(round(float(r.target), 12))
        if t is None:
            continue
        test, ph = t["test"], t.get("p_holm", np.nan)
        if ph is not None and np.isfinite(ph) and ph >= 0.05:
            verdict = "inconcl."
        else:
            verdict = "HPX" if test["median"] < 1 else "UTM"
        inc = ("--" if not np.isfinite(r.inconclusive_lo)
               else f"${r.inconclusive_lo:.2f}$--${r.inconclusive_hi:.2f}$")
        rows.append(
            f"{_label(r.target)} & {int(r.n)} & {r.gamma_star_median:.2f} & "
            f"[{r.gamma_star_q25:.2f}, {r.gamma_star_q75:.2f}] & "
            f"{r.recovered_fraction:.2f} & "
            f"{100*r['win_at_1']:.0f} & {100*r['win_at_1.3']:.0f} & "
            f"{100*r['win_at_2']:.0f} & {inc} & "
            f"{test['median']:.2f} & {_p(ph)} & {verdict} \\\\")
    if not rows:
        return ""
    return (
        "\\begin{table*}[t]\n\\centering\n"
        + _caption("The verdict per RMSE target: overhead-free ratio, critical "
                   "overhead and the illustrative operating point",
                   "tab:synthesis")
        + "\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n"
        "\\begin{tabular}{@{}lrrlrrrrlrrl@{}}\n\\toprule\n"
        " & & \\multicolumn{3}{c}{$\\gamma=1$} & "
        "\\multicolumn{3}{c}{HPX smaller (\\% of $n$)} & & "
        f"\\multicolumn{{3}}{{c}}{{$\\gamma_{{\\mathrm{{tile}}}}={gamma_utm:g}$}} \\\\\n"
        "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){10-12}\n"
        "Target RMSE & $n$ & $\\gamma_\\star$ & IQR & $f(e)$ & "
        "$\\gamma{=}1$ & $\\gamma{=}1.3$ & $\\gamma{=}2$ & Inconclusive $\\gamma$ & "
        "Ratio & $p_{\\mathrm{Holm}}$ & Verdict \\\\\n\\midrule\n"
        + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n\\end{table*}\n")


# --- 2. floors ---------------------------------------------------------------
def _sci(x):
    return f"{x:.2e}".replace("e-0", "e-")


def table_floors():
    rows = []
    fl = pd.read_csv(cfg.RESULTS_DIR / "common_reference_floor.csv")
    fl = fl[fl["mode"] == "nearest"] if "mode" in fl else fl
    fl = fl.dropna(subset=["representation_floor_rmse"])
    std_col = "utm_std_on_mask" if "utm_std_on_mask" in fl else "utm_std"
    for cls, g in fl.groupby("scene_class"):
        pct = 100 * np.median(g.representation_floor_rmse / g[std_col])
        rows.append(f"level 20 & C--T & nearest & {cls} & {len(g)} & "
                    f"{_sci(g.representation_floor_rmse.median())} & {pct:.1f} & "
                    f"[{_sci(g.representation_floor_rmse.min())}, "
                    f"{_sci(g.representation_floor_rmse.max())}] \\\\")
    pct = 100 * np.median(fl.representation_floor_rmse / fl[std_col])
    rows.append(f"level 20 & C--T & nearest & all & {len(fl)} & "
                f"{_sci(fl.representation_floor_rmse.median())} & {pct:.1f} & "
                f"[{_sci(fl.representation_floor_rmse.min())}, "
                f"{_sci(fl.representation_floor_rmse.max())}] \\\\")
    rows.append("\\midrule")
    rf = pd.read_csv(cfg.RESULTS_DIR / "reader_floor.csv")
    rf = rf[rf.error.isna()] if "error" in rf else rf
    names = {"clough_tocher": "C--T", "psf": "PSF"}
    order = [("clough_tocher", "nearest"), ("clough_tocher", "bilinear"),
             ("psf", "nearest"), ("psf", "bilinear"), ("psf", "psf")]
    for k, rd in order:
        g = rf[(rf.kernel == k) & (rf.reader == rd)]
        if not len(g):
            continue
        if k == "clough_tocher" and rd == "nearest":
            continue            # the class block above already holds this row
        med = float(g.floor.median())
        note = f"{100*g.support_fraction.median():.1f}\\% of cells read"
        if rd != "nearest":
            # paired per-region ratio to the nearest read of the same tiles,
            # as reader_floor.py reports it (median ratio, then 1 - r^2)
            near = rf[(rf.kernel == k) & (rf.reader == "nearest")].set_index("scene").floor
            ratio = (g.set_index("scene").floor / near).dropna()
            if len(ratio):
                note += f"; MSE $-{100*(1 - float(ratio.median()) ** 2):.0f}\\%$ vs nearest"
        rows.append(f"level 20 & {names[k]} & {rd} & all & {len(g)} & {_sci(med)} & "
                    f"{100*g.floor_frac_of_std.median():.1f} & {note} \\\\")
    rows.append("\\midrule")
    d19 = json.loads((cfg.RESULTS_DIR / "level19_direct.json").read_text())
    rows.append(f"level 19, direct & C--T & nearest & all & {d19['n']} & "
                f"{_sci(d19['floor19_direct'])} & -- & "
                f"${d19['ratio_direct']:.2f}\\times$ level 20; lower on {d19['direct_better_on']} \\\\")
    rows.append(f"level 19, aggregate & C--T & nearest & all & {d19['n']} & "
                f"{_sci(d19['floor19_agg'])} & -- & "
                f"${d19['ratio_agg']:.2f}\\times$ level 20; lower on "
                f"{d19['n'] - d19['direct_better_on']} \\\\")
    return (
        "\\begin{table*}[t]\n\\centering\n"
        + _caption("Transfer-and-read floors against the native observation",
                   "tab:floors")
        + "\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n"
        "\\begin{tabular}{@{}llllrrrl@{}}\n\\toprule\n"
        "Carrier & Transfer & Read & Class & $n$ & Floor (RMSE) & \\% of std & "
        "Range or note \\\\\n\\midrule\n"
        + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n\\end{table*}\n")


# --- 3. transfer -------------------------------------------------------------
def table_transfer():
    tc = json.loads((cfg.RESULTS_DIR / "transfer_cost_digest.json").read_text())
    # One population for the whole block: the regions bracketed on all four
    # sweeps (superres_cost.py), so the three carriers are compared on the
    # same regions. The digest stores X/UTM at gamma_tile; gamma = 1 is that
    # times gamma_tile.
    rows = []
    for t in cfg.TARGET_RMSE:
        b = tc.get(f"{t:g}")
        if not b:
            continue
        g = cfg.GAMMA_UTM
        rows.append(f"{_label(t)} & {b['n']} & {g*b['clough_tocher_over_utm']:.2f} & "
                    f"{g*b['psf_over_utm']:.2f} & {g*b['linear_over_utm']:.2f} & "
                    f"{b['psf_over_base']:.2f} & {b['linear_over_base']:.2f} & "
                    f"{b['deconv_over_control']:.2f} \\\\")
    part_a = (
        "\\begin{tabular}{@{}lrrrrrrr@{}}\n\\toprule\n"
        " & & \\multicolumn{3}{c}{HPX/UTM at $\\gamma=1$} & "
        "\\multicolumn{3}{c}{Relative to C--T} \\\\\n"
        "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\n"
        "Target RMSE & $n$ & C--T & PSF & IDW & PSF & IDW & PSF/IDW \\\\\n"
        "\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}")

    comp = pd.read_csv(cfg.RESULTS_DIR / "compactability.csv")
    names = {"clough_tocher": "C--T", "psf": "PSF", "linear": "IDW",
             "native UTM": "native UTM"}
    crow = []
    for _, r in comp.iterrows():
        vals = " & ".join(f"{r[f't{t:g}']:.3f}" if np.isfinite(r[f"t{t:g}"]) else "--"
                          for t in cfg.TARGET_RMSE)
        crow.append(f"{names.get(r.transfer, r.transfer)} & "
                    f"{r.lossless_bytes/1e6:.2f} & {vals} \\\\")
    heads = " & ".join(_label(t) for t in cfg.TARGET_RMSE)
    part_b = (
        "\\begin{tabular}{@{}lrrrrrr@{}}\n\\toprule\n"
        " & Lossless & \\multicolumn{5}{c}{Fraction of own lossless size at target RMSE} \\\\\n"
        "\\cmidrule(lr){3-7}\n"
        f"Transfer & (MB) & {heads} \\\\\n\\midrule\n"
        + "\n".join(crow) + "\n\\bottomrule\n\\end{tabular}")

    ka = pd.read_csv(cfg.RESULTS_DIR / "kernel_ablation.csv")
    from kernel_ablation import paired          # same pairing as the source table
    kname = {"clough_tocher": "C--T", "psf": "PSF", "linear": "IDW",
             "nearest": "nearest", "bicubic": "bicubic"}
    fp = paired(ka, "floor").set_index("kernel")
    sizes = {}
    for t in cfg.TARGET_RMSE[:3]:
        col = f"size_at_{t:g}"
        if col in ka:
            sizes[t] = paired(ka, col).set_index("kernel")
    krow = []
    for k in ("clough_tocher", "psf", "linear", "nearest", "bicubic"):
        g = ka[ka.kernel == k]
        if not len(g):
            continue
        tf = (f"{100*g.transfer_frac_of_std.median():.2f}"
              if np.isfinite(g.transfer_frac_of_std.median()) else "--")
        ovr = (f"{100*g.cells_outside_source_range.median():.2f}"
               if "cells_outside_source_range" in g else "--")
        rel = f"{fp.loc[k, 'median']:.2f}" if k in fp.index else "--"
        cells = "".join(f" & {sizes[t].loc[k, 'median']:.2f}"
                        if t in sizes and k in sizes[t].index else " & --"
                        for t in cfg.TARGET_RMSE[:3])
        krow.append(f"{kname[k]} & {len(g)} & {tf} & {ovr} & "
                    f"{100*g.floor_frac_of_std.median():.1f} & {rel} & "
                    f"{g.logpsd_rms.median():.2f} & {g.logpsd_top_octave.median():.2f}"
                    + cells + " \\\\")
    heads3 = " & ".join(_label(t) for t in cfg.TARGET_RMSE[:3])
    part_c = (
        "\\begin{tabular}{@{}lrrrrrrrrrr@{}}\n\\toprule\n"
        " & & Transfer & Over & Floor & Floor & & Top & "
        "\\multicolumn{3}{c}{Size / C--T at target RMSE} \\\\\n"
        "\\cmidrule(lr){9-11}\n"
        f"Kernel & $n$ & (\\% std) & (\\%) & (\\% std) & rel.\\ C--T & log-PSD & oct. & {heads3} \\\\\n"
        "\\midrule\n" + "\n".join(krow) + "\n\\bottomrule\n\\end{tabular}")

    return (
        "\\begin{table*}[t]\n\\centering\n"
        + _caption("What the transfer costs: storage, compactability and accuracy",
                   "tab:transfer")
        + "\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n"
        + part_a + "\n\\vspace{4pt}\n\n" + part_b + "\n\\vspace{4pt}\n\n" + part_c
        + "\n\\end{table*}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", choices=["synthesis", "floors", "transfer"], default=None)
    args = ap.parse_args(argv)
    todo = {"synthesis": table_synthesis, "floors": table_floors,
            "transfer": table_transfer}
    if args.only:
        todo = {args.only: todo[args.only]}
    cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
    for name, fn in todo.items():
        tex = fn()
        if not tex:
            print(f"  {name}: nothing to write")
            continue
        out = cfg.TEX_TABLES / f"table_{name}.tex"
        out.write_text(tex, encoding="utf-8")
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
