"""The one number a reader wants: how big is the HEALPix tile, relative to UTM?

Everything in this study reduces, for a practitioner, to a single question:
once compressed to the same reconstruction error, how many bytes does the
level-20 HEALPix tile take compared with the native UTM raster of the same
ground? The answer is spread over several sections of the paper because it
depends on three things -- the tolerated error, the transfer that populated
the carrier, and the overhead charged to the native product -- and each is
established separately there. This script assembles the answer in one table
and one figure, with no assumption about overhead (gamma = 1), so that the
abstract and the introduction can quote it and a reader can find it without
reading the rest.

Sources (nothing is measured here):
  paper_numbers_self.json     per-region ratios at gamma_tile, hence at gamma=1
  transfer_cost_digest.json   the PSF-aware and inverse-distance carriers, on
                              the regions all three sweeps bracket
  level_choice.csv            the tile's standard deviation, to express each
                              RMSE target as a fraction of scene contrast
  storage_accounting          the raw, uncompressed ratio (a geometric constant)

Usage
-----
    python headline.py
"""

from __future__ import annotations

import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402

import paper_config as cfg                # noqa: E402

RMSE_LABEL = {1e-4: r"$10^{-4}$", 3e-4: r"$3\times10^{-4}$", 1e-3: r"$10^{-3}$",
              3e-3: r"$3\times10^{-3}$", 1e-2: r"$10^{-2}$"}
CLASSES = ("urban", "agriculture", "forest", "water")


def load():
    doc = json.loads((cfg.RESULTS_DIR / "paper_numbers_self.json").read_text())
    gamma = float(doc["gamma_utm"])
    per_target = {}
    for t in doc["targets"]:
        if t["kind"] != "rmse" or not t.get("reported", True):
            continue
        r = pd.Series(t["ratios"], dtype=float) * gamma       # gamma = 1
        cls = r.index.map(lambda s: s.split("__")[1])
        per_target[float(t["target"])] = dict(
            median=float(r.median()), n=int(r.notna().sum()),
            by_class={c: float(r[cls == c].median()) for c in CLASSES})
    raw = float(doc["storage_accounting"]["ratio_gamma1"]["median"])
    p = cfg.RESULTS_DIR / "transfer_cost_digest.json"
    digest = json.loads(p.read_text()) if p.exists() else {}
    transfers = {}
    for k, v in digest.items():
        # The digest also carries diagnostics that are not per-target blocks
        # (superres_cost writes the lossless-size correlation there); only the
        # numeric keys are RMSE targets.
        try:
            target = float(k)
        except ValueError:
            continue
        if not isinstance(v, dict) or "n" not in v:
            continue
        transfers[target] = dict(
            n=int(v["n"]),
            ct=float(v["clough_tocher_over_utm"]) * gamma,
            psf=float(v["psf_over_utm"]) * gamma,
            idw=float(v["linear_over_utm"]) * gamma)
    lc = pd.read_csv(cfg.RESULTS_DIR / "level_choice.csv")
    std = float(lc[lc.level == cfg.LEVEL].tile_std.median())
    return per_target, raw, transfers, std, gamma


def latex_table(per_target, raw, transfers, std, gamma):
    targets = sorted(per_target)
    rows = []
    rows.append("uncompressed & -- & "
                + " & ".join([f"{raw:.2f}"] * 3) + r" & -- \\")
    rows.append(r"\midrule")
    for t in targets:
        d = per_target[t]
        bc = d["by_class"]
        lo, hi = min(bc.values()), max(bc.values())
        tr = transfers.get(t)
        psf = f"{tr['psf']:.2f}" if tr else "--"
        idw = f"{tr['idw']:.2f}" if tr else "--"
        n_tr = f"{tr['n']}" if tr else "--"
        rows.append(f"{RMSE_LABEL.get(t, f'${t:g}$')} & {100*t/std:.2g} & "
                    f"{d['median']:.2f} [{lo:.2f}--{hi:.2f}] & {psf} & {idw} & "
                    f"{n_tr} \\\\")
    note = (
        "Median over the 40 regions of the compressed size of the HEALPix "
        "tile divided by that of the native UTM crop of the same nominal "
        "ground area, at the same reconstruction error and with no storage "
        "overhead charged to either side ($\\gamma=1$); a value above one "
        "means the HEALPix tile is larger. The second column expresses each "
        "target as a fraction of a typical scene's contrast (the median "
        f"standard deviation of the 40 tiles, {std:.3f} in reflectance), so "
        "that $10^{-2}$ is an error of about a quarter of the contrast and "
        "$10^{-4}$ a quarter of a percent. The C--T (Clough--Tocher) column is "
        "the carrier the paper uses, with the range of the four land-cover "
        "class medians in brackets; the PSF (PSF-aware) and IDW "
        "(inverse-distance) columns are the same tile rebuilt with the other "
        "two transfers of "
        "Section~\\ref{sec:transfer_cost}, on the $n$ regions all three "
        "sweeps bracket. Under an operational overhead $\\gamma$ on the native "
        f"product every entry divides by $\\gamma$ ({gamma:g} in the "
        "illustrative scenario of Section~\\ref{sec:rd_illustrative}).")
    return "\n".join([
        r"\begin{table}[t]", r"\centering",
        *cfg.tex_caption("How much larger is the HEALPix tile, once compressed?",
                         "tab:headline"),
        r"\scriptsize", r"\setlength{\tabcolsep}{4pt}",
        r"\begin{tabular}{@{}lrlrrr@{}}", r"\toprule",
        r"Target & Error & \multicolumn{3}{c}{HEALPix / UTM size, by transfer} & \\",
        r"\cmidrule(lr){3-5}",
        r"RMSE & (\% contrast) & C--T [class range] & PSF & IDW & $n$ \\",
        r"\midrule", *rows, r"\bottomrule", r"\end{tabular}",
        *cfg.tex_note(note), r"\end{table}"]) + "\n"


def figure(per_target, raw, transfers, std, gamma, publish=True):
    targets = sorted(per_target)
    x = np.array(targets)
    med = np.array([per_target[t]["median"] for t in targets])
    lo = np.array([min(per_target[t]["by_class"].values()) for t in targets])
    hi = np.array([max(per_target[t]["by_class"].values()) for t in targets])

    fig, ax = plt.subplots(figsize=(5.0, 3.3))
    ax.axhline(raw, color="0.35", linestyle="--", linewidth=1)
    ax.text(x[0], raw * 1.04, f"uncompressed: {raw:.2f}", fontsize=8, color="0.35",
            va="bottom")
    ax.axhline(1.0, color="0.6", linestyle=":", linewidth=1)
    ax.text(x[2], 1.02, "same size as UTM", fontsize=8, color="0.45", ha="center",
            va="bottom")
    ax.axhline(1.0 / gamma, color="0.6", linestyle=":", linewidth=1)
    ax.text(x[-1], 1.0 / gamma * 1.02, f"same size once UTM carries a "
            f"{100*(gamma-1):.0f}% overhead", fontsize=7.5, color="0.45",
            ha="right", va="bottom")

    ax.fill_between(x, lo, hi, color="#0072B2", alpha=0.18, linewidth=0,
                    label="Clough–Tocher, range of class medians")
    ax.plot(x, med, "o-", color="#0072B2", label="Clough–Tocher (the paper's carrier)")
    if transfers:
        tt = sorted(t for t in transfers if t in per_target)
        ax.plot(tt, [transfers[t]["psf"] for t in tt], "s-", color="#D55E00",
                label="PSF-aware transfer")
        ax.plot(tt, [transfers[t]["idw"] for t in tt], "^-", color="#009E73",
                label="inverse-distance transfer")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("tolerated reconstruction error\n"
                  "(RMSE in reflectance; below, as a share of scene contrast)",
                  fontsize=8)
    ax.set_ylabel("HEALPix tile size / native UTM size")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{t:g}\n{100*t/std:.2g}%" for t in x], fontsize=7.5)
    ax.tick_params(axis="x", which="minor", bottom=False, labelbottom=False)
    ax.set_ylim(0.42, 3.2)
    ax.set_yticks([0.5, 0.7, 1, 1.5, 2, 2.58, 3])
    ax.set_yticklabels(["0.5", "0.7", "1", "1.5", "2", "2.58", "3"])
    ax.tick_params(axis="y", which="minor", left=False, labelleft=False)
    ax.grid(alpha=0.25, which="both", linewidth=0.4)
    ax.legend(fontsize=7, frameon=False, loc="lower left", borderaxespad=0.3)
    fig.tight_layout()
    cfg.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(cfg.FIGURES_DIR / f"fig_headline.{ext}", dpi=160,
                    bbox_inches="tight")
    if publish:
        import shutil
        cfg.TEX_FIGURES.mkdir(parents=True, exist_ok=True)
        src, dst = cfg.FIGURES_DIR / "fig_headline.pdf", cfg.TEX_FIGURES / "fig_headline.pdf"
        if src.resolve() != dst.resolve():
            shutil.copyfile(src, dst)
    plt.close(fig)


def digest(per_target, raw, transfers, std, gamma):
    """The sentences the abstract and the introduction quote."""
    ts = sorted(per_target)
    tight, loose = per_target[ts[0]], per_target[ts[-1]]
    out = {
        "raw_ratio": raw,
        "scene_std_median": std,
        "tightest": dict(target=ts[0], ratio=tight["median"],
                         pct_of_std=100 * ts[0] / std,
                         class_lo=min(tight["by_class"].values()),
                         class_hi=max(tight["by_class"].values())),
        "loosest": dict(target=ts[-1], ratio=loose["median"],
                        pct_of_std=100 * ts[-1] / std,
                        class_lo=min(loose["by_class"].values()),
                        class_hi=max(loose["by_class"].values())),
        "transfers_gamma1": {f"{t:g}": v for t, v in transfers.items()},
    }
    if transfers:
        vals = [transfers[t] for t in sorted(transfers)]
        out["psf_gamma1_range"] = [min(v["psf"] for v in vals), max(v["psf"] for v in vals)]
        out["idw_gamma1_range"] = [min(v["idw"] for v in vals), max(v["idw"] for v in vals)]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-publish", action="store_true")
    args = ap.parse_args(argv)
    per_target, raw, transfers, std, gamma = load()
    if not transfers:
        print("note: transfer_cost_digest.json missing; the PSF-aware and "
              "inverse-distance columns will be blank (run superres_cost.py --tex)")
    cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
    (cfg.TEX_TABLES / "table_headline.tex").write_text(
        latex_table(per_target, raw, transfers, std, gamma), encoding="utf-8")
    figure(per_target, raw, transfers, std, gamma, publish=not args.no_publish)
    d = digest(per_target, raw, transfers, std, gamma)
    (cfg.RESULTS_DIR / "headline.json").write_text(json.dumps(d, indent=1))
    print(json.dumps(d, indent=1))
    print(f"wrote {cfg.TEX_TABLES / 'table_headline.tex'} and fig_headline")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
