"""How much of the transfer-and-read floor is the read? (revision 3, item 3.2)

Section V-F of the paper measures the nearest-cell HEALPix -> UTM read as a
sub-cell positional jitter of about 2.1 m, estimates that the read alone would
produce an error of 14 % of the field's standard deviation on a field with the
imagery's spectrum, and then says: "What would lower it is ... a reconstructor
that interpolates within the carrier raster instead of reading the containing
cell ... We have not implemented one, and we do not claim to know how much of
the floor it would remove." The reviewer asked for exactly that reconstructor,
because the common-reference verdict of Section V-H rests on this floor.

`common_reference.CommonReference(mode="bilinear")` is that reconstructor. This
script measures, on every region and with no codec in the loop, the floor of
the uncompressed carrier under each reader:

    nearest    the containing cell -- the paper's reader
    bilinear   bilinear within the carrier's own lattice, at sub-cell position
    psf        the companion's forward operator (on the PSF tiles only, where
               it is the matched reader)

on the Clough--Tocher tiles the paper's results use, and on the PSF tiles of
Section V-J. Everything is paired per region, so the headline is a ratio:

    floor_bilinear / floor_nearest  on the same tile, same mask, same observation

and its square is the share of the floor's variance the read was responsible
for. If that share is large, the reversal under the common reference is a
property of the reader and Section V-H must be re-scored with the interpolating
one (run_codec_bench.py --common-reference --cr-mode bilinear, under a fresh
part tag). If it is small, the reader is exonerated and the floor is the
transfer's.

Usage
-----
    python reader_floor.py --limit 4
    python reader_floor.py
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

import paper_config as cfg
import build_tiles
import common_reference
import region_sites

READERS = {"clough_tocher": ("nearest", "bilinear"),
           "psf": ("nearest", "bilinear", "psf")}


def region_rows(scene, kernels, verbose=True):
    rows = []
    for kernel in kernels:
        path = build_tiles.tile_path(scene, kernel)
        if not path.exists():
            print(f"    no {kernel} tile ({path.name}); skipped")
            continue
        cell_ids, cell_data, utm_crop, meta, npz = build_tiles.load_tile(scene, kernel)
        x = np.asarray(cell_data, dtype=np.float64).ravel()
        obs = np.asarray(utm_crop, dtype=np.float64)
        for reader in READERS[kernel]:
            try:
                cr = common_reference.build_for_tile(npz, mode=reader)
            except Exception as exc:            # noqa: BLE001
                print(f"    {kernel}/{reader} FAILED: {type(exc).__name__}: {exc}")
                rows.append(dict(scene=scene, cls=scene.split("__")[1],
                                 kernel=kernel, reader=reader,
                                 error=f"{type(exc).__name__}: {exc}"))
                continue
            floor = cr.rmse_on_mask(cr.predict(x), obs)
            sd = float(obs[cr.mask].std())
            rows.append(dict(
                scene=scene, cls=scene.split("__")[1], kernel=kernel,
                reader=reader, mode_used=cr.mode, floor=floor,
                floor_frac_of_std=floor / sd if sd > 0 else np.nan,
                obs_std_on_mask=sd, mask_fraction=cr.mask_fraction,
                support_fraction=float(cr.read_support().mean()),
                constant_field=cr.support_median,
                offset_q99=getattr(cr, "offset_q99", np.nan), error=""))
            if verbose:
                print(f"    {kernel:14s} {reader:9s} floor {floor:.3e} "
                      f"({100*floor/sd:5.1f}% of std)  support "
                      f"{100*rows[-1]['support_fraction']:5.1f}%  "
                      f"{cr.describe()}")
    return rows


def summarize(df):
    ok = df[df.error.fillna("").astype(str).str.len() == 0]
    out = ["# How much of the floor is the read?", "",
           f"- regions: {ok.scene.nunique()}; no codec involved, the "
           f"uncompressed carrier read on the observation over the "
           f"intersection mask", ""]
    out += ["| tiles | reader | floor (median) | % of std | cells reached | "
            "regions |", "|---|---|---|---|---|---|"]
    for (k, r), g in ok.groupby(["kernel", "reader"], sort=False):
        out.append(f"| {k} | {r} | {g.floor.median():.3e} | "
                   f"{100*g.floor_frac_of_std.median():.1f}% | "
                   f"{100*g.support_fraction.median():.1f}% | {len(g)} |")
    out.append("")
    for k in ok.kernel.unique():
        piv = ok[ok.kernel == k].pivot_table(index="scene", columns="reader",
                                             values="floor", aggfunc="median")
        if "nearest" not in piv or "bilinear" not in piv:
            continue
        r = (piv.bilinear / piv.nearest).dropna()
        share = 1.0 - r ** 2
        out.append(f"- **{k} tiles**: bilinear/nearest floor = "
                   f"{r.median():.3f} (paired median, q10 {r.quantile(.1):.3f}, "
                   f"q90 {r.quantile(.9):.3f}; lower on {(r < 1).sum()}/{len(r)}"
                   f"), so the read accounts for **{100*share.median():.0f}%** "
                   f"of the floor's variance on the median region")
        if "psf" in piv:
            rp = (piv.psf / piv.nearest).dropna()
            rb = (piv.psf / piv.bilinear).dropna()
            out.append(f"  - psf/nearest {rp.median():.3f}, psf/bilinear "
                       f"{rb.median():.3f}: the matched operator against the "
                       f"best non-deconvolving read")
    out.append("")
    ct = ok[(ok.kernel == "clough_tocher")]
    if len(ct):
        piv = ct.pivot_table(index="scene", columns="reader", values="floor")
        if "nearest" in piv and "bilinear" in piv:
            share = float((1.0 - (piv.bilinear / piv.nearest) ** 2).median())
            if share >= 0.5:
                out.append("- **The floor is mostly the read.** Section V-H's "
                           "common-reference verdict was scored through the "
                           "nearest read and must be re-scored with the "
                           "interpolating one before it is interpreted as a "
                           "statement about the carrier.")
            elif share >= 0.2:
                out.append("- The read is a substantial but not dominant part "
                           "of the floor; the common-reference verdict should "
                           "be reported under both readers.")
            else:
                out.append("- **The read is exonerated**: an interpolating "
                           "reconstructor removes little of the floor, which "
                           "is therefore the transfer's, and the "
                           "common-reference verdict stands as scored.")
    return "\n".join(out)


def latex_table(df, path):
    ok = df[df.error.fillna("").astype(str).str.len() == 0]
    n = ok.scene.nunique()
    per = ok.groupby("reader").scene.nunique()
    short = {r: int(v) for r, v in per.items() if v < n}
    cover = f"{n} regions" + (
        " (" + ", ".join(f"\\emph{{{r}}} on {v}" for r, v in sorted(short.items()))
        + ", the others on all)" if short else "")
    note = (f"On {cover}: the RMSE of the \\emph{{uncompressed}} carrier, read "
            "on the native pixels of the footprints' intersection, against the "
            "Sentinel-2 observation. \\emph{nearest} is the containing cell, "
            "the read used throughout the paper; \\emph{bilinear} interpolates "
            "within the carrier's own lattice at sub-cell position; "
            "\\emph{psf} is the companion's forward operator, reported only on "
            "the tiles it is matched to. Ratios are paired per region, then "
            "the median is taken.")
    lines = ["\\begin{table}[t]", "\\centering",
             *cfg.tex_caption("The transfer-and-read floor under three readers",
                              "tab:reader_floor"), "\\footnotesize",
             "\\begin{tabular}{@{}llrrr@{}}", "\\toprule",
             "Tiles & Reader & Floor & \\% of std & Cells read \\\\", "\\midrule"]
    for (k, r), g in ok.groupby(["kernel", "reader"], sort=False):
        kname = {"clough_tocher": "Clough--Tocher", "psf": "PSF"}.get(k, k)
        lines.append(f"{kname} & {r} & {g.floor.median():.2e} & "
                     f"{100*g.floor_frac_of_std.median():.1f} & "
                     f"{100*g.support_fraction.median():.1f}\\% \\\\")
    lines += ["\\midrule"]
    for k in ok.kernel.unique():
        piv = ok[ok.kernel == k].pivot_table(index="scene", columns="reader",
                                             values="floor")
        if "nearest" in piv and "bilinear" in piv:
            r = (piv.bilinear / piv.nearest).dropna()
            kname = {"clough_tocher": "Clough--Tocher", "psf": "PSF"}.get(k, k)
            lines.append(f"\\multicolumn{{2}}{{@{{}}l}}{{{kname}: bilinear / "
                         f"nearest}} & {r.median():.3f} & "
                         # Not a variance decomposition: it is the measured
                         # drop in mean squared error between two complete
                         # chains that differ only in the read.
                         f"\\multicolumn{{2}}{{l@{{}}}}{{MSE reduction, "
                         f"{100*(1 - r.median()**2):.0f}\\%}} \\\\")
    lines += ["\\bottomrule", "\\end{tabular}", *cfg.tex_note(note),
              "\\end{table}"]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  wrote {path}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--kernels", nargs="*", default=list(READERS))
    ap.add_argument("--no-tex", action="store_true")
    args = ap.parse_args(argv)

    scenes = args.regions or list(region_sites.load_sites()["scene"])
    scenes = [s for s in scenes if build_tiles.tile_path(s).exists()]
    if args.limit:
        scenes = scenes[:args.limit]
    if not scenes:
        print(f"no tile under {cfg.TILE_DIR}; run build_tiles.py first")
        return 1

    rows = []
    for s in scenes:
        print(f"  {s}")
        rows += region_rows(s, args.kernels)
    df = pd.DataFrame(rows)
    if "error" not in df:
        df["error"] = ""
    ok = df[df.error.fillna("").astype(str).str.len() == 0]
    if not len(ok):
        errs = sorted(set(df.error.dropna()))
        raise SystemExit("no reader could be built on any tile:\n  "
                         + "\n  ".join(errs))
    df.to_csv(cfg.RESULTS_DIR / "reader_floor.csv", index=False)
    md = summarize(df)
    (cfg.RESULTS_DIR / "reader_floor.md").write_text(md, encoding="utf-8")
    if not args.no_tex:
        latex_table(df, cfg.REPO_ROOT / "tex" / "tables" / "table_reader_floor.tex")
    print("\n" + md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
