# -*- coding: utf-8 -*-
"""Adds the paper's two transfer-cost tables to superres_cost.py.

Every column of the main table is computed on the SAME region set at a given
target -- the regions bracketed on the baseline sweep, on every tested
transfer's sweep, and on the native sweep -- so the columns are comparable to
each other and the prose can quote them without a caveat about which subset
each came from.
"""
import sys
from pathlib import Path

p = Path(sys.argv[1] if len(sys.argv) > 1 else "superres_cost.py")
s = p.read_text(encoding="utf-8")
anchor = "def main(argv=None) -> int:"
assert s.count(anchor) == 1

BLOCK = r'''SHORT = {"clough_tocher": "C--T", "psf": "PSF", "linear": "bilin.",
         "nearest": "near.", "bicubic": "bicub.", "utm": "UTM",
         "native UTM": "native UTM"}


def _short(k):
    """Column labels have to fit: a repeated \\texttt{clough\\_tocher} in a
    header overruns even a full-width float, so the header carries short names
    and the caption spells them out."""
    return SHORT.get(k, k.replace("_", chr(92) + "_"))


def lossless_size(df):
    """Smallest lossless size per region: that raster's entropy, model-free."""
    d = df[df.codec.isin(cfg.LOSSLESS_CODECS)]
    return d.groupby("scene").total_bytes.min()


def compactability(sources, codecs=None):
    """Compressed size as a fraction of that raster's own lossless size.

    This is what separates the transfers. Two rasters of nearly equal entropy
    can sit far apart in bytes at a loose target if one concentrates its energy
    in coarse coefficients and the other spreads it into fine ones, because a
    lossy coder discards small coefficients for free and must keep large ones.
    Dividing by each raster's own lossless size removes the entropy difference
    and leaves the compaction.
    """
    codecs = codecs or cfg.COMPRESSOR_CODECS
    rows = []
    for label, df in sources.items():
        ll = lossless_size(df)
        row = dict(transfer=label, lossless_bytes=float(ll.median()),
                   n_lossless=int(len(ll)))
        for t in cfg.TARGET_RMSE:
            v = []
            for sc in ll.index:
                if not ps.any_bracketed(df, sc, codecs, t, "rmse"):
                    continue
                sz = ps.best_of_codecs_size(df, sc, codecs, t)
                if np.isfinite(sz):
                    v.append(sz / ll[sc])
            row[f"t{t:g}"] = float(np.median(v)) if len(v) >= 20 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def _common_regions(sources, target, codecs):
    """Regions bracketed on every sweep at this target."""
    sets = []
    for df in sources.values():
        sets.append({sc for sc in df.scene.unique()
                     if ps.any_bracketed(df, sc, codecs, target, "rmse")
                     and np.isfinite(ps.best_of_codecs_size(df, sc, codecs,
                                                            target))})
    out = sets[0]
    for x in sets[1:]:
        out = out & x
    return sorted(out)


def _median_ratio(num_df, den_df, scenes, target, codecs, scale=1.0):
    """Paired median of num/den over the given regions, and its p-value."""
    r = []
    for sc in scenes:
        a = ps.best_of_codecs_size(den_df, sc, codecs, target)
        b = ps.best_of_codecs_size(num_df, sc, codecs, target)
        if np.isfinite(a) and np.isfinite(b) and a > 0:
            r.append(b / (a * scale))
    if not r:
        return np.nan, np.nan
    st = ps.paired_ratio_test(r)
    return st["median"], st["p"]


def latex_tables(baseline="clough_tocher", tested=("psf", "linear"),
                 codecs=None, publish=True):
    """The two tables Section V-I needs, plus a digest of what the prose quotes."""
    codecs = codecs or cfg.COMPRESSOR_CODECS
    have = {}
    for k in (baseline,) + tuple(tested):
        try:
            have[k] = sweep(k)
        except FileNotFoundError as exc:
            print(f"note: {exc}")
    if baseline not in have or len(have) < 2:
        print("not enough sweeps for the transfer-cost tables")
        return None
    kn = [k for k in tested if k in have]
    utm = utm_sweep()
    sources = {baseline: have[baseline]}
    sources.update({k: have[k] for k in kn})
    if utm is not None:
        sources["utm"] = utm

    digest, rows = {}, []
    for t in cfg.TARGET_RMSE:
        scenes = _common_regions(sources, t, codecs)
        cells, d = [], dict(n=len(scenes))
        for k in kn:
            m, pv = _median_ratio(have[k], have[baseline], scenes, t, codecs)
            cells.append(f"{m:.2f}")
            d[f"{k}_over_base"] = m
            d[f"{k}_over_base_p"] = pv
        if len(kn) == 2:
            m, pv = _median_ratio(have[kn[0]], have[kn[1]], scenes, t, codecs)
            cells.append(f"{m:.2f}")
            d["deconv_over_control"] = m
            d["deconv_over_control_p"] = pv
        if utm is not None:
            for k in [baseline] + kn:
                m, _ = _median_ratio(have[k], utm, scenes, t, codecs,
                                     scale=cfg.GAMMA_UTM)
                cells.append(f"{m:.2f}")
                d[f"{k}_over_utm"] = m
        digest[f"{t:g}"] = d
        rows.append(f"${t:g}$ & {len(scenes)} & " + " & ".join(cells) + r" \\")

    over = " & ".join(_short(k) for k in kn)
    if len(kn) == 2:
        over += f" & {_short(kn[0])}/{_short(kn[1])}"
    n_over = len(kn) + (1 if len(kn) == 2 else 0)
    utm_head = " & ".join(_short(k) for k in [baseline] + kn) \
        if utm is not None else ""
    n_utm = (1 + len(kn)) if utm is not None else 0
    ncol = 2 + n_over + n_utm
    t1 = (
        "\\begin{table*}[t]\n\\centering\n"
        "\\caption{What the transfer costs to store. C--T is "
        "\\texttt{clough\\_tocher}, PSF the PSF-aware transfer, bilin.\\ "
        "bilinear interpolation. Every ratio at a given "
        "target is a paired median over the same regions --- those bracketed "
        "on all sweeps at that target --- so the columns are comparable to one "
        "another. Left block: compressed size relative to the "
        f"\\texttt{{{baseline.replace('_', chr(92)+'_')}}} tile of the same "
        "region at the same absolute RMSE; its last column is the deconvolving "
        "transfer over the non-deconvolving one, which is the cost of "
        "deconvolving separated from the cost of leaving Clough--Tocher. Right "
        "block: each product against the equal-area native crop scaled by "
        f"$\\gamma_{{\\mathrm{{tile}}}}={cfg.GAMMA_UTM}$, composed per region "
        "rather than as a product of medians; below one favours HEALPix.}\n"
        "\\label{tab:transfer_cost}\n\\footnotesize\n"
        f"\\begin{{tabular}}{{@{{}}lr{'r'*(ncol-2)}@{{}}}}\n\\toprule\n"
        f" & & \\multicolumn{{{n_over}}}{{c}}{{Relative to {_short(baseline)}}}"
        + (f" & \\multicolumn{{{n_utm}}}{{c}}{{HPX/UTM}}" if utm is not None
           else "") + " \\\\\n"
        f"\\cmidrule(lr){{3-{2+n_over}}}"
        + (f"\\cmidrule(lr){{{3+n_over}-{ncol}}}" if utm is not None else "")
        + "\n"
        f"Target RMSE & $n$ & {over}"
        + (f" & {utm_head}" if utm is not None else "") + " \\\\\n\\midrule\n"
        + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n\\end{table*}\n")

    csrc = {k: have[k] for k in [baseline] + kn}
    if utm is not None:
        csrc["native UTM"] = utm
    comp = compactability(csrc, codecs)
    body = []
    for _, r in comp.iterrows():
        name = _short(r.transfer)
        cells = [f"{r.lossless_bytes/1e6:.2f}"]
        for t in cfg.TARGET_RMSE:
            v = r[f"t{t:g}"]
            cells.append("--" if v != v else f"{v:.3f}")
        body.append(f"{name} & " + " & ".join(cells) + r" \\")
    t2 = (
        "\\begin{table}[t]\n\\centering\n"
        "\\caption{Why the transfers differ: not in how much they carry, but "
        "in how compactly. Transfers are labelled as in "
        "Table~\\ref{tab:transfer_cost}. The second column is each raster's own lossless "
        "size, which separates the HEALPix rasters by at most a few percent. "
        "The rest is the compressed size at each target as a fraction of that "
        "same lossless size, so the entropy difference divides out and what "
        "remains is how much of the content a lossy coder can discard. A "
        "curvature-minimising $C^1$ interpolant concentrates its energy in "
        "coarse coefficients; a $C^0$ transfer leaves gradient discontinuities "
        "at cell boundaries and a deconvolving one restores fine structure and "
        "amplifies noise with it, and both put energy where a coder targeting "
        "an RMSE has to keep it.}\n"
        "\\label{tab:compactability}\n\\footnotesize\n"
        f"\\begin{{tabular}}{{@{{}}l r {'r'*len(cfg.TARGET_RMSE)}@{{}}}}\n"
        "\\toprule\n"
        f" & Lossless & \\multicolumn{{{len(cfg.TARGET_RMSE)}}}{{c}}"
        "{Fraction of that lossless size} \\\\\n"
        f"\\cmidrule(lr){{3-{2+len(cfg.TARGET_RMSE)}}}\n"
        "Transfer & (MB) & "
        + " & ".join(f"${t:g}$" for t in cfg.TARGET_RMSE) + " \\\\\n"
        "\\midrule\n" + "\n".join(body) + "\n\\bottomrule\n\\end{tabular}\n"
        "\\end{table}\n")

    if publish:
        cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
        (cfg.TEX_TABLES / "table_transfer_cost.tex").write_text(
            t1, encoding="utf-8")
        (cfg.TEX_TABLES / "table_compactability.tex").write_text(
            t2, encoding="utf-8")
        print(f"wrote {(cfg.TEX_TABLES/'table_transfer_cost.tex').name}, "
              f"{(cfg.TEX_TABLES/'table_compactability.tex').name}")
    comp.to_csv(cfg.RESULTS_DIR / "compactability.csv", index=False)
    # Everything the prose quotes, so a number in the text can be checked
    # against the run that produced the table rather than against a memory.
    (cfg.RESULTS_DIR / "transfer_cost_digest.json").write_text(
        json.dumps(digest, indent=2, default=float), encoding="utf-8")
    print(f"wrote {(cfg.RESULTS_DIR/'transfer_cost_digest.json').name}")
    return t1, t2


'''
s = s.replace(anchor, BLOCK + anchor)

old = """import argparse
import shutil"""
new = """import argparse
import json
import shutil"""
assert s.count(old) == 1
s = s.replace(old, new)

old = '''    ap.add_argument("--no-publish", action="store_true")
    args = ap.parse_args(argv)'''
new = '''    ap.add_argument("--no-publish", action="store_true")
    ap.add_argument("--tex", action="store_true",
                    help="also write the paper's transfer-cost tables and the "
                         "digest the prose is checked against; needs the "
                         "baseline and every tested transfer swept")
    args = ap.parse_args(argv)'''
assert s.count(old) == 1
s = s.replace(old, new)

old = '''    if args.control:
        print("\\nRead the two together: a surcharge that the control also "
              "shows is a cost of changing transfer, not of deconvolving.")
    return 0'''
new = '''    if args.control:
        print("\\nRead the two together: a surcharge that the control also "
              "shows is a cost of changing transfer, not of deconvolving.")
    if args.tex:
        latex_tables(args.baseline, tuple(todo), args.codecs,
                     publish=not args.no_publish)
    return 0'''
assert s.count(old) == 1
s = s.replace(old, new)

p.write_text(s, encoding="utf-8")
print("superres_cost.py patched")
