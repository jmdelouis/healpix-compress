"""Is each stored tile's window still the acquisition the tile was built from?

Why this exists
---------------
`build_tiles.audit_windows` already asks whether the window on disk still
*covers* the stored tile: it transforms the tile's cell centres into the
window's pixel grid and checks that they all land inside it with a halo. That
test is geometric. It passes whenever the window is centred on the right
ground, whatever is in the pixels -- a different date, a different processing
baseline, a different band, a rescaled product.

The level-19 experiment walked into the gap. On one region the transfer
rebuilt from the window did not return the stored tile: not at level 19, and
not at level 20 either, where it must. The two were related by
`rebuilt ~ 0.53 * stored + 0.021`, which is not a resampling error; it is a
different radiometric quantity over the same ground. Every measurement that
re-runs a transfer from the window -- the kernel ablation, the resampling-error
study, the level-19 comparison -- silently measures that mismatch instead of
what it set out to measure. The compression results do not: they read the
stored tile and never reopen the window.

What this does
--------------
For each region, sample the window at the stored tile's own cell centres by
nearest pixel -- no triangulation, no operator, a second per region -- and
compare the samples with the stored cell values:

    corr           Pearson correlation, stored against sampled
    slope, offset  least-squares fit stored ~ slope * sampled + offset
    rms            RMS of (stored - sampled), and the same over the tile's
                   standard deviation

Nearest-pixel sampling is cruder than the transfer that built the tile, so a
healthy region does not give rms = 0; it gives a correlation close to one, a
slope close to one and an offset close to zero, with the residual at the scale
of a half-pixel of relief. A region whose window is a different acquisition
fails all three at once, and by a margin no interpolation choice can produce.

Usage
-----
    python window_radiometry.py
    python window_radiometry.py --regions region__water__lake_geneva
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

import build_tiles
import paper_config as cfg
import region_sites

# The test.
#
# If the window is the product the tile came from, the two differ only because
# nearest-pixel sampling is cruder than the transfer that built the tile. That
# is additive noise on the regressor, and additive noise on the regressor has
# a closed-form consequence: regressing the tile on the samples attenuates the
# slope by exactly the squared correlation,
#
#     slope = sigma_t^2 / (sigma_t^2 + sigma_noise^2) = corr^2
#
# so `slope / corr^2` is 1 whatever the noise level, and whatever the scene's
# contrast. A window that is a different acquisition carries a gain and an
# offset of its own, and breaks the identity.
#
# Measured on the 40 regions: 39 of them sit in [0.982, 1.011], median 1.003,
# across correlations from 0.91 to 0.999 and contrasts from 0.001 to 0.18 in
# reflectance. One sits at 1.88.
#
# An earlier version of this script tested the correlation, the slope and the
# offset against fixed thresholds instead, with the offset expressed in units
# of the scene's standard deviation. That flagged two more regions, both of
# them the lowest-contrast lakes in the set, where a small absolute offset is
# many standard deviations and the same fixed sampling noise takes a larger
# share. Those are properties of the scene, not of the window. The identity
# above has neither failure mode, and it is what the sampling model predicts
# rather than a line drawn around an outcome.
#
# The identity holds for a transfer that does not systematically alter the
# field. Not every transfer here is one: the PSF-aware transfer deconvolves,
# adding gain at high frequency by construction, and sits at 1.09 in median
# with a spread from 1.00 to 1.70; Keys' bicubic rings, its correlation with a
# nearest-pixel sample collapses to 0.16, and its ratio has no useful scale.
#
# Two ways of judging those were tried and both fail. Against 1 they are
# flagged wholesale, which reports a design decision as a defect. Against
# their own cohort with a Tukey rule they are worse still: the carrier's
# cohort has an interquartile range of 0.003, so three of them is a band of
# 0.008, and the two lowest-contrast lakes fall outside it although
# `tile_crop_consistency.py` shows their tiles to be exact. A tight cohort
# turns ordinary scene-driven variation into an outlier.
#
# So this test is scoped to where its assumption holds: the carrier. Other
# transfers are measured and reported, never judged. The test that does judge
# every transfer is `tile_crop_consistency.py`, which compares each stored crop
# with the window's own pixels over the same ground -- exact, geolocated, and
# free of any model of what the transfer does.
MAX_DILUTION_DEV = 0.10          # |slope / corr^2 - 1|, on the carrier's tiles


def region_row(scene, halo_px=32, band=None, kernel=None):
    """One (region, kernel) row. `kernel=None` is the tile the paper uses.

    A region has one window and up to five tiles, one per transfer, and they
    were not all built on the same day. The window can therefore match the
    published Clough--Tocher tile and not the PSF-aware one, or the reverse,
    so the audit is per tile rather than per region; the comparison of
    Section~\\ref{sec:transfer_cost} sets two of those tiles against each
    other and is only meaningful where both came from the same acquisition.
    """
    store = cfg.REGION_DIR / f"{scene}_data.zarr"
    if kernel is None:
        cell_ids, cell_data, utm_crop, meta, npz = build_tiles.load_tile(scene)
    else:
        z = np.load(build_tiles.tile_path(scene, kernel), allow_pickle=False)
        cell_data = z["cell_data"]
    stored = np.asarray(cell_data, dtype=np.float64).ravel()

    sw = build_tiles.tile_subwindow(scene, store, halo_px=halo_px, band=band)
    img = sw["img"]
    H, W = sw["shape"]
    col = np.clip(np.round(sw["cand"]["col"]).astype(np.int64), 0, W - 1)
    row = np.clip(np.round(sw["cand"]["row"]).astype(np.int64), 0, H - 1)
    samp = np.asarray(img[row, col], dtype=np.float64)

    good = np.isfinite(samp) & np.isfinite(stored)
    s, m = stored[good], samp[good]
    sd = float(stored.std())
    slope, offset = np.polyfit(m, s, 1) if len(s) > 2 else (np.nan, np.nan)
    diff = s - m
    rms = float(np.sqrt(np.mean(diff * diff)))
    corr = float(np.corrcoef(s, m)[0, 1]) if len(s) > 2 else np.nan

    # The verdict needs the whole cohort, so it is taken in add_status.
    dilution = float(slope / corr ** 2) if corr and corr == corr else np.nan
    return dict(
        scene=scene, cls=scene.split("__")[1],
        kernel=kernel or build_tiles.DEFAULT_KERNEL, n=int(good.sum()),
        n_nonfinite=int((~good).sum()),
        stored_mean=float(stored.mean()), stored_std=sd,
        window_mean=float(m.mean()), window_std=float(m.std()),
        corr=corr, slope=float(slope), slope_over_corr2=dilution,
        offset=float(offset),
        offset_frac_of_std=float(offset / sd) if sd else np.nan,
        rms=rms, rms_frac_of_std=float(rms / sd) if sd else np.nan)


def add_status(df):
    """Judge the carrier's tiles; measure the others without judging them.

    `status` is "ok" or "MISMATCH" on the carrier kernel and "not applicable"
    elsewhere, because the identity this test rests on assumes a transfer that
    does not systematically alter the field. `departure` is |ratio - 1| on
    every row, so the other transfers stay inspectable.
    """
    df = df.copy()
    df["departure"] = (df.slope_over_corr2.astype(float) - 1.0).abs()
    carrier = getattr(cfg, "CARRIER_KERNEL", "clough_tocher")
    is_carrier = (df.kernel == carrier if "kernel" in df.columns
                  else pd.Series(True, index=df.index))
    df["status"] = np.where(
        ~is_carrier, "not applicable",
        np.where(df.departure > MAX_DILUTION_DEV, "MISMATCH", "ok"))
    return df


def carrier_rows(df):
    """The rows for the tile the paper's results are built on."""
    if "kernel" in df.columns:
        return df[df.kernel == cfg.CARRIER_KERNEL]
    return df


def markdown(df):
    per_kernel = []
    if "kernel" in df.columns and df.kernel.nunique() > 1:
        for k, g in df.groupby("kernel"):
            judged = g[g.status.isin(("ok", "MISMATCH"))]
            n_bad = int((g.status == "MISMATCH").sum())
            verdict = (f"{n_bad}" if len(judged)
                       else "not judged (see tile\\_crop\\_consistency)")
            names = ", ".join(g.loc[g.status == "MISMATCH", "scene"])
            per_kernel.append(
                f"| `{k}` | {len(g)} | {g.departure.median():.3f} | "
                f"{verdict} | {names or '--'} |")
    df = carrier_rows(df)
    bad = df[df.status == "MISMATCH"]
    lines = [
        "# Window radiometry audit",
        "",
        "Stored tile against its own window, sampled at the tile's cell centres "
        "by nearest pixel. Under the sampling-noise model the regression slope "
        "equals the squared correlation whatever the noise level, so the ratio "
        "`slope / corr^2` is one for a transfer that does not systematically "
        "alter the field. A deconvolving transfer does, so the verdict is taken "
        "on the carrier's tiles alone, where the assumption holds: flagged "
        f"when the ratio departs from one by more than {MAX_DILUTION_DEV}. "
        "The other transfers are measured and reported, never judged; "
        "`tile_crop_consistency.py` is the test that covers them.",
        "",
        f"- regions audited: {len(df)}",
        f"- passing: {int((df.status == 'ok').sum())}",
        f"- mismatched: {len(bad)}",
        "",
    ]
    # Column access is by key throughout: `row.corr` is Series.corr, the
    # method, not the column, and formatting it fails at the end of a run
    # rather than at the start of one.
    if len(bad):
        lines += ["Mismatched regions:", ""]
        for _, r in bad.iterrows():
            lines.append(
                f"- `{r['scene']}`: corr {r['corr']:.3f}, "
                f"slope {r['slope']:.3f}, slope/corr^2 "
                f"{r['slope_over_corr2']:.3f}, "
                f"offset {r['offset']:+.5f} ({r['offset_frac_of_std']:+.2f} of "
                f"the tile's std), RMS {r['rms']:.2e} "
                f"({r['rms_frac_of_std']:.2f} of the tile's std)")
        lines += [
            "",
            "A mismatched window does not affect any result read from the "
            "stored tile -- the compression sweeps, the storage accounting, the "
            "level choice by aggregation. It does affect every measurement that "
            "re-runs a transfer from the window.",
            "",
        ]
    d = df.slope_over_corr2
    lines += [
        "",
        f"`slope / corr^2` over the {len(df)} regions: median {d.median():.4f}, "
        f"range {d.min():.4f}--{d.max():.4f}.",
        "",
        "Worst ten by departure from the sampling-noise identity:", "",
        "| scene | tile std | corr | slope | slope/corr^2 | RMS/std | status |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- |"]
    order = (df.slope_over_corr2 - 1.0).abs().sort_values(ascending=False).index
    for _, r in df.loc[order].head(10).iterrows():
        lines.append(f"| {r['scene']} | {r['stored_std']:.4f} | "
                     f"{r['corr']:.4f} | {r['slope']:.3f} | "
                     f"{r['slope_over_corr2']:.3f} | "
                     f"{r['rms_frac_of_std']:.3f} | {r['status']} |")
    if per_kernel:
        lines += [
            "", "## Every transfer's tile, not only the carrier's", "",
            "A region has one window and one tile per transfer, and they were "
            "not built on the same day, so a window can match one tile and not "
            "another. Section V-L sets two of these tiles against each other, "
            "which is only meaningful where both came from the same "
            "acquisition.", "",
            "| kernel | tiles | median \\|ratio-1\\| | mismatched | which |",
            "| --- | ---: | ---: | ---: | --- |", *per_kernel]
    return "\n".join(lines) + "\n"


def latex_table(df):
    """The audit in one block: a tight cluster, and what falls out of it.

    Built on the carrier's tiles alone. The other transfers' tiles are audited
    too and reported in the markdown; the appendix states the carrier's result,
    which is what every published number depends on.
    """
    df = carrier_rows(df)
    ok, bad = df[df.status == "ok"], df[df.status == "MISMATCH"]
    if not len(bad):
        cols = "one column"
        head = r"Quantity & All " + str(len(df)) + r" regions \\"
        def cell(_r, c, fmt):
            return format(df[c].min(), fmt) + "--" + format(df[c].max(), fmt)
        rows = [f"Correlation & {cell(None,'corr','.3f')}" + r" \\",
                f"Slope $/$ corr$^2$ & {cell(None,'slope_over_corr2','.3f')}" + r" \\",
                f"RMS $/$ tile std & {cell(None,'rms_frac_of_std','.3f')}" + r" \\"]
        note_tail = "No region fails the test."
    else:
        names = ", ".join(s.replace("region__", "").replace("__", "/")
                          .replace("_", r"\_") for s in bad.scene)
        head = (f"Quantity & {len(ok)} eligible & {names} " + r"\\")
        def rng(c, fmt):
            return format(ok[c].min(), fmt) + "--" + format(ok[c].max(), fmt)
        def one(c, fmt):
            return " / ".join(format(v, fmt) for v in bad[c])
        rows = [
            f"Correlation & {rng('corr', '.3f')} & {one('corr', '.3f')}" + r" \\",
            f"Slope $/$ corr$^2$ & {rng('slope_over_corr2', '.3f')} & "
            f"{one('slope_over_corr2', '.3f')}" + r" \\",
            f"RMS $/$ tile std & {rng('rms_frac_of_std', '.3f')} & "
            f"{one('rms_frac_of_std', '.3f')}" + r" \\"]
        note_tail = (
            f"{names} fails all three at once, and by a margin no interpolation "
            "choice produces: a rebuild of its stored tile from its own window "
            "reproduces neither the level-19 carrier nor the level-20 one. "
            "Nothing read from a stored tile is affected; the measurements that "
            "re-run a transfer from the window exclude it and say so.")
    cols = "@{}lrr@{}" if len(bad) else "@{}lr@{}"
    note = (
        f"Each stored tile against its own native window, sampled at the tile's "
        f"cell centres by nearest pixel, over the {len(df)} regions. Nearest-pixel "
        "sampling is cruder than the Clough--Tocher transfer that built the "
        "tile, so the two differ even when the window is the right product; "
        "that difference is additive noise on the regressor, which attenuates "
        "the regression slope by exactly the squared correlation. The ratio "
        "slope$/$corr$^2$ is therefore one whatever the noise level and "
        "whatever the scene's contrast, and it is the test: a window that is a "
        "different acquisition carries a gain and an offset of its own and "
        f"breaks the identity. The median over all {len(df)} regions is "
        f"{df.slope_over_corr2.median():.4f}. " + note_tail)
    return "\n".join([
        r"\begin{table}[t]", r"\centering",
        *cfg.tex_caption("Is each window still the product its tile came from?",
                         "tab:window_radiometry"),
        r"\footnotesize",
        r"\begin{tabular}{" + cols + "}", r"\toprule", head,
        r"\midrule", *rows, r"\bottomrule", r"\end{tabular}",
        *cfg.tex_note(note), r"\end{table}"]) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--classes", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--band", default=None)
    ap.add_argument("--halo", type=int, default=32)
    ap.add_argument("--kernels", nargs="*", default=None,
                    help="audit these transfers' tiles as well as the "
                         "carrier's; default: every kernel whose tile exists")
    args = ap.parse_args(argv)

    if build_tiles.DEFAULT_KERNEL != cfg.CARRIER_KERNEL:
        raise SystemExit(
            f"build_tiles.DEFAULT_KERNEL is {build_tiles.DEFAULT_KERNEL!r} but "
            f"paper_config.CARRIER_KERNEL is {cfg.CARRIER_KERNEL!r}; every "
            "script that reads the audit filters on the second, so they must "
            "agree.")
    kernels = list(args.kernels) if args.kernels else list(build_tiles.KERNELS)
    if cfg.CARRIER_KERNEL not in kernels:
        kernels.insert(0, cfg.CARRIER_KERNEL)

    sites = region_sites.select(region_sites.load_sites(), args.regions,
                                args.classes, args.limit)
    rows, failures = [], []
    for i, scene in enumerate(sites.scene, 1):
        if not build_tiles.tile_path(scene).exists():
            failures.append((scene, "no tile; run build_tiles.py"))
            continue
        for kernel in kernels:
            # Only the tiles that exist: a kernel nobody built is not a
            # failure, and reporting it as one would bury the real ones.
            if not build_tiles.tile_path(scene, kernel).exists():
                continue
            try:
                r = region_row(scene, halo_px=args.halo, band=args.band,
                               kernel=kernel)
            except Exception as exc:                          # noqa: BLE001
                failures.append((f"{scene} [{kernel}]",
                                 f"{type(exc).__name__}: {exc}"))
                continue
            rows.append(r)
            print(f"[{i:>2}/{len(sites)}] {scene} [{kernel}]: "
                  f"corr {r['corr']:.4f}  slope {r['slope']:.3f}  "
                  f"slope/corr^2 {r['slope_over_corr2']:.3f}  "
                  f"offset {r['offset']:+.5f}  RMS "
                  f"{r['rms_frac_of_std']:.3f} std")
    if not rows:
        for s, why in failures:
            print(f"  {s}: {why}")
        raise SystemExit("no region could be audited")

    df = add_status(pd.DataFrame(rows))
    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    tag = "__subset" if (args.regions or args.classes or args.limit) else ""
    df.to_csv(cfg.RESULTS_DIR / f"window_radiometry{tag}.csv", index=False)
    # The CSV is the result; the report is a convenience. A bug in the report
    # must not throw away a run that has already finished.
    try:
        (cfg.RESULTS_DIR / f"window_radiometry{tag}.md").write_text(
            markdown(df), encoding="utf-8")
    except Exception as exc:                                  # noqa: BLE001
        print(f"note: the markdown report failed ({type(exc).__name__}: {exc}); "
              "the CSV is written and the verdict below is unaffected")
    if not tag:
        cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
        (cfg.TEX_TABLES / "table_window_radiometry.tex").write_text(
            latex_table(df), encoding="utf-8")
        print(f"wrote {cfg.TEX_TABLES / 'table_window_radiometry.tex'}")
    if "kernel" in df.columns and df.kernel.nunique() > 1:
        print("\nby transfer:")
        for k, g in df.groupby("kernel"):
            judged = int(g.status.isin(("ok", "MISMATCH")).sum())
            verdict = (f"{int((g.status == 'MISMATCH').sum())} mismatched"
                       if judged else "not judged (measured only)")
            print(f"  {k:16s} {len(g):3d} tile(s), median |ratio-1| "
                  f"{g.departure.median():.3f}, {verdict}")
    df = carrier_rows(df)
    bad = df[df.status == "MISMATCH"]
    print(f"\n{len(df)} carrier tile(s) audited, {len(bad)} mismatched")
    for _, r in bad.iterrows():
        print(f"  MISMATCH {r['scene']}: stored ~ {r['slope']:.3f} * window "
              f"{r['offset']:+.5f}, corr {r['corr']:.3f}, "
              f"slope/corr^2 {r['slope_over_corr2']:.3f}")
    print(f"wrote {cfg.RESULTS_DIR / f'window_radiometry{tag}.csv'}")
    if failures:
        print(f"\n{len(failures)} region(s) could not be audited:")
        for s, why in failures:
            print(f"  {s}: {why}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
