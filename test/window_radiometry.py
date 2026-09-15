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
MAX_DILUTION_DEV = 0.10          # |slope / corr^2 - 1|
MIN_CORR = 0.80                  # a backstop, an order of magnitude from the data


def region_row(scene, halo_px=32, band=None):
    store = cfg.REGION_DIR / f"{scene}_data.zarr"
    cell_ids, cell_data, utm_crop, meta, npz = build_tiles.load_tile(scene)
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

    dilution = float(slope / corr ** 2) if corr and corr == corr else np.nan
    ok = (corr > MIN_CORR and abs(dilution - 1.0) < MAX_DILUTION_DEV)
    return dict(
        scene=scene, cls=scene.split("__")[1], n=int(good.sum()),
        n_nonfinite=int((~good).sum()),
        stored_mean=float(stored.mean()), stored_std=sd,
        window_mean=float(m.mean()), window_std=float(m.std()),
        corr=corr, slope=float(slope), slope_over_corr2=dilution,
        offset=float(offset),
        offset_frac_of_std=float(offset / sd) if sd else np.nan,
        rms=rms, rms_frac_of_std=float(rms / sd) if sd else np.nan,
        status="ok" if ok else "MISMATCH")


def markdown(df):
    bad = df[df.status != "ok"]
    lines = [
        "# Window radiometry audit",
        "",
        "Stored tile against its own window, sampled at the tile's cell centres "
        "by nearest pixel. Under the sampling-noise model the regression slope "
        "equals the squared correlation whatever the noise level, so a region "
        f"passes when `slope / corr^2` is within {MAX_DILUTION_DEV} of one "
        f"(and the correlation clears {MIN_CORR} as a backstop).",
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
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--classes", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--band", default=None)
    ap.add_argument("--halo", type=int, default=32)
    args = ap.parse_args(argv)

    sites = region_sites.select(region_sites.load_sites(), args.regions,
                                args.classes, args.limit)
    rows, failures = [], []
    for i, scene in enumerate(sites.scene, 1):
        if not build_tiles.tile_path(scene).exists():
            failures.append((scene, "no tile; run build_tiles.py"))
            continue
        try:
            r = region_row(scene, halo_px=args.halo, band=args.band)
        except Exception as exc:                              # noqa: BLE001
            failures.append((scene, f"{type(exc).__name__}: {exc}"))
            continue
        rows.append(r)
        print(f"[{i:>2}/{len(sites)}] {scene}: corr {r['corr']:.4f}  slope "
              f"{r['slope']:.3f}  slope/corr^2 {r['slope_over_corr2']:.3f}  "
              f"offset {r['offset']:+.5f}  RMS "
              f"{r['rms_frac_of_std']:.3f} std  -> {r['status']}")
    if not rows:
        for s, why in failures:
            print(f"  {s}: {why}")
        raise SystemExit("no region could be audited")

    df = pd.DataFrame(rows)
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
    bad = df[df.status != "ok"]
    print(f"\n{len(df)} region(s) audited, {len(bad)} mismatched")
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
