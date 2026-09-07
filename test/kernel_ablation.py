"""Ablation of the transfer kernel that populates the carrier grid (V3).

Section~\\ref{sec:resampler} of the paper reports, qualitatively, that an
inverse-distance-weighted resampler left excess high-frequency power in the
carrier tile and that replacing it with Clough--Tocher removed most of it. That
observation decided the pipeline but was never measured: no spectral distance
was computed, and no rate--distortion consequence was quantified. Everything
the paper says about the kernel therefore rests on a visual comparison, which
is why the review asked for either the ablation or the removal of the claims.

This runs the ablation. For each kernel in

    nearest, linear, bicubic, clough_tocher, psf

the carrier tile is rebuilt from the same native window, and three quantities
are measured against the same observation:

    floor        the RMSE of the uncompressed tile read on the native pixels,
                 over the two footprints' intersection. This is the paper's
                 representation floor, now as a function of the kernel. It is
                 the transfer error and the reading error together, but the
                 read is identical across kernels, so differences between rows
                 are the transfer alone -- which is the decomposition the
                 operator characterisation of `operator_tests.py` makes
                 possible;

    log-PSD      the radially averaged power spectrum of the projected tile
                 against that of the observed crop, as an RMS distance in
                 log10 over the band, plus the ratio in the top octave alone.
                 This is the quantitative spectral metric the paper says it
                 never computed, and it is what tells excess high-frequency
                 content from a deficit;

    rate         with --rd, the compressed size at each target RMSE, so the
                 kernel's effect on the headline ratios is bounded rather than
                 assumed. This is the expensive half and is off by default.

Only the default kernel writes to the unsuffixed tile path, so nothing the
paper's results are built on is overwritten.

Usage
-----
    python kernel_ablation.py --limit 3                # floor + spectrum
    python kernel_ablation.py --limit 3 --rd           # + rate-distortion
    python kernel_ablation.py --kernels nearest psf --limit 2

Writes results/kernel_ablation.csv, results/kernel_ablation.md,
figures/fig_kernel_psd.pdf and ../tex/tables/table_kernel_ablation.tex.
"""

from __future__ import annotations

import argparse
import math
import shutil
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt           # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402

import paper_config as cfg                # noqa: E402
import paper_stats as st                  # noqa: E402
import build_tiles                        # noqa: E402
import common_reference                   # noqa: E402
import region_sites                       # noqa: E402
import support_controls as sup            # noqa: E402

REFERENCE_KERNEL = build_tiles.DEFAULT_KERNEL


def ensure_tile(scene: str, kernel: str, store: Path, build: bool, force=False):
    """Return the tile path, building it if it is not there yet."""
    path = build_tiles.tile_path(scene, kernel)
    if path.exists() and not force:
        return path, 0.0
    if not build:
        raise FileNotFoundError(f"{path} is missing and --no-build was given")
    if not store.exists():
        raise FileNotFoundError(f"missing window store {store}")
    t0 = time.perf_counter()
    build_tiles.build_tile(scene, store, kernel=kernel, verbose=False)
    return path, time.perf_counter() - t0


def spectral_distance(projected, observed, mask):
    """Log-PSD distance between a projected tile and the observation.

    The comparison is made on the largest rectangle inscribed in the two
    footprints' intersection: a spectrum taken over a masked octagon would
    measure the mask's edges rather than the field, and the whole point of this
    metric is to be about content.
    """
    rect = sup.largest_true_rectangle(mask)
    if rect is None:
        return dict(logpsd_rms=np.nan, logpsd_top_octave=np.nan,
                    psd_rect_px=0)
    r0, r1, c0, c1 = rect
    n = min(r1 - r0, c1 - c0)
    p = np.asarray(projected)[r0:r0 + n, c0:c0 + n]
    o = np.asarray(observed)[r0:r0 + n, c0:c0 + n]
    if not (np.all(np.isfinite(p)) and np.all(np.isfinite(o))):
        p = np.nan_to_num(p, nan=float(np.nanmean(p)))
        o = np.nan_to_num(o, nan=float(np.nanmean(o)))
    k, psd_p = sup.radial_psd(p, cfg.NATIVE_PIXEL_M)
    _k, psd_o = sup.radial_psd(o, cfg.NATIVE_PIXEL_M)
    k_nyq = math.pi / cfg.NATIVE_PIXEL_M
    band = (k > 0) & (k <= k_nyq) & (psd_p > 0) & (psd_o > 0)
    if band.sum() < 6:
        return dict(logpsd_rms=np.nan, logpsd_top_octave=np.nan, psd_rect_px=n)
    d = np.log10(psd_p[band] / psd_o[band])
    top = (k >= 0.5 * k_nyq) & band
    return dict(
        logpsd_rms=float(np.sqrt(np.mean(d ** 2))),
        # >1: the transfer leaves more high-frequency power than the
        # observation has, which is the failure mode the paper describes
        logpsd_top_octave=(float(psd_p[top].sum() / psd_o[top].sum())
                           if top.any() else np.nan),
        psd_rect_px=int(n))



def transfer_fidelity(scene: str, kernel: str, store: Path, seed: int = 17,
                      halo_px: int = 32):
    """Transfer error against an analytic truth, measured at the cell centres.

    Scoring a transfer by comparing the tile back to the source imagery --
    which is what the floor does -- rewards any transfer that is close to the
    identity on the source lattice, whatever it does to the continuous field.
    A nearest transfer copies native pixel values, and the nearest read then
    hands most of them straight back, so the two nearly cancel: `nearest` wins
    that comparison by construction of the metric and not by fidelity. The
    log-PSD distance has the same weakness, since a copy reproduces the source
    spectrum exactly.

    The truth here is an analytic field, known at every position, with the
    slope and instrument width fitted to the imagery in the support controls.
    It is evaluated at the native pixel centres to play the observation,
    transferred to the cells by the kernel under test, and compared to the same
    field evaluated at the cell centres. Neither lattice is privileged and the
    identity has no advantage, so this is the number that says whether a kernel
    represents the field or merely copies samples.
    """
    import pyproj
    ds, img, lon, lat, x, y, crs_wkt = build_tiles._open_window(store, cfg.BAND)
    H, W = img.shape
    crs = pyproj.CRS.from_wkt(crs_wkt)
    to_crs = pyproj.Transformer.from_crs(pyproj.CRS.from_epsg(4326), crs,
                                         always_xy=True)
    # The parent is whichever one the paper's own tile used, so the geometry is
    # identical to the real-data rows above.
    _ids, _cd, _crop, _meta, npz = build_tiles.load_tile(scene)
    parent = int(npz["parent_id"])
    cand = build_tiles._evaluate_parent(parent, to_crs, x, y, img, H, W, halo_px)

    c0 = max(int(np.floor(cand["col"].min())) - halo_px, 0)
    c1 = min(int(np.ceil(cand["col"].max())) + halo_px + 1, W)
    r0 = max(int(np.floor(cand["row"].min())) - halo_px, 0)
    r1 = min(int(np.ceil(cand["row"].max())) + halo_px + 1, H)
    sub_lon = lon[r0:r1, c0:c1].ravel()
    sub_lat = lat[r0:r1, c0:c1].ravel()
    sx, sy = to_crs.transform(sub_lon, sub_lat)

    x0 = float(np.mean(cand["cx"]))
    y0 = float(np.mean(cand["cy"]))
    k_min = 2.0 * math.pi / cfg.TILE_SIDE_M
    k_max = math.pi / cfg.DELTA_EQ
    sigma = sup.psf_sigma_from_fwhm(11.9)
    f_src = sup.synthetic_field(np.asarray(sx) - x0, np.asarray(sy) - y0,
                                beta=1.9, k_min=k_min, k_max=k_max, seed=seed,
                                psf_sigma_m=sigma)
    f_truth = sup.synthetic_field(np.asarray(cand["cx"]) - x0,
                                  np.asarray(cand["cy"]) - y0,
                                  beta=1.9, k_min=k_min, k_max=k_max, seed=seed,
                                  psf_sigma_m=sigma)
    op = build_tiles.make_resampler(kernel, sub_lon, sub_lat, cand["cell_ids"],
                                    verbose=False)
    res = op.resample(np.asarray(f_src, dtype=np.float64))
    ids = np.asarray(res.cell_ids, dtype=np.int64).ravel()
    val = np.asarray(res.cell_data, dtype=np.float64).ravel()
    # Same completion as the real-data tiles, so a kernel is not judged on the
    # handful of hull-edge cells the library drops rather than on its transfer.
    try:
        val, n_filled = build_tiles.complete_block(ids, val, cand["cell_ids"])
    except RuntimeError:
        return dict(transfer_rmse=np.nan, transfer_frac_of_std=np.nan,
                    transfer_p999=np.nan, transfer_overshoot=np.nan,
                    transfer_cells_kept=int(ids.size), transfer_cells_filled=-1)
    d = val - np.asarray(f_truth, dtype=np.float64)
    sd = float(np.std(f_truth)) or 1.0
    lo, hi = float(np.min(f_src)), float(np.max(f_src))
    return dict(
        transfer_rmse=float(np.sqrt(np.mean(d ** 2))),
        transfer_frac_of_std=float(np.sqrt(np.mean(d ** 2)) / sd),
        # a diverging interpolant shows up in the tail long before the RMSE
        transfer_p999=float(np.percentile(np.abs(d), 99.9) / sd),
        # values the transfer invented outside the range of what it was given
        transfer_overshoot=float(np.mean((val < lo) | (val > hi))),
        transfer_cells_kept=int(ids.size),
        transfer_cells_filled=int(n_filled))

def run_region(scene: str, kernel: str, store: Path, codecs, targets,
               build=True, force=False, do_rd=False, truth=True, verbose=True):
    path, build_s = ensure_tile(scene, kernel, store, build, force)
    cell_ids, cell_data, utm_crop, meta, npz = build_tiles.load_tile(scene, kernel)
    cr = common_reference.CommonReference(cell_ids, np.asarray(npz["utm_lon"]),
                                          np.asarray(npz["utm_lat"]),
                                          mode="nearest")
    proj = cr.predict(cell_data)
    floor = cr.rmse_on_mask(proj, utm_crop)
    obs = np.asarray(utm_crop, dtype=np.float64)
    scene_std = float(np.std(obs[cr.mask]))

    # A kernel that diverges shows up in the tail and in the range long before
    # the RMSE makes it obvious, and the constant-field check passes it: a
    # bicubic interpolant is exact on constants and can still overshoot wildly
    # on real data. Both diagnostics travel with every row.
    resid = (np.asarray(proj, dtype=np.float64) - obs)[cr.mask]
    resid = resid[np.isfinite(resid)]
    lo, hi = float(np.min(obs[cr.mask])), float(np.max(obs[cr.mask]))
    cd = np.asarray(cell_data, dtype=np.float64)
    row = dict(scene=scene, kernel=kernel, floor=floor,
               floor_frac_of_std=floor / scene_std if scene_std else np.nan,
               floor_p999=(float(np.percentile(np.abs(resid), 99.9) / scene_std)
                           if resid.size and scene_std else np.nan),
               cells_outside_source_range=float(np.mean((cd < lo) | (cd > hi))),
               scene_std=scene_std, build_seconds=build_s,
               constant_check=meta.get("resampler_constant_check", np.nan),
               cells_filled=int(meta.get("resampler_cells_filled", 0) or 0),
               hpx_std=float(np.std(cell_data)))
    row.update(spectral_distance(proj, obs, cr.mask))
    if truth:
        try:
            row.update(transfer_fidelity(scene, kernel, store))
        except Exception as exc:                           # noqa: BLE001
            print(f"      transfer-fidelity failed: {type(exc).__name__}: {exc}")
            row.update(transfer_rmse=np.nan, transfer_frac_of_std=np.nan,
                       transfer_p999=np.nan, transfer_overshoot=np.nan,
                       transfer_cells_kept=-1, transfer_cells_filled=-1,
                       transfer_fail=f"{type(exc).__name__}: {exc}"[:200])

    if do_rd:
        tile = sup.expanded(cell_data)
        sweep = sup.sweep_array(tile, scene, codecs, label=f"hpx_{kernel}")
        for t in targets:
            row[f"size_at_{t:g}"] = st.best_of_codecs_size(sweep, scene, codecs, t)
            row[f"bracketed_at_{t:g}"] = st.any_bracketed(sweep, scene, codecs,
                                                          t, "rmse")
        sweep.to_csv(cfg.RESULTS_DIR / "kernel_ablation_sweeps.csv", mode="a",
                     header=not (cfg.RESULTS_DIR /
                                 "kernel_ablation_sweeps.csv").exists(),
                     index=False)
    if verbose:
        extra = ""
        if do_rd and np.isfinite(row.get(f"size_at_{targets[0]:g}", np.nan)):
            extra = (f" | size@{targets[0]:g} "
                     f"{row[f'size_at_{targets[0]:g}']/cfg.SHARED_RAW:.3f}")
        tf = row.get("transfer_frac_of_std", np.nan)
        flag = ""
        if row["cells_outside_source_range"] > 0.01:
            flag = (f"  <-- DIVERGES: "
                    f"{100*row['cells_outside_source_range']:.1f}% of cells "
                    f"outside the source range")
        print(f"    {kernel:14s} floor {floor:.3e} "
              f"({100*row['floor_frac_of_std']:.1f}% of std) | "
              f"log-PSD {row['logpsd_rms']:.3f} | "
              f"top octave {row['logpsd_top_octave']:.2f}x | "
              + (f"transfer {100*tf:.1f}% of std" if np.isfinite(tf)
                 else "transfer --")
              + extra + flag)
    return row


# --- reporting ---------------------------------------------------------------
def paired(df, col, reference=REFERENCE_KERNEL):
    """Per-region ratio of each kernel's column to the reference kernel's."""
    ref = df[df.kernel == reference].set_index("scene")[col]
    out = []
    for kernel, g in df.groupby("kernel"):
        v = g.set_index("scene")[col]
        common = v.index.intersection(ref.index)
        if not len(common):
            continue
        r = (v.loc[common] / ref.loc[common]).to_numpy(dtype=float)
        r = r[np.isfinite(r)]
        out.append(dict(kernel=kernel, n=len(r),
                        median=float(np.median(r)) if r.size else np.nan,
                        lo=float(np.min(r)) if r.size else np.nan,
                        hi=float(np.max(r)) if r.size else np.nan))
    return pd.DataFrame(out)


def summarize(df, targets, do_rd):
    out = ["# Transfer-kernel ablation (V3)", "",
           f"- regions: {df.scene.nunique()}, kernels: "
           f"{', '.join(sorted(df.kernel.unique()))}",
           f"- reference kernel: `{REFERENCE_KERNEL}` (the one the paper uses)",
           ""]
    has_truth = "transfer_frac_of_std" in df
    out += ["## Against an analytic truth (the metric that decides)", ""]
    if has_truth:
        out += ["One analytic field, evaluated at the native pixel centres to "
                "play the observation, transferred to the cells by each kernel, "
                "and compared to the same field evaluated at the cell centres. "
                "Neither lattice is privileged, so a transfer gets no credit "
                "for merely copying samples.", "",
                "| kernel | n | transfer RMSE / std | 99.9th pct / std | "
                "overshoot outside source range | cells kept |",
                "|---|---|---|---|---|---|"]
        # A kernel measured on a different set of regions from the others is
        # not comparable with them, and averaging over whichever regions
        # happened to succeed is the attrition this project has been removing
        # everywhere else. The support is printed, and the headline comparison
        # below is restricted to the regions every kernel reached.
        ok = df.dropna(subset=["transfer_frac_of_std"])
        common = None
        for _k, g in ok.groupby("kernel"):
            sc = set(g.scene)
            common = sc if common is None else (common & sc)
        common = common or set()
        for kernel, g in ok.groupby("kernel"):
            kept = g.transfer_cells_kept.median()
            out.append(f"| `{kernel}` | {g.scene.nunique()} | "
                       f"{100*g.transfer_frac_of_std.median():.2f}% | "
                       f"{g.transfer_p999.median():.3f} | "
                       f"{100*g.transfer_overshoot.median():.3f}% | "
                       f"{'all' if kept == cfg.N_CHILDREN else f'{kept:.0f}'}"
                       + (f" (+{g.transfer_cells_filled.median():.0f} filled)"
                          if "transfer_cells_filled" in g
                          and g.transfer_cells_filled.median() > 0 else "")
                       + " |")
        out.append("")
        n_reg = df.scene.nunique()
        if len(common) < n_reg:
            miss = sorted(set(df.scene) - common)
            out.append(f"- **support is not common**: only {len(common)} of "
                       f"{n_reg} regions were measured for every kernel "
                       f"(missing: {', '.join(m.split('__')[-1] for m in miss)}). "
                       f"The medians above are each over that kernel's own "
                       f"regions; the paired comparison below uses the common "
                       f"support only.")
            if "transfer_fail" in df:
                for _, r in df[df.transfer_fail.notna()].iterrows():
                    out.append(f"  - `{r.kernel}` on "
                               f"{r.scene.split('__')[-1]}: {r.transfer_fail}")
            out.append("")
        if common:
            sub = ok[ok.scene.isin(common)]
            out += ["Restricted to the common support:", "",
                    "| kernel | n | transfer RMSE / std |", "|---|---|---|"]
            for kernel, g in sub.groupby("kernel"):
                out.append(f"| `{kernel}` | {g.scene.nunique()} | "
                           f"{100*g.transfer_frac_of_std.median():.2f}% |")
            out.append("")
    else:
        out.append("_Not run (`--no-truth`)._\n")

    out += ["## Against the imagery (read back through the same operator)", "",
            "**These two columns reward a transfer for being close to the "
            "identity on the native lattice.** A nearest transfer copies pixel "
            "values and the nearest read hands most of them straight back, so "
            "the two nearly cancel; the log-PSD distance has the same weakness, "
            "since a copy reproduces the source spectrum exactly. Read them as "
            "a description of the transfer-and-read pair, not as a ranking of "
            "transfers.", "",
            "| kernel | floor | floor / scene std | 99.9th pct / std | log-PSD "
            "RMS | top-octave power | cells outside source range | build (s) |",
            "|---|---|---|---|---|---|---|---|"]
    for kernel, g in df.groupby("kernel"):
        out.append(f"| `{kernel}` | {g.floor.median():.3e} | "
                   f"{100*g.floor_frac_of_std.median():.1f}% | "
                   f"{g.floor_p999.median():.2f} | "
                   f"{g.logpsd_rms.median():.3f} | "
                   f"{g.logpsd_top_octave.median():.2f} | "
                   f"{100*g.cells_outside_source_range.median():.3f}% | "
                   # cached tiles record no build time; an all-cached kernel
                   # must print "cached", not the NaN of an empty median
                   + (lambda b: f"{b.median():.0f} |" if b.notna().any()
                      else "cached |")(g.build_seconds.replace(0, np.nan)))
    out.append("")
    out.append("A top-octave ratio above 1 means the transfer leaves more "
               "high-frequency power than the observation carries, which is "
               "the failure the paper describes qualitatively; below 1 means "
               "it smooths.")
    # What separates a diverging interpolant from one that merely overshoots is
    # the size of the excursion, not the number of cells involved: a per-cent of
    # cells a hair outside the source range is ordinary for any non-monotone
    # kernel. The tail is the discriminator.
    for kernel, g in df.groupby("kernel"):
        p999 = g.floor_p999.median()
        over = 100 * g.cells_outside_source_range.median()
        if np.isfinite(p999) and p999 > 5.0:
            out.append(f"- **`{kernel}` diverges**: its 99.9th-percentile error "
                       f"reaches {p999:.0f} times the scene's standard "
                       f"deviation and its floor {100*g.floor_frac_of_std.median():.0f}%"
                       f" of it, with {over:.1f}% of cells outside the range of "
                       f"the data it was given (median over all "
                       f"{g.scene.nunique()} regions). It passes the "
                       f"constant-field check, which any interpolant exact on "
                       f"polynomials will: that check is necessary and not "
                       f"sufficient.")
        elif over > 0.5:
            out.append(f"- `{kernel}` overshoots mildly: {over:.1f}% of cells "
                       f"outside the source range, but its 99.9th-percentile "
                       f"error stays at {p999:.2f} times the scene's standard "
                       f"deviation, so this is non-monotonicity and not "
                       f"divergence.")
    out.append("")

    out += ["## Paired against the reference kernel", "",
            "| kernel | n | floor ratio | range |", "|---|---|---|---|"]
    for _, r in paired(df, "floor").iterrows():
        out.append(f"| `{r.kernel}` | {int(r.n)} | {r['median']:.3f} | "
                   f"[{r.lo:.3f}, {r.hi:.3f}] |")
    out.append("")

    if do_rd:
        out += ["## Compressed size at matched RMSE, paired against the "
                "reference kernel", "",
                "| kernel | " + " | ".join(f"{t:g}" for t in targets) + " |",
                "|---" * (len(targets) + 1) + "|"]
        for kernel in sorted(df.kernel.unique()):
            cells = []
            for t in targets:
                col = f"size_at_{t:g}"
                if col not in df:
                    cells.append("--")
                    continue
                p = paired(df, col)
                v = p[p.kernel == kernel]
                cells.append(f"{v['median'].iloc[0]:.3f}" if len(v) else "--")
            out.append(f"| `{kernel}` | " + " | ".join(cells) + " |")
        out.append("")
        out.append("Values are per-region ratios of compressed size to the "
                   "reference kernel's, so 1.00 means the kernel changes "
                   "nothing at that target.")
        out.append("")

    ref = df[df.kernel == REFERENCE_KERNEL]
    others = df[df.kernel != REFERENCE_KERNEL]
    col = "transfer_frac_of_std" if has_truth else "floor"
    label = "transfer error against the analytic truth" if has_truth else "floor"
    if len(ref) and len(others) and col in df:
        allk = df.groupby("kernel")[col].median().dropna()
        m = others.groupby("kernel")[col].median().dropna()
        r = float(ref[col].median())
        if len(allk) and np.isfinite(r) and r > 0:
            best = allk.idxmin()
            if best == REFERENCE_KERNEL:
                runner = m.idxmin()
                out.append(f"- **`{REFERENCE_KERNEL}` has the lowest {label}**, "
                           f"{allk.min():.3e}; the next best is `{runner}` at "
                           f"{m.min():.3e}, {m.min()/r:.1f} times worse. The "
                           f"kernel the paper uses is the right one on this "
                           f"measure.")
            else:
                out.append(f"- **`{best}` has the lowest {label}**, "
                           f"{allk.min():.3e}, against {r:.3e} for "
                           f"`{REFERENCE_KERNEL}`: the kernel the paper uses is "
                           f"{r/allk.min():.1f} times worse on this measure.")
    return "\n".join(out)


def figure(df, publish=True, name="fig_kernel_psd"):
    """Three panels, in the order the argument runs.

    The analytic-truth panel comes first because it is the only one of the
    three that does not reward a transfer for staying close to the identity on
    the native lattice; the other two describe the transfer-and-read pair.
    """
    if not len(df):
        return
    has_truth = "transfer_frac_of_std" in df and df.transfer_frac_of_std.notna().any()
    panels = ([("transfer_frac_of_std", "transfer error\n(% of field std)", 100.0,
                "log")] if has_truth else [])
    panels += [("floor_frac_of_std", "floor\n(% of scene std)", 100.0, "log"),
               ("logpsd_top_octave", "top-octave power\nvs observation", 1.0, "log")]
    fig, axes = plt.subplots(1, len(panels), figsize=(2.4 * len(panels) + 0.6, 3.0))
    axes = np.atleast_1d(axes)
    kernels = sorted(df.kernel.unique())
    cmap = plt.get_cmap("viridis")
    colors = {k: cmap(i / max(len(kernels) - 1, 1)) for i, k in enumerate(kernels)}
    for ax, (col, ylab, scale, yscale) in zip(axes, panels):
        for i, k in enumerate(kernels):
            g = df[df.kernel == k]
            v = g[col].to_numpy(dtype=float) * scale
            v = v[np.isfinite(v)]
            if not v.size:
                continue
            ax.scatter(np.full(v.size, i), v, color=colors[k], s=18,
                       alpha=0.75, linewidths=0)
            ax.plot([i], [np.median(v)], marker="_", ms=18, mew=2.2,
                    color=colors[k])
        ax.set_xticks(range(len(kernels)))
        ax.set_xticklabels([k.replace("clough_tocher", "C--T") for k in kernels],
                           rotation=35, ha="right", fontsize=8)
        ax.set_ylabel(ylab, fontsize=9)
        ax.set_yscale(yscale)
        ax.tick_params(labelsize=8)
        ax.set_xlim(-0.6, len(kernels) - 0.4)
    if has_truth:
        axes[0].set_title("against an analytic truth", fontsize=9)
        axes[1].set_title("against the imagery", fontsize=9)
    axes[-1].axhline(1.0, color="0.6", ls="--", lw=1.0)
    fig.suptitle("Transfer kernel: fidelity, representation floor and spectrum",
                 fontsize=10)
    fig.tight_layout()
    cfg.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(cfg.FIGURES_DIR / f"{name}.{ext}", bbox_inches="tight", dpi=160)
    if publish:
        cfg.TEX_FIGURES.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(cfg.FIGURES_DIR / f"{name}.pdf",
                        cfg.TEX_FIGURES / f"{name}.pdf")
    plt.close(fig)
    print(f"  {name}")


def latex_table(df, targets, do_rd):
    rows = []
    fp = paired(df, "floor").set_index("kernel")
    for kernel, g in df.groupby("kernel"):
        ratio = (f"{fp.loc[kernel,'median']:.2f}" if kernel in fp.index else "--")
        cells = ""
        if do_rd:
            for t in targets[:3]:
                col = f"size_at_{t:g}"
                p = paired(df, col).set_index("kernel") if col in df else None
                cells += (f" & {p.loc[kernel,'median']:.2f}"
                          if p is not None and kernel in p.index else " & --")
        tf = (f"{100*g.transfer_frac_of_std.median():.2f}"
              if "transfer_frac_of_std" in g
              and np.isfinite(g.transfer_frac_of_std.median()) else "--")
        ovr = (f"{100*g.cells_outside_source_range.median():.2f}"
               if "cells_outside_source_range" in g else "--")
        rows.append(
            f"\\texttt{{{kernel.replace('_', chr(92)+'_')}}} & {tf} & {ovr} & "
            f"{100*g.floor_frac_of_std.median():.1f} & "
            f"{ratio} & {g.logpsd_rms.median():.2f} & "
            f"{g.logpsd_top_octave.median():.2f}" + cells + " \\\\")
    extra_h = "".join(f" & ${t:g}$" for t in targets[:3]) if do_rd else ""
    spec = "lrrrrrr" + ("r" * (3 if do_rd else 0))
    return (
        "\\begin{table}[t]\n\\centering\n"
        "\\caption{Transfer-kernel ablation. Each row rebuilds the carrier tile "
        "from the same native window with a different kernel, with the same "
        "reading operator throughout, so differences between rows are the "
        "transfer alone. \\emph{Transfer} is the error against an analytic "
        "field known at every position --- evaluated at the native pixels to "
        "play the observation, transferred to the cells, compared to itself at "
        "the cell centres --- as a percentage of the field's standard "
        "deviation; it is the only column here that gives a transfer no credit "
        "for merely copying samples. \\emph{Over} is the percentage of cells "
        "the transfer places outside the range of the data it was given, which "
        "detects a diverging interpolant that the constant-field check passes. "
        "\\emph{Floor} is the RMSE of the uncompressed tile read on the native "
        "pixels over the two footprints' intersection, as a fraction of the "
        "scene's standard deviation, and \\emph{rel.} its per-region ratio to "
        "\\texttt{clough\\_tocher}, the kernel used throughout this paper; both "
        "reward a transfer for staying close to the identity on the native "
        "lattice and are reported for that reason, not as a ranking. "
        "\\emph{log-PSD} is the RMS distance in $\\log_{10}$ between the "
        "projected tile's radial power spectrum and the observation's over the "
        "band up to the native Nyquist, and \\emph{top oct.} the ratio of "
        "their power in the highest octave alone: above 1 the transfer leaves "
        "more high-frequency power than the observation carries, which is the "
        "failure this paper reports qualitatively for its first resampler."
        + (" The last columns are per-region compressed-size ratios to the "
           "same reference at three RMSE targets." if do_rd else "") + "}\n"
        "\\label{tab:kernel_ablation}\n\\footnotesize\n"
        f"\\begin{{tabular}}{{@{{}}{spec}@{{}}}}\n\\toprule\n"
        "Kernel & Transfer & Over & Floor & rel. & log-PSD & top oct." + extra_h
        + " \\\\\n\\midrule\n" + "\n".join(rows)
        + "\n\\bottomrule\n\\end{tabular}\n\\end{table}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kernels", nargs="*", default=list(build_tiles.KERNELS))
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=3)
    ap.add_argument("--codecs", nargs="*", default=cfg.COMPRESSOR_CODECS)
    ap.add_argument("--targets", type=float, nargs="*", default=cfg.TARGET_RMSE)
    ap.add_argument("--rd", action="store_true",
                    help="also run the codec sweep per kernel (the expensive "
                         "half; without it only the floor and the spectrum are "
                         "measured)")
    ap.add_argument("--no-truth", action="store_true",
                    help="skip the analytic-truth transfer metric (which is "
                         "the only one of the three that does not reward a "
                         "transfer for being close to the identity)")
    ap.add_argument("--no-build", action="store_true",
                    help="fail instead of building a missing tile")
    ap.add_argument("--force", action="store_true",
                    help="rebuild tiles even when they exist")
    ap.add_argument("--no-publish", action="store_true")
    ap.add_argument("--no-tex", action="store_true")
    args = ap.parse_args(argv)

    sites = region_sites.load_sites()
    scenes = args.regions or list(sites["scene"])
    scenes = [s for s in scenes
              if (cfg.REGION_DIR / f"{s}_data.zarr").exists()
              or build_tiles.tile_path(s).exists()]
    if not scenes:
        print(f"no region windows under {cfg.REGION_DIR} and no tiles under "
              f"{cfg.TILE_DIR}")
        return 1
    if args.limit:
        scenes = scenes[:args.limit]

    print(f"{len(scenes)} region(s), kernels {args.kernels}, "
          f"rate-distortion {'on' if args.rd else 'off'}")
    rows, failures = [], []
    for scene in scenes:
        print(f"  {scene}")
        store = cfg.REGION_DIR / f"{scene}_data.zarr"
        for kernel in args.kernels:
            try:
                rows.append(run_region(scene, kernel, store, args.codecs,
                                       args.targets, build=not args.no_build,
                                       force=args.force, do_rd=args.rd,
                                       truth=not args.no_truth))
            except Exception as exc:                       # noqa: BLE001
                print(f"    {kernel:14s} FAILED {type(exc).__name__}: {exc}")
                failures.append((scene, kernel, f"{type(exc).__name__}: {exc}"))
    if not rows:
        print("no (region, kernel) pair produced a row")
        return 1

    df = pd.DataFrame(rows)
    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    csv = cfg.RESULTS_DIR / "kernel_ablation.csv"
    df.to_csv(csv, index=False)
    md = summarize(df, args.targets, args.rd)
    (cfg.RESULTS_DIR / "kernel_ablation.md").write_text(md, encoding="utf-8")
    print("\n" + md)
    print(f"wrote {csv} ({len(df)} rows)")

    figure(df, publish=not args.no_publish)
    if not args.no_tex:
        out = cfg.TEX_TABLES / "table_kernel_ablation.tex"
        cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
        out.write_text(latex_table(df, args.targets, args.rd), encoding="utf-8")
        print(f"wrote {out}")

    if failures:
        print(f"\n{len(failures)} (region, kernel) pair(s) failed:")
        for s, k, why in failures:
            print(f"  {s} {k}: {why}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
