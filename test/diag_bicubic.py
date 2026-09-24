"""Diagnose the bicubic row of the kernel ablation (revision 7, point 7).

Table XIII(c) reports a transfer error of about 392% of the analytic field's
standard deviation for `healpix_resample.BicubicResampler`, with only ~1.3% of
cells outside the source range and a 99.9th percentile near 47 standard
deviations. Those three numbers together mean that a small fraction of cells
carries enormous values, which is not what the ringing of a Keys kernel does
(it overshoots by a fraction of the local contrast). The hypothesis to test,
not to assume: the radial signed kernel's per-cell normalisation divides by a
signed weight sum that is close to zero at some spatial phases, and the
`_floor_signed` guard (10^-3 of the absolute weight) turns that into an
amplification of up to a thousand.

The script never modifies the operator. It runs, per region:

  1. an inventory of the operator's attributes (normalisers, Npt, sigma_m,
     floor), whatever the installed version exposes;
  2. four probe fields on the same geometry -- a constant, an affine plane,
     a single-frequency wave, and the ablation's own analytic field -- with
     RMSE / std, the 99.9th percentile, the overshoot fraction and the output
     extrema;
  3. the location of the worst cells: their share, whether they cluster near
     the hull edge, and (when the operator exposes it) the raw normaliser at
     those cells against the others;
  4. a small parameter sweep over sigma_m and Npt where the constructor
     accepts them;
  5. the same analytic field transferred by a classical separable bicubic on
     the regular source raster (scipy.ndimage.map_coordinates, order 3), so
     that the radial variant is compared with the interpolation it stands for.

Usage (on the machine where healpix_resample and the region stores are):

    python diag_bicubic.py                       # three regions, one per class
    python diag_bicubic.py --scenes region__urban__paris ...
    python diag_bicubic.py --all

Writes results/diag_bicubic.md and results/diag_bicubic.csv. Nothing else in
the pipeline reads them; the paper's prose on the bicubic row is to be
rewritten from this report once it exists.
"""
from __future__ import annotations

import argparse
import inspect
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import paper_config as cfg
import build_tiles
import region_sites
import support_controls as sup

DEFAULT_SCENES = None   # resolved at run time: the first region of each class


# --- helpers ---------------------------------------------------------------
def _attr_inventory(op) -> dict:
    """Everything the operator exposes that could explain a normalisation."""
    keys = ("norm", "weight", "sigma", "Npt", "npt", "floor", "knn", "q",
            "kernel", "mode", "radial", "eps")
    out = {}
    for name in dir(op):
        if name.startswith("__"):
            continue
        if any(k.lower() in name.lower() for k in keys):
            try:
                v = getattr(op, name)
            except Exception as e:  # noqa: BLE001
                out[name] = f"<error {e}>"
                continue
            if callable(v):
                try:
                    out[name] = f"<callable {inspect.signature(v)}>"
                except (TypeError, ValueError):
                    out[name] = "<callable>"
            elif hasattr(v, "shape"):
                arr = np.asarray(v.detach().cpu() if hasattr(v, "detach") else v)
                out[name] = (f"array shape={arr.shape} dtype={arr.dtype} "
                             f"min={np.nanmin(arr):.3e} med={np.nanmedian(arr):.3e} "
                             f"max={np.nanmax(arr):.3e}")
            else:
                out[name] = repr(v)[:120]
    return out


def _normaliser(op):
    """Return (raw, floored) per-cell signed normalisers as numpy, or None."""
    for raw_name, fl_name in (("norm_col_raw", "norm_col"),
                              ("_norm_col_raw", "_norm_col"),
                              ("norm_raw", "norm")):
        raw = getattr(op, raw_name, None)
        fl = getattr(op, fl_name, None)
        if raw is not None and fl is not None:
            to_np = lambda v: np.asarray(v.detach().cpu() if hasattr(v, "detach") else v, dtype=np.float64).ravel()  # noqa: E731
            return to_np(raw), to_np(fl)
    return None


def _score(val, truth, src):
    d = val - truth
    sd = float(np.std(truth)) or 1.0
    lo, hi = float(np.min(src)), float(np.max(src))
    return dict(rmse_over_std=float(np.sqrt(np.mean(d ** 2)) / sd),
                p999_over_std=float(np.percentile(np.abs(d), 99.9) / sd),
                overshoot_frac=float(np.mean((val < lo) | (val > hi))),
                out_min=float(np.min(val)), out_max=float(np.max(val)),
                src_min=lo, src_max=hi,
                bad_frac_5std=float(np.mean(np.abs(d) > 5 * sd)))


def _run(op, src):
    res = op.resample(np.asarray(src, dtype=np.float64))
    return (np.asarray(res.cell_ids, dtype=np.int64).ravel(),
            np.asarray(res.cell_data, dtype=np.float64).ravel())


def _complete(ids, val, cell_ids):
    val, _n = build_tiles.complete_block(ids, val, cell_ids)
    return val


# --- one region ------------------------------------------------------------
def diagnose(scene: str, sweep: bool = True, seed: int = 17, halo_px: int = 32):
    import healpix_resample as hr
    store = cfg.REGION_DIR / f"{scene}_data.zarr"
    w = build_tiles.tile_subwindow(scene, store, halo_px=halo_px)
    cand, to_crs = w["cand"], w["to_crs"]
    sub_lon, sub_lat = w["sub_lon"], w["sub_lat"]
    r0, r1, c0, c1 = w["bounds"]
    sx, sy = to_crs.transform(sub_lon, sub_lat)
    x0, y0 = float(np.mean(cand["cx"])), float(np.mean(cand["cy"]))
    src_x, src_y = np.asarray(sx) - x0, np.asarray(sy) - y0
    cell_x, cell_y = np.asarray(cand["cx"]) - x0, np.asarray(cand["cy"]) - y0
    cell_ids = cand["cell_ids"]
    lines = [f"## {scene}", ""]
    rows = []

    op = build_tiles.make_resampler("bicubic", sub_lon, sub_lat, cell_ids, verbose=False)
    lines += ["### Operator inventory", "", "```"]
    for k, v in _attr_inventory(op).items():
        lines.append(f"{k}: {v}")
    lines += ["```", ""]
    try:
        lines.append(f"constructor signature: `{inspect.signature(hr.BicubicResampler)}`")
    except (TypeError, ValueError):
        pass
    lines.append(f"module: `{getattr(hr, '__file__', '?')}` version "
                 f"`{getattr(hr, '__version__', '?')}`")
    lines.append("")

    # --- probe fields --------------------------------------------------------
    k_min = 2.0 * math.pi / cfg.TILE_SIDE_M
    k_max = math.pi / cfg.DELTA_EQ
    sigma = sup.psf_sigma_from_fwhm(11.9)
    probes = {}
    probes["constant"] = (np.full(src_x.shape, 0.3), np.full(cell_x.shape, 0.3))
    probes["affine"] = (0.2 + 1e-4 * src_x + 5e-5 * src_y,
                        0.2 + 1e-4 * cell_x + 5e-5 * cell_y)
    lam_m = 60.0
    probes["wave_60m"] = (0.2 + 0.05 * np.sin(2 * math.pi * src_x / lam_m),
                          0.2 + 0.05 * np.sin(2 * math.pi * cell_x / lam_m))
    f_src = sup.synthetic_field(src_x, src_y, beta=1.9, k_min=k_min, k_max=k_max,
                                seed=seed, psf_sigma_m=sigma)
    f_truth = sup.synthetic_field(cell_x, cell_y, beta=1.9, k_min=k_min,
                                  k_max=k_max, seed=seed, psf_sigma_m=sigma)
    probes["analytic"] = (f_src, f_truth)

    lines += ["### Probe fields (radial bicubic, default parameters)", "",
              "| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |",
              "|---|---|---|---|---|---|---|---|---|"]
    worst = None
    for name, (src, truth) in probes.items():
        ids, val = _run(op, src)
        val = _complete(ids, val, cell_ids)
        s = _score(val, truth, src)
        rows.append(dict(scene=scene, probe=name, variant="radial_default", **s))
        lines.append(f"| {name} | {s['rmse_over_std']:.3f} | {s['p999_over_std']:.2f} | "
                     f"{100*s['overshoot_frac']:.2f}% | {100*s['bad_frac_5std']:.3f}% | "
                     f"{s['out_min']:.3g} | {s['out_max']:.3g} | {s['src_min']:.3g} | {s['src_max']:.3g} |")
        if name == "analytic":
            worst = (ids, val, truth)
    lines.append("")

    # --- where are the worst cells? ------------------------------------------
    ids, val, truth = worst
    d = np.abs(val - truth)
    sd = float(np.std(truth)) or 1.0
    bad = d > 5 * sd
    lines += ["### Worst cells on the analytic field", ""]
    lines.append(f"- cells with |error| > 5 std: {bad.sum()} of {bad.size} "
                 f"({100*bad.mean():.3f}%); their RMSE alone accounts for "
                 f"{100*np.sum(d[bad]**2)/max(np.sum(d**2),1e-30):.1f}% of the total squared error")
    # distance to the hull edge, in cells, from the row/col box of the tile
    col, row = np.asarray(cand["col"]), np.asarray(cand["row"])
    edge = np.minimum.reduce([col - col.min(), col.max() - col,
                              row - row.min(), row.max() - row])
    if bad.any():
        lines.append(f"- distance to the tile's pixel box (native pixels): bad cells "
                     f"median {np.median(edge[bad]):.1f}, all cells median {np.median(edge):.1f}")
        # spatial clustering: share of bad cells whose 4 NESTED siblings are also bad
        sib = (ids // 4)
        bad_by_parent = pd.Series(bad).groupby(sib).mean()
        lines.append(f"- share of level-19 parents with all four children bad: "
                     f"{100*np.mean(bad_by_parent.values == 1.0):.2f}% "
                     f"(random placement would give ~{100*bad.mean()**4:.4f}%)")
    nrm = _normaliser(op)
    if nrm is not None:
        raw, fl = nrm
        if raw.size == val.size:
            lines.append(f"- raw signed normaliser |norm_col_raw|: bad cells median "
                         f"{np.median(np.abs(raw[bad])):.3e}, others median "
                         f"{np.median(np.abs(raw[~bad])):.3e}")
            floored = np.abs(fl) > np.abs(raw) * (1 + 1e-9)
            lines.append(f"- cells where the floor was applied: {floored.sum()} "
                         f"({100*floored.mean():.3f}%); among bad cells {100*floored[bad].mean():.1f}%")
            lines.append(f"- correlation of log|error| with -log|norm_col_raw|: "
                         f"{np.corrcoef(np.log(d + 1e-12), -np.log(np.abs(raw) + 1e-12))[0,1]:.3f}")
        else:
            lines.append(f"- normaliser has {raw.size} entries, tile has {val.size}: not aligned, skipped")
    else:
        lines.append("- operator exposes no norm_col / norm_col_raw: inspect "
                     "healpix_resample/bicubic.py by hand for the guard value")
    lines.append("")

    # --- parameter sweep -------------------------------------------------------
    if sweep:
        lines += ["### Parameter sweep on the analytic field", "",
                  "| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |", "|---|---|---|---|---|"]
        try:
            sig = inspect.signature(hr.BicubicResampler)
            params = sig.parameters
        except (TypeError, ValueError):
            params = {}
        sig_key = "sigma_m" if "sigma_m" in params else None
        npt_key = next((k for k in ("Npt", "npt", "q", "n_neighbours") if k in params), None)
        common = dict(lon_deg=sub_lon, lat_deg=sub_lat, level=cfg.LEVEL,
                      out_cell_ids=cell_ids, verbose=False)
        for sg in ([None, 10.0, 15.0, 20.0] if sig_key else [None]):
            for npt in ([None, 16, 36] if npt_key else [None]):
                kw = dict(common)
                if sg is not None:
                    kw[sig_key] = sg
                if npt is not None:
                    kw[npt_key] = npt
                try:
                    op2 = hr.BicubicResampler(**kw)
                    ids2, val2 = _run(op2, f_src)
                    val2 = _complete(ids2, val2, cell_ids)
                    s = _score(val2, f_truth, f_src)
                except Exception as e:  # noqa: BLE001
                    lines.append(f"| {sg} | {npt} | error: {str(e)[:60]} | | |")
                    continue
                rows.append(dict(scene=scene, probe="analytic",
                                 variant=f"radial_sigma{sg}_npt{npt}", **s))
                lines.append(f"| {sg or 'default'} | {npt or 'default'} | "
                             f"{s['rmse_over_std']:.3f} | {s['p999_over_std']:.2f} | "
                             f"{100*s['overshoot_frac']:.2f}% |")
        lines.append("")

    # --- classical separable bicubic on the regular raster ---------------------
    try:
        from scipy.ndimage import map_coordinates
        H = r1 - r0
        W = c1 - c0
        grid = np.asarray(f_src, dtype=np.float64).reshape(H, W)
        coords = np.vstack([row - r0, col - c0])
        val3 = map_coordinates(grid, coords, order=3, mode="nearest")
        s = _score(val3, f_truth, f_src)
        rows.append(dict(scene=scene, probe="analytic", variant="separable_scipy_order3", **s))
        lines += ["### Classical separable bicubic on the source raster (scipy order 3)", "",
                  f"RMSE/std {s['rmse_over_std']:.3f}, p99.9/std {s['p999_over_std']:.2f}, "
                  f"overshoot {100*s['overshoot_frac']:.2f}%, >5 std {100*s['bad_frac_5std']:.3f}%", ""]
    except Exception as e:  # noqa: BLE001
        lines += [f"separable bicubic skipped: {e}", ""]
    return lines, rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--scenes", nargs="*", default=None)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--no-sweep", action="store_true")
    args = ap.parse_args(argv)
    if args.all:
        scenes = list(region_sites.load_sites()["scene"])
    elif args.scenes:
        scenes = args.scenes
    else:
        sites = region_sites.load_sites()
        scenes = list(sites.groupby("scene_class").head(1)["scene"])
    report = ["# Bicubic transfer diagnosis", "",
              "Hypothesis under test: near-zero signed weight sums, floored and "
              "then divided, on a small fraction of cells. Facts are in the "
              "tables; the hypothesis is confirmed only if the worst cells are "
              "the floored ones and the separable bicubic shows no such tail.", ""]
    rows = []
    for scene in scenes:
        try:
            lines, r = diagnose(scene, sweep=not args.no_sweep)
        except Exception as e:  # noqa: BLE001
            lines, r = [f"## {scene}", "", f"failed: {e}", ""], []
        report += lines
        rows += r
        print("\n".join(lines))
    out_md = cfg.RESULTS_DIR / "diag_bicubic.md"
    out_md.write_text("\n".join(report), encoding="utf-8")
    if rows:
        pd.DataFrame(rows).to_csv(cfg.RESULTS_DIR / "diag_bicubic.csv", index=False)
    print(f"wrote {out_md}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
