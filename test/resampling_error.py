"""Reconstruction error of each transfer, against its own estimand (V3, step 1).

A compression target has no meaning on its own. It is only interpretable
against the error the product already carries from being transferred onto the
carrier grid, and that error is not the same quantity for the two transfers
this study cares about:

    clough_tocher  claims to represent the *observed* field -- what the
                   instrument, with its response, actually measured. Its
                   estimand is the PSF-convolved scene, and it neither adds nor
                   removes spatial structure by design.

    psf            claims to represent the *latent* field -- the scene before
                   the instrument blurred it. It deconvolves, so it restores
                   small-scale structure by construction, and its estimand is
                   the unconvolved scene.

Scoring both against the same reference, or comparing their compression at the
same absolute RMSE, compares errors on two different quantities. This module
measures each against its own estimand, which is only possible synthetically:
the latent field is analytic, the observation is that field convolved with the
instrument response, sampled at the native pixel centres, plus sensor noise.

The noise term is not optional. A deconvolving transfer amplifies noise exactly
as it amplifies structure, and it does so at the high frequencies a codec pays
most for, so a noiseless synthetic would flatter the PSF transfer twice over --
once on reconstruction error and again on compressibility. The level is
estimated from the imagery itself with Immerkaer's high-pass estimator, whose
mask annihilates local quadratics and therefore responds to noise rather than
to scene structure.

Usage
-----
    python resampling_error.py --limit 8          # ~15 min
    python resampling_error.py                    # all 40 regions
    python resampling_error.py --noise 0          # noiseless control

Writes results/resampling_error.csv and results/resampling_error.md.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

import paper_config as cfg
import build_tiles
import region_sites
import support_controls as sup

# Which reference each transfer is answerable to.
ESTIMAND = {"clough_tocher": "observed", "linear": "observed",
            "nearest": "observed", "bicubic": "observed", "psf": "latent"}


def estimate_noise(img) -> float:
    """Per-pixel noise standard deviation, by Immerkaer's estimator.

    The 3x3 mask [[1,-2,1],[-2,4,-2],[1,-2,1]] annihilates any local quadratic,
    so its response over a smooth scene is dominated by noise rather than by
    structure. The constant makes the result an unbiased standard deviation for
    white Gaussian noise. Preferred here over fitting the flat tail of the
    radial spectrum, which is entangled with the response width the same fit
    has to estimate.
    """
    a = np.asarray(img, dtype=np.float64)
    a = a[np.isfinite(a).all(axis=1)] if a.ndim == 2 else a
    if a.ndim != 2 or min(a.shape) < 5:
        return np.nan
    m = np.array([[1.0, -2.0, 1.0], [-2.0, 4.0, -2.0], [1.0, -2.0, 1.0]])
    h, w = a.shape
    acc = np.zeros((h - 2, w - 2), dtype=np.float64)
    for i in range(3):
        for j in range(3):
            acc += m[i, j] * a[i:h - 2 + i, j:w - 2 + j]
    return float(np.sum(np.abs(acc)) * math.sqrt(0.5 * math.pi)
                 / (6.0 * (w - 2) * (h - 2)))


def scene_parameters(scene: str, spectrum: pd.DataFrame | None):
    """Per-region spectral slope and response width, fitted or median."""
    beta, fwhm = 1.9, 11.9
    if spectrum is not None and len(spectrum):
        row = spectrum[spectrum.scene == scene]
        if len(row):
            b = float(row.beta_fit.iloc[0])
            f = float(row.psf_fwhm_m.iloc[0])
            # A fit that lands outside the physically sensible range is a
            # noise-dominated crop, not a measurement: fall back rather than
            # propagate it (Lake Victoria returns beta = -0.19).
            if np.isfinite(b) and 0.5 <= b <= 4.0:
                beta = b
            if np.isfinite(f) and 6.0 <= f <= 20.0:
                fwhm = f
    return beta, fwhm


def measure(scene: str, store: Path, kernels, spectrum, noise_scale=1.0,
            seed=23, halo_px=32, verbose=True):
    # Checked once, with a diagnostic: a window that is not the one the tile
    # was built from produces a partial sub-window and a downstream "43% of
    # cells missing" that points at the halo instead of at the geometry.
    w = build_tiles.tile_subwindow(scene, store, halo_px=halo_px)
    cand, to_crs = w["cand"], w["to_crs"]
    sub_lon, sub_lat = w["sub_lon"], w["sub_lat"]
    _i, _c, utm_crop, _m, npz = build_tiles.load_tile(scene)
    sx, sy = to_crs.transform(sub_lon, sub_lat)
    x0, y0 = float(np.mean(cand["cx"])), float(np.mean(cand["cy"]))

    beta, fwhm = scene_parameters(scene, spectrum)
    sigma = sup.psf_sigma_from_fwhm(fwhm)
    k_min = 2.0 * math.pi / cfg.TILE_SIDE_M
    k_max = math.pi / cfg.DELTA_EQ

    def field(px, py, psf):
        return sup.synthetic_field(px, py, beta=beta, k_min=k_min, k_max=k_max,
                                   seed=seed, psf_sigma_m=psf)

    # The latent scene, and what the instrument makes of it. Both are analytic,
    # so each is available at whatever positions a transfer needs.
    src_x, src_y = np.asarray(sx) - x0, np.asarray(sy) - y0
    cell_x, cell_y = np.asarray(cand["cx"]) - x0, np.asarray(cand["cy"]) - y0
    obs_clean = field(src_x, src_y, sigma)
    latent_cells = field(cell_x, cell_y, 0.0)
    observed_cells = field(cell_x, cell_y, sigma)

    # Put the synthetic on the same scale as the region so the noise estimated
    # from the imagery is the right size relative to it.
    crop = np.asarray(utm_crop, dtype=np.float64)
    sd_obs = float(np.std(obs_clean)) or 1.0
    gain = float(np.std(crop[np.isfinite(crop)])) / sd_obs
    obs_clean = obs_clean * gain
    latent_cells = latent_cells * gain
    observed_cells = observed_cells * gain

    noise_std = estimate_noise(crop) * float(noise_scale)
    rng = np.random.default_rng(seed + 1)
    obs = obs_clean + (rng.normal(0.0, noise_std, obs_clean.shape)
                       if noise_std > 0 else 0.0)

    rows = []
    for kernel in kernels:
        estimand = ESTIMAND.get(kernel, "observed")
        truth = latent_cells if estimand == "latent" else observed_cells
        try:
            op = build_tiles.make_resampler(kernel, sub_lon, sub_lat,
                                            cand["cell_ids"], verbose=False)
            res = op.resample(np.asarray(obs, dtype=np.float64))
            val, n_filled = build_tiles.complete_block(
                np.asarray(res.cell_ids, dtype=np.int64).ravel(),
                np.asarray(res.cell_data, dtype=np.float64).ravel(),
                cand["cell_ids"])
        except Exception as exc:                           # noqa: BLE001
            print(f"    {kernel:14s} FAILED {type(exc).__name__}: {exc}")
            rows.append(dict(scene=scene, kernel=kernel, estimand=estimand,
                             e_resample=np.nan,
                             fail=f"{type(exc).__name__}: {exc}"[:200]))
            continue
        d = val - truth
        sd_truth = float(np.std(truth)) or 1.0
        # What the transfer would score against the other estimand, so the
        # asymmetry between the two products is visible rather than asserted.
        other = observed_cells if estimand == "latent" else latent_cells
        rows.append(dict(
            scene=scene, kernel=kernel, estimand=estimand,
            beta=beta, psf_fwhm_m=fwhm, gain=gain,
            noise_std=noise_std,
            noise_frac_of_std=noise_std / (float(np.std(obs_clean)) or 1.0),
            e_resample=float(np.sqrt(np.mean(d ** 2))),
            e_resample_frac_of_std=float(np.sqrt(np.mean(d ** 2)) / sd_truth),
            e_against_other_estimand=float(
                np.sqrt(np.mean((val - other) ** 2))),
            blur_error=float(np.sqrt(np.mean(
                (observed_cells - latent_cells) ** 2))),
            truth_std=sd_truth, cells_filled=int(n_filled)))
        if verbose:
            r = rows[-1]
            print(f"    {kernel:14s} [{estimand:8s}] e_resample "
                  f"{r['e_resample']:.4e} "
                  f"({100*r['e_resample_frac_of_std']:.1f}% of std) | "
                  f"vs other estimand {r['e_against_other_estimand']:.4e} | "
                  f"blur {r['blur_error']:.4e}")
    return rows


def summarize(df):
    out = ["# Reconstruction error per transfer, against its own estimand", ""]
    ok = df.dropna(subset=["e_resample"])
    if not len(ok):
        return "\n".join(out + ["_no measurement succeeded_"])
    out.append(f"- regions: {ok.scene.nunique()}, "
               f"noise: median {ok.noise_std.median():.4e} "
               f"({100*ok.noise_frac_of_std.median():.1f}% of the field's std)")
    out.append("")
    ok = ok.assign(cls=ok.scene.str.split("__").str[1])
    out += ["## By transfer and scene class", "",
            "| kernel | estimand | class | n | e_resample | % of estimand std | "
            "blur error |", "|---|---|---|---|---|---|---|"]
    for (kernel, cls), g in ok.groupby(["kernel", "cls"]):
        out.append(f"| `{kernel}` | {g.estimand.iloc[0]} | {cls} | {len(g)} | "
                   f"{g.e_resample.median():.4e} | "
                   f"{100*g.e_resample_frac_of_std.median():.1f}% | "
                   f"{g.blur_error.median():.4e} |")
    out.append("")
    out += ["## Medians over all regions", "",
            "| kernel | estimand | e_resample | % of estimand std | against the "
            "other estimand |", "|---|---|---|---|---|"]
    for kernel, g in ok.groupby("kernel"):
        out.append(f"| `{kernel}` | {g.estimand.iloc[0]} | "
                   f"{g.e_resample.median():.4e} | "
                   f"{100*g.e_resample_frac_of_std.median():.1f}% | "
                   f"{g.e_against_other_estimand.median():.4e} |")
    out.append("")
    out.append("The last column is what each transfer would score if it were "
               "held to the other product's reference. A large gap between the "
               "two columns is the whole reason these products cannot be "
               "compared at the same absolute RMSE.")
    out.append("")
    ct = ok[ok.kernel == "clough_tocher"]
    ps = ok[ok.kernel == "psf"]
    if len(ct) and len(ps):
        out.append(f"- `clough_tocher` against the observed field: "
                   f"{ct.e_resample.median():.4e}")
        out.append(f"- `psf` against the latent field: "
                   f"{ps.e_resample.median():.4e}")
        out.append(f"- blur alone (observed minus latent): "
                   f"{ok.blur_error.median():.4e}; `psf` recovers "
                   f"{100*(1 - ps.e_resample.median()/ok.blur_error.median()):.1f}% "
                   f"of it")
        out.append("")
        out.append("These are the two numbers a compression target has to be "
                   "read against, one per product. They are written to "
                   "results/resampling_error.csv for `alpha_budget.py`.")
    return "\n".join(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kernels", nargs="*", default=["clough_tocher", "psf"])
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--noise", type=float, default=1.0,
                    help="multiplier on the noise estimated from the imagery; "
                         "0 gives the noiseless control, which flatters a "
                         "deconvolving transfer and is reported as a control "
                         "only")
    ap.add_argument("--seed", type=int, default=23)
    args = ap.parse_args(argv)

    spec_path = cfg.RESULTS_DIR / "support_controls_spectrum.csv"
    spectrum = pd.read_csv(spec_path) if spec_path.exists() else None
    if spectrum is None:
        print("note: no support_controls_spectrum.csv; using the median slope "
              "and response width for every region")

    scenes = args.regions or list(region_sites.load_sites()["scene"])
    scenes = [s for s in scenes
              if (cfg.REGION_DIR / f"{s}_data.zarr").exists()
              and build_tiles.tile_path(s).exists()]
    if not scenes:
        print(f"no region with both a window under {cfg.REGION_DIR} and a tile "
              f"under {cfg.TILE_DIR}")
        return 1
    if args.limit:
        scenes = scenes[:args.limit]

    print(f"{len(scenes)} region(s), kernels {args.kernels}, "
          f"noise x{args.noise:g}")
    rows = []
    for scene in scenes:
        print(f"  {scene}")
        store = cfg.REGION_DIR / f"{scene}_data.zarr"
        try:
            rows.extend(measure(scene, store, args.kernels, spectrum,
                                noise_scale=args.noise, seed=args.seed))
        except Exception as exc:                           # noqa: BLE001
            print(f"    FAILED {type(exc).__name__}: {exc}")
    if not rows:
        print("no region produced a row")
        return 1

    df = pd.DataFrame(rows)
    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    csv = cfg.RESULTS_DIR / "resampling_error.csv"
    df.to_csv(csv, index=False)
    md = summarize(df)
    (cfg.RESULTS_DIR / "resampling_error.md").write_text(md, encoding="utf-8")
    print("\n" + md)
    print(f"wrote {csv} ({len(df)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
