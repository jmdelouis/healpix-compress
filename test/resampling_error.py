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


class _Case(dict):
    """Attribute access over a plain dict, so a case reads like a record."""
    __getattr__ = dict.__getitem__


def synthetic_case(scene: str, store: Path, spectrum, noise_scale=1.0,
                   seed=23, halo_px=32, children=16):
    """Build one region's synthetic scene: latent field, observation, truth.

    Factored out of `measure` so that anything else needing the same field --
    the cross-term measurement of cross_term.py, in particular -- gets it from
    here rather than from a second copy of the construction. Two copies of a
    forward model drift, and a drift between the reconstruction error and the
    compression error measured against it would be invisible in the result.
    """
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

    def cell_field(psf, n_children=children):
        """The estimand is a cell *average*, not the value at the centre.

        The companion study makes this explicit and measures the difference:
        its reconstruction agrees 2.3--2.5 times better with a child-sampled
        cell average than with a point sample at the cell centre. Scoring
        against the centre therefore understates every method, and understates
        a deconvolving one most, since that is the estimator whose output is a
        cell average by construction. The average is taken over the cell's
        children one or two levels down, which is exactly how the companion
        builds it.
        """
        import healpix_geo
        if n_children <= 1:
            return field(cell_x, cell_y, psf)
        depth = int(round(math.log(n_children, 4)))
        kids = (cand["cell_ids"][:, None] * (4 ** depth)
                + np.arange(4 ** depth, dtype=np.int64)[None, :]).ravel()
        klon, klat = healpix_geo.nested.healpix_to_lonlat(
            kids, cfg.LEVEL + depth, ellipsoid="WGS84")
        kx, ky = to_crs.transform(np.asarray(klon, dtype=np.float64),
                                  np.asarray(klat, dtype=np.float64))
        v = field(np.asarray(kx) - x0, np.asarray(ky) - y0, psf)
        return v.reshape(-1, 4 ** depth).mean(axis=1)

    # The latent scene, and what the instrument makes of it. Both are analytic,
    # so each is available at whatever positions a transfer needs.
    src_x, src_y = np.asarray(sx) - x0, np.asarray(sy) - y0
    cell_x, cell_y = np.asarray(cand["cx"]) - x0, np.asarray(cand["cy"]) - y0
    obs_clean = field(src_x, src_y, sigma)
    latent_cells = cell_field(0.0)
    observed_cells = cell_field(sigma)

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

    return _Case(scene=scene, cand=cand, sub_lon=sub_lon, sub_lat=sub_lat,
                 obs=obs, obs_clean=obs_clean, latent_cells=latent_cells,
                 observed_cells=observed_cells, beta=beta, psf_fwhm_m=fwhm,
                 sigma=sigma, gain=gain, noise_std=noise_std,
                 children=int(children), utm_crop=crop)


def transfer_cells(case, kernel: str, lam: float = 1e-3):
    """Run one transfer on the case's observation; return (cells, n_filled)."""
    # `case.sigma` is the Gaussian standard deviation of the response (FWHM /
    # 2 sqrt(2 ln 2)); PSFResampler's `sigma_m` is the operator *scale* s of the
    # companion's Eq. (2), exp(-2 d^2 / s^2), i.e. s = 2 sigma, which is what
    # build_tiles passes for the real tiles through fwhm_to_scale. Passing
    # sigma here handed the inverse a response half as wide as the one that
    # blurred the field (revision 8 diagnosis, diag_psf_matched.py): under the
    # library convention in force before 17 September the two coincided by
    # accident, because psf_sigma_from_fwhm then returned s. The stored
    # results of 16 September must be regenerated with this line in place.
    op = build_tiles.make_resampler(kernel, case.sub_lon, case.sub_lat,
                                    case.cand["cell_ids"], verbose=False,
                                    sigma_m=2.0 * case.sigma)
    # PSFResampler defaults to lam=0, i.e. an undamped solve; the companion
    # study's reported recovery uses lam=1e-3, and running the two at different
    # damping compares different estimators.
    obs = np.asarray(case.obs, dtype=np.float64)
    res = op.resample(obs, lam=lam) if kernel == "psf" else op.resample(obs)
    return build_tiles.complete_block(
        np.asarray(res.cell_ids, dtype=np.int64).ravel(),
        np.asarray(res.cell_data, dtype=np.float64).ravel(),
        case.cand["cell_ids"])


def measure(scene: str, store: Path, kernels, spectrum, noise_scale=1.0,
            seed=23, halo_px=32, lam=1e-3, children=16, verbose=True):
    case = synthetic_case(scene, store, spectrum, noise_scale=noise_scale,
                          seed=seed, halo_px=halo_px, children=children)
    latent_cells, observed_cells = case.latent_cells, case.observed_cells
    obs_clean, beta, fwhm = case.obs_clean, case.beta, case.psf_fwhm_m
    gain, noise_std = case.gain, case.noise_std

    rows = []
    for kernel in kernels:
        estimand = ESTIMAND.get(kernel, "observed")
        truth = latent_cells if estimand == "latent" else observed_cells
        try:
            val, n_filled = transfer_cells(case, kernel, lam=lam)
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
            children=int(children), lam=float(lam),
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


def latex_table(df, control=None):
    """Per-class transfer error for both estimands, plus the noiseless control.

    The recovery fraction is the companion study's own definition,
    1 - RMSE(reconstruction)/RMSE(blurred observation), so the two numbers are
    directly comparable. It is reported here at the measured sensor noise and
    again at zero noise, because the published figure is a noise-free one and
    quoting only one of ours would make the comparison look like agreement or
    like failure depending on which we chose.
    """
    ok = df.dropna(subset=["e_resample"])
    ok = ok.assign(cls=ok.scene.str.split("__").str[1])
    ctrl = None
    if control is not None and len(control):
        ctrl = control.dropna(subset=["e_resample"])
        ctrl = ctrl.assign(cls=ctrl.scene.str.split("__").str[1])
    rows = []
    for cls, g in sorted(ok.groupby("cls")):
        ct = g[g.kernel == "clough_tocher"]
        ps = g[g.kernel == "psf"]
        if not len(ct) or not len(ps):
            continue
        blur = float(ps.blur_error.median())
        # The recovery fraction is formed per region, as the companion's
        # Eq. (24) does, and the class value is the median of those fractions.
        # A ratio of class medians pairs the median error of one region with
        # the median blur of another; over water it read -16% where every
        # region but two recovers (revision 8).
        rec = 100.0 * float((1.0 - ps.e_resample / ps.blur_error).median())
        rec0 = float("nan")
        if ctrl is not None:
            c = ctrl[(ctrl.cls == cls) & (ctrl.kernel == "psf")]
            if len(c):
                rec0 = 100.0 * float((1.0 - c.e_resample / c.blur_error).median())
        rows.append(
            f"{cls.capitalize()} & {float(g.beta.median()):.2f} & "
            f"{float(g.psf_fwhm_m.median()):.1f} & "
            f"{100*float(g.noise_frac_of_std.median()):.1f} & "
            f"{blur:.4f} & {float(ct.e_resample.median()):.4f} & "
            f"{float(ps.e_resample.median()):.4f} & "
            f"{rec:+.1f} & "
            + ("---" if rec0 != rec0 else f"{rec0:+.1f}") + " \\\\")
    n = int(ok.scene.nunique())
    have_control = ctrl is not None and len(ctrl) > 0
    # The estimand's depth is read from the rows, never written by hand: a
    # caption that said "four level-21 children" outlived the run that used
    # sixteen, and a reviewer caught it before we did.
    depths = sorted(set(int(c) for c in ok.children.dropna().unique()))
    if len(depths) != 1:
        raise ValueError(f"one table cannot carry several cell-average depths: "
                         f"{depths}")
    n_ch = depths[0]
    words = {1: "one", 4: "four", 16: "sixteen", 64: "sixty-four"}
    sub_level = cfg.LEVEL + int(round(math.log(n_ch, 4))) if n_ch > 1 else cfg.LEVEL
    estimand = (f"the truth at each cell is the average over its "
                f"{words.get(n_ch, n_ch)} level-{sub_level} children, matching "
                f"the companion study's cell-average estimand"
                if n_ch > 1 else
                "the truth at each cell is a point sample at the cell centre")
    note = (
        "On a synthetic scene built per region from that region's own fitted "
        f"spectral slope $\\beta$ and response width, over {n} regions. The "
        "latent field is analytic, the observation is that field convolved "
        "exactly with the fitted Gaussian response, sampled at the native "
        "pixel positions and given additive noise at the level Immerkaer's "
        f"estimator reads off the region's own imagery; {estimand}. "
        "\\emph{Blur} is the RMSE between the observed and latent fields, i.e.\\ "
        "what a non-deconvolving transfer leaves on the table by construction. "
        "\\texttt{clough\\_tocher} is scored against the observed field it "
        "claims to represent, \\texttt{psf} against the latent field it claims "
        "to recover; the last two columns are the median over regions of the "
        "companion study's per-region recovery fraction "
        "$1-e_{\\mathrm{psf}}/e_{\\mathrm{blur}}$, at the measured noise "
        "and in a noiseless control"
        + ("" if have_control else
           ", the latter left blank because no control at the same "
           "cell-average depth is available")
        + ".")
    return (
        "\\begin{table*}[t]\n\\centering\n"
        + "\n".join(cfg.tex_caption(
            "What each transfer costs against the reference it answers to",
            "tab:resampling_error")) + "\n\\footnotesize\n"
        "\\begin{tabular}{@{}lrrr rrr rr@{}}\n\\toprule\n"
        " & \\multicolumn{3}{c}{Scene} & \\multicolumn{3}{c}{RMSE} & "
        "\\multicolumn{2}{c}{\\texttt{psf} recovery (\\%)} \\\\\n"
        "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-9}\n"
        "Class & $\\beta$ & FWHM (m) & Noise (\\%) & Blur & "
        "\\texttt{clough\\_tocher} & \\texttt{psf} & measured & noiseless "
        "\\\\\n\\midrule\n"
        + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n"
        + "\n".join(cfg.tex_note(note, wide=True)) + "\n\\end{table*}\n")


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
    ap.add_argument("--lam", type=float, default=1e-3,
                    help="Tikhonov damping of the PSF resampler's CG solve; "
                         "the companion study reports its recovery at 1e-3, "
                         "and the library defaults to 0")
    ap.add_argument("--children", type=int, default=16, choices=[1, 4, 16],
                    help="samples per cell used to build the cell-average "
                         "estimand; 1 reverts to a point sample at the cell "
                         "centre, which understates every method")
    ap.add_argument("--seed", type=int, default=23)
    ap.add_argument("--keep-mismatched-windows", action="store_true",
                    help="keep the regions window_radiometry.py rejects; this "
                         "measurement re-runs the transfer from the native "
                         "window, so on those regions it prices a difference "
                         "between two products, not a transfer error")
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
    if not args.keep_mismatched_windows:
        drop = set(cfg.window_mismatched_scenes(verbose=True))
        excluded = [s for s in scenes if s in drop]
        if excluded:
            print(f"window audit excludes {len(excluded)} region(s): "
                  + ", ".join(excluded))
            scenes = [s for s in scenes if s not in drop]
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
                                noise_scale=args.noise, seed=args.seed,
                                lam=args.lam, children=args.children))
        except Exception as exc:                           # noqa: BLE001
            print(f"    FAILED {type(exc).__name__}: {exc}")
    if not rows:
        print("no region produced a row")
        return 1

    df = pd.DataFrame(rows)
    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    # The noiseless run is a control, not the measurement, and it used to write
    # the same two filenames -- so running it second silently replaced the
    # result it was meant to be compared against.
    tag = "" if args.noise == 1.0 else f"_noise{args.noise:g}"
    tag += "" if args.children == 16 else f"_c{args.children}"
    csv = cfg.RESULTS_DIR / f"resampling_error{tag}.csv"
    df.to_csv(csv, index=False)
    md = summarize(df)
    (cfg.RESULTS_DIR / f"resampling_error{tag}.md").write_text(md,
                                                               encoding="utf-8")
    print("\n" + md)
    print(f"wrote {csv} ({len(df)} rows)")

    # Only the measurement publishes a table; the control feeds its last two
    # columns, so a control-only run cannot overwrite it with half a result.
    if not tag:
        c = cfg.RESULTS_DIR / f"resampling_error_noise0{tag}.csv"
        control = pd.read_csv(c) if c.exists() else None
        if control is None:
            print(f"note: no {c.name}; the noiseless recovery column will be "
                  f"blank. Run --noise 0 --children {args.children} as well.")
        elif "children" not in control:
            # A control written before the depth was recorded cannot be shown
            # to match, and an unverifiable control is worse than none: it
            # would be published in the same table as a matched one.
            print(f"note: {c.name} carries no cell-average depth, so it "
                  f"predates this script's provenance columns and cannot be "
                  f"matched to the measurement; dropping it. Re-run "
                  f"--noise 0 --children {args.children} with this version.")
            control = None
        elif set(control.children) != set(df.children):
            print(f"note: {c.name} was measured at a different cell-average "
                  f"depth ({sorted(set(control.children))} against "
                  f"{sorted(set(df.children))}); dropping it rather than "
                  f"publishing two depths in one table.")
            control = None
        cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
        out = cfg.TEX_TABLES / "table_resampling_error.tex"
        out.write_text(latex_table(df, control), encoding="utf-8")
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
