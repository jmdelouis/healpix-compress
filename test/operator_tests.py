"""Characterise the HEALPix -> UTM reconstruction operator B (V5).

The common-reference evaluation of the paper rests on one operator: the map
that reads a stored HEALPix field on the native UTM pixels, so that both grids
can be scored against the same observation. Until now the only thing verified
about it was that it preserves a constant field. That check is necessary --- an
operator that loses a constant is worthless, and the failure is invisible in
the results because every rate-distortion point then reports the same "error"
--- but it is very far from sufficient: a badly registered operator, a
smoothing one, and an aliasing one all pass it.

This module measures what the operator actually does, on fields whose exact
value at every position is known analytically. That is the point of using
analytic fields rather than imagery: the ground truth on the UTM side is the
same function evaluated at the UTM pixel centres, so nothing is resampled
anywhere and the residual is the operator's own error and the unavoidable
discretisation, and nothing else.

The battery, in order of what it isolates:

    constant     the existing sanity check, kept so its value is reported
                 alongside the rest rather than in isolation;
    ramp         linear in x and in y: an operator reading the containing
                 cell's centre must err by about half a cell times the
                 gradient, so this checks that the error is the geometric one
                 and not a registration bug;
    sinusoid     the main measurement. At each wavenumber and orientation the
                 projected field is regressed on the same wave to recover an
                 amplitude response and a phase, giving a transfer function and
                 a registration error in metres, per direction. HEALPix cells
                 are equal-area but not equal-shape, so the response is
                 expected to be anisotropic and that anisotropy is measurable
                 here rather than argued;
    spectrum     one band-limited field with the slope and instrument width
                 fitted to the imagery (Section on support controls), which
                 gives a single number directly comparable with the
                 representation floor the paper reports: how much of that floor
                 is the operator, and how much is the carrier grid itself;
    edge         a step across the tile at four orientations, where a
                 nearest-cell read is at its worst;
    impulse      a single cell set to one, to see where the operator puts it.

Usage
-----
    python operator_tests.py --limit 3
    python operator_tests.py --modes nearest knn      # if the weighted
                                                      # operator is usable here
Writes results/operator_tests.csv, results/operator_tests.md,
figures/fig_operator_response.pdf and ../tex/tables/table_operator_tests.tex.
"""

from __future__ import annotations

import argparse
import math
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt           # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402

import paper_config as cfg                # noqa: E402
import build_tiles                        # noqa: E402
import common_reference                   # noqa: E402
import region_sites                       # noqa: E402
import support_controls as sup            # noqa: E402

# Orientations in the projected frame. 0 and 90 degrees are the UTM raster
# axes; 45 and 135 are where the HEALPix raster axes sit for a tile near the
# equator, which is why both pairs are needed to see anisotropy.
ORIENTATIONS_DEG = [0.0, 45.0, 90.0, 135.0]


def geometry(scene: str):
    """Cell centres and UTM pixel centres of one region, in one common frame."""
    import pyproj
    reg = sup.load_region(scene)
    crs = pyproj.CRS.from_wkt(reg["crs_wkt"])
    to_crs = pyproj.Transformer.from_crs(pyproj.CRS.from_epsg(4326), crs,
                                         always_xy=True)
    lon, lat = sup.cell_lonlat(reg["cell_ids"])
    cx, cy = to_crs.transform(lon, lat)
    x0, y0 = float(np.mean(cx)), float(np.mean(cy))
    gx, gy = np.meshgrid(reg["utm_x"] - x0, reg["utm_y"] - y0)
    return reg, (np.asarray(cx) - x0, np.asarray(cy) - y0), (gx, gy)


def wave(px, py, k, theta, phase=0.0):
    return np.cos(k * (np.cos(theta) * px + np.sin(theta) * py) + phase)


def fit_response(projected, px, py, k, theta, mask):
    """Amplitude response and phase shift of a projected single-frequency wave.

    Regressing the projected field on the same wave in quadrature recovers
    both at once: the amplitude is the operator's transfer at this wavenumber
    and direction, and the phase converts to a displacement in metres, which is
    the quantity a registration error actually is.
    """
    c = wave(px, py, k, theta)[mask]
    s = np.sin(k * (np.cos(theta) * px + np.sin(theta) * py))[mask]
    p = np.asarray(projected)[mask]
    ok = np.isfinite(p) & np.isfinite(c) & np.isfinite(s)
    if ok.sum() < 16:
        return np.nan, np.nan, np.nan
    A = np.column_stack([c[ok], s[ok]])
    coef, *_ = np.linalg.lstsq(A, p[ok], rcond=None)
    amp = float(np.hypot(coef[0], coef[1]))
    # a cos(u) + b sin(u) = R cos(u - phi) with phi = atan2(b, a), while a field
    # displaced by d along the wave vector is cos(u + k d). Matching the two
    # gives d = -phi/k: the sign matters, because it is the difference between
    # saying the operator reads the field ahead of where it is and behind it.
    phi = float(np.arctan2(coef[1], coef[0]))
    shift_m = float(-phi / k)   # positive: the operator places the field
                                # displaced by +shift along the wave vector
    resid = p[ok] - A @ coef
    rms_resid = float(np.sqrt(np.mean(resid ** 2)))
    return amp, shift_m, rms_resid


def run_region(scene: str, mode: str, k_points: int = 10, verbose=True):
    reg, (cx, cy), (gx, gy) = geometry(scene)
    cr = common_reference.CommonReference(reg["cell_ids"], reg["utm_lon"],
                                          reg["utm_lat"], mode=mode)
    mask = cr.mask
    rows = []

    def record(**kw):
        rows.append(dict(scene=scene, mode=cr.mode, **kw))

    def project(cell_values):
        return np.asarray(cr.predict(np.asarray(cell_values, dtype=np.float64)))

    def err(pred, truth):
        d = (np.asarray(pred) - np.asarray(truth))[mask]
        d = d[np.isfinite(d)]
        return (float(np.sqrt(np.mean(d ** 2))) if d.size else np.nan,
                float(np.max(np.abs(d))) if d.size else np.nan)

    # --- constant -----------------------------------------------------------
    r, m = err(project(np.ones(cfg.N_CHILDREN)), np.ones_like(gx))
    record(family="constant", parameter="", rmse=r, max_abs=m,
           amplitude=np.nan, shift_m=np.nan)

    # --- ramps --------------------------------------------------------------
    # Slope of one unit per tile side, so the error is directly readable as a
    # fraction of a cell: a nearest-cell read should give ~ half a cell.
    for name, (pc, pg) in (("ramp_x", (cx, gx)), ("ramp_y", (cy, gy))):
        slope = 1.0 / cfg.TILE_SIDE_M
        r, m = err(project(pc * slope), pg * slope)
        record(family="ramp", parameter=name, rmse=r, max_abs=m,
               rmse_cells=r / (slope * cfg.DELTA_EQ),
               max_cells=m / (slope * cfg.DELTA_EQ),
               amplitude=np.nan, shift_m=np.nan)

    # --- sinusoids: the transfer function -----------------------------------
    k_min = 2.0 * math.pi / cfg.TILE_SIDE_M
    k_max = math.pi / cfg.DELTA_EQ
    for k in np.geomspace(k_min * 4, k_max, k_points):
        for deg in ORIENTATIONS_DEG:
            th = math.radians(deg)
            proj = project(wave(cx, cy, k, th))
            truth = wave(gx, gy, k, th)
            r, m = err(proj, truth)
            amp, shift, resid = fit_response(proj, gx, gy, k, th, mask)
            record(family="sinusoid", parameter=f"{deg:g}deg",
                   k_rad_per_m=float(k), wavelength_m=float(2 * math.pi / k),
                   orientation_deg=deg, rmse=r, max_abs=m,
                   amplitude=amp, shift_m=shift, fit_residual_rms=resid)

    # --- one band-limited field with the imagery's own parameters -----------
    sigma = sup.psf_sigma_from_fwhm(11.9)
    f_cells = sup.synthetic_field(cx, cy, beta=1.9, k_min=k_min, k_max=k_max,
                                  seed=11, psf_sigma_m=sigma)
    f_pix = sup.synthetic_field(gx, gy, beta=1.9, k_min=k_min, k_max=k_max,
                                seed=11, psf_sigma_m=sigma)
    sd = float(np.std(f_pix)) or 1.0
    f_cells = (f_cells - np.mean(f_pix)) / sd
    f_pix = (f_pix - np.mean(f_pix)) / sd
    r, m = err(project(f_cells), f_pix)
    record(family="spectrum", parameter="beta=1.9, fwhm=11.9m",
           rmse=r, max_abs=m, amplitude=np.nan, shift_m=np.nan,
           note="unit-variance field, so rmse is a fraction of the field's std")

    # --- step edges ---------------------------------------------------------
    for deg in ORIENTATIONS_DEG:
        th = math.radians(deg)
        proj_c = np.cos(th) * cx + np.sin(th) * cy
        proj_g = np.cos(th) * gx + np.sin(th) * gy
        r, m = err(project((proj_c > 0).astype(float)),
                   (proj_g > 0).astype(float))
        record(family="edge", parameter=f"{deg:g}deg", orientation_deg=deg,
               rmse=r, max_abs=m, amplitude=np.nan, shift_m=np.nan)

    # --- impulse ------------------------------------------------------------
    # A single cell at the tile centre: where does the operator put it, and
    # over how many native pixels does it spread?
    imp = np.zeros(cfg.N_CHILDREN)
    centre = int(np.argmin(cx ** 2 + cy ** 2))
    imp[centre] = 1.0
    p = project(imp)
    hit = np.asarray(p)[mask] > 0.5
    n_hit = int(hit.sum())
    if n_hit:
        sel = np.zeros_like(np.asarray(p), dtype=bool)
        sel[mask] = hit
        dx = float(np.mean(gx[sel])) - float(cx[centre])
        dy = float(np.mean(gy[sel])) - float(cy[centre])
        offset = math.hypot(dx, dy)
    else:
        offset = np.nan
    record(family="impulse", parameter="one cell", rmse=np.nan, max_abs=np.nan,
           amplitude=np.nan, shift_m=np.nan, n_pixels_hit=n_hit,
           centroid_offset_m=offset,
           note=f"one level-20 cell covers about "
                f"{(cfg.DELTA_EQ / cfg.NATIVE_PIXEL_M) ** 2:.2f} native pixels")

    if verbose:
        sin_ = [r for r in rows if r["family"] == "sinusoid"]
        hi = [r for r in sin_ if r["wavelength_m"] < 3 * cfg.NATIVE_PIXEL_M]
        print(f"    constant {rows[0]['rmse']:.2e} | "
              f"ramp {rows[1]['rmse_cells']:.2f} cells | "
              f"spectrum {[r for r in rows if r['family']=='spectrum'][0]['rmse']:.3e} | "
              f"|shift| max {np.nanmax([abs(r['shift_m']) for r in sin_]):.2f} m"
              + (f" | response at 3 px {np.nanmedian([r['amplitude'] for r in hi]):.2f}"
                 if hi else ""))
    return rows


# --- reporting ---------------------------------------------------------------
def figure(df, publish=True, name="fig_operator_response"):
    d = df[df.family == "sinusoid"].dropna(subset=["amplitude"])
    if not len(d):
        print("  no sinusoid rows; figure skipped")
        return
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.0, 5.4), sharex=True,
                                   gridspec_kw=dict(hspace=0.12))
    cmap = plt.get_cmap("cividis")
    degs = sorted(d.orientation_deg.unique())
    for i, deg in enumerate(degs):
        g = d[d.orientation_deg == deg].groupby("wavelength_m")
        lam = np.array(sorted(g.groups))
        amp = g.amplitude.median().reindex(lam).to_numpy()
        sh = g.shift_m.median().reindex(lam).to_numpy()
        c = cmap(i / max(len(degs) - 1, 1))
        ax1.plot(lam, amp, "o-", ms=3, lw=1.4, color=c, label=f"{deg:g}$^\\circ$")
        ax2.plot(lam, sh, "o-", ms=3, lw=1.4, color=c)
    for ax in (ax1, ax2):
        ax.set_xscale("log")
        ax.axvline(2 * cfg.NATIVE_PIXEL_M, color="0.6", ls=":", lw=1.0)
        ax.axvline(2 * cfg.DELTA_EQ, color="0.6", ls="--", lw=1.0)
    ax1.axhline(1.0, color="0.75", lw=0.8, zorder=0)
    ax1.set_ylabel("amplitude response")
    ax1.legend(fontsize=8, frameon=False, ncol=4, title="wave direction",
               title_fontsize=8, loc="lower right")
    ax1.set_title("HEALPix $\\to$ UTM operator: transfer and registration",
                  fontsize=10)
    ax2.axhline(0.0, color="0.75", lw=0.8, zorder=0)
    ax2.set_ylabel("registration shift (m)")
    ax2.set_xlabel("wavelength (m); dotted: UTM Nyquist, dashed: HEALPix Nyquist")
    cfg.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(cfg.FIGURES_DIR / f"{name}.{ext}", bbox_inches="tight", dpi=160)
    if publish:
        cfg.TEX_FIGURES.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(cfg.FIGURES_DIR / f"{name}.pdf",
                        cfg.TEX_FIGURES / f"{name}.pdf")
    plt.close(fig)
    print(f"  {name}")


def _med(df, family, col, **eq):
    d = df[df.family == family]
    for k, v in eq.items():
        d = d[d[k] == v]
    v = d[col].to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    return float(np.median(v)) if v.size else np.nan


def summarize(df, floor=None):
    out = ["# Reconstruction operator B: what it does beyond a constant (V5)", ""]
    out.append(f"- regions: {df.scene.nunique()}, "
               f"modes: {', '.join(sorted(df['mode'].unique()))}")
    out.append("")
    out.append("| test | quantity | median |")
    out.append("|---|---|---|")
    out.append(f"| constant | RMSE (should be 0) | "
               f"{_med(df,'constant','rmse'):.2e} |")
    for p in ("ramp_x", "ramp_y"):
        out.append(f"| {p} | RMSE, in cells of gradient | "
                   f"{_med(df,'ramp','rmse_cells',parameter=p):.3f} |")
    sp = _med(df, "spectrum", "rmse")
    out.append(f"| band-limited field | RMSE, fraction of the field's std | "
               f"{sp:.3e} |")
    for deg in ORIENTATIONS_DEG:
        out.append(f"| edge {deg:g}deg | RMSE across a step | "
                   f"{_med(df,'edge','rmse',parameter=f'{deg:g}deg'):.3f} |")
    out.append(f"| impulse | native pixels lit by one cell | "
               f"{_med(df,'impulse','n_pixels_hit'):.1f} |")
    out.append(f"| impulse | centroid offset (m) | "
               f"{_med(df,'impulse','centroid_offset_m'):.2f} |")
    out.append("")

    d = df[df.family == "sinusoid"].dropna(subset=["amplitude"])
    if len(d):
        out += ["## Transfer function", "",
                "| wavelength (m) | response 0deg | 45deg | 90deg | 135deg | "
                "max &#124;shift&#124; (m) |", "|---|---|---|---|---|---|"]
        for lam, g in d.groupby("wavelength_m"):
            cells = [f"{g[g.orientation_deg==deg].amplitude.median():.3f}"
                     if (g.orientation_deg == deg).any() else "--"
                     for deg in ORIENTATIONS_DEG]
            out.append(f"| {lam:.1f} | " + " | ".join(cells)
                       + f" | {np.nanmax(np.abs(g.shift_m)):.2f} |")
        out.append("")
        aniso = []
        for lam, g in d.groupby("wavelength_m"):
            a = g.groupby("orientation_deg").amplitude.median()
            if len(a) > 1 and np.isfinite(a).all() and a.max() > 0:
                aniso.append((lam, float(a.max() / a.min())))
        if aniso:
            lam, worst = max(aniso, key=lambda t: t[1])
            out.append(f"- strongest directional anisotropy: "
                       f"**{worst:.2f}x** between directions at "
                       f"{lam:.1f} m wavelength")
        out.append(f"- largest registration shift over all wavelengths and "
                   f"directions: **{np.nanmax(np.abs(d.shift_m)):.2f} m** "
                   f"({np.nanmax(np.abs(d.shift_m))/cfg.DELTA_EQ:.2f} cells)")
        out.append("")
    if floor is not None and np.isfinite(sp):
        out += ["## Against the representation floor", "",
                f"The paper reports a representation floor of {floor:.3e} in "
                f"median, about 10% of a scene's standard deviation. On a field "
                f"of the imagery's own spectrum this operator contributes "
                f"{sp:.3e} of the field's standard deviation, i.e. "
                f"{100*sp:.2f}% against the floor's ~10%. The floor is "
                + ("therefore dominated by the carrier grid's own transfer "
                   "step and not by the reading operator."
                   if sp < 0.03 else
                   "NOT clearly dominated by the carrier grid: the operator "
                   "accounts for a comparable share, and Section V-E's "
                   "conclusion would need restating."), ""]
    return "\n".join(out)


def latex_table(df):
    d = df[df.family == "sinusoid"].dropna(subset=["amplitude"])
    if not len(d):
        return ""
    rows = []
    for lam, g in d.groupby("wavelength_m"):
        cells = [(f"{g[g.orientation_deg==deg].amplitude.median():.3f}"
                  if (g.orientation_deg == deg).any() else "--")
                 for deg in ORIENTATIONS_DEG]
        rows.append(f"{lam:.0f} & " + " & ".join(cells)
                    + f" & {np.nanmax(np.abs(g.shift_m)):.2f} \\\\")
    sp = _med(df, "spectrum", "rmse")
    imp = _med(df, "impulse", "n_pixels_hit")
    ramp = _med(df, "ramp", "rmse_cells", parameter="ramp_x")
    return (
        "\\begin{table}[t]\n\\centering\n"
        "\\caption{What the HEALPix$\\to$UTM reading operator does beyond "
        "preserving a constant. Each row is a single-frequency wave placed "
        "analytically on the level-20 cell centres, read on the native UTM "
        "pixels, and regressed against the same wave evaluated at those pixels; "
        "the columns are the recovered amplitude response by wave direction in "
        "the projected frame, and the largest phase shift converted to metres. "
        "Nothing is resampled: the ground truth is the same analytic function, "
        "so the residual is the operator and the unavoidable discretisation of "
        "reading a cell centre. On a band-limited field with the imagery's own "
        f"fitted spectrum the operator contributes {sp:.1e} of the field's "
        f"standard deviation; one cell lights {imp:.1f} native pixels; a unit "
        f"gradient is reproduced to {ramp:.2f} cells.}}\n"
        "\\label{tab:operator_tests}\n\\footnotesize\n"
        "\\begin{tabular}{@{}rccccr@{}}\n\\toprule\n"
        "Wavelength & \\multicolumn{4}{c}{Amplitude response} & Max shift \\\\\n"
        "\\cmidrule(lr){2-5}\n"
        "(m) & $0^\\circ$ & $45^\\circ$ & $90^\\circ$ & $135^\\circ$ & (m) \\\\\n"
        "\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\n"
        "\\end{tabular}\n\\end{table}\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=3)
    ap.add_argument("--modes", nargs="*", default=["nearest"],
                    choices=["nearest", "knn", "psf"])
    ap.add_argument("--k-points", type=int, default=10)
    ap.add_argument("--floor", type=float, default=None,
                    help="representation floor to compare against; read from "
                         "results/common_reference_floor.csv when omitted")
    ap.add_argument("--no-publish", action="store_true")
    ap.add_argument("--no-tex", action="store_true")
    args = ap.parse_args(argv)

    scenes = args.regions or list(region_sites.load_sites()["scene"])
    scenes = [s for s in scenes if (cfg.TILE_DIR / f"{s}.npz").exists()]
    if not scenes:
        print(f"no tiles under {cfg.TILE_DIR}; run build_tiles.py first")
        return 1
    if args.limit:
        scenes = scenes[:args.limit]

    floor = args.floor
    if floor is None:
        f = cfg.RESULTS_DIR / "common_reference_floor.csv"
        if f.exists():
            col = [c for c in pd.read_csv(f).columns if "floor" in c.lower()]
            if col:
                floor = float(pd.read_csv(f)[col[0]].median())

    print(f"{len(scenes)} region(s), modes {args.modes}")
    rows, failures = [], []
    for mode in args.modes:
        for scene in scenes:
            print(f"  [{mode}] {scene}")
            try:
                rows.extend(run_region(scene, mode, k_points=args.k_points))
            except Exception as exc:                       # noqa: BLE001
                print(f"    FAILED {type(exc).__name__}: {exc}")
                failures.append((mode, scene, f"{type(exc).__name__}: {exc}"))
    if not rows:
        print("no region produced any row")
        return 1

    df = pd.DataFrame(rows)
    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    csv = cfg.RESULTS_DIR / "operator_tests.csv"
    df.to_csv(csv, index=False)
    md = summarize(df, floor)
    (cfg.RESULTS_DIR / "operator_tests.md").write_text(md, encoding="utf-8")
    print("\n" + md)
    print(f"wrote {csv} ({len(df)} rows)")

    figure(df, publish=not args.no_publish)
    if not args.no_tex:
        tex = latex_table(df)
        if tex:
            out = cfg.TEX_TABLES / "table_operator_tests.tex"
            cfg.TEX_TABLES.mkdir(parents=True, exist_ok=True)
            out.write_text(tex, encoding="utf-8")
            print(f"wrote {out}")

    if failures:
        print(f"\n{len(failures)} region/mode pair(s) failed:")
        for m, s, why in failures:
            print(f"  {m} {s}: {why}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
