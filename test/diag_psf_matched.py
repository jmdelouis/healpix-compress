"""Matched synthetic test for the PSF-aware recovery gap (revision 7, point 4).

The companion study (healpix-resample, Sect. IV-B and IV-G, Table I) reports a
recovery of 15.4% (forest) to 19.6% (agriculture, urban) of the blur-driven
RMSE on its synthetic detector experiment. This paper's harness
(resampling_error.py, Appendix D) finds +3.0% to +6.5% without noise and
nothing at the measured noise level. Noise does not explain the gap, since it
persists in the noiseless control; the aggregation of the recovery fraction
does not either (mean of per-region ratios 3.5-6.3% against ratio of class
medians 3.0-6.7%, both from resampling_error_noise0.csv).

The two protocols differ in three things this script can change one at a time
on the SAME geometry with the SAME transfer code (healpix_resample.PSFResampler,
lambda = 1e-3, cell-average estimand over 16 level-22 children):

  latent field    : analytic power law with the region's fitted slope, band-
                    limited at the level-20 Nyquist (ours)  vs  a real high-
                    resolution texture at 2.5 m (Esri World Imagery, theirs);
  acquisition     : exact Gaussian convolution point-sampled at the native
                    pixel centres (ours)  vs  Gaussian blur on the 2.5 m grid
                    followed by a 4x4 block average, the detector aperture
                    (theirs; their inverse still assumes the Gaussian alone);
  response width  : the region's fitted FWHM, 10.1-13.7 m (ours)  vs  the
                    nominal 12.5 m for both generation and inversion (theirs).

Ladder of variants, each adding one change to the previous one:

  V0  our harness, noiseless (should reproduce Appendix D's noiseless column)
  V1  V0 with FWHM fixed at 12.5 m for generation and inversion
  V2  V1 with the latent field replaced by a texture on a 2.5 m grid,
      Gaussian-blurred and POINT-sampled at the pixel centres
  V3  V2 with the 4x4 block average (detector aperture) instead of point
      sampling  -- this is the companion's protocol
  V4  V3 restricted to a 256 x 256 native patch, as the companion's patches

The first variant whose recovery moves from ~5% towards ~15-20% names the
cause. Recovery follows the companion's Eq. (24):
    recovery = 1 - RMSE(reconstruction, latent) / RMSE(blurred, latent)
with both RMSE on the cell-average estimand, and is reported per region and
as the mean over regions (the companion's aggregation).

The texture: pass --texture PATH to a 2-D array (npz with key 'img', npy, or
a GeoTIFF readable by rasterio) sampled at 2.5 m, at least 4x the sub-window
in pixels; the companion's Esri patches, exported at zoom 17 and converted to
luminance, are the intended input. Without --texture the script builds a
PROXY: piecewise-constant parcels plus a power-law field with power up to the
2.5 m Nyquist. The proxy has sharp edges and high-frequency content like a
texture, but it is not Esri imagery and the report says so on every line.

Usage (where healpix_resample and the region stores are):

    python diag_psf_matched.py                       # one region per class
    python diag_psf_matched.py --texture esri_paris.npz --scenes region__urban__paris
    python diag_psf_matched.py --all --variants V0 V1 V3

Writes results/diag_psf_matched.md and results/diag_psf_matched.csv.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import paper_config as cfg
import build_tiles
import region_sites
import resampling_error as rerr
import support_controls as sup

LAM = 1e-3
NOMINAL_FWHM_M = 12.5
FINE_M = 2.5                      # the texture grid (native / 4)
OVERSAMPLE = int(round(cfg.NATIVE_PIXEL_M / FINE_M))   # 4
CHILDREN = 16                     # level-22 children per level-20 cell


# --- geometry -----------------------------------------------------------------
def _geometry(scene: str, halo_px: int = 32):
    store = cfg.REGION_DIR / f"{scene}_data.zarr"
    w = build_tiles.tile_subwindow(scene, store, halo_px=halo_px)
    cand, to_crs = w["cand"], w["to_crs"]
    r0, r1, c0, c1 = w["bounds"]
    x, y = np.asarray(w["x"], dtype=np.float64), np.asarray(w["y"], dtype=np.float64)
    dx = float(np.median(np.diff(x)))
    dy = float(np.median(np.diff(y)))
    # fractional (row, col) in the sub-window frame for any projected (X, Y)
    def index(X, Y):
        col = (np.asarray(X) - x[c0]) / dx
        row = (np.asarray(Y) - y[r0]) / dy
        return row, col
    import healpix_geo
    depth = int(round(math.log(CHILDREN, 4)))
    kids = (cand["cell_ids"][:, None] * (4 ** depth)
            + np.arange(4 ** depth, dtype=np.int64)[None, :]).ravel()
    klon, klat = healpix_geo.nested.healpix_to_lonlat(kids, cfg.LEVEL + depth,
                                                      ellipsoid="WGS84")
    kx, ky = to_crs.transform(np.asarray(klon, dtype=np.float64),
                              np.asarray(klat, dtype=np.float64))
    krow, kcol = index(kx, ky)
    return dict(scene=scene, store=store, w=w, cand=cand, to_crs=to_crs,
                shape=(r1 - r0, c1 - c0), bounds=(r0, r1, c0, c1),
                sub_lon=w["sub_lon"], sub_lat=w["sub_lat"],
                kid_rc=(krow, kcol), index=index)


def _sigma_pix_fine(fwhm_m):
    return sup.psf_sigma_from_fwhm(fwhm_m) / FINE_M


# --- latent textures ---------------------------------------------------------------
def load_texture(path: Path, shape_fine):
    p = Path(path)
    if p.suffix == ".npz":
        z = np.load(p)
        img = z["img"] if "img" in z else z[list(z.keys())[0]]
    elif p.suffix == ".npy":
        img = np.load(p)
    else:
        import rasterio
        with rasterio.open(p) as ds:
            img = ds.read(1)
    img = np.asarray(img, dtype=np.float64)
    if img.ndim == 3:
        img = img.mean(axis=0)
    H, W = shape_fine
    if img.shape[0] < H or img.shape[1] < W:
        raise ValueError(f"texture {img.shape} smaller than the sub-window at 2.5 m {shape_fine}")
    r = (img.shape[0] - H) // 2
    c = (img.shape[1] - W) // 2
    return img[r:r + H, c:c + W]


def proxy_texture(shape_fine, seed=5, n_parcels=400):
    """Piecewise-constant parcels plus a power-law field with power to 2.5 m.

    A stand-in for a high-resolution image: sharp edges everywhere and real
    energy above the native Nyquist. It is NOT Esri imagery.
    """
    H, W = shape_fine
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:H, 0:W]
    seeds = rng.uniform(0, [H, W], size=(n_parcels, 2))
    vals = rng.uniform(0.05, 0.35, size=n_parcels)
    from scipy.spatial import cKDTree
    _, idx = cKDTree(seeds).query(np.column_stack([yy.ravel(), xx.ravel()]))
    parcels = vals[idx].reshape(H, W)
    k_min = 2 * math.pi / (min(H, W) * FINE_M)
    k_max = math.pi / FINE_M
    px = (xx.ravel() - W / 2) * FINE_M
    py = (yy.ravel() - H / 2) * FINE_M
    tex = sup.synthetic_field(px, py, beta=1.9, k_min=k_min, k_max=k_max,
                              seed=seed, psf_sigma_m=0.0).reshape(H, W)
    tex = 0.03 * tex / (np.std(tex) or 1.0)
    return parcels + tex


# --- one variant ------------------------------------------------------------------
def _cell_average_from_fine(fine, geo):
    from scipy.ndimage import map_coordinates
    krow, kcol = geo["kid_rc"]
    v = map_coordinates(fine, [krow * OVERSAMPLE + (OVERSAMPLE - 1) / 2,
                               kcol * OVERSAMPLE + (OVERSAMPLE - 1) / 2],
                        order=1, mode="nearest")
    return v.reshape(-1, CHILDREN).mean(axis=1)


def _run_psf(geo, obs_flat, sigma, mask=None, cell_sel=None):
    """Run the PSF-aware transfer; `sigma` is the response's standard deviation.

    PSFResampler's `sigma_m` is the operator scale s = 2 sigma of the
    companion's Eq. (2), the quantity build_tiles passes for the real tiles
    (fwhm_to_scale). The first run of this script passed sigma itself, which
    gave the inverse a response half as wide as the blur, dropped 3-22% of the
    cells for lack of kernel weight and returned negative recoveries.
    """
    sub_lon, sub_lat = geo["sub_lon"], geo["sub_lat"]
    cell_ids = geo["cand"]["cell_ids"]
    if cell_sel is not None:
        cell_ids = cell_ids[cell_sel]
    if mask is not None:
        sub_lon, sub_lat = sub_lon[mask], sub_lat[mask]
        obs_flat = obs_flat[mask]
    op = build_tiles.make_resampler("psf", sub_lon, sub_lat, cell_ids,
                                    verbose=False, sigma_m=2.0 * float(sigma))
    res = op.resample(np.asarray(obs_flat, dtype=np.float64), lam=LAM)
    val, _n = build_tiles.complete_block(
        np.asarray(res.cell_ids, dtype=np.int64).ravel(),
        np.asarray(res.cell_data, dtype=np.float64).ravel(), cell_ids)
    return val


def _recovery(rec, blurred_cells, latent_cells, sel=None):
    if sel is not None:
        rec, blurred_cells, latent_cells = rec[sel], blurred_cells[sel], latent_cells[sel]
    e_psf = float(np.sqrt(np.mean((rec - latent_cells) ** 2)))
    e_blur = float(np.sqrt(np.mean((blurred_cells - latent_cells) ** 2)))
    return dict(rmse_blurred=e_blur, rmse_matched=e_psf,
                recovery=1.0 - e_psf / e_blur if e_blur > 0 else np.nan)


def variant_analytic(scene, geo, fwhm_override=None):
    """V0 / V1: this paper's harness, noiseless, via resampling_error.synthetic_case."""
    spectrum = None
    p = cfg.RESULTS_DIR / "support_controls_spectrum.csv"
    if p.exists():
        spectrum = pd.read_csv(p)
    if fwhm_override is not None:
        orig = rerr.scene_parameters
        rerr.scene_parameters = lambda sc, sp: (orig(sc, sp)[0], float(fwhm_override))
    try:
        case = rerr.synthetic_case(scene, geo["store"], spectrum, noise_scale=0.0,
                                   children=CHILDREN)
    finally:
        if fwhm_override is not None:
            rerr.scene_parameters = orig
    rec = _run_psf(geo, case.obs, case.sigma)
    out = _recovery(rec, case.observed_cells, case.latent_cells)
    out.update(fwhm_m=case.psf_fwhm_m, beta=case.beta)
    return out


def variant_texture(geo, fine, aperture: bool, fwhm_m=NOMINAL_FWHM_M, patch=None):
    """V2 / V3 / V4: texture latent, Gaussian blur on the 2.5 m grid, point or
    block-averaged acquisition, inversion with the Gaussian alone."""
    from scipy.ndimage import gaussian_filter, map_coordinates
    H, W = geo["shape"]
    blurred = gaussian_filter(fine, _sigma_pix_fine(fwhm_m), mode="reflect")
    if aperture:
        obs = blurred.reshape(H, OVERSAMPLE, W, OVERSAMPLE).mean(axis=(1, 3))
    else:
        rr, cc = np.mgrid[0:H, 0:W]
        obs = map_coordinates(blurred, [rr.ravel() * OVERSAMPLE + (OVERSAMPLE - 1) / 2,
                                        cc.ravel() * OVERSAMPLE + (OVERSAMPLE - 1) / 2],
                              order=1, mode="nearest").reshape(H, W)
    latent_cells = _cell_average_from_fine(fine, geo)
    blurred_cells = _cell_average_from_fine(blurred, geo)
    mask = None
    sel = None
    if patch is not None:
        # a 256 x 256 native patch centred on the tile: observations outside it
        # are withheld, and only cells whose footprint lies inside are scored
        rc, cc0 = H // 2, W // 2
        rr, cc = np.mgrid[0:H, 0:W]
        mask = ((np.abs(rr - rc) < patch // 2) & (np.abs(cc - cc0) < patch // 2)).ravel()
        krow, kcol = geo["kid_rc"]
        kr = krow.reshape(-1, CHILDREN).mean(axis=1)
        kc = kcol.reshape(-1, CHILDREN).mean(axis=1)
        margin = 20  # native pixels kept clear of the patch edge, as the companion does
        sel = (np.abs(kr - rc) < patch // 2 - margin) & (np.abs(kc - cc0) < patch // 2 - margin)
        if sel.sum() < 1000:
            raise RuntimeError(f"only {int(sel.sum())} cells inside the {patch}-pixel patch")
    rec = _run_psf(geo, obs.ravel(), sup.psf_sigma_from_fwhm(fwhm_m), mask=mask,
                   cell_sel=sel)
    if sel is not None:
        out = _recovery(rec, blurred_cells[sel], latent_cells[sel])
    else:
        out = _recovery(rec, blurred_cells, latent_cells)
    out.update(fwhm_m=fwhm_m, n_cells=int(sel.sum()) if sel is not None else int(rec.size))
    return out


VARIANTS = ("V0", "V1", "V2", "V3", "V4")


def run_scene(scene, variants, texture_path=None):
    geo = _geometry(scene)
    H, W = geo["shape"]
    fine = None
    if any(v in variants for v in ("V2", "V3", "V4")):
        if texture_path:
            fine, source = load_texture(texture_path, (H * OVERSAMPLE, W * OVERSAMPLE)), Path(texture_path).name
        else:
            fine, source = proxy_texture((H * OVERSAMPLE, W * OVERSAMPLE)), "PROXY (parcels + power law)"
    rows = []
    for v in variants:
        try:
            if v == "V0":
                r = variant_analytic(scene, geo); desc = "analytic field, fitted FWHM, point-sampled, noiseless"
            elif v == "V1":
                r = variant_analytic(scene, geo, fwhm_override=NOMINAL_FWHM_M); desc = "analytic field, FWHM 12.5 m, point-sampled"
            elif v == "V2":
                r = variant_texture(geo, fine, aperture=False); desc = f"texture [{source}], Gaussian 12.5 m, point-sampled"
            elif v == "V3":
                r = variant_texture(geo, fine, aperture=True); desc = f"texture [{source}], Gaussian 12.5 m + 4x4 aperture (companion protocol)"
            elif v == "V4":
                r = variant_texture(geo, fine, aperture=True, patch=256); desc = f"texture [{source}], companion protocol, 256 x 256 patch"
            else:
                continue
            r.update(scene=scene, variant=v, description=desc)
        except Exception as e:  # noqa: BLE001
            r = dict(scene=scene, variant=v, description=f"failed: {e}", recovery=np.nan,
                     rmse_blurred=np.nan, rmse_matched=np.nan)
        rows.append(r)
        print(f"{scene} {v}: recovery {100*r.get('recovery', np.nan):.1f}%  "
              f"blur {r.get('rmse_blurred', np.nan):.4g}  matched {r.get('rmse_matched', np.nan):.4g}  ({r['description']})")
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--scenes", nargs="*", default=None)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--variants", nargs="*", default=list(VARIANTS))
    ap.add_argument("--texture", default=None, help="2-D array at 2.5 m (npz/npy/GeoTIFF)")
    args = ap.parse_args(argv)
    sites = region_sites.load_sites()
    if args.all:
        scenes = list(sites["scene"])
    elif args.scenes:
        scenes = args.scenes
    else:
        scenes = list(sites.groupby("scene_class").head(1)["scene"])
    rows = []
    for sc in scenes:
        rows += run_scene(sc, args.variants, args.texture)
    df = pd.DataFrame(rows)
    df["scene_class"] = df.scene.str.split("__").str[1]
    out_csv = cfg.RESULTS_DIR / "diag_psf_matched.csv"
    df.to_csv(out_csv, index=False)
    lines = ["# Matched PSF-aware recovery test", "",
             "recovery = 1 - RMSE(matched, latent) / RMSE(blurred, latent), cell-average "
             "estimand (16 level-22 children), lambda = 1e-3, healpix_resample.PSFResampler.",
             "The companion reports 15.4-19.6% (Table I). This paper's Appendix D "
             "reports +3.0 to +6.5% noiseless (V0 should reproduce it).", ""]
    if args.texture is None and any(v in args.variants for v in ("V2", "V3", "V4")):
        lines += ["**V2-V4 use a PROXY texture (parcels + power law), not Esri imagery.** "
                  "Rerun with --texture on an exported Esri patch before drawing a conclusion.", ""]
    lines += ["| variant | description | mean recovery over regions | per class |", "|---|---|---|---|"]
    for v in args.variants:
        g = df[df.variant == v]
        if not len(g):
            continue
        per = ", ".join(f"{c} {100*x:.1f}%" for c, x in g.groupby("scene_class").recovery.mean().items())
        lines.append(f"| {v} | {g.description.iloc[0]} | {100*g.recovery.mean():.1f}% | {per} |")
    lines += ["", "Per region:", "", df.to_string(index=False), ""]
    out_md = cfg.RESULTS_DIR / "diag_psf_matched.md"
    out_md.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {out_md} and {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
