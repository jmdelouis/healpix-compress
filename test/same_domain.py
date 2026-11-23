"""Same-ground control of the support mismatch (external review, 2026-11-23).

The carrier tile and the native crop have equal areas but not the same
ground: they overlap over 74-83% of it. Every ratio of the main comparison is
therefore taken between two rasters that partly see different scenes. The
controls of Table VII (tab:support_controls) bound that effect or remove the content altogether; the
common-support row (Appendix B until 2026-11-23) compared a 256 x 256 cell block with a
native rectangle five to eight times larger, which is biased at the loose
end where fixed costs weigh. None compresses exactly the same ground on both
grids at full raster size, which is what this control does.

The domain D is the intersection of the two footprints, sampled by each grid:

    native crop   pixels whose centre lies in the carrier tile, i.e. whose
                  level-20 cell has the tile as its level-10 parent;
    carrier       cells whose centre lies inside the crop square.

Outside D each raster is replaced by a fill computed from its own values on
D, and the error is measured on D only. Both rasters keep their lattice and
their full shape, so the codecs see the layout of the main comparison, and
since tile and crop have equal areas the filled fraction is about the same
on both sides (the same ground, the complement of the intersection). Three
variants are swept per region, all with the paper's codec ladders, the
oracle-best codec and the bracketing rule on both grids:

    full          the unmodified rasters, error over the whole raster: the
                  main comparison, recomputed here with the same ladders so
                  that each region's same-ground ratio has its own baseline;
    on_domain     the unmodified rasters, error over D only (no fill);
    harmonic      fill = harmonic (Laplace) extension of the values on D's
                  boundary: smooth, cheap to code, no edge along D;
    mean          fill = mean over D: a step along D's boundary, as a check
                  that the fill does not decide the ratio.

Usage (the paper's run, 40 regions, about ten minutes on 16 cores)
-----
    python same_domain.py --workers 12
    python same_domain.py --regions region__agriculture__ukraine
    python same_domain.py --summary-only

Writes results/same_domain_parts/<scene>.csv (one per region, resumable),
results/same_domain_sweeps.csv, results/same_domain.csv and
results/same_domain.json.
"""

from __future__ import annotations

import argparse
import json
import os
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

import paper_config as cfg

PARTS = cfg.RESULTS_DIR / "same_domain_parts"
VARIANTS = ("full", "on_domain", "harmonic", "mean")


# --- domain -------------------------------------------------------------------
def domain_masks(reg):
    """(utm_mask, hpx_mask_2d): the intersection D sampled by each grid."""
    import healpix_geo
    import support_controls as sc

    lon, lat = reg["utm_lon"], reg["utm_lat"]
    cells = np.asarray(healpix_geo.nested.lonlat_to_healpix(
        lon.ravel(), lat.ravel(), cfg.LEVEL, ellipsoid="WGS84"), dtype=np.int64)
    utm_mask = (cells // cfg.N_CHILDREN == reg["parent"]).reshape(lon.shape)

    clon, clat = sc.cell_lonlat(reg["cell_ids"])
    col, row = sc.crop_pixel_index(clon, clat, reg["crs_wkt"],
                                   reg["utm_x"], reg["utm_y"])
    n = reg["utm_crop"].shape[0]
    inside = ((col >= -0.5) & (col < n - 0.5) & (row >= -0.5) & (row < n - 0.5))
    hpx_mask = sc.expand_bool(inside)
    return utm_mask, hpx_mask


def harmonic_fill(img, mask):
    """Solve Laplace's equation outside `mask`, with `img` on `mask` as data.

    Five-point stencil, reflecting (Neumann) at the raster border. Every
    unknown region touches D, so the system is non-singular.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import spsolve

    img = np.asarray(img, dtype=np.float64)
    H, W = img.shape
    unk = ~np.asarray(mask, dtype=bool)
    n = int(unk.sum())
    if n == 0:
        return img.astype(np.float32)
    idx = -np.ones((H, W), dtype=np.int64)
    idx[unk] = np.arange(n)
    ys, xs = np.nonzero(unk)
    p_all = idx[ys, xs]
    deg = np.zeros(n)
    b = np.zeros(n)
    off_r, off_c = [], []
    for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        ny, nx = ys + dy, xs + dx
        ok = (ny >= 0) & (ny < H) & (nx >= 0) & (nx < W)
        p = p_all[ok]
        ny, nx = ny[ok], nx[ok]
        np.add.at(deg, p, 1.0)
        q_unk = unk[ny, nx]
        off_r.append(p[q_unk])
        off_c.append(idx[ny[q_unk], nx[q_unk]])
        np.add.at(b, p[~q_unk], img[ny[~q_unk], nx[~q_unk]])
    r = np.concatenate([p_all] + off_r)
    c = np.concatenate([p_all] + off_c)
    v = np.concatenate([deg] + [-np.ones(len(x)) for x in off_r])
    A = coo_matrix((v, (r, c)), shape=(n, n)).tocsc()
    out = img.copy()
    out[unk] = spsolve(A, b)
    return out.astype(np.float32)


def filled(img, mask, how):
    img = np.asarray(img, dtype=np.float32)
    if how == "harmonic":
        return harmonic_fill(img, mask)
    if how == "mean":
        out = img.copy()
        out[~mask] = float(img[mask].mean())
        return out
    raise ValueError(how)


# --- sweeps -------------------------------------------------------------------
def sweep(arr, mask, scene, grid, variant, codecs):
    """One rate-distortion sweep; error over the whole raster and over D."""
    import codec_lib
    a = np.asarray(arr, dtype=np.float32)
    rows = []
    for run in codec_lib.iter_runs(a, codecs=codecs, on_error="warn"):
        d = np.asarray(run.reconstruction, dtype=np.float64) - a
        rows.append(dict(scene=scene, grid=grid, variant=variant,
                         codec=run.codec,
                         config=json.dumps(run.config, sort_keys=True),
                         total_bytes=int(run.total_bytes),
                         rmse_full=float(np.sqrt(np.mean(d * d))),
                         rmse_domain=float(np.sqrt(np.mean(d[mask] ** 2)))))
    return rows


def run_region(scene, codecs, force=False):
    PARTS.mkdir(parents=True, exist_ok=True)
    out = PARTS / f"{scene}.csv"
    if out.exists() and not force:
        return f"{scene}: cached"
    import support_controls as sc
    import build_tiles  # noqa: F401  (load_region needs it)
    t0 = time.time()
    reg = sc.load_region(scene)
    utm = np.asarray(reg["utm_crop"], dtype=np.float32)
    hpx = sc.expanded(reg["cell_data"])
    utm_mask, hpx_mask = domain_masks(reg)
    rows = []
    rows += sweep(hpx, hpx_mask, scene, "hpx", "orig", codecs)
    rows += sweep(utm, utm_mask, scene, "utm", "orig", codecs)
    for how in ("harmonic", "mean"):
        rows += sweep(filled(hpx, hpx_mask, how), hpx_mask, scene, "hpx", how, codecs)
        rows += sweep(filled(utm, utm_mask, how), utm_mask, scene, "utm", how, codecs)
    df = pd.DataFrame(rows)
    df["utm_domain_frac"] = float(utm_mask.mean())
    df["hpx_domain_frac"] = float(hpx_mask.mean())
    df["utm_domain_km2"] = float(utm_mask.sum() * cfg.NATIVE_PIXEL_M ** 2 / 1e6)
    df["hpx_domain_km2"] = float(hpx_mask.sum() * cfg.DELTA_EQ ** 2 / 1e6)
    df["utm_std_domain"] = float(utm[utm_mask].std())
    df["hpx_std_domain"] = float(hpx[hpx_mask].std())
    tmp = out.with_suffix(".tmp")
    df.to_csv(tmp, index=False)
    tmp.replace(out)
    return (f"{scene}: {len(df)} points, D = {100 * utm_mask.mean():.1f}% of the "
            f"crop, {100 * hpx_mask.mean():.1f}% of the tile, "
            f"{time.time() - t0:.0f} s")


def _worker(args):
    scene, codecs, force = args
    try:
        return run_region(scene, codecs, force)
    except Exception:                                         # noqa: BLE001
        return f"{scene}: FAILED\n{traceback.format_exc()[-1500:]}"


# --- ratios -------------------------------------------------------------------
def variant_frame(sw, variant):
    """Sweep rows of one variant, with the error column the variant scores."""
    src = "orig" if variant in ("full", "on_domain") else variant
    col = "rmse_full" if variant == "full" else "rmse_domain"
    d = sw[sw.variant == src].copy()
    d["rmse"] = d[col]
    return d


def ratios(sw, codecs, targets):
    import paper_stats as st
    rows = []
    for variant in VARIANTS:
        d = variant_frame(sw, variant)
        for scene, g in d.groupby("scene"):
            h, u = g[g.grid == "hpx"], g[g.grid == "utm"]
            for t in targets:
                hb = st.best_of_codecs_size(h, scene, codecs, t)
                ub = st.best_of_codecs_size(u, scene, codecs, t)
                bh = st.oracle_bracketed(h, scene, codecs, t, "rmse")
                bu = st.oracle_bracketed(u, scene, codecs, t, "rmse")
                rows.append(dict(
                    scene=scene, variant=variant, target=t,
                    hpx_bytes=hb, utm_bytes=ub,
                    ratio=(hb / ub if np.isfinite(hb) and np.isfinite(ub) and ub > 0
                           else np.nan),
                    hpx_codec=st.argbest_codec_size(h, scene, codecs, t)[0],
                    utm_codec=st.argbest_codec_size(u, scene, codecs, t)[0],
                    bracketed_hpx=bh, bracketed_utm=bu,
                    both_bracketed=bool(bh and bu),
                    utm_domain_frac=float(g.utm_domain_frac.iloc[0]),
                    hpx_domain_frac=float(g.hpx_domain_frac.iloc[0])))
    return pd.DataFrame(rows)


def summarize(r, targets):
    import paper_stats as st
    out = {"targets": [float(t) for t in targets], "variants": {}}
    for variant in VARIANTS:
        rows = []
        for t in targets:
            d = r[(r.variant == variant) & (r.target == t) & r.both_bracketed]
            test = st.paired_ratio_test(d.ratio.values)
            rows.append(dict(target=float(t), n=test["n"], median=test["median"],
                             q25=test["q25"], q75=test["q75"], p=test["p"],
                             p_sign=test["p_sign"], ci_lo=test["ci_lo"],
                             ci_hi=test["ci_hi"],
                             n_below_one=int((d.ratio < 1).sum())))
        out["variants"][variant] = rows
    # region by region: how far does restricting to the same ground move the ratio?
    shifts = {}
    for variant in ("on_domain", "harmonic", "mean"):
        rr = []
        for t in targets:
            a = r[(r.variant == "full") & (r.target == t) & r.both_bracketed]
            b = r[(r.variant == variant) & (r.target == t) & r.both_bracketed]
            m = a.merge(b, on="scene", suffixes=("_full", "_v"))
            q = (m.ratio_v / m.ratio_full).values
            test = st.paired_ratio_test(q)
            rr.append(dict(target=float(t), n=test["n"], median=test["median"],
                           q25=test["q25"], q75=test["q75"], p=test["p"]))
        shifts[variant] = rr
    out["over_full"] = shifts
    d = r[(r.variant == "harmonic") & (r.target == targets[0])]
    out["domain_frac_utm"] = [float(d.utm_domain_frac.min()),
                              float(d.utm_domain_frac.median()),
                              float(d.utm_domain_frac.max())]
    out["domain_frac_hpx"] = [float(d.hpx_domain_frac.min()),
                              float(d.hpx_domain_frac.median()),
                              float(d.hpx_domain_frac.max())]
    out["n_regions"] = int(r.scene.nunique())
    return out


def compare_with_paper(r, targets):
    """Region-by-region agreement of the `full` variant with the main bench."""
    import paper_stats as st
    path = cfg.RESULTS_DIR / "paper_numbers_self.json"
    if not path.exists():
        return None
    doc = json.loads(path.read_text())
    out = []
    for t in targets:
        tt = [x for x in doc["targets"] if x["kind"] == "rmse"
              and abs(float(x["target"]) - t) < 1e-12]
        if not tt:
            continue
        # stored at the illustrative overhead: ratio = HPX / (gamma * UTM)
        g = float(doc.get("gamma_utm", 1.0))
        paper = {k: float(v) * g for k, v in tt[0].get("ratios", {}).items()}
        d = r[(r.variant == "full") & (r.target == t) & r.both_bracketed]
        common = [s for s in d.scene if s in paper]
        if not common:
            continue
        mine = d.set_index("scene").ratio.loc[common].values
        theirs = np.array([paper[s] for s in common], dtype=float)
        out.append(dict(target=t, n_common=len(common), n_paper=len(paper),
                        n_here=int(len(d)),
                        median_here=float(np.median(mine)),
                        median_paper=float(np.median(theirs)),
                        max_abs_rel_diff=float(np.max(np.abs(mine / theirs - 1)))))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--summary-only", action="store_true")
    ap.add_argument("--codecs", nargs="*", default=cfg.COMPRESSOR_CODECS)
    ap.add_argument("--targets", type=float, nargs="*", default=cfg.TARGET_RMSE)
    args = ap.parse_args(argv)

    import region_sites
    scenes = list(region_sites.select(region_sites.load_sites(), args.regions,
                                      None, None).scene)
    if not args.summary_only:
        jobs = [(s, args.codecs, args.force) for s in scenes]
        if args.workers > 1:
            from multiprocessing import get_context
            with get_context("spawn").Pool(args.workers) as pool:
                for msg in pool.imap_unordered(_worker, jobs):
                    print(msg, flush=True)
        else:
            for j in jobs:
                print(_worker(j), flush=True)

    parts = [PARTS / f"{s}.csv" for s in scenes if (PARTS / f"{s}.csv").exists()]
    if parts:
        sw = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
        sw.to_csv(cfg.RESULTS_DIR / "same_domain_sweeps.csv", index=False)
    elif (cfg.RESULTS_DIR / "same_domain_sweeps.csv").exists():
        # a checkout with the assembled sweeps only
        sw = pd.read_csv(cfg.RESULTS_DIR / "same_domain_sweeps.csv")
    else:
        print("no sweeps")
        return 1
    r = ratios(sw, args.codecs, args.targets)
    r.to_csv(cfg.RESULTS_DIR / "same_domain.csv", index=False)
    summary = summarize(r, args.targets)
    summary["paper_full_agreement"] = compare_with_paper(r, args.targets)
    summary["codecs"] = list(args.codecs)
    (cfg.RESULTS_DIR / "same_domain.json").write_text(json.dumps(summary, indent=1))
    for variant, rows in summary["variants"].items():
        print(f"{variant:10s} " + "  ".join(
            f"{x['target']:g}: {x['median']:.3f} ({x['n']}, p={x['p']:.3g})"
            for x in rows))
    for variant, rows in summary["over_full"].items():
        print(f"{variant:10s}/full " + "  ".join(
            f"{x['target']:g}: {x['median']:.3f} ({x['n']})" for x in rows))
    print("D as a fraction of the crop (min, median, max):", summary["domain_frac_utm"])
    print("D as a fraction of the tile:", summary["domain_frac_hpx"])
    print("agreement with the paper:", summary["paper_full_agreement"])
    return 0


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    raise SystemExit(main())
