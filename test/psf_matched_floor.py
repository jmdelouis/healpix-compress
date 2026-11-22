#!/usr/bin/env python
# Usage (from test/, with healpix-resample at tag jstars-submission on PYTHONPATH and
# the four regions' tiles and windows from the data deposit under data/):
#     python psf_matched_floor.py
"""PSF-aware tiles of four regions: floor against the observation under
several reads, at the deposited damping (1e-3, rebuilt as a check) and at
the trial's damping (0.1). healpix-compress review M5.3, 2026-11-22.

healpix-resample at tag jstars-submission: PSFResampler(sigma_m=s) uses the
kernel exp(-2 d^2 / s^2), so s = 2 sigma and fwhm_to_scale(12.5) = 10.62 m.
The run behind results/psf_matched_floor.json used common_reference.py as it
stood before 2026-11-22, which passed sigma = FWHM / 2.355 as sigma_m:
fwhm_m=12.5 gave a 6.25 m kernel (the published psf/psf block) and fwhm_m=25
the 12.5 m kernel the tiles were built with. common_reference.py now uses the
package's convention, so the same two kernels are fwhm_m=6.25 and 12.5 below.
"""
import json, os, shutil, sys, time, traceback
import numpy as np
sys.path.insert(0, os.getcwd())
import paper_config as cfg
import build_tiles, common_reference
import healpix_resample as hr

OUT = os.path.join(str(cfg.RESULTS_DIR), "psf_matched_floor.json")
SCENES = ["region__agriculture__ukraine", "region__forest__black_forest",
          "region__urban__montreal", "region__water__lake_taupo"]
res = {"healpix_resample_file": hr.__file__, "rows": [], "errors": []}
try:
    from healpix_resample.psf_geometry import fwhm_to_scale
    res["fwhm_to_scale_12.5"] = float(fwhm_to_scale(12.5))
except Exception as e:
    res["fwhm_to_scale_12.5"] = repr(e)


def save():
    with open(OUT, "w") as fh:
        json.dump(res, fh, indent=1)


READS = (("nearest", "nearest", 12.5), ("bilinear", "bilinear", 12.5),
         ("psf_fwhm6.25_published", "psf", 6.25), ("psf_fwhm12.5_matched", "psf", 12.5))
for s in SCENES:
    try:
        store = cfg.REGION_DIR / f"{s}_data.zarr"
        dep = build_tiles.tile_path(s, "psf")
        keep = dep.with_name(f"{s}__psf__deposit.npz")
        if not keep.exists():
            shutil.copy(dep, keep)
        variants = {"deposit_lam1e-3": keep}
        for lam in (1e-3, 0.1):
            p = dep.with_name(f"{s}__psf__lam{lam:g}.npz")
            if not p.exists():
                t0 = time.time()
                build_tiles.build_tile(s, store, kernel="psf", lam=lam, verbose=False)
                shutil.move(dep, p)
                res["rows"].append(dict(scene=s, step="build", lam=lam,
                                        seconds=round(time.time() - t0, 1)))
                save()
            variants[f"rebuilt_lam{lam:g}"] = p
        shutil.copy(keep, dep)
        x_dep = np.load(keep, allow_pickle=True)["cell_data"].astype(float).ravel()
        for tag, p in variants.items():
            npz = np.load(p, allow_pickle=True)
            x = npz["cell_data"].astype(float).ravel()
            obs = np.asarray(npz["utm_crop"], dtype=float)
            row = dict(scene=s, variant=tag,
                       rms_diff_vs_deposit=float(np.sqrt(np.mean((x - x_dep) ** 2))))
            for name, mode, fw in READS:
                t0 = time.time()
                try:
                    cr = common_reference.build_for_tile(npz, mode=mode, fwhm_m=fw)
                    fl = float(cr.rmse_on_mask(cr.predict(x), obs))
                    sd = float(obs[cr.mask].std())
                    row[name] = fl
                    row[name + "_frac_std"] = fl / sd
                    row[name + "_scale_m"] = float(getattr(cr, "scale_m", float("nan")))
                except Exception as e:
                    row[name] = repr(e)[:300]
                row[name + "_s"] = round(time.time() - t0, 1)
            res["rows"].append(row)
            save()
    except Exception:
        res["errors"].append(dict(scene=s, tb=traceback.format_exc()[-1500:]))
        save()
res["done"] = True
save()

