"""Common reconstruction reference for both grids (paper Limitation 5).

The paper measures HEALPix RMSE against its own already-resampled tile and UTM
RMSE against the raw UTM crop: two different notions of distortion, and the one
place where the two grids' different pixel orientations really bite. This module
removes that asymmetry by bringing the decoded HEALPix field onto the native UTM
pixels, so both grids are scored against the same Sentinel-2 observation, on the
same pixels, in the same frame.

Two numbers come out of it, and the paper needs both:

  * ``rmse_common``  -- decoded HEALPix, read on the UTM pixels, against the
                        observed crop, restricted to the intersection of the two
                        footprints (the tile is a rotated square, the crop is
                        axis-aligned: they enclose the same area but not the
                        same ground, and a UTM pixel outside the tile has no
                        cell behind it);
  * ``floor``        -- the same projection applied to the *uncompressed* tile.
                        This is the carrier grid's representation error: no
                        codec setting can go below it, and quoting
                        ``rmse_common`` without it would charge the resampling
                        error to the codec.

Projection operator
-------------------
The default (``mode="nearest"``) reads each UTM pixel from the level-20 cell
that contains it: ``cell_of(lon, lat)`` is one ``healpix_geo`` call, exact, and
introduces no smoothing of its own, so the floor it measures is the carrier
grid's own error and nothing else. It answers "what does the stored HEALPix
field say at this location", which is the question a data cube's reader asks.

``mode="knn"`` / ``mode="psf"`` instead apply healpix-resample's Gaussian
observation operator, answering "what would the instrument have measured given
this stored field". That path is kept for the sensitivity check the paper
mentions, but it is not the default: on the environment used here the sparse
operator returns an all-zero projection, which is invisible in the results (every
sweep point reports the same "error", the RMS of the observations themselves) and
is guarded against below.
"""

from __future__ import annotations

import numpy as np

import paper_config as cfg

DEFAULT_FWHM_M = 12.5     # healpix-resample's calibrated effective B04 response

# The weighted-operator fallback is a property of the environment, not of a
# region: once it has been reported, repeating it 39 more times adds nothing.
# The mode actually used is recorded per region in common_reference_floor.csv.
_FALLBACK_REPORTED = False


class CommonReference:
    """Read the HEALPix tile on the UTM crop's pixels, where both overlap."""

    def __init__(self, cell_ids, utm_lon, utm_lat, *, mode="nearest",
                 fwhm_m=DEFAULT_FWHM_M, npt=None, erode=3, verbose=False):
        import healpix_geo

        self.shape = tuple(np.shape(utm_lon))
        self.cell_ids = np.asarray(cell_ids, dtype=np.int64)
        self.base_id = int(self.cell_ids[0])
        self.parent = self.base_id // cfg.N_CHILDREN
        self.mode = mode

        # --- which level-20 cell does each UTM pixel fall in? ----------------
        ipix = healpix_geo.nested.lonlat_to_healpix(
            np.asarray(utm_lon, dtype=np.float64).ravel(),
            np.asarray(utm_lat, dtype=np.float64).ravel(),
            cfg.LEVEL, ellipsoid="WGS84").astype(np.int64)
        self._child_index = (ipix - self.base_id).reshape(self.shape)

        # --- where can the two grids be compared at all? ---------------------
        inside = (ipix // cfg.N_CHILDREN == self.parent).reshape(self.shape)
        self.raw_inside_fraction = float(inside.mean())
        if erode:
            from scipy.ndimage import binary_erosion
            inside = binary_erosion(inside, iterations=int(erode), border_value=0)
        self.mask = inside
        self.mask_fraction = float(inside.mean())
        if self.mask_fraction < 0.2:
            raise RuntimeError(
                f"only {100*self.mask_fraction:.1f}% of the UTM crop lies inside "
                f"the HEALPix tile footprint (before erosion "
                f"{100*self.raw_inside_fraction:.1f}%); the two are not comparable")

        # Indices are only meaningful on the mask; clamp elsewhere so the gather
        # is always in range, and never read those pixels.
        self._safe_index = np.clip(self._child_index, 0, cfg.N_CHILDREN - 1)

        self.op = None
        self.scale_m = np.nan
        self.npt = npt
        if mode in ("knn", "psf"):
            self._build_operator(utm_lon, utm_lat, fwhm_m, npt, verbose)

        # --- the invariant every projection must satisfy ---------------------
        # A projection that loses a constant field is worthless and, worse,
        # invisible: every rate-distortion point then reports the same "error",
        # the RMS of the observations themselves. The weighted operator does
        # exactly that on some torch builds, so a failure there falls back to
        # the nearest-cell read -- which is verified, not a guess -- rather than
        # dropping the region. Only a nearest-cell failure is fatal, because it
        # would mean the tile and the crop are genuinely misregistered.
        self.support_median = self._constant_field_check()
        if not self._support_ok(self.support_median) and self.op is not None:
            global _FALLBACK_REPORTED
            if not _FALLBACK_REPORTED:
                print(f"    NOTE: the {mode!r} operator projects a constant field "
                      f"to {self.support_median:.4f} instead of 1.0 and is not "
                      f"usable in this environment; falling back to the "
                      f"nearest-cell read for every region (recorded as mode "
                      f"'nearest' in the results). Pass --cr-mode nearest to "
                      f"skip building it at all.")
                _FALLBACK_REPORTED = True
            self.op = None
            self.mode = f"nearest (fallback from {mode})"
            self.support_median = self._constant_field_check()
        if not self._support_ok(self.support_median):
            raise RuntimeError(
                f"the HEALPix -> UTM projection does not preserve a constant "
                f"field (median {self.support_median} on the comparable area, "
                f"expected 1.0) even with a nearest-cell read: the tile and the "
                f"crop must be misregistered.")

    @staticmethod
    def _support_ok(med) -> bool:
        return bool(np.isfinite(med)) and abs(float(med) - 1.0) <= 0.05

    def _constant_field_check(self) -> float:
        support = self.predict(np.ones(cfg.N_CHILDREN, dtype=np.float64))
        if not self.mask.any():
            return float("nan")
        return float(np.nanmedian(np.asarray(support)[self.mask]))

    # -- optional PSF-weighted operator --------------------------------------
    def _build_operator(self, utm_lon, utm_lat, fwhm_m, npt, verbose):
        """Build the reading operator the same way build_tiles.py does.

        healpix-resample's constructor arguments are keyword-only and have
        changed across versions; this path had drifted into passing four
        positionals and stayed broken unnoticed, because no published number
        used a mode other than `nearest`. It is now built from the same
        keyword call `build_tiles.make_resampler` uses, and the optional
        arguments are passed only if the installed class accepts them.
        """
        import inspect

        import healpix_resample as hr
        try:
            fwhm_to_scale = hr.fwhm_to_scale
        except AttributeError:
            from healpix_resample.psf_geometry import fwhm_to_scale
        self.scale_m = float(fwhm_to_scale(fwhm_m))
        cls = hr.PSFResampler if self.mode == "psf" else hr.KNeighborsResampler
        accepted = set(inspect.signature(cls.__init__).parameters)
        kw = dict(lon_deg=np.asarray(utm_lon, dtype=np.float64).ravel(),
                  lat_deg=np.asarray(utm_lat, dtype=np.float64).ravel(),
                  level=cfg.LEVEL, out_cell_ids=self.cell_ids,
                  verbose=verbose)
        if "sigma_m" in accepted:
            kw["sigma_m"] = self.scale_m
        if "npt" in accepted:
            if npt is None:
                from healpix_resample.psf_geometry import recommend_npt
                npt = int(recommend_npt(self.scale_m, cfg.LEVEL)["npt"])
            kw["npt"] = int(npt)
            self.npt = int(npt)
        else:
            self.npt = None
        missing = [k for k in kw if k not in accepted]
        if missing and not any(
                pr.kind is inspect.Parameter.VAR_KEYWORD
                for pr in inspect.signature(cls.__init__).parameters.values()):
            raise TypeError(
                f"{cls.__name__} does not accept {missing}; the installed "
                f"healpix-resample has a different constructor than this "
                f"code expects")
        self.op = cls(**kw)
        kept = np.asarray(self.op.get_cell_ids(), dtype=np.int64).ravel()
        self._kept_index = kept - self.base_id
        self.kept_fraction = kept.size / float(self.cell_ids.size)

    # -- projection -----------------------------------------------------------
    def predict(self, cell_values):
        """Read a full-block cell array on the UTM crop grid."""
        vals = np.asarray(cell_values, dtype=np.float64).ravel()
        if vals.size != cfg.N_CHILDREN:
            raise ValueError(f"expected {cfg.N_CHILDREN} cells, got {vals.size}")
        if self.op is None:
            return vals[self._safe_index]
        out = self.op.invert(vals[self._kept_index])
        return np.asarray(out, dtype=np.float64).reshape(self.shape)

    def read_at_level(self, cell_values, depth: int = 0):
        """Read the tile on the UTM grid after aggregating it `depth` levels up.

        HEALPix NESTED makes the aggregation exact: cells are equal-area and a
        contiguous NESTED block lays siblings out adjacently, so
        `reshape(-1, 4**depth).mean(axis=1)` is exactly the cell-average of the
        stored carrier over each level-(L-depth) cell.

        It is NOT the carrier the same transfer would have built at that level.
        The carrier holds point evaluations of an interpolant at cell centres,
        and for a general interpolant the value at the parent centre is not the
        mean of the values at the four child centres; the two differ by the
        curvature of the interpolant over the cell. What this returns is
        therefore a level-(L-depth) carrier built by interpolation followed by
        averaging, which bounds rather than reproduces the directly evaluated
        one. `level19_direct.py` builds the other object and measures the gap;
        `level_choice.py` reports this one and says so.

        The read stays the nearest containing cell at every depth, so a coarser
        level is charged twice: once for what the aggregation discards, once for
        a coarser read. `level_choice.py` separates the two.
        """
        d = int(depth)
        vals = np.asarray(cell_values, dtype=np.float64).ravel()
        if vals.size != cfg.N_CHILDREN:
            raise ValueError(f"expected {cfg.N_CHILDREN} cells, got {vals.size}")
        if d < 0 or 4 ** d > cfg.N_CHILDREN:
            raise ValueError(f"depth {d} out of range for this block")
        if d == 0:
            return vals[self._safe_index]
        agg = vals.reshape(-1, 4 ** d).mean(axis=1)
        idx = np.clip(self._safe_index // (4 ** d), 0, agg.size - 1)
        return agg[idx]

    @staticmethod
    def aggregate(cell_values, depth: int = 0):
        """The level-(L-depth) carrier, and the level-L field it prolongs to."""
        d = int(depth)
        vals = np.asarray(cell_values, dtype=np.float64).ravel()
        if d == 0:
            return vals, vals
        agg = vals.reshape(-1, 4 ** d).mean(axis=1)
        return agg, np.repeat(agg, 4 ** d)

    # -- the reader's kernel (idea 3) -----------------------------------------
    def read_weights(self):
        """How many mask pixels read each carrier cell. Nearest read only.

        For the nearest read the projection is a gather, so this vector *is*
        the diagonal of R^T R: a cell with weight zero cannot influence the
        product at all, whatever value it carries. The oversampling the paper
        charges for is exactly what makes that set large.
        """
        if self.op is not None:
            raise ValueError("read_weights is defined for the nearest read; "
                             "use read_support() for an operator read")
        idx = np.asarray(self._safe_index)[np.asarray(self.mask)].ravel()
        return np.bincount(idx, minlength=cfg.N_CHILDREN)

    def read_support(self):
        """Boolean per carrier cell: can it influence the masked read at all?

        The reader's kernel is the complement of this set. For the nearest
        read it is exact. For an operator read only the kept cells enter the
        projection, so the kept set is returned -- an *over*-estimate of the
        support, hence a conservative (under-stated) kernel.
        """
        if self.op is None:
            return self.read_weights() > 0
        sup = np.zeros(cfg.N_CHILDREN, dtype=bool)
        sup[np.asarray(self._kept_index, dtype=np.int64)] = True
        return sup

    def rmse_on_mask(self, predicted, observed) -> float:
        """RMSE over the two footprints' intersection only."""
        p = np.asarray(predicted, dtype=np.float64)[self.mask]
        o = np.asarray(observed, dtype=np.float64)[self.mask]
        ok = np.isfinite(p) & np.isfinite(o)
        return float(np.sqrt(np.mean((p[ok] - o[ok]) ** 2))) if ok.any() else np.nan

    def describe(self) -> str:
        extra = "" if self.op is None else (
            f", npt={self.npt}, scale={self.scale_m:.2f} m, cells kept "
            f"{100*self.kept_fraction:.1f}%")
        return (f"mode={self.mode} footprint overlap "
                f"{100*self.raw_inside_fraction:.1f}% -> comparable area "
                f"{100*self.mask_fraction:.1f}% of the crop, constant-field "
                f"check {self.support_median:.4f}{extra}")


def build_for_tile(npz, *, mode="nearest", fwhm_m=DEFAULT_FWHM_M, erode=3,
                   verbose=False):
    """Build a CommonReference from a tile .npz written by build_tiles.py."""
    parent = int(npz["parent_id"])
    cell_ids = parent * cfg.N_CHILDREN + np.arange(cfg.N_CHILDREN, dtype=np.int64)
    return CommonReference(cell_ids, npz["utm_lon"], npz["utm_lat"], mode=mode,
                           fwhm_m=fwhm_m, erode=erode, verbose=verbose)
