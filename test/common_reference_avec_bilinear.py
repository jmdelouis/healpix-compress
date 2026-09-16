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
        self._bil_idx = None
        self.offset_q99 = np.nan
        if mode in ("knn", "psf"):
            self._build_operator(utm_lon, utm_lat, fwhm_m, npt, verbose)
        elif mode == "bilinear":
            self._build_bilinear(utm_lon, utm_lat)
        elif mode != "nearest":
            raise ValueError(f"unknown reader mode {mode!r}")

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
        # The companion package is under active development and its helper
        # module has already moved once between two runs of the same day. The
        # width conversion is the plain Gaussian one (checked against the
        # companion's own values: 11.921 m FWHM -> 5.0625 m, ratio 2.3548), so
        # it is computed here when neither import resolves, rather than
        # letting the reader disappear from a run.
        fwhm_to_scale = getattr(hr, "fwhm_to_scale", None)
        if fwhm_to_scale is None:
            try:
                from healpix_resample.psf_geometry import fwhm_to_scale
            except ImportError:
                def fwhm_to_scale(f):
                    return float(f) / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        self.scale_m = float(fwhm_to_scale(fwhm_m))
        cls = hr.PSFResampler if self.mode == "psf" else hr.KNeighborsResampler
        accepted = set(inspect.signature(cls.__init__).parameters)
        kw = dict(lon_deg=np.asarray(utm_lon, dtype=np.float64).ravel(),
                  lat_deg=np.asarray(utm_lat, dtype=np.float64).ravel(),
                  level=cfg.LEVEL, out_cell_ids=self.cell_ids,
                  verbose=verbose)
        if "sigma_m" in accepted:
            kw["sigma_m"] = self.scale_m
        # npt is passed only when the caller chose one; otherwise the class's
        # own default applies, exactly as in build_tiles.make_resampler. The
        # old code derived it from a helper that no longer imports.
        if npt is not None and "npt" in accepted:
            kw["npt"] = int(npt)
        self.npt = int(npt) if (npt is not None and "npt" in accepted) else None
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

    def _build_bilinear(self, utm_lon, utm_lat):
        """An interpolating read: bilinear within the carrier's own lattice.

        A NESTED block's local index is the Morton code of the cell's integer
        (x, y) position on its base face, so the tile is a 1024 x 1024 lattice
        in which a cell's lattice neighbours are its spatial neighbours. The
        read places each UTM pixel in that lattice at sub-cell precision by
        inverting, at the containing cell, the local affine map from lattice
        steps to metres -- measured on the cell centres themselves, so cell
        anisotropy and orientation enter through the Jacobian rather than by
        assumption -- and takes the bilinear combination of the four bracketing
        centres. The four weights sum to one, so a constant field is preserved
        exactly, and the read is linear, so everything `transmission` and
        `read_support` say about a reader applies to it.

        This is the reconstructor Section V-F of the paper says it has not
        implemented: one that interpolates within the carrier raster instead of
        taking the containing cell, whose sub-cell jitter that section measures
        at about 2.1 m. Comparing the two floors on the same tiles says how much
        of the floor the read was responsible for.
        """
        import healpix_geo
        delta, side = cfg.DELTA, cfg.TILE_SIDE
        code = np.arange(cfg.N_CHILDREN, dtype=np.int64)
        ix = np.zeros_like(code)
        iy = np.zeros_like(code)
        for b in range(delta):                      # zorder._decode_zorder
            ix |= ((code >> (2 * b)) & 1) << b
            iy |= ((code >> (2 * b + 1)) & 1) << b

        def enc(x, y):
            out = np.zeros_like(x)
            for b in range(delta):
                out |= ((x >> b) & 1) << (2 * b)
                out |= ((y >> b) & 1) << (2 * b + 1)
            return out

        clon, clat = healpix_geo.nested.healpix_to_lonlat(
            self.cell_ids, cfg.LEVEL, ellipsoid="WGS84")
        clon = np.asarray(clon, dtype=np.float64)
        clat = np.asarray(clat, dtype=np.float64)
        # Longitudes arrive in whichever convention each library uses --
        # healpix_geo hands back [0, 360), the UTM crop carries [-180, 180) --
        # and the tile's mean is no safer than either. Every difference is
        # therefore wrapped to (-180, 180]; without this, the western
        # hemisphere reads as an offset of 2^22 cells and the guard below
        # rejects the region (it did, on 15 of 40).
        def wrap(d):
            return (np.asarray(d, dtype=np.float64) + 180.0) % 360.0 - 180.0

        lon0 = float(clon[0])
        lon0 = lon0 + float(np.mean(wrap(clon - lon0)))
        lat0 = float(clat.mean())
        k = cfg.R_EARTH * np.cos(np.deg2rad(lat0))

        def to_xy(lon, lat):                        # local metric frame
            return (np.deg2rad(wrap(lon - lon0)) * k,
                    np.deg2rad(lat - lat0) * cfg.R_EARTH)

        ex, ey = to_xy(clon, clat)
        E = np.empty((side, side))
        N = np.empty((side, side))
        E[iy, ix] = ex
        N[iy, ix] = ey

        px, py = to_xy(np.asarray(utm_lon, dtype=np.float64).ravel(),
                       np.asarray(utm_lat, dtype=np.float64).ravel())
        c = self._safe_index.ravel()
        cx, cy = ix[c], iy[c]
        xm, xp = np.clip(cx - 1, 0, side - 1), np.clip(cx + 1, 0, side - 1)
        ym, yp = np.clip(cy - 1, 0, side - 1), np.clip(cy + 1, 0, side - 1)
        dEdx = (E[cy, xp] - E[cy, xm]) / (xp - xm)
        dNdx = (N[cy, xp] - N[cy, xm]) / (xp - xm)
        dEdy = (E[yp, cx] - E[ym, cx]) / (yp - ym)
        dNdy = (N[yp, cx] - N[ym, cx]) / (yp - ym)
        det = dEdx * dNdy - dEdy * dNdx
        rx, ry = px - E[cy, cx], py - N[cy, cx]
        u = (dNdy * rx - dEdy * ry) / det
        v = (-dNdx * rx + dEdx * ry) / det
        fx, fy = cx + u, cy + v
        x0 = np.clip(np.floor(fx).astype(np.int64), 0, side - 2)
        y0 = np.clip(np.floor(fy).astype(np.int64), 0, side - 2)
        tx = np.clip(fx - x0, 0.0, 1.0)
        ty = np.clip(fy - y0, 0.0, 1.0)
        idx = np.stack([enc(x0, y0), enc(x0 + 1, y0),
                        enc(x0, y0 + 1), enc(x0 + 1, y0 + 1)], axis=-1)
        w = np.stack([(1 - tx) * (1 - ty), tx * (1 - ty),
                      (1 - tx) * ty, tx * ty], axis=-1)
        self._bil_idx = idx.reshape(self.shape + (4,))
        self._bil_w = w.reshape(self.shape + (4,))
        # A pixel lies inside its containing cell, which is a unit square of
        # the lattice, so |u|, |v| <= 1/2 up to the affine approximation. If
        # this quantile is not close to 1/2 the local frame is wrong.
        off = np.maximum(np.abs(u), np.abs(v))[np.asarray(self.mask).ravel()]
        self.offset_q99 = float(np.quantile(off, 0.99)) if off.size else np.nan
        if not (self.offset_q99 < 0.65):
            raise RuntimeError(
                f"bilinear read: 99th-percentile lattice offset "
                f"{self.offset_q99:.3f} cells, expected <= 0.5; the local "
                f"frame does not place pixels inside their containing cell")

    # -- projection -----------------------------------------------------------
    def predict(self, cell_values):
        """Read a full-block cell array on the UTM crop grid."""
        vals = np.asarray(cell_values, dtype=np.float64).ravel()
        if vals.size != cfg.N_CHILDREN:
            raise ValueError(f"expected {cfg.N_CHILDREN} cells, got {vals.size}")
        if self._bil_idx is not None:
            return np.einsum("...k,...k->...", vals[self._bil_idx], self._bil_w)
        if self.op is None:
            return vals[self._safe_index]
        out = self.op.invert(vals[self._kept_index])
        return np.asarray(out, dtype=np.float64).reshape(self.shape)

    def read_at_level(self, cell_values, depth: int = 0):
        """Read the tile on the UTM grid after aggregating it `depth` levels up.

        HEALPix NESTED makes this exact rather than approximate: cells are
        equal-area, so the value of a parent cell *is* the mean of its four
        children, and a contiguous NESTED block is laid out so that siblings
        are adjacent. Aggregating by `reshape(-1, 4**depth).mean(axis=1)` is
        therefore the level-(L-depth) carrier that would be built from the same
        transfer, not an approximation of it -- which is what lets the question
        "is a coarser level enough" be answered from the stored tiles instead of
        from a second pass over the imagery.

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
        if self.op is not None or self._bil_idx is not None:
            raise ValueError("read_weights is defined for the nearest read; "
                             "use read_support() for any other")
        idx = np.asarray(self._safe_index)[np.asarray(self.mask)].ravel()
        return np.bincount(idx, minlength=cfg.N_CHILDREN)

    def read_support(self):
        """Boolean per carrier cell: can it influence the masked read at all?

        The reader's kernel is the complement of this set. For the nearest
        read it is exact. For an operator read only the kept cells enter the
        projection, so the kept set is returned -- an *over*-estimate of the
        support, hence a conservative (under-stated) kernel.
        """
        if self._bil_idx is not None:
            m = np.asarray(self.mask)
            idx = self._bil_idx[m]
            w = self._bil_w[m]
            sup = np.zeros(cfg.N_CHILDREN, dtype=bool)
            sup[idx[w > 0]] = True
            return sup
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
        if self._bil_idx is not None:
            return (f"mode=bilinear footprint overlap "
                    f"{100*self.mask_fraction:.1f}%, constant-field check "
                    f"{self.support_median:.4f}, lattice offset q99 "
                    f"{self.offset_q99:.3f} cells")
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
