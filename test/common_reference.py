"""Common reconstruction reference for both grids (paper Limitation 5).

Today the paper measures HEALPix RMSE against its own already-resampled tile
and UTM RMSE against the raw UTM crop: two different notions of distortion.
This module closes that gap by pushing the decoded HEALPix field back onto the
native UTM pixels with healpix-resample's forward observation operator
(``KNeighborsResampler.invert``), so both grids are scored against the *same*
Sentinel-2 observation.

Two numbers come out of it, and the paper needs both:

  * ``rmse_common``  -- decoded HEALPix -> UTM pixels vs. the observed crop,
                        **restricted to the intersection of the two
                        footprints**: both grids are then scored on exactly the
                        same pixels, in the native UTM frame, so the fact that
                        HEALPix cells are rotated with respect to UTM ones
                        stops being something the reader has to take on trust;
  * ``floor``        -- the same projection applied to the *uncompressed*
                        tile. This is the HEALPix representation error: no
                        compression setting can go below it, and quoting
                        ``rmse_common`` without it would charge the carrier
                        grid's resampling error to the codec.

Note the deliberate model asymmetry, which the paper must state: tiles are
populated with Clough-Tocher interpolation (the paper's reference operator)
while the projection back uses the Gaussian observation operator. The tile is
therefore not the inverse of the operator used to predict, and the floor is
correspondingly higher than it would be with ``PSFResampler`` on both ends.
``--mode psf`` runs the consistent variant as a sensitivity check.
"""

from __future__ import annotations

import numpy as np

import paper_config as cfg

DEFAULT_FWHM_M = 12.5     # healpix-resample's calibrated effective B04 response


class CommonReference:
    """Forward operator from the tile's cells to its 637x637 UTM pixels."""

    def __init__(self, cell_ids, utm_lon, utm_lat, *, fwhm_m=DEFAULT_FWHM_M,
                 npt=None, mode="knn", erode=3, verbose=False):
        from healpix_resample import KNeighborsResampler, PSFResampler
        from healpix_resample.psf_geometry import fwhm_to_scale, recommend_npt

        self.shape = utm_lon.shape
        self.cell_ids = np.asarray(cell_ids, dtype=np.int64)
        self.base_id = int(self.cell_ids[0])
        scale_m = fwhm_to_scale(fwhm_m)
        if npt is None:
            npt = int(recommend_npt(scale_m, cfg.LEVEL)["npt"])
        self.npt = npt
        self.scale_m = scale_m
        self.mode = mode

        cls = PSFResampler if mode == "psf" else KNeighborsResampler
        self.op = cls(np.asarray(utm_lon, dtype=np.float64).ravel(),
                      np.asarray(utm_lat, dtype=np.float64).ravel(),
                      npt, cfg.LEVEL, sigma_m=scale_m,
                      out_cell_ids=self.cell_ids, verbose=verbose)

        kept = np.asarray(self.op.get_cell_ids(), dtype=np.int64).ravel()
        self.kept_index = kept - self.base_id
        if self.kept_index.min() < 0 or self.kept_index.max() >= self.cell_ids.size:
            raise RuntimeError("operator returned cells outside the requested block")
        self.kept_fraction = kept.size / float(self.cell_ids.size)

        # --- where can the two grids actually be compared? -------------------
        # The HEALPix tile is a square rotated by roughly 45 degrees in UTM;
        # the equal-area crop is axis-aligned. They enclose the same area but
        # not the same ground: each one's corners stick out of the other. A UTM
        # pixel in such a corner has no HEALPix cell behind it, and projecting
        # the tile there returns 0 -- comparing that against real reflectance
        # charges HEALPix for ground it never stored, and swamps the codec's
        # own error completely (it is what made rmse_common come out constant
        # across every sweep point).
        #
        # The intersection is computed by definition rather than inferred from
        # the operator: a crop pixel is inside the tile exactly when the
        # level-20 cell containing it belongs to the tile's child block. That
        # is one healpix_geo call, it cannot be thrown off by how the operator
        # normalises or thresholds its weights, and it is checkable by hand.
        import healpix_geo
        parent = self.base_id // cfg.N_CHILDREN
        ipix = healpix_geo.nested.lonlat_to_healpix(
            np.asarray(utm_lon, dtype=np.float64).ravel(),
            np.asarray(utm_lat, dtype=np.float64).ravel(),
            cfg.LEVEL, ellipsoid="WGS84").astype(np.int64)
        inside = (ipix // cfg.N_CHILDREN == parent).reshape(self.shape)
        self.raw_inside_fraction = float(inside.mean())
        if erode:
            from scipy.ndimage import binary_erosion
            # A pixel on the boundary sees a truncated kernel, so its
            # prediction is biased even though it is nominally inside.
            inside = binary_erosion(inside, iterations=int(erode), border_value=0)
        self.mask = inside
        self.mask_fraction = float(inside.mean())

        # --- the projection itself -------------------------------------------
        # `KNeighborsResampler.invert` multiplies a dense vector by a sparse CSR
        # tensor from the left. That dispatch silently returns zeros in some
        # torch builds, and a projection that is identically zero is invisible:
        # every rate-distortion point then reports the same "error", namely the
        # RMS of the observations themselves. We therefore do the multiplication
        # explicitly, with the operator's own coefficients, and check it.
        MT = self.op.MT
        coo = (MT.to_sparse_coo() if hasattr(MT, "to_sparse_coo") else MT).coalesce()
        idx = coo.indices().cpu().numpy()
        self._cell_of_link = idx[0].astype(np.int64)     # row of MT  = cell
        self._sample_of_link = idx[1].astype(np.int64)   # col of MT  = sample
        self._weight = coo.values().detach().cpu().numpy().astype(np.float64)
        self._n_samples = int(MT.shape[1])
        self._n_cells = int(MT.shape[0])
        if self._n_cells != kept.size:
            raise RuntimeError(f"operator has {self._n_cells} rows but reports "
                               f"{kept.size} cells")

        support = self._project(np.ones(kept.size, dtype=np.float64))
        self.support_stats = dict(
            min=float(np.nanmin(support)), median=float(np.nanmedian(support)),
            max=float(np.nanmax(support)),
            median_on_mask=float(np.nanmedian(support[self.mask]))
            if self.mask.any() else float("nan"))
        # Row-normalised weights: a fully supported sample must sum to one.
        med = self.support_stats["median_on_mask"]
        if not np.isfinite(med) or abs(med - 1.0) > 0.05:
            raise RuntimeError(
                "the HEALPix -> UTM projection does not preserve a constant "
                f"field (median support {med:.4f} on the comparable area, "
                f"expected 1.000). Operator statistics: {self.support_stats}")

        if self.mask_fraction < 0.2:
            raise RuntimeError(
                f"only {100*self.mask_fraction:.1f}% of the UTM crop lies inside "
                f"the HEALPix tile footprint (before erosion: "
                f"{100*self.raw_inside_fraction:.1f}%); operator support stats "
                f"{self.support_stats}")

    def _project(self, kept_values):
        """Sparse mat-vec MT^T @ v, written out rather than dispatched.

        Each non-zero of MT links one retained cell to one UTM sample with a
        row-normalised weight, so accumulating weight * value into the sample
        it belongs to is exactly the operator's forward prediction.
        """
        v = np.asarray(kept_values, dtype=np.float64).ravel()
        if v.size != self._n_cells:
            raise ValueError(f"expected {self._n_cells} retained cells, got {v.size}")
        out = np.bincount(self._sample_of_link,
                          weights=self._weight * v[self._cell_of_link],
                          minlength=self._n_samples)
        return out.reshape(self.shape)

    def predict(self, cell_values):
        """Project a full-block cell array onto the UTM crop grid."""
        vals = np.asarray(cell_values, dtype=np.float64).ravel()
        if vals.size != self.cell_ids.size:
            raise ValueError(f"expected {self.cell_ids.size} cells, got {vals.size}")
        return self._project(vals[self.kept_index])

    def rmse_on_mask(self, predicted, observed) -> float:
        """RMSE over the two footprints' intersection only."""
        p = np.asarray(predicted, dtype=np.float64)[self.mask]
        o = np.asarray(observed, dtype=np.float64)[self.mask]
        ok = np.isfinite(p) & np.isfinite(o)
        return float(np.sqrt(np.mean((p[ok] - o[ok]) ** 2))) if ok.any() else np.nan

    def describe(self) -> str:
        return (f"mode={self.mode} npt={self.npt} scale={self.scale_m:.2f} m "
                f"(FWHM {self.scale_m*1.1774:.2f} m) cells kept "
                f"{100*self.kept_fraction:.2f}%, footprint overlap "
                f"{100*self.raw_inside_fraction:.2f}% -> comparable area "
                f"{100*self.mask_fraction:.2f}% of the crop | support median "
                f"{self.support_stats.get('median_on_mask', float('nan')):.4f}")


def build_for_tile(npz, *, fwhm_m=DEFAULT_FWHM_M, mode="knn", erode=3,
                   verbose=False):
    """Build a CommonReference from a tile .npz written by build_tiles.py."""
    parent = int(npz["parent_id"])
    cell_ids = parent * cfg.N_CHILDREN + np.arange(cfg.N_CHILDREN, dtype=np.int64)
    return CommonReference(cell_ids, npz["utm_lon"], npz["utm_lat"],
                           fwhm_m=fwhm_m, mode=mode, erode=erode, verbose=verbose)
