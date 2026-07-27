"""
A wavelet transform computed **natively on HEALPix's NESTED (Morton)
ordering** -- no ``expand_tile``/``flatten_tile`` detour through a 2-D
raster, no :mod:`pywt`.

Motivation
----------
:func:`hcc.compression.grad_dwt.dwt_forward` works by first calling
:func:`hcc.compression.zorder.expand_tile` to de-interleave a HEALPix
NESTED-ordered 1-D array into a square 2-D raster (undoing the Morton bit
interleaving), then running a generic separable 2-D DWT (:mod:`pywt`) on
that raster. This is correct, but it is not aware that the array it is
given is a quadtree: it treats HEALPix's dyadic pixel hierarchy as an
implementation detail to be unpacked away before the "real" transform
starts.

But HEALPix NESTED indexing already *is* a quadtree by construction: pixel
``i`` at some resolution has *exactly* its 4 children at indices ``4*i,
4*i+1, 4*i+2, 4*i+3`` at the next-finer resolution (this is exactly what
Morton/Z-order numbering means). :func:`hcc.compression.zorder._decode_zorder`
shows the precise bit convention this codebase uses: within one group of 4
consecutive NESTED indices (local code ``k in {0,1,2,3}``), bit 0 of ``k``
selects the x-offset and bit 1 selects the y-offset, i.e.
``k=0 -> (dx,dy)=(0,0)``, ``k=1 -> (1,0)``, ``k=2 -> (0,1)``, ``k=3 -> (1,1)``.

A single-level 2-D Haar transform on that 2x2 block is *entirely local* to
one group of 4 -- it needs no neighbouring information -- so it can be
computed directly from 4 consecutive array elements, with no 2-D expansion
at all. And because NESTED numbering is self-similar across scales (the
sequence of "LL" approximations produced by one grouping pass is itself a
valid NESTED-ordered array at the next-coarser HEALPix resolution), the
same grouping step can simply be repeated to get a full multi-level
pyramid -- each level corresponding *exactly* to one HEALPix resolution
level (nside -> nside/2), by construction, not by later reshaping.

Scope of this first version (Haar)
------------------------------------
The :func:`qt_dwt_forward`/:func:`qt_dwt_inverse` pair implements the
**Haar** case (2-tap filter): the natural, exact fit for quadtree
grouping, since a 2x2 Haar butterfly touches precisely the 4 children of
one parent and nothing else -- no cross-group neighbour lookups needed.

Longer filters: CDF 5/3, via lifting with neighbouring quadtree groups
------------------------------------------------------------------------
:func:`qt_dwt53_forward`/:func:`qt_dwt53_inverse` extend this to the CDF
5/3 filter (the JPEG2000-lossless wavelet), which needs *neighbouring*
quadtree groups, not just one parent's 4 children. This is done with a
lifting scheme computed **directly on the NESTED array**, using
:func:`hcc.compression.zorder.encode_zorder`/``_decode_zorder`` to look up
a pixel's left/right (or up/down) neighbour by index arithmetic -- with
periodic wraparound at each level's own resolution -- instead of forming a
2-D raster. Concretely, for one level of resolution ``2**delta`` per side:

1. Decode every NESTED index in the level to (x, y) via ``_decode_zorder``.
   Because of the bit-interleaving convention, the low bit of the NESTED
   index directly equals ``x & 1``, and the next bit equals ``y & 1`` --
   i.e. the same 4-way grouping used by the Haar step (``arr[0::4]``,
   ``arr[1::4]``, ...) already sorts pixels by (x parity, y parity).
2. Run a standard CDF 5/3 lifting pass *in place* along x (predict odd-x
   pixels from their even-x neighbours, then update even-x pixels from
   the new odd-x details), using neighbour lookups of the form
   ``encode_zorder((x + dx) % side, y, delta)``.
3. Repeat along y on the result.
4. Because of the parity alignment from step 1, subsampling to extract
   LL/LH/HL/HL is then just ``arr[0::4], arr[2::4], arr[1::4], arr[3::4]``
   -- exactly the same slicing the Haar step uses, no separate gather.

The periodic wraparound at each level's own side length mirrors the
``mode="periodization"`` convention already used by
:func:`hcc.compression.grad_dwt.dwt_forward` for the generic
:mod:`pywt`-based path, so results are directly comparable.

Validated properties (see ``hcc/tests/test_quadtree_wavelet.py``): exact
perfect reconstruction (float64, multi-level, various tile sizes); on a
realistic 1024x1024 test image, quantized rate-distortion (ratio and
RMSE) matching ``bior4.4`` via the existing ``expand_tile`` + :mod:`pywt`
path almost exactly -- markedly better than the native Haar step, which is
expected given CDF 5/3's longer support. The current implementation is,
however, substantially *slower* than the generic pywt path (the per-level
neighbour-gather via ``encode_zorder`` is pure NumPy fancy-indexing, not
vectorized C convolution) -- a known, unoptimized cost of doing the
lookup natively rather than via a 2-D raster; see module-level TODO if
revisiting for production use.

Base-pixel boundary caveat
----------------------------
Both the Haar and the CDF 5/3 native transforms operate on a single
HEALPix base-pixel's quadtree (or an equivalent square tile). Periodic
wraparound is applied at each level's own tile boundary -- correct for a
single tile treated in isolation (matching how the existing residual/DWT
codecs already process per-chunk tiles), but the pixel-adjacency quirks
*between* different HEALPix base pixels (which are not simple
index-arithmetic neighbours) are out of scope here, exactly as for the
existing ``expand_tile``-based path.

Coefficient layout
-------------------
Unlike :func:`hcc.compression.grad_dwt.dwt_forward` (which returns one flat
2-D array + a separate ``slices`` structure describing how to carve it into
subbands), this module returns subbands as **plain 1-D arrays directly**:
there is no shared raster to slice into, since nothing is ever expanded to
2-D. See :func:`qt_dwt_forward` / :func:`qt_dwt_inverse` (Haar) and
:func:`qt_dwt53_forward` / :func:`qt_dwt53_inverse` (CDF 5/3).

Public API
----------
qt_dwt_forward : Haar forward transform, NESTED array -> (LL, per-level details).
qt_dwt_inverse : Exact inverse of :func:`qt_dwt_forward`.
qt_dwt53_forward : CDF 5/3 forward transform, NESTED array -> (LL, per-level details).
qt_dwt53_inverse : Exact inverse of :func:`qt_dwt53_forward`.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from hcc.compression.zorder import _decode_zorder, encode_zorder

__all__ = [
    "qt_dwt_forward",
    "qt_dwt_inverse",
    "qt_dwt53_forward",
    "qt_dwt53_inverse",
]


# ─────────────────────────────────────────────────────────────────
#  Single-level 2x2 Haar butterfly, applied to groups of 4 consecutive
#  NESTED-ordered elements
# ─────────────────────────────────────────────────────────────────


def _haar_analysis_step(x: npt.NDArray) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray]:
    """One level of quadtree-grouped 2-D Haar analysis.

    Groups *x* (length divisible by 4) into consecutive runs of 4 --
    exactly the 4 children of one parent pixel in NESTED order -- and
    applies an orthonormal 2x2 Haar transform to each group.

    Using the local-code -> (dx, dy) convention from
    :func:`hcc.compression.zorder._decode_zorder`
    (``k=0,1,2,3 -> (dx,dy) = (0,0),(1,0),(0,1),(1,1)``), writing
    ``a0..a3`` for the 4 children in NESTED order, the separable 2-D Haar
    transform (1-D Haar along x, then along y) works out to::

        LL = ( a0 + a1 + a2 + a3) / 2
        LH = ( a0 + a1 - a2 - a3) / 2   -- low-x,  high-y detail
        HL = ( a0 - a1 + a2 - a3) / 2   -- high-x, low-y  detail
        HH = ( a0 - a1 - a2 + a3) / 2   -- high-x, high-y detail

    (an orthonormal, energy-preserving 4x4 Hadamard-like transform; see
    :func:`_haar_synthesis_step` for the matching, exact inverse).

    Returns
    -------
    ll, lh, hl, hh : numpy.ndarray, each of length ``len(x) // 4``
    """
    a0 = x[0::4]
    a1 = x[1::4]
    a2 = x[2::4]
    a3 = x[3::4]
    ll = (a0 + a1 + a2 + a3) * 0.5
    lh = (a0 + a1 - a2 - a3) * 0.5
    hl = (a0 - a1 + a2 - a3) * 0.5
    hh = (a0 - a1 - a2 + a3) * 0.5
    return ll, lh, hl, hh


def _haar_synthesis_step(
    ll: npt.NDArray, lh: npt.NDArray, hl: npt.NDArray, hh: npt.NDArray
) -> npt.NDArray:
    """Exact inverse of :func:`_haar_analysis_step`.

    The forward step's 4x4 matrix ``M`` (mapping ``[a0,a1,a2,a3]`` to
    ``2*[LL,LH,HL,HH]``) is symmetric with ``M @ M == 4*I`` (each row
    orthogonal to every other, each with squared norm 4), so
    ``M^-1 == M / 4`` and the inverse has the same coefficients as the
    forward step, just re-summed::

        a0 = ( LL + LH + HL + HH) / 2
        a1 = ( LL + LH - HL - HH) / 2
        a2 = ( LL - LH + HL - HH) / 2
        a3 = ( LL - LH - HL + HH) / 2

    Returns
    -------
    numpy.ndarray of length ``4 * len(ll)``, NESTED order (interleaved
    a0, a1, a2, a3 per group, matching :func:`_haar_analysis_step`'s
    input layout).
    """
    a0 = (ll + lh + hl + hh) * 0.5
    a1 = (ll + lh - hl - hh) * 0.5
    a2 = (ll - lh + hl - hh) * 0.5
    a3 = (ll - lh - hl + hh) * 0.5

    out = np.empty(ll.size * 4, dtype=ll.dtype)
    out[0::4] = a0
    out[1::4] = a1
    out[2::4] = a2
    out[3::4] = a3
    return out


# ─────────────────────────────────────────────────────────────────
#  Multi-level pyramid
# ─────────────────────────────────────────────────────────────────


def qt_dwt_forward(
    data: npt.NDArray[np.floating], levels: int
) -> tuple[npt.NDArray[np.float64], list[tuple[npt.NDArray, npt.NDArray, npt.NDArray]]]:
    """Multi-level quadtree-native Haar analysis of a HEALPix NESTED tile.

    Each level exactly halves the tile's linear HEALPix resolution
    (quarters its pixel count) -- by construction, not by later
    resampling -- since it is computed by directly regrouping the NESTED
    array 4 elements at a time, recursively.

    Parameters
    ----------
    data : numpy.ndarray of shape (n,), n a power of 4
        HEALPix NESTED-ordered 1-D tile (the same layout
        :func:`hcc.compression.grad_dwt.compress`/``compress_residuals``'s
        callers work with before calling
        :func:`hcc.compression.zorder.expand_tile` -- this module skips
        that step entirely).
    levels : int
        Number of decomposition levels. Must satisfy
        ``4**levels <= data.size``.

    Returns
    -------
    ll : numpy.ndarray of shape (data.size // 4**levels,), dtype float64
        Coarsest-level approximation -- itself a valid NESTED-ordered
        array at the corresponding coarser HEALPix resolution.
    details : list of (lh, hl, hh) tuples, length *levels*
        One entry per level, **coarsest first** (index 0 is the level
        closest to *ll*, matching the ordering convention of
        :func:`hcc.compression.grad_dwt._subband_slices`), each a
        3-tuple of 1-D arrays of matching length at that level's
        resolution.

    See Also
    --------
    qt_dwt_inverse : Exact inverse of this function.
    """
    n = data.shape[-1]
    if 4**levels > n:
        raise ValueError(
            f"levels={levels} requires 4**levels={4**levels} <= data.size={n}."
        )

    current = np.asarray(data, dtype=np.float64)
    details_finest_first: list[tuple[npt.NDArray, npt.NDArray, npt.NDArray]] = []
    for _ in range(levels):
        ll, lh, hl, hh = _haar_analysis_step(current)
        details_finest_first.append((lh, hl, hh))
        current = ll

    return current, list(reversed(details_finest_first))


def qt_dwt_inverse(
    ll: npt.NDArray[np.floating],
    details: list[tuple[npt.NDArray, npt.NDArray, npt.NDArray]],
) -> npt.NDArray[np.float64]:
    """Exact inverse of :func:`qt_dwt_forward`.

    Parameters
    ----------
    ll : numpy.ndarray
        Coarsest-level approximation, as returned by
        :func:`qt_dwt_forward`.
    details : list of (lh, hl, hh) tuples
        Coarsest-first, as returned by :func:`qt_dwt_forward`.

    Returns
    -------
    numpy.ndarray, dtype float64
        Reconstructed HEALPix NESTED-ordered 1-D tile.
    """
    current = np.asarray(ll, dtype=np.float64)
    for lh, hl, hh in details:  # coarsest-first, matching qt_dwt_forward's return order
        current = _haar_synthesis_step(current, lh, hl, hh)
    return current


# ─────────────────────────────────────────────────────────────────
#  CDF 5/3 lifting, computed natively on the NESTED array via
#  encode_zorder-based neighbour lookups (no 2-D raster)
# ─────────────────────────────────────────────────────────────────


def _shift(
    values: npt.NDArray,
    x: npt.NDArray[np.uint64],
    y: npt.NDArray[np.uint64],
    side: int,
    delta: int,
    axis: str,
    step: int,
) -> npt.NDArray:
    """Gather ``values`` at the neighbour ``step`` pixels away along
    ``axis`` ("x" or "y"), with periodic wraparound at this level's own
    side length -- i.e. ``values[i]`` is replaced by the value that lives
    at ``(x[i] + step, y[i])`` (axis="x") or ``(x[i], y[i] + step)``
    (axis="y"), looked up via :func:`~hcc.compression.zorder.encode_zorder`.
    """
    if axis == "x":
        new_x = (x.astype(np.int64) + step) % side
        idx = encode_zorder(new_x.astype(np.uint64), y, delta)
    else:
        new_y = (y.astype(np.int64) + step) % side
        idx = encode_zorder(x, new_y.astype(np.uint64), delta)
    return values[idx]


def _lift_axis(
    values: npt.NDArray,
    coords_axis: npt.NDArray[np.uint64],
    x: npt.NDArray[np.uint64],
    y: npt.NDArray[np.uint64],
    side: int,
    delta: int,
    axis: str,
) -> npt.NDArray:
    """One in-place CDF 5/3 analysis lifting pass along *axis*.

    Odd-``coords_axis`` positions are overwritten with their predict
    residual (detail); even-``coords_axis`` positions are overwritten with
    their update result (smoothed approximation) -- both computed from
    their immediate +-1 neighbours along *axis* (periodic wraparound via
    :func:`_shift`). This is the standard CDF 5/3 (JPEG2000-lossless)
    lifting scheme::

        d[odd]  = o[odd]  - 0.5 * (e[odd-1]  + e[odd+1])     (predict)
        s[even] = e[even] + 0.25 * (d[even-1] + d[even+1])   (update)

    where ``e``/``o`` denote the even-/odd-indexed samples of the original
    signal along *axis*; here both are read from (and written back into)
    the same array, since even and odd positions never overwrite each
    other's *inputs* within one pass (the update step reads the *already
    predicted* detail values at the odd neighbours, exactly as the
    classic in-place lifting scheme requires).
    """
    is_odd = (coords_axis & np.uint64(1)) == np.uint64(1)

    left = _shift(values, x, y, side, delta, axis, -1)
    right = _shift(values, x, y, side, delta, axis, +1)
    predicted = values - 0.5 * (left + right)
    d_full = np.where(is_odd, predicted, values)

    left_d = _shift(d_full, x, y, side, delta, axis, -1)
    right_d = _shift(d_full, x, y, side, delta, axis, +1)
    updated = values + 0.25 * (left_d + right_d)

    return np.where(is_odd, d_full, updated)


def _unlift_axis(
    values: npt.NDArray,
    coords_axis: npt.NDArray[np.uint64],
    x: npt.NDArray[np.uint64],
    y: npt.NDArray[np.uint64],
    side: int,
    delta: int,
    axis: str,
) -> npt.NDArray:
    """Exact inverse of :func:`_lift_axis` (reversed lifting steps)."""
    is_odd = (coords_axis & np.uint64(1)) == np.uint64(1)

    left_d = _shift(values, x, y, side, delta, axis, -1)
    right_d = _shift(values, x, y, side, delta, axis, +1)
    e_recovered = values - 0.25 * (left_d + right_d)
    tmp = np.where(is_odd, values, e_recovered)

    left_e = _shift(tmp, x, y, side, delta, axis, -1)
    right_e = _shift(tmp, x, y, side, delta, axis, +1)
    o_recovered = tmp + 0.5 * (left_e + right_e)

    return np.where(is_odd, o_recovered, tmp)


def _cdf53_analysis_step(
    current: npt.NDArray[np.float64], delta: int
) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray]:
    """One level of quadtree-native CDF 5/3 analysis (separable: x then y).

    Because the NESTED bit-interleaving convention puts ``x``'s parity in
    the low index bit and ``y``'s parity in the next bit (see module
    docstring), after lifting along x then y in place, the four subbands
    are recovered by exactly the same grouping the Haar step uses::

        ll = lifted[0::4]   (x even, y even -- approx/approx)
        lh = lifted[2::4]   (x even, y odd  -- approx-x, detail-y)
        hl = lifted[1::4]   (x odd,  y even -- detail-x, approx-y)
        hh = lifted[3::4]   (x odd,  y odd  -- detail/detail)
    """
    side = 2**delta
    x, y = _decode_zorder(delta)
    lifted_x = _lift_axis(current, x, x, y, side, delta, "x")
    lifted_xy = _lift_axis(lifted_x, y, x, y, side, delta, "y")
    return lifted_xy[0::4], lifted_xy[2::4], lifted_xy[1::4], lifted_xy[3::4]


def _cdf53_synthesis_step(
    ll: npt.NDArray, lh: npt.NDArray, hl: npt.NDArray, hh: npt.NDArray, delta: int
) -> npt.NDArray:
    """Exact inverse of :func:`_cdf53_analysis_step`.

    ``delta`` is the level of the *output* (upsampled) resolution, i.e.
    ``4 * ll.size == 4**delta``.
    """
    side = 2**delta
    x, y = _decode_zorder(delta)

    full = np.empty(ll.size * 4, dtype=np.float64)
    full[0::4] = ll
    full[2::4] = lh
    full[1::4] = hl
    full[3::4] = hh

    unlifted_y = _unlift_axis(full, y, x, y, side, delta, "y")
    return _unlift_axis(unlifted_y, x, x, y, side, delta, "x")


def qt_dwt53_forward(
    data: npt.NDArray[np.floating], levels: int
) -> tuple[npt.NDArray[np.float64], list[tuple[npt.NDArray, npt.NDArray, npt.NDArray]]]:
    """Multi-level quadtree-native CDF 5/3 analysis of a HEALPix NESTED tile.

    Longer-filter counterpart of :func:`qt_dwt_forward`: same recursive
    per-level structure (each level exactly halves the tile's linear
    HEALPix resolution), but each level's analysis step is the CDF 5/3
    lifting pass of :func:`_cdf53_analysis_step` instead of the local 2x2
    Haar butterfly -- i.e. it looks beyond one parent's 4 children to
    neighbouring quadtree groups, via periodic-wraparound neighbour
    lookups at each level's own resolution.

    Parameters
    ----------
    data : numpy.ndarray of shape (n,), n a power of 4
        HEALPix NESTED-ordered 1-D tile.
    levels : int
        Number of decomposition levels. Must satisfy
        ``4**levels <= data.size``.

    Returns
    -------
    ll : numpy.ndarray of shape (data.size // 4**levels,), dtype float64
        Coarsest-level approximation.
    details : list of (lh, hl, hh) tuples, length *levels*, coarsest first
        (matching :func:`qt_dwt_forward`'s convention).

    See Also
    --------
    qt_dwt53_inverse : Exact inverse of this function.
    qt_dwt_forward : The Haar (2-tap) counterpart.
    """
    n = data.shape[-1]
    if 4**levels > n:
        raise ValueError(
            f"levels={levels} requires 4**levels={4**levels} <= data.size={n}."
        )
    delta_total = round(np.log(n) / np.log(4))
    if 4**delta_total != n:
        raise ValueError(f"data.size={n} must be a power of 4.")

    current = np.asarray(data, dtype=np.float64)
    details_finest_first: list[tuple[npt.NDArray, npt.NDArray, npt.NDArray]] = []
    delta = delta_total
    for _ in range(levels):
        ll, lh, hl, hh = _cdf53_analysis_step(current, delta)
        details_finest_first.append((lh, hl, hh))
        current = ll
        delta -= 1

    return current, list(reversed(details_finest_first))


def qt_dwt53_inverse(
    ll: npt.NDArray[np.floating],
    details: list[tuple[npt.NDArray, npt.NDArray, npt.NDArray]],
) -> npt.NDArray[np.float64]:
    """Exact inverse of :func:`qt_dwt53_forward`.

    Parameters
    ----------
    ll : numpy.ndarray
        Coarsest-level approximation, as returned by
        :func:`qt_dwt53_forward`.
    details : list of (lh, hl, hh) tuples
        Coarsest-first, as returned by :func:`qt_dwt53_forward`.

    Returns
    -------
    numpy.ndarray, dtype float64
        Reconstructed HEALPix NESTED-ordered 1-D tile.
    """
    current = np.asarray(ll, dtype=np.float64)
    total_size = ll.size * 4 ** len(details)
    delta = round(np.log(total_size) / np.log(4)) - len(details)
    for lh, hl, hh in details:  # coarsest-first, matching qt_dwt53_forward's return order
        delta += 1
        current = _cdf53_synthesis_step(current, lh, hl, hh, delta)
    return current
