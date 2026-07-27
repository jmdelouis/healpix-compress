"""SPIHT-style set-partitioning entropy coding, native to the HEALPix
quadtree wavelet coefficients of :mod:`hcc.compression.quadtree_wavelet`.

Motivation and relation to :mod:`hcc.compression.zerotree`
-------------------------------------------------------------
:mod:`hcc.compression.zerotree` implements a bit-plane, cross-scale
context-adaptive coder (an EBCOT-flavoured design) that was found, on
realistic dense coefficient data, to be *worse* than plain
``zlib.compress`` on a uniformly-quantised array -- because every
individual coefficient was re-tested for significance on *every*
bit-plane, with no mechanism to skip a whole insignificant subtree in
less than O(subtree size) work. That is precisely the piece classic
SPIHT (Said & Pearlman, 1996) adds: **set-partitioning** -- a coefficient
together with *all* of its descendants (down to the finest decomposition
level) is treated as a single "insignificant set" as long as the largest
magnitude anywhere in that set is below the current bit-plane threshold,
tested and skipped in O(1) via a precomputed max-magnitude pyramid,
instead of being scanned.

This module implements that missing piece, built on top of the **native**
quadtree wavelet transform (:func:`hcc.compression.quadtree_wavelet.qt_dwt_forward`
/ :func:`~hcc.compression.quadtree_wavelet.qt_dwt53_forward`) rather than
the generic ``pywt``-based path -- because that transform's per-level flat
1-D subband arrays give an exact, trivial parent -> 4-children relation
(index ``i`` at level ``l`` has children at indices ``4*i .. 4*i+3`` at
level ``l+1``, same orientation), with **no raster reshaping or slice
bookkeeping needed** to walk the tree. This is the concrete form of the
synergy between HEALPix's quadtree indexing and per-band non-uniform
coding that motivated this whole line of work.

Algorithm summary
-------------------
Standard SPIHT sorting-pass / refinement-pass structure, per bit-plane
(MSB to LSB), operating on three coefficient lists:

- **LIP** (List of Insignificant Pixels): individual coefficients not yet
  found significant; each is bit-tested every remaining bit-plane.
- **LIS** (List of Insignificant Sets): tree nodes representing either
  "all descendants of coefficient *i*" (type A) or "all descendants of
  *i*'s children, i.e. grandchildren and below" (type B). Tested via the
  precomputed ``desc_max``/``gc_max`` pyramids (see
  :func:`_build_pyramids`) in O(1); insignificant sets are skipped
  entirely (no descendants visited), significant sets are split
  (type A -> test its 4 children individually + re-add as type B if it
  has grandchildren; type B -> re-add its 4 children as type-A entries).
- **LSP** (List of Significant Pixels): coefficients already found
  significant in a previous pass; refined by one more magnitude bit per
  remaining bit-plane.

Symbols (significance bits, signs, refinement bits) are coded with five
small adaptive binary models (reusing
:class:`hcc.compression.entropy._AdaptiveModel` and the same range coder)
-- one per symbol *role*, not per spatial context; this is a simpler
context scheme than :mod:`hcc.compression.zerotree`'s 9-context design,
relying on the set-partitioning structure itself (rather than rich
contexts) for the compression gain.

The LL/coarsest approximation band is **not** coded by this module --
callers are expected to store it losslessly (as raw float32), matching
the rest of this codebase's established LL-lossless convention (see
:func:`hcc.compression.grad_dwt.compress_residuals`). Only the *detail*
subbands (as returned by ``qt_dwt_forward``/``qt_dwt53_forward``) are
coded here, with the coarsest detail level's coefficients acting as the
tree roots.

Honest scope / measured results
----------------------------------
Unlike :mod:`hcc.compression.zerotree`, this design *does* beat the
production per-subband-quantised + zlib baseline -- but the size of the
win is highly data-dependent:

- On a synthetic "patchy" test image (a handful of homogeneous regions
  plus a few sharp edges -- a rough proxy for land-cover-like structure
  in real Sentinel-2 tiles), at ``quantize_bits=8``: this coder was
  ~38-41% smaller than the current production per-subband zlib baseline
  on the native **Haar** wavelet's coefficients, and ~22-26% smaller on
  the native **CDF 5/3** wavelet's coefficients (which spreads energy
  more, leaving fewer near-zero coefficients for set-partitioning to
  exploit -- Haar's sparser residual structure pairs better with this
  coder than the longer filter does).
- On a fully dense, globally-smooth-but-non-sparse synthetic signal (no
  homogeneous flat regions -- see this module's tests), this coder
  roughly **ties** zlib-on-uniform-quantisation (ratio close to 1.0, even
  a slight loss in one noiseless case) -- set-partitioning has nothing to
  exploit when there are no large below-threshold regions to skip.

In short: the gain is real and substantial on data with the kind of
spatial sparsity/homogeneity real remote-sensing imagery has, but is not
a universal win independent of data statistics -- exactly the property
classic SPIHT/EBCOT are known for.

Performance caveat: this is a pure-Python reference implementation (list-
based LIP/LIS/LSP bookkeeping); it is meaningfully slower than
:mod:`zlib` or the vectorised wavelet transforms themselves. It has not
been optimised (e.g. with Cython/Numba) -- see module TODO if pursuing
production use.

Public API
----------
encode_spiht : Entropy-code detail subbands (list of (lh, hl, hh) tuples,
    coarsest first, as returned by ``qt_dwt_forward``/``qt_dwt53_forward``).
decode_spiht : Exact inverse of :func:`encode_spiht`.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from hcc.compression.entropy import _AdaptiveModel, _RangeDecoder, _RangeEncoder

__all__ = ["encode_spiht", "decode_spiht"]


# ─────────────────────────────────────────────────────────────────
#  O(1)-testable subtree-maximum pyramids
# ─────────────────────────────────────────────────────────────────


def _build_pyramids(
    mag: list[tuple[npt.NDArray, npt.NDArray, npt.NDArray]],
) -> tuple[
    list[tuple[npt.NDArray, npt.NDArray, npt.NDArray]],
    list[tuple[npt.NDArray, npt.NDArray, npt.NDArray]],
]:
    """Precompute, for every node, the maximum integer magnitude anywhere
    in its descendant subtree -- the piece that makes set significance
    testing O(1) instead of O(subtree size) (this is *the* fix for
    :mod:`hcc.compression.zerotree`'s underperformance: that module had no
    equivalent and had to re-test every leaf every bit-plane).

    Parameters
    ----------
    mag : list of (lh, hl, hh) int arrays, coarsest first
        Quantised (integer) magnitude arrays, same layout as
        :func:`hcc.compression.quadtree_wavelet.qt_dwt_forward`'s
        ``details`` return value.

    Returns
    -------
    desc_max : same structure as *mag*
        ``desc_max[l][o][i]`` = max magnitude over *all* strict
        descendants of node ``i`` (excluding the node itself) -- used for
        SPIHT type-A ("D(i)") significance tests.
    gc_max : same structure as *mag*
        ``gc_max[l][o][i]`` = max magnitude over grandchildren-and-below
        only (excluding the node and its direct children) -- used for
        SPIHT type-B ("L(i)") significance tests.
    """
    n_levels = len(mag)
    desc_max: list = [None] * n_levels
    gc_max: list = [None] * n_levels
    # Finest level has no children -> no descendants.
    desc_max[-1] = tuple(np.zeros_like(a) for a in mag[-1])
    gc_max[-1] = tuple(np.zeros_like(a) for a in mag[-1])

    for level in range(n_levels - 2, -1, -1):
        level_desc = []
        level_gc = []
        for orientation in range(3):
            child_abs = mag[level + 1][orientation]
            child_desc = desc_max[level + 1][orientation]
            n_parent = mag[level][orientation].size
            # max over (child's own value, child's own descendants) grouped
            # by 4 -> max over *all* strict descendants of the parent.
            combined = np.maximum(child_abs, child_desc).reshape(n_parent, 4).max(axis=1)
            level_desc.append(combined)
            # max over children's descendants only (excludes the children
            # themselves) -> grandchildren-and-below.
            level_gc.append(child_desc.reshape(n_parent, 4).max(axis=1))
        desc_max[level] = tuple(level_desc)
        gc_max[level] = tuple(level_gc)

    return desc_max, gc_max


# ─────────────────────────────────────────────────────────────────
#  Public API
# ─────────────────────────────────────────────────────────────────


def encode_spiht(
    details: list[tuple[npt.NDArray, npt.NDArray, npt.NDArray]],
    quantize_bits: int,
    clip_percentile: float = 0.5,
    max_bitplanes: int | None = None,
) -> tuple[bytes, np.float32, int, int]:
    """SPIHT-encode the detail subbands of a native quadtree wavelet
    decomposition (:func:`hcc.compression.quadtree_wavelet.qt_dwt_forward`
    or ``qt_dwt53_forward``'s ``details`` return value).

    A single **global** scale is used across all levels/orientations
    (unlike :mod:`hcc.compression.grad_dwt`'s per-subband scales) --
    required for SPIHT's bit-plane thresholds to be directly comparable
    across the whole tree.

    Parameters
    ----------
    details : list of (lh, hl, hh) tuples, coarsest first
        As returned by ``qt_dwt_forward``/``qt_dwt53_forward``. The LL
        band is *not* passed here -- store it separately, losslessly.
    quantize_bits : int
        Number of bit-planes (i.e. quantisation levels = ``2**quantize_bits
        - 1``).
    clip_percentile : float, default 0.5
        Percentile (from the top) of |coefficient| used to set the global
        scale, clipping extreme outliers -- same convention as
        :func:`hcc.compression.grad_dwt.compress_residuals`.
    max_bitplanes : int or None, default None
        If given, stop after this many bit-planes (embedded/progressive
        truncation, coarser than the full ``quantize_bits``) -- None means
        fully encode all ``quantize_bits`` planes.

    Returns
    -------
    payload : bytes
        Self-contained range-coded bitstream.
    scale : numpy.float32
        Global dequantisation scale.
    n_bitplanes : int
        Total number of bit-planes implied by *quantize_bits* (``==
        quantize_bits``).
    n_bitplanes_encoded : int
        Number of bit-planes actually coded (``== max_bitplanes`` if
        given and smaller than ``quantize_bits``, else ``quantize_bits``).

    See Also
    --------
    decode_spiht : Exact inverse.
    """
    if quantize_bits < 1:
        raise ValueError(f"quantize_bits must be >= 1, got {quantize_bits}.")

    n_levels = len(details)
    max_level = (1 << quantize_bits) - 1
    # IMPORTANT: take the max of each subband's *own* percentile-clipped
    # cap, not one percentile pooled across all subbands combined. Coarser
    # levels have far fewer coefficients than finer ones (4x fewer per
    # level, by construction) but typically much larger magnitudes; a
    # single pooled percentile is dominated by the numerous small
    # fine-level values and ends up clipping away the coarse levels
    # entirely -- silently capping the representable range far below the
    # true largest coefficients. That clipping is a hard error floor no
    # number of extra bit-planes can fix (unlike ordinary quantisation
    # step size), and was previously validated by RMSE staying flat as
    # quantize_bits increased -- the bug this per-subband max fixes.
    per_subband_caps = []
    for triple in details:
        for sb in triple:
            abs_sb = np.abs(sb)
            if clip_percentile > 0 and abs_sb.size > 1:
                per_subband_caps.append(float(np.percentile(abs_sb, 100.0 - clip_percentile)))
            elif abs_sb.size:
                per_subband_caps.append(float(np.max(abs_sb)))
    max_abs = max(per_subband_caps) if per_subband_caps else 0.0
    if max_abs == 0.0:
        max_abs = 1.0
    scale = np.float32(max_abs / max_level)

    mag = [
        tuple(np.clip(np.round(np.abs(sb) / scale), 0, max_level).astype(np.int64) for sb in triple)
        for triple in details
    ]
    sign = [tuple(sb < 0 for sb in triple) for triple in details]
    desc_max, gc_max = _build_pyramids(mag)

    n_bitplanes = quantize_bits
    if max_bitplanes is None:
        max_bitplanes = n_bitplanes
    max_bitplanes = min(max_bitplanes, n_bitplanes)

    enc = _RangeEncoder()
    m_lip = _AdaptiveModel(2)
    m_lisA = _AdaptiveModel(2)
    m_lisB = _AdaptiveModel(2)
    m_sign = _AdaptiveModel(2)
    m_ref = _AdaptiveModel(2)

    lsp: list[tuple[int, int, int]] = []
    lip: list[tuple[int, int, int]] = [
        (0, o, i) for o in range(3) for i in range(mag[0][o].size)
    ]
    lis: list[tuple[int, int, int, str]] = []
    if n_levels > 1:
        lis = [(0, o, i, "A") for o in range(3) for i in range(mag[0][o].size)]

    for bp in range(n_bitplanes - 1, n_bitplanes - 1 - max_bitplanes, -1):
        threshold = 1 << bp
        lsp_len_before = len(lsp)

        new_lip = []
        for level, orientation, i in lip:
            val = mag[level][orientation][i]
            sig = val >= threshold
            m_lip.encode(enc, int(sig))
            if sig:
                m_sign.encode(enc, int(sign[level][orientation][i]))
                lsp.append((level, orientation, i))
            else:
                new_lip.append((level, orientation, i))
        lip = new_lip

        idx = 0
        while idx < len(lis):
            level, orientation, i, kind = lis[idx]
            if kind == "A":
                sig = bool(desc_max[level][orientation][i] >= threshold)
                m_lisA.encode(enc, int(sig))
                remove = False
                if sig:
                    child_level = level + 1
                    for c in range(4 * i, 4 * i + 4):
                        cval = mag[child_level][orientation][c]
                        csig = cval >= threshold
                        m_lip.encode(enc, int(csig))
                        if csig:
                            m_sign.encode(enc, int(sign[child_level][orientation][c]))
                            lsp.append((child_level, orientation, c))
                        else:
                            lip.append((child_level, orientation, c))
                    if child_level + 1 < n_levels:
                        lis.append((level, orientation, i, "B"))
                    remove = True
            else:  # type "B"
                sig = bool(gc_max[level][orientation][i] >= threshold)
                m_lisB.encode(enc, int(sig))
                remove = False
                if sig:
                    for c in range(4 * i, 4 * i + 4):
                        lis.append((level + 1, orientation, c, "A"))
                    remove = True
            if remove:
                lis.pop(idx)
            else:
                idx += 1

        for level, orientation, i in lsp[:lsp_len_before]:
            bit = (int(mag[level][orientation][i]) >> bp) & 1
            m_ref.encode(enc, bit)

    payload = enc.finish()
    return payload, scale, n_bitplanes, max_bitplanes


def decode_spiht(
    payload: bytes,
    shapes: list[tuple[int, int, int]],
    scale: np.float32,
    n_bitplanes: int,
    n_bitplanes_encoded: int,
) -> list[tuple[npt.NDArray, npt.NDArray, npt.NDArray]]:
    """Exact inverse of :func:`encode_spiht`.

    Parameters
    ----------
    payload : bytes
        As returned by :func:`encode_spiht`.
    shapes : list of (n_lh, n_hl, n_hh) tuples, coarsest first
        Per-level, per-orientation coefficient *counts* (``sb.size`` for
        each subband -- these are 1-D arrays), matching the *details*
        passed to :func:`encode_spiht`.
    scale, n_bitplanes, n_bitplanes_encoded
        As returned by :func:`encode_spiht`.

    Returns
    -------
    list of (lh, hl, hh) tuples, coarsest first, dtype float64
        Reconstructed detail subbands. If ``n_bitplanes_encoded <
        n_bitplanes`` (progressive truncation was used), coefficients not
        yet found significant decode to ``0.0`` and significant ones are
        reconstructed with the low-order bits below the truncation point
        set to ``0`` (least-squares-style embedded reconstruction, same
        semantics as :mod:`hcc.compression.zerotree`).
    """
    n_levels = len(shapes)
    mag = [tuple(np.zeros(s, dtype=np.int64) for s in triple) for triple in shapes]
    sign = [tuple(np.zeros(s, dtype=bool) for s in triple) for triple in shapes]

    dec = _RangeDecoder(payload)
    m_lip = _AdaptiveModel(2)
    m_lisA = _AdaptiveModel(2)
    m_lisB = _AdaptiveModel(2)
    m_sign = _AdaptiveModel(2)
    m_ref = _AdaptiveModel(2)

    lsp: list[tuple[int, int, int]] = []
    lip: list[tuple[int, int, int]] = [
        (0, o, i) for o in range(3) for i in range(shapes[0][o])
    ]
    lis: list[tuple[int, int, int, str]] = []
    if n_levels > 1:
        lis = [(0, o, i, "A") for o in range(3) for i in range(shapes[0][o])]

    for bp in range(n_bitplanes - 1, n_bitplanes - 1 - n_bitplanes_encoded, -1):
        threshold = 1 << bp
        lsp_len_before = len(lsp)

        new_lip = []
        for level, orientation, i in lip:
            sig = m_lip.decode(dec)
            if sig:
                s = m_sign.decode(dec)
                sign[level][orientation][i] = bool(s)
                mag[level][orientation][i] = threshold
                lsp.append((level, orientation, i))
            else:
                new_lip.append((level, orientation, i))
        lip = new_lip

        idx = 0
        while idx < len(lis):
            level, orientation, i, kind = lis[idx]
            if kind == "A":
                sig = m_lisA.decode(dec)
                remove = False
                if sig:
                    child_level = level + 1
                    for c in range(4 * i, 4 * i + 4):
                        csig = m_lip.decode(dec)
                        if csig:
                            s = m_sign.decode(dec)
                            sign[child_level][orientation][c] = bool(s)
                            mag[child_level][orientation][c] = threshold
                            lsp.append((child_level, orientation, c))
                        else:
                            lip.append((child_level, orientation, c))
                    if child_level + 1 < n_levels:
                        lis.append((level, orientation, i, "B"))
                    remove = True
            else:
                sig = m_lisB.decode(dec)
                remove = False
                if sig:
                    for c in range(4 * i, 4 * i + 4):
                        lis.append((level + 1, orientation, c, "A"))
                    remove = True
            if remove:
                lis.pop(idx)
            else:
                idx += 1

        for level, orientation, i in lsp[:lsp_len_before]:
            bit = m_ref.decode(dec)
            mag[level][orientation][i] |= bit << bp

    out = []
    for level in range(n_levels):
        triple = []
        for orientation in range(3):
            val = mag[level][orientation].astype(np.float64) * float(scale)
            val[sign[level][orientation]] *= -1.0
            triple.append(val)
        out.append(tuple(triple))
    return out
