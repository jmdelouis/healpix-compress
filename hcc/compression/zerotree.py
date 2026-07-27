"""
Context-adaptive, embedded (bit-plane) zero-tree coding for quantised DWT
detail coefficients -- an EBCOT-equivalent entropy back-end that adds the
cross-scale ("zero-tree") correlation missing from
:mod:`hcc.compression.entropy`'s per-subband order-0 model.

Why this module exists
-----------------------
A controlled measurement (see this project's paper Discussion) found that
``hcc.compression.entropy``'s adaptive range coder -- order-0 per subband,
with zero-run collapsing -- only gains ~5% over plain zlib on realistic
detail-coefficient data. The reason: an order-0 model has no way to exploit
that a coefficient's magnitude is *predictable from context* -- specifically
from (a) whether its immediate spatial neighbours in the same subband are
already significant, and (b) whether its "parent" coefficient -- the
co-located coefficient one scale coarser, same orientation -- is already
significant. This second, cross-scale kind of correlation is exactly what
JPEG 2000's EBCOT and the classic EZW/SPIHT coders exploit, and it is a
particularly natural fit for this codebase: HEALPix's NESTED indexing is
itself a quadtree, and a dyadic 2-D DWT's per-orientation subbands across
levels form an *identical* quadtree (each coefficient at level L has
exactly 4 co-located children at level L-1, same orientation) -- so the
parent-lookup this module needs is a single integer division, no extra
bookkeeping.

Algorithm (SPIHT/EZW-style zero-tree, with EBCOT-style context-adaptive
arithmetic coding in place of SPIHT's plain bit flags)
--------------------------------------------------------------------------
1. Per detail subband, derive an integer magnitude array and a sign array
   from the float coefficients, using the same percentile-clip scale
   derivation as :func:`hcc.compression.grad_dwt.compress_residuals`
   (``quantize_bits`` here means "number of bit-planes", not a fixed
   dtype -- see below).
2. Coefficients are coded **bit-plane by bit-plane, MSB first**, and
   *within* each bit-plane, subband by subband from the coarsest detail
   level to the finest. This traversal order guarantees that whenever a
   coefficient is coded, its parent's current significance state (in a
   coarser, already-fully-processed subband) is known.
3. For each not-yet-significant coefficient, the significance bit for the
   current plane is coded using an adaptively-modelled context combining:
   whether its parent is already significant (3 states: no parent /
   parent insignificant / parent significant), and how many of its causal
   (already-visited, left/up) same-subband neighbours are already
   significant (3 buckets: 0 / 1 / 2). 9 contexts per subband, each with
   its own adaptive binary model.
4. The first time a coefficient becomes significant, its sign is coded
   (one shared adaptive binary model).
5. Once significant, a coefficient's later bit-planes are coded as
   "refinement" bits (one adaptive binary model per subband).

This produces a genuinely **embedded** bitstream: stopping after fewer
bit-planes than were encoded (``max_bitplanes`` at encode time, or simply
truncating the returned payload and adjusting ``n_bitplanes_encoded`` at
decode time) yields a smaller, coarser -- but still validly decodable --
reconstruction, exactly like a truncated JPEG 2000 codestream.

Honest scope / what this is *not*
----------------------------------
This is a coefficient-level realisation of EBCOT/SPIHT's core mechanisms
(bit-plane embedding + cross-scale + within-scale context modelling), not
a literal reimplementation of either:

* No code-block partitioning or independent per-block decoding regions
  (real EBCOT).
* No post-compression rate-distortion-optimal truncation search
  (PCRD-opt, real EBCOT) -- truncation here is simply "stop after N
  bit-planes", not the Lagrangian-optimal byte-exact truncation point.
* No SPIHT-style set-partitioning (LIS/LIP/LSP lists) to skip whole
  all-zero quadtree subtrees in O(1); every coefficient's significance is
  tested individually at every bit-plane until it becomes significant.
  This is the main scope simplification to be aware of: it is *correct*
  and *does* capture the cross-scale statistical dependency (via context
  modelling), but it does not get SPIHT's extra speed/rate advantage from
  skipping whole insignificant subtrees in a single symbol. Revisit this
  if profiling shows it is the bottleneck for production use.

Public API
----------
encode_zerotree : Encode a full DWT coefficient array's detail subbands.
decode_zerotree : Exact (or, if truncated, best-effort progressive)
    inverse of :func:`encode_zerotree`.
"""

from __future__ import annotations

import struct

import numpy as np

from hcc.compression.entropy import _AdaptiveModel, _RangeDecoder, _RangeEncoder

__all__ = ["encode_zerotree", "decode_zerotree"]


# ─────────────────────────────────────────────────────────────────
#  Subband quadtree geometry
# ─────────────────────────────────────────────────────────────────


def _detail_subband_geometry(n_detail_subbands: int) -> list[tuple[int, int]]:
    """Return ``(level_index, parent_index)`` for each detail subband.

    Detail subbands are ordered coarsest-level-first, 3 orientations per
    level (the order produced by
    :func:`hcc.compression.grad_dwt._subband_slices` with the LL entry
    dropped). ``level_index`` is 0 for the coarsest level (closest to LL)
    up to ``levels - 1`` for the finest. ``parent_index`` is the index (in
    the same list) of the co-located, same-orientation subband one level
    coarser, or ``-1`` if this is a coarsest-level (root) subband -- the
    classic EZW/SPIHT convention of rooting zero-trees at the coarsest
    detail level, not at LL.
    """
    if n_detail_subbands % 3 != 0:
        raise ValueError(
            f"Expected a multiple of 3 detail subbands (3 orientations per "
            f"level), got {n_detail_subbands}."
        )
    geometry = []
    for i in range(n_detail_subbands):
        level_index = i // 3
        parent_index = i - 3 if level_index > 0 else -1
        geometry.append((level_index, parent_index))
    return geometry


# ─────────────────────────────────────────────────────────────────
#  Contexts
# ─────────────────────────────────────────────────────────────────

# parent bucket: 0 = no parent (root subband), 1 = parent insignificant,
# 2 = parent significant.
_N_PARENT_BUCKETS = 3
# neighbour bucket: 0, 1, or 2-plus already-significant causal (left/up)
# same-subband neighbours.
_N_NEIGHBOR_BUCKETS = 3
_N_CONTEXTS = _N_PARENT_BUCKETS * _N_NEIGHBOR_BUCKETS


def _context_index(parent_bucket: int, neighbor_bucket: int) -> int:
    return parent_bucket * _N_NEIGHBOR_BUCKETS + neighbor_bucket


# ─────────────────────────────────────────────────────────────────
#  Encode
# ─────────────────────────────────────────────────────────────────


def encode_zerotree(
    coeffs_arr: np.ndarray,
    subbands: list,
    quantize_bits: int,
    clip_percentile: float,
    max_bitplanes: int | None = None,
) -> tuple[bytes, np.ndarray, int, int]:
    """Encode a DWT coefficient array's detail subbands with context-adaptive
    zero-tree bit-plane coding.

    Parameters
    ----------
    coeffs_arr : numpy.ndarray of shape (n, n), dtype float
        Full DWT coefficient array, as returned by
        :func:`hcc.compression.grad_dwt.dwt_forward` (LL region included
        but ignored here -- LL is handled separately, losslessly, exactly
        as in :func:`hcc.compression.grad_dwt.compress_residuals`).
    subbands : list
        Subband slice list from
        :func:`hcc.compression.grad_dwt._subband_slices` (LL first, then
        detail subbands coarsest-level-first).
    quantize_bits : int
        Number of bit-planes to derive the per-subband quantisation scale
        from (``max_level = 2**quantize_bits - 1``), analogous to
        ``quantize_bits`` in
        :func:`hcc.compression.grad_dwt.compress_residuals` but not
        restricted to ``{8, 16}`` here -- any positive integer is valid
        since bit-plane coding has no fixed-width dtype to match.
    clip_percentile : float
        Same meaning as in :func:`hcc.compression.grad_dwt.compress_residuals`.
    max_bitplanes : int, optional
        If given, stop after coding this many bit-planes (MSB first) --
        produces a smaller, lossier, but still validly decodable payload
        (the embedded/progressive property). Defaults to ``quantize_bits``
        (i.e. full, lossless-at-the-chosen-quantisation coding).

    Returns
    -------
    payload : bytes
        The entropy-coded bitstream.
    scales : numpy.ndarray of shape (n_detail_subbands,), dtype float32
        Per-subband quantisation step size (same role as ``q_scale`` in
        :func:`hcc.compression.grad_dwt.compress_residuals`).
    n_bitplanes : int
        Total bit-planes the scale derivation was calibrated for
        (``quantize_bits``).
    n_bitplanes_encoded : int
        How many bit-planes were actually coded (``<= n_bitplanes``).
        Required, together with *scales* and *n_bitplanes*, to decode.

    See Also
    --------
    decode_zerotree : Inverse of this function.
    """
    if quantize_bits < 1:
        raise ValueError(f"quantize_bits must be >= 1, got {quantize_bits!r}")
    detail_subbands = subbands[1:]
    n_sub = len(detail_subbands)
    geometry = _detail_subband_geometry(n_sub)

    n_bitplanes = quantize_bits
    n_bitplanes_encoded = n_bitplanes if max_bitplanes is None else min(max_bitplanes, n_bitplanes)
    if n_bitplanes_encoded < 0:
        raise ValueError("max_bitplanes must be >= 0")

    max_level = (1 << quantize_bits) - 1
    mags: list[np.ndarray] = []
    signs: list[np.ndarray] = []
    scales = np.empty(n_sub, dtype=np.float32)

    for i, sl in enumerate(detail_subbands):
        block = coeffs_arr[sl]
        abs_block = np.abs(block)
        if clip_percentile > 0 and abs_block.size > 1:
            max_abs = float(np.percentile(abs_block, 100.0 - clip_percentile))
        else:
            max_abs = float(np.max(abs_block)) if abs_block.size else 0.0
        if max_abs == 0.0:
            max_abs = 1.0
        scale = np.float32(max_abs / max_level)
        mag = np.clip(np.round(np.abs(block) / scale), 0, max_level).astype(np.int64)
        mags.append(mag)
        signs.append(block < 0)
        scales[i] = scale

    enc = _RangeEncoder()
    sig_models = [[_AdaptiveModel(2) for _ in range(_N_CONTEXTS)] for _ in range(n_sub)]
    refine_models = [_AdaptiveModel(2) for _ in range(n_sub)]
    sign_model = _AdaptiveModel(2)
    significance = [np.zeros(mag.shape, dtype=bool) for mag in mags]

    for b in range(n_bitplanes - 1, n_bitplanes - 1 - n_bitplanes_encoded, -1):
        for i in range(n_sub):
            _, parent_idx = geometry[i]
            mag = mags[i]
            sig = significance[i]
            sign_arr = signs[i]
            parent_sig = significance[parent_idx] if parent_idx >= 0 else None
            rows, cols = mag.shape
            sig_models_i = sig_models[i]
            refine_model_i = refine_models[i]

            mag_rows = mag.tolist()
            sig_rows = sig.tolist()

            for r in range(rows):
                mag_row = mag_rows[r]
                sig_row = sig_rows[r]
                sig_row_up = sig_rows[r - 1] if r > 0 else None
                parent_row = parent_sig[r // 2].tolist() if parent_sig is not None else None
                for c in range(cols):
                    if sig_row[c]:
                        bit = (mag_row[c] >> b) & 1
                        refine_model_i.encode(enc, bit)
                        continue

                    if parent_row is None:
                        parent_bucket = 0
                    else:
                        parent_bucket = 2 if parent_row[c // 2] else 1

                    n_neighbors_sig = 0
                    if c > 0 and sig_row[c - 1]:
                        n_neighbors_sig += 1
                    if sig_row_up is not None and sig_row_up[c]:
                        n_neighbors_sig += 1
                    neighbor_bucket = min(n_neighbors_sig, 2)

                    ctx = _context_index(parent_bucket, neighbor_bucket)
                    bit = (mag_row[c] >> b) & 1
                    sig_models_i[ctx].encode(enc, bit)

                    if bit:
                        sig[r, c] = True
                        sig_row[c] = True
                        sign_model.encode(enc, 1 if sign_arr[r, c] else 0)

    payload = enc.finish()
    return payload, scales, n_bitplanes, n_bitplanes_encoded


# ─────────────────────────────────────────────────────────────────
#  Decode
# ─────────────────────────────────────────────────────────────────


def decode_zerotree(
    payload: bytes,
    subbands: list,
    shape: tuple,
    scales: np.ndarray,
    n_bitplanes: int,
    n_bitplanes_encoded: int,
) -> np.ndarray:
    """Exact (or, if truncated, best-effort progressive) inverse of
    :func:`encode_zerotree`.

    Parameters
    ----------
    payload : bytes
        Bitstream from :func:`encode_zerotree`.
    subbands : list
        Same subband slice list used at encode time.
    shape : tuple[int, int]
        Shape of the full coefficient array (``coeffs_arr.shape`` at
        encode time).
    scales : numpy.ndarray of shape (n_detail_subbands,), dtype float32
        Per-subband quantisation step sizes from :func:`encode_zerotree`.
    n_bitplanes : int
        Total bit-planes the scale derivation was calibrated for, from
        :func:`encode_zerotree`.
    n_bitplanes_encoded : int
        How many bit-planes were actually coded, from
        :func:`encode_zerotree`. May be less than *n_bitplanes* (a
        truncated/progressive payload); undecoded low bit-planes are
        reconstructed as 0 (i.e. the coefficient magnitude is the best
        estimate given only the coded high bit-planes).

    Returns
    -------
    numpy.ndarray of shape *shape*, dtype float32
        Reconstructed *detail* coefficients (LL region left as 0 -- the
        caller is responsible for placing the separately-stored, lossless
        LL subband, exactly as with
        :func:`hcc.compression.grad_dwt.decompress_residuals`).
    """
    detail_subbands = subbands[1:]
    n_sub = len(detail_subbands)
    geometry = _detail_subband_geometry(n_sub)

    scales = np.ravel(np.asarray(scales, dtype=np.float32))
    if scales.size != n_sub:
        raise ValueError(
            f"scales has {scales.size} element(s) but {n_sub} detail "
            "subbands are expected. scales must come from the matching "
            "encode_zerotree() call."
        )

    shapes = [tuple(s.stop - (s.start or 0) for s in sl) for sl in detail_subbands]
    mags = [np.zeros(shape, dtype=np.int64) for shape in shapes]
    signs = [np.zeros(shape, dtype=bool) for shape in shapes]
    significance = [np.zeros(shape, dtype=bool) for shape in shapes]

    dec = _RangeDecoder(payload)
    sig_models = [[_AdaptiveModel(2) for _ in range(_N_CONTEXTS)] for _ in range(n_sub)]
    refine_models = [_AdaptiveModel(2) for _ in range(n_sub)]
    sign_model = _AdaptiveModel(2)

    for b in range(n_bitplanes - 1, n_bitplanes - 1 - n_bitplanes_encoded, -1):
        for i in range(n_sub):
            _, parent_idx = geometry[i]
            mag = mags[i]
            sig = significance[i]
            sign_arr = signs[i]
            parent_sig = significance[parent_idx] if parent_idx >= 0 else None
            rows, cols = mag.shape
            sig_models_i = sig_models[i]
            refine_model_i = refine_models[i]

            for r in range(rows):
                parent_row = parent_sig[r // 2] if parent_sig is not None else None
                for c in range(cols):
                    if sig[r, c]:
                        bit = refine_model_i.decode(dec)
                        if bit:
                            mag[r, c] |= (1 << b)
                        continue

                    if parent_row is None:
                        parent_bucket = 0
                    else:
                        parent_bucket = 2 if parent_row[c // 2] else 1

                    n_neighbors_sig = 0
                    if c > 0 and sig[r, c - 1]:
                        n_neighbors_sig += 1
                    if r > 0 and sig[r - 1, c]:
                        n_neighbors_sig += 1
                    neighbor_bucket = min(n_neighbors_sig, 2)

                    ctx = _context_index(parent_bucket, neighbor_bucket)
                    bit = sig_models_i[ctx].decode(dec)
                    if bit:
                        mag[r, c] |= (1 << b)
                        sig[r, c] = True
                        sign_arr[r, c] = bool(sign_model.decode(dec))

    out = np.zeros(shape, dtype=np.float32)
    for sl, mag, sign_arr, scale in zip(detail_subbands, mags, signs, scales):
        val = mag.astype(np.float32) * float(scale)
        val[sign_arr] *= -1.0
        out[sl] = val
    return out
