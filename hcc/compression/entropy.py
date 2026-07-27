"""
Adaptive entropy coding for quantised DWT/DCT coefficients.

Motivation
----------
:func:`hcc.compression.grad_dwt.compress_residuals` produces, per subband, a
grid of uniformly-quantised int8/int16 coefficients. Historically these were
serialised with ``zlib.compress(coeffs_q.tobytes())`` (DEFLATE: LZ77 +
static/canonical Huffman). A controlled measurement (see the project's paper
Discussion and this module's tests) showed that on real Sentinel-2 tiles this
already sits close to the *order-0* Shannon entropy of the quantised
coefficient stream for the densely-textured subbands (zlib is not "leaving
bits on the table" there) -- so the ~2-4x compression ceiling observed for
``healpix-compress`` versus 5-30x for JPEG2000/ZFP/SZ3 at comparable RMSE is
not a zlib-inefficiency bug. The real gap is architectural: JPEG2000's EBCOT
does context-adaptive, bit-plane arithmetic coding that exploits the strong
skew of wavelet coefficient magnitudes (most detail coefficients are exactly
or nearly zero) essentially for free, whereas quantise-then-zlib pays a
roughly uniform per-symbol cost regardless of magnitude.

This module closes part of that gap with a from-scratch, dependency-free
entropy back-end purpose-built for quantised wavelet residuals:

1. **Zero run-length coding.** Quantised detail coefficients are typically
   dominated by exact zeros (values that clip to zero after per-subband
   uniform quantisation). Runs of zeros are coded as a single (class,
   extra-bits) event instead of one symbol per zero, so sparse subbands
   collapse to very few coding events.
2. **Magnitude-class + sign + extra-bits coding for non-zero values**
   (the same scheme JPEG uses for AC/DC coefficients, generalised with
   arithmetic instead of static Huffman coding): each non-zero value is
   coded as its bit-length "class" (an adaptively-modelled small-alphabet
   symbol) plus raw extra bits pinpointing its exact value within that
   class, plus a sign bit.
3. **Adaptive, per-subband, order-0 arithmetic (range) coding** of the
   event type / run-length class / magnitude class symbols. Coding is
   restarted fresh per subband (LL vs. detail subbands, and different
   decomposition levels, have very different statistics -- see
   :func:`hcc.compression.grad_dwt.compress_residuals`'s per-subband
   quantisation rationale for the same argument applied to scale
   calibration).

The range coder itself is the classic "carryless" design (Subbotin-style):
32-bit low/range state, renormalisation by one byte at a time, no explicit
carry propagation (traded for a negligible -- a small fraction of a bit per
symbol -- compression loss in exchange for a simple, easily-audited
implementation). Frequency tables are tiny (at most ~24 symbols), so the
adaptive model update is a fast linear scan rather than a Fenwick tree.

Public API
----------
encode_coeffs : Entropy-code a full quantised coefficient array, subband by
    subband, into a single ``bytes`` payload.
decode_coeffs : Exact inverse of :func:`encode_coeffs`.

Both are drop-in replacements for
``zlib.compress(coeffs_q.tobytes())`` / ``zlib.decompress(...)`` wherever
*coeffs_q* is the standard-layout quantised DWT coefficient array together
with the ``subbands`` slice list from
:func:`hcc.compression.grad_dwt._subband_slices`.
"""

import struct

import numpy as np

__all__ = ["encode_coeffs", "decode_coeffs"]


# ─────────────────────────────────────────────────────────────────
#  Range coder (Subbotin-style carryless range coder)
# ─────────────────────────────────────────────────────────────────

_TOP = 1 << 24
_BOT = 1 << 16
_MASK = 0xFFFFFFFF

# Raw/uniform bits and adaptive-model total frequencies must each stay well
# below _BOT (65536) so that `range // tot` never underflows to 0 during
# renormalisation. 12 bits (4096) leaves a large safety margin.
_RAW_CHUNK_BITS = 12


class _RangeEncoder:
    __slots__ = ("low", "range", "out")

    def __init__(self):
        self.low = 0
        self.range = _MASK
        self.out = bytearray()

    def encode_freq(self, cum: int, freq: int, tot: int) -> None:
        self.range //= tot
        self.low = (self.low + cum * self.range) & _MASK
        self.range = (self.range * freq) & _MASK
        while True:
            cond_a = ((self.low ^ ((self.low + self.range) & _MASK)) & _MASK) < _TOP
            cond_b = False
            if not cond_a and self.range < _BOT:
                self.range = (-self.low) & (_BOT - 1)
                cond_b = True
            if not (cond_a or cond_b):
                break
            self.out.append((self.low >> 24) & 0xFF)
            self.low = (self.low << 8) & _MASK
            self.range = (self.range << 8) & _MASK

    def finish(self) -> bytes:
        for _ in range(4):
            self.out.append((self.low >> 24) & 0xFF)
            self.low = (self.low << 8) & _MASK
        return bytes(self.out)


class _RangeDecoder:
    __slots__ = ("data", "pos", "low", "range", "code")

    def __init__(self, data: bytes):
        self.data = data
        self.pos = 0
        self.low = 0
        self.range = _MASK
        self.code = 0
        for _ in range(4):
            self.code = ((self.code << 8) | self._next_byte()) & _MASK

    def _next_byte(self) -> int:
        if self.pos < len(self.data):
            b = self.data[self.pos]
            self.pos += 1
            return b
        return 0  # past-the-end padding; matches encoder's 4-byte flush

    def get_freq(self, tot: int) -> int:
        self.range //= tot
        return ((self.code - self.low) & _MASK) // self.range

    def decode_update(self, cum: int, freq: int, tot: int) -> None:
        self.low = (self.low + cum * self.range) & _MASK
        self.range = (self.range * freq) & _MASK
        while True:
            cond_a = ((self.low ^ ((self.low + self.range) & _MASK)) & _MASK) < _TOP
            cond_b = False
            if not cond_a and self.range < _BOT:
                self.range = (-self.low) & (_BOT - 1)
                cond_b = True
            if not (cond_a or cond_b):
                break
            self.code = ((self.code << 8) | self._next_byte()) & _MASK
            self.low = (self.low << 8) & _MASK
            self.range = (self.range << 8) & _MASK


def _encode_raw_bits(enc: _RangeEncoder, value: int, nbits: int) -> None:
    """Encode *value* (0 <= value < 2**nbits) as uniformly-distributed bits,
    split into <=_RAW_CHUNK_BITS chunks to keep every range-coder `tot`
    comfortably below _BOT."""
    remaining = nbits
    while remaining > 0:
        take = min(remaining, _RAW_CHUNK_BITS)
        remaining -= take
        chunk = (value >> remaining) & ((1 << take) - 1)
        enc.encode_freq(chunk, 1, 1 << take)


def _decode_raw_bits(dec: _RangeDecoder, nbits: int) -> int:
    value = 0
    remaining = nbits
    while remaining > 0:
        take = min(remaining, _RAW_CHUNK_BITS)
        remaining -= take
        tot = 1 << take
        chunk = dec.get_freq(tot)
        dec.decode_update(chunk, 1, tot)
        value = (value << take) | chunk
    return value


# ─────────────────────────────────────────────────────────────────
#  Adaptive order-0 model (small alphabet, linear-scan frequency table)
# ─────────────────────────────────────────────────────────────────


class _AdaptiveModel:
    """Adaptive order-0 frequency model over a small alphabet {0, ..., k-1}.

    Starts from a uniform prior (every symbol has frequency 1, so nothing is
    ever "impossible") and reinforces symbols as they are observed, with
    periodic halving to bound the total frequency (keeps `tot` well under
    `_BOT` and lets the model track locally-varying statistics instead of
    converging to a fixed long-run average).
    """

    __slots__ = ("freq", "total", "k", "increment", "max_total")

    def __init__(self, k: int, increment: int = 24, max_total: int = 1 << 13):
        self.k = k
        self.freq = [1] * k
        self.total = k
        self.increment = increment
        self.max_total = max_total

    def encode(self, enc: _RangeEncoder, sym: int) -> None:
        cum = 0
        freq = self.freq
        for i in range(sym):
            cum += freq[i]
        enc.encode_freq(cum, freq[sym], self.total)
        self._update(sym)

    def decode(self, dec: _RangeDecoder) -> int:
        target = dec.get_freq(self.total)
        cum = 0
        sym = 0
        freq = self.freq
        while cum + freq[sym] <= target:
            cum += freq[sym]
            sym += 1
        dec.decode_update(cum, freq[sym], self.total)
        self._update(sym)
        return sym

    def _update(self, sym: int) -> None:
        self.freq[sym] += self.increment
        self.total += self.increment
        if self.total > self.max_total:
            new_total = 0
            freq = self.freq
            for i in range(self.k):
                freq[i] = (freq[i] + 1) // 2
                new_total += freq[i]
            self.total = new_total


# Alphabet sizes for the bit-length "class" models. Class k (1-indexed, via
# int.bit_length()) covers values in [2**(k-1), 2**k - 1]; class 0 is unused
# (values are always >= 1 by construction -- 0 is handled by the separate
# zero-run path). Sized generously above the largest values we expect:
# run lengths up to a few million coefficients, magnitudes up to 32767
# (int16 quantisation ceiling).
_RUN_CLASSES = 24
_MAG_CLASSES = 17


def _encode_varlen(enc: _RangeEncoder, class_model: _AdaptiveModel, value: int) -> None:
    """Encode an integer value >= 1 as (bit-length class, extra bits)."""
    k = value.bit_length()
    class_model.encode(enc, k)
    extra_bits = k - 1
    if extra_bits:
        _encode_raw_bits(enc, value - (1 << extra_bits), extra_bits)


def _decode_varlen(dec: _RangeDecoder, class_model: _AdaptiveModel) -> int:
    k = class_model.decode(dec)
    extra_bits = k - 1
    extra = _decode_raw_bits(dec, extra_bits) if extra_bits else 0
    return (1 << extra_bits) + extra


# ─────────────────────────────────────────────────────────────────
#  Per-subband stream codec
# ─────────────────────────────────────────────────────────────────


def _encode_subband(values) -> bytes:
    """Entropy-code a 1-D sequence of signed integers (a flattened subband).

    Scheme: scan left to right; runs of exact zeros are coded as a single
    "zero-run" event (bit-length-classed run length); non-zero values are
    coded as sign + bit-length-classed magnitude. Fresh adaptive models are
    used per call (i.e. per subband), matching the per-subband quantisation
    calibration already used upstream in
    :func:`hcc.compression.grad_dwt.compress_residuals`.
    """
    n = len(values)
    enc = _RangeEncoder()
    if n == 0:
        return enc.finish()

    event_model = _AdaptiveModel(2)       # 0 = zero-run, 1 = non-zero
    run_class_model = _AdaptiveModel(_RUN_CLASSES)
    mag_class_model = _AdaptiveModel(_MAG_CLASSES)
    sign_model = _AdaptiveModel(2)        # 0 = positive, 1 = negative

    i = 0
    while i < n:
        v = values[i]
        if v == 0:
            run = 1
            while i + run < n and values[i + run] == 0:
                run += 1
            event_model.encode(enc, 0)
            _encode_varlen(enc, run_class_model, run)
            i += run
        else:
            event_model.encode(enc, 1)
            sign_model.encode(enc, 1 if v < 0 else 0)
            _encode_varlen(enc, mag_class_model, abs(v))
            i += 1

    return enc.finish()


def _decode_subband(data: bytes, n: int) -> list:
    out = [0] * n
    if n == 0:
        return out
    dec = _RangeDecoder(data)

    event_model = _AdaptiveModel(2)
    run_class_model = _AdaptiveModel(_RUN_CLASSES)
    mag_class_model = _AdaptiveModel(_MAG_CLASSES)
    sign_model = _AdaptiveModel(2)

    i = 0
    while i < n:
        event = event_model.decode(dec)
        if event == 0:
            run = _decode_varlen(dec, run_class_model)
            i += run  # `out` slots already 0-initialised
        else:
            sign = sign_model.decode(dec)
            mag = _decode_varlen(dec, mag_class_model)
            out[i] = -mag if sign else mag
            i += 1

    return out


# ─────────────────────────────────────────────────────────────────
#  Public array-level API
# ─────────────────────────────────────────────────────────────────


def _slice_shape(sl: tuple) -> tuple:
    # pywt's coeffs_to_array slices use `start=None` (meaning 0) for the
    # first/LL block, e.g. slice(None, 64, None) rather than slice(0, 64,
    # None) -- normalise before subtracting.
    return tuple(s.stop - (s.start or 0) for s in sl)


def encode_coeffs(coeffs_q: np.ndarray, subbands: list) -> bytes:
    """Entropy-code a quantised DWT coefficient array, subband by subband.

    Parameters
    ----------
    coeffs_q : numpy.ndarray of shape (n, n), dtype int8 or int16
        Quantised DWT coefficients in the standard quad-tree layout, as
        returned by :func:`hcc.compression.grad_dwt.compress_residuals`.
    subbands : list of tuple[slice, slice]
        Per-subband 2-D slices, as returned by
        :func:`hcc.compression.grad_dwt._subband_slices`. Must be the exact
        same list (same order) used to produce *coeffs_q*.

    Returns
    -------
    bytes
        Self-contained payload: a small header (subband count + one
        uint32 length per subband) followed by the concatenated per-subband
        entropy-coded streams. Pass to :func:`decode_coeffs` together with
        the same *subbands* list, the array *shape*, and *dtype* to invert.
    """
    parts = []
    for sl in subbands:
        block = coeffs_q[sl].ravel()
        parts.append(_encode_subband(block.tolist()))

    header = struct.pack("<I", len(parts))
    header += struct.pack(f"<{len(parts)}I", *(len(p) for p in parts))
    return header + b"".join(parts)


def decode_coeffs(
    data: bytes, subbands: list, shape: tuple, dtype: np.dtype
) -> np.ndarray:
    """Exact inverse of :func:`encode_coeffs`.

    Parameters
    ----------
    data : bytes
        Payload produced by :func:`encode_coeffs`.
    subbands : list of tuple[slice, slice]
        Same subband slice list used at encode time.
    shape : tuple[int, int]
        Shape of the full coefficient array (``coeffs_q.shape`` at encode
        time).
    dtype : numpy.dtype
        Output dtype (``np.int8`` or ``np.int16``, matching
        ``quantize_bits`` at encode time).

    Returns
    -------
    numpy.ndarray of shape *shape*, dtype *dtype*
        Reconstructed quantised coefficient array, bit-for-bit identical to
        the *coeffs_q* passed to :func:`encode_coeffs`.
    """
    (n_sub,) = struct.unpack_from("<I", data, 0)
    if n_sub != len(subbands):
        raise ValueError(
            f"Payload has {n_sub} subband(s) but {len(subbands)} were "
            "expected -- `subbands` must match the list used at encode time."
        )
    offset = 4
    lengths = struct.unpack_from(f"<{n_sub}I", data, offset)
    offset += 4 * n_sub

    out = np.zeros(shape, dtype=dtype)
    for sl, length in zip(subbands, lengths):
        chunk = data[offset : offset + length]
        offset += length
        block_shape = _slice_shape(sl)
        n_elems = int(np.prod(block_shape))
        values = _decode_subband(chunk, n_elems)
        out[sl] = np.asarray(values, dtype=dtype).reshape(block_shape)

    return out
