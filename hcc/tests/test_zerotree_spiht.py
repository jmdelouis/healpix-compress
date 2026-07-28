import zlib

import numpy as np
import pytest

from hcc.compression.grad_dwt import _subband_slices, dwt_forward
from hcc.compression.zerotree import (
    decode_zerotree,
    decode_zerotree_spiht,
    encode_zerotree,
    encode_zerotree_spiht,
)


def _patchy_tile(n=64, seed=0):
    """A handful of homogeneous regions plus sharp edges -- a rough proxy
    for land-cover-like spatial structure, as opposed to fully dense
    texture/noise (which has no large below-threshold regions for
    set-partitioning to exploit -- see this module's data-dependence
    tests below)."""
    rng = np.random.default_rng(seed)
    img = np.zeros((n, n), dtype=np.float32)
    for _ in range(5):
        cx, cy = rng.integers(0, n, size=2)
        w, h = rng.integers(max(2, n // 8), max(3, n // 3), size=2)
        val = rng.uniform(0.2, 0.9)
        x0, x1 = max(0, cx - w // 2), min(n, cx + w // 2)
        y0, y1 = max(0, cy - h // 2), min(n, cy + h // 2)
        img[y0:y1, x0:x1] = val
    img[n // 2 - 1 : n // 2 + 1, :] = 0.95
    img += 0.01 * rng.standard_normal((n, n)).astype(np.float32)
    return img


def _direct_global_quantize(coeffs_arr, subbands, scale, quantize_bits):
    max_level = (1 << quantize_bits) - 1
    out = np.zeros_like(coeffs_arr)
    for sl in subbands[1:]:
        block = coeffs_arr[sl]
        mag = np.clip(np.round(np.abs(block) / scale), 0, max_level)
        val = mag.astype(np.float32) * scale
        val[block < 0] *= -1.0
        out[sl] = val
    return out


@pytest.mark.parametrize("quantize_bits", [4, 8])
@pytest.mark.parametrize("levels", [2, 3])
def test_full_bitplane_roundtrip_matches_global_quantization(quantize_bits, levels):
    tile = _patchy_tile(64, seed=quantize_bits + levels)
    coeffs_arr, slices = dwt_forward(tile, "db4", levels)
    subbands = _subband_slices(slices)

    payload, scales, n_bp, n_bp_enc = encode_zerotree_spiht(
        coeffs_arr, subbands, quantize_bits=quantize_bits, clip_percentile=0.5
    )
    assert n_bp_enc == n_bp == quantize_bits
    # global scale -> every entry identical
    assert np.all(scales == scales[0])

    decoded = decode_zerotree_spiht(payload, subbands, coeffs_arr.shape, scales, n_bp, n_bp_enc)
    expected = _direct_global_quantize(coeffs_arr, subbands, float(scales[0]), quantize_bits)
    np.testing.assert_allclose(decoded, expected, atol=1e-4)


def test_all_zero_tile():
    tile = np.zeros((32, 32), dtype=np.float32)
    coeffs_arr, slices = dwt_forward(tile, "db2", 2)
    subbands = _subband_slices(slices)
    payload, scales, n_bp, n_bp_enc = encode_zerotree_spiht(
        coeffs_arr, subbands, quantize_bits=8, clip_percentile=0.5
    )
    decoded = decode_zerotree_spiht(payload, subbands, coeffs_arr.shape, scales, n_bp, n_bp_enc)
    assert np.all(decoded == 0.0)


def test_levels_one_no_lis_ever_used():
    tile = _patchy_tile(64, seed=7)
    coeffs_arr, slices = dwt_forward(tile, "db2", 1)
    subbands = _subband_slices(slices)
    payload, scales, n_bp, n_bp_enc = encode_zerotree_spiht(
        coeffs_arr, subbands, quantize_bits=8, clip_percentile=0.5
    )
    decoded = decode_zerotree_spiht(payload, subbands, coeffs_arr.shape, scales, n_bp, n_bp_enc)
    expected = _direct_global_quantize(coeffs_arr, subbands, float(scales[0]), 8)
    np.testing.assert_allclose(decoded, expected, atol=1e-4)


def test_truncation_is_monotonically_better_with_more_bitplanes():
    tile = _patchy_tile(96, seed=3)
    coeffs_arr, slices = dwt_forward(tile, "db4", 3)
    subbands = _subband_slices(slices)

    full_payload, full_scales, n_bp, n_bp_full = encode_zerotree_spiht(
        coeffs_arr, subbands, quantize_bits=10, clip_percentile=0.5
    )
    reference = decode_zerotree_spiht(
        full_payload, subbands, coeffs_arr.shape, full_scales, n_bp, n_bp_full
    )

    prev_size = -1
    prev_rmse = np.inf
    for max_bp in (1, 3, 5, 7, 10):
        payload, scales, n_bp2, n_bp_enc = encode_zerotree_spiht(
            coeffs_arr, subbands, quantize_bits=10, clip_percentile=0.5, max_bitplanes=max_bp
        )
        assert n_bp_enc == max_bp
        decoded = decode_zerotree_spiht(payload, subbands, coeffs_arr.shape, scales, n_bp2, n_bp_enc)
        rmse = float(np.sqrt(np.mean((decoded - reference) ** 2)))

        assert len(payload) >= prev_size
        assert rmse <= prev_rmse + 1e-9
        prev_size = len(payload)
        prev_rmse = rmse

    assert prev_rmse == pytest.approx(0.0, abs=1e-6)


def test_invalid_quantize_bits_raises():
    tile = _patchy_tile(32, seed=0)
    coeffs_arr, slices = dwt_forward(tile, "db2", 2)
    subbands = _subband_slices(slices)
    with pytest.raises(ValueError):
        encode_zerotree_spiht(coeffs_arr, subbands, quantize_bits=0, clip_percentile=0.5)


def test_wrong_scales_size_raises():
    tile = _patchy_tile(32, seed=0)
    coeffs_arr, slices = dwt_forward(tile, "db2", 2)
    subbands = _subband_slices(slices)
    payload, scales, n_bp, n_bp_enc = encode_zerotree_spiht(
        coeffs_arr, subbands, quantize_bits=8, clip_percentile=0.5
    )
    with pytest.raises(ValueError):
        decode_zerotree_spiht(payload, subbands, coeffs_arr.shape, scales[:-1], n_bp, n_bp_enc)


def _zlib_baseline_bytes(coeffs_arr, subbands, bits=8, clip_pct=0.5):
    raw = bytearray()
    max_level = (1 << bits) - 1
    for sl in subbands[1:]:
        block = coeffs_arr[sl]
        abs_block = np.abs(block)
        max_abs = float(np.percentile(abs_block, 100 - clip_pct)) if abs_block.size > 1 else float(np.max(abs_block))
        if max_abs == 0.0:
            max_abs = 1.0
        scale = max_abs / max_level
        mag = np.clip(np.round(np.abs(block) / scale), 0, max_level).astype(np.uint8)
        raw += mag.tobytes()
        raw += np.packbits(block < 0).tobytes()
    return len(zlib.compress(bytes(raw), level=9))


def test_beats_zlib_and_plain_zerotree_on_patchy_data():
    """The key claim: on data with genuine spatial sparsity (unlike fully
    dense texture), real set-partitioning should beat both the production
    per-subband zlib baseline and the earlier context-only zerotree
    coder."""
    tile = _patchy_tile(128, seed=42)
    coeffs_arr, slices = dwt_forward(tile, "db4", 4)
    subbands = _subband_slices(slices)

    spiht_payload, _, _, _ = encode_zerotree_spiht(
        coeffs_arr, subbands, quantize_bits=8, clip_percentile=0.5
    )
    zt_payload, _, _, _ = encode_zerotree(coeffs_arr, subbands, quantize_bits=8, clip_percentile=0.5)
    zlib_bytes = _zlib_baseline_bytes(coeffs_arr, subbands, bits=8, clip_pct=0.5)

    assert len(spiht_payload) < zlib_bytes
    assert len(spiht_payload) < len(zt_payload)


def test_dense_textured_data_gain_is_small_or_neutral():
    """Data-dependence sanity check: without large homogeneous regions,
    set-partitioning has little to skip, so the gain over zlib should be
    small (this guards against a future regression making the module
    silently worse than zlib on this kind of data, as an earlier
    per-subband-scale prototype did)."""
    n = 96
    rng = np.random.default_rng(5)
    xx, yy = np.meshgrid(np.linspace(0, 6, n), np.linspace(0, 6, n))
    tile = (0.5 + 0.3 * np.sin(xx) * np.cos(yy) + 0.1 * np.sin(3 * xx + yy)).astype(np.float32)
    tile += (0.03 * rng.standard_normal((n, n))).astype(np.float32)

    coeffs_arr, slices = dwt_forward(tile, "db4", 3)
    subbands = _subband_slices(slices)

    spiht_payload, _, _, _ = encode_zerotree_spiht(
        coeffs_arr, subbands, quantize_bits=8, clip_percentile=0.5
    )
    zlib_bytes = _zlib_baseline_bytes(coeffs_arr, subbands, bits=8, clip_pct=0.5)

    assert len(spiht_payload) < zlib_bytes * 1.1  # not meaningfully worse
