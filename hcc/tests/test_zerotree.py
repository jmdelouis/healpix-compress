import numpy as np
import pytest

from hcc.compression.grad_dwt import _coeff_slices, _subband_slices, dwt_forward
from hcc.compression.zerotree import decode_zerotree, encode_zerotree


def _make_tile(n=64, seed=0, noise=0.05):
    rng = np.random.default_rng(seed)
    x = np.linspace(0, 4, n, dtype="float32")
    y = np.linspace(0, 4, n, dtype="float32")
    xx, yy = np.meshgrid(x, y)
    tile = (0.5 + 0.3 * np.sin(xx) * np.cos(yy)).astype(np.float32)
    tile += (noise * rng.standard_normal((n, n))).astype(np.float32)
    return tile.astype(np.float32)


def _direct_quantize(coeffs_arr, subbands, quantize_bits, clip_percentile):
    """Reference implementation of encode_zerotree's per-subband scale +
    round-to-nearest quantisation, used to check bit-exact full-bitplane
    round-trips independently of the bit-plane machinery itself."""
    detail_subbands = subbands[1:]
    max_level = (1 << quantize_bits) - 1
    out = np.zeros_like(coeffs_arr)
    for sl in detail_subbands:
        block = coeffs_arr[sl]
        abs_block = np.abs(block)
        if clip_percentile > 0 and abs_block.size > 1:
            max_abs = float(np.percentile(abs_block, 100.0 - clip_percentile))
        else:
            max_abs = float(np.max(abs_block)) if abs_block.size else 0.0
        if max_abs == 0.0:
            max_abs = 1.0
        scale = np.float32(max_abs / max_level)
        mag = np.clip(np.round(np.abs(block) / scale), 0, max_level)
        val = mag.astype(np.float32) * scale
        val[block < 0] *= -1.0
        out[sl] = val
    return out


@pytest.mark.parametrize("quantize_bits", [4, 8, 10])
@pytest.mark.parametrize("levels", [2, 3])
def test_full_bitplane_roundtrip_matches_direct_quantization(quantize_bits, levels):
    tile = _make_tile(64, seed=quantize_bits + levels)
    coeffs_arr, slices = dwt_forward(tile, "db4", levels)
    subbands = _subband_slices(slices)

    expected = _direct_quantize(coeffs_arr, subbands, quantize_bits, clip_percentile=0.5)

    payload, scales, n_bp, n_bp_enc = encode_zerotree(
        coeffs_arr, subbands, quantize_bits=quantize_bits, clip_percentile=0.5
    )
    assert n_bp_enc == n_bp == quantize_bits

    decoded = decode_zerotree(payload, subbands, coeffs_arr.shape, scales, n_bp, n_bp_enc)
    np.testing.assert_allclose(decoded, expected, atol=1e-6)


def test_all_zero_tile():
    tile = np.zeros((32, 32), dtype=np.float32)
    coeffs_arr, slices = dwt_forward(tile, "db2", 2)
    subbands = _subband_slices(slices)

    payload, scales, n_bp, n_bp_enc = encode_zerotree(
        coeffs_arr, subbands, quantize_bits=8, clip_percentile=0.5
    )
    decoded = decode_zerotree(payload, subbands, coeffs_arr.shape, scales, n_bp, n_bp_enc)
    assert np.all(decoded == 0.0)


def test_levels_one_no_parent_anywhere():
    """levels=1 means every detail subband is a zero-tree root (no
    parent) -- exercises the parent_bucket == 0 path exclusively."""
    tile = _make_tile(64, seed=7)
    coeffs_arr, slices = dwt_forward(tile, "db2", 1)
    subbands = _subband_slices(slices)

    payload, scales, n_bp, n_bp_enc = encode_zerotree(
        coeffs_arr, subbands, quantize_bits=8, clip_percentile=0.5
    )
    decoded = decode_zerotree(payload, subbands, coeffs_arr.shape, scales, n_bp, n_bp_enc)
    expected = _direct_quantize(coeffs_arr, subbands, 8, clip_percentile=0.5)
    np.testing.assert_allclose(decoded, expected, atol=1e-6)


def test_truncation_is_monotonically_better_with_more_bitplanes():
    """The embedded/progressive property: decoding with more bitplanes
    should never increase RMSE against the fully-decoded reference, and
    the payload should grow monotonically with bitplane count."""
    tile = _make_tile(96, seed=3, noise=0.02)
    coeffs_arr, slices = dwt_forward(tile, "db4", 3)
    subbands = _subband_slices(slices)

    full_payload, full_scales, n_bp, n_bp_full = encode_zerotree(
        coeffs_arr, subbands, quantize_bits=10, clip_percentile=0.5
    )
    reference = decode_zerotree(
        full_payload, subbands, coeffs_arr.shape, full_scales, n_bp, n_bp_full
    )

    prev_size = -1
    prev_rmse = np.inf
    for max_bp in (1, 3, 5, 7, 10):
        payload, scales, n_bp2, n_bp_enc = encode_zerotree(
            coeffs_arr, subbands, quantize_bits=10, clip_percentile=0.5, max_bitplanes=max_bp
        )
        assert n_bp_enc == max_bp
        decoded = decode_zerotree(payload, subbands, coeffs_arr.shape, scales, n_bp2, n_bp_enc)
        rmse = float(np.sqrt(np.mean((decoded - reference) ** 2)))

        assert len(payload) >= prev_size
        assert rmse <= prev_rmse + 1e-9
        prev_size = len(payload)
        prev_rmse = rmse

    assert prev_rmse == pytest.approx(0.0, abs=1e-6)  # max_bp == n_bp -> exact


def test_max_bitplanes_zero_gives_all_zero_output():
    tile = _make_tile(64, seed=11)
    coeffs_arr, slices = dwt_forward(tile, "db2", 2)
    subbands = _subband_slices(slices)

    payload, scales, n_bp, n_bp_enc = encode_zerotree(
        coeffs_arr, subbands, quantize_bits=8, clip_percentile=0.5, max_bitplanes=0
    )
    assert n_bp_enc == 0
    decoded = decode_zerotree(payload, subbands, coeffs_arr.shape, scales, n_bp, n_bp_enc)
    assert np.all(decoded == 0.0)


def test_invalid_quantize_bits_raises():
    tile = _make_tile(32, seed=0)
    coeffs_arr, slices = dwt_forward(tile, "db2", 2)
    subbands = _subband_slices(slices)
    with pytest.raises(ValueError):
        encode_zerotree(coeffs_arr, subbands, quantize_bits=0, clip_percentile=0.5)


def test_wrong_scales_size_raises():
    tile = _make_tile(32, seed=0)
    coeffs_arr, slices = dwt_forward(tile, "db2", 2)
    subbands = _subband_slices(slices)
    payload, scales, n_bp, n_bp_enc = encode_zerotree(
        coeffs_arr, subbands, quantize_bits=8, clip_percentile=0.5
    )
    with pytest.raises(ValueError):
        decode_zerotree(payload, subbands, coeffs_arr.shape, scales[:-1], n_bp, n_bp_enc)
