import numpy as np
import pytest

from hcc.compression.entropy import decode_coeffs, encode_coeffs
from hcc.compression.grad_dwt import _coeff_slices, _subband_slices


def _whole_array_subbands(shape):
    return [tuple(slice(0, s) for s in shape)]


def _roundtrip(arr, subbands):
    data = encode_coeffs(arr, subbands)
    out = decode_coeffs(data, subbands, arr.shape, arr.dtype)
    assert out.dtype == arr.dtype
    assert out.shape == arr.shape
    np.testing.assert_array_equal(out, arr)
    return data


def test_all_zeros_int8():
    arr = np.zeros((16, 16), dtype=np.int8)
    _roundtrip(arr, _whole_array_subbands(arr.shape))


def test_all_zeros_int16():
    arr = np.zeros((16, 16), dtype=np.int16)
    _roundtrip(arr, _whole_array_subbands(arr.shape))


def test_dense_random_int8_roundtrip():
    rng = np.random.default_rng(0)
    arr = rng.integers(-127, 128, size=(64, 64)).astype(np.int8)
    _roundtrip(arr, _whole_array_subbands(arr.shape))


def test_dense_random_int16_roundtrip():
    rng = np.random.default_rng(0)
    arr = rng.integers(-32767, 32768, size=(32, 32)).astype(np.int16)
    _roundtrip(arr, _whole_array_subbands(arr.shape))


def test_sparse_laplace_int8_roundtrip_and_compresses():
    rng = np.random.default_rng(1)
    lap = rng.laplace(0, 3, size=(128, 128))
    arr = np.clip(np.round(lap), -127, 127).astype(np.int8)
    subbands = _whole_array_subbands(arr.shape)
    data = _roundtrip(arr, subbands)
    # A sparse, peaked distribution should compress well below the raw
    # int8 size (this is the whole point of run-length + adaptive coding).
    assert len(data) < arr.nbytes


def test_extreme_values_int8():
    arr = np.full((8, 8), 127, dtype=np.int8)
    arr[0, 0] = -127
    _roundtrip(arr, _whole_array_subbands(arr.shape))


def test_extreme_values_int16():
    arr = np.full((8, 8), 32767, dtype=np.int16)
    arr[0, 0] = -32767
    _roundtrip(arr, _whole_array_subbands(arr.shape))


def test_single_element():
    for val, dtype in ((0, np.int8), (5, np.int8), (-5, np.int16)):
        arr = np.array([[val]], dtype=dtype)
        _roundtrip(arr, _whole_array_subbands(arr.shape))


def test_multi_subband_matches_real_dwt_layout():
    """Round-trip against the actual quad-tree subband layout produced for
    a real (wavelet, levels, shape) combination, as used by
    hcc.compression.grad_dwt.compress_residuals -- not just a single
    whole-array subband."""
    rng = np.random.default_rng(2)
    shape = (256, 256)
    slices = _coeff_slices(shape, "db4", 3)
    subbands = _subband_slices(slices)
    assert len(subbands) > 1

    arr = np.zeros(shape, dtype=np.int8)
    for sl in subbands:
        block_shape = tuple(s.stop - (s.start or 0) for s in sl)
        arr[sl] = np.clip(
            np.round(rng.laplace(0, 2, size=block_shape)), -127, 127
        ).astype(np.int8)

    _roundtrip(arr, subbands)


def test_wrong_subband_count_raises():
    arr = np.zeros((8, 8), dtype=np.int8)
    data = encode_coeffs(arr, _whole_array_subbands(arr.shape))
    bad_subbands = [
        (slice(0, 4), slice(0, 8)),
        (slice(4, 8), slice(0, 8)),
    ]
    with pytest.raises(ValueError):
        decode_coeffs(data, bad_subbands, arr.shape, arr.dtype)
