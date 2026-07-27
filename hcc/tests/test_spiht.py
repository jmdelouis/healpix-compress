import numpy as np
import pytest

from hcc.compression.quadtree_wavelet import qt_dwt_forward, qt_dwt53_forward
from hcc.compression.spiht import decode_spiht, encode_spiht
from hcc.compression.zorder import flatten_tile


def _patchy_image(n=64, seed=0):
    """A handful of homogeneous regions plus a couple of sharp edges -- a
    rough proxy for land-cover-like spatial structure in real imagery
    (as opposed to fully dense noise, which has no large below-threshold
    regions for set-partitioning to exploit)."""
    rng = np.random.default_rng(seed)
    img = np.zeros((n, n))
    for _ in range(4):
        cx, cy = rng.integers(0, n, size=2)
        w, h = rng.integers(max(2, n // 8), max(3, n // 3), size=2)
        val = rng.uniform(0.2, 0.9)
        x0, x1 = max(0, cx - w // 2), min(n, cx + w // 2)
        y0, y1 = max(0, cy - h // 2), min(n, cy + h // 2)
        img[y0:y1, x0:x1] = val
    img[n // 2 - 1 : n // 2 + 1, :] = 0.95
    img += 0.01 * rng.standard_normal((n, n))
    return img


def _global_quantize(d, scale, bits):
    max_level = (1 << bits) - 1
    mag = np.clip(np.round(np.abs(d) / scale), 0, max_level)
    val = mag * scale
    val[d < 0] *= -1.0
    return val


@pytest.mark.parametrize("forward_fn", [qt_dwt_forward, qt_dwt53_forward])
@pytest.mark.parametrize("levels", [1, 2, 3])
def test_full_bitplane_roundtrip_matches_global_quantization(forward_fn, levels):
    img = _patchy_image(64, seed=levels)
    tile = flatten_tile(img)
    ll, details = forward_fn(tile, levels)

    payload, scale, n_bp, n_bp_enc = encode_spiht(details, quantize_bits=8, clip_percentile=0.5)
    assert n_bp_enc == n_bp == 8

    shapes = [tuple(sb.size for sb in triple) for triple in details]
    decoded = decode_spiht(payload, shapes, scale, n_bp, n_bp_enc)

    for (lh, hl, hh), (dlh, dhl, dhh) in zip(details, decoded):
        for original, out in zip((lh, hl, hh), (dlh, dhl, dhh)):
            expected = _global_quantize(original, scale, 8)
            np.testing.assert_allclose(out, expected, atol=1e-6)


def test_all_zero_details():
    tile = np.zeros(4**3)
    ll, details = qt_dwt_forward(tile, 3)
    payload, scale, n_bp, n_bp_enc = encode_spiht(details, quantize_bits=8, clip_percentile=0.5)
    shapes = [tuple(sb.size for sb in triple) for triple in details]
    decoded = decode_spiht(payload, shapes, scale, n_bp, n_bp_enc)
    for lh, hl, hh in decoded:
        assert np.all(lh == 0.0) and np.all(hl == 0.0) and np.all(hh == 0.0)


def test_single_level_no_lis_ever_used():
    """levels=1 -> the single detail level is also the finest -> no
    descendants exist anywhere -> LIS must start (and stay) empty, only
    LIP/LSP are exercised."""
    img = _patchy_image(16, seed=3)
    tile = flatten_tile(img)
    ll, details = qt_dwt_forward(tile, 1)
    payload, scale, n_bp, n_bp_enc = encode_spiht(details, quantize_bits=6, clip_percentile=0.5)
    shapes = [tuple(sb.size for sb in triple) for triple in details]
    decoded = decode_spiht(payload, shapes, scale, n_bp, n_bp_enc)
    for (lh, hl, hh), (dlh, dhl, dhh) in zip(details, decoded):
        for original, out in zip((lh, hl, hh), (dlh, dhl, dhh)):
            np.testing.assert_allclose(out, _global_quantize(original, scale, 6), atol=1e-6)


def test_truncation_is_monotonically_better_with_more_bitplanes():
    img = _patchy_image(64, seed=11)
    tile = flatten_tile(img)
    ll, details = qt_dwt_forward(tile, 3)

    full_payload, full_scale, n_bp, n_bp_full = encode_spiht(details, quantize_bits=10, clip_percentile=0.5)
    shapes = [tuple(sb.size for sb in triple) for triple in details]
    reference = decode_spiht(full_payload, shapes, full_scale, n_bp, n_bp_full)

    prev_size = -1
    prev_rmse = np.inf
    for max_bp in (1, 3, 5, 7, 10):
        payload, scale, n_bp2, n_bp_enc = encode_spiht(
            details, quantize_bits=10, clip_percentile=0.5, max_bitplanes=max_bp
        )
        assert n_bp_enc == max_bp
        decoded = decode_spiht(payload, shapes, scale, n_bp2, n_bp_enc)

        sq_err = 0.0
        n_total = 0
        for (rlh, rhl, rhh), (dlh, dhl, dhh) in zip(reference, decoded):
            for r, d in zip((rlh, rhl, rhh), (dlh, dhl, dhh)):
                sq_err += float(np.sum((r - d) ** 2))
                n_total += r.size
        rmse = float(np.sqrt(sq_err / n_total))

        assert len(payload) >= prev_size
        assert rmse <= prev_rmse + 1e-9
        prev_size = len(payload)
        prev_rmse = rmse

    assert prev_rmse == pytest.approx(0.0, abs=1e-6)


def test_invalid_quantize_bits_raises():
    tile = np.zeros(4**2)
    ll, details = qt_dwt_forward(tile, 2)
    with pytest.raises(ValueError):
        encode_spiht(details, quantize_bits=0, clip_percentile=0.5)


def test_beats_persubband_zlib_baseline_on_patchy_data():
    """The key claim motivating this module: on data with genuine spatial
    sparsity/homogeneity (unlike fully dense noise), SPIHT's set-
    partitioning should meaningfully beat the production
    per-subband-quantised + zlib baseline used elsewhere in this
    codebase."""
    import zlib

    img = _patchy_image(128, seed=42)
    tile = flatten_tile(img)
    ll, details = qt_dwt_forward(tile, 3)

    payload, scale, n_bp, n_bp_enc = encode_spiht(details, quantize_bits=8, clip_percentile=0.5)

    # Production-style baseline: per-subband scale + uniform quantisation + zlib.
    raw = bytearray()
    max_level = 255
    for triple in details:
        for d in triple:
            abs_d = np.abs(d)
            max_abs = float(np.percentile(abs_d, 99.5)) if abs_d.size > 1 else float(np.max(abs_d))
            if max_abs == 0.0:
                max_abs = 1.0
            sc = max_abs / max_level
            mag = np.clip(np.round(np.abs(d) / sc), 0, max_level).astype(np.uint8)
            raw += mag.tobytes()
            raw += np.packbits(d < 0).tobytes()
    zlib_bytes = len(zlib.compress(bytes(raw), level=9))

    assert len(payload) < zlib_bytes
