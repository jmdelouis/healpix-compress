import numpy as np
import pytest

from hcc.compression.quadtree_wavelet import qt_dwt_forward, qt_dwt_inverse
from hcc.compression.zorder import expand_tile, flatten_tile


@pytest.mark.parametrize("delta", [2, 3, 4, 5])
def test_perfect_reconstruction_all_level_counts(delta):
    n = 4**delta
    rng = np.random.default_rng(delta)
    x = rng.standard_normal(n)
    for levels in range(1, delta + 1):
        ll, details = qt_dwt_forward(x, levels)
        rec = qt_dwt_inverse(ll, details)
        np.testing.assert_allclose(rec, x, atol=1e-10)


def test_energy_is_preserved_orthonormal():
    rng = np.random.default_rng(0)
    n = 4**4
    x = rng.standard_normal(n)
    ll, details = qt_dwt_forward(x, 4)
    energy_in = np.sum(x**2)
    energy_out = np.sum(ll**2) + sum(
        np.sum(lh**2) + np.sum(hl**2) + np.sum(hh**2) for lh, hl, hh in details
    )
    assert energy_out == pytest.approx(energy_in, rel=1e-10)


def test_matches_pywt_haar_on_expanded_raster():
    """The native quadtree transform must be mathematically identical to
    the existing expand_tile -> pywt 2-D Haar path (just computed without
    ever forming the 2-D raster) -- this is the key correctness property
    tying it back to the rest of the codebase."""
    pywt = pytest.importorskip("pywt")
    rng = np.random.default_rng(1)
    n = 16
    x = rng.standard_normal(n)

    ll, details = qt_dwt_forward(x, 1)
    lh, hl, hh = details[0]

    x2d = expand_tile(x)
    cA, (cH, cV, cD) = pywt.dwt2(x2d, "haar", mode="periodization")

    np.testing.assert_allclose(ll, flatten_tile(cA), atol=1e-10)
    # orientation labels differ between this module's (lh, hl, hh) and
    # pywt's (cH, cV, cD) -- both are valid, self-consistent conventions;
    # what matters is that the *set* of subbands matches exactly.
    np.testing.assert_allclose(lh, flatten_tile(cV), atol=1e-10)
    np.testing.assert_allclose(hl, flatten_tile(cH), atol=1e-10)
    np.testing.assert_allclose(hh, flatten_tile(cD), atol=1e-10)


def test_levels_exceeding_capacity_raises():
    x = np.zeros(16)
    with pytest.raises(ValueError):
        qt_dwt_forward(x, levels=3)  # 4**3 = 64 > 16


def test_details_ordered_coarsest_first():
    """details[0] must correspond to the level closest to ll (coarsest),
    matching hcc.compression.grad_dwt._subband_slices' convention -- each
    successive level should have 4x more coefficients than the previous."""
    rng = np.random.default_rng(2)
    x = rng.standard_normal(4**3)
    ll, details = qt_dwt_forward(x, 3)
    sizes = [lh.size for lh, hl, hh in details]
    assert sizes == sorted(sizes)  # strictly increasing (coarsest = smallest, first)
    assert sizes[0] == ll.size
    assert sizes[-1] == x.size // 4


def test_single_level_matches_direct_formula():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    ll, details = qt_dwt_forward(x, 1)
    lh, hl, hh = details[0]
    a0, a1, a2, a3 = x
    assert ll[0] == pytest.approx((a0 + a1 + a2 + a3) * 0.5)
    assert lh[0] == pytest.approx((a0 + a1 - a2 - a3) * 0.5)
    assert hl[0] == pytest.approx((a0 - a1 + a2 - a3) * 0.5)
    assert hh[0] == pytest.approx((a0 - a1 - a2 + a3) * 0.5)
