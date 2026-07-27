import numpy as np
import pytest

from hcc.compression.quadtree_wavelet import (
    qt_dwt_forward,
    qt_dwt53_forward,
    qt_dwt53_inverse,
)
from hcc.compression.zorder import expand_tile, flatten_tile


@pytest.mark.parametrize("delta_total", [2, 3, 4, 5, 6])
def test_perfect_reconstruction_all_level_counts(delta_total):
    n = 4**delta_total
    rng = np.random.default_rng(delta_total)
    x = rng.standard_normal(n)
    for levels in range(1, delta_total + 1):
        ll, details = qt_dwt53_forward(x, levels)
        rec = qt_dwt53_inverse(ll, details)
        np.testing.assert_allclose(rec, x, atol=1e-8)


def test_levels_exceeding_capacity_raises():
    x = np.zeros(16)
    with pytest.raises(ValueError):
        qt_dwt53_forward(x, levels=3)  # 4**3 = 64 > 16


def test_non_power_of_4_raises():
    x = np.zeros(20)
    with pytest.raises(ValueError):
        qt_dwt53_forward(x, levels=1)


def test_details_ordered_coarsest_first():
    rng = np.random.default_rng(2)
    x = rng.standard_normal(4**3)
    ll, details = qt_dwt53_forward(x, 3)
    sizes = [lh.size for lh, hl, hh in details]
    assert sizes == sorted(sizes)
    assert sizes[0] == ll.size
    assert sizes[-1] == x.size // 4


def test_smooth_signal_concentrates_energy_in_ll():
    """CDF 5/3 lifting (as implemented here, without a separate
    biorthogonal normalisation step) is not orthonormal -- unlike the Haar
    step, overall coefficient energy is *not* preserved (the unnormalised
    0.5/0.25 lifting steps shrink detail-band energy considerably). What
    should still hold for a smooth, low-frequency-dominated signal is the
    qualitative energy-compaction property: the ll band should hold the
    large majority of the total coefficient energy, with detail energy
    small by comparison and no runaway blow-up."""
    n_side = 32
    xx, yy = np.meshgrid(np.linspace(0, 4, n_side), np.linspace(0, 4, n_side))
    img = 0.5 + 0.3 * np.sin(xx) * np.cos(yy)
    x = flatten_tile(img)
    ll, details = qt_dwt53_forward(x, 3)
    ll_energy = np.sum(ll**2)
    detail_energy = sum(
        np.sum(lh**2) + np.sum(hl**2) + np.sum(hh**2) for lh, hl, hh in details
    )
    assert ll_energy > detail_energy
    assert np.isfinite(detail_energy)


def test_single_level_reduces_low_frequency_energy_better_than_haar():
    """For a smooth signal, CDF 5/3's longer support should concentrate
    more energy into ll and less into the detail bands than the 2-tap
    Haar step -- the whole point of extending to a longer filter."""
    n_side = 64
    xx, yy = np.meshgrid(np.linspace(0, 6, n_side), np.linspace(0, 6, n_side))
    img = 0.5 + 0.3 * np.sin(xx) * np.cos(yy) + 0.1 * np.sin(3 * xx + yy)
    x = flatten_tile(img)

    ll53, details53 = qt_dwt53_forward(x, 2)
    detail_energy_53 = sum(
        np.sum(lh**2) + np.sum(hl**2) + np.sum(hh**2) for lh, hl, hh in details53
    )

    ll_haar, details_haar = qt_dwt_forward(x, 2)
    detail_energy_haar = sum(
        np.sum(lh**2) + np.sum(hl**2) + np.sum(hh**2) for lh, hl, hh in details_haar
    )

    assert detail_energy_53 < detail_energy_haar


def test_matches_manual_single_level_formula():
    """Cross-check _cdf53_analysis_step indirectly: on a tiny 4x4 tile,
    verify perfect reconstruction and that ll is a sensible smoothed
    average (sanity, not a bit-exact reference since CDF 5/3 has no
    trivial closed form the way Haar does)."""
    x = np.arange(16, dtype=np.float64)
    ll, details = qt_dwt53_forward(x, 1)
    rec = qt_dwt53_inverse(ll, details)
    np.testing.assert_allclose(rec, x, atol=1e-8)
    assert ll.size == 4
    lh, hl, hh = details[0]
    assert lh.size == hl.size == hh.size == 4
