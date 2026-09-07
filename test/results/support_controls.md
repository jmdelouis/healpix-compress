# Support-mismatch controls (V1)

## What the real imagery's spectrum says

Radial power spectrum of each UTM crop, fitted as a power law times a Gaussian instrument response. The fitted FWHM is what the synthetic control should use: the nominal 12.5 m is a specification, this is what the data show.

- fitted slope beta: median 1.90, range [-0.19, 3.33]
- fitted PSF FWHM: median 11.9 m, range [4.3, 14.3] m
- variance above half the UTM Nyquist: median 0.21%

The same estimator is applied to the synthetic field, and it carries a small known bias from the Hann window and the radial binning (about +0.2 in slope and -1 m in FWHM on a field of known parameters). Matching fitted-to-fitted cancels it to first order, which is what feeding these medians straight back in does:

    python support_controls.py --control synthetic --betas 1.9 --psf-fwhm 11.9

## One analytic field on both geometries

Same function sampled by both lattices: no content difference at all, so the ratio is the grid effect alone. Rows where the target is bracketed on both sides for fewer than 75% of regions are marked withheld: there the sweeps did not reach the target and the ratio is a quotient of bounds, not an estimate.

| beta | target | n | bracketed | median ratio | range | reportable |
|---|---|---|---|---|---|---|
| 1.9 | 0.0001 | 40 | 100% | 2.221 | [1.115, 2.363] | yes |
| 1.9 | 0.0003 | 40 | 100% | 2.151 | [1.001, 2.295] | yes |
| 1.9 | 0.001 | 40 | 100% | 1.899 | [0.602, 2.258] | yes |
| 1.9 | 0.003 | 40 | 92% | 1.494 | [0.959, 2.134] | yes |
| 1.9 | 0.01 | 40 | 90% | 1.183 | [0.396, 1.941] | yes |
