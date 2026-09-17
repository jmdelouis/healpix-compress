# Support-mismatch controls (V1)

## Shifted-crop ensemble

How far the ratio moves when the UTM crop is re-cut at offsets of up to 16 native pixels, with the HEALPix side fixed. Spread is across offsets, within a region.

| target | regions | offsets | median ratio | median spread (max-min) | spread / (ratio-1) |
|---|---|---|---|---|---|
| 0.0001 | 4 | 9 | 2.123 | 0.008 | 0.01 |
| 0.0003 | 4 | 9 | 1.989 | 0.011 | 0.01 |
| 0.001 | 4 | 9 | 1.582 | 0.018 | 0.03 |
| 0.003 | 4 | 9 | 1.159 | 0.020 | 0.13 |
| 0.01 | 4 | 9 | 1.278 | 0.030 | 0.11 |

## Common central support

Largest common-ground region on each native lattice, normalised by the ground area each covers.

- UTM rectangle: 405x401 px median (17.55 km2)
- HEALPix block: 256 cells per side (2.53 km2)

| target | regions | HPX B/km2 | UTM B/km2 | ratio |
|---|---|---|---|---|
| 0.0001 | 10 | 1.67e+04 | 7.69e+03 | 2.134 |
| 0.0003 | 10 | 1.13e+04 | 5.6e+03 | 1.975 |
| 0.001 | 10 | 5.58e+03 | 3.41e+03 | 1.617 |
| 0.003 | 10 | 1.89e+03 | 1.55e+03 | 1.378 |
| 0.01 | 10 | 512 | 270 | 1.355 |

## One analytic field on both geometries

Same function sampled by both lattices: no content difference at all, so the ratio is the grid effect alone. Rows where the target is bracketed on both sides for fewer than 75% of regions are marked withheld: there the sweeps did not reach the target and the ratio is a quotient of bounds, not an estimate. `transfer` = none is the field evaluated at the cell centres; otherwise the field was sampled on the native lattice and carried to the cells by that resampler, the one step of the real chain the direct rows omit.

| transfer | beta | target | n | bracketed | median ratio | range | reportable |
|---|---|---|---|---|---|---|---|
| clough_tocher | 1.9 | 0.0001 | 40 | 100% | 2.173 | [0.927, 2.313] | yes |
| clough_tocher | 1.9 | 0.0003 | 40 | 100% | 2.099 | [0.817, 2.238] | yes |
| clough_tocher | 1.9 | 0.001 | 40 | 98% | 1.900 | [0.751, 2.211] | yes |
| clough_tocher | 1.9 | 0.003 | 40 | 95% | 1.387 | [0.758, 2.370] | yes |
| clough_tocher | 1.9 | 0.01 | 40 | 80% | 0.944 | [0.817, 1.648] | yes |
| none | 1.9 | 0.0001 | 40 | 100% | 2.167 | [1.095, 2.371] | yes |
| none | 1.9 | 0.0003 | 40 | 100% | 2.076 | [0.989, 2.246] | yes |
| none | 1.9 | 0.001 | 40 | 98% | 1.885 | [0.911, 2.203] | yes |
| none | 1.9 | 0.003 | 40 | 95% | 1.449 | [0.927, 2.391] | yes |
| none | 1.9 | 0.01 | 40 | 80% | 1.155 | [0.973, 1.691] | yes |

Transfer error on the analytic field (RMSE against the exact cell values, in units of the scene std): median 8.80%, range [8.66%, 9.26%] over 40 regions.

## The transfer alone, on the native lattice

The native crop read through the same Clough-Tocher operator at points shifted by a fraction of a pixel, at the native density (operator: healpix_resample). Smoothed/native is what the transfer alone buys at each target; HPX/smoothed is what remains of the carrier's cost once both sides carry the same smoothing. Sizes are per sample (the shifted lattice has one row and one column fewer).

- std ratio smoothed/native: median 0.9854

| phase, target | n | smoothed/native | range | n | HPX/smoothed | range | HPX/native |
|---|---|---|---|---|---|---|---|
| (0.5,0.5) 0.0001 | 40 | 0.943 | [0.855, 0.967] | 40 | 2.228 | [0.971, 2.394] | 2.063 |
| (0.5,0.5) 0.0003 | 40 | 0.916 | [0.512, 0.957] | 40 | 2.110 | [0.780, 2.403] | 1.936 |
| (0.5,0.5) 0.001 | 40 | 0.890 | [0.037, 1.002] | 39 | 1.837 | [0.471, 17.974] | 1.595 |
| (0.5,0.5) 0.003 | 39 | 0.842 | [0.205, 1.039] | 37 | 1.470 | [0.283, 2.896] | 1.195 |
| (0.5,0.5) 0.01 | 34 | 0.778 | [0.638, 0.992] | 32 | 1.237 | [0.327, 1.589] | 0.933 |
