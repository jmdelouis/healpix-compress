# Support-mismatch controls (V1)

## Shifted-crop ensemble

How far the ratio moves when the UTM crop is re-cut at offsets of up to 16 native pixels, with the HEALPix side fixed. Spread is across offsets, within a region.

| target | regions | offsets | median ratio | median spread (max-min) | spread / (ratio-1) |
|---|---|---|---|---|---|
| 0.0001 | 40 | 9 | 2.074 | 0.012 | 0.01 |
| 0.0003 | 40 | 9 | 1.952 | 0.015 | 0.02 |
| 0.001 | 40 | 9 | 1.591 | 0.022 | 0.04 |
| 0.003 | 40 | 9 | 1.221 | 0.040 | 0.18 |
| 0.01 | 40 | 9 | 0.968 | 0.057 | 1.79 |

## Common central support

Largest common-ground region on each native lattice, normalised by the ground area each covers.

- UTM rectangle: 417x391 px median (19.74 km2)
- HEALPix block: 256 cells per side (2.53 km2)

| target | regions | HPX B/km2 | UTM B/km2 | ratio |
|---|---|---|---|---|
| 0.0001 | 40 | 1.61e+04 | 7.36e+03 | 2.103 |
| 0.0003 | 40 | 1.08e+04 | 5.24e+03 | 1.973 |
| 0.001 | 40 | 5.16e+03 | 3.06e+03 | 1.617 |
| 0.003 | 40 | 1.89e+03 | 1.27e+03 | 1.569 |
| 0.01 | 40 | 435 | 258 | 1.486 |

## One analytic field on both geometries

Same function sampled by both lattices: no content difference at all, so the ratio is the grid effect alone. Rows where the target is bracketed on both sides for fewer than 75% of regions are marked withheld: there the sweeps did not reach the target and the ratio is a quotient of bounds, not an estimate. `transfer` = none is the field evaluated at the cell centres; otherwise the field was sampled on the native lattice and carried to the cells by that resampler, the one step of the real chain the direct rows omit.

| transfer | beta | target | n | bracketed | median ratio | range | reportable |
|---|---|---|---|---|---|---|---|
| clough_tocher | 3 | 0.0001 | 40 | 100% | 1.953 | [0.925, 2.214] | yes |
| clough_tocher | 3 | 0.0003 | 40 | 98% | 1.611 | [0.745, 2.085] | yes |
| clough_tocher | 3 | 0.001 | 40 | 90% | 1.025 | [0.740, 1.856] | yes |
| clough_tocher | 3 | 0.003 | 40 | 88% | 1.028 | [0.753, 2.478] | yes |
| clough_tocher | 3 | 0.01 | 40 | 60% | 1.027 | [0.968, 1.327] | **withheld** |
| clough_tocher | 4 | 0.0001 | 40 | 100% | 1.111 | [0.825, 1.841] | yes |
| clough_tocher | 4 | 0.0003 | 40 | 90% | 1.098 | [0.854, 1.372] | yes |
| clough_tocher | 4 | 0.001 | 40 | 88% | 1.113 | [0.969, 1.589] | yes |
| clough_tocher | 4 | 0.003 | 40 | 72% | 1.117 | [0.979, 2.075] | **withheld** |
| clough_tocher | 4 | 0.01 | 40 | 18% | 1.299 | [1.060, 1.977] | **withheld** |
| clough_tocher | 5 | 0.0001 | 40 | 95% | 1.159 | [0.874, 1.363] | yes |
| clough_tocher | 5 | 0.0003 | 40 | 90% | 1.204 | [1.031, 1.442] | yes |
| clough_tocher | 5 | 0.001 | 40 | 82% | 1.246 | [1.123, 1.750] | yes |
| clough_tocher | 5 | 0.003 | 40 | 55% | 1.329 | [1.222, 1.730] | **withheld** |
| clough_tocher | 5 | 0.01 | 40 | 30% | 1.661 | [1.333, 1.868] | **withheld** |
| none | 3 | 0.0001 | 40 | 100% | 1.793 | [1.018, 2.241] | yes |
| none | 3 | 0.0003 | 40 | 98% | 1.545 | [0.837, 1.998] | yes |
| none | 3 | 0.001 | 40 | 90% | 1.127 | [0.834, 1.751] | yes |
| none | 3 | 0.003 | 40 | 85% | 1.109 | [0.842, 2.478] | yes |
| none | 3 | 0.01 | 40 | 60% | 1.049 | [0.979, 1.405] | **withheld** |
| none | 4 | 0.0001 | 40 | 100% | 1.175 | [0.886, 1.804] | yes |
| none | 4 | 0.0003 | 40 | 90% | 1.165 | [0.926, 1.403] | yes |
| none | 4 | 0.001 | 40 | 88% | 1.114 | [0.979, 1.589] | yes |
| none | 4 | 0.003 | 40 | 72% | 1.123 | [0.990, 2.075] | **withheld** |
| none | 4 | 0.01 | 40 | 18% | 1.300 | [1.059, 1.977] | **withheld** |
| none | 5 | 0.0001 | 40 | 95% | 1.160 | [0.898, 1.369] | yes |
| none | 5 | 0.0003 | 40 | 90% | 1.208 | [1.033, 1.442] | yes |
| none | 5 | 0.001 | 40 | 82% | 1.247 | [1.124, 1.750] | yes |
| none | 5 | 0.003 | 40 | 55% | 1.329 | [1.221, 1.730] | **withheld** |
| none | 5 | 0.01 | 40 | 30% | 1.661 | [1.335, 1.868] | **withheld** |

Transfer error on the analytic field (RMSE against the exact cell values, in units of the scene std): median 0.99%, range [0.94%, 1.06%] over 40 regions.

## The transfer alone, on the native lattice

The native crop read through the same Clough-Tocher operator at points shifted by a fraction of a pixel, at the native density (operator: healpix_resample). Smoothed/native is what the transfer alone buys at each target; HPX/smoothed is what remains of the carrier's cost once both sides carry the same smoothing. Sizes are per sample (the shifted lattice has one row and one column fewer).

- std ratio smoothed/native: median 0.9868

| phase, target | n | smoothed/native | range | n | HPX/smoothed | range | HPX/native |
|---|---|---|---|---|---|---|---|
| (0.25,0.25) 0.0001 | 40 | 0.957 | [0.853, 0.986] | 40 | 2.187 | [0.977, 2.347] | 2.072 |
| (0.25,0.25) 0.0003 | 40 | 0.945 | [0.579, 0.983] | 40 | 2.058 | [0.779, 2.320] | 1.954 |
| (0.25,0.25) 0.001 | 40 | 0.918 | [0.285, 1.050] | 39 | 1.780 | [0.458, 3.563] | 1.578 |
| (0.25,0.25) 0.003 | 39 | 0.888 | [0.238, 0.964] | 37 | 1.356 | [0.261, 2.398] | 1.226 |
| (0.25,0.25) 0.01 | 35 | 0.863 | [0.790, 0.994] | 32 | 1.101 | [0.316, 1.500] | 0.933 |
| (0.5,0.5) 0.0001 | 40 | 0.944 | [0.855, 0.967] | 40 | 2.223 | [0.971, 2.394] | 2.072 |
| (0.5,0.5) 0.0003 | 40 | 0.919 | [0.513, 0.957] | 40 | 2.105 | [0.780, 2.402] | 1.954 |
| (0.5,0.5) 0.001 | 40 | 0.886 | [0.037, 1.002] | 39 | 1.828 | [0.471, 17.974] | 1.578 |
| (0.5,0.5) 0.003 | 39 | 0.845 | [0.205, 1.030] | 37 | 1.484 | [0.283, 2.246] | 1.226 |
| (0.5,0.5) 0.01 | 35 | 0.778 | [0.638, 0.992] | 32 | 1.237 | [0.327, 1.589] | 0.933 |
