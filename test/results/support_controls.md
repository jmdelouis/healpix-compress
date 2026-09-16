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

Same function sampled by both lattices: no content difference at all, so the ratio is the grid effect alone. Rows where the target is bracketed on both sides for fewer than 75% of regions are marked withheld: there the sweeps did not reach the target and the ratio is a quotient of bounds, not an estimate.

| beta | target | n | bracketed | median ratio | range | reportable |
|---|---|---|---|---|---|---|
| 1.9 | 0.0001 | 40 | 100% | 1.793 | [0.803, 2.223] | yes |
| 1.9 | 0.0003 | 40 | 100% | 1.650 | [0.821, 2.025] | yes |
| 1.9 | 0.001 | 40 | 95% | 1.438 | [0.811, 1.858] | yes |
| 1.9 | 0.003 | 40 | 95% | 1.142 | [0.810, 1.822] | yes |
| 1.9 | 0.01 | 40 | 78% | 1.114 | [0.936, 1.344] | yes |
