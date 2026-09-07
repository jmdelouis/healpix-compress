# Transfer-kernel ablation (V3)

- regions: 3, kernels: bicubic, clough_tocher, linear, nearest, psf
- reference kernel: `clough_tocher` (the one the paper uses)

## Against an analytic truth (the metric that decides)

One analytic field, evaluated at the native pixel centres to play the observation, transferred to the cells by each kernel, and compared to the same field evaluated at the cell centres. Neither lattice is privileged, so a transfer gets no credit for merely copying samples.

| kernel | n | transfer RMSE / std | 99.9th pct / std | overshoot outside source range | cells kept |
|---|---|---|---|---|---|
| `bicubic` | 2 | 381.75% | 45.508 | 1.309% | all |
| `clough_tocher` | 2 | 1.19% | 0.050 | 0.000% | all |
| `linear` | 2 | 12.07% | 1.656 | 0.000% | all |
| `nearest` | 3 | 9.84% | 0.398 | 0.000% | all |
| `psf` | 2 | 8.78% | 1.145 | 0.005% | all |

- **support is not common**: only 2 of 3 regions were measured for every kernel (missing: iowa). The medians above are each over that kernel's own regions; the paired comparison below uses the common support only.

Restricted to the common support:

| kernel | n | transfer RMSE / std |
|---|---|---|
| `bicubic` | 2 | 381.75% |
| `clough_tocher` | 2 | 1.19% |
| `linear` | 2 | 12.07% |
| `nearest` | 2 | 9.81% |
| `psf` | 2 | 8.78% |

## Against the imagery (read back through the same operator)

**These two columns reward a transfer for being close to the identity on the native lattice.** A nearest transfer copies pixel values and the nearest read hands most of them straight back, so the two nearly cancel; the log-PSD distance has the same weakness, since a copy reproduces the source spectrum exactly. Read them as a description of the transfer-and-read pair, not as a ranking of transfers.

| kernel | floor | floor / scene std | 99.9th pct / std | log-PSD RMS | top-octave power | cells outside source range | build (s) |
|---|---|---|---|---|---|---|---|
| `bicubic` | 1.763e-01 | 285.8% | 21.31 | 1.737 | 108.16 | 1.710% | 30 |
| `clough_tocher` | 4.090e-03 | 8.3% | 0.64 | 0.055 | 0.89 | 0.000% | cached |
| `linear` | 2.834e-03 | 6.9% | 0.63 | 0.043 | 0.90 | 0.948% | 3 |
| `nearest` | 9.254e-04 | 2.4% | 0.22 | 0.010 | 1.00 | 0.000% | cached |
| `psf` | 6.704e-03 | 18.9% | 1.43 | 0.402 | 2.40 | 0.097% | cached |

A top-octave ratio above 1 means the transfer leaves more high-frequency power than the observation carries, which is the failure the paper describes qualitatively; below 1 means it smooths.
- **`bicubic` diverges**: its 99.9th-percentile error reaches 21 times the scene's standard deviation and its floor 286% of it, with 1.7% of cells outside the range of the data it was given (median over all 3 regions). It passes the constant-field check, which any interpolant exact on polynomials will: that check is necessary and not sufficient.
- `linear` overshoots mildly: 0.9% of cells outside the source range, but its 99.9th-percentile error stays at 0.63 times the scene's standard deviation, so this is non-monotonicity and not divergence.

## Paired against the reference kernel

| kernel | n | floor ratio | range |
|---|---|---|---|
| `bicubic` | 3 | 37.389 | [34.623, 44.473] |
| `clough_tocher` | 3 | 1.000 | [1.000, 1.000] |
| `linear` | 3 | 0.693 | [0.642, 0.840] |
| `nearest` | 3 | 0.232 | [0.226, 0.292] |
| `psf` | 3 | 1.698 | [1.639, 2.295] |

- **`clough_tocher` has the lowest transfer error against the analytic truth**, 1.190e-02; the next best is `psf` at 8.780e-02, 7.4 times worse. The kernel the paper uses is the right one on this measure.