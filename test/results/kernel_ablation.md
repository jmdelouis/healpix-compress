# Transfer-kernel ablation (V3)

- regions: 40, kernels: bicubic, clough_tocher, linear, nearest, psf
- reference kernel: `clough_tocher` (the one the paper uses)

## Against an analytic truth (the metric that decides)

One analytic field, evaluated at the native pixel centres to play the observation, transferred to the cells by each kernel, and compared to the same field evaluated at the cell centres. Neither lattice is privileged, so a transfer gets no credit for merely copying samples.

| kernel | n | transfer RMSE / std | 99.9th pct / std | overshoot outside source range | cells kept |
|---|---|---|---|---|---|
| `bicubic` | 39 | 634.28% | 77.759 | 2.037% | all |
| `clough_tocher` | 40 | 7.97% | 0.320 | 0.000% | all |
| `linear` | 39 | 16.20% | 1.576 | 0.000% | all |
| `nearest` | 40 | 18.30% | 0.729 | 0.000% | all |
| `psf` | 40 | 24.46% | 2.019 | 0.015% | all |

- **support is not common**: only 39 of 40 regions were measured for every kernel (missing: lake_taupo). The medians above are each over that kernel's own regions; the paired comparison below uses the common support only.

Restricted to the common support:

| kernel | n | transfer RMSE / std |
|---|---|---|
| `bicubic` | 39 | 634.28% |
| `clough_tocher` | 39 | 7.97% |
| `linear` | 39 | 16.20% |
| `nearest` | 39 | 18.30% |
| `psf` | 39 | 24.49% |

## Against the imagery (read back through the same operator)

**These two columns reward a transfer for being close to the identity on the native lattice.** A nearest transfer copies pixel values and the nearest read hands most of them straight back, so the two nearly cancel; the log-PSD distance has the same weakness, since a copy reproduces the source spectrum exactly. Read them as a description of the transfer-and-read pair, not as a ranking of transfers.

| kernel | floor | floor / scene std | 99.9th pct / std | log-PSD RMS | top-octave power | cells outside source range | build (s) |
|---|---|---|---|---|---|---|---|
| `bicubic` | 5.078e-03 | 11.9% | 0.87 | 0.162 | 1.51 | 2.454% | cached |
| `clough_tocher` | 4.038e-03 | 10.0% | 0.85 | 0.048 | 0.89 | 0.000% | cached |
| `linear` | 2.834e-03 | 6.9% | 0.74 | 0.030 | 0.92 | 0.854% | cached |
| `nearest` | 6.215e-05 | 0.7% | 0.00 | 0.001 | 1.00 | 0.000% | cached |
| `psf` | 7.835e-03 | 23.6% | 1.79 | 0.443 | 2.80 | 0.194% | cached |

A top-octave ratio above 1 means the transfer leaves more high-frequency power than the observation carries, which is the failure the paper describes qualitatively; below 1 means it smooths.
- `bicubic` overshoots mildly: 2.5% of cells outside the source range, but its 99.9th-percentile error stays at 0.87 times the scene's standard deviation, so this is non-monotonicity and not divergence.
- `linear` overshoots mildly: 0.9% of cells outside the source range, but its 99.9th-percentile error stays at 0.74 times the scene's standard deviation, so this is non-monotonicity and not divergence.

## Paired against the reference kernel

| kernel | n | floor ratio | range |
|---|---|---|---|
| `bicubic` | 39 | 0.957 | [0.687, 54.023] |
| `clough_tocher` | 40 | 1.000 | [1.000, 1.000] |
| `linear` | 39 | 0.742 | [0.576, 0.946] |
| `nearest` | 40 | 0.118 | [0.000, 0.361] |
| `psf` | 40 | 2.546 | [1.565, 4.774] |

## Compressed size at matched RMSE, paired against the reference kernel

| kernel | 0.0001 | 0.0003 | 0.001 | 0.003 | 0.01 |
|---|---|---|---|---|---|
| `bicubic` | 1.654 | 1.924 | 2.872 | 5.685 | 12.198 |
| `clough_tocher` | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| `linear` | 1.205 | 1.306 | 1.603 | 1.924 | 2.008 |
| `nearest` | 1.208 | 1.312 | 1.514 | 1.712 | 1.462 |
| `psf` | 1.267 | 1.366 | 1.642 | 2.078 | 2.392 |

Values are per-region ratios of compressed size to the reference kernel's, so 1.00 means the kernel changes nothing at that target.

- **`clough_tocher` has the lowest transfer error against the analytic truth**, 7.973e-02; the next best is `linear` at 1.620e-01, 2.0 times worse. The kernel the paper uses is the right one on this measure.