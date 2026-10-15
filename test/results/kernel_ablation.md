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
| `bicubic` | 5.507e-03 | 11.9% | 0.90 | 0.162 | 1.51 | 2.676% | cached |
| `clough_tocher` | 4.038e-03 | 10.1% | 0.89 | 0.049 | 0.89 | 0.001% | cached |
| `linear` | 2.834e-03 | 7.2% | 0.79 | 0.030 | 0.92 | 0.854% | cached |
| `nearest` | 6.215e-05 | 0.9% | 0.00 | 0.001 | 1.00 | 0.000% | cached |
| `psf` | 7.809e-03 | 24.4% | 1.88 | 0.443 | 2.80 | 0.204% | cached |

A top-octave ratio above 1 means the transfer leaves more high-frequency power than the observation carries, which is the failure the paper describes qualitatively; below 1 means it smooths.
- `bicubic` overshoots mildly: 2.7% of cells outside the source range, but its 99.9th-percentile error stays at 0.90 times the scene's standard deviation, so this is non-monotonicity and not divergence.
- `linear` overshoots mildly: 0.9% of cells outside the source range, but its 99.9th-percentile error stays at 0.79 times the scene's standard deviation, so this is non-monotonicity and not divergence.

## Paired against the reference kernel

| kernel | n | floor ratio | range |
|---|---|---|---|
| `bicubic` | 39 | 0.957 | [0.687, 49.772] |
| `clough_tocher` | 40 | 1.000 | [1.000, 1.000] |
| `linear` | 39 | 0.741 | [0.576, 0.881] |
| `nearest` | 40 | 0.107 | [0.000, 0.361] |
| `psf` | 40 | 2.546 | [1.565, 4.774] |

## Compressed size at matched RMSE, paired against the reference kernel

| kernel | 0.0001 | 0.0003 | 0.001 | 0.003 | 0.01 |
|---|---|---|---|---|---|
| `bicubic` | 1.647 | 1.915 | 2.870 | 5.700 | 12.211 |
| `clough_tocher` | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| `linear` | 1.203 | 1.306 | 1.603 | 1.924 | 1.976 |
| `nearest` | 1.208 | 1.312 | 1.499 | 1.674 | 1.439 |
| `psf` | 1.265 | 1.349 | 1.639 | 2.041 | 2.336 |

Values are per-region ratios of compressed size to the reference kernel's, so 1.00 means the kernel changes nothing at that target.

- **`clough_tocher` has the lowest transfer error against the analytic truth**, 7.973e-02; the next best is `linear` at 1.620e-01, 2.0 times worse. The kernel the paper uses is the right one on this measure.