# Transfer-kernel ablation (V3)

- regions: 39, kernels: bicubic, clough_tocher, linear, nearest, psf
- reference kernel: `clough_tocher` (the one the paper uses)

## Against an analytic truth (the metric that decides)

One analytic field, evaluated at the native pixel centres to play the observation, transferred to the cells by each kernel, and compared to the same field evaluated at the cell centres. Neither lattice is privileged, so a transfer gets no credit for merely copying samples.

| kernel | n | transfer RMSE / std | 99.9th pct / std | overshoot outside source range | cells kept |
|---|---|---|---|---|---|
| `bicubic` | 38 | 392.34% | 47.585 | 1.307% | all |
| `clough_tocher` | 39 | 1.20% | 0.050 | 0.000% | all |
| `linear` | 38 | 11.48% | 1.584 | 0.000% | all |
| `nearest` | 39 | 10.04% | 0.405 | 0.000% | all |
| `psf` | 39 | 8.15% | 1.064 | 0.002% | all |

- **support is not common**: only 38 of 39 regions were measured for every kernel (missing: lake_taupo). The medians above are each over that kernel's own regions; the paired comparison below uses the common support only.

Restricted to the common support:

| kernel | n | transfer RMSE / std |
|---|---|---|
| `bicubic` | 38 | 392.34% |
| `clough_tocher` | 38 | 1.20% |
| `linear` | 38 | 11.48% |
| `nearest` | 38 | 10.04% |
| `psf` | 38 | 8.15% |

## Against the imagery (read back through the same operator)

**These two columns reward a transfer for being close to the identity on the native lattice.** A nearest transfer copies pixel values and the nearest read hands most of them straight back, so the two nearly cancel; the log-PSD distance has the same weakness, since a copy reproduces the source spectrum exactly. Read them as a description of the transfer-and-read pair, not as a ranking of transfers.

| kernel | floor | floor / scene std | 99.9th pct / std | log-PSD RMS | top-octave power | cells outside source range | build (s) |
|---|---|---|---|---|---|---|---|
| `bicubic` | 4.677e-03 | 11.7% | 0.83 | 0.162 | 1.51 | 2.565% | cached |
| `clough_tocher` | 4.090e-03 | 10.0% | 0.86 | 0.048 | 0.89 | 0.000% | cached |
| `linear` | 2.924e-03 | 7.0% | 0.76 | 0.030 | 0.92 | 0.857% | cached |
| `nearest` | 8.233e-05 | 0.8% | 0.00 | 0.001 | 1.00 | 0.000% | cached |
| `psf` | 8.008e-03 | 24.0% | 1.82 | 0.443 | 2.80 | 0.183% | cached |

A top-octave ratio above 1 means the transfer leaves more high-frequency power than the observation carries, which is the failure the paper describes qualitatively; below 1 means it smooths.
- `bicubic` overshoots mildly: 2.6% of cells outside the source range, but its 99.9th-percentile error stays at 0.83 times the scene's standard deviation, so this is non-monotonicity and not divergence.
- `linear` overshoots mildly: 0.9% of cells outside the source range, but its 99.9th-percentile error stays at 0.76 times the scene's standard deviation, so this is non-monotonicity and not divergence.

## Paired against the reference kernel

| kernel | n | floor ratio | range |
|---|---|---|---|
| `bicubic` | 38 | 0.947 | [0.687, 49.772] |
| `clough_tocher` | 39 | 1.000 | [1.000, 1.000] |
| `linear` | 38 | 0.741 | [0.576, 0.868] |
| `nearest` | 39 | 0.106 | [0.000, 0.361] |
| `psf` | 39 | 2.538 | [1.565, 4.774] |

## Compressed size at matched RMSE, paired against the reference kernel

| kernel | 0.0001 | 0.0003 | 0.001 | 0.003 | 0.01 |
|---|---|---|---|---|---|
| `bicubic` | 1.653 | 1.920 | 2.838 | 5.583 | 12.117 |
| `clough_tocher` | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| `linear` | 1.204 | 1.306 | 1.601 | 1.937 | 2.006 |
| `nearest` | 1.208 | 1.311 | 1.498 | 1.748 | 1.458 |
| `psf` | 1.265 | 1.360 | 1.639 | 2.108 | 2.376 |

Values are per-region ratios of compressed size to the reference kernel's, so 1.00 means the kernel changes nothing at that target.

- **`clough_tocher` has the lowest transfer error against the analytic truth**, 1.203e-02; the next best is `psf` at 8.153e-02, 6.8 times worse. The kernel the paper uses is the right one on this measure.