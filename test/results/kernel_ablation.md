# Transfer-kernel ablation (V3)

- regions: 39, kernels: bicubic
- reference kernel: `clough_tocher` (the one the paper uses)

## Against an analytic truth (the metric that decides)

One analytic field, evaluated at the native pixel centres to play the observation, transferred to the cells by each kernel, and compared to the same field evaluated at the cell centres. Neither lattice is privileged, so a transfer gets no credit for merely copying samples.

| kernel | n | transfer RMSE / std | 99.9th pct / std | overshoot outside source range | cells kept |
|---|---|---|---|---|---|
| `bicubic` | 39 | 634.28% | 77.759 | 2.037% | all |

Restricted to the common support:

| kernel | n | transfer RMSE / std |
|---|---|---|
| `bicubic` | 39 | 634.28% |

## Against the imagery (read back through the same operator)

**These two columns reward a transfer for being close to the identity on the native lattice.** A nearest transfer copies pixel values and the nearest read hands most of them straight back, so the two nearly cancel; the log-PSD distance has the same weakness, since a copy reproduces the source spectrum exactly. Read them as a description of the transfer-and-read pair, not as a ranking of transfers.

| kernel | floor | floor / scene std | 99.9th pct / std | log-PSD RMS | top-octave power | cells outside source range | build (s) |
|---|---|---|---|---|---|---|---|
| `bicubic` | 5.078e-03 | 11.9% | 0.87 | 0.162 | 1.51 | 2.454% | cached |

A top-octave ratio above 1 means the transfer leaves more high-frequency power than the observation carries, which is the failure the paper describes qualitatively; below 1 means it smooths.
- `bicubic` overshoots mildly: 2.5% of cells outside the source range, but its 99.9th-percentile error stays at 0.87 times the scene's standard deviation, so this is non-monotonicity and not divergence.

## Paired against the reference kernel

| kernel | n | floor ratio | range |
|---|---|---|---|

## Compressed size at matched RMSE, paired against the reference kernel

| kernel | 0.0001 | 0.0003 | 0.001 | 0.003 | 0.01 |
|---|---|---|---|---|---|
| `bicubic` | -- | -- | -- | -- | -- |

Values are per-region ratios of compressed size to the reference kernel's, so 1.00 means the kernel changes nothing at that target.
