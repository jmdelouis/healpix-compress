# Transfer-kernel ablation (V3)

- regions: 3, kernels: bicubic, clough_tocher, linear, nearest, psf
- reference kernel: `clough_tocher` (the one the paper uses)

## Against an analytic truth (the metric that decides)

One analytic field, evaluated at the native pixel centres to play the observation, transferred to the cells by each kernel, and compared to the same field evaluated at the cell centres. Neither lattice is privileged, so a transfer gets no credit for merely copying samples.

| kernel | transfer RMSE / std | 99.9th pct / std | overshoot outside source range | cells kept |
|---|---|---|---|---|
| `bicubic` | 381.75% | 45.508 | 1.309% | all |
| `clough_tocher` | 1.19% | 0.050 | 0.000% | all |
| `linear` | 12.07% | 1.656 | 0.000% | all |
| `nearest` | 9.84% | 0.398 | 0.000% | all |
| `psf` | 8.78% | 1.145 | 0.005% | all |

## Against the imagery (read back through the same operator)

**These two columns reward a transfer for being close to the identity on the native lattice.** A nearest transfer copies pixel values and the nearest read hands most of them straight back, so the two nearly cancel; the log-PSD distance has the same weakness, since a copy reproduces the source spectrum exactly. Read them as a description of the transfer-and-read pair, not as a ranking of transfers.

| kernel | floor | floor / scene std | 99.9th pct / std | log-PSD RMS | top-octave power | cells outside source range | build (s) |
|---|---|---|---|---|---|---|---|
| `bicubic` | 1.763e-01 | 285.8% | 21.31 | 1.737 | 108.16 | 1.710% | 30 |
| `clough_tocher` | 4.090e-03 | 8.3% | 0.64 | 0.055 | 0.89 | 0.000% | nan |
| `linear` | 2.834e-03 | 6.9% | 0.63 | 0.043 | 0.90 | 0.948% | 3 |
| `nearest` | 9.254e-04 | 2.4% | 0.22 | 0.010 | 1.00 | 0.000% | nan |
| `psf` | 6.704e-03 | 18.9% | 1.43 | 0.402 | 2.40 | 0.097% | nan |

A top-octave ratio above 1 means the transfer leaves more high-frequency power than the observation carries, which is the failure the paper describes qualitatively; below 1 means it smooths.
- **`bicubic` diverges**: 1.7% of cells fall outside the range of the data it was given, and its floor reaches 286% of the scene's standard deviation. Note that it passes the constant-field check, which an interpolant exact on polynomials always will: that check is necessary and not sufficient.
- **`linear` diverges**: 1.1% of cells fall outside the range of the data it was given, and its floor reaches 4% of the scene's standard deviation. Note that it passes the constant-field check, which an interpolant exact on polynomials always will: that check is necessary and not sufficient.
- **`psf` diverges**: 2.3% of cells fall outside the range of the data it was given, and its floor reaches 19% of the scene's standard deviation. Note that it passes the constant-field check, which an interpolant exact on polynomials always will: that check is necessary and not sufficient.

## Paired against the reference kernel

| kernel | n | floor ratio | range |
|---|---|---|---|
| `bicubic` | 3 | 37.389 | [34.623, 44.473] |
| `clough_tocher` | 3 | 1.000 | [1.000, 1.000] |
| `linear` | 3 | 0.693 | [0.642, 0.840] |
| `nearest` | 3 | 0.232 | [0.226, 0.292] |
| `psf` | 3 | 1.698 | [1.639, 2.295] |

- lowest transfer error against the analytic truth: `psf` at 8.780e-02, against 1.190e-02 for `clough_tocher` (+637.6%)