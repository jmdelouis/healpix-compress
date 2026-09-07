# Reconstruction error per transfer, against its own estimand

- regions: 40, noise: median 0.0000e+00 (0.0% of the field's std)

## By transfer and scene class

| kernel | estimand | class | n | e_resample | % of estimand std | blur error |
|---|---|---|---|---|---|---|
| `clough_tocher` | observed | agriculture | 10 | 4.5446e-04 | 1.3% | 1.8008e-02 |
| `clough_tocher` | observed | forest | 10 | 3.5820e-04 | 1.0% | 6.7250e-03 |
| `clough_tocher` | observed | urban | 10 | 1.1019e-03 | 1.8% | 8.8076e-02 |
| `clough_tocher` | observed | water | 10 | 1.3465e-04 | 1.5% | 4.9626e-03 |
| `psf` | latent | agriculture | 10 | 1.7383e-02 | 36.4% | 1.8008e-02 |
| `psf` | latent | forest | 10 | 6.3270e-03 | 24.0% | 6.7250e-03 |
| `psf` | latent | urban | 10 | 8.6551e-02 | 77.2% | 8.8076e-02 |
| `psf` | latent | water | 10 | 4.6710e-03 | 44.1% | 4.9626e-03 |

## Medians over all regions

| kernel | estimand | e_resample | % of estimand std | against the other estimand |
|---|---|---|---|---|
| `clough_tocher` | observed | 4.5411e-04 | 1.4% | 1.8224e-02 |
| `psf` | latent | 1.7383e-02 | 44.5% | 3.2000e-03 |

The last column is what each transfer would score if it were held to the other product's reference. A large gap between the two columns is the whole reason these products cannot be compared at the same absolute RMSE.

- `clough_tocher` against the observed field: 4.5411e-04
- `psf` against the latent field: 1.7383e-02
- blur alone (observed minus latent): 1.8008e-02; `psf` recovers 3.5% of it

These are the two numbers a compression target has to be read against, one per product. They are written to results/resampling_error.csv for `alpha_budget.py`.