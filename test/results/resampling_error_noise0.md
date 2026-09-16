# Reconstruction error per transfer, against its own estimand

- regions: 39, noise: median 0.0000e+00 (0.0% of the field's std)

## By transfer and scene class

| kernel | estimand | class | n | e_resample | % of estimand std | blur error |
|---|---|---|---|---|---|---|
| `clough_tocher` | observed | agriculture | 10 | 3.8201e-04 | 1.1% | 1.5717e-02 |
| `clough_tocher` | observed | forest | 10 | 2.8891e-04 | 0.8% | 5.6887e-03 |
| `clough_tocher` | observed | urban | 10 | 8.9860e-04 | 1.5% | 7.4814e-02 |
| `clough_tocher` | observed | water | 9 | 1.8631e-04 | 1.4% | 6.8073e-03 |
| `psf` | latent | agriculture | 10 | 1.4721e-02 | 31.3% | 1.5717e-02 |
| `psf` | latent | forest | 10 | 5.3170e-03 | 20.6% | 5.6887e-03 |
| `psf` | latent | urban | 10 | 7.2556e-02 | 71.8% | 7.4814e-02 |
| `psf` | latent | water | 9 | 6.3543e-03 | 41.9% | 6.8073e-03 |

## Medians over all regions

| kernel | estimand | e_resample | % of estimand std | against the other estimand |
|---|---|---|---|---|
| `clough_tocher` | observed | 3.9380e-04 | 1.3% | 1.6496e-02 |
| `psf` | latent | 1.5422e-02 | 39.5% | 2.1100e-03 |

The last column is what each transfer would score if it were held to the other product's reference. A large gap between the two columns is the whole reason these products cannot be compared at the same absolute RMSE.

- `clough_tocher` against the observed field: 3.9380e-04
- `psf` against the latent field: 1.5422e-02
- blur alone (observed minus latent): 1.6430e-02; `psf` recovers 6.1% of it

These are the two numbers a compression target has to be read against, one per product. They are written to results/resampling_error.csv for `alpha_budget.py`.