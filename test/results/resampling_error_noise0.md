# Reconstruction error per transfer, against its own estimand

- regions: 39, noise: median 0.0000e+00 (0.0% of the field's std)

## By transfer and scene class

| kernel | estimand | class | n | e_resample | % of estimand std | blur error |
|---|---|---|---|---|---|---|
| `clough_tocher` | observed | agriculture | 10 | 1.9089e-03 | 5.7% | 8.8204e-03 |
| `clough_tocher` | observed | forest | 10 | 1.3100e-03 | 3.9% | 2.8265e-03 |
| `clough_tocher` | observed | urban | 10 | 5.3992e-03 | 10.3% | 4.0441e-02 |
| `clough_tocher` | observed | water | 9 | 1.1064e-03 | 6.9% | 3.4545e-03 |
| `psf` | latent | agriculture | 10 | 7.5061e-03 | 17.6% | 8.8204e-03 |
| `psf` | latent | forest | 10 | 3.2897e-03 | 11.4% | 2.8265e-03 |
| `psf` | latent | urban | 10 | 3.3979e-02 | 45.1% | 4.0441e-02 |
| `psf` | latent | water | 9 | 4.0077e-03 | 24.5% | 3.4545e-03 |

## Medians over all regions

| kernel | estimand | e_resample | % of estimand std | against the other estimand |
|---|---|---|---|---|
| `clough_tocher` | observed | 1.9426e-03 | 6.7% | 1.0364e-02 |
| `psf` | latent | 7.6073e-03 | 23.0% | 6.6941e-03 |

The last column is what each transfer would score if it were held to the other product's reference. A large gap between the two columns is the whole reason these products cannot be compared at the same absolute RMSE.

- `clough_tocher` against the observed field: 1.9426e-03
- `psf` against the latent field: 7.6073e-03
- blur alone (observed minus latent): 9.3521e-03; `psf` recovers 18.7% of it

These are the two numbers a compression target has to be read against, one per product. They are written to results/resampling_error.csv for `alpha_budget.py`.