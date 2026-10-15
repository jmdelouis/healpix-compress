# Reconstruction error per transfer, against its own estimand

- regions: 40, noise: median 0.0000e+00 (0.0% of the field's std)

## By transfer and scene class

| kernel | estimand | class | n | e_resample | % of estimand std | blur error |
|---|---|---|---|---|---|---|
| `clough_tocher` | observed | agriculture | 10 | 1.9089e-03 | 5.7% | 8.8204e-03 |
| `clough_tocher` | observed | forest | 10 | 1.0849e-03 | 3.9% | 2.1256e-03 |
| `clough_tocher` | observed | urban | 10 | 5.3992e-03 | 10.3% | 4.0441e-02 |
| `clough_tocher` | observed | water | 10 | 7.2177e-04 | 6.7% | 2.4350e-03 |
| `psf` | latent | agriculture | 10 | 7.5061e-03 | 17.6% | 8.8204e-03 |
| `psf` | latent | forest | 10 | 2.5265e-03 | 11.4% | 2.1256e-03 |
| `psf` | latent | urban | 10 | 3.3979e-02 | 45.1% | 4.0441e-02 |
| `psf` | latent | water | 10 | 2.6097e-03 | 24.0% | 2.4350e-03 |

## Medians over all regions

| kernel | estimand | e_resample | % of estimand std | against the other estimand |
|---|---|---|---|---|
| `clough_tocher` | observed | 1.8926e-03 | 6.6% | 9.7947e-03 |
| `psf` | latent | 7.5061e-03 | 23.0% | 6.5737e-03 |

The last column is what each transfer would score if it were held to the other product's reference. A large gap between the two columns is the whole reason these products cannot be compared at the same absolute RMSE.

- `clough_tocher` against the observed field: 1.8926e-03
- `psf` against the latent field: 7.5061e-03
- blur alone (observed minus latent): 8.8204e-03; `psf` recovers 14.9% of it

These are the two numbers a compression target has to be read against, one per product. They are written to results/resampling_error.csv for `alpha_budget.py`.