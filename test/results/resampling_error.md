# Reconstruction error per transfer, against its own estimand

- regions: 39, noise: median 2.3061e-03 (6.6% of the field's std)

## By transfer and scene class

| kernel | estimand | class | n | e_resample | % of estimand std | blur error |
|---|---|---|---|---|---|---|
| `clough_tocher` | observed | agriculture | 10 | 1.9444e-03 | 5.1% | 1.5717e-02 |
| `clough_tocher` | observed | forest | 10 | 1.3159e-03 | 5.0% | 5.6887e-03 |
| `clough_tocher` | observed | urban | 10 | 6.2967e-03 | 11.0% | 7.4814e-02 |
| `clough_tocher` | observed | water | 9 | 8.2038e-04 | 6.8% | 6.8073e-03 |
| `psf` | latent | agriculture | 10 | 1.5902e-02 | 33.2% | 1.5717e-02 |
| `psf` | latent | forest | 10 | 6.0265e-03 | 23.1% | 5.6887e-03 |
| `psf` | latent | urban | 10 | 7.5831e-02 | 74.4% | 7.4814e-02 |
| `psf` | latent | water | 9 | 6.4011e-03 | 55.2% | 6.8073e-03 |

## Medians over all regions

| kernel | estimand | e_resample | % of estimand std | against the other estimand |
|---|---|---|---|---|
| `clough_tocher` | observed | 1.9263e-03 | 5.3% | 1.6722e-02 |
| `psf` | latent | 1.7179e-02 | 40.8% | 5.4127e-03 |

The last column is what each transfer would score if it were held to the other product's reference. A large gap between the two columns is the whole reason these products cannot be compared at the same absolute RMSE.

- `clough_tocher` against the observed field: 1.9263e-03
- `psf` against the latent field: 1.7179e-02
- blur alone (observed minus latent): 1.6430e-02; `psf` recovers -4.6% of it

These are the two numbers a compression target has to be read against, one per product. They are written to results/resampling_error.csv for `alpha_budget.py`.