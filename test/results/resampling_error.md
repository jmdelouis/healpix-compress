# Reconstruction error per transfer, against its own estimand

- regions: 40, noise: median 2.2719e-03 (6.8% of the field's std)

## By transfer and scene class

| kernel | estimand | class | n | e_resample | % of estimand std | blur error |
|---|---|---|---|---|---|---|
| `clough_tocher` | observed | agriculture | 10 | 1.9444e-03 | 5.1% | 1.5717e-02 |
| `clough_tocher` | observed | forest | 10 | 1.3159e-03 | 5.0% | 5.6887e-03 |
| `clough_tocher` | observed | urban | 10 | 6.2967e-03 | 11.0% | 7.4814e-02 |
| `clough_tocher` | observed | water | 10 | 6.9155e-04 | 9.8% | 4.2650e-03 |
| `psf` | latent | agriculture | 10 | 1.5902e-02 | 33.2% | 1.5717e-02 |
| `psf` | latent | forest | 10 | 6.0265e-03 | 23.1% | 5.6887e-03 |
| `psf` | latent | urban | 10 | 7.5831e-02 | 74.4% | 7.4814e-02 |
| `psf` | latent | water | 10 | 4.5914e-03 | 54.8% | 4.2650e-03 |

## Medians over all regions

| kernel | estimand | e_resample | % of estimand std | against the other estimand |
|---|---|---|---|---|
| `clough_tocher` | observed | 1.9140e-03 | 5.5% | 1.5946e-02 |
| `psf` | latent | 1.5902e-02 | 41.4% | 5.2959e-03 |

The last column is what each transfer would score if it were held to the other product's reference. A large gap between the two columns is the whole reason these products cannot be compared at the same absolute RMSE.

- `clough_tocher` against the observed field: 1.9140e-03
- `psf` against the latent field: 1.5902e-02
- blur alone (observed minus latent): 1.5717e-02; `psf` recovers -1.2% of it

These are the two numbers a compression target has to be read against, one per product. They are written to results/resampling_error.csv for `alpha_budget.py`.