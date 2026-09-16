# Reconstruction error per transfer, against its own estimand

- regions: 40, noise: median 0.0000e+00 (0.0% of the field's std)

## By transfer and scene class

| kernel | estimand | class | n | e_resample | % of estimand std | blur error |
|---|---|---|---|---|---|---|
| `clough_tocher` | observed | agriculture | 10 | 3.8675e-04 | 1.1% | 1.6090e-02 |
| `clough_tocher` | observed | forest | 10 | 2.9658e-04 | 0.8% | 5.8559e-03 |
| `clough_tocher` | observed | urban | 10 | 9.1693e-04 | 1.6% | 7.7056e-02 |
| `clough_tocher` | observed | water | 10 | 1.1539e-04 | 1.3% | 4.3831e-03 |
| `psf` | latent | agriculture | 10 | 1.5138e-02 | 32.1% | 1.6090e-02 |
| `psf` | latent | forest | 10 | 5.5099e-03 | 21.2% | 5.8559e-03 |
| `psf` | latent | urban | 10 | 7.4896e-02 | 72.7% | 7.7056e-02 |
| `psf` | latent | water | 10 | 4.1286e-03 | 39.5% | 4.3831e-03 |

## Medians over all regions

| kernel | estimand | e_resample | % of estimand std | against the other estimand |
|---|---|---|---|---|
| `clough_tocher` | observed | 3.8994e-04 | 1.2% | 1.6189e-02 |
| `psf` | latent | 1.5138e-02 | 40.0% | 2.0128e-03 |

The last column is what each transfer would score if it were held to the other product's reference. A large gap between the two columns is the whole reason these products cannot be compared at the same absolute RMSE.

- `clough_tocher` against the observed field: 3.8994e-04
- `psf` against the latent field: 1.5138e-02
- blur alone (observed minus latent): 1.6090e-02; `psf` recovers 5.9% of it

These are the two numbers a compression target has to be read against, one per product. They are written to results/resampling_error.csv for `alpha_budget.py`.