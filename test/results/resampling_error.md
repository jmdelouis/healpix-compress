# Reconstruction error per transfer, against its own estimand

- regions: 40, noise: median 2.2759e-03 (7.2% of the field's std)

## By transfer and scene class

| kernel | estimand | class | n | e_resample | % of estimand std | blur error |
|---|---|---|---|---|---|---|
| `clough_tocher` | observed | agriculture | 10 | 2.7954e-03 | 7.5% | 8.8204e-03 |
| `clough_tocher` | observed | forest | 10 | 1.9111e-03 | 8.0% | 2.1256e-03 |
| `clough_tocher` | observed | urban | 10 | 8.1397e-03 | 14.8% | 4.0441e-02 |
| `clough_tocher` | observed | water | 10 | 1.1380e-03 | 9.4% | 2.4350e-03 |
| `psf` | latent | agriculture | 10 | 9.2047e-03 | 21.2% | 8.8204e-03 |
| `psf` | latent | forest | 10 | 4.3720e-03 | 17.2% | 2.1256e-03 |
| `psf` | latent | urban | 10 | 4.0501e-02 | 52.0% | 4.0441e-02 |
| `psf` | latent | water | 10 | 3.2386e-03 | 36.4% | 2.4350e-03 |

## Medians over all regions

| kernel | estimand | e_resample | % of estimand std | against the other estimand |
|---|---|---|---|---|
| `clough_tocher` | observed | 2.7954e-03 | 8.4% | 9.9238e-03 |
| `psf` | latent | 8.6732e-03 | 26.0% | 7.8630e-03 |

The last column is what each transfer would score if it were held to the other product's reference. A large gap between the two columns is the whole reason these products cannot be compared at the same absolute RMSE.

- `clough_tocher` against the observed field: 2.7954e-03
- `psf` against the latent field: 8.6732e-03
- blur alone (observed minus latent): 8.8204e-03; `psf` recovers 1.7% of it

These are the two numbers a compression target has to be read against, one per product. They are written to results/resampling_error.csv for `alpha_budget.py`.