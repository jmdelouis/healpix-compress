# Reconstruction error per transfer, against its own estimand

- regions: 39, noise: median 2.3061e-03 (6.6% of the field's std)

## By transfer and scene class

| kernel | estimand | class | n | e_resample | % of estimand std | blur error |
|---|---|---|---|---|---|---|
| `clough_tocher` | observed | agriculture | 10 | 2.7954e-03 | 7.5% | 8.8204e-03 |
| `clough_tocher` | observed | forest | 10 | 1.9111e-03 | 7.1% | 2.8265e-03 |
| `clough_tocher` | observed | urban | 10 | 8.1397e-03 | 14.8% | 4.0441e-02 |
| `clough_tocher` | observed | water | 9 | 1.4335e-03 | 10.7% | 3.4545e-03 |
| `psf` | latent | agriculture | 10 | 9.2047e-03 | 21.2% | 8.8204e-03 |
| `psf` | latent | forest | 10 | 4.3047e-03 | 15.7% | 2.8265e-03 |
| `psf` | latent | urban | 10 | 4.0501e-02 | 52.0% | 4.0441e-02 |
| `psf` | latent | water | 9 | 4.0820e-03 | 40.0% | 3.4545e-03 |

## Medians over all regions

| kernel | estimand | e_resample | % of estimand std | against the other estimand |
|---|---|---|---|---|
| `clough_tocher` | observed | 2.9534e-03 | 8.7% | 1.0716e-02 |
| `psf` | latent | 9.9157e-03 | 26.6% | 8.7980e-03 |

The last column is what each transfer would score if it were held to the other product's reference. A large gap between the two columns is the whole reason these products cannot be compared at the same absolute RMSE.

- `clough_tocher` against the observed field: 2.9534e-03
- `psf` against the latent field: 9.9157e-03
- blur alone (observed minus latent): 9.3521e-03; `psf` recovers -6.0% of it

These are the two numbers a compression target has to be read against, one per product. They are written to results/resampling_error.csv for `alpha_budget.py`.