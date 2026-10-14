# The verdict in units of the product's own error

- transfer: `clough_tocher`; reconstruction error measured here on 40 regions in 4 class(es), clough_tocher against its own estimand, the observed field
- operational overhead applied to the ratios: $\gamma$ = 1.3

| class | e_resample |
|---|---|
| agriculture | 2.7954e-03 |
| forest | 1.9111e-03 |
| urban | 8.1397e-03 |
| water | 1.1380e-03 |

## What each target costs, per class

alpha is the compression error as a fraction of the reconstruction error; inflation is how much the product's total error grows by compressing at that target, the two errors adding in quadrature.

| class | target | alpha | total-error inflation | HPX/UTM ratio | n |
|---|---|---|---|---|---|
| agriculture | 0.0001 | 0.04 | +0.1% | 1.63 | 10 |
| agriculture | 0.0003 | 0.11 | +0.6% | 1.51 | 10 |
| agriculture | 0.001 | 0.36 | +6.2% | 1.21 | 10 |
| agriculture | 0.003 | 1.07 | +46.7% | 0.88 | 10 |
| agriculture | 0.01 | 3.58 | +271.4% | 0.69 | 10 |
| forest | 0.0001 | 0.05 | +0.1% | 1.59 | 10 |
| forest | 0.0003 | 0.16 | +1.2% | 1.46 | 10 |
| forest | 0.001 | 0.52 | +12.9% | 1.07 | 10 |
| forest | 0.003 | 1.57 | +86.1% | 0.84 | 10 |
| forest | 0.01 | 5.23 | +432.7% | 0.85 | 8 |
| urban | 0.0001 | 0.01 | +0.0% | 1.67 | 10 |
| urban | 0.0003 | 0.04 | +0.1% | 1.63 | 10 |
| urban | 0.001 | 0.12 | +0.8% | 1.48 | 10 |
| urban | 0.003 | 0.37 | +6.6% | 1.21 | 10 |
| urban | 0.01 | 1.23 | +58.4% | 0.70 | 10 |
| water | 0.0001 | 0.09 | +0.4% | 1.53 | 10 |
| water | 0.0003 | 0.26 | +3.4% | 1.12 | 10 |
| water | 0.001 | 0.88 | +33.1% | 0.98 | 9 |
| water | 0.003 | 2.64 | +181.9% | 0.97 | 7 |
| water | 0.01 | 8.79 | +784.4% | 0.75 | 4 |

## Where HEALPix becomes the cheaper representation

| class | first target with ratio < 1 | alpha there | inflation | ratio |
|---|---|---|---|---|
| agriculture | 0.003 | 1.07 | +46.7% | 0.88 |
| forest | 0.003 | 1.57 | +86.1% | 0.84 |
| urban | 0.01 | 1.23 | +58.4% | 0.70 |
| water | 0.001 | 0.88 | +33.1% | 0.98 |

- **Over water scenes, HEALPix becomes 2% cheaper to store at a compression error that inflates the product's total error by 33.1%** (target 0.001, alpha 0.88).

- **7 of 20 (class, target) pairs have alpha > 1**, i.e. the compression error there exceeds the reconstruction error itself: agriculture at 0.003, agriculture at 0.01, forest at 0.003, forest at 0.01, urban at 0.01, water at 0.003, water at 0.01. Under this denominator those targets are not conservative, and the inflation figures beside them say by how much.