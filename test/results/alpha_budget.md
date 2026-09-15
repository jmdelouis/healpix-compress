# The verdict in units of the product's own error

- transfer: `clough_tocher`; reconstruction error measured here on 40 regions in 4 class(es), clough_tocher against its own estimand, the observed field
- operational overhead applied to the ratios: $\gamma$ = 1.3

| class | e_resample |
|---|---|
| agriculture | 1.9444e-03 |
| forest | 1.3159e-03 |
| urban | 6.2967e-03 |
| water | 6.9155e-04 |

## What each target costs, per class

alpha is the compression error as a fraction of the reconstruction error; inflation is how much the product's total error grows by compressing at that target, the two errors adding in quadrature.

| class | target | alpha | total-error inflation | HPX/UTM ratio | n |
|---|---|---|---|---|---|
| agriculture | 0.0001 | 0.05 | +0.1% | 1.76 | 10 |
| agriculture | 0.0003 | 0.15 | +1.2% | 1.65 | 10 |
| agriculture | 0.001 | 0.51 | +12.5% | 1.36 | 10 |
| agriculture | 0.003 | 1.54 | +83.9% | 1.11 | 10 |
| agriculture | 0.01 | 5.14 | +423.9% | 1.20 | 10 |
| forest | 0.0001 | 0.08 | +0.3% | 1.69 | 10 |
| forest | 0.0003 | 0.23 | +2.6% | 1.71 | 10 |
| forest | 0.001 | 0.76 | +25.6% | 1.22 | 10 |
| forest | 0.003 | 2.28 | +148.9% | 1.05 | 10 |
| forest | 0.01 | 7.60 | +666.5% | 1.23 | 9 |
| urban | 0.0001 | 0.02 | +0.0% | 1.79 | 10 |
| urban | 0.0003 | 0.05 | +0.1% | 1.71 | 10 |
| urban | 0.001 | 0.16 | +1.3% | 1.56 | 10 |
| urban | 0.003 | 0.48 | +10.8% | 1.38 | 10 |
| urban | 0.01 | 1.59 | +87.7% | 0.98 | 10 |
| water | 0.0001 | 0.14 | +1.0% | 1.68 | 10 |
| water | 0.0003 | 0.43 | +9.0% | 1.51 | 10 |
| water | 0.001 | 1.45 | +75.8% | 1.21 | 9 |
| water | 0.003 | 4.34 | +345.2% | 1.46 | 7 |
| water | 0.01 | 14.46 | +1349.5% | 1.02 | 4 |

## Where HEALPix becomes the cheaper representation

| class | first target with ratio < 1 | alpha there | inflation | ratio |
|---|---|---|---|---|
| agriculture | none tested | (at 0.01: 5.14) | +423.9% | 1.20 |
| forest | none tested | (at 0.01: 7.60) | +666.5% | 1.23 |
| urban | 0.01 | 1.59 | +87.7% | 0.98 |
| water | none tested | (at 0.01: 14.46) | +1349.5% | 1.02 |

- **Over urban scenes, HEALPix becomes 2% cheaper to store at a compression error that inflates the product's total error by 87.7%** (target 0.01, alpha 1.59).

- **8 of 20 (class, target) pairs have alpha > 1**, i.e. the compression error there exceeds the reconstruction error itself: agriculture at 0.003, agriculture at 0.01, forest at 0.003, forest at 0.01, urban at 0.01, water at 0.001, water at 0.003, water at 0.01. Under this denominator those targets are not conservative, and the inflation figures beside them say by how much.