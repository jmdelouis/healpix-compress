# The verdict in units of the product's own error

- transfer: `clough_tocher`; reconstruction error measured here on 39 regions in 4 class(es), clough_tocher against its own estimand, the observed field
- operational overhead applied to the ratios: $\gamma$ = 1.3

| class | e_resample |
|---|---|
| agriculture | 1.9444e-03 |
| forest | 1.3159e-03 |
| urban | 6.2967e-03 |
| water | 8.2038e-04 |

## What each target costs, per class

alpha is the compression error as a fraction of the reconstruction error; inflation is how much the product's total error grows by compressing at that target, the two errors adding in quadrature.

| class | target | alpha | total-error inflation | HPX/UTM ratio | n |
|---|---|---|---|---|---|
| agriculture | 0.0001 | 0.05 | +0.1% | 1.63 | 10 |
| agriculture | 0.0003 | 0.15 | +1.2% | 1.51 | 10 |
| agriculture | 0.001 | 0.51 | +12.5% | 1.21 | 10 |
| agriculture | 0.003 | 1.54 | +83.9% | 0.88 | 10 |
| agriculture | 0.01 | 5.14 | +423.9% | 0.69 | 10 |
| forest | 0.0001 | 0.08 | +0.3% | 1.57 | 10 |
| forest | 0.0003 | 0.23 | +2.6% | 1.42 | 10 |
| forest | 0.001 | 0.76 | +25.6% | 1.00 | 10 |
| forest | 0.003 | 2.28 | +148.9% | 0.83 | 10 |
| forest | 0.01 | 7.60 | +666.5% | 0.85 | 8 |
| urban | 0.0001 | 0.02 | +0.0% | 1.67 | 10 |
| urban | 0.0003 | 0.05 | +0.1% | 1.63 | 10 |
| urban | 0.001 | 0.16 | +1.3% | 1.48 | 10 |
| urban | 0.003 | 0.48 | +10.8% | 1.21 | 10 |
| urban | 0.01 | 1.59 | +87.7% | 0.70 | 10 |
| water | 0.0001 | 0.12 | +0.7% | 1.55 | 10 |
| water | 0.0003 | 0.37 | +6.5% | 1.12 | 10 |
| water | 0.001 | 1.22 | +57.7% | 0.98 | 9 |
| water | 0.003 | 3.66 | +279.1% | 1.03 | 7 |
| water | 0.01 | 12.19 | +1123.0% | 0.75 | 4 |

## Where HEALPix becomes the cheaper representation

| class | first target with ratio < 1 | alpha there | inflation | ratio |
|---|---|---|---|---|
| agriculture | 0.003 | 1.54 | +83.9% | 0.88 |
| forest | 0.003 | 2.28 | +148.9% | 0.83 |
| urban | 0.01 | 1.59 | +87.7% | 0.70 |
| water | 0.001 | 1.22 | +57.7% | 0.98 |

- **Over water scenes, HEALPix becomes 2% cheaper to store at a compression error that inflates the product's total error by 57.7%** (target 0.001, alpha 1.22).

- **8 of 20 (class, target) pairs have alpha > 1**, i.e. the compression error there exceeds the reconstruction error itself: agriculture at 0.003, agriculture at 0.01, forest at 0.003, forest at 0.01, urban at 0.01, water at 0.001, water at 0.003, water at 0.01. Under this denominator those targets are not conservative, and the inflation figures beside them say by how much.