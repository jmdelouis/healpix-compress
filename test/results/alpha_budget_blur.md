# The verdict in units of the product's own error

- transfer: `clough_tocher`; reconstruction error measured here on 40 regions in 4 class(es), clough_tocher against the latent scene (i.e. the blur)
- operational overhead applied to the ratios: $\gamma$ = 1.3

| class | e_resample |
|---|---|
| agriculture | 1.5946e-02 |
| forest | 5.9653e-03 |
| urban | 7.5302e-02 |
| water | 4.4331e-03 |

## What each target costs, per class

alpha is the compression error as a fraction of the reconstruction error; inflation is how much the product's total error grows by compressing at that target, the two errors adding in quadrature.

| class | target | alpha | total-error inflation | HPX/UTM ratio | n |
|---|---|---|---|---|---|
| agriculture | 0.0001 | 0.01 | +0.0% | 1.63 | 10 |
| agriculture | 0.0003 | 0.02 | +0.0% | 1.51 | 10 |
| agriculture | 0.001 | 0.06 | +0.2% | 1.21 | 10 |
| agriculture | 0.003 | 0.19 | +1.8% | 0.88 | 10 |
| agriculture | 0.01 | 0.63 | +18.0% | 0.81 | 10 |
| forest | 0.0001 | 0.02 | +0.0% | 1.57 | 10 |
| forest | 0.0003 | 0.05 | +0.1% | 1.42 | 10 |
| forest | 0.001 | 0.17 | +1.4% | 1.00 | 10 |
| forest | 0.003 | 0.50 | +11.9% | 0.83 | 10 |
| forest | 0.01 | 1.68 | +95.2% | 0.96 | 10 |
| urban | 0.0001 | 0.00 | +0.0% | 1.67 | 10 |
| urban | 0.0003 | 0.00 | +0.0% | 1.63 | 10 |
| urban | 0.001 | 0.01 | +0.0% | 1.48 | 10 |
| urban | 0.003 | 0.04 | +0.1% | 1.21 | 10 |
| urban | 0.01 | 0.13 | +0.9% | 0.70 | 10 |
| water | 0.0001 | 0.02 | +0.0% | 1.55 | 10 |
| water | 0.0003 | 0.07 | +0.2% | 1.12 | 10 |
| water | 0.001 | 0.23 | +2.5% | 1.16 | 10 |
| water | 0.003 | 0.68 | +20.7% | 1.14 | 10 |
| water | 0.01 | 2.26 | +146.7% | 0.94 | 10 |

## Where HEALPix becomes the cheaper representation

| class | first target with ratio < 1 | alpha there | inflation | ratio |
|---|---|---|---|---|
| agriculture | 0.003 | 0.19 | +1.8% | 0.88 |
| forest | 0.003 | 0.50 | +11.9% | 0.83 |
| urban | 0.01 | 0.13 | +0.9% | 0.70 |
| water | 0.01 | 2.26 | +146.7% | 0.94 |

- **Over urban scenes, HEALPix becomes 30% cheaper to store at a compression error that inflates the product's total error by 0.9%** (target 0.01, alpha 0.13).

- **2 of 20 (class, target) pairs have alpha > 1**, i.e. the compression error there exceeds the reconstruction error itself: forest at 0.01, water at 0.01. Under this denominator those targets are not conservative, and the inflation figures beside them say by how much.