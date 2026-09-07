# The verdict in units of the product's own error

- transfer: `clough_tocher`; reconstruction error measured here on 3 regions, clough_tocher against the observed field
- operational overhead applied to the ratios: $\gamma$ = 1.3

| class | e_resample |
|---|---|
| agriculture | 1.7745e-03 |

## What each target costs, per class

alpha is the compression error as a fraction of the reconstruction error; inflation is how much the product's total error grows by compressing at that target, the two errors adding in quadrature.

| class | target | alpha | total-error inflation | HPX/UTM ratio | n |
|---|---|---|---|---|---|
| agriculture | 0.0001 | 0.06 | +0.2% | 1.63 | 10 |
| agriculture | 0.0003 | 0.17 | +1.4% | 1.51 | 10 |
| agriculture | 0.001 | 0.56 | +14.8% | 1.21 | 10 |
| agriculture | 0.003 | 1.69 | +96.4% | 0.88 | 10 |
| agriculture | 0.01 | 5.64 | +472.3% | 0.81 | 10 |

## Where HEALPix becomes the cheaper representation

| class | first target with ratio < 1 | alpha there | inflation | ratio |
|---|---|---|---|---|
| agriculture | 0.003 | 1.69 | +96.4% | 0.88 |

- **Over agriculture scenes, HEALPix becomes 12% cheaper to store at a compression error that inflates the product's total error by 96.4%** (target 0.003, alpha 1.69).

Every reported target sits below the reconstruction error of its own class, so no result in this paper is obtained by compressing more finely than the transfer step is accurate.