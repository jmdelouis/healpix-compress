# The verdict in units of the product's own error

- transfer: `clough_tocher`; reconstruction error healpix-resample, controlled synthetic (Esri texture)
- operational overhead applied to the ratios: $\gamma$ = 1.3

| class | e_resample |
|---|---|
| agriculture | 2.2300e-02 |
| forest | 2.5900e-02 |
| urban | 5.1100e-02 |
| water | 1.1400e-02 |

## What each target costs, per class

alpha is the compression error as a fraction of the reconstruction error; inflation is how much the product's total error grows by compressing at that target, the two errors adding in quadrature.

| class | target | alpha | total-error inflation | HPX/UTM ratio | n |
|---|---|---|---|---|---|
| agriculture | 0.0001 | 0.00 | +0.0% | 1.63 | 10 |
| agriculture | 0.0003 | 0.01 | +0.0% | 1.51 | 10 |
| agriculture | 0.001 | 0.04 | +0.1% | 1.22 | 10 |
| agriculture | 0.003 | 0.13 | +0.9% | 0.89 | 10 |
| agriculture | 0.01 | 0.45 | +9.6% | 0.85 | 10 |
| forest | 0.0001 | 0.00 | +0.0% | 1.57 | 10 |
| forest | 0.0003 | 0.01 | +0.0% | 1.42 | 10 |
| forest | 0.001 | 0.04 | +0.1% | 1.01 | 10 |
| forest | 0.003 | 0.12 | +0.7% | 0.83 | 10 |
| forest | 0.01 | 0.39 | +7.2% | 0.96 | 10 |
| urban | 0.0001 | 0.00 | +0.0% | 1.67 | 10 |
| urban | 0.0003 | 0.01 | +0.0% | 1.63 | 10 |
| urban | 0.001 | 0.02 | +0.0% | 1.48 | 10 |
| urban | 0.003 | 0.06 | +0.2% | 1.21 | 10 |
| urban | 0.01 | 0.20 | +1.9% | 0.70 | 10 |
| water | 0.0001 | 0.01 | +0.0% | 1.55 | 10 |
| water | 0.0003 | 0.03 | +0.0% | 1.12 | 10 |
| water | 0.001 | 0.09 | +0.4% | 1.16 | 10 |
| water | 0.003 | 0.26 | +3.4% | 1.18 | 10 |
| water | 0.01 | 0.88 | +33.0% | 0.98 | 10 |

## Where HEALPix becomes the cheaper representation

| class | first target with ratio < 1 | alpha there | inflation | ratio |
|---|---|---|---|---|
| agriculture | 0.003 | 0.13 | +0.9% | 0.89 |
| forest | 0.003 | 0.12 | +0.7% | 0.83 |
| urban | 0.01 | 0.20 | +1.9% | 0.70 |
| water | 0.01 | 0.88 | +33.0% | 0.98 |

- **Over forest scenes, HEALPix becomes 17% cheaper to store at a compression error that inflates the product's total error by 0.7%** (target 0.003, alpha 0.12).

Every reported target sits below the reconstruction error of its own class, so no result in this paper is obtained by compressing more finely than the transfer step is accurate.