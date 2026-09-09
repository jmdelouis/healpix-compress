# Is the compression error orthogonal to the transfer error?

- regions: 40, kernels: clough_tocher, codecs: jpeg2000_lossy, zfp
- operating points: 720

- transfer error reproduces the published measurement to 3.4e-11% at worst over 720 (region, kernel) pairs

## rho by kernel, over all operating points

| kernel | n | median rho | 10th | 90th | share |rho| > 0.1 |
|---|---|---|---|---|---|
| `clough_tocher` | 720 | -0.004 | -0.369 | +0.006 | 38% |

## what the quadrature assumption costs, near the reported kappa values

| kernel | kappa | n | measured inflation | quadrature | difference |
|---|---|---|---|---|---|
| `clough_tocher` | ~0.1 | 30 | +0.4% | +0.5% | -0.1 pt |
| `clough_tocher` | ~0.2 | 66 | +1.4% | +1.9% | -0.5 pt |
| `clough_tocher` | ~0.5 | 95 | +3.1% | +14.2% | -11.2 pt |
| `clough_tocher` | ~1 | 119 | +16.7% | +40.4% | -23.7 pt |

- **The cross term is not negligible**: median $\rho$ = -0.004, and only 62% of operating points have $|\rho| \le 0.1$. Section V-H must either carry $\rho$ explicitly or be declared a scenario under an orthogonality assumption.