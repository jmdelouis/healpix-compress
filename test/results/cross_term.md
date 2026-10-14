# Is the compression error orthogonal to the transfer error?

- regions: 40, kernels: clough_tocher, codecs: jpeg2000_lossy, zfp
- operating points: 2240

- transfer error reproduces the published measurement to 2.2e-11% at worst over 1120 (region, kernel) pairs

## rho by transfer and reference

| kernel | reference | n | median rho | 10th | 90th | share |rho| > 0.1 |
|---|---|---|---|---|---|---|
| `clough_tocher` | latent | 1120 | +0.047 | -0.011 | +0.316 | 48% |
| `clough_tocher` | observed | 1120 | -0.041 | -0.204 | +0.001 | 28% |

## what the quadrature assumption costs, near the reported kappa values

| kernel | reference | kappa | n | measured | quadrature | difference |
|---|---|---|---|---|---|---|
| `clough_tocher` | latent | ~0.1 | 66 | +0.6% | +0.5% | +0.0 pt |
| `clough_tocher` | latent | ~0.2 | 93 | +3.4% | +2.3% | +1.1 pt |
| `clough_tocher` | latent | ~0.5 | 148 | +20.4% | +12.9% | +7.5 pt |
| `clough_tocher` | latent | ~1 | 248 | +57.5% | +43.5% | +14.0 pt |
| `clough_tocher` | observed | ~0.1 | 37 | +0.4% | +0.4% | +0.0 pt |
| `clough_tocher` | observed | ~0.2 | 59 | +1.9% | +2.4% | -0.5 pt |
| `clough_tocher` | observed | ~0.5 | 77 | +10.5% | +16.2% | -5.6 pt |
| `clough_tocher` | observed | ~1 | 118 | +26.3% | +40.7% | -14.5 pt |

- `clough_tocher` against the latent field: median $\rho$ = +0.047; quadrature **understates** the inflation, by at most 14.0 percentage points (at $\kappa\approx1$). Quadrature therefore understates the damage and Section V-H cannot use it as written.
- `clough_tocher` against the observed field: median $\rho$ = -0.041; quadrature **overstates** the inflation, by at most 14.5 percentage points (at $\kappa\approx1$). Quadrature is therefore a conservative bound on the damage, and Section V-H can say so.