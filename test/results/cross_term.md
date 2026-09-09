# Is the compression error orthogonal to the transfer error?

- regions: 40, kernels: clough_tocher, codecs: jpeg2000_lossy, zfp
- operating points: 1440

- transfer error reproduces the published measurement to 3.4e-11% at worst over 720 (region, kernel) pairs

## rho by transfer and reference

| kernel | reference | n | median rho | 10th | 90th | share |rho| > 0.1 |
|---|---|---|---|---|---|---|
| `clough_tocher` | latent | 720 | +0.001 | -0.023 | +0.237 | 23% |
| `clough_tocher` | observed | 720 | -0.004 | -0.369 | +0.006 | 38% |

## what the quadrature assumption costs, near the reported kappa values

| kernel | reference | kappa | n | measured | quadrature | difference |
|---|---|---|---|---|---|---|
| `clough_tocher` | latent | ~0.1 | 99 | +1.2% | +0.5% | +0.6 pt |
| `clough_tocher` | latent | ~0.2 | 97 | +4.5% | +2.0% | +2.5 pt |
| `clough_tocher` | latent | ~0.5 | 59 | +16.7% | +10.1% | +6.6 pt |
| `clough_tocher` | latent | ~1 | 13 | +23.9% | +29.4% | -5.5 pt |
| `clough_tocher` | observed | ~0.1 | 30 | +0.4% | +0.5% | -0.1 pt |
| `clough_tocher` | observed | ~0.2 | 66 | +1.4% | +1.9% | -0.5 pt |
| `clough_tocher` | observed | ~0.5 | 95 | +3.1% | +14.2% | -11.2 pt |
| `clough_tocher` | observed | ~1 | 119 | +16.7% | +40.4% | -23.7 pt |

- `clough_tocher` against the latent field: median $\rho$ = +0.001; quadrature **understates** the inflation, by at most 6.6 percentage points (at $\kappa\approx0.5$). Quadrature therefore understates the damage and Section V-H cannot use it as written.
- `clough_tocher` against the observed field: median $\rho$ = -0.004; quadrature **overstates** the inflation, by at most 23.7 percentage points (at $\kappa\approx1$). Quadrature is therefore a conservative bound on the damage, and Section V-H can say so.