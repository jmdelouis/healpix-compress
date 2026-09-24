# Bicubic transfer diagnosis

Hypothesis under test: near-zero signed weight sums, floored and then divided, on a small fraction of cells. Facts are in the tables; the hypothesis is confirmed only if the worst cells are the floored ones and the separable bicubic shows no such tail.

## region__agriculture__beauce

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/srv/conda/envs/notebook/lib/python3.13/site-packages/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 70294697415971.766 | 83.00 | 56.36% | 1.621% | 0.000847 | 0.3 | 0.3 | 0.3 |
| affine | 0.087 | 1.18 | 0.00% | 0.001% | -0.628 | 0.959 | -0.496 | 0.896 |
| wave_60m | 11.426 | 156.29 | 37.86% | 3.524% | -25 | 25.2 | 0.15 | 0.25 |
| analytic | 6.223 | 73.87 | 2.04% | 1.642% | -6.86e+03 | 7.53e+03 | -52.7 | 51 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17214 of 1048576 (1.642%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 186.1, all cells median 187.5
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.223 | 73.87 | 2.04% |
| default | 16 | 6.223 | 73.87 | 2.04% |
| default | 36 | 6.223 | 73.87 | 2.04% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.057, p99.9/std 0.21, overshoot 0.00%, >5 std 0.000%

## region__forest__amazon

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/srv/conda/envs/notebook/lib/python3.13/site-packages/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 73332161904367.672 | 109.00 | 56.15% | 1.609% | 0.00134 | 0.3 | 0.3 | 0.3 |
| affine | 0.062 | 0.84 | 0.00% | 0.000% | -0.703 | 1.05 | -0.545 | 0.945 |
| wave_60m | 11.549 | 145.79 | 40.50% | 3.234% | -23.3 | 23.6 | 0.153 | 0.247 |
| analytic | 6.285 | 77.76 | 1.99% | 1.687% | -7.76e+03 | 6.99e+03 | -52.6 | 53.2 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17686 of 1048576 (1.687%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 223.0, all cells median 223.4
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.285 | 77.76 | 1.99% |
| default | 16 | 6.285 | 77.76 | 1.99% |
| default | 36 | 6.285 | 77.76 | 1.99% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__urban__cape_town

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/srv/conda/envs/notebook/lib/python3.13/site-packages/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 74960747002066.953 | 115.00 | 56.32% | 1.623% | 0.00188 | 0.3 | 0.3 | 0.3 |
| affine | 0.071 | 0.94 | 0.01% | 0.000% | -0.7 | 1.03 | -0.504 | 0.904 |
| wave_60m | 12.156 | 149.05 | 39.86% | 3.279% | -23.7 | 23.8 | 0.153 | 0.247 |
| analytic | 6.391 | 78.77 | 2.05% | 1.701% | -6.25e+03 | 7.52e+03 | -50.7 | 52.8 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17841 of 1048576 (1.701%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 223.2, all cells median 222.3
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.391 | 78.77 | 2.05% |
| default | 16 | 6.391 | 78.77 | 2.05% |
| default | 36 | 6.391 | 78.77 | 2.05% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__water__caspian_sea

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/srv/conda/envs/notebook/lib/python3.13/site-packages/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 75007740790705.141 | 115.00 | 56.28% | 1.622% | 0.000407 | 0.3 | 0.3 | 0.3 |
| affine | 0.073 | 0.97 | 0.01% | 0.001% | -0.677 | 1.06 | -0.491 | 0.891 |
| wave_60m | 11.927 | 150.26 | 39.06% | 3.318% | -24.1 | 24.3 | 0.152 | 0.248 |
| analytic | 6.468 | 78.05 | 2.00% | 1.705% | -6.24e+03 | 6.69e+03 | -52.3 | 52.8 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17878 of 1048576 (1.705%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 215.9, all cells median 215.2
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.468 | 78.05 | 2.00% |
| default | 16 | 6.468 | 78.05 | 2.00% |
| default | 36 | 6.468 | 78.05 | 2.00% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%
