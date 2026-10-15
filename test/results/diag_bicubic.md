# Bicubic transfer diagnosis

Hypothesis under test: near-zero signed weight sums, floored and then divided, on a small fraction of cells. Facts are in the tables; the hypothesis is confirmed only if the worst cells are the floored ones and the separable bicubic shows no such tail.

## region__agriculture__beauce

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 70286443412830.781 | 93.00 | 55.32% | 1.633% | 0.000847 | 0.3 | 0.3 | 0.3 |
| affine | 0.087 | 1.18 | 0.00% | 0.001% | -0.628 | 0.959 | -0.496 | 0.896 |
| wave_60m | 11.429 | 156.29 | 37.86% | 3.524% | -25 | 25.2 | 0.15 | 0.25 |
| analytic | 6.226 | 73.80 | 2.04% | 1.642% | -6.86e+03 | 7.53e+03 | -52.7 | 51 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17214 of 1048576 (1.642%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 186.1, all cells median 187.5
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.226 | 73.80 | 2.04% |
| default | 16 | 6.226 | 73.80 | 2.04% |
| default | 36 | 6.226 | 73.80 | 2.04% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.057, p99.9/std 0.21, overshoot 0.00%, >5 std 0.000%

## region__agriculture__canterbury

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 71569014890663.328 | 109.00 | 55.42% | 1.608% | 0.00192 | 0.3 | 0.3 | 0.3 |
| affine | 0.088 | 1.17 | 0.00% | 0.001% | -0.672 | 0.931 | -0.493 | 0.893 |
| wave_60m | 11.584 | 140.49 | 45.18% | 3.126% | -22.6 | 22.8 | 0.155 | 0.245 |
| analytic | 6.185 | 74.43 | 2.02% | 1.633% | -8.41e+03 | 7.13e+03 | -51.8 | 51.7 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17124 of 1048576 (1.633%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 191.3, all cells median 192.9
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.185 | 74.43 | 2.02% |
| default | 16 | 6.185 | 74.43 | 2.02% |
| default | 36 | 6.185 | 74.43 | 2.02% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.057, p99.9/std 0.21, overshoot 0.00%, >5 std 0.000%

## region__agriculture__iowa

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 74064456448809.234 | 115.00 | 55.32% | 1.634% | 0.000281 | 0.3 | 0.3 | 0.3 |
| affine | 0.048 | 0.62 | 0.01% | 0.000% | -0.835 | 1.13 | -0.489 | 0.89 |
| wave_60m | 11.815 | 153.51 | 38.80% | 3.352% | -24.2 | 24.4 | 0.152 | 0.248 |
| analytic | 6.521 | 82.24 | 2.00% | 1.732% | -7.03e+03 | 6.26e+03 | -52.3 | 52.9 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 18166 of 1048576 (1.732%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 184.3, all cells median 184.2
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.521 | 82.24 | 2.00% |
| default | 16 | 6.521 | 82.24 | 2.00% |
| default | 36 | 6.521 | 82.24 | 2.00% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.059, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__agriculture__mato_grosso

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 72346999361300.656 | 115.00 | 55.12% | 1.636% | 0.000183 | 0.3 | 0.3 | 0.3 |
| affine | 0.064 | 0.86 | 0.00% | 0.000% | -0.762 | 1.17 | -0.538 | 0.938 |
| wave_60m | 11.580 | 141.94 | 42.66% | 3.179% | -22.9 | 23.1 | 0.154 | 0.246 |
| analytic | 6.499 | 77.37 | 2.03% | 1.677% | -7.05e+03 | 6.62e+03 | -51.3 | 52.4 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17589 of 1048576 (1.677%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 224.8, all cells median 223.9
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.499 | 77.37 | 2.03% |
| default | 16 | 6.499 | 77.37 | 2.03% |
| default | 36 | 6.499 | 77.37 | 2.03% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__agriculture__nile_delta

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 75137218245837.516 | 115.00 | 55.22% | 1.623% | 0.000662 | 0.3 | 0.3 | 0.3 |
| affine | 0.069 | 0.92 | 0.01% | 0.000% | -0.69 | 1.02 | -0.51 | 0.91 |
| wave_60m | 11.683 | 145.53 | 51.49% | 2.997% | -21.8 | 22 | 0.156 | 0.244 |
| analytic | 6.330 | 78.41 | 2.06% | 1.679% | -5.57e+03 | 5.94e+03 | -49.8 | 52.3 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17609 of 1048576 (1.679%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 223.9, all cells median 223.6
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.330 | 78.41 | 2.06% |
| default | 16 | 6.330 | 78.41 | 2.06% |
| default | 36 | 6.330 | 78.41 | 2.06% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__agriculture__pampas

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 75925337951870.516 | 115.00 | 55.29% | 1.631% | 0.000122 | 0.3 | 0.3 | 0.3 |
| affine | 0.072 | 0.99 | 0.01% | 0.000% | -0.704 | 1.06 | -0.503 | 0.903 |
| wave_60m | 12.225 | 151.85 | 39.91% | 3.271% | -23.6 | 23.8 | 0.153 | 0.247 |
| analytic | 6.413 | 80.49 | 2.00% | 1.707% | -5.74e+03 | 6.18e+03 | -51.9 | 54.1 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17901 of 1048576 (1.707%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 222.5, all cells median 221.9
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.413 | 80.49 | 2.00% |
| default | 16 | 6.413 | 80.49 | 2.00% |
| default | 36 | 6.413 | 80.49 | 2.00% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__agriculture__po_valley

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 71274928681372.594 | 115.00 | 55.20% | 1.637% | 4.54e-05 | 0.3 | 0.3 | 0.3 |
| affine | 0.089 | 1.20 | 0.00% | 0.001% | -0.742 | 0.995 | -0.489 | 0.89 |
| wave_60m | 11.597 | 144.04 | 39.67% | 3.287% | -23.8 | 23.9 | 0.152 | 0.248 |
| analytic | 6.366 | 77.06 | 2.05% | 1.659% | -6.7e+03 | 6.12e+03 | -51.2 | 52.5 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17397 of 1048576 (1.659%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 193.6, all cells median 193.3
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.366 | 77.06 | 2.05% |
| default | 16 | 6.366 | 77.06 | 2.05% |
| default | 36 | 6.366 | 77.06 | 2.05% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.057, p99.9/std 0.21, overshoot 0.00%, >5 std 0.000%

## region__agriculture__punjab

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 71386966389422.922 | 115.00 | 55.16% | 1.626% | 0.000325 | 0.3 | 0.3 | 0.3 |
| affine | 0.069 | 0.88 | 0.01% | 0.000% | -0.644 | 1.06 | -0.51 | 0.91 |
| wave_60m | 11.433 | 145.97 | 39.35% | 3.308% | -23.9 | 24.1 | 0.152 | 0.248 |
| analytic | 6.322 | 76.33 | 1.98% | 1.687% | -6.94e+03 | 7.25e+03 | -52.4 | 53.4 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17691 of 1048576 (1.687%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 222.6, all cells median 223.6
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.322 | 76.33 | 1.98% |
| default | 16 | 6.322 | 76.33 | 1.98% |
| default | 36 | 6.322 | 76.33 | 1.98% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__agriculture__ukraine

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 71297020308885.578 | 115.00 | 55.35% | 1.629% | 0.000569 | 0.3 | 0.3 | 0.3 |
| affine | 0.084 | 1.10 | 0.01% | 0.001% | -0.669 | 1.02 | -0.49 | 0.889 |
| wave_60m | 12.049 | 147.99 | 39.18% | 3.346% | -24 | 24.2 | 0.152 | 0.248 |
| analytic | 6.527 | 82.12 | 2.06% | 1.695% | -7e+03 | 6.91e+03 | -52.8 | 50.4 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17778 of 1048576 (1.695%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 211.4, all cells median 210.0
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.527 | 82.12 | 2.06% |
| default | 16 | 6.527 | 82.12 | 2.06% |
| default | 36 | 6.527 | 82.12 | 2.06% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__agriculture__western_australia

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 73586624536996.234 | 115.00 | 55.27% | 1.641% | 0.000506 | 0.3 | 0.3 | 0.3 |
| affine | 0.069 | 0.92 | 0.00% | 0.000% | -0.697 | 1.03 | -0.508 | 0.908 |
| wave_60m | 11.506 | 147.06 | 38.73% | 3.350% | -24.2 | 24.4 | 0.152 | 0.248 |
| analytic | 6.142 | 74.42 | 2.00% | 1.690% | -6.19e+03 | 8.92e+03 | -51.6 | 53.7 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17716 of 1048576 (1.690%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 222.4, all cells median 223.3
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.142 | 74.42 | 2.00% |
| default | 16 | 6.142 | 74.42 | 2.00% |
| default | 36 | 6.142 | 74.42 | 2.00% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__forest__amazon

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 73332161904367.719 | 115.00 | 55.15% | 1.627% | 0.00134 | 0.3 | 0.3 | 0.3 |
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

## region__forest__appalachian

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 71373067412806.953 | 115.00 | 55.31% | 1.648% | 0.00118 | 0.3 | 0.3 | 0.3 |
| affine | 0.073 | 0.98 | 0.01% | 0.001% | -0.693 | 1.08 | -0.501 | 0.901 |
| wave_60m | 11.785 | 153.32 | 38.49% | 3.391% | -24.4 | 24.6 | 0.151 | 0.249 |
| analytic | 6.677 | 82.71 | 2.04% | 1.727% | -6.23e+03 | 8.36e+03 | -52.1 | 53 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 18111 of 1048576 (1.727%); their RMSE alone accounts for 99.1% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 221.5, all cells median 221.3
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.677 | 82.71 | 2.04% |
| default | 16 | 6.677 | 82.71 | 2.04% |
| default | 36 | 6.677 | 82.71 | 2.04% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__forest__black_forest

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 71193540660544.328 | 115.00 | 55.31% | 1.637% | 5.22e-05 | 0.3 | 0.3 | 0.3 |
| affine | 0.088 | 1.16 | 0.00% | 0.002% | -0.651 | 1.04 | -0.494 | 0.894 |
| wave_60m | 11.534 | 142.52 | 41.07% | 3.228% | -23.2 | 23.5 | 0.153 | 0.247 |
| analytic | 6.091 | 72.19 | 2.09% | 1.667% | -7.43e+03 | 6.22e+03 | -51 | 50.9 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17485 of 1048576 (1.667%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 194.7, all cells median 195.0
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.091 | 72.19 | 2.09% |
| default | 16 | 6.091 | 72.19 | 2.09% |
| default | 36 | 6.091 | 72.19 | 2.09% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.057, p99.9/std 0.21, overshoot 0.00%, >5 std 0.000%

## region__forest__borneo

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 74257571171200.078 | 115.00 | 55.25% | 1.641% | 0.000175 | 0.3 | 0.3 | 0.3 |
| affine | 0.062 | 0.82 | 0.00% | 0.000% | -0.705 | 1.02 | -0.545 | 0.945 |
| wave_60m | 11.773 | 144.65 | 42.69% | 3.195% | -22.9 | 23.1 | 0.154 | 0.246 |
| analytic | 6.343 | 76.22 | 2.03% | 1.689% | -1.04e+04 | 7.03e+03 | -51.7 | 52.7 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17706 of 1048576 (1.689%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 223.1, all cells median 223.4
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.343 | 76.22 | 2.03% |
| default | 16 | 6.343 | 76.22 | 2.03% |
| default | 36 | 6.343 | 76.22 | 2.03% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__forest__canadian_rockies

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 71508846257037.578 | 109.00 | 55.24% | 1.635% | 0.00025 | 0.3 | 0.3 | 0.3 |
| affine | 0.062 | 0.80 | 0.01% | 0.000% | -0.664 | 1.11 | -0.495 | 0.895 |
| wave_60m | 11.690 | 161.48 | 37.96% | 3.446% | -24.8 | 25 | 0.15 | 0.25 |
| analytic | 6.228 | 78.47 | 2.05% | 1.718% | -7.97e+03 | 7.58e+03 | -52.8 | 50.8 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 18013 of 1048576 (1.718%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 218.0, all cells median 215.4
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.228 | 78.47 | 2.05% |
| default | 16 | 6.228 | 78.47 | 2.05% |
| default | 36 | 6.228 | 78.47 | 2.05% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__forest__congo

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 73846233904675.094 | 115.00 | 55.20% | 1.639% | 0.000181 | 0.3 | 0.3 | 0.3 |
| affine | 0.063 | 0.81 | 0.00% | 0.000% | -0.713 | 1.1 | -0.544 | 0.945 |
| wave_60m | 11.867 | 156.18 | 38.29% | 3.382% | -24.5 | 24.6 | 0.151 | 0.249 |
| analytic | 6.553 | 78.24 | 2.01% | 1.693% | -6.08e+03 | 6.84e+03 | -52.4 | 52.9 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17751 of 1048576 (1.693%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 223.0, all cells median 223.0
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.553 | 78.24 | 2.01% |
| default | 16 | 6.553 | 78.24 | 2.01% |
| default | 36 | 6.553 | 78.24 | 2.01% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__forest__finland_boreal

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 73279731301557.344 | 115.00 | 55.30% | 1.628% | 0.000337 | 0.3 | 0.3 | 0.3 |
| affine | 0.084 | 1.10 | 0.00% | 0.001% | -0.657 | 0.95 | -0.499 | 0.9 |
| wave_60m | 11.801 | 166.70 | 37.87% | 3.463% | -25 | 25.2 | 0.15 | 0.25 |
| analytic | 6.215 | 74.81 | 2.01% | 1.690% | -6.23e+03 | 7.35e+03 | -52.3 | 52.1 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17721 of 1048576 (1.690%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 218.3, all cells median 214.9
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.215 | 74.81 | 2.01% |
| default | 16 | 6.215 | 74.81 | 2.01% |
| default | 36 | 6.215 | 74.81 | 2.01% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__forest__new_zealand

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 72968256594765.641 | 115.00 | 55.27% | 1.625% | 0.000378 | 0.3 | 0.3 | 0.3 |
| affine | 0.074 | 0.97 | 0.01% | 0.000% | -0.73 | 1.04 | -0.491 | 0.891 |
| wave_60m | 12.030 | 147.96 | 40.13% | 3.288% | -23.5 | 23.7 | 0.153 | 0.247 |
| analytic | 6.252 | 78.40 | 2.08% | 1.712% | -7.59e+03 | 7.11e+03 | -50.8 | 50.7 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17951 of 1048576 (1.712%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 215.9, all cells median 215.6
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.252 | 78.40 | 2.08% |
| default | 16 | 6.252 | 78.40 | 2.08% |
| default | 36 | 6.252 | 78.40 | 2.08% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__forest__patagonia

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 68701561153882.062 | 115.00 | 55.25% | 1.615% | 0.00156 | 0.3 | 0.3 | 0.3 |
| affine | 0.052 | 0.67 | 0.01% | 0.000% | -0.751 | 1.09 | -0.486 | 0.887 |
| wave_60m | 11.390 | 143.09 | 50.35% | 3.038% | -21.9 | 22.1 | 0.156 | 0.244 |
| analytic | 6.443 | 76.77 | 2.04% | 1.715% | -7.03e+03 | 7.86e+03 | -51.9 | 51.6 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17986 of 1048576 (1.715%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 200.0, all cells median 198.5
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.443 | 76.77 | 2.04% |
| default | 16 | 6.443 | 76.77 | 2.04% |
| default | 36 | 6.443 | 76.77 | 2.04% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.059, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__forest__tasmania

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 72590475982728.281 | 115.00 | 55.30% | 1.648% | 0.000464 | 0.3 | 0.3 | 0.3 |
| affine | 0.081 | 1.08 | 0.01% | 0.001% | -0.717 | 1.06 | -0.49 | 0.889 |
| wave_60m | 11.513 | 143.54 | 40.43% | 3.238% | -23.4 | 23.6 | 0.153 | 0.247 |
| analytic | 6.475 | 77.74 | 2.07% | 1.677% | -7.18e+03 | 7.29e+03 | -52.5 | 49.9 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17588 of 1048576 (1.677%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 212.2, all cells median 211.8
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.475 | 77.74 | 2.07% |
| default | 16 | 6.475 | 77.74 | 2.07% |
| default | 36 | 6.475 | 77.74 | 2.07% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__urban__cape_town

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 74976162165670.391 | 115.00 | 55.32% | 1.642% | 0.00188 | 0.3 | 0.3 | 0.3 |
| affine | 0.071 | 0.94 | 0.01% | 0.000% | -0.7 | 1.03 | -0.504 | 0.904 |
| wave_60m | 12.156 | 149.05 | 39.86% | 3.279% | -23.7 | 23.8 | 0.153 | 0.247 |
| analytic | 6.392 | 78.77 | 2.05% | 1.702% | -6.25e+03 | 7.52e+03 | -50.7 | 52.8 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17842 of 1048576 (1.702%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 223.2, all cells median 222.3
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.392 | 78.77 | 2.05% |
| default | 16 | 6.392 | 78.77 | 2.05% |
| default | 36 | 6.392 | 78.77 | 2.05% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__urban__mexico_city

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 73727018029008.156 | 115.00 | 55.19% | 1.643% | 0.000539 | 0.3 | 0.3 | 0.3 |
| affine | 0.066 | 0.90 | 0.00% | 0.000% | -0.764 | 0.986 | -0.529 | 0.929 |
| wave_60m | 11.826 | 158.58 | 38.40% | 3.393% | -24.4 | 24.6 | 0.151 | 0.249 |
| analytic | 6.527 | 77.28 | 2.03% | 1.707% | -8.21e+03 | 8.27e+03 | -51.3 | 53.4 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17896 of 1048576 (1.707%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 223.9, all cells median 224.4
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.527 | 77.28 | 2.03% |
| default | 16 | 6.527 | 77.28 | 2.03% |
| default | 36 | 6.527 | 77.28 | 2.03% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__urban__montreal

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 75967750853793.469 | 115.00 | 55.32% | 1.646% | 0.000163 | 0.3 | 0.3 | 0.3 |
| affine | 0.090 | 1.21 | 0.01% | 0.002% | -0.628 | 1.03 | -0.489 | 0.889 |
| wave_60m | 11.944 | 146.64 | 42.46% | 3.175% | -22.9 | 23.2 | 0.154 | 0.246 |
| analytic | 6.193 | 74.99 | 2.05% | 1.666% | -6.55e+03 | 8.19e+03 | -52.2 | 51.4 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17467 of 1048576 (1.666%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 201.0, all cells median 199.9
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.193 | 74.99 | 2.05% |
| default | 16 | 6.193 | 74.99 | 2.05% |
| default | 36 | 6.193 | 74.99 | 2.05% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.057, p99.9/std 0.21, overshoot 0.00%, >5 std 0.000%

## region__urban__nairobi

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 72204573543821.328 | 109.00 | 55.16% | 1.634% | 0.000432 | 0.3 | 0.3 | 0.3 |
| affine | 0.062 | 0.82 | 0.00% | 0.000% | -0.76 | 1.15 | -0.545 | 0.945 |
| wave_60m | 11.709 | 154.30 | 38.28% | 3.391% | -24.5 | 24.7 | 0.151 | 0.249 |
| analytic | 6.201 | 77.92 | 2.02% | 1.684% | -6.73e+03 | 5.79e+03 | -52 | 52.7 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17660 of 1048576 (1.684%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 223.9, all cells median 223.2
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.201 | 77.92 | 2.02% |
| default | 16 | 6.201 | 77.92 | 2.02% |
| default | 36 | 6.201 | 77.92 | 2.02% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.21, overshoot 0.00%, >5 std 0.000%

## region__urban__paris

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 75421423355250.812 | 115.00 | 55.26% | 1.627% | 0.000295 | 0.3 | 0.3 | 0.3 |
| affine | 0.088 | 1.17 | 0.00% | 0.001% | -0.664 | 1.05 | -0.494 | 0.894 |
| wave_60m | 11.531 | 161.42 | 37.83% | 3.464% | -24.9 | 25.1 | 0.15 | 0.25 |
| analytic | 6.192 | 74.45 | 2.06% | 1.666% | -6.12e+03 | 6e+03 | -52.7 | 50.7 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17468 of 1048576 (1.666%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 187.0, all cells median 187.6
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.192 | 74.45 | 2.06% |
| default | 16 | 6.192 | 74.45 | 2.06% |
| default | 36 | 6.192 | 74.45 | 2.06% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.057, p99.9/std 0.21, overshoot 0.00%, >5 std 0.000%

## region__urban__reykjavik

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 74358145249685.156 | 109.00 | 55.39% | 1.621% | 0.000947 | 0.3 | 0.3 | 0.3 |
| affine | 0.055 | 0.70 | 0.01% | 0.000% | -0.76 | 1.1 | -0.496 | 0.896 |
| wave_60m | 11.497 | 154.30 | 39.00% | 3.352% | -24.1 | 24.3 | 0.152 | 0.248 |
| analytic | 6.456 | 79.75 | 2.05% | 1.720% | -7.02e+03 | 7.92e+03 | -52.5 | 50.8 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 18038 of 1048576 (1.720%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 209.0, all cells median 209.9
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.456 | 79.75 | 2.05% |
| default | 16 | 6.456 | 79.75 | 2.05% |
| default | 36 | 6.456 | 79.75 | 2.05% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.059, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__urban__sao_paulo

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 74964603749530.141 | 115.00 | 55.16% | 1.614% | 0.000398 | 0.3 | 0.3 | 0.3 |
| affine | 0.066 | 0.90 | 0.00% | 0.000% | -0.695 | 1.11 | -0.523 | 0.923 |
| wave_60m | 11.874 | 149.01 | 45.71% | 3.108% | -22.5 | 22.7 | 0.155 | 0.245 |
| analytic | 6.332 | 78.86 | 2.03% | 1.670% | -6.41e+03 | 6.87e+03 | -50.1 | 52.7 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17516 of 1048576 (1.670%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 223.8, all cells median 224.6
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.332 | 78.86 | 2.03% |
| default | 16 | 6.332 | 78.86 | 2.03% |
| default | 36 | 6.332 | 78.86 | 2.03% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__urban__singapore

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 72555828029304.625 | 109.00 | 55.21% | 1.629% | 0.00091 | 0.3 | 0.3 | 0.3 |
| affine | 0.064 | 0.84 | 0.00% | 0.000% | -0.658 | 1.11 | -0.544 | 0.944 |
| wave_60m | 11.721 | 147.92 | 38.74% | 3.337% | -24.2 | 24.4 | 0.152 | 0.248 |
| analytic | 6.166 | 74.02 | 2.06% | 1.688% | -8.88e+03 | 5.92e+03 | -50.9 | 52.1 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17702 of 1048576 (1.688%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 222.5, all cells median 223.0
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.166 | 74.02 | 2.06% |
| default | 16 | 6.166 | 74.02 | 2.06% |
| default | 36 | 6.166 | 74.02 | 2.06% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__urban__sydney

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 72195297451090.812 | 115.00 | 55.22% | 1.637% | 0.000107 | 0.3 | 0.3 | 0.3 |
| affine | 0.070 | 0.93 | 0.01% | 0.000% | -0.71 | 1.06 | -0.504 | 0.904 |
| wave_60m | 11.795 | 148.60 | 39.15% | 3.323% | -24 | 24.2 | 0.152 | 0.248 |
| analytic | 6.485 | 77.96 | 2.04% | 1.693% | -6.37e+03 | 5.99e+03 | -50.5 | 52.8 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17749 of 1048576 (1.693%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 223.5, all cells median 222.3
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.485 | 77.96 | 2.04% |
| default | 16 | 6.485 | 77.96 | 2.04% |
| default | 36 | 6.485 | 77.96 | 2.04% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__urban__tokyo

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 71400340436468.562 | 109.00 | 55.24% | 1.616% | 0.00228 | 0.3 | 0.3 | 0.3 |
| affine | 0.069 | 0.90 | 0.00% | 0.000% | -0.645 | 1.03 | -0.5 | 0.9 |
| wave_60m | 11.674 | 165.66 | 37.68% | 3.413% | -25 | 25.2 | 0.15 | 0.25 |
| analytic | 6.454 | 73.16 | 2.00% | 1.693% | -7.5e+03 | 8.22e+03 | -51.8 | 52.5 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17753 of 1048576 (1.693%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 221.3, all cells median 221.2
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.454 | 73.16 | 2.00% |
| default | 16 | 6.454 | 73.16 | 2.00% |
| default | 36 | 6.454 | 73.16 | 2.00% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__water__caspian_sea

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 74919130685109.344 | 115.00 | 55.30% | 1.626% | 0.000407 | 0.3 | 0.3 | 0.3 |
| affine | 0.073 | 0.97 | 0.01% | 0.001% | -0.677 | 1.06 | -0.491 | 0.891 |
| wave_60m | 11.925 | 150.26 | 39.06% | 3.318% | -24.1 | 24.3 | 0.152 | 0.248 |
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

## region__water__great_salt_lake

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 69840851987873.891 | 109.00 | 55.23% | 1.635% | 0.000146 | 0.3 | 0.3 | 0.3 |
| affine | 0.070 | 0.92 | 0.01% | 0.000% | -0.652 | 1.05 | -0.492 | 0.892 |
| wave_60m | 11.820 | 154.59 | 37.95% | 3.451% | -24.8 | 25 | 0.15 | 0.25 |
| analytic | 6.202 | 75.85 | 2.03% | 1.703% | -6.61e+03 | 6.07e+03 | -52.4 | 51.8 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17853 of 1048576 (1.703%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 215.5, all cells median 215.8
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.202 | 75.85 | 2.03% |
| default | 16 | 6.202 | 75.85 | 2.03% |
| default | 36 | 6.202 | 75.85 | 2.03% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__water__lake_baikal

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 71825157948511.500 | 115.00 | 55.36% | 1.629% | 0.00082 | 0.3 | 0.3 | 0.3 |
| affine | 0.085 | 1.18 | 0.00% | 0.001% | -0.617 | 1.05 | -0.489 | 0.889 |
| wave_60m | 11.965 | 158.49 | 38.28% | 3.409% | -24.6 | 24.8 | 0.151 | 0.249 |
| analytic | 6.363 | 80.14 | 2.04% | 1.677% | -7.11e+03 | 6.92e+03 | -52.3 | 51.5 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17583 of 1048576 (1.677%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 201.6, all cells median 201.6
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.363 | 80.14 | 2.04% |
| default | 16 | 6.363 | 80.14 | 2.04% |
| default | 36 | 6.363 | 80.14 | 2.04% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.057, p99.9/std 0.21, overshoot 0.00%, >5 std 0.000%

## region__water__lake_constance

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 72630833790735.328 | 115.00 | 55.22% | 1.640% | 0.000361 | 0.3 | 0.3 | 0.3 |
| affine | 0.090 | 1.20 | 0.00% | 0.001% | -0.638 | 0.954 | -0.492 | 0.891 |
| wave_60m | 11.787 | 158.88 | 38.07% | 3.427% | -24.7 | 24.9 | 0.151 | 0.249 |
| analytic | 6.288 | 75.06 | 2.02% | 1.670% | -6.36e+03 | 7.2e+03 | -53 | 52.2 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17509 of 1048576 (1.670%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 193.7, all cells median 194.5
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.288 | 75.06 | 2.02% |
| default | 16 | 6.288 | 75.06 | 2.02% |
| default | 36 | 6.288 | 75.06 | 2.02% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.057, p99.9/std 0.21, overshoot 0.00%, >5 std 0.000%

## region__water__lake_geneva

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 74917893402128.719 | 115.00 | 55.34% | 1.638% | 0.000846 | 0.3 | 0.3 | 0.3 |
| affine | 0.087 | 1.18 | 0.01% | 0.001% | -0.684 | 1.06 | -0.485 | 0.885 |
| wave_60m | 11.961 | 154.15 | 51.20% | 3.015% | -21.8 | 22 | 0.156 | 0.244 |
| analytic | 6.293 | 78.00 | 2.08% | 1.660% | -5.7e+03 | 7.38e+03 | -51.7 | 50.9 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17402 of 1048576 (1.660%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 184.6, all cells median 185.8
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.293 | 78.00 | 2.08% |
| default | 16 | 6.293 | 78.00 | 2.08% |
| default | 36 | 6.293 | 78.00 | 2.08% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.057, p99.9/std 0.21, overshoot 0.00%, >5 std 0.000%

## region__water__lake_michigan

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 74208575633056.891 | 115.00 | 55.27% | 1.637% | 0.00015 | 0.3 | 0.3 | 0.3 |
| affine | 0.092 | 1.27 | 0.00% | 0.001% | -0.659 | 1.05 | -0.492 | 0.893 |
| wave_60m | 11.880 | 149.50 | 44.70% | 3.141% | -22.6 | 22.9 | 0.155 | 0.245 |
| analytic | 6.354 | 79.63 | 2.06% | 1.651% | -6.61e+03 | 6.93e+03 | -51.1 | 51.5 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17314 of 1048576 (1.651%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 184.8, all cells median 184.7
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.354 | 79.63 | 2.06% |
| default | 16 | 6.354 | 79.63 | 2.06% |
| default | 36 | 6.354 | 79.63 | 2.06% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.056, p99.9/std 0.21, overshoot 0.00%, >5 std 0.000%

## region__water__lake_tahoe

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 76541270577400.297 | 115.00 | 55.22% | 1.627% | 0.000552 | 0.3 | 0.3 | 0.3 |
| affine | 0.072 | 0.96 | 0.00% | 0.000% | -0.627 | 1.03 | -0.494 | 0.895 |
| wave_60m | 11.950 | 152.21 | 51.65% | 3.014% | -21.8 | 22 | 0.156 | 0.244 |
| analytic | 6.625 | 80.53 | 2.04% | 1.700% | -6.2e+03 | 8.19e+03 | -51.3 | 52 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17824 of 1048576 (1.700%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 218.1, all cells median 218.3
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.625 | 80.53 | 2.04% |
| default | 16 | 6.625 | 80.53 | 2.04% |
| default | 36 | 6.625 | 80.53 | 2.04% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__water__lake_taupo

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 77313920801550.922 | 115.00 | 55.15% | 1.635% | 0.000146 | 0.3 | 0.3 | 0.3 |
| affine | 0.073 | 0.96 | 0.00% | 0.000% | -0.629 | 0.966 | -0.496 | 0.895 |
| wave_60m | 12.224 | 162.08 | 38.25% | 3.369% | -24.5 | 24.8 | 0.151 | 0.249 |
| analytic | 6.461 | 77.80 | 1.99% | 1.696% | -7.26e+03 | 6.42e+03 | -52.4 | 52.9 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17780 of 1048576 (1.696%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 219.5, all cells median 218.7
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.461 | 77.80 | 1.99% |
| default | 16 | 6.461 | 77.80 | 1.99% |
| default | 36 | 6.461 | 77.80 | 1.99% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__water__lake_titicaca

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 73261062198439.359 | 115.00 | 55.25% | 1.632% | 0.000443 | 0.3 | 0.3 | 0.3 |
| affine | 0.066 | 0.86 | 0.00% | 0.000% | -0.737 | 1.11 | -0.534 | 0.934 |
| wave_60m | 11.746 | 149.93 | 47.21% | 3.088% | -22.3 | 22.4 | 0.155 | 0.245 |
| analytic | 6.203 | 75.24 | 2.03% | 1.694% | -5.19e+03 | 6.67e+03 | -51 | 53.4 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17760 of 1048576 (1.694%); their RMSE alone accounts for 98.9% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 223.6, all cells median 224.1
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.203 | 75.24 | 2.03% |
| default | 16 | 6.203 | 75.24 | 2.03% |
| default | 36 | 6.203 | 75.24 | 2.03% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%

## region__water__lake_victoria

### Operator inventory

```
Npt: 16
sigma_m: 6.217588859870031
```

constructor signature: `(*args, Npt: int = 16, area: Optional[~T_Array] = None, **kwargs)`
module: `/mnt/c/Users/jeanm/Claude/healpix-resample/healpix_resample/__init__.py` version `?`

### Probe fields (radial bicubic, default parameters)

| field | RMSE/std | p99.9/std | overshoot | >5 std | out min | out max | src min | src max |
|---|---|---|---|---|---|---|---|---|
| constant | 72654063761372.453 | 115.00 | 55.14% | 1.643% | 0.000837 | 0.3 | 0.3 | 0.3 |
| affine | 0.066 | 0.88 | 0.00% | 0.000% | -0.653 | 1.14 | -0.544 | 0.944 |
| wave_60m | 12.181 | 153.43 | 40.09% | 3.259% | -23.5 | 23.7 | 0.153 | 0.247 |
| analytic | 6.493 | 80.35 | 2.05% | 1.688% | -6.67e+03 | 5.62e+03 | -50.3 | 52.7 |

### Worst cells on the analytic field

- cells with |error| > 5 std: 17700 of 1048576 (1.688%); their RMSE alone accounts for 99.0% of the total squared error
- distance to the tile's pixel box (native pixels): bad cells median 223.3, all cells median 223.0
- share of level-19 parents with all four children bad: 0.00% (random placement would give ~0.0000%)
- operator exposes no norm_col / norm_col_raw: inspect healpix_resample/bicubic.py by hand for the guard value

### Parameter sweep on the analytic field

| sigma_m | Npt | RMSE/std | p99.9/std | overshoot |
|---|---|---|---|---|
| default | default | 6.493 | 80.35 | 2.05% |
| default | 16 | 6.493 | 80.35 | 2.05% |
| default | 36 | 6.493 | 80.35 | 2.05% |

### Classical separable bicubic on the source raster (scipy order 3)

RMSE/std 0.058, p99.9/std 0.22, overshoot 0.00%, >5 std 0.000%
