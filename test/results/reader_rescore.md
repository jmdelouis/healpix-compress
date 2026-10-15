# The common-reference verdict under two readers

- floors: nearest 4.038e-03, bilinear 1.988e-03
- a target is reported when bracketed on both grids for >= 75% of regions

| target | nearest n | nearest median | p_Holm | bilinear n | bilinear median | p_Holm | gamma* (bilinear) |
|---|---|---|---|---|---|---|---|
| budget 0.02 | 40 | 1.500 | 7.3e-12 | 40 | 1.119 | 4.9e-06 | 1.47 |
| budget 0.05 | 39 | 4.046 | 1.1e-11 | 39 | 2.496 | 1.5e-11 | 5.12 |
| budget 0.1 | 38 | 32.517 | 1.5e-11 | 38 | 19.309 | 2.2e-11 | 62.69 |
| budget 0.2 | 24 | 1040.681 (censored) | nan | 25 | 631.193 (censored) | nan | --- |
| budget 0.4 | 0 | nan (censored) | nan | 0 | nan (censored) | nan | --- |
| rmse 0.0001 | 0 | nan (censored) | nan | 0 | nan (censored) | nan | --- |
| rmse 0.0003 | 2 | 2.014 (censored) | nan | 3 | 1.496 (censored) | nan | --- |
| rmse 0.001 | 7 | 1.495 (censored) | nan | 10 | 1.253 (censored) | nan | --- |
| rmse 0.003 | 12 | 1.001 (censored) | nan | 23 | 1.096 (censored) | nan | --- |
| rmse 0.01 | 31 | 1.094 | 0.18 | 32 | 0.926 | 0.0026 | 1.20 |

- at the 2% budget the interpolating read moves the ratio from 1.50 to 1.12; native UTM still leads (p_Holm 4.9e-06)
- at RMSE 0.01 the two grids reach parity under the interpolating read: 0.926, p_Holm 0.00, n=32; intrinsic ratio 1.20, so gamma* = 1.20 sits at the illustrative 1.30