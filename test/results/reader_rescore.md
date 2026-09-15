# The common-reference verdict under two readers

- floors: nearest 4.038e-03, bilinear 1.940e-03
- a target is reported when bracketed on both grids for >= 75% of regions

| target | nearest n | nearest median | p_Holm | bilinear n | bilinear median | p_Holm | gamma* (bilinear) |
|---|---|---|---|---|---|---|---|
| budget 0.02 | 37 | 1.437 | 4.4e-11 | 40 | 1.119 | 4.5e-06 | 1.47 |
| budget 0.05 | 36 | 3.855 | 5.8e-11 | 39 | 2.496 | 1.5e-11 | 5.12 |
| budget 0.1 | 34 | 32.180 | 1.2e-10 | 37 | 19.000 | 4.4e-11 | 64.44 |
| budget 0.2 | 22 | 979.195 (censored) | nan | 26 | 618.471 (censored) | nan | --- |
| budget 0.4 | 0 | nan (censored) | nan | 0 | nan (censored) | nan | --- |
| rmse 0.0001 | 0 | nan (censored) | nan | 0 | nan (censored) | nan | --- |
| rmse 0.0003 | 3 | 2.082 (censored) | nan | 3 | 1.272 (censored) | nan | --- |
| rmse 0.001 | 6 | 1.538 (censored) | nan | 10 | 1.253 (censored) | nan | --- |
| rmse 0.003 | 12 | 1.001 (censored) | nan | 24 | 1.129 (censored) | nan | --- |
| rmse 0.01 | 28 | 1.083 (censored) | nan | 32 | 0.926 | 0.0011 | 1.20 |

- at the 2% budget the interpolating read moves the ratio from 1.44 to 1.12; native UTM still leads (p_Holm 4.5e-06)
- at RMSE 0.01 the two grids reach parity under the interpolating read: 0.926, p_Holm 0.00, n=32; intrinsic ratio 1.20, so gamma* = 1.20 sits at the illustrative 1.30