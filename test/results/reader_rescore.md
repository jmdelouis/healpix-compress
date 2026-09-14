# The common-reference verdict under two readers

- floors: nearest 4.038e-03, bilinear 1.940e-03
- a target is reported when bracketed on both grids for >= 75% of regions

| target | nearest n | nearest median | p_Holm | bilinear n | bilinear median | p_Holm | gamma* (bilinear) |
|---|---|---|---|---|---|---|---|
| budget 0.02 | 37 | 1.437 | 4.4e-11 | 40 | 1.119 | 4.5e-06 | 1.47 |
| budget 0.05 | 37 | 3.882 | 4.4e-11 | 40 | 2.501 | 7.3e-12 | 5.16 |
| budget 0.1 | 37 | 32.297 | 4.4e-11 | 40 | 19.309 | 7.3e-12 | 62.86 |
| budget 0.2 | 37 | 1212.270 (censored) | nan | 40 | 698.804 (censored) | nan | --- |
| budget 0.4 | 37 | 6088.537 (censored) | nan | 40 | 2961.469 (censored) | nan | --- |
| rmse 0.0001 | 0 | nan (censored) | nan | 0 | nan (censored) | nan | --- |
| rmse 0.0003 | 3 | 2.082 (censored) | nan | 3 | 1.272 (censored) | nan | --- |
| rmse 0.001 | 8 | 1.534 (censored) | nan | 12 | 1.377 (censored) | nan | --- |
| rmse 0.003 | 15 | 1.324 (censored) | nan | 27 | 1.256 (censored) | nan | --- |
| rmse 0.01 | 36 | 1.163 (censored) | nan | 40 | 1.015 | 0.55 | 1.32 |

- at the 2% budget the interpolating read moves the ratio from 1.44 to 1.12; native UTM still leads (p_Holm 4.5e-06)
- at RMSE 0.01 the two grids reach parity under the interpolating read: 1.015, p_Holm 0.55, n=40; intrinsic ratio 1.32, so gamma* = 1.32 sits at the illustrative 1.30