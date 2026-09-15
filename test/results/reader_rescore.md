# The common-reference verdict under two readers

- floors: nearest 4.038e-03, bilinear 1.940e-03
- a target is reported when bracketed on both grids for >= 75% of regions

| target | nearest n | nearest median | p_Holm | bilinear n | bilinear median | p_Holm | gamma* (bilinear) |
|---|---|---|---|---|---|---|---|
| budget 0.02 | 40 | 1.278 | 5.5e-12 | 40 | 1.119 | 2.2e-06 | 1.47 |
| budget 0.05 | 39 | 3.172 | 7.3e-12 | 39 | 2.496 | 1.1e-11 | 5.12 |
| budget 0.1 | 37 | 24.879 | 1.5e-11 | 37 | 19.000 | 2.9e-11 | 64.44 |
| budget 0.2 | 26 | 806.865 (censored) | nan | 26 | 618.471 (censored) | nan | --- |
| budget 0.4 | 0 | nan (censored) | nan | 0 | nan (censored) | nan | --- |
| rmse 0.0001 | 0 | nan (censored) | nan | 0 | nan (censored) | nan | --- |
| rmse 0.0003 | 3 | 1.580 (censored) | nan | 3 | 1.272 (censored) | nan | --- |
| rmse 0.001 | 10 | 1.670 (censored) | nan | 10 | 1.253 (censored) | nan | --- |
| rmse 0.003 | 23 | 1.419 (censored) | nan | 24 | 1.129 (censored) | nan | --- |
| rmse 0.01 | 29 | 1.024 (censored) | nan | 29 | 0.926 (censored) | nan | --- |

- at the 2% budget the interpolating read moves the ratio from 1.28 to 1.12; native UTM still leads (p_Holm 2.2e-06)