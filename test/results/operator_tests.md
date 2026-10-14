# Reconstruction operator B: what it does beyond a constant (V5)

- regions: 3, modes: nearest

| test | quantity | median |
|---|---|---|
| constant | RMSE (should be 0) | 0.00e+00 |
| ramp_x | RMSE, in cells of gradient | 0.319 |
| ramp_y | RMSE, in cells of gradient | 0.355 |
| band-limited field | RMSE, fraction of the field's std | 1.439e-01 |
| edge 0deg | RMSE across a step | 0.015 |
| edge 45deg | RMSE across a step | 0.021 |
| edge 90deg | RMSE across a step | 0.023 |
| edge 135deg | RMSE across a step | 0.017 |
| impulse | native pixels lit by one cell | 0.0 |
| impulse | centroid offset (m) | 2.72 |

## Transfer function

| wavelength (m) | response 0deg | 45deg | 90deg | 135deg | max &#124;shift&#124; (m) |
|---|---|---|---|---|---|
| 12.4 | 0.568 | 0.810 | 0.514 | 0.353 | 0.00 |
| 21.3 | 0.837 | 0.933 | 0.805 | 0.732 | 0.00 |
| 36.6 | 0.943 | 0.977 | 0.930 | 0.902 | 0.00 |
| 62.7 | 0.980 | 0.992 | 0.976 | 0.966 | 0.00 |
| 107.4 | 0.993 | 0.997 | 0.992 | 0.988 | 0.01 |
| 184.2 | 0.998 | 0.999 | 0.997 | 0.996 | 0.00 |
| 315.8 | 0.999 | 1.000 | 0.999 | 0.999 | 0.00 |
| 541.5 | 1.000 | 1.000 | 1.000 | 1.000 | 0.00 |
| 928.4 | 1.000 | 1.000 | 1.000 | 1.000 | 0.00 |
| 1591.7 | 1.000 | 1.000 | 1.000 | 1.000 | 0.00 |

- strongest directional anisotropy: **2.29x** between directions at 12.4 m wavelength
- largest registration shift over all wavelengths and directions: **0.01 m** (0.00 cells)

## Against the representation floor

The paper reports a representation floor of 4.038e-03 in median, about 10% of a scene's standard deviation. On a field of the imagery's own spectrum this operator contributes 1.439e-01 of the field's standard deviation, i.e. 14.39% against the floor's ~10%. The floor is NOT clearly dominated by the carrier grid: the operator accounts for a comparable share, and Section V-E's conclusion would need restating.
