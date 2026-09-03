# Paper numbers (self RMSE reference)

- regions: **7**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 9 targets

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| rmse | 0.0001 | 7 | 1.588 | [1.258, 1.715] | 0.0156 | 0.1406 |
| rmse | 0.0003 | 7 | 1.467 | [0.613, 1.676] | 0.2969 | 0.5938 |
| rmse | 0.001 | 7 | 1.239 | [0.938, 1.893] | 0.0312 | 0.1562 |
| rmse | 0.003 | 5 | 0.801 | [0.540, 1.237] | 0.1875 | 0.5625 |
| rmse | 0.01 | 5 | 0.910 | [0.683, 2.019] | 1.0000 | 1.0000 |
| budget | 0.2 | 6 | 13.445 | [8.652, 33.507] | 0.0312 | 0.1562 |
| budget | 0.1 | 7 | 4.167 | [2.884, 4.425] | 0.0156 | 0.1406 |
| budget | 0.05 | 7 | 1.315 | [1.258, 1.440] | 0.0156 | 0.1406 |
| budget | 0.02 | 7 | 0.856 | [0.728, 0.956] | 0.0156 | 0.1406 |

## Intrinsic gap (gamma=1) and critical overhead

- RMSE 0.0001: gap +106%, gamma* = 2.06
- RMSE 0.0003: gap +91%, gamma* = 1.91
- RMSE 0.001: gap +61%, gamma* = 1.61
- RMSE 0.003: gap +4%, gamma* = 1.04
- RMSE 0.01: gap +18%, gamma* = 1.18

## HEALPix regime (latitude)


## Codec ranking (HEALPix side)

- 65 operating points, wins sum to the same number
- jpeg2000_lossy: 26 wins / 54 evaluable, mean rank 1.56
- jpegxl_lossy: 22 wins / 28 evaluable, mean rank 1.21
- sz3: 8 wins / 56 evaluable, mean rank 2.57
- zfp: 6 wins / 48 evaluable, mean rank 3.81
- dwt_zstd: 3 wins / 54 evaluable, mean rank 3.06
