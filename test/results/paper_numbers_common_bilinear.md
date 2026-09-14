# Paper numbers (common RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 4 targets
- 16-bit requantization floor: HPX 5.79e-06, UTM 6.23e-06 (the image codecs are fed a 16-bit tile, so their `lossless` points are lossless after that step, not bit-exact)

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| rmse | 0.01 | 40 | 1.015 | [0.402, 2.582] | 0.5537 | 0.5537 |
| budget | 0.1 | 40 | 19.309 | [12.993, 1862.679] | 0.0000 | 0.0000 |
| budget | 0.05 | 40 | 2.501 | [1.797, 293.043] | 0.0000 | 0.0000 |
| budget | 0.02 | 40 | 1.119 | [0.896, 15.142] | 0.0000 | 0.0000 |

## Censoring (fraction of regions where the target is bracketed
on both grids; the rest are conservative bounds, not estimates)

- rmse 0.0003: 8% bracketed
- rmse 0.001: 30% bracketed
- rmse 0.003: 57% bracketed
- rmse 0.01: 82% bracketed
- budget 0.4: 0% bracketed
- budget 0.2: 18% bracketed
- budget 0.1: 98% bracketed
- budget 0.05: 100% bracketed
- budget 0.02: 100% bracketed

## Intrinsic gap (gamma=1) and critical overhead

- RMSE 0.01: gap +32%, gamma* = 1.32

## HEALPix regime (latitude)


## Codec ranking (HEALPix side)

- 282 operating points, wins sum to the same number
- jpeg2000_lossy: 139 wins / 282 evaluable, mean rank 1.97
- jpegxl_lossy: 61 wins / 282 evaluable, mean rank 2.82
- sz3: 51 wins / 282 evaluable, mean rank 2.46
- dwt_zstd: 17 wins / 282 evaluable, mean rank 3.29
- zfp: 14 wins / 215 evaluable, mean rank 4.29

## HEALPix representation floor (common reference)

- median 1.940e-03, range [1.211e-04, 6.738e-03] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
