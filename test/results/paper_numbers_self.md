# Paper numbers (self RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 8 targets
- 16-bit requantization floor: HPX 6.29e-06, UTM 6.23e-06 (the image codecs are fed a 16-bit tile, so their `lossless` points are lossless after that step, not bit-exact)

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| rmse | 0.0001 | 40 | 1.594 | [0.705, 1.857] | 0.0000 | 0.0000 |
| rmse | 0.0003 | 40 | 1.489 | [0.537, 1.811] | 0.0000 | 0.0000 |
| rmse | 0.001 | 40 | 1.249 | [0.345, 6.128] | 0.0026 | 0.0079 |
| rmse | 0.003 | 40 | 0.996 | [0.069, 5.034] | 0.7346 | 0.7346 |
| rmse | 0.01 | 40 | 0.828 | [0.536, 2.582] | 0.0031 | 0.0079 |
| budget | 0.1 | 40 | 4.060 | [0.867, 8.961] | 0.0000 | 0.0000 |
| budget | 0.05 | 40 | 1.333 | [0.906, 2.316] | 0.0000 | 0.0000 |
| budget | 0.02 | 40 | 0.829 | [0.510, 1.418] | 0.0000 | 0.0001 |

## Censoring (fraction of regions where the target is bracketed
on both grids; the rest are conservative bounds, not estimates)

- rmse 0.0001: 100% bracketed
- rmse 0.0003: 100% bracketed
- rmse 0.001: 100% bracketed
- rmse 0.003: 90% bracketed
- rmse 0.01: 85% bracketed
- budget 0.4: 0% bracketed
- budget 0.2: 18% bracketed
- budget 0.1: 98% bracketed
- budget 0.05: 100% bracketed
- budget 0.02: 100% bracketed

## Intrinsic gap (gamma=1) and critical overhead

- RMSE 0.0001: gap +107%, gamma* = 2.07
- RMSE 0.0003: gap +94%, gamma* = 1.94
- RMSE 0.001: gap +62%, gamma* = 1.62
- RMSE 0.003: gap +29%, gamma* = 1.29
- RMSE 0.01: gap +8%, gamma* = 1.08

## HEALPix regime (latitude)

- rmse:0.0001: Spearman rho=-0.050 (p=0.760, Holm=1.000), belt 2.092 vs cap 2.062 (p=0.455, Holm=1.000)
- rmse:0.0003: Spearman rho=-0.040 (p=0.807, Holm=1.000), belt 1.904 vs cap 1.960 (p=0.407, Holm=1.000)
- rmse:0.001: Spearman rho=0.025 (p=0.879, Holm=1.000), belt 1.683 vs cap 1.596 (p=0.438, Holm=1.000)
- rmse:0.003: Spearman rho=0.139 (p=0.392, Holm=1.000), belt 1.216 vs cap 1.327 (p=0.187, Holm=0.749)
- rmse:0.01: Spearman rho=0.318 (p=0.046, Holm=0.275), belt 0.967 vs cap 1.275 (p=0.004, Holm=0.031)
- budget:0.1: Spearman rho=0.062 (p=0.703, Holm=1.000), belt 12.365 vs cap 14.223 (p=0.146, Holm=0.729)
- budget:0.05: Spearman rho=0.372 (p=0.018, Holm=0.127), belt 2.678 vs cap 2.836 (p=0.007, Holm=0.041)
- budget:0.02: Spearman rho=0.499 (p=0.001, Holm=0.008), belt 0.947 vs cap 1.228 (p=0.001, Holm=0.007)

## Codec ranking (HEALPix side)

- 400 operating points, wins sum to the same number
- jpeg2000_lossy: 165 wins / 400 evaluable, mean rank 1.86
- jpegxl_lossy: 150 wins / 400 evaluable, mean rank 2.23
- sz3: 37 wins / 400 evaluable, mean rank 2.72
- zfp: 35 wins / 333 evaluable, mean rank 4.49
- dwt_zstd: 13 wins / 396 evaluable, mean rank 3.60

## HEALPix representation floor (common reference)

- median 4.038e-03, range [1.340e-04, 1.370e-02] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
