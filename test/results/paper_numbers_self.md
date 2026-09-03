# Paper numbers (self RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 9 targets
- 16-bit requantization floor: HPX 6.29e-06, UTM 6.23e-06 (the image codecs are fed a 16-bit tile, so their `lossless` points are lossless after that step, not bit-exact)

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| rmse | 0.0001 | 40 | 1.594 | [0.705, 1.857] | 0.0000 | 0.0000 |
| rmse | 0.0003 | 40 | 1.489 | [0.537, 1.811] | 0.0000 | 0.0000 |
| rmse | 0.001 | 40 | 1.249 | [0.345, 6.128] | 0.0009 | 0.0019 |
| rmse | 0.003 | 36 | 0.944 | [0.065, 1.587] | 0.2325 | 0.2325 |
| rmse | 0.01 | 34 | 0.787 | [0.235, 1.543] | 0.0002 | 0.0005 |
| budget | 0.2 | 38 | 14.041 | [10.178, 29.227] | 0.0000 | 0.0000 |
| budget | 0.1 | 39 | 4.061 | [2.554, 7.161] | 0.0000 | 0.0000 |
| budget | 0.05 | 40 | 1.335 | [1.019, 7.589] | 0.0000 | 0.0000 |
| budget | 0.02 | 40 | 0.829 | [0.510, 1.418] | 0.0000 | 0.0001 |

## Intrinsic gap (gamma=1) and critical overhead

- RMSE 0.0001: gap +107%, gamma* = 2.07
- RMSE 0.0003: gap +94%, gamma* = 1.94
- RMSE 0.001: gap +62%, gamma* = 1.62
- RMSE 0.003: gap +23%, gamma* = 1.23
- RMSE 0.01: gap +2%, gamma* = 1.02

## HEALPix regime (latitude)

- rmse:0.0001: Spearman rho=-0.050 (p=0.758), belt 2.092 vs cap 2.062 (p=0.455)
- rmse:0.0003: Spearman rho=-0.041 (p=0.802), belt 1.904 vs cap 1.960 (p=0.407)
- rmse:0.001: Spearman rho=0.019 (p=0.906), belt 1.683 vs cap 1.596 (p=0.596)
- rmse:0.003: Spearman rho=0.088 (p=0.609), belt 1.130 vs cap 1.273 (p=0.220)
- rmse:0.01: Spearman rho=0.209 (p=0.234), belt 0.925 vs cap 1.191 (p=0.015)
- budget:0.2: Spearman rho=0.893 (p=0.007), belt nan vs cap nan (p=nan)
- budget:0.1: Spearman rho=0.186 (p=0.256), belt 12.341 vs cap 14.223 (p=0.088)
- budget:0.05: Spearman rho=0.362 (p=0.022), belt 2.678 vs cap 2.836 (p=0.010)
- budget:0.02: Spearman rho=0.497 (p=0.001), belt 0.947 vs cap 1.228 (p=0.001)

## Codec ranking (HEALPix side)

- 383 operating points, wins sum to the same number
- jpeg2000_lossy: 162 wins / 317 evaluable, mean rank 1.56
- jpegxl_lossy: 146 wins / 327 evaluable, mean rank 1.93
- zfp: 33 wins / 269 evaluable, mean rank 4.42
- sz3: 30 wins / 346 evaluable, mean rank 2.71
- dwt_zstd: 12 wins / 333 evaluable, mean rank 3.51

## HEALPix representation floor (common reference)

- median 8.116e-02, range [1.992e-03, 2.093e-01] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
