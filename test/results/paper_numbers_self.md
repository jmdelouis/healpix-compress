# Paper numbers (self RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 8 targets
- 16-bit requantization floor: HPX 5.79e-06, UTM 6.23e-06 (the image codecs are fed a 16-bit tile, so their `lossless` points are lossless after that step, not bit-exact)

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| rmse | 0.0001 | 40 | 1.587 | [0.705, 1.740] | 0.0000 | 0.0000 |
| rmse | 0.0003 | 40 | 1.489 | [0.537, 1.704] | 0.0000 | 0.0000 |
| rmse | 0.001 | 40 | 1.227 | [0.345, 6.128] | 0.0039 | 0.0077 |
| rmse | 0.003 | 40 | 0.962 | [0.069, 5.034] | 0.4929 | 0.4929 |
| rmse | 0.01 | 40 | 0.812 | [0.536, 2.582] | 0.0007 | 0.0021 |
| budget | 0.1 | 40 | 4.016 | [0.867, 8.961] | 0.0000 | 0.0000 |
| budget | 0.05 | 40 | 1.333 | [0.906, 2.251] | 0.0000 | 0.0000 |
| budget | 0.02 | 40 | 0.829 | [0.510, 1.418] | 0.0000 | 0.0000 |

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

- RMSE 0.0001: gap +106%, gamma* = 2.06
- RMSE 0.0003: gap +94%, gamma* = 1.94
- RMSE 0.001: gap +60%, gamma* = 1.60
- RMSE 0.003: gap +25%, gamma* = 1.25
- RMSE 0.01: gap +6%, gamma* = 1.06

## HEALPix regime (latitude)

- rmse:0.0001: Spearman rho=nan (p=nan, Holm=nan), belt nan vs cap nan (p=nan, Holm=nan)
- rmse:0.0003: Spearman rho=nan (p=nan, Holm=nan), belt nan vs cap nan (p=nan, Holm=nan)
- rmse:0.001: Spearman rho=nan (p=nan, Holm=nan), belt nan vs cap nan (p=nan, Holm=nan)
- rmse:0.003: Spearman rho=nan (p=nan, Holm=nan), belt nan vs cap nan (p=nan, Holm=nan)
- rmse:0.01: Spearman rho=nan (p=nan, Holm=nan), belt nan vs cap nan (p=nan, Holm=nan)
- budget:0.1: Spearman rho=nan (p=nan, Holm=nan), belt nan vs cap nan (p=nan, Holm=nan)
- budget:0.05: Spearman rho=nan (p=nan, Holm=nan), belt nan vs cap nan (p=nan, Holm=nan)
- budget:0.02: Spearman rho=nan (p=nan, Holm=nan), belt nan vs cap nan (p=nan, Holm=nan)

## Codec ranking (HEALPix side)

- 400 operating points, wins sum to the same number
- jpeg2000_lossy: 166 wins / 400 evaluable, mean rank 1.85
- jpegxl_lossy: 149 wins / 400 evaluable, mean rank 2.25
- sz3: 37 wins / 400 evaluable, mean rank 2.71
- zfp: 35 wins / 333 evaluable, mean rank 4.48
- dwt_zstd: 13 wins / 395 evaluable, mean rank 3.60

## HEALPix representation floor (common reference)

- median 4.038e-03, range [1.340e-04, 1.370e-02] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
