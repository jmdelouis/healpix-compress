# Paper numbers (common RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 3 targets
- 16-bit requantization floor: HPX 5.79e-06, UTM 6.23e-06 (the image codecs are fed a 16-bit tile, so their `lossless` points are lossless after that step, not bit-exact)

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| budget | 0.1 | 34 | 32.180 | [18.391, 162.914] | 0.0000 | 0.0000 |
| budget | 0.05 | 36 | 3.855 | [2.401, 23.338] | 0.0000 | 0.0000 |
| budget | 0.02 | 37 | 1.437 | [1.084, 24.153] | 0.0000 | 0.0000 |

## Censoring (fraction of regions where the target is bracketed
on both grids; the rest are conservative bounds, not estimates)

- rmse 0.0003: 8% bracketed
- rmse 0.001: 15% bracketed
- rmse 0.003: 30% bracketed
- rmse 0.01: 70% bracketed
- budget 0.2: 55% bracketed
- budget 0.1: 85% bracketed
- budget 0.05: 90% bracketed
- budget 0.02: 92% bracketed

## Intrinsic gap (gamma=1) and critical overhead


## HEALPix regime (latitude)

- budget:0.1: Spearman rho=0.241 (p=0.170, Holm=0.180), belt 100.421 vs cap 126.380 (p=0.077, Holm=0.184)
- budget:0.05: Spearman rho=0.287 (p=0.090, Holm=0.180), belt 7.176 vs cap 8.514 (p=0.067, Holm=0.184)
- budget:0.02: Spearman rho=0.328 (p=0.048, Holm=0.143), belt 1.631 vs cap 2.060 (p=0.061, Holm=0.184)

## Codec ranking (HEALPix side)

- 247 operating points, wins sum to the same number
- jpeg2000_lossy: 138 wins / 247 evaluable, mean rank 1.61
- jpegxl_lossy: 49 wins / 247 evaluable, mean rank 3.01
- sz3: 33 wins / 247 evaluable, mean rank 2.54
- zfp: 16 wins / 227 evaluable, mean rank 4.36
- dwt_zstd: 11 wins / 247 evaluable, mean rank 3.43

## HEALPix representation floor (common reference)

- median 4.038e-03, range [1.340e-04, 1.370e-02] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
