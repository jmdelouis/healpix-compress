# Paper numbers (common RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 3 targets
- 16-bit requantization floor: HPX 6.75e-06, UTM 6.23e-06 (the image codecs are fed a 16-bit tile, so their `lossless` points are lossless after that step, not bit-exact)

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| budget | 0.1 | 37 | 24.879 | [15.692, 134.764] | 0.0000 | 0.0000 |
| budget | 0.05 | 39 | 3.172 | [2.231, 19.315] | 0.0000 | 0.0000 |
| budget | 0.02 | 40 | 1.278 | [1.010, 19.648] | 0.0000 | 0.0000 |

## Censoring (fraction of regions where the target is bracketed
on both grids; the rest are conservative bounds, not estimates)

- rmse 0.0003: 8% bracketed
- rmse 0.001: 25% bracketed
- rmse 0.003: 57% bracketed
- rmse 0.01: 72% bracketed
- budget 0.2: 65% bracketed
- budget 0.1: 92% bracketed
- budget 0.05: 98% bracketed
- budget 0.02: 100% bracketed

## Intrinsic gap (gamma=1) and critical overhead


## HEALPix regime (latitude)

- budget:0.1: Spearman rho=0.236 (p=0.160, Holm=0.160), belt 77.252 vs cap 98.393 (p=0.037, Holm=0.040)
- budget:0.05: Spearman rho=0.325 (p=0.043, Holm=0.087), belt 5.735 vs cap 7.337 (p=0.020, Holm=0.040)
- budget:0.02: Spearman rho=0.398 (p=0.011, Holm=0.033), belt 1.456 vs cap 1.827 (p=0.012, Holm=0.036)

## Codec ranking (HEALPix side)

- 282 operating points, wins sum to the same number
- jpeg2000_lossy: 164 wins / 282 evaluable, mean rank 1.68
- jpegxl_lossy: 49 wins / 282 evaluable, mean rank 3.00
- sz3: 32 wins / 282 evaluable, mean rank 2.51
- zfp: 21 wins / 261 evaluable, mean rank 4.36
- dwt_zstd: 16 wins / 282 evaluable, mean rank 3.40

## HEALPix representation floor (common reference)

- median 4.038e-03, range [1.340e-04, 1.370e-02] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
