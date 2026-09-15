# Paper numbers (common RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 3 targets
- 16-bit requantization floor: HPX 5.79e-06, UTM 6.23e-06 (the image codecs are fed a 16-bit tile, so their `lossless` points are lossless after that step, not bit-exact)

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| budget | 0.1 | 37 | 19.000 | [12.993, 106.614] | 0.0000 | 0.0000 |
| budget | 0.05 | 39 | 2.496 | [1.797, 15.293] | 0.0000 | 0.0000 |
| budget | 0.02 | 40 | 1.119 | [0.896, 15.142] | 0.0000 | 0.0000 |

## Censoring (fraction of regions where the target is bracketed
on both grids; the rest are conservative bounds, not estimates)

- rmse 0.0003: 8% bracketed
- rmse 0.001: 25% bracketed
- rmse 0.003: 60% bracketed
- rmse 0.01: 72% bracketed
- budget 0.2: 65% bracketed
- budget 0.1: 92% bracketed
- budget 0.05: 98% bracketed
- budget 0.02: 100% bracketed

## Intrinsic gap (gamma=1) and critical overhead


## HEALPix regime (latitude)

- budget:0.1: Spearman rho=0.468 (p=0.004, Holm=0.004), belt 52.106 vs cap 73.660 (p=0.001, Holm=0.002)
- budget:0.05: Spearman rho=0.538 (p=0.000, Holm=0.001), belt 4.459 vs cap 5.867 (p=0.001, Holm=0.002)
- budget:0.02: Spearman rho=0.520 (p=0.001, Holm=0.001), belt 1.305 vs cap 1.602 (p=0.002, Holm=0.002)

## Codec ranking (HEALPix side)

- 282 operating points, wins sum to the same number
- jpeg2000_lossy: 158 wins / 282 evaluable, mean rank 1.72
- jpegxl_lossy: 48 wins / 282 evaluable, mean rank 2.95
- sz3: 45 wins / 282 evaluable, mean rank 2.53
- dwt_zstd: 17 wins / 282 evaluable, mean rank 3.33
- zfp: 14 wins / 261 evaluable, mean rank 4.42

## HEALPix representation floor (common reference)

- median 1.940e-03, range [1.211e-04, 6.738e-03] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
