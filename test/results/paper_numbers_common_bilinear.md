# Paper numbers (common RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 4 targets
- 16-bit requantization floor: HPX 5.79e-06, UTM 6.23e-06 (the image codecs are fed a 16-bit tile, so their `lossless` points are lossless after that step, not bit-exact)

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| rmse | 0.01 | 32 | 0.926 | [0.203, 1.449] | 0.0011 | 0.0011 |
| budget | 0.1 | 37 | 19.000 | [12.993, 106.614] | 0.0000 | 0.0000 |
| budget | 0.05 | 39 | 2.496 | [1.797, 15.293] | 0.0000 | 0.0000 |
| budget | 0.02 | 40 | 1.119 | [0.896, 15.142] | 0.0000 | 0.0000 |

## Censoring (fraction of regions where the target is bracketed
on both grids; the rest are conservative bounds, not estimates)

- rmse 0.0003: 8% bracketed
- rmse 0.001: 25% bracketed
- rmse 0.003: 60% bracketed
- rmse 0.01: 80% bracketed
- budget 0.2: 65% bracketed
- budget 0.1: 92% bracketed
- budget 0.05: 98% bracketed
- budget 0.02: 100% bracketed

## Intrinsic gap (gamma=1) and critical overhead

- RMSE 0.01: gap +20%, gamma* = 1.20

## HEALPix regime (latitude)

- rmse:0.01: Spearman rho=0.093 (p=0.611, Holm=0.611), belt 1.139 vs cap 1.211 (p=0.326, Holm=0.326)
- budget:0.1: Spearman rho=0.468 (p=0.004, Holm=0.007), belt 52.106 vs cap 73.660 (p=0.001, Holm=0.002)
- budget:0.05: Spearman rho=0.538 (p=0.000, Holm=0.002), belt 4.459 vs cap 5.867 (p=0.001, Holm=0.002)
- budget:0.02: Spearman rho=0.520 (p=0.001, Holm=0.002), belt 1.305 vs cap 1.602 (p=0.002, Holm=0.003)

## Codec ranking (HEALPix side)

- 282 operating points, wins sum to the same number
- jpeg2000_lossy: 162 wins / 282 evaluable, mean rank 1.68
- jpegxl_lossy: 48 wins / 282 evaluable, mean rank 2.96
- sz3: 44 wins / 282 evaluable, mean rank 2.54
- dwt_zstd: 14 wins / 282 evaluable, mean rank 3.36
- zfp: 14 wins / 261 evaluable, mean rank 4.42

## HEALPix representation floor (common reference)

- median 1.940e-03, range [1.211e-04, 6.738e-03] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
