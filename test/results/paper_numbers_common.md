# Paper numbers (common RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 3 targets
- 16-bit requantization floor: HPX 5.79e-06, UTM 6.23e-06 (the image codecs are fed a 16-bit tile, so their `lossless` points are lossless after that step, not bit-exact)

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| budget | 0.1 | 37 | 32.297 | [18.391, 2971.193] | 0.0000 | 0.0000 |
| budget | 0.05 | 37 | 3.882 | [2.401, 467.438] | 0.0000 | 0.0000 |
| budget | 0.02 | 37 | 1.437 | [1.084, 24.153] | 0.0000 | 0.0000 |

## Censoring (fraction of regions where the target is bracketed
on both grids; the rest are conservative bounds, not estimates)

- rmse 0.0003: 8% bracketed
- rmse 0.001: 20% bracketed
- rmse 0.003: 28% bracketed
- rmse 0.01: 72% bracketed
- budget 0.4: 0% bracketed
- budget 0.2: 18% bracketed
- budget 0.1: 90% bracketed
- budget 0.05: 92% bracketed
- budget 0.02: 92% bracketed

## Intrinsic gap (gamma=1) and critical overhead


## HEALPix regime (latitude)

- budget:0.1: Spearman rho=nan (p=nan, Holm=nan), belt nan vs cap nan (p=nan, Holm=nan)
- budget:0.05: Spearman rho=nan (p=nan, Holm=nan), belt nan vs cap nan (p=nan, Holm=nan)
- budget:0.02: Spearman rho=nan (p=nan, Holm=nan), belt nan vs cap nan (p=nan, Holm=nan)

## Codec ranking (HEALPix side)

- 247 operating points, wins sum to the same number
- jpeg2000_lossy: 121 wins / 247 evaluable, mean rank 1.87
- jpegxl_lossy: 59 wins / 247 evaluable, mean rank 2.89
- sz3: 40 wins / 247 evaluable, mean rank 2.46
- zfp: 16 wins / 185 evaluable, mean rank 4.21
- dwt_zstd: 11 wins / 247 evaluable, mean rank 3.38

## HEALPix representation floor (common reference)

- median 4.038e-03, range [1.340e-04, 1.370e-02] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
