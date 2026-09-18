# Paper numbers (common RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 3 targets
- 16-bit requantization floor: HPX 9.24e-06, UTM 6.23e-06 (the image codecs are fed a 16-bit tile, so their `lossless` points are lossless after that step, not bit-exact)

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| budget | 0.1 | 37 | 43.488 | [32.061, 220.952] | 0.0000 | 0.0000 |
| budget | 0.05 | 39 | 4.584 | [3.104, 31.540] | 0.0000 | 0.0000 |
| budget | 0.02 | 40 | 1.341 | [0.936, 39.898] | 0.0000 | 0.0000 |

## Censoring (fraction of regions where the target is bracketed
on both grids; the rest are conservative bounds, not estimates)

- rmse 0.001: 10% bracketed
- rmse 0.003: 30% bracketed
- rmse 0.01: 62% bracketed
- budget 0.2: 50% bracketed
- budget 0.1: 92% bracketed
- budget 0.05: 98% bracketed
- budget 0.02: 100% bracketed

## Intrinsic gap (gamma=1) and critical overhead


## HEALPix regime (latitude)

- budget:0.1: Spearman rho=0.246 (p=0.142, Holm=0.142), belt 126.608 vs cap 159.029 (p=0.091, Holm=0.091)
    given the spectral slope beta: partial rho=0.246 (p=0.149, n=37); rho(beta, |lat|)=0.030 (p=0.862); rho(beta, ratio)=0.214
- budget:0.05: Spearman rho=0.440 (p=0.005, Holm=0.010), belt 7.939 vs cap 11.223 (p=0.016, Holm=0.032)
    given the spectral slope beta: partial rho=0.451 (p=0.005, n=39); rho(beta, |lat|)=0.050 (p=0.762); rho(beta, ratio)=0.343
- budget:0.02: Spearman rho=0.502 (p=0.001, Holm=0.003), belt 1.464 vs cap 2.078 (p=0.006, Holm=0.019)
    given the spectral slope beta: partial rho=0.522 (p=0.001, n=40); rho(beta, |lat|)=0.049 (p=0.763); rho(beta, ratio)=0.377

## Codec ranking (HEALPix side)

- 253 operating points, wins sum to the same number
- jpegxl_lossy: 223 wins / 253 evaluable, mean rank 1.25
- jpeg2000_lossy: 16 wins / 249 evaluable, mean rank 2.67
- dwt_zstd: 9 wins / 249 evaluable, mean rank 3.28
- sz3: 5 wins / 249 evaluable, mean rank 3.25

## HEALPix representation floor (common reference)

- median 4.444e-03, range [3.773e-04, 1.496e-02] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
