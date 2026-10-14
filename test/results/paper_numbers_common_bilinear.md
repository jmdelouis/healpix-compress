# Paper numbers (common RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 4 targets
- 16-bit requantization floor: HPX 5.79e-06, UTM 6.47e-06 (the image codecs are fed a 16-bit tile, so their `lossless` points are lossless after that step, not bit-exact)

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| rmse | 0.01 | 32 | 0.926 | [0.203, 1.449] | 0.0026 | 0.0026 |
| budget | 0.1 | 38 | 19.309 | [12.993, 106.614] | 0.0000 | 0.0000 |
| budget | 0.05 | 39 | 2.496 | [1.797, 15.293] | 0.0000 | 0.0000 |
| budget | 0.02 | 40 | 1.119 | [0.896, 15.142] | 0.0000 | 0.0000 |

## Censoring (fraction of regions where the target is bracketed
on both grids; the rest are conservative bounds, not estimates)

- rmse 0.0003: 8% bracketed
- rmse 0.001: 25% bracketed
- rmse 0.003: 57% bracketed
- rmse 0.01: 80% bracketed
- budget 0.2: 62% bracketed
- budget 0.1: 95% bracketed
- budget 0.05: 98% bracketed
- budget 0.02: 100% bracketed

## Intrinsic gap (gamma=1) and critical overhead

- RMSE 0.01: gap +20%, gamma* = 1.20

## HEALPix regime (latitude)

- rmse:0.01: Spearman rho=0.126 (p=0.493, Holm=0.493), belt 1.139 vs cap 1.219 (p=0.257, Holm=0.257)
    given the spectral slope beta: partial rho=0.134 (p=0.474, n=32); rho(beta, |lat|)=0.000 (p=0.998); rho(beta, ratio)=-0.336
- budget:0.1: Spearman rho=0.445 (p=0.005, Holm=0.010), belt 52.106 vs cap 73.176 (p=0.001, Holm=0.003)
    given the spectral slope beta: partial rho=0.446 (p=0.006, n=38); rho(beta, |lat|)=0.043 (p=0.800); rho(beta, ratio)=0.203
- budget:0.05: Spearman rho=0.528 (p=0.001, Holm=0.002), belt 4.459 vs cap 5.867 (p=0.001, Holm=0.003)
    given the spectral slope beta: partial rho=0.546 (p=0.000, n=39); rho(beta, |lat|)=0.049 (p=0.767); rho(beta, ratio)=0.349
- budget:0.02: Spearman rho=0.516 (p=0.001, Holm=0.002), belt 1.305 vs cap 1.602 (p=0.002, Holm=0.003)
    given the spectral slope beta: partial rho=0.539 (p=0.000, n=40); rho(beta, |lat|)=0.048 (p=0.768); rho(beta, ratio)=0.379

## Codec ranking (HEALPix side)

- 280 operating points, wins sum to the same number
- jpeg2000_lossy: 156 wins / 280 evaluable, mean rank 1.70
- sz3: 50 wins / 280 evaluable, mean rank 2.52
- jpegxl_lossy: 48 wins / 280 evaluable, mean rank 2.97
- dwt_zstd: 13 wins / 280 evaluable, mean rank 3.35
- zfp: 13 wins / 258 evaluable, mean rank 4.41

## HEALPix representation floor (common reference)

- median 1.988e-03, range [1.211e-04, 6.738e-03] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
