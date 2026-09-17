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
| rmse | 0.001 | 39 | 1.227 | [0.308, 5.929] | 0.0078 | 0.0157 |
| rmse | 0.003 | 37 | 0.919 | [0.167, 2.313] | 0.1366 | 0.1366 |
| rmse | 0.01 | 32 | 0.718 | [0.192, 1.107] | 0.0000 | 0.0000 |
| budget | 0.1 | 37 | 4.006 | [2.554, 5.749] | 0.0000 | 0.0000 |
| budget | 0.05 | 39 | 1.333 | [1.019, 2.251] | 0.0000 | 0.0000 |
| budget | 0.02 | 40 | 0.829 | [0.510, 1.418] | 0.0000 | 0.0000 |

## Censoring (fraction of regions where the target is bracketed
on both grids; the rest are conservative bounds, not estimates)

- rmse 0.0001: 100% bracketed
- rmse 0.0003: 100% bracketed
- rmse 0.001: 98% bracketed
- rmse 0.003: 92% bracketed
- rmse 0.01: 80% bracketed
- budget 0.2: 65% bracketed
- budget 0.1: 92% bracketed
- budget 0.05: 98% bracketed
- budget 0.02: 100% bracketed

## Intrinsic gap (gamma=1) and critical overhead

- RMSE 0.0001: gap +106%, gamma* = 2.06
- RMSE 0.0003: gap +94%, gamma* = 1.94
- RMSE 0.001: gap +59%, gamma* = 1.59
- RMSE 0.003: gap +19%, gamma* = 1.19
- RMSE 0.01: gap -7%, gamma* = 0.93

## HEALPix regime (latitude)

- rmse:0.0001: Spearman rho=-0.106 (p=0.515, Holm=1.000), belt 2.092 vs cap 2.056 (p=0.881, Holm=1.000)
    given the spectral slope beta: partial rho=-0.099 (p=0.550, n=40); rho(beta, |lat|)=0.049 (p=0.763); rho(beta, ratio)=-0.189
- rmse:0.0003: Spearman rho=-0.074 (p=0.652, Holm=1.000), belt 1.904 vs cap 1.960 (p=0.654, Holm=1.000)
    given the spectral slope beta: partial rho=-0.068 (p=0.681, n=40); rho(beta, |lat|)=0.049 (p=0.763); rho(beta, ratio)=-0.124
- rmse:0.001: Spearman rho=-0.022 (p=0.892, Holm=1.000), belt 1.651 vs cap 1.586 (p=0.767, Holm=1.000)
    given the spectral slope beta: partial rho=-0.014 (p=0.933, n=39); rho(beta, |lat|)=0.046 (p=0.781); rho(beta, ratio)=-0.190
- rmse:0.003: Spearman rho=0.077 (p=0.649, Holm=1.000), belt 1.130 vs cap 1.226 (p=0.353, Holm=1.000)
    given the spectral slope beta: partial rho=0.107 (p=0.534, n=37); rho(beta, |lat|)=0.108 (p=0.523); rho(beta, ratio)=-0.241
- rmse:0.01: Spearman rho=0.216 (p=0.236, Holm=1.000), belt 0.895 vs cap 1.022 (p=0.031, Holm=0.188)
    given the spectral slope beta: partial rho=0.225 (p=0.223, n=32); rho(beta, |lat|)=0.000 (p=0.998); rho(beta, ratio)=0.292
- budget:0.1: Spearman rho=-0.036 (p=0.832, Holm=1.000), belt 12.365 vs cap 12.662 (p=0.659, Holm=1.000)
    given the spectral slope beta: partial rho=-0.034 (p=0.845, n=37); rho(beta, |lat|)=0.030 (p=0.862); rho(beta, ratio)=-0.078
- budget:0.05: Spearman rho=0.376 (p=0.018, Holm=0.128), belt 2.667 vs cap 2.811 (p=0.006, Holm=0.039)
    given the spectral slope beta: partial rho=0.403 (p=0.012, n=39); rho(beta, |lat|)=0.050 (p=0.762); rho(beta, ratio)=0.488
- budget:0.02: Spearman rho=0.496 (p=0.001, Holm=0.009), belt 0.947 vs cap 1.209 (p=0.001, Holm=0.008)
    given the spectral slope beta: partial rho=0.522 (p=0.001, n=40); rho(beta, |lat|)=0.049 (p=0.763); rho(beta, ratio)=0.408

## Codec ranking (HEALPix side)

- 400 operating points, wins sum to the same number
- jpeg2000_lossy: 196 wins / 400 evaluable, mean rank 1.63
- jpegxl_lossy: 131 wins / 400 evaluable, mean rank 2.35
- zfp: 35 wins / 379 evaluable, mean rank 4.55
- sz3: 28 wins / 400 evaluable, mean rank 2.78
- dwt_zstd: 10 wins / 395 evaluable, mean rank 3.66

## HEALPix representation floor (common reference)

- median 4.038e-03, range [1.340e-04, 1.370e-02] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
