# Paper numbers (self RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 8 targets
- 16-bit requantization floor: HPX 5.79e-06, UTM 6.47e-06 (the image codecs are fed a 16-bit tile, so their `lossless` points are lossless after that step, not bit-exact)

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| rmse | 0.0001 | 40 | 1.594 | [0.705, 1.740] | 0.0000 | 0.0000 |
| rmse | 0.0003 | 40 | 1.503 | [0.537, 1.704] | 0.0000 | 0.0000 |
| rmse | 0.001 | 39 | 1.214 | [0.308, 1.616] | 0.0082 | 0.0164 |
| rmse | 0.003 | 37 | 0.943 | [0.167, 1.779] | 0.1866 | 0.1866 |
| rmse | 0.01 | 32 | 0.718 | [0.192, 1.107] | 0.0000 | 0.0000 |
| budget | 0.1 | 38 | 4.024 | [2.554, 7.132] | 0.0000 | 0.0000 |
| budget | 0.05 | 39 | 1.333 | [1.019, 2.251] | 0.0000 | 0.0000 |
| budget | 0.02 | 40 | 0.824 | [0.510, 1.418] | 0.0000 | 0.0000 |

## Censoring (fraction of regions where the target is bracketed
on both grids; the rest are conservative bounds, not estimates)

- rmse 0.0001: 100% bracketed
- rmse 0.0003: 100% bracketed
- rmse 0.001: 98% bracketed
- rmse 0.003: 92% bracketed
- rmse 0.01: 80% bracketed
- budget 0.2: 62% bracketed
- budget 0.1: 95% bracketed
- budget 0.05: 98% bracketed
- budget 0.02: 100% bracketed

## Intrinsic gap (gamma=1) and critical overhead

- RMSE 0.0001: gap +107%, gamma* = 2.07
- RMSE 0.0003: gap +95%, gamma* = 1.95
- RMSE 0.001: gap +58%, gamma* = 1.58
- RMSE 0.003: gap +23%, gamma* = 1.23
- RMSE 0.01: gap -7%, gamma* = 0.93

## HEALPix regime (latitude)

- rmse:0.0001: Spearman rho=-0.055 (p=0.735, Holm=1.000), belt 2.092 vs cap 2.062 (p=0.796, Holm=1.000)
    given the spectral slope beta: partial rho=-0.049 (p=0.767, n=40); rho(beta, |lat|)=0.048 (p=0.768); rho(beta, ratio)=-0.141
- rmse:0.0003: Spearman rho=-0.034 (p=0.833, Holm=1.000), belt 1.904 vs cap 1.967 (p=0.596, Holm=1.000)
    given the spectral slope beta: partial rho=-0.030 (p=0.854, n=40); rho(beta, |lat|)=0.048 (p=0.768); rho(beta, ratio)=-0.084
- rmse:0.001: Spearman rho=0.007 (p=0.965, Holm=1.000), belt 1.651 vs cap 1.572 (p=0.725, Holm=1.000)
    given the spectral slope beta: partial rho=0.014 (p=0.932, n=39); rho(beta, |lat|)=0.045 (p=0.786); rho(beta, ratio)=-0.151
- rmse:0.003: Spearman rho=0.125 (p=0.463, Holm=1.000), belt 1.130 vs cap 1.229 (p=0.266, Holm=1.000)
    given the spectral slope beta: partial rho=0.151 (p=0.380, n=37); rho(beta, |lat|)=0.107 (p=0.528); rho(beta, ratio)=-0.205
- rmse:0.01: Spearman rho=0.212 (p=0.245, Holm=1.000), belt 0.895 vs cap 1.031 (p=0.031, Holm=0.188)
    given the spectral slope beta: partial rho=0.221 (p=0.232, n=32); rho(beta, |lat|)=0.000 (p=0.998); rho(beta, ratio)=0.295
- budget:0.1: Spearman rho=-0.017 (p=0.919, Holm=1.000), belt 12.365 vs cap 12.923 (p=0.589, Holm=1.000)
    given the spectral slope beta: partial rho=-0.014 (p=0.935, n=38); rho(beta, |lat|)=0.043 (p=0.800); rho(beta, ratio)=-0.078
- budget:0.05: Spearman rho=0.363 (p=0.023, Holm=0.163), belt 2.667 vs cap 2.811 (p=0.007, Holm=0.050)
    given the spectral slope beta: partial rho=0.387 (p=0.016, n=39); rho(beta, |lat|)=0.049 (p=0.767); rho(beta, ratio)=0.481
- budget:0.02: Spearman rho=0.490 (p=0.001, Holm=0.011), belt 0.947 vs cap 1.194 (p=0.001, Holm=0.008)
    given the spectral slope beta: partial rho=0.515 (p=0.001, n=40); rho(beta, |lat|)=0.048 (p=0.768); rho(beta, ratio)=0.406

## Codec ranking (HEALPix side)

- 400 operating points, wins sum to the same number
- jpeg2000_lossy: 199 wins / 400 evaluable, mean rank 1.62
- jpegxl_lossy: 131 wins / 400 evaluable, mean rank 2.33
- zfp: 35 wins / 378 evaluable, mean rank 4.54
- sz3: 26 wins / 400 evaluable, mean rank 2.79
- dwt_zstd: 9 wins / 395 evaluable, mean rank 3.68

## HEALPix representation floor (common reference)

- median 4.038e-03, range [1.340e-04, 1.370e-02] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
