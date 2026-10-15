# Paper numbers (common RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 4 targets
- 16-bit requantization floor: HPX 5.79e-06, UTM 6.47e-06 (the image codecs are fed a 16-bit tile, so their `lossless` points are lossless after that step, not bit-exact)

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| rmse | 0.01 | 31 | 1.094 | [0.204, 3.609] | 0.1757 | 0.1757 |
| budget | 0.1 | 38 | 32.517 | [18.391, 162.914] | 0.0000 | 0.0000 |
| budget | 0.05 | 39 | 4.046 | [2.401, 23.338] | 0.0000 | 0.0000 |
| budget | 0.02 | 40 | 1.500 | [1.084, 24.153] | 0.0000 | 0.0000 |

## Censoring (fraction of regions where the target is bracketed
on both grids; the rest are conservative bounds, not estimates)

- rmse 0.001: 18% bracketed
- rmse 0.003: 30% bracketed
- rmse 0.01: 78% bracketed
- budget 0.2: 60% bracketed
- budget 0.1: 95% bracketed
- budget 0.05: 98% bracketed
- budget 0.02: 100% bracketed

## Intrinsic gap (gamma=1) and critical overhead

- RMSE 0.01: gap +42%, gamma* = 1.42

## HEALPix regime (latitude)

- rmse:0.01: Spearman rho=-0.076 (p=0.684, Holm=0.684), belt 1.454 vs cap 1.300 (p=0.858, Holm=0.858)
    given the spectral slope beta: partial rho=-0.100 (p=0.599, n=31); rho(beta, |lat|)=-0.016 (p=0.931); rho(beta, ratio)=-0.531
- budget:0.1: Spearman rho=0.224 (p=0.176, Holm=0.352), belt 100.416 vs cap 126.380 (p=0.068, Holm=0.135)
    given the spectral slope beta: partial rho=0.232 (p=0.167, n=38); rho(beta, |lat|)=0.043 (p=0.800); rho(beta, ratio)=0.478
- budget:0.05: Spearman rho=0.317 (p=0.050, Holm=0.166), belt 7.176 vs cap 9.825 (p=0.023, Holm=0.093)
    given the spectral slope beta: partial rho=0.344 (p=0.034, n=39); rho(beta, |lat|)=0.049 (p=0.767); rho(beta, ratio)=0.535
- budget:0.02: Spearman rho=0.324 (p=0.042, Holm=0.166), belt 1.631 vs cap 2.148 (p=0.035, Holm=0.105)
    given the spectral slope beta: partial rho=0.350 (p=0.029, n=40); rho(beta, |lat|)=0.048 (p=0.768); rho(beta, ratio)=0.522

## Codec ranking (HEALPix side)

- 264 operating points, wins sum to the same number
- jpeg2000_lossy: 151 wins / 264 evaluable, mean rank 1.58
- jpegxl_lossy: 49 wins / 264 evaluable, mean rank 3.05
- sz3: 38 wins / 264 evaluable, mean rank 2.52
- zfp: 16 wins / 242 evaluable, mean rank 4.34
- dwt_zstd: 10 wins / 264 evaluable, mean rank 3.45

## HEALPix representation floor (common reference)

- median 4.038e-03, range [1.340e-04, 1.370e-02] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
