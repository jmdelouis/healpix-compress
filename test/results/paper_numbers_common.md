# Paper numbers (common RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 4 targets
- 16-bit requantization floor: HPX 6.29e-06, UTM 6.23e-06 (the image codecs are fed a 16-bit tile, so their `lossless` points are lossless after that step, not bit-exact)

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| rmse | 0.01 | 39 | 1.182 | [0.402, 3.843] | 0.0002 | 0.0002 |
| budget | 0.1 | 40 | 33.386 | [18.391, 2971.193] | 0.0000 | 0.0000 |
| budget | 0.05 | 40 | 4.081 | [2.401, 467.438] | 0.0000 | 0.0000 |
| budget | 0.02 | 40 | 1.502 | [1.084, 24.153] | 0.0000 | 0.0000 |

## Censoring (fraction of regions where the target is bracketed
on both grids; the rest are conservative bounds, not estimates)

- rmse 0.0003: 8% bracketed
- rmse 0.001: 20% bracketed
- rmse 0.003: 30% bracketed
- rmse 0.01: 80% bracketed
- budget 0.4: 0% bracketed
- budget 0.2: 18% bracketed
- budget 0.1: 98% bracketed
- budget 0.05: 100% bracketed
- budget 0.02: 100% bracketed

## Intrinsic gap (gamma=1) and critical overhead

- RMSE 0.01: gap +54%, gamma* = 1.54

## HEALPix regime (latitude)

- rmse:0.01: Spearman rho=0.034 (p=0.835, Holm=0.835), belt 1.470 vs cap 1.623 (p=0.506, Holm=0.506)
- budget:0.1: Spearman rho=0.251 (p=0.118, Holm=0.237), belt 100.421 vs cap 128.380 (p=0.055, Holm=0.120)
- budget:0.05: Spearman rho=0.319 (p=0.045, Holm=0.134), belt 7.401 vs cap 10.479 (p=0.040, Holm=0.120)
- budget:0.02: Spearman rho=0.359 (p=0.023, Holm=0.092), belt 1.631 vs cap 2.284 (p=0.020, Holm=0.080)

## Codec ranking (HEALPix side)

- 266 operating points, wins sum to the same number
- jpeg2000_lossy: 134 wins / 266 evaluable, mean rank 1.83
- jpegxl_lossy: 61 wins / 266 evaluable, mean rank 2.90
- sz3: 43 wins / 266 evaluable, mean rank 2.46
- zfp: 16 wins / 199 evaluable, mean rank 4.22
- dwt_zstd: 12 wins / 266 evaluable, mean rank 3.39

## HEALPix representation floor (common reference)

- median 4.038e-03, range [1.340e-04, 1.370e-02] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
