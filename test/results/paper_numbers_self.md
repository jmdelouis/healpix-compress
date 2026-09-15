# Paper numbers (self RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 8 targets
- 16-bit requantization floor: HPX 6.75e-06, UTM 6.23e-06 (the image codecs are fed a 16-bit tile, so their `lossless` points are lossless after that step, not bit-exact)

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| rmse | 0.0001 | 40 | 1.747 | [0.939, 2.022] | 0.0000 | 0.0000 |
| rmse | 0.0003 | 40 | 1.677 | [0.537, 2.134] | 0.0000 | 0.0000 |
| rmse | 0.001 | 39 | 1.404 | [0.444, 3.723] | 0.0002 | 0.0009 |
| rmse | 0.003 | 37 | 1.241 | [0.413, 5.219] | 0.0076 | 0.0151 |
| rmse | 0.01 | 33 | 1.174 | [0.208, 2.791] | 0.0230 | 0.0230 |
| budget | 0.1 | 37 | 5.681 | [2.555, 9.744] | 0.0000 | 0.0000 |
| budget | 0.05 | 39 | 1.926 | [1.218, 3.434] | 0.0000 | 0.0000 |
| budget | 0.02 | 40 | 1.213 | [0.714, 2.611] | 0.0003 | 0.0010 |

## Censoring (fraction of regions where the target is bracketed
on both grids; the rest are conservative bounds, not estimates)

- rmse 0.0001: 100% bracketed
- rmse 0.0003: 100% bracketed
- rmse 0.001: 98% bracketed
- rmse 0.003: 92% bracketed
- rmse 0.01: 82% bracketed
- budget 0.2: 65% bracketed
- budget 0.1: 92% bracketed
- budget 0.05: 98% bracketed
- budget 0.02: 100% bracketed

## Intrinsic gap (gamma=1) and critical overhead

- RMSE 0.0001: gap +127%, gamma* = 2.27
- RMSE 0.0003: gap +118%, gamma* = 2.18
- RMSE 0.001: gap +82%, gamma* = 1.82
- RMSE 0.003: gap +61%, gamma* = 1.61
- RMSE 0.01: gap +53%, gamma* = 1.53

## HEALPix regime (latitude)

- rmse:0.0001: Spearman rho=-0.049 (p=0.762, Holm=1.000), belt 2.267 vs cap 2.275 (p=0.860, Holm=1.000)
- rmse:0.0003: Spearman rho=-0.113 (p=0.488, Holm=1.000), belt 2.180 vs cap 2.176 (p=0.615, Holm=1.000)
- rmse:0.001: Spearman rho=-0.072 (p=0.664, Holm=1.000), belt 1.905 vs cap 1.769 (p=0.745, Holm=1.000)
- rmse:0.003: Spearman rho=-0.027 (p=0.873, Holm=1.000), belt 1.624 vs cap 1.531 (p=0.964, Holm=1.000)
- rmse:0.01: Spearman rho=0.259 (p=0.146, Holm=1.000), belt 1.412 vs cap 1.560 (p=0.124, Holm=0.871)
- budget:0.1: Spearman rho=0.123 (p=0.470, Holm=1.000), belt 15.150 vs cap 19.089 (p=0.131, Holm=0.871)
- budget:0.05: Spearman rho=0.137 (p=0.405, Holm=1.000), belt 3.794 vs cap 4.343 (p=0.242, Holm=1.000)
- budget:0.02: Spearman rho=0.171 (p=0.291, Holm=1.000), belt 1.523 vs cap 1.737 (p=0.080, Holm=0.636)

## Codec ranking (HEALPix side)

- 400 operating points, wins sum to the same number
- jpeg2000_lossy: 152 wins / 400 evaluable, mean rank 2.04
- jpegxl_lossy: 97 wins / 400 evaluable, mean rank 2.60
- sz3: 87 wins / 400 evaluable, mean rank 2.47
- dwt_zstd: 33 wins / 395 evaluable, mean rank 3.44
- zfp: 31 wins / 379 evaluable, mean rank 4.40

## HEALPix representation floor (common reference)

- median 4.038e-03, range [1.340e-04, 1.370e-02] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
