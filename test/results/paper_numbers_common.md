# Paper numbers (common RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 8 targets
- 16-bit requantization floor: HPX 6.29e-06, UTM 6.23e-06 (the image codecs are fed a 16-bit tile, so their `lossless` points are lossless after that step, not bit-exact)

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| rmse | 0.0003 | 3 | 2.082 | [1.626, 2.403] | nan | nan |
| rmse | 0.001 | 8 | 1.791 | [0.827, 3.865] | 0.0234 | 0.0469 |
| rmse | 0.003 | 12 | 1.255 | [0.239, 4.312] | 0.3394 | 0.3394 |
| rmse | 0.01 | 32 | 1.219 | [0.229, 3.843] | 0.0044 | 0.0132 |
| budget | 0.2 | 38 | 961.698 | [593.039, 5632.986] | 0.0000 | 0.0000 |
| budget | 0.1 | 39 | 32.297 | [16.322, 162.914] | 0.0000 | 0.0000 |
| budget | 0.05 | 40 | 4.081 | [2.401, 214.183] | 0.0000 | 0.0000 |
| budget | 0.02 | 40 | 1.502 | [1.084, 24.153] | 0.0000 | 0.0000 |

## Intrinsic gap (gamma=1) and critical overhead

- RMSE 0.0003: gap +171%, gamma* = 2.71
- RMSE 0.001: gap +133%, gamma* = 2.33
- RMSE 0.003: gap +63%, gamma* = 1.63
- RMSE 0.01: gap +58%, gamma* = 1.58

## HEALPix regime (latitude)

- rmse:0.01: Spearman rho=-0.192 (p=0.293), belt 1.661 vs cap 1.436 (p=0.556)
- budget:0.1: Spearman rho=0.330 (p=0.040), belt 94.248 vs cap 128.380 (p=0.014)
- budget:0.05: Spearman rho=0.319 (p=0.045), belt 7.401 vs cap 10.479 (p=0.040)
- budget:0.02: Spearman rho=0.359 (p=0.023), belt 1.631 vs cap 2.284 (p=0.020)

## Codec ranking (HEALPix side)

- 248 operating points, wins sum to the same number
- jpeg2000_lossy: 125 wins / 189 evaluable, mean rank 1.43
- jpegxl_lossy: 63 wins / 198 evaluable, mean rank 2.45
- zfp: 33 wins / 146 evaluable, mean rank 3.85
- sz3: 19 wins / 210 evaluable, mean rank 2.51
- dwt_zstd: 8 wins / 203 evaluable, mean rank 3.35

## HEALPix representation floor (common reference)

- median 4.038e-03, range [1.340e-04, 1.370e-02] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
