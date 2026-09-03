# Paper numbers (common RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 4 targets
- 16-bit requantization floor: HPX 6.29e-06, UTM 6.23e-06 (the image codecs are fed a 16-bit tile, so their `lossless` points are lossless after that step, not bit-exact)

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| budget | 0.2 | 38 | 25287.890 | [7079.545, 76297.534] | 0.0000 | 0.0000 |
| budget | 0.1 | 39 | 850.727 | [230.622, 3498.911] | 0.0000 | 0.0000 |
| budget | 0.05 | 40 | 94.587 | [25.844, 3695.193] | 0.0000 | 0.0000 |
| budget | 0.02 | 40 | 24.227 | [8.038, 416.820] | 0.0000 | 0.0000 |

## Intrinsic gap (gamma=1) and critical overhead


## HEALPix regime (latitude)

- budget:0.2: Spearman rho=0.179 (p=0.702), belt nan vs cap nan (p=nan)
- budget:0.1: Spearman rho=-0.063 (p=0.703), belt 2157.918 vs cap 2846.540 (p=0.526)
- budget:0.05: Spearman rho=-0.002 (p=0.991), belt 167.708 vs cap 202.010 (p=0.559)
- budget:0.02: Spearman rho=0.016 (p=0.920), belt 29.422 vs cap 33.718 (p=0.541)

## Codec ranking (HEALPix side)

- 192 operating points, wins sum to the same number
- zfp: 123 wins / 123 evaluable, mean rank 1.00
- sz3: 67 wins / 156 evaluable, mean rank 1.57
- jpeg2000_lossy: 2 wins / 160 evaluable, mean rank 2.54

## HEALPix representation floor (common reference)

- median 8.116e-02, range [1.992e-03, 2.093e-01] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
