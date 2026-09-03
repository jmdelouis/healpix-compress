# Paper numbers (common RMSE reference)

- regions: **40**, codecs: zfp, sz3, jpeg2000_lossy, jpegxl_lossy, dwt_zstd
- raw ratio HPX / (UTM x gamma=1.3): **1.988** (gamma=1: 2.584)
- multiple-testing family: 4 targets

## Ratios (oracle codec, per region)

| protocol | target | n | median | range | p | p_Holm |
|---|---|---|---|---|---|---|
| budget | 0.2 | 38 | 24442.904 | [7042.258, 76093.506] | 0.0000 | 0.0000 |
| budget | 0.1 | 39 | 851.753 | [231.072, 3489.555] | 0.0000 | 0.0000 |
| budget | 0.05 | 40 | 91.873 | [25.894, 5771.246] | 0.0000 | 0.0000 |
| budget | 0.02 | 40 | 24.060 | [8.053, 651.001] | 0.0000 | 0.0000 |

## Intrinsic gap (gamma=1) and critical overhead


## HEALPix regime (latitude)

- budget:0.2: Spearman rho=0.179 (p=0.702), belt nan vs cap nan (p=nan)
- budget:0.1: Spearman rho=-0.039 (p=0.813), belt 2130.084 vs cap 2775.396 (p=0.455)
- budget:0.05: Spearman rho=-0.001 (p=0.994), belt 166.851 vs cap 196.334 (p=0.596)
- budget:0.02: Spearman rho=0.025 (p=0.876), belt 29.580 vs cap 33.747 (p=0.559)

## Codec ranking (HEALPix side)

- 192 operating points, wins sum to the same number
- zfp: 123 wins / 123 evaluable, mean rank 1.00
- sz3: 67 wins / 156 evaluable, mean rank 1.57
- jpeg2000_lossy: 2 wins / 160 evaluable, mean rank 2.54

## HEALPix representation floor (common reference)

- median 8.229e-02, range [1.996e-03, 2.127e-01] over 40 regions
- every target RMSE below this floor is unreachable on the HEALPix side by construction and must be reported as such.
