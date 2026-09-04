# Byte and type audit

- rows: **424**, regions: 2, grids: hpx, utm, scenarios: dn, float32
- files whose re-read bytes differ from what was written: **0**
- points an independent decoder could not read back: **0**
- points where the independent decoder disagrees with the adapter: **0**
- reported size minus payload length: min 0 B, max 4 B (the DWT baseline charges 4 B for a scale it does not write)

## Per codec

| codec | library | version | self-contained | uncharged sidecar | 16-bit floor (RMSE) | lossless claimed | bit-exact vs input | bit-exact vs 16-bit input |
|---|---|---|---|---|---|---|---|---|
| `dwt_zstd` | pywt | 1.8.0 | no | 0 B | -- | no | -- | -- |
| `jpeg2000_lossless` | imagecodecs | 2026.5.10 | no | 8 B | 1.81e-02 | yes | no | yes |
| `jpeg2000_lossy` | PIL | 12.2.0 | no | 8 B | 1.81e-02 | no | -- | -- |
| `jpegxl_lossless` | imagecodecs | 2026.5.10 | no | 8 B | 1.81e-02 | yes | no | yes |
| `jpegxl_lossy` | imagecodecs | 2026.5.10 | no | 8 B | 1.81e-02 | no | -- | -- |
| `sz3` | pysz | unknown | no | 0 B | -- | no | -- | -- |
| `uncompressed` | numpy | 2.4.6 | no | 0 B | -- | yes | yes | yes |
| `zfp` | zfpy | unknown | yes | 0 B | -- | no | -- | -- |

## Per scenario

- **dn**: declared dtype uint16, unit DN (DN = 10000 x rho), conversion error before any codec 2.887e-05 (reflectance units), 0.000% of samples clipped, 0.000% non-finite
  - lossless points: 4/12 bit-exact against the array fed in, 12/12 bit-exact against its 16-bit renormalization
- **float32**: declared dtype float32, unit reflectance, conversion error before any codec 0.000e+00 (reflectance units), 0.000% of samples clipped, 0.000% non-finite
  - lossless points: 4/12 bit-exact against the array fed in, 12/12 bit-exact against its 16-bit renormalization
