# Block-size and odd-side bias

5 regions; bytes per pixel, median over regions of the per-region ratio.

| quantity | pair | median | range |
|---|---|---|---|
| lossless_bytes_per_px | utm_odd_over_utm_512 | 0.9982 | [0.9913, 1.0089] |
| lossless_bytes_per_px | utm_640_over_utm_odd | 0.9979 | [0.9972, 0.9996] |
| lossless_bytes_per_px | hpx_odd_over_hpx_512 | 1.0045 | [0.9953, 1.0156] |
| lossless_bytes_per_px | hpx_full_over_hpx_512 | 0.9579 | [0.9407, 1.0133] |
| lossless_bytes_per_px | hpx_full_over_hpx_odd | 0.9624 | [0.9365, 1.0003] |
| bytes_per_px@0.001 | utm_odd_over_utm_512 | 0.9685 | [0.9431, 0.9991] |
| bytes_per_px@0.001 | utm_640_over_utm_odd | 0.9935 | [0.9915, 0.9972] |
| bytes_per_px@0.001 | hpx_odd_over_hpx_512 | 0.9713 | [0.9323, 1.0214] |
| bytes_per_px@0.001 | hpx_full_over_hpx_512 | 0.8927 | [0.8187, 1.0305] |
| bytes_per_px@0.001 | hpx_full_over_hpx_odd | 0.9213 | [0.8781, 1.0198] |
| bytes_per_px@0.01 | utm_odd_over_utm_512 | 0.9332 | [0.8481, 0.9814] |
| bytes_per_px@0.01 | utm_640_over_utm_odd | 0.9943 | [0.9870, 0.9963] |
| bytes_per_px@0.01 | hpx_odd_over_hpx_512 | 0.9278 | [0.8835, 0.9919] |
| bytes_per_px@0.01 | hpx_full_over_hpx_512 | 0.7899 | [0.7192, 0.9053] |
| bytes_per_px@0.01 | hpx_full_over_hpx_odd | 0.8249 | [0.7960, 0.9659] |

Odd-side excess (UTM 637/512 divided by HPX 637/512, lossless): 0.9938, i.e. -0.62% per pixel.
