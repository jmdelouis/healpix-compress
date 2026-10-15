# Block-size and odd-side bias

40 regions; bytes per pixel, median over regions of the per-region ratio.

| quantity | pair | median | range |
|---|---|---|---|
| lossless_bytes_per_px | utm_odd_over_utm_512 | 1.0015 | [0.9724, 2.5343] |
| lossless_bytes_per_px | utm_640_over_utm_odd | 0.9978 | [0.9962, 1.0071] |
| lossless_bytes_per_px | hpx_odd_over_hpx_512 | 1.0010 | [0.8687, 1.0318] |
| lossless_bytes_per_px | hpx_full_over_hpx_512 | 0.9638 | [0.7536, 2.5886] |
| lossless_bytes_per_px | hpx_full_over_hpx_odd | 0.9795 | [0.8383, 2.5699] |
| bytes_per_px@0.001 | utm_odd_over_utm_512 | 0.9961 | [0.7985, 4.1834] |
| bytes_per_px@0.001 | utm_640_over_utm_odd | 0.9933 | [0.9866, 1.0009] |
| bytes_per_px@0.001 | hpx_odd_over_hpx_512 | 0.9975 | [0.7851, 3.9922] |
| bytes_per_px@0.001 | hpx_full_over_hpx_512 | 0.9741 | [0.5366, 43.0912] |
| bytes_per_px@0.001 | hpx_full_over_hpx_odd | 0.9728 | [0.6131, 46.4520] |
| bytes_per_px@0.01 | utm_odd_over_utm_512 | 0.9683 | [0.5942, 5.9338] |
| bytes_per_px@0.01 | utm_640_over_utm_odd | 0.9899 | [0.9782, 1.0080] |
| bytes_per_px@0.01 | hpx_odd_over_hpx_512 | 0.9567 | [0.6855, 3.9922] |
| bytes_per_px@0.01 | hpx_full_over_hpx_512 | 0.8412 | [0.3281, 1.5258] |
| bytes_per_px@0.01 | hpx_full_over_hpx_odd | 0.8579 | [0.2087, 1.2564] |

Odd-side excess (UTM 637/512 divided by HPX 637/512, lossless): 1.0004, i.e. +0.04% per pixel.
