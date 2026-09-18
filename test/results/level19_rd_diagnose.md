# Level-19 diagnostic: aggregate against direct

| target | n direct | n agg | direct/UTM | agg/UTM | agg/direct (paired, same regions) | agg < direct | compact. UTM | compact. direct | compact. agg |
|---|---|---|---|---|---|---|---|---|---|
| 0.0001 | 39 | 40 | 0.666 | 0.636 | 0.951 (n=39) | 1.00 | 0.187 | 0.191 | 0.182 |
| 0.0003 | 39 | 40 | 0.673 | 0.634 | 0.926 (n=39) | 1.00 | 0.132 | 0.139 | 0.128 |
| 0.001 | 38 | 40 | 0.695 | 0.634 | 0.894 (n=38) | 1.00 | 0.080 | 0.088 | 0.079 |
| 0.003 | 36 | 38 | 0.721 | 0.641 | 0.877 (n=36) | 1.00 | 0.035 | 0.043 | 0.037 |
| 0.01 | 34 | 34 | 0.783 | 0.642 | 0.829 (n=34) | 0.97 | 0.008 | 0.011 | 0.009 |

## Oracle codec per carrier (share of regions)

- aggregate, target 0.0001: jpeg2000_lossy 85%, sz3 10%, jpegxl_lossy 5%
- aggregate, target 0.0003: jpeg2000_lossy 92%, jpegxl_lossy 8%
- aggregate, target 0.001: jpeg2000_lossy 90%, jpegxl_lossy 8%, sz3 2%
- aggregate, target 0.003: jpeg2000_lossy 85%, dwt_zstd 10%, jpegxl_lossy 2%, sz3 2%
- aggregate, target 0.01: jpeg2000_lossy 57%, jpegxl_lossy 22%, dwt_zstd 18%, sz3 2%
- direct, target 0.0001: jpeg2000_lossy 82%, sz3 13%, jpegxl_lossy 5%
- direct, target 0.0003: jpeg2000_lossy 95%, jpegxl_lossy 5%
- direct, target 0.001: jpeg2000_lossy 82%, sz3 8%, jpegxl_lossy 8%, dwt_zstd 3%
- direct, target 0.003: jpeg2000_lossy 85%, dwt_zstd 8%, sz3 5%, jpegxl_lossy 3%
- direct, target 0.01: jpeg2000_lossy 62%, jpegxl_lossy 21%, dwt_zstd 15%, sz3 3%

## Reference check
- aggregate sweep rows: 2720, raw_bytes per row [np.int64(1048576)] (expected 1048576 for a 512^2 float32 raster: the RMSE column is against that raster)
- compactability columns are compressed / raw bytes of each raster; the aggregate is a 2x2 mean of the Clough--Tocher tile, hence smoother than the directly evaluated carrier, and a flat agg/UTM column with a falling agg/direct column means the aggregate gains compressibility with the target at the same rate as the native crop