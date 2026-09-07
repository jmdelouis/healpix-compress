# Verdict as a function of the overhead (V4)

Every column below is a function of the per-region critical overhead, so nothing here required re-running a sweep.

| target | n | gamma* median | IQR | HPX wins at g=1 | at 1.3 | at 2.0 | inconclusive range |
|---|---|---|---|---|---|---|---|
| 0.0001 | 40 | 2.07 | [2.00, 2.17] | 2% | 2% | 25% | 2.04-2.12 |
| 0.0003 | 40 | 1.94 | [1.74, 2.07] | 5% | 8% | 68% | 1.76-1.98 |
| 0.001 | 40 | 1.62 | [1.25, 1.92] | 12% | 28% | 88% | 1.39-1.75 |
| 0.003 | 40 | 1.29 | [1.04, 1.63] | 22% | 50% | 95% | 1.15-1.40 |
| 0.01 | 40 | 1.08 | [0.92, 1.31] | 35% | 72% | 98% | 1.04-1.24 |

- median critical overhead spans **1.08 to 2.07** across the reported targets
- targets where HEALPix is still behind for a majority of regions even at an overhead ratio of 2.0: 0.0001

## Measured container cost on each side

- **HPX**: 737 B, 0.0176% of the raw tile — zarr v3 group and array metadata for one level-20 chunk; the cell identifiers are implicit for a contiguous NESTED block, so no index is stored
- **UTM**: 8461 B, 0.5213% of the raw tile — complete single-tile deflate GeoTIFF of an all-zero 637x637 float32 plane, i.e. header, tag directory and a near-empty payload
Container format therefore accounts for at most 0.5213% on either side, a negligible part of any plausible gamma. The illustrative 30% is granule geometry -- swath padding, inter-granule overlap, per-granule metadata -- which is a property of a product layout and not of a format, and is why it is reported as a parameter rather than measured here.
