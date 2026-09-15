# Verdict as a function of the overhead (V4)

Every column below is a function of the per-region critical overhead, so nothing here required re-running a sweep.

| target | n | gamma* median | IQR | HPX wins at g=1 | at 1.3 | at 2.0 | inconclusive range |
|---|---|---|---|---|---|---|---|
| 0.0001 | 40 | 2.27 | [2.18, 2.37] | 0% | 2% | 10% | 2.21-2.33 |
| 0.0003 | 40 | 2.18 | [2.00, 2.28] | 2% | 2% | 25% | 2.06-2.21 |
| 0.001 | 39 | 1.82 | [1.48, 2.02] | 8% | 21% | 72% | 1.54-1.92 |
| 0.003 | 37 | 1.61 | [1.18, 1.84] | 11% | 32% | 86% | 1.37-1.71 |
| 0.01 | 33 | 1.53 | [1.25, 1.98] | 12% | 33% | 76% | 1.33-1.72 |

- median critical overhead spans **1.53 to 2.27** across the reported targets
- targets where HEALPix is still behind for a majority of regions even at an overhead ratio of 2.0: 0.0001, 0.0003
