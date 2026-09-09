# Verdict as a function of the overhead (V4)

Every column below is a function of the per-region critical overhead, so nothing here required re-running a sweep.

| target | n | gamma* median | IQR | HPX wins at g=1 | at 1.3 | at 2.0 | inconclusive range |
|---|---|---|---|---|---|---|---|
| 0.0001 | 40 | 2.06 | [1.99, 2.16] | 2% | 2% | 28% | 2.03-2.10 |
| 0.0003 | 40 | 1.94 | [1.74, 2.01] | 5% | 8% | 72% | 1.76-1.97 |
| 0.001 | 40 | 1.60 | [1.25, 1.88] | 12% | 30% | 92% | 1.38-1.73 |
| 0.003 | 40 | 1.25 | [1.04, 1.59] | 22% | 52% | 98% | 1.14-1.37 |
| 0.01 | 40 | 1.06 | [0.92, 1.27] | 35% | 78% | 98% | 1.02-1.20 |

- median critical overhead spans **1.06 to 2.06** across the reported targets
- targets where HEALPix is still behind for a majority of regions even at an overhead ratio of 2.0: 0.0001
