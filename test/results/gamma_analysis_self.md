# Verdict as a function of the overhead (V4)

Every column below is a function of the per-region critical overhead, so nothing here required re-running a sweep.

| target | n | gamma* median | IQR | HPX wins at g=1 | at 1.3 | at 2.0 | inconclusive range |
|---|---|---|---|---|---|---|---|
| 0.0001 | 40 | 2.07 | [1.99, 2.16] | 2% | 2% | 28% | 2.03-2.10 |
| 0.0003 | 40 | 1.95 | [1.74, 2.05] | 5% | 8% | 70% | 1.76-1.97 |
| 0.001 | 39 | 1.58 | [1.22, 1.90] | 13% | 28% | 92% | 1.35-1.72 |
| 0.003 | 37 | 1.23 | [0.97, 1.56] | 27% | 59% | 97% | 1.08-1.34 |
| 0.01 | 32 | 0.93 | [0.89, 1.03] | 69% | 94% | 100% | 0.90-0.99 |

- median critical overhead spans **0.93 to 2.07** across the reported targets
- targets where HEALPix is still behind for a majority of regions even at an overhead ratio of 2.0: 0.0001
