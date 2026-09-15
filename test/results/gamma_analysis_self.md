# Verdict as a function of the overhead (V4)

Every column below is a function of the per-region critical overhead, so nothing here required re-running a sweep.

| target | n | gamma* median | IQR | HPX wins at g=1 | at 1.3 | at 2.0 | inconclusive range |
|---|---|---|---|---|---|---|---|
| 0.0001 | 40 | 2.06 | [1.99, 2.16] | 2% | 2% | 28% | 2.03-2.10 |
| 0.0003 | 40 | 1.94 | [1.74, 2.01] | 5% | 8% | 72% | 1.76-1.97 |
| 0.001 | 39 | 1.59 | [1.22, 1.89] | 13% | 31% | 92% | 1.35-1.73 |
| 0.003 | 37 | 1.19 | [0.97, 1.56] | 27% | 59% | 97% | 1.08-1.33 |
| 0.01 | 32 | 0.93 | [0.89, 1.02] | 69% | 94% | 100% | 0.90-0.99 |

- median critical overhead spans **0.93 to 2.06** across the reported targets
- targets where HEALPix is still behind for a majority of regions even at an overhead ratio of 2.0: 0.0001
