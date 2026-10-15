# Which carrier level? Level 20 against coarser ones

- regions: 40; the read is the nearest containing cell at every level

| level | cell (m) | cells | raw HPX/UTM | floor | floor % of std | discarded % of variance |
|---|---|---|---|---|---|---|
| 20 | 6.22 | 1048576 | 2.584 | 4.038e-03 | 10.0% | 0.0000% |
| 19 | 12.44 | 262144 | 0.646 | 8.565e-03 | 21.9% | 3.4653% |
| 18 | 24.87 | 65536 | 0.162 | 1.475e-02 | 35.8% | 10.8384% |
| 17 | 49.74 | 16384 | 0.040 | 2.082e-02 | 48.2% | 22.2391% |

- level 20 floor 4.038e-03, level 19 floor 8.565e-03: a factor of 2.12
- aggregating to level 19 discards 3.4653% of the tile's variance. This is *not* the power above the level-19 Nyquist: a $2\times2$ mean is a high-pass whose response $1-\cos^2(k_xd/2)\cos^2(k_yd/2)$ grows like $(kd)^4$, so it feeds mainly on the shoulder below the native Nyquist
- raw storage: level 20 is 2.584 of the equal-area UTM crop, level 19 is 0.646

- **The level-19 floor is 2.12 times level 20's**, so a coarser carrier is not free against the observation. That excess is the information the aggregation removes, not a worse read: the level-19 floor is the level-20 floor and the discarded intra-cell energy in quadrature to within 7% (median measured/quadrature 1.074), and the discarded energy alone accounts for 84% of the excess variance.

## By scene class, level 19 against level 20

| class | floor 20 | floor 19 | ratio | discarded % of variance |
|---|---|---|---|---|
| agriculture | 3.821e-03 | 8.350e-03 | 2.19 | 2.6766% |
| forest | 1.988e-03 | 4.394e-03 | 2.21 | 3.1789% |
| urban | 7.880e-03 | 1.828e-02 | 2.32 | 7.1957% |
| water | 8.052e-04 | 1.948e-03 | 2.42 | 3.7090% |