# Which carrier level? Level 20 against coarser ones

- regions: 40; the read is the nearest containing cell at every level

| level | cell (m) | cells | raw HPX/UTM | floor | floor % of std | discarded % of variance |
|---|---|---|---|---|---|---|
| 20 | 6.22 | 1048576 | 2.584 | 4.038e-03 | 9.9% | 0.0000% |
| 19 | 12.44 | 262144 | 0.646 | 8.565e-03 | 21.7% | 3.0685% |
| 18 | 24.87 | 65536 | 0.162 | 1.427e-02 | 34.9% | 10.6063% |
| 17 | 49.74 | 16384 | 0.040 | 2.004e-02 | 48.1% | 20.4015% |

- level 20 floor 4.038e-03, level 19 floor 8.565e-03: a factor of 2.12
- aggregating to level 19 discards 3.0685% of the tile's variance, against the 0.21% the imagery carries above half the native Nyquist
- raw storage: level 20 is 2.584 of the equal-area UTM crop, level 19 is 0.646

- **The level-19 floor is 2.12 times level 20's**, so a coarser carrier is not free against the observation. How much of that is information and how much is the coarser read is the next question: the aggregation itself discards only 3.0685% of the variance, so most of the gap is the read, not the level.

## By scene class, level 19 against level 20

| class | floor 20 | floor 19 | ratio | discarded % of variance |
|---|---|---|---|---|
| agriculture | 3.821e-03 | 8.350e-03 | 2.19 | 2.6766% |
| forest | 1.857e-03 | 4.047e-03 | 2.18 | 2.1696% |
| urban | 7.880e-03 | 1.828e-02 | 2.32 | 7.1957% |
| water | 8.052e-04 | 1.948e-03 | 2.42 | 3.7065% |