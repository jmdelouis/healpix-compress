# Window radiometry audit

Stored tile against its own window, sampled at the tile's cell centres by nearest pixel. Under the sampling-noise model the regression slope equals the squared correlation whatever the noise level, so the ratio `slope / corr^2` is one for a transfer that does not systematically alter the field. A deconvolving transfer does, so the verdict is taken on the carrier's tiles alone, where the assumption holds: flagged when the ratio departs from one by more than 0.1. The other transfers are measured and reported, never judged; `tile_crop_consistency.py` is the test that covers them.

- regions audited: 40
- passing: 39
- mismatched: 1

Mismatched regions:

- `region__water__lake_geneva`: corr 0.580, slope 0.634, slope/corr^2 1.884, offset -0.00888 (-0.87 of the tile's std), RMS 2.00e-02 (1.97 of the tile's std)

A mismatched window does not affect any result read from the stored tile -- the compression sweeps, the storage accounting, the level choice by aggregation. It does affect every measurement that re-runs a transfer from the window.


`slope / corr^2` over the 40 regions: median 1.0032, range 0.9820--1.8845.

Worst ten by departure from the sampling-noise identity:

| scene | tile std | corr | slope | slope/corr^2 | RMS/std | status |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| region__water__lake_geneva | 0.0102 | 0.5803 | 0.634 | 1.884 | 1.969 | MISMATCH |
| region__water__lake_baikal | 0.0014 | 0.9393 | 0.866 | 0.982 | 0.372 | ok |
| region__urban__sao_paulo | 0.0627 | 0.9645 | 0.940 | 1.011 | 0.271 | ok |
| region__water__lake_victoria | 0.0028 | 0.9102 | 0.820 | 0.990 | 0.460 | ok |
| region__urban__nairobi | 0.0367 | 0.9722 | 0.954 | 1.009 | 0.239 | ok |
| region__urban__sydney | 0.0670 | 0.9738 | 0.955 | 1.007 | 0.232 | ok |
| region__urban__mexico_city | 0.0593 | 0.9686 | 0.945 | 1.007 | 0.255 | ok |
| region__urban__montreal | 0.0511 | 0.9638 | 0.935 | 1.007 | 0.275 | ok |
| region__agriculture__pampas | 0.0528 | 0.9781 | 0.963 | 1.006 | 0.212 | ok |
| region__water__great_salt_lake | 0.0122 | 0.9777 | 0.962 | 1.006 | 0.214 | ok |

## Every transfer's tile, not only the carrier's

A region has one window and one tile per transfer, and they were not built on the same day, so a window can match one tile and not another. Section V-L sets two of these tiles against each other, which is only meaningful where both came from the same acquisition.

| kernel | tiles | median \|ratio-1\| | mismatched | which |
| --- | ---: | ---: | ---: | --- |
| `bicubic` | 39 | 37.798 | not judged (see tile\_crop\_consistency) | -- |
| `clough_tocher` | 40 | 0.004 | 1 | region__water__lake_geneva |
| `linear` | 39 | 0.023 | not judged (see tile\_crop\_consistency) | -- |
| `nearest` | 40 | 0.000 | not judged (see tile\_crop\_consistency) | -- |
| `psf` | 40 | 0.092 | not judged (see tile\_crop\_consistency) | -- |
