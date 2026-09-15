# Window radiometry audit

Stored tile against its own window, sampled at the tile's cell centres by nearest pixel. A region passes when the correlation exceeds 0.95, the slope is within 0.15 of one and the offset is under 0.5 of the tile's standard deviation.

- regions audited: 40
- passing: 37
- mismatched: 3

Mismatched regions:

- `region__water__lake_baikal`: corr 0.939, slope 0.866, offset +0.00182 (+1.28 of the tile's std), RMS 5.27e-04 (0.37 of the tile's std)
- `region__water__lake_geneva`: corr 0.580, slope 0.634, offset -0.00888 (-0.87 of the tile's std), RMS 2.00e-02 (1.97 of the tile's std)
- `region__water__lake_victoria`: corr 0.910, slope 0.820, offset +0.00602 (+2.14 of the tile's std), RMS 1.29e-03 (0.46 of the tile's std)

A mismatched window does not affect any result read from the stored tile -- the compression sweeps, the storage accounting, the level choice by aggregation. It does affect every measurement that re-runs a transfer from the window.


Worst ten by RMS over the tile's standard deviation:

| scene | corr | slope | offset/std | RMS/std | status |
| --- | ---: | ---: | ---: | ---: | --- |
| region__water__lake_geneva | 0.5803 | 0.634 | -0.874 | 1.969 | MISMATCH |
| region__water__lake_victoria | 0.9102 | 0.820 | +2.141 | 0.460 | MISMATCH |
| region__water__lake_baikal | 0.9393 | 0.866 | +1.283 | 0.372 | MISMATCH |
| region__urban__tokyo | 0.9504 | 0.907 | +0.185 | 0.325 | ok |
| region__urban__montreal | 0.9638 | 0.935 | +0.116 | 0.275 | ok |
| region__urban__sao_paulo | 0.9645 | 0.940 | +0.096 | 0.271 | ok |
| region__forest__appalachian | 0.9666 | 0.937 | +0.203 | 0.265 | ok |
| region__water__lake_tahoe | 0.9674 | 0.934 | +0.140 | 0.262 | ok |
| region__urban__paris | 0.9669 | 0.940 | +0.109 | 0.262 | ok |
| region__urban__mexico_city | 0.9686 | 0.945 | +0.111 | 0.255 | ok |
