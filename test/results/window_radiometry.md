# Window radiometry audit

Stored tile against its own window, sampled at the tile's cell centres by nearest pixel. Under the sampling-noise model the regression slope equals the squared correlation whatever the noise level, so a region passes when `slope / corr^2` is within 0.1 of one (and the correlation clears 0.8 as a backstop).

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

| kernel | tiles audited | mismatched | which |
| --- | ---: | ---: | --- |
| `bicubic` | 39 | 39 | region__agriculture__beauce, region__agriculture__canterbury, region__agriculture__iowa, region__agriculture__mato_grosso, region__agriculture__nile_delta, region__agriculture__pampas, region__agriculture__po_valley, region__agriculture__punjab, region__agriculture__ukraine, region__agriculture__western_australia, region__forest__amazon, region__forest__appalachian, region__forest__black_forest, region__forest__borneo, region__forest__canadian_rockies, region__forest__congo, region__forest__finland_boreal, region__forest__new_zealand, region__forest__patagonia, region__forest__tasmania, region__urban__cape_town, region__urban__mexico_city, region__urban__montreal, region__urban__nairobi, region__urban__paris, region__urban__reykjavik, region__urban__sao_paulo, region__urban__singapore, region__urban__sydney, region__urban__tokyo, region__water__caspian_sea, region__water__great_salt_lake, region__water__lake_baikal, region__water__lake_constance, region__water__lake_geneva, region__water__lake_michigan, region__water__lake_tahoe, region__water__lake_titicaca, region__water__lake_victoria |
| `clough_tocher` | 40 | 1 | region__water__lake_geneva |
| `linear` | 39 | 6 | region__agriculture__mato_grosso, region__forest__congo, region__water__great_salt_lake, region__water__lake_baikal, region__water__lake_titicaca, region__water__lake_victoria |
| `nearest` | 40 | 0 | -- |
| `psf` | 40 | 19 | region__agriculture__iowa, region__agriculture__nile_delta, region__agriculture__pampas, region__forest__amazon, region__forest__appalachian, region__forest__tasmania, region__urban__mexico_city, region__urban__montreal, region__urban__nairobi, region__urban__paris, region__urban__reykjavik, region__urban__sao_paulo, region__urban__sydney, region__urban__tokyo, region__water__great_salt_lake, region__water__lake_baikal, region__water__lake_tahoe, region__water__lake_titicaca, region__water__lake_victoria |
