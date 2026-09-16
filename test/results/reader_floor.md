# How much of the floor is the read?

- regions: 40; no codec involved, the uncompressed carrier read on the observation over the intersection mask

| tiles | reader | floor (median) | % of std | cells reached | regions |
|---|---|---|---|---|---|
| clough_tocher | nearest | 4.038e-03 | 10.0% | 30.8% | 40 |
| clough_tocher | bilinear | 4.038e-03 | 10.0% | 30.8% | 40 |
| psf | nearest | 7.835e-03 | 23.6% | 30.8% | 40 |
| psf | bilinear | 7.835e-03 | 23.6% | 30.8% | 40 |
| psf | psf | 1.787e-05 | 0.1% | 81.2% | 40 |

- **clough_tocher tiles**: bilinear/nearest floor = 1.000 (paired median, q10 1.000, q90 1.000; lower on 0/40), so the read accounts for **0%** of the floor's variance on the median region
- **psf tiles**: bilinear/nearest floor = 1.000 (paired median, q10 1.000, q90 1.000; lower on 0/40), so the read accounts for **0%** of the floor's variance on the median region
  - psf/nearest 0.002, psf/bilinear 0.002: the matched operator against the best non-deconvolving read

- **The read is exonerated**: an interpolating reconstructor removes little of the floor, which is therefore the transfer's, and the common-reference verdict stands as scored.