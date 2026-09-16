# How much of the floor is the read?

- regions: 40; no codec involved, the uncompressed carrier read on the observation over the intersection mask

| tiles | reader | floor (median) | % of std | cells reached | regions |
|---|---|---|---|---|---|
| clough_tocher | nearest | 4.038e-03 | 10.0% | 30.8% | 40 |
| clough_tocher | bilinear | 1.940e-03 | 5.9% | 78.4% | 40 |
| psf | nearest | 7.835e-03 | 23.6% | 30.8% | 40 |
| psf | bilinear | 4.288e-03 | 13.8% | 78.4% | 40 |
| psf | psf | 1.787e-05 | 0.1% | 81.2% | 40 |

- **clough_tocher tiles**: bilinear/nearest floor = 0.582 (paired median, q10 0.494, q90 0.750; lower on 40/40), so the read accounts for **66%** of the floor's variance on the median region
- **psf tiles**: bilinear/nearest floor = 0.581 (paired median, q10 0.512, q90 0.611; lower on 40/40), so the read accounts for **66%** of the floor's variance on the median region
  - psf/nearest 0.002, psf/bilinear 0.004: the matched operator against the best non-deconvolving read

- **The floor is mostly the read.** Section V-H's common-reference verdict was scored through the nearest read and must be re-scored with the interpolating one before it is interpreted as a statement about the carrier.