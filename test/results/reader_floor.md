# How much of the floor is the read?

- regions: 40; no codec involved, the uncompressed carrier read on the observation over the intersection mask

| tiles | reader | floor (median) | % of std | cells reached | regions |
|---|---|---|---|---|---|
| clough_tocher | nearest | 4.038e-03 | 10.0% | 30.8% | 40 |
| clough_tocher | bilinear | 1.869e-03 | 5.0% | 78.4% | 25 |
| psf | nearest | 7.835e-03 | 23.6% | 30.8% | 40 |
| psf | bilinear | 3.495e-03 | 11.2% | 78.4% | 25 |
| psf | psf | 4.444e-03 | 13.2% | 81.1% | 40 |

- **clough_tocher tiles**: bilinear/nearest floor = 0.579 (paired median, q10 0.482, q90 0.729; lower on 25/25), so the read accounts for **66%** of the floor's variance on the median region
- **psf tiles**: bilinear/nearest floor = 0.574 (paired median, q10 0.511, q90 0.605; lower on 25/25), so the read accounts for **67%** of the floor's variance on the median region
  - psf/nearest 0.558, psf/bilinear 0.960: the matched operator against the best non-deconvolving read

- **The floor is mostly the read.** Section V-H's common-reference verdict was scored through the nearest read and must be re-scored with the interpolating one before it is interpreted as a statement about the carrier.