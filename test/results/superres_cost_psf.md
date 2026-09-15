# Storage surcharge of `psf` over `clough_tocher`

- denominator: latent; regions: 40

| class | e_resample base | e_resample test | ratio |
|---|---|---|---|
| agriculture | 1.5946e-02 | 1.5902e-02 | 0.997 |
| forest | 5.9653e-03 | 6.0265e-03 | 1.010 |
| urban | 7.5302e-02 | 7.5831e-02 | 1.007 |
| water | 4.4331e-03 | 4.5914e-03 | 1.036 |

The two denominators being within a few percent of each other is not an assumption: it is what a recovery fraction near zero implies, and it means the kappa comparison below is close to the absolute one.

## At equal kappa

| target | n | bracketed | cens. base | cens. test | median size ratio | IQR | p |
|---|---|---|---|---|---|---|---|
| 0.1 (+0.5% total error) | 40 | 100% | 0% | 0% | 2.005 | [1.679, 2.531] | 1.8e-12 |
| 0.2 (+2.0% total error) | 40 | 100% | 0% | 0% | 2.435 | [2.114, 2.787] | 1.8e-12 |
| 0.5 (+11.8% total error) | 40 | 100% | 0% | 0% | 2.540 | [1.751, 3.366] | 5.3e-08 |
| 1 (+41.4% total error) | 40 | 95% | 5% | 5% | 2.026 | [1.149, 2.538] | 1.7e-07 |

- The ratio falls from 2.54 at 0.5 to 2.03 at 1, where only 95% of regions are bracketed on both sides. Censoring inflates the smoother product's size more than the deconvolved one's, so it biases the ratio towards 1: the fall is a property of the sweep's reach, not a turn in the surcharge, and the affected values are floors.

## At equal absolute RMSE

| target | n | bracketed | cens. base | cens. test | median size ratio | IQR | p |
|---|---|---|---|---|---|---|---|
| 0.0001 | 40 | 100% | 0% | 0% | 1.267 | [1.208, 1.330] | 1.8e-12 |
| 0.0003 | 40 | 100% | 0% | 0% | 1.366 | [1.274, 1.517] | 1.8e-12 |
| 0.001 | 40 | 100% | 0% | 0% | 1.642 | [1.381, 2.262] | 1.8e-12 |
| 0.003 | 40 | 98% | 3% | 3% | 2.078 | [1.661, 2.689] | 5.7e-08 |
| 0.01 | 40 | 95% | 5% | 3% | 2.392 | [1.732, 2.642] | 1.2e-07 |

## Against native UTM (gamma = 1.3)

Each HEALPix product's compressed size over the equal-area native crop's, scaled by the assumed overhead, at that product's own target. Composed per region, not as a product of medians. Values below 1 favour HEALPix. Restricted to the regions bracketed on both HEALPix sweeps, so these differ slightly from the paper's own HPX/UTM table, which admits censored bounds over all 40 regions; the comparison that matters here is between the two columns, which are computed on the same regions.

| target | n | `clough_tocher`/UTM | `psf`/UTM | still favours HEALPix |
|---|---|---|---|---|
| kappa 0.1 | 40 | 1.05 | 2.05 | neither |
| kappa 0.2 | 40 | 0.82 | 2.08 | baseline only |
| kappa 0.5 | 40 | 0.74 | 2.02 | baseline only |
| kappa 1 | 38 | 0.86 | 1.68 | baseline only |
| RMSE 0.0001 | 40 | 1.59 | 2.02 | neither |
| RMSE 0.0003 | 40 | 1.49 | 2.03 | neither |
| RMSE 0.001 | 40 | 1.23 | 2.03 | neither |
| RMSE 0.003 | 39 | 0.92 | 2.00 | baseline only |
| RMSE 0.01 | 38 | 0.73 | 1.71 | baseline only |

- **At 0.003, 0.01 the non-deconvolving product is the cheaper representation and the deconvolved one is not**: the surcharge spends the whole margin. This is the storage price of super-resolution, stated as the thing a product designer has to trade away.

- **At kappa = 0.1, the deconvolved product costs 2.01 times the bytes of the non-deconvolving one** for the same relative quality (+0.5% total error each).
- Across the reportable kappa grid the surcharge runs 2.01 to 2.54.