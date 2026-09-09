# Storage surcharge of `linear` over `clough_tocher`

- denominator: latent; regions: 39

| class | e_resample base | e_resample test | ratio |
|---|---|---|---|
| agriculture | 1.5946e-02 | 1.7189e-02 | 1.078 |
| forest | 5.9653e-03 | 8.8406e-03 | 1.482 |
| urban | 7.5302e-02 | 7.6176e-02 | 1.012 |
| water | 4.4331e-03 | 4.8187e-03 | 1.087 |

The two denominators being within a few percent of each other is not an assumption: it is what a recovery fraction near zero implies, and it means the kappa comparison below is close to the absolute one.

## At equal kappa

| target | n | bracketed | cens. base | cens. test | median size ratio | IQR | p |
|---|---|---|---|---|---|---|---|
| 0.1 (+0.5% total error) | 39 | 100% | 0% | 0% | 1.706 | [1.278, 2.219] | 3.6e-12 |
| 0.2 (+2.0% total error) | 39 | 100% | 0% | 0% | 1.862 | [1.308, 2.410] | 2.5e-11 |
| 0.5 (+11.8% total error) | 39 | 82% | 18% | 15% | 1.414 | [1.056, 2.627] | 3.3e-06 |
| 1 (+41.4% total error) *(withheld)* | 39 | 64% | 36% | 33% | 1.011 | [0.969, 1.410] | 1.7e-01 |

- The ratio falls from 1.41 at 0.5 to 1.01 at 1, where only 64% of regions are bracketed on both sides. Censoring inflates the smoother product's size more than the deconvolved one's, so it biases the ratio towards 1: the fall is a property of the sweep's reach, not a turn in the surcharge, and the affected values are floors.

## At equal absolute RMSE

| target | n | bracketed | cens. base | cens. test | median size ratio | IQR | p |
|---|---|---|---|---|---|---|---|
| 0.0001 | 39 | 100% | 0% | 0% | 1.205 | [1.159, 1.257] | 3.6e-12 |
| 0.0003 | 39 | 100% | 0% | 0% | 1.306 | [1.246, 1.406] | 3.6e-12 |
| 0.001 | 39 | 100% | 0% | 0% | 1.603 | [1.330, 1.866] | 7.3e-12 |
| 0.003 | 39 | 92% | 8% | 5% | 1.926 | [1.530, 2.454] | 1.8e-11 |
| 0.01 | 39 | 85% | 15% | 15% | 2.005 | [1.566, 2.577] | 5.4e-09 |

## Against native UTM (gamma = 1.3)

Each HEALPix product's compressed size over the equal-area native crop's, scaled by the assumed overhead, at that product's own target. Composed per region, not as a product of medians. Values below 1 favour HEALPix. Restricted to the regions bracketed on both HEALPix sweeps, so these differ slightly from the paper's own HPX/UTM table, which admits censored bounds over all 40 regions; the comparison that matters here is between the two columns, which are computed on the same regions.

| target | n | `clough_tocher`/UTM | `linear`/UTM | still favours HEALPix |
|---|---|---|---|---|
| kappa 0.1 | 39 | 1.06 | 1.94 | neither |
| kappa 0.2 | 39 | 0.82 | 1.89 | baseline only |
| kappa 0.5 | 32 | 0.81 | 1.83 | baseline only |
| RMSE 0.0001 | 39 | 1.59 | 1.91 | neither |
| RMSE 0.0003 | 39 | 1.50 | 1.93 | neither |
| RMSE 0.001 | 39 | 1.23 | 1.94 | neither |
| RMSE 0.003 | 36 | 0.94 | 1.90 | baseline only |
| RMSE 0.01 | 33 | 0.78 | 1.77 | baseline only |

- **At 0.003, 0.01 the non-deconvolving product is the cheaper representation and the deconvolved one is not**: the surcharge spends the whole margin. This is the storage price of super-resolution, stated as the thing a product designer has to trade away.

- **At kappa = 0.1, the deconvolved product costs 1.71 times the bytes of the non-deconvolving one** for the same relative quality (+0.5% total error each).
- Across the reportable kappa grid the surcharge runs 1.41 to 1.86.