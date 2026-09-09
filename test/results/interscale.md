# Does the link between scales pay?

- regions: 40; rates are ideal entropy bounds, so they flatter this coder against the measured bytes of real codecs

## Subband energy of the exact HEALPix Haar pyramid

| level | coefficients | % of cells | % of variance |
|---|---|---|---|
| 20 | 786432 | 75.00% | 3.09% |
| 19 | 196608 | 18.75% | 7.38% |
| 18 | 49152 | 4.69% | 10.18% |
| 17 | 12288 | 1.17% | 11.20% |
| 16 | 3072 | 0.29% | 10.58% |
| 15 | 768 | 0.07% | 9.96% |
| 14 | 256 | 0.02% | 39.52% |

## Bits per coefficient by context model

| model | b/coef (median) | saving vs order-0 |
|---|---|---|
| no context | 2.745 | 0.0% |
| parent (inter-scale) | 2.466 | 5.9% |
| neighbour (intra-band) | 2.554 | 3.0% |
| parent + neighbour | 2.373 | 7.9% |

- the inter-scale link alone saves 5.9% of the order-0 rate; the intra-band link alone saves 3.0%
- **the parent adds 5.3% on top of the intra-band context**, and that is the number the design lives on, because deployed coders already charge for the intra-band link

- A context-modelled inter-scale coder has 5.3% of rate to trade against Haar's compaction deficit relative to a 9/7 wavelet. The head-to-head below says whether that is enough.

## The zerotree premise

| orientation | P(child = 0) | given parent = 0 | given parent != 0 |
|---|---|---|---|
| d | 0.296 | 0.507 | 0.230 |
| h | 0.239 | 0.444 | 0.184 |
| v | 0.208 | 0.400 | 0.166 |

- a zerotree is worth signalling only if an insignificant parent really predicts insignificant children; here it lifts P(child = 0) by 20.5 points

## By distortion regime

| quartile | RMSE | parent adds over neighbour | bound / codec bytes |
|---|---|---|---|
| tightest | 6.56e-07 | 0.4% | 0.91 |
| tight | 2.39e-04 | 3.0% | 1.13 |
| loose | 1.94e-03 | 9.2% | 1.45 |
| loosest | 1.23e-02 | 11.1% | 1.06 |

- the parent earns its keep where most coefficients are zero, which is the low-rate end; that is the regime a zerotree coder is built for

## Head to head against the measured codec front

| model | entropy bytes / codec bytes (median) |
|---|---|
| no context | 1.20 |
| parent (inter-scale) | 1.11 |
| neighbour (intra-band) | 1.15 |
| parent + neighbour | 1.08 |

- At equal distortion the pyramid's *entropy bound* is 1.08 times the bytes real codecs actually spend on the same tiles. A bound above one is decisive against the design: no implementation can reach its own bound, and this one is already losing.