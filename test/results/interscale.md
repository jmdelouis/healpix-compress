# Does the link between scales pay?

- regions: 40; rates are ideal entropy bounds, so they flatter this coder against the measured bytes of real codecs
- transform haar, tiles clough_tocher

## Subband energy of the exact HEALPix Haar pyramid

| level | coefficients | % of cells | % of variance |
|---|---|---|---|
| 20 | 786432 | 75.00% | 3.52% |
| 19 | 196608 | 18.75% | 8.05% |
| 18 | 49152 | 4.69% | 11.09% |
| 17 | 12288 | 1.17% | 11.91% |
| 16 | 3072 | 0.29% | 10.83% |
| 15 | 768 | 0.07% | 9.68% |
| 14 | 256 | 0.02% | 39.52% |

## Bits per coefficient by context model

| model | b/coef (median) | saving vs order-0 |
|---|---|---|
| no context | 1.473 | 0.0% |
| parent (inter-scale) | 1.230 | 7.0% |
| neighbour (intra-band) | 1.309 | 3.2% |
| parent + neighbour | 1.183 | 9.6% |

- the inter-scale link alone saves 7.0% of the order-0 rate; the intra-band link alone saves 3.2%
- **the parent adds 5.7% on top of the intra-band context**, and that is the number the design lives on, because deployed coders already charge for the intra-band link

- A context-modelled inter-scale coder has 5.7% of rate to trade against Haar's compaction deficit relative to a 9/7 wavelet. The head-to-head below says whether that is enough.

## The zerotree premise

| orientation | P(child = 0) | given parent = 0 | given parent != 0 |
|---|---|---|---|
| d | 0.559 | 0.774 | 0.410 |
| h | 0.479 | 0.730 | 0.333 |
| v | 0.449 | 0.673 | 0.320 |

- a zerotree is worth signalling only if an insignificant parent really predicts insignificant children; here it lifts P(child = 0) by 22.4 points

## By distortion regime

| quartile | RMSE | parent adds over neighbour | bound / codec bytes |
|---|---|---|---|
| tightest | 6.56e-07 | 0.4% | 0.91 |
| tight | 3.30e-04 | 3.9% | 1.17 |
| loose | 3.11e-03 | 12.1% | 1.53 |
| loosest | 2.21e-02 | 9.7% | 2.09 |

- the parent earns its keep where most coefficients are zero, which is the low-rate end; that is the regime a zerotree coder is built for

## Head to head against the measured codec front

| model | entropy bytes / codec bytes (median) |
|---|---|
| no context | 1.46 |
| parent (inter-scale) | 1.32 |
| neighbour (intra-band) | 1.38 |
| parent + neighbour | 1.29 |

- At equal distortion the pyramid's *entropy bound* is 1.29 times the bytes real codecs actually spend on the same tiles. A bound above one is decisive against the design: no implementation can reach its own bound, and this one is already losing.

## Idea 3: what the reader never sees

| reader | floor | transmission ||Re||/||e|| | error energy in ker(R) | cells R can reach | rate saved by a product target |
|---|---|---|---|---|---|
| nearest | 4.038e-03 | 1.005 | 68.8% | 30.8% | 0.1% |
| psf | 1.952e-03 | 0.702 | 18.8% | 81.1% | 24.1% |

- paired floors, psf over nearest: median 0.589, worse on 0/40 regions

- **nearest**: transmission 1.01 with 69% of the error energy in the kernel. The kernel does not make the product error smaller -- it reveals wasted rate. The product depends on 30.8% of the cells, so at high rate a reader-aware coder has a bound of about 3.2x on rate, and it is a bound on *allocation*, not something a coarser step buys.
- **psf**: transmission 0.70. This reader *attenuates* a white coding error rather than ignoring particular cells, so the same product distortion is reached with a coarser step -- a saving that needs no new codec, only the target expressed on the product.

- The gain is a property of the *reader*, not of HEALPix. A gather read has a large kernel exactly because the carrier is oversampled relative to it, which is what Section V-F calls the dominant error; fixing the reader closes the kernel and removes this gain with it. Shaping coding noise into a smoothing reader's stopband is the safer half, and it is safe only for a product consumed through that operator: the companion's point is that the small scales there are real, so a viewer reading the carrier as an image would see the noise put in.
