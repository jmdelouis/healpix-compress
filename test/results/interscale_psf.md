# Does the link between scales pay?

- regions: 40; rates are ideal entropy bounds, so they flatter this coder against the measured bytes of real codecs
- transform haar, tiles psf

## Subband energy of the exact HEALPix Haar pyramid

| level | coefficients | % of cells | % of variance |
|---|---|---|---|
| 20 | 786432 | 75.00% | 8.17% |
| 19 | 196608 | 18.75% | 11.31% |
| 18 | 49152 | 4.69% | 10.34% |
| 17 | 12288 | 1.17% | 9.97% |
| 16 | 3072 | 0.29% | 9.00% |
| 15 | 768 | 0.07% | 8.33% |
| 14 | 256 | 0.02% | 32.40% |

## Bits per coefficient by context model

| model | b/coef (median) | saving vs order-0 |
|---|---|---|
| no context | 2.413 | 0.0% |
| parent (inter-scale) | 2.254 | 5.2% |
| neighbour (intra-band) | 2.294 | 2.6% |
| parent + neighbour | 2.208 | 7.2% |

- the inter-scale link alone saves 5.2% of the order-0 rate; the intra-band link alone saves 2.6%
- **the parent adds 4.2% on top of the intra-band context**, and that is the number the design lives on, because deployed coders already charge for the intra-band link

- A context-modelled inter-scale coder has 4.2% of rate to trade against Haar's compaction deficit relative to a 9/7 wavelet. The head-to-head below says whether that is enough.

## The zerotree premise

| orientation | P(child = 0) | given parent = 0 | given parent != 0 |
|---|---|---|---|
| d | 0.422 | 0.618 | 0.297 |
| h | 0.332 | 0.571 | 0.251 |
| v | 0.325 | 0.545 | 0.243 |

- a zerotree is worth signalling only if an insignificant parent really predicts insignificant children; here it lifts P(child = 0) by 22.0 points

## By distortion regime

| quartile | RMSE | parent adds over neighbour | bound / codec bytes |
|---|---|---|---|
| tightest | 6.56e-07 | 0.4% | 0.89 |
| tight | 3.86e-04 | 1.9% | 1.01 |
| loose | 3.49e-03 | 8.4% | 1.18 |
| loosest | 2.14e-02 | 9.7% | 0.98 |

- the parent earns its keep where most coefficients are zero, which is the low-rate end; that is the regime a zerotree coder is built for

## Head to head against the measured codec front

| model | entropy bytes / codec bytes (median) |
|---|---|
| no context | 1.05 |
| parent (inter-scale) | 1.02 |
| neighbour (intra-band) | 1.03 |
| parent + neighbour | 0.99 |

- At equal distortion the pyramid's entropy bound is 0.99 of the bytes real codecs spend on the same tiles. That is headroom, and it is necessary but not sufficient -- the bound charges nothing for coding overhead, signalling, or the loss from a per-tile adaptive model.

## Idea 3: what the reader never sees

| reader | floor | transmission ||Re||/||e|| | error energy in ker(R) | cells R can reach | rate saved by a product target |
|---|---|---|---|---|---|
| nearest | 7.835e-03 | 1.007 | 68.7% | 30.8% | 0.1% |
| psf | 4.444e-03 | 0.672 | 19.0% | 81.1% | 20.4% |

- paired floors, psf over nearest: median 0.558, worse on 0/40 regions

- **nearest**: transmission 1.01 with 69% of the error energy in the kernel. The kernel does not make the product error smaller -- it reveals wasted rate. The product depends on 30.8% of the cells, so at high rate a reader-aware coder has a bound of about 3.2x on rate, and it is a bound on *allocation*, not something a coarser step buys.
- **psf**: transmission 0.67. This reader *attenuates* a white coding error rather than ignoring particular cells, so the same product distortion is reached with a coarser step -- a saving that needs no new codec, only the target expressed on the product.

- The gain is a property of the *reader*, not of HEALPix. A gather read has a large kernel exactly because the carrier is oversampled relative to it, which is what Section V-F calls the dominant error; fixing the reader closes the kernel and removes this gain with it. Shaping coding noise into a smoothing reader's stopband is the safer half, and it is safe only for a product consumed through that operator: the companion's point is that the small scales there are real, so a viewer reading the carrier as an image would see the noise put in.
