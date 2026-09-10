# Does the link between scales pay?

- regions: 40; rates are ideal entropy bounds, so they flatter this coder against the measured bytes of real codecs

## Subband energy of the exact HEALPix Haar pyramid

| level | coefficients | % of cells | % of variance |
|---|---|---|---|
| 20 | 786432 | 75.00% | 0.95% |
| 19 | 196608 | 18.75% | 4.14% |
| 18 | 49152 | 4.69% | 10.12% |
| 17 | 12288 | 1.17% | 10.97% |
| 16 | 3072 | 0.29% | 12.49% |
| 15 | 768 | 0.07% | 12.42% |
| 14 | 256 | 0.02% | 30.12% |

## Bits per coefficient by context model

| model | b/coef (median) | saving vs order-0 |
|---|---|---|
| no context | 2.045 | 0.0% |
| parent (inter-scale) | 1.988 | 3.6% |
| neighbour (intra-band) | 2.002 | 2.7% |
| parent + neighbour | 1.932 | 6.2% |

- the inter-scale link alone saves 3.6% of the order-0 rate; the intra-band link alone saves 2.7%
- **the parent adds 3.4% on top of the intra-band context**, and that is the number the design lives on, because deployed coders already charge for the intra-band link

- A context-modelled inter-scale coder has 3.4% of rate to trade against Haar's compaction deficit relative to a 9/7 wavelet. The head-to-head below says whether that is enough.

## The zerotree premise

| orientation | P(child = 0) | given parent = 0 | given parent != 0 |
|---|---|---|---|
| d | 0.323 | 0.504 | 0.256 |
| h | 0.232 | 0.366 | 0.197 |
| v | 0.254 | 0.425 | 0.210 |

- a zerotree is worth signalling only if an insignificant parent really predicts insignificant children; here it lifts P(child = 0) by 17.1 points

## By distortion regime

| quartile | RMSE | parent adds over neighbour | bound / codec bytes |
|---|---|---|---|
| tightest | 6.56e-07 | 0.3% | 0.87 |
| tight | 2.39e-04 | 2.4% | 0.96 |
| loose | 1.94e-03 | 7.0% | 0.99 |
| loosest | 1.23e-02 | 6.4% | 0.89 |

- the parent earns its keep where most coefficients are zero, which is the low-rate end; that is the regime a zerotree coder is built for

## Head to head against the measured codec front

| model | entropy bytes / codec bytes (median) |
|---|---|
| no context | 1.00 |
| parent (inter-scale) | 0.96 |
| neighbour (intra-band) | 0.97 |
| parent + neighbour | 0.93 |

- At equal distortion the pyramid's entropy bound is 0.93 of the bytes real codecs spend on the same tiles. That is headroom, and it is necessary but not sufficient -- the bound charges nothing for coding overhead, signalling, or the loss from a per-tile adaptive model.

## Idea 3: what the reader never sees

| reader | floor | transmission ||Re||/||e|| | error energy in ker(R) | cells R can reach | rate saved by a product target |
|---|---|---|---|---|---|
| nearest | 4.038e-03 | 1.016 | 68.7% | 30.8% | 0.1% |


- **nearest**: transmission 1.02 with 69% of the error energy in the kernel. The kernel does not make the product error smaller -- it reveals wasted rate. The product depends on 30.8% of the cells, so at high rate a reader-aware coder has a bound of about 3.2x on rate, and it is a bound on *allocation*, not something a coarser step buys.

- The gain is a property of the *reader*, not of HEALPix. A gather read has a large kernel exactly because the carrier is oversampled relative to it, which is what Section V-F calls the dominant error; fixing the reader closes the kernel and removes this gain with it. Shaping coding noise into a smoothing reader's stopband is the safer half, and it is safe only for a product consumed through that operator: the companion's point is that the small scales there are real, so a viewer reading the carrier as an image would see the noise put in.
