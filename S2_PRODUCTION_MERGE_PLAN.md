# Sentinel-2 production on merged packages: review and merge order

Status on 2026-10-15. Goal: run the Sentinel-2 L2A to HEALPix production with released versions of the GRID4EARTH packages only. The ten-product production test of healpix-compress (paper appendix "Large-Scale Trial on Whole Sentinel-2 Products") ran with unmerged branches and with workarounds in its driver script (`test/s2_full_products.py`). Each workaround now has a pull request. This document lists them in the order they should be reviewed and merged, and why.

## 1. Already done

| Package | What | Why it matters |
|---|---|---|
| healpix-resample | Every open PR merged; **v2026.10.0** on PyPI | Contains the fixes the production depended on: no all-NaN band from one NaN input pixel (#50, v2026.07.0), no crash on a chunk without support (#71), no memory blow-up in the PSF fallback |
| healpix-compress | #30, #33, #31 merged; **v0.2.0** released on prefix.dev/grid4earth | The `hcc-sentinel2` and `hcc-sentinel2-pyramid` serializers. Streams of the production decode bit-identically with this release |

## 2. healpix-convert (reviewer: @benbovy)

Merge in this order, then cut a release.

| Order | PR | Author | What | Why this place |
|---|---|---|---|---|
| 1 | #24 | tinaok | Top-left-corner transform from x/y coordinates (half-pixel shift) | **Correctness.** Without it every EOPF Sentinel-2 L2A product is placed half a pixel off. The production ran with it. Touches `core/utils.py` only. |
| 2 | #25 | tinaok | Leave a chunk empty when its input points cannot fill it | **Robustness.** Replaces the driver's `_skip_chunks_without_support`: a tile-edge chunk no longer stops a whole conversion. |
| 3 | #40 | jmdelouis | Write `object` string arrays as variable-length strings (closes #33) | **Blocking for L2A.** The `conditions/geometry` group of an EOPF product cannot be written with zarr 3 without it. Same behaviour as the production workaround. |
| 4 | #41 | jmdelouis | Warn when a variable with valid input comes back all-NaN (closes #37) | **Visibility.** Silent chunk loss is how 4 of 9 products lost 21 to 106 chunks per band with healpix-resample v2026.06.0. Logging only, no change to the output. |
| 5 | #42 | jmdelouis | Delta-code the cell ids before zstd (closes #32) | **Size.** 938 bytes instead of 267 kB per level-20 chunk of cell ids; the mask, footprint and classification groups were almost entirely cell ids. **Needs a decision:** `numcodecs.delta` is not a core zarr v3 codec (readers need numcodecs). The production used it. |
| 6 | #45 | jmdelouis | NaN fill value for floating resampled variables (closes #44) | **Correctness.** With zarr's default fill value 0.0, the cells without data are written as zero reflectances, which `hcc-sentinel2-pyramid` codes as values: on the Paris product, ~131 000 cells of one level-20 chunk at the edge of the data read back as near-zero data (a black band in a viewer). Changes the `fill_value` metadata of every floating output. After #25: the author adapts ten of #25's tests, which compare with `==` (never true for NaN). |
| 7 | #27 | jmdelouis | Check the group codecs when the settings are created; document the healpix-compress serializers | Turns a missing codec package into an error at settings time instead of on a worker mid-conversion. Useful before legacy-converters#100 puts these codecs in the presets. |
| 8 | #31 | jmdelouis | Resample the variables of a chunk in one batched call | **Speed** (one PSF operator per chunk instead of one per variable). Raises the `healpix-resample` floor to `>=2026.7.0`, below which one NaN pixel turns the whole batch NaN. Conflicts with #41 (same loop): the author rebases it once #41 is in. |
| 9 | #39 | jmdelouis | Write the output arrays of a chunk concurrently (closes #38) | **Speed** (with the serializers' worker threads, 830 to 440 ms per 10 m chunk). Built on the same function as #31: the author rebases it once #31 is in. |
| 10 | #29 | tinaok | Multi-resolution (pyramid) stage, for Gridlook | Not needed to convert, needed to view a product at every zoom. Tested on the whole Paris product transcoded to the default codecs: 2 352 comparisons with an independent block mean, no difference above 8e-15, 2 h 07. On its own it refuses the production format (`hcc-sentinel2-pyramid`). |
| 11 | #43 | jmdelouis | Parent-aligned dense chunks for the Sentinel-2 serializers (stacked on #29) | Makes #29 work on the production format: dense chunks per parent cell, native levels copied as stored, derived levels coded with `hcc-sentinel2-pyramid`. Tested on the Paris production product: 750/750 native chunks byte-identical, 2 352 derived comparisons within half the quantisation step, 1.39 GB for the whole pyramid against 1.17 GB for the product (18.4 GB with #29 on the default codecs), 3 h 48. |

Items 1 to 7 and 10, 11 merge cleanly on top of each other (checked with local trial merges); with #25 and #45 together, ten of #25's tests need the NaN-aware comparison above. Items 8 and 9 are last among the converter changes because they rewrite the per-chunk loop the others touch, and their author rebases them.

Before the review by @benbovy, the author's PRs (#27, #31, #39, #40, #41, #42) got a self-review against the GRID4EARTH rubric (review-rules#3), as tinaok did on #24 and #25 and as she asked on #27.

Not needed for the Sentinel-2 production (review independently): #4 (`cell_methods`, has conflicts).

## 3. legacy-converters

| Order | PR | What | Why |
|---|---|---|---|
| 1 | #100 | Sentinel-2 L2A settings: reflectance under `hcc-sentinel2-pyramid`, atmosphere, scene classification and probabilities under lossless `hcc-sentinel2` (closes #96) | Puts the production's codecs in the preset; today the driver overrides the settings (`make_settings`). The codecs equal the production's on all 19 groups. **Needs a decision:** whether `pixi.toml` gets the private `grid4earth` channel and `healpix-compress`, or the package stays a requirement of the environment that converts. |

Then a release, depending on the healpix-convert release of section 2.

Deliberately left out: #98 (bounded PSF fit `clip=(0, None)`), which would change the reflectance relative to the validated production; #99 (`units`, `dtype`, `fill_value` attributes of the reflectance), metadata only, which needs healpix-convert#36 first.

## 4. g4e-desp-argo

| PR | What | Change needed |
|---|---|---|
| #29 | Pin healpix-resample in the `s2msi` images (closes #28) | Today it pins commit `612c6bd`. Once sections 2 and 3 are released, pin the released versions instead: healpix-resample `>=2026.10.0`, healpix-convert and legacy-converters at their new releases, healpix-compress `>=0.2.0` with `imagecodecs`, then rebuild the images. |

## 5. Validation

Convert one tile (T31UDQ, 2026-08-12) with the merged and released stack only, no driver workaround, and compare with the production store:

- the serializers are bit-identical, so any difference comes from the transfer;
- expected difference: healpix-resample v2026.10.0 includes #76, which changes the cells at the tile edges whose neighbourhood is partial; the production used commit `4786405`;
- check sizes, the error on the native pixels (no-data pixels left out) and that no chunk is reported lost (#41).

## Driver workarounds and where they go

| Workaround in `s2_full_products.py` | Replaced by |
|---|---|
| `_string_coordinates_as_zarr_strings` (object strings) | healpix-convert#40 |
| same function, delta filter on `cell_ids` | healpix-convert#42 |
| `_skip_chunks_without_support` | healpix-convert#25 and healpix-resample#71 (v2026.10.0) |
| float64 retry of the nearest-neighbour search | healpix-resample v2026.07.0 and later |
| `make_settings` (serializers per group) | legacy-converters#100 |
| healpix-resample from a checkout (`healpix_resample_path`) | healpix-resample v2026.10.0, pinned by g4e-desp-argo#29 |
| `_share_neighbour_search` (disk cache of neighbour searches) | Speed only; not needed for correctness |
