# Reproducing the paper

*When Does an Oversampled HEALPix Carrier Cost Less Than Native UTM?
A Rate–Distortion and Critical-Overhead Study on Sentinel-2*

This document is the recipe that rebuilds, from scratch, every number,
table and figure of the manuscript, and the recipe that built the dataset
deposited on Zenodo. Every figure quoted in the text is produced by a
script from the stored results, and `test/check_paper_numbers.py` refuses
to let the manuscript state a value the current results do not support.
If the last step of this document ends with *every checked figure from
the current results appears in the manuscript*, the manuscript and the
results agree.

All scripts live in `test/` and are run from that directory. Every long
step is resumable: a region already on disk is skipped unless `--force`
is given.

`test/reproduce.sh` runs steps 2.3 to 6 below in order and stops at the
first failure; it starts from the native windows on disk, so step 2.2
(the only one that needs the network) is run by hand first.

```bash
cd test
REGIONS="paris amazon" bash reproduce.sh    # smoke test on two regions
bash reproduce.sh                           # the 40 regions; log in results/reproduce.log
SKIP_TEX=1 bash reproduce.sh                # everything but the pdflatex run
FORCE=1 bash reproduce.sh                   # recompute every cached tile and sweep from the windows on disk
```

`FORCE=1` is the run that produces the deposited state: after a new
acquisition the cached tiles and sweeps belong to the previous windows,
and a chain run without it would mix the two.

---

## 0. Layout

| Path | Content |
|---|---|
| `tex/` | manuscript: `tgrs_healpix_compress.tex`, `sections/`, `tables/` (generated), `figures/` (generated) |
| `test/paper_config.py` | the single source of constants: level 20, δ = 10 (1024×1024 tiles), 637×637 equal-area crop, RMSE targets, γ = 1.30, paths |
| `test/data/regions/` | native Sentinel-2 B04 windows, one Zarr store per region, plus `acquisition_log.csv` |
| `test/data/tiles/` | derived tiles and crops, one `.npz` per region, plus `tiles_index.csv` |
| `test/results/` | every CSV / JSON / Markdown digest the sweeps and controls write; `bench_parts/` holds the per-region parts of the sweeps |
| `test/figures/` | figures as drawn; `make_figures.py`, `headline.py`, `gamma_analysis.py` and `fig_workflow.py` copy them into `tex/figures/` |

Three levels of reproduction, from the cheapest to the most complete:

1. **Tables, figures and the checker only** (minutes): unpack the Zenodo
   *results* archive into `test/results/` and run step 6. Nothing is
   recomputed.
2. **Full recomputation from the derived tiles** (about half a machine
   day): unpack the *tiles* archive into `test/data/tiles/` and run
   steps 3 to 6. No network is needed.
3. **Everything, including the acquisition** (adds ~15 min of download):
   run steps 2 to 6. Only step 2 needs the network.

---

## 1. Environment

```bash
git clone https://github.com/jmdelouis/healpix-compress.git
cd healpix-compress
git checkout tgrs-submission            # the tag the manuscript was built from
pip install -e ".[paper]"               # the evaluation stack (pyproject.toml, extra "paper")
```

The transfer operators come from the companion package `healpix-resample`
and the pipeline reads its frozen `data_manifest.csv` to get the 40
regions, so a checkout of that repository is required, at the tag the
sweeps ran against. The pin matters: the PSF width convention of
`healpix-resample` changed once, and every sweep records the checkout it
used in `results/bench_environment.json`.

```bash
# next to this repository, or anywhere with HEALPIX_RESAMPLE_ROOT pointing to it
git clone https://github.com/GRID4EARTH/healpix-resample.git ../healpix-resample
git -C ../healpix-resample checkout jstars-submission
pip install -e ../healpix-resample
export HEALPIX_RESAMPLE_ROOT=$(realpath ../healpix-resample)
```

Codecs and their Python bindings (all pulled in by the `paper` extra; the
compiled ones are easier to get from conda-forge if a wheel is missing):

| Codec in the paper | Package | Ladder swept |
|---|---|---|
| ZFP | `zfpy` | fixed-accuracy tolerance |
| SZ3 | `pysz` | absolute error bound |
| JPEG 2000 | `imagecodecs` (OpenJPEG) | compression ratio |
| JPEG XL | `imagecodecs` (libjxl) | quality level, then Butteraugli distance |
| DWT + Zstandard | `PyWavelets` + `zstandard` (`codec_lib.py`) | quantisation bits |

The HEALPix suite used is GRID4EARTH only (`healpix-geo`,
`healpix-resample`; `healpix-plot` for the maps); `healpy` is not used.

Check the installation before anything else:

```bash
cd test
python check_install.py          # every pipeline module importable, every third-party dependency present
python check_install.py --sha    # per-file sha256, to diff two machines
python paper_config.py           # prints the geometry: 637 px crop, 901 px bounding box, 1801 px safe window
python region_sites.py           # the 40 regions, classes, latitudes, pinned products
```

`region_sites.py` should print `40 regions | classes {agriculture: 10,
forest: 10, urban: 10, water: 10} | HEALPix regime {equatorial: 22, cap:
18} | lat -43.5..64.1 | pinned products 36/40`, followed by the note that
the companion study's own acquisition fell back to EOPF for four regions
(see 2.2).

---

## 2. Building the dataset (what the Zenodo deposit contains)

### 2.1 The 40 regions

The regions are those of the companion study: 40 Sentinel-2 L2A sites,
ten per land-cover class, latitudes from 43.5° S to 64.1° N. Their
anchors, classes and the Sentinel-2 product each patch came from are not
stored in this repository; they are read from `healpix-resample`'s frozen
`data_manifest.csv` (Zenodo DOI 10.5281/zenodo.22210945), so this study
cannot drift away from the published dataset. `region_sites.py` is the
only reader of that manifest.

### 2.2 Acquiring the native windows (network)

`build_region_windows.py` re-extracts, for each region, one B04 (10 m)
window from the same Sentinel-2 product as the companion study. The
window is centred on the level-10 HEALPix tile that contains the anchor
(not on the anchor itself) and sized per region from the tile's bounding
box plus a 32-pixel halo, 984 to 1194 native pixels across (median 1135);
`window_geometry.py` prints what each anchor requires without reading an
image.

The source is the `earth-search` STAC catalogue (Element 84; Sentinel-2
L2A cloud-optimised GeoTIFFs on AWS, read anonymously), queried over the
collections `sentinel-2-l2a`, `sentinel-2-c1-l2a` and
`sentinel-2-pre-c1-l2a`. For each region the script builds an ordered
list of candidates: the product pinned by the companion manifest, then
the same acquisition as the companion's patch when that patch came from
EOPF, then the search sorted by cloud cover. A candidate whose `red`
asset is a JP2 in the requester-pays archive is discarded before any
read, and a `RasterioIOError` moves to the next candidate. A product is
accepted only if the window has at least 99.9 % finite pixels
(`--min-finite`).

Four regions are a known case: `black_forest`, `finland_boreal`,
`lake_constance` and `lake_geneva`. For these, `earth-search` served the
asset from the requester-pays JP2 archive, the companion study's own
acquisition failed on the same four, and it fell back to ESA's EOPF Zarr
Sample Service. This pipeline resolves them through the candidate list
above (same acquisition, matched on MGRS tile and day, since EOPF and
`earth-search` identifiers differ); `acquisition_log.csv` records, per
region, which product and collection served it, whether it was pinned or
resolved, the cloud cover and the finite fraction.

```bash
unset AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY   # GDAL must read the public bucket anonymously
python build_region_windows.py --regions paris amazon   # two regions first (~10-30 s each)
python build_region_windows.py                          # the 40
python build_region_windows.py --check                  # offline inventory of what is on disk
```

Output: `data/regions/<scene>_data.zarr` (one per region, float32
reflectance, UTM zone of the product) and `data/regions/acquisition_log.csv`.
**This directory is the *windows* part of the Zenodo deposit**
(`make_data_archive.py --windows`). It is the only input that depends on
an external service, which is the reason to deposit it: a reader who
downloads it never touches `earth-search` or EOPF, and the four fallback
regions are the same bytes for everyone.

### 2.3 Deriving the tiles and crops

`build_tiles.py` turns each window into the arrays the sweeps compress,
and rejects a region (never patches it) if no level-10 parent tile lies
entirely inside the window with its halo, if the no-data fraction inside
the tile footprint exceeds `--max-nodata` (0.1 %), if the 637×637 crop
does not fit or contains no-data, or if the resampler does not return
exactly the 1 048 576 cells requested.

```bash
python build_tiles.py --regions paris amazon        # two regions first (a few minutes each: Delaunay)
python build_tiles.py                               # the 40, Clough--Tocher transfer (the reference carrier)
python build_tiles.py --resampler psf               # the PSF-aware transfer (lambda = 1e-3, 12.5 m FWHM)
python build_tiles.py --resampler linear            # the inverse-distance control
python build_tiles.py --audit                       # geometric audit: does each window still cover its tile?
```

Each `data/tiles/<scene>[__<kernel>].npz` holds the 1024×1024 level-20
tile (Morton de-interleaved from the nested index), the 637×637 equal-area
UTM crop cut from the same window, and the metadata (parent cell, UTM
zone, no-data counts); `tiles_index.csv` lists them. The level-19
carriers and the smoothed-native controls are not stored: `level19_rd.py`,
`level19_direct.py` and `support_controls.py` derive them from these
files. **This directory is the *tiles* part of the Zenodo deposit**
(`make_data_archive.py --tiles`, about 8 MB per region) and the starting
point of reproduction level 2.

One region is documented not to be covered by the inverse-distance
transfer: on Lake Taupo that kernel leaves 0.0154 % of the cells
unreachable, above the 0.01 % ceiling, so the region was left out rather
than filled, and Section V-L reports 39 regions for it. `paper_config.KNOWN_REJECTIONS`
lists this case; `build_tiles.py` and `run_codec_bench.py` report it as
expected and exit 0, while any other rejection still fails the run.

### 2.4 Audits that validate the dataset

```bash
python window_radiometry.py        # is each window still the product its tile came from (spectrum, response, MISMATCH flags)
python tile_crop_consistency.py    # tile/crop pairs against the window on disk (198 pairs over the transfers)
python codec_audit.py --grid both  # bytes written, closed, re-read by an independent decoder (two regions)
python operator_tests.py           # the HEALPix -> UTM read operator: constant, registration, smoothing, aliasing
```

`window_radiometry.csv` is read by the kernel-ablation table, which
excludes a region flagged `MISMATCH` (one lake in the current dataset),
and by the checker. The four files these audits write are the per-region
tables the manuscript says are "in the data repository", so
`make_data_archive.py` refuses to build a deposit without them.

### 2.5 Packaging for Zenodo

```bash
python make_data_archive.py --tag tgrs-submission --tiles --windows --zip
```

This writes `dist/healpix-compress-data-tgrs-submission/` and the same
as one `.zip` beside it, containing `regions.csv` (one row per region
with its pinned product), `results/` (every CSV, JSON and Markdown digest,
including `bench_environment.json` with the codec library versions and
the `healpix-resample` checkout of each run), `data/tiles/*.npz`, the
windows zipped one archive per region, `MANIFEST.csv` (path, bytes,
sha256 and the paper element each file feeds, so a reader can go from a
table to its data), a `README.md` with the provenance, the licences and
the command that rebuilds each table, and `zenodo_data.json`. The
creators of `zenodo_data.json` are read from `.zenodo.json` at the
repository root, so the code and data deposits cannot drift apart; add
that file before building the archive. Nothing is uploaded: attach the
`.zip` to the Zenodo deposit, then put its DOI in the Reproducibility
subsection of the paper and the tag's URL in `\repourl`.

Licences: results are CC-BY-4.0; the Sentinel-2 extracts are Copernicus
Sentinel data under the Copernicus licence; the code is Apache-2.0.

---

## 3. The compression sweeps

Each sweep writes one part file per (grid, region) into
`results/bench_parts/` and assembles them into `results/<grid>_codec_bench[__<kernel>][__<reader>].csv`.
A run that dies at region 27 resumes without redoing the first 26;
`--assemble` only rebuilds the CSVs. Each grid prints `N/M sweep points`
and writes a `.skips.json` per region; a `WARNING` that a codec produced
no point on one grid must be fixed before the statistics, because
`update_paper_numbers.py` refuses an asymmetric codec set.

```bash
python run_codec_bench.py --limit 2 --grid both      # dry run on two regions
python run_codec_bench.py --grid both                # 5 codecs x 71 operating points, both grids, 40 regions
python extend_sweep.py                               # adds the coarse end of the ladders the loosest target needs, reassembles
```

The same tiles scored against the native observation instead of the
carrier (Section V-I), with the two reads of Section V-G:

```bash
python run_codec_bench.py --grid hpx --common-reference --force                      # nearest-cell read -> rmse_common, common_reference_floor.csv
python run_codec_bench.py --grid hpx --common-reference --cr-mode bilinear --force   # interpolating read -> *__bilinear.csv
```

The other transfers (Section V-L):

```bash
python run_codec_bench.py --grid hpx --resampler psf        # -> hpx_codec_bench__psf.csv
python run_codec_bench.py --grid hpx --resampler linear     # -> hpx_codec_bench__linear.csv
python superres_cost.py --kernel psf --control linear --tex   # PSF-aware vs C--T at matched error, IDW as control -> superres_cost_{psf,linear}.csv, compactability.csv, transfer_cost_digest.json
python kernel_ablation.py --rd --limit 0                    # five kernels from the same window: transfer error, floors, spectra, sizes
```

The level-19 carriers (Section V-H, Appendix):

```bash
python level19_rd.py              # the aggregated level-19 carrier through the same codec sweep
python level19_direct.py          # the directly filled level-19 carrier and its floor
python level_choice.py            # what level 20 buys over level 19
```

The controls (Section V-E, Appendix):

```bash
python support_controls.py --control synthetic --fitted --transfer clough_tocher   # analytic field at the fitted slope and response (beta = 1.9), direct and through the paper's transfer
python support_controls.py --control smoothed                         # smoothed-native control at half a pixel (default --smooth-shift 0.5 0.5)
python support_controls.py --control smoothed --smooth-shift 0.25 0.25    # at a quarter of a pixel
python support_controls.py --control common shift                     # common-support and shifted-crop controls, all regions
python block_size_bias.py --limit 0                                   # 512^2 / 637^2 / 1024^2 per-pixel cost
python granule_overhead.py                                            # the tiling overhead from the published MGRS geometry
python censoring_sensitivity.py                                       # share of censored (bound) observations per target
```

The error budget of the transfers (Appendix D) and its decomposition:

```bash
python reader_floor.py                      # nearest, bilinear and PSF reads of the same uncompressed carrier
python resampling_error.py                  # C--T and PSF error against the observed and latent fields, at the measured noise
python resampling_error.py --noise 0        # the noiseless control -> resampling_error_noise0.csv
python alpha_budget.py --no-publish         # the targets in units of the transfer error
python cross_term.py                        # the cross term between transfer and coding errors
python diag_bicubic.py --all                # why the bicubic row of Table XIII(c) is marked invalid
python interscale.py                        # the quadtree-native (Haar / CDF 5/3) coder test of the appendix
```

Order of magnitude on one Pangeo-EOSC node: acquisition 10–30 s per
region, tiling a few minutes per region, the main sweep a few minutes per
region and grid, the common reference one to two minutes more per region;
half a machine day for the whole chain on 40 regions, unattended thanks to
the resume.

---

## 4. Statistics and digests

The log–log interpolation at each target, the oracle-codec selection, the
Wilcoxon signed-rank tests on log ratios, the Holm correction and the
latitude tests have exactly one implementation, `paper_stats.py`. The
digests are written by:

```bash
python update_paper_numbers.py                                  # paper_numbers_self.json/.md, the within-representation tables
python update_paper_numbers.py --metric common                  # paper_numbers_common.json, nearest-cell read
python update_paper_numbers.py --metric common --cr-mode bilinear   # paper_numbers_common_bilinear.json, interpolating read
python gamma_analysis.py            # the critical overhead gamma*(e) as a continuous curve, bootstrap CI; writes fig_verdict_vs_gamma
python headline.py                  # the headline medians at gamma = 1; writes table_headline and fig_headline
python reader_rescore.py            # the two reads side by side -> table_common_two_readers
```

The percentile bootstrap uses 10 000 resamples of regions with the seed
set in `paper_config.py`, so the intervals reproduce to the last digit.

---

## 5. Tables and figures

Tables are generated, never hand-edited. The three merged tables of the
paper are assembled from the digests above:

```bash
python merged_tables.py --only synthesis    # tex/tables/table_synthesis.tex   (the verdict per target)
python merged_tables.py --only floors       # tex/tables/table_floors.tex      (transfer-and-read floors)
python merged_tables.py --only transfer     # tex/tables/table_transfer.tex    (what the transfer costs, with the dagger notes)
python rebuild_tables.py                    # every other table from the stored results, without re-measuring
python rebuild_tables.py --list             # what each table is built from and by which script
```

The manuscript inputs fifteen tables: `table_storage_accounting`,
`table_synthesis`, `table_error_at_constant_size`, `table_fixed_codec`,
`table_codec_ranking`, `table_latitude_effect`, `table_support_controls`,
`table_floors`, `table_common_two_readers`, `table_level19_rd`,
`table_transfer`, `table_alpha_budget`, `table_resampling_error`,
`table_cross_term` and `table_censoring`.

Figures:

```bash
python fig_workflow.py        # Fig. 1, the workflow diagram (no data enters it)
python make_figures.py        # fig_utm_crop_footprints, fig_region_map, fig_storage_at_constant_rmse (and the diagnostic figures)
```

The manuscript includes six figures: `fig_workflow`,
`fig_utm_crop_footprints`, `fig_region_map`, `fig_headline` (from
`headline.py`), `fig_verdict_vs_gamma` (from `gamma_analysis.py`) and
`fig_storage_at_constant_rmse`. HEALPix data are drawn with
`healpix_plot` (`HealpixGrid` on cartopy axes).

---

## 6. The checker and the manuscript

```bash
python bold_glossary_terms.py ../tex/sections/*.tex   # idempotent: sets the glossary terms in italics
python check_paper_numbers.py                         # recomputes every quoted figure and finds it in the text
cd ../tex
pdflatex tgrs_healpix_compress && bibtex tgrs_healpix_compress && pdflatex tgrs_healpix_compress && pdflatex tgrs_healpix_compress
```

`check_paper_numbers.py` computes, from the files in `results/`, the
exact string each claim should contain (305 at the time of writing) and
asserts that it appears in the section its label names. A `MISSING` line
means that the prose is stale or that the run differs from the stored
results; it never means that the prose is right. The manuscript is never
edited by hand to make a check pass: the script that produces the number
is fixed, the table regenerated, the check re-run.

One part of the paper is outside this chain: the large-scale trial on
ten whole Sentinel-2 products (Appendix on the large-scale trial). It
ran the operational conversion (EOPF reader, converter,
`healpix-resample`, the `hcc-sentinel2` codecs) on an object store, and
the HEALPix products stay there. Two scripts measured its two tables:
`test/s2_full_products.py` converted the products and wrote their sizes
and errors to the store, and `test/s2_ten_tiles_tables.py` fetches these
results (`--fetch`) and writes the tables. `test/validation_s2_ten_tiles.ipynb`
draws its two figures. They need access to that store and the production
environment, so `reproduce.sh` does not run them, their inputs are not in
the Zenodo deposit, and the figures quoted in the appendix are not
asserted by `check_paper_numbers.py`. The paper says so in the appendix
and in its limitations.

---

## 7. Checklist before tagging

- [x] `check_install.py` reports the versions recorded in `results/bench_environment.json`
- [x] `check_paper_numbers.py` ends with *every checked figure from the current results appears in the manuscript*
- [x] `pdflatex` reports no undefined reference and no overfull box
- [x] `.zenodo.json` is at the repository root and `make_data_archive.py` builds without a missing required file
- [x] the `.zip` checksums in `MANIFEST.csv` match the Zenodo upload
- [x] `\repourl` in `tgrs_healpix_compress.tex` points to the tag, and the data DOI is in the Reproducibility subsection
