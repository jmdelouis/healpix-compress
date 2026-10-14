# Reproducing the HEALPix Sentinel-2 compression runs

This note covers the T31UDQ test and the T31UDQ+T31UDR mosaic.

It lists what you need to reproduce the HEALPix production of Sentinel-2 L2A 10 m reflectance used in our tests. The pipeline converts to HEALPix level 20, nested indexing on the WGS84 ellipsoid, chunked by level-11 cells, with the PSF resampler and the `hcc-sentinel2-pyramid` codec. Several pieces are still in open PRs, so use the branches below rather than `main`.

> The repositories `GRID4EARTH/healpix-compress`, `GRID4EARTH/legacy-converters` and `jmdelouis/healpix-compress-paper` are private: you need access to them.

## 1. Codec: `healpix-compress`, PR GRID4EARTH/healpix-compress#30

- Branch `sentinel2-jpegxl` on the fork `jmdelouis/healpix-compress`, head `abe831e`.
- It provides the `hcc-sentinel2` and `hcc-sentinel2-pyramid` Zarr codecs, plus the single-tile driver `test/eopf_s2_to_healpix.py` and its configurations.
- Install it with the JPEG XL extra:

```bash
git clone -b sentinel2-jpegxl git@github.com:jmdelouis/healpix-compress.git healpix-compress-s2
pip install -e "healpix-compress-s2[jpegxl]"
```

## 2. Conversion: `healpix-convert`, `main` plus three open PRs

- **PR #24, `fix-center-transform`.** Required. It handles x/y as pixel centres and fixes a half-pixel shift.
- **PR #25, `fix-empty-chunk-input`.** Required. It leaves a chunk empty when no input point can fill it, which happens at tile edges.
- **PR #31, `perf/batch-resample`.** Optional. It resamples all bands of a chunk in one call. The mosaic driver detects it and runs faster with it.
- PR #27 is not needed: the drivers set the codecs themselves.

Merging the three together has not been tested, so expect possible conflicts:

```bash
git clone git@github.com:GRID4EARTH/healpix-convert.git && cd healpix-convert
git fetch origin pull/24/head:pr24 pull/25/head:pr25 pull/31/head:pr31
git checkout -b s2-production origin/main
git merge --no-edit pr24 pr25 pr31
pip install -e .
```

## 3. Other dependencies

- **`healpix-resample`**: GRID4EARTH `main`, commit `533ceab`. It provides `PSFResampler` and needs torch. The CPU build is enough.
- **`legacy-converters`**: GRID4EARTH, commit `a2b67a4` (PR #92, 9 Oct 2026). It holds the Sentinel-2 L2A preset.
- **From PyPI**: `healpix-geo`, `pystac-client`, `s3fs`, `zarr>=3`, `numcodecs`, `requests`.

```bash
pip install -e "git+ssh://git@github.com/GRID4EARTH/healpix-resample.git@533ceab#egg=healpix-resample"
pip install "git+ssh://git@github.com/GRID4EARTH/legacy-converters.git@a2b67a4"
```

## 4. Production settings

These are already in the configuration files, so you should not need to change them:

- **Resampler**: PSF with `Npt=9`, `threshold=0.1`, `lam=0.1`, `max_iter=100`, `tol=1e-8`.
- **Codec**: `hcc-sentinel2-pyramid` in its default configuration. That means `relative_step=0.57`, `min_step=1`, `predictor_size=5` and `effort=7`. The quantification value (10000) and offset (−1000) are derived from the band's `scale_factor` and `add_offset`.
- **Bands**: `b02`, `b03`, `b04`, `b08` from `measurements/reflectance/r10m`.

## 5. Drivers and configurations

**Single tile**, with the uncompressed reference written alongside. Use `test/eopf_s2_to_healpix.py` from branch `sentinel2-jpegxl`, with `test/eopf_s2_t31udq_10m.json` (T31UDQ, 13 June 2025):

```bash
python test/eopf_s2_to_healpix.py test/eopf_s2_t31udq_10m.json
```

**Multi-tile mosaic.** Use `benchmarks/mosaic_convert_profile.py` from `jmdelouis/healpix-compress-paper`, `main`, commit `28402e3`, with `benchmarks/mosaic_20260812_R051.json` (T31UDQ and T31UDR, 12 August 2026, read from the EOPF service). Chunks that straddle two tiles get the points of both tiles. The run is resumable.

```bash
python benchmarks/mosaic_convert_profile.py benchmarks/mosaic_20260812_R051.json --workers 4 --threads 4
```

## 6. Output and credentials

- Both configurations write to `s3://grid4earth/test-compression/...` on the OVH GRA endpoint (`https://s3.gra.io.cloud.ovh.net`). They expect an AWS profile named `grid4earth` with write access. Change `output.s3` or `dest` and the profile if you write elsewhere.
- The EOPF STAC catalogue withdraws products over time. The 13 June 2025 product is no longer listed there, which is why its configuration reads the copy in `s3://grid4earth/public/eopf-mirror/`.

## 7. Validation

The notebooks comparing UTM, uncompressed HEALPix and compressed HEALPix are on branch `Makepaper` of `jmdelouis/healpix-compress`:

- `test/validation_T31UDQ_20250613.ipynb`
- `test/validation_mosaic_20260812_R051.ipynb`

They show four 512×512 patches, the reprojection to UTM with `PSFResampler`, and residual maps and spectra in reflectance. On the 13 June tile, the codec error is about 1e-3 reflectance RMS on B02–B04 and 2e-3 on B08 (`sizes.json`, whole tile). That is roughly 3.7–4.4 bits per cell.

> The 12 August mosaic currently on S3 is incomplete, with 909 of 2397 chunks written. Do not use it as a reference until it has been rerun.
