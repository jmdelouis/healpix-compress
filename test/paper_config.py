"""Shared configuration for the HEALPix-vs-UTM compression benchmark.

Single source of truth for every geometric and protocol constant that used to
be duplicated across compress_test.ipynb, hpx_codec_bench.ipynb,
utm_codec_bench.ipynb and stats_significance.ipynb.

Nothing here touches the network or the filesystem at import time.
"""

from __future__ import annotations

import math
import os
from pathlib import Path

# --- HEALPix carrier geometry ------------------------------------------------
LEVEL = 20                    # carrier level (6.22 m equal-area cell)
DELTA = 10                    # refinement levels below the parent tile
TILE_SIDE = 2 ** DELTA        # 1024 cells per side
N_CHILDREN = TILE_SIDE * TILE_SIDE
PARENT_LEVEL = LEVEL - DELTA  # 10

R_EARTH = 6371008.8           # m, matches eq. (1) of the paper
DELTA_EQ = math.sqrt(4 * math.pi * R_EARTH ** 2 / (12 * 4 ** LEVEL))  # ~6.2176 m
TILE_SIDE_M = TILE_SIDE * DELTA_EQ                                    # ~6366.8 m

# --- native Sentinel-2 grid --------------------------------------------------
NATIVE_PIXEL_M = 10.0
BAND = "b04"

# Equal-area native square compared against the HEALPix tile (Section III-B).
UTM_SIDE_PX = int(round(TILE_SIDE * DELTA_EQ / NATIVE_PIXEL_M))        # 637

# Smallest native window that can contain the tile whatever its rotation:
# a square rotated by 45 deg has a bounding box sqrt(2) times its own side.
MIN_WINDOW_PX = int(math.ceil(TILE_SIDE_M * math.sqrt(2) / NATIVE_PIXEL_M))  # 901

# Window actually acquired per region. Kept at the legacy 1512 px so the new
# 40-region bundle and the seven original scenes stay directly comparable, and
# so a footprint-exact (rather than equal-area) crop stays possible later.
WINDOW_PX = int(os.environ.get("HPXC_WINDOW_PX", 1512))

# --- storage accounting ------------------------------------------------------
GAMMA_UTM = 1.30              # assumed Sentinel-2/MGRS granule overhead
SHARED_RAW = N_CHILDREN * 4   # 4_194_304 B, fixed HEALPix tile raw size
UTM_RAW = UTM_SIDE_PX * UTM_SIDE_PX * 4

# --- rate-distortion protocol ------------------------------------------------
TARGET_RMSE = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2]
BUDGETS = [0.40, 0.20, 0.10, 0.05, 0.02]
MIN_SCENES = 3                # minimum evaluable regions for a reported target

COMPRESSOR_CODECS = ["zfp", "sz3", "jpeg2000_lossy", "jpegxl_lossy", "dwt_zstd"]
LOSSLESS_CODECS = ["jpeg2000_lossless", "jpegxl_lossless"]

# HEALPix cell shape regime: the equatorial belt is |z| < 2/3 on the sphere.
EQUATORIAL_LAT_DEG = math.degrees(math.asin(2.0 / 3.0))   # 41.81 deg

# --- paths -------------------------------------------------------------------
TEST_DIR = Path(__file__).resolve().parent
REPO_ROOT = TEST_DIR.parent
DATA_DIR = TEST_DIR / "data"
REGION_DIR = DATA_DIR / "regions"          # 1512 px windows, 40 regions
TILE_DIR = DATA_DIR / "tiles"              # derived HEALPix tiles + UTM crops
RESULTS_DIR = TEST_DIR / "results"
PARTS_DIR = RESULTS_DIR / "bench_parts"    # per-region partial sweeps (resume)
FIGURES_DIR = TEST_DIR / "figures"
TEX_FIGURES = REPO_ROOT / "tex" / "figures"
TEX_TABLES = REPO_ROOT / "tex" / "tables"

LEGACY_SCENES = ["urban", "water", "forest", "agriculture",
                 "snow_ice", "clouds", "ocean"]


def resample_root() -> Path:
    """Locate the healpix-resample checkout this study is synchronised with.

    Order: $HEALPIX_RESAMPLE_ROOT, then a sibling directory of this repo.
    """
    env = os.environ.get("HEALPIX_RESAMPLE_ROOT")
    candidates = []
    if env:
        candidates.append(Path(env))
    candidates += [REPO_ROOT.parent / "healpix-resample",
                   REPO_ROOT.parent / "healpix_resample"]
    for c in candidates:
        if (c / "healpix_resample").is_dir() and (c / "data_manifest.csv").exists():
            return c.resolve()
    raise FileNotFoundError(
        "healpix-resample checkout not found. Set HEALPIX_RESAMPLE_ROOT to it, "
        f"or place it next to this repository. Tried: "
        + ", ".join(str(c) for c in candidates)
    )


def ensure_dirs() -> None:
    for d in (DATA_DIR, REGION_DIR, TILE_DIR, RESULTS_DIR, PARTS_DIR, FIGURES_DIR):
        d.mkdir(parents=True, exist_ok=True)


def describe() -> str:
    return "\n".join([
        f"HEALPix level {LEVEL}, delta {DELTA} -> tile {TILE_SIDE}x{TILE_SIDE} cells",
        f"  delta_eq({LEVEL})   = {DELTA_EQ:.4f} m",
        f"  tile side           = {TILE_SIDE_M:.1f} m",
        f"  equal-area UTM crop = {UTM_SIDE_PX} x {UTM_SIDE_PX} native px "
        f"({UTM_RAW/1e6:.3f} MB raw)",
        f"  shared raw ref      = {SHARED_RAW} B ({SHARED_RAW/1e6:.3f} MB)",
        f"  min native window   = {MIN_WINDOW_PX} px (rotation-proof)",
        f"  acquired window     = {WINDOW_PX} px "
        f"({'OK' if WINDOW_PX >= MIN_WINDOW_PX else 'TOO SMALL'})",
        f"  gamma_UTM           = {GAMMA_UTM}",
        f"  equatorial belt     = |lat| < {EQUATORIAL_LAT_DEG:.2f} deg",
    ])


if __name__ == "__main__":
    print(describe())
