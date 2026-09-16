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

# Two different bounds, and confusing them is what makes a window look big
# enough when it is not.
#
# TILE_BBOX_PX -- the tile's own bounding box, i.e. the smallest window that can
# hold it *if the window happens to be centred on the tile*. A square rotated by
# 45 deg has a bounding box sqrt(2) times its own side.
TILE_BBOX_PX = int(math.ceil(TILE_SIDE_M * math.sqrt(2) / NATIVE_PIXEL_M))   # 901

# SAFE_WINDOW_PX -- the smallest window that is guaranteed to contain some whole
# level-10 tile whatever the anchor. Windows are centred on the region anchor,
# which is an arbitrary point of the HEALPix grid: in the worst case it lands on
# a four-tile junction, so the nearest tile centre is a half-diagonal away
# (4502 m) and that tile's far corner another half-diagonal beyond it. The
# window half-side must therefore reach 2 half-diagonals, i.e. the window must
# be twice the tile's bounding box.
TILE_HALF_DIAGONAL_M = TILE_SIDE_M * math.sqrt(2) / 2.0                      # 4502 m
SAFE_WINDOW_PX = int(math.ceil(4.0 * TILE_HALF_DIAGONAL_M / NATIVE_PIXEL_M))  # 1801
MIN_WINDOW_PX = SAFE_WINDOW_PX

# Window actually acquired per region: the safe bound plus room for the
# resampling halo on both sides.
WINDOW_PX = int(os.environ.get("HPXC_WINDOW_PX", 1920))

# --- storage accounting ------------------------------------------------------
GAMMA_UTM = 1.30              # assumed Sentinel-2/MGRS granule overhead
SHARED_RAW = N_CHILDREN * 4   # 4_194_304 B, fixed HEALPix tile raw size
UTM_RAW = UTM_SIDE_PX * UTM_SIDE_PX * 4

# --- rate-distortion protocol ------------------------------------------------
TARGET_RMSE = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2]
BUDGETS = [0.40, 0.20, 0.10, 0.05, 0.02]
MIN_SCENES = 3                # minimum evaluable regions for a reported target
# A target is only reported when this fraction of regions has it bracketed by
# two measured sweep points on both grids. Below it the row is a set of
# conservative bounds rather than an estimate: the 40% and 20% budgets are
# bracketed for 0% and 18% of regions respectively, which is what produced the
# uninterpretable ratios of the first version.
MIN_BRACKETED = 0.75

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


# The transfer the paper's carrier is built with. `build_tiles.DEFAULT_KERNEL`
# is the same string; window_radiometry.py checks the two agree at startup,
# because this module must not import build_tiles (build_tiles imports it).
CARRIER_KERNEL = "clough_tocher"


def window_mismatched_scenes(kernel: str | None = None, verbose: bool = False):
    """Regions whose native window is no longer the product their tile holds.

    Most of this study reads the stored tile and never reopens the window it
    was transferred from. The measurements that re-run a transfer do reopen
    it -- the kernel ablation, the resampling-error study and the level-19
    comparison -- and for those the window has to be the same product, which
    the geometric fit check in `build_tiles.tile_subwindow` does not
    establish: it passes whenever the window covers the tile, whatever the
    pixels contain. `window_radiometry.py` tests the radiometry and writes its
    verdict to results/window_radiometry.csv; every script that rebuilds a
    transfer reads it through this one function, so the three cannot drift
    apart on which regions are eligible.

    Returns an empty list when the audit has not been run, with a note rather
    than an error: a first pass of the pipeline has no audit yet, and refusing
    to proceed would be worse than proceeding with a stated caveat.
    """
    p = RESULTS_DIR / "window_radiometry.csv"
    if not p.exists():
        if verbose:
            print(f"note: {p.name} not found, so no region is excluded. This "
                  "measurement rebuilds a transfer from the native window and "
                  "is only valid where that window is the product the tile "
                  "came from; run window_radiometry.py.")
        return []
    import pandas as pd
    a = pd.read_csv(p)
    if "status" not in a.columns:
        if verbose:
            print(f"note: {p.name} predates the status column; re-run "
                  "window_radiometry.py")
        return []
    # A region has one window and one tile per transfer, built on different
    # days, so the audit is per (region, kernel): a window can match the
    # published Clough--Tocher tile and not the PSF-aware one. Each caller
    # asks about the kernel it actually rebuilds, and the default is the
    # carrier the paper uses.
    if "kernel" in a.columns:
        a = a[a.kernel == (kernel or CARRIER_KERNEL)]
    elif kernel not in (None, CARRIER_KERNEL) and verbose:
        print(f"note: {p.name} has no kernel column, so the verdict for "
              f"{kernel!r} is the carrier kernel's; re-run "
              "window_radiometry.py to audit every tile.")
    return sorted(a.loc[a.status != "ok", "scene"].astype(str))


def ensure_dirs() -> None:
    for d in (DATA_DIR, REGION_DIR, TILE_DIR, RESULTS_DIR, PARTS_DIR, FIGURES_DIR):
        d.mkdir(parents=True, exist_ok=True)


# --- LaTeX table helpers ----------------------------------------------------
# IEEEtran sets table captions in small caps, so a caption that explains the
# columns becomes a block of capitals above the data and is read by nobody.
# Every generated table therefore emits a short title as the caption and puts
# what the columns mean in a note under the table, in normal case.

def tex_caption(title: str, label: str) -> list[str]:
    """The two caption lines of a table: a short title, then its label."""
    title = " ".join(str(title).split())
    return [f"\\caption{{{title}}}", f"\\label{{{label}}}"]


def tex_note(text: str | None, wide: bool = False) -> list[str]:
    """What the columns mean -- which no longer goes under the table.

    These notes were set as a minipage under each table. IEEE's style manual
    asks for a short title above the table, centred and in small caps, set as
    an inverted pyramid; it does not provide for a paragraph of prose below
    one, and 42 such paragraphs made the paper hard to size and harder to
    read. The explanation belongs in the running text that introduces the
    table.

    The signature and the call sites are kept: each script still composes its
    note, which documents in one place what that table's columns are, and
    keeping it means the text is available to whoever rewrites the prose and
    the table generators need not be touched. Nothing is emitted.

    `text` and `wide` are accepted and ignored.
    """
    return []


def describe() -> str:
    return "\n".join([
        f"HEALPix level {LEVEL}, delta {DELTA} -> tile {TILE_SIDE}x{TILE_SIDE} cells",
        f"  delta_eq({LEVEL})   = {DELTA_EQ:.4f} m",
        f"  tile side           = {TILE_SIDE_M:.1f} m",
        f"  equal-area UTM crop = {UTM_SIDE_PX} x {UTM_SIDE_PX} native px "
        f"({UTM_RAW/1e6:.3f} MB raw)",
        f"  shared raw ref      = {SHARED_RAW} B ({SHARED_RAW/1e6:.3f} MB)",
        f"  tile bounding box   = {TILE_BBOX_PX} px (window centred on the tile)",
        f"  safe native window  = {SAFE_WINDOW_PX} px (any anchor, incl. a "
        f"four-tile junction)",
        f"  acquired window     = {WINDOW_PX} px "
        f"({'OK' if WINDOW_PX >= MIN_WINDOW_PX else 'TOO SMALL'})",
        f"  gamma_UTM           = {GAMMA_UTM}",
        f"  equatorial belt     = |lat| < {EQUATORIAL_LAT_DEG:.2f} deg",
    ])


if __name__ == "__main__":
    print(describe())
