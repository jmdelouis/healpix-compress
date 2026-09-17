"""The native product's tiling overhead, from the tiling geometry (revision 5, run 2.3).

The paper carries the overhead of the native layout as a parameter gamma and
uses gamma_tile = 1.30 as an illustrative value. Two reviews asked for a
measured one. This script derives what the Sentinel-2 tiling itself imposes,
from its published geometry, and can check the numbers on a real granule.

Geometry (ESA, SentiWiki "Sentinel-2 Products": tiles are "110 km x 110 km
ortho-images" on a grid "using a 100 km step ... in order to provide large
overlap with the neighbouring"; a 10 m band of a tile is 10980 x 10980 pixels,
i.e. 109.8 km). Three overheads follow, from the one that certainly applies
to the one that applies to raw bytes only:

    overlap     Every tile stores 109.8 km of ground on a 100 km step in both
                axes, so the archive holds each pixel of a fully covered area
                (109.8/100)^2 = 1.206 times. These are real pixels with real
                content, stored twice; the overhead survives compression
                unchanged. This is a floor for gamma_tile.

    swath       A 290 km swath is cut into tiles of 109.8 km on a 100 km
                step across track. Depending on where the swath edge falls
                on the grid, 3 or 4 tiles are needed; the edge tiles are
                partially covered and the rest of their raster is no-data
                (DN = 0). Averaged over the offset, the tiles store
                E[n] x 109.8 km of raster for 290 km of swath across track,
                and 1.098 along track. This is the raw-byte overhead of the
                tiled product, no-data included.

    compressed  No-data compresses to almost nothing, so under compression the
                swath overhead collapses towards the overlap overhead. How
                far is what the --granule measurement says: the fraction of
                no-data pixels in a real edge granule, and the compressed size
                of that granule's B04 relative to a fully covered one, if two
                paths are given.

Usage
-----
    python granule_overhead.py
    python granule_overhead.py --granule /path/to/T31TCJ_..._B04_10m.jp2 [--full /path/to/full_tile_B04.jp2]

Writes results/granule_overhead.json and results/granule_overhead.md.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import paper_config as cfg

TILE_PX = 10980            # pixels per side of a 10 m band
PIXEL_KM = 0.010
TILE_KM = TILE_PX * PIXEL_KM   # 109.8 km
STEP_KM = 100.0            # MGRS 100 km grid step
SWATH_KM = 290.0           # nominal MSI swath

SOURCES = {
    "tile_and_step": "ESA SentiWiki, Sentinel-2 Products (MSI overview): "
                     "'110 km x 110 km ortho-images ... using a 100 km step'",
    "swath": "ESA, Sentinel-2 mission description: 290 km swath",
    "tile_px": "product raster dimensions: 10980 x 10980 pixels at 10 m",
}


def overlap_overhead(tile_km=TILE_KM, step_km=STEP_KM):
    """Stored ground per unique ground, fully covered area, both axes."""
    return (tile_km / step_km) ** 2


def swath_overhead(tile_km=TILE_KM, step_km=STEP_KM, swath_km=SWATH_KM,
                   n_offsets=10000):
    """Raw-byte overhead of tiling a swath across track, averaged over where
    the swath edge falls on the grid, times the along-track overlap."""
    u = np.linspace(0.0, step_km, n_offsets, endpoint=False)
    # tile i covers [step*i, step*i + tile); it intersects [u, u + swath) iff
    # step*i < u + swath and step*i + tile > u
    i_lo = np.ceil((u - tile_km) / step_km + 1e-12)
    i_hi = np.ceil((u + swath_km) / step_km - 1e-12) - 1
    n_tiles = i_hi - i_lo + 1
    stored = n_tiles * tile_km
    across = stored / swath_km
    along = tile_km / step_km
    nodata_fraction = 1.0 - swath_km / stored          # per across-track row of tiles
    return dict(n_tiles_min=int(n_tiles.min()), n_tiles_max=int(n_tiles.max()),
                n_tiles_mean=float(n_tiles.mean()),
                across_track=float(across.mean()),
                across_track_min=float(across.min()),
                across_track_max=float(across.max()),
                along_track=float(along),
                total=float(across.mean() * along),
                nodata_fraction_mean=float(nodata_fraction.mean()))


def measure_granule(path: Path):
    """Dimensions and no-data fraction of one 10 m band file; compressed size."""
    import rasterio
    with rasterio.open(path) as src:
        arr = src.read(1)
        meta = dict(width=src.width, height=src.height, dtype=str(arr.dtype),
                    crs=str(src.crs), driver=src.driver,
                    res_m=[float(r) for r in src.res])
    nodata = int((arr == 0).sum())
    return dict(path=str(path), file_bytes=int(Path(path).stat().st_size),
                n_pixels=int(arr.size), nodata_pixels=nodata,
                nodata_fraction=nodata / arr.size,
                raw_bytes_uint16=int(arr.size * 2), **meta)


def summarize(d):
    ov, sw = d["overlap"], d["swath"]
    out = ["# Tiling overhead of the native Sentinel-2 product", "",
           f"Tile {d['tile_km']:.1f} km ({d['tile_px']} px at 10 m) on a "
           f"{d['step_km']:.0f} km step; swath {d['swath_km']:.0f} km.", "",
           f"- overlap alone (both axes, fully covered area): gamma = "
           f"{ov:.3f} -- real pixels stored twice, survives compression; a "
           f"floor for gamma_tile",
           f"- swath tiling, raw bytes: {sw['n_tiles_min']}-{sw['n_tiles_max']} "
           f"tiles across track (mean {sw['n_tiles_mean']:.2f}), across-track "
           f"factor {sw['across_track']:.3f} [{sw['across_track_min']:.3f}, "
           f"{sw['across_track_max']:.3f}], along-track {sw['along_track']:.3f}, "
           f"total gamma = {sw['total']:.3f}; no-data is "
           f"{100*sw['nodata_fraction_mean']:.1f}% of the stored raster",
           f"- the illustrative gamma_tile = {cfg.GAMMA_UTM} of the paper lies "
           f"between the two", ""]
    for g in d.get("granules", []):
        out.append(f"- granule {Path(g['path']).name}: {g['width']}x{g['height']} "
                   f"{g['dtype']}, no-data {100*g['nodata_fraction']:.1f}%, file "
                   f"{g['file_bytes']/1e6:.2f} MB for {g['raw_bytes_uint16']/1e6:.1f} "
                   f"MB raw uint16")
    if d.get("compressed_edge_over_full") is not None:
        out.append(f"- compressed size, edge granule / fully covered granule: "
                   f"{d['compressed_edge_over_full']:.3f} (no-data fraction "
                   f"{100*d['edge_nodata_fraction']:.1f}%)")
    return "\n".join(out) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tile-px", type=int, default=TILE_PX)
    ap.add_argument("--step-km", type=float, default=STEP_KM)
    ap.add_argument("--swath-km", type=float, default=SWATH_KM)
    ap.add_argument("--granule", type=Path, default=None,
                    help="a 10 m band file (JP2/COG) of any granule; reports its "
                         "dimensions and no-data fraction")
    ap.add_argument("--full", type=Path, default=None,
                    help="a fully covered granule's 10 m band, to compare "
                         "compressed sizes with --granule")
    args = ap.parse_args(argv)

    tile_km = args.tile_px * PIXEL_KM
    d = dict(tile_px=args.tile_px, tile_km=tile_km, step_km=args.step_km,
             swath_km=args.swath_km, sources=SOURCES,
             overlap=overlap_overhead(tile_km, args.step_km),
             swath=swath_overhead(tile_km, args.step_km, args.swath_km),
             gamma_illustrative=cfg.GAMMA_UTM, granules=[])
    for p in (args.granule, args.full):
        if p is not None:
            d["granules"].append(measure_granule(p))
    if args.granule is not None and args.full is not None:
        g, f = d["granules"]
        d["compressed_edge_over_full"] = g["file_bytes"] / f["file_bytes"]
        d["edge_nodata_fraction"] = g["nodata_fraction"]
    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (cfg.RESULTS_DIR / "granule_overhead.json").write_text(json.dumps(d, indent=1))
    md = summarize(d)
    (cfg.RESULTS_DIR / "granule_overhead.md").write_text(md, encoding="utf-8")
    print(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
