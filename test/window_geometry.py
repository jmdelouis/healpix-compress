"""How large a native window does *this* anchor actually need?

`paper_config.SAFE_WINDOW_PX` (1801 px) is a worst-case bound: it holds for any
anchor, including one landing exactly on a four-tile junction. Most anchors need
far less, and asking for more than necessary is not free -- a Sentinel-2 granule
only extends a few kilometres past its nominal MGRS square, so an over-large
window can fail to be covered by any single product even where a smaller one
would have been (this is what rejected Mato Grosso and Sao Paulo on 408 and 255
candidates respectively).

This module computes, for one anchor, the smallest window that contains a whole
level-10 tile plus the resampling halo. It needs no imagery: HEALPix cell
positions are a property of the grid.

Method: the tile's extent is set by its four corner cells. In NESTED indexing a
tile's children are Morton-ordered, so the corners are at fixed offsets inside
the child block -- (x=0,y=0) is offset 0, (1023,0) is 0x55555, (0,1023) is
0xAAAAA and (1023,1023) is 0xFFFFF. Candidate tiles are those covering a disc of
about one tile diagonal around the anchor, and everything is measured in an
azimuthal-equidistant frame centred on the anchor, so no scene CRS is needed.

    python window_geometry.py              # the window every region needs
    python window_geometry.py --regions sao_paulo mato_grosso
"""

from __future__ import annotations

import argparse

import numpy as np

import paper_config as cfg

# Morton corner offsets inside a delta=10 child block.
_CORNER_OFFSETS = np.array([
    0,                                        # (x=0,     y=0)
    int("01" * cfg.DELTA, 2),                 # (x=1023,  y=0)     0x55555
    int("10" * cfg.DELTA, 2),                 # (x=0,     y=1023)  0xAAAAA
    (1 << (2 * cfg.DELTA)) - 1,               # (x=1023,  y=1023)  0xFFFFF
], dtype=np.int64)


def _local_frame(lon0, lat0):
    import pyproj
    crs = pyproj.CRS.from_proj4(
        f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +x_0=0 +y_0=0 "
        f"+ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    fwd = pyproj.Transformer.from_crs(pyproj.CRS.from_epsg(4326), crs, always_xy=True)
    inv = pyproj.Transformer.from_crs(crs, pyproj.CRS.from_epsg(4326), always_xy=True)
    return fwd, inv


def candidate_parents(lon0, lat0, n_rings=6, n_angles=64):
    """Level-10 tiles covering a disc of ~one tile diagonal around the anchor."""
    import healpix_geo
    _fwd, inv = _local_frame(lon0, lat0)
    radii = np.linspace(0.0, 1.2 * cfg.TILE_SIDE_M * np.sqrt(2), n_rings)
    ang = np.linspace(0.0, 2 * np.pi, n_angles, endpoint=False)
    xs = np.concatenate([r * np.cos(ang) for r in radii])
    ys = np.concatenate([r * np.sin(ang) for r in radii])
    lon, lat = inv.transform(xs, ys)
    ipix = healpix_geo.nested.lonlat_to_healpix(
        np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64),
        cfg.LEVEL, ellipsoid="WGS84").astype(np.int64)
    return np.unique(ipix // cfg.N_CHILDREN)


def tile_half_extent_m(parent, lon0, lat0):
    """Half-extent (x, y) of one tile in the anchor-centred frame, in metres."""
    import healpix_geo
    fwd, _inv = _local_frame(lon0, lat0)
    ids = int(parent) * cfg.N_CHILDREN + _CORNER_OFFSETS
    clon, clat = healpix_geo.nested.healpix_to_lonlat(ids, cfg.LEVEL,
                                                      ellipsoid="WGS84")
    x, y = fwd.transform(np.asarray(clon, dtype=np.float64),
                         np.asarray(clat, dtype=np.float64))
    # Corner cell centres sit half a cell inside the tile; add that back.
    pad = cfg.DELTA_EQ
    return float(np.abs(x).max() + pad), float(np.abs(y).max() + pad)


def required_window_px(lon0, lat0, halo_px=32, safety_px=8, max_window=None):
    """Smallest odd-ish window (in native pixels) holding a whole tile + halo.

    Returns (window_px, best_parent, diagnostics). `window_px` is None when no
    candidate fits under `max_window`.
    """
    max_window = int(max_window or cfg.SAFE_WINDOW_PX + 4 * halo_px)
    rows = []
    for parent in candidate_parents(lon0, lat0):
        hx, hy = tile_half_extent_m(parent, lon0, lat0)
        half_px = max(hx, hy) / cfg.NATIVE_PIXEL_M + halo_px + safety_px
        need = int(2 * np.ceil(half_px))
        rows.append(dict(parent=int(parent), half_x_m=hx, half_y_m=hy,
                         window_px=need))
    rows.sort(key=lambda r: r["window_px"])
    if not rows:
        return None, None, rows
    best = rows[0]
    if best["window_px"] > max_window:
        return None, best["parent"], rows
    return best["window_px"], best["parent"], rows


# --- centring the acquisition on the tile rather than on the anchor ----------
def tile_for_anchor(lon0, lat0) -> int:
    """The level-10 tile containing the region anchor."""
    import healpix_geo
    ipix = healpix_geo.nested.lonlat_to_healpix(
        np.array([float(lon0)]), np.array([float(lat0)]), cfg.LEVEL,
        ellipsoid="WGS84").astype(np.int64)[0]
    return int(ipix // cfg.N_CHILDREN)


def plan_acquisition(lon0, lat0, halo_px=32, safety_px=8):
    """Where to centre the native window, and how big it must be.

    Centring on the *anchor* forces a window twice the tile's bounding box,
    because the anchor can sit anywhere inside its tile (worst case: a corner).
    Centring on the tile that contains the anchor needs only the bounding box
    plus the halo -- roughly half as wide, which is the difference between a
    window a single Sentinel-2 granule covers and one no granule covers.

    The anchor still selects the site and the tile; nothing here looks at any
    imagery, so the choice cannot favour either grid.
    """
    parent = tile_for_anchor(lon0, lat0)
    fwd, inv = _local_frame(lon0, lat0)
    import healpix_geo
    ids = parent * cfg.N_CHILDREN + _CORNER_OFFSETS
    clon, clat = healpix_geo.nested.healpix_to_lonlat(ids, cfg.LEVEL,
                                                      ellipsoid="WGS84")
    cx, cy = fwd.transform(np.asarray(clon, dtype=np.float64),
                           np.asarray(clat, dtype=np.float64))
    mx, my = float(np.mean(cx)), float(np.mean(cy))
    centre_lon, centre_lat = inv.transform(mx, my)
    half_m = max(float(np.abs(cx - mx).max()), float(np.abs(cy - my).max())) + cfg.DELTA_EQ
    window_px = int(2 * np.ceil(half_m / cfg.NATIVE_PIXEL_M + halo_px + safety_px))
    anchor_offset_m = float(np.hypot(mx, my))
    return dict(parent=parent, centre_lon=float(centre_lon),
                centre_lat=float(centre_lat), window_px=window_px,
                anchor_offset_m=anchor_offset_m, half_extent_m=half_m)


def main(argv=None) -> int:
    import pandas as pd
    import region_sites

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--halo", type=int, default=32)
    ap.add_argument("--max-window", type=int, default=None)
    args = ap.parse_args(argv)

    sites = region_sites.select(region_sites.load_sites(), args.regions, None, None)
    rows = []
    for site in sites.itertuples():
        need, parent, diag = required_window_px(site.lon, site.lat,
                                                halo_px=args.halo,
                                                max_window=args.max_window)
        plan = plan_acquisition(site.lon, site.lat, halo_px=args.halo)
        rows.append(dict(scene=site.scene, lat=site.lat, lon=site.lon,
                         anchor_centred_px=need if need else np.nan,
                         tile_centred_px=plan["window_px"],
                         tile_parent=plan["parent"],
                         anchor_offset_km=plan["anchor_offset_m"] / 1000.0,
                         n_candidates=len(diag)))
    df = pd.DataFrame(rows)
    print(f"safe worst-case bound: {cfg.SAFE_WINDOW_PX} px | "
          f"default acquisition: {cfg.WINDOW_PX} px")
    print(df.to_string(index=False))
    ok = df.anchor_centred_px.dropna()
    if len(ok):
        print(f"\nanchor-centred window: min {ok.min():.0f}, median {ok.median():.0f}, "
              f"max {ok.max():.0f} px")
    print(f"tile-centred window:   min {df.tile_centred_px.min():.0f}, "
          f"median {df.tile_centred_px.median():.0f}, "
          f"max {df.tile_centred_px.max():.0f} px  <-- what is acquired")

    # Cross-check: tiles already built must be the ones the anchors select.
    idx = cfg.TILE_DIR / "tiles_index.csv"
    if idx.exists():
        built = pd.read_csv(idx)[["scene", "parent_id"]]
        chk = df.merge(built, on="scene", how="inner")
        if len(chk):
            same = (chk.tile_parent == chk.parent_id)
            print(f"\nalready-built tiles: {int(same.sum())}/{len(chk)} use the "
                  f"tile their anchor selects")
            if not same.all():
                print(chk.loc[~same, ["scene", "tile_parent", "parent_id"]]
                      .to_string(index=False))
                print("  -> those were built from a differently centred window; "
                      "rebuild them with --force for a uniform protocol")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
