# Tiling overhead of the native Sentinel-2 product

Tile 109.8 km (10980 px at 10 m) on a 100 km step; swath 290 km.

- overlap alone (both axes, fully covered area): gamma = 1.206 -- real pixels stored twice, survives compression; a floor for gamma_tile
- swath tiling, raw bytes: 3-4 tiles across track (mean 4.00), across-track factor 1.514 [1.136, 1.514], along-track 1.098, total gamma = 1.662; no-data is 33.9% of the stored raster
- the illustrative gamma_tile = 1.3 of the paper lies between the two

