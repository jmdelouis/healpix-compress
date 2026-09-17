"""Does the raster's size or its non-power-of-two side bias the codecs? (revision 5, run 2.7)

The HEALPix tile is 1024 x 1024, the equal-area UTM crop 637 x 637. Codecs
work in blocks (ZFP 4 x 4, JPEG 2000 code-blocks, wavelet levels that pad an
odd side, see codec_lib.match_shape), so a reviewer asks whether part of the
measured HEALPix/UTM ratio is the crop's awkward size rather than its content.

The test compares each raster with its own central sub-blocks, so the texture
is the same and only the raster size and its divisibility change:

    UTM  637^2 (the crop)   vs  UTM 512^2 (central)   vs  UTM 640^2 (edge-replicated)
    HPX 1024^2 (the tile)   vs  HPX 637^2 (central)   vs  HPX 512^2 (central)

For each, the smallest lossless size per pixel and the oracle-best size per
pixel at two RMSE targets. If bytes per pixel of the 637^2 rasters match
those of the 512^2 and 1024^2 ones to within a percent, no block-size bias is
at work; if the 637^2 UTM crop costs more per pixel than its 512^2 sub-block
while the HEALPix 637^2 sub-block does not, the odd side has a price and the
digest gives it.

Usage
-----
    python block_size_bias.py --limit 5

Writes results/block_size_bias.csv, results/block_size_bias.json and
results/block_size_bias.md. No table: the result is one sentence.
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

import paper_config as cfg
import paper_stats as st
import codec_lib
import build_tiles
import region_sites

TARGETS = (1e-3, 1e-2)


def central(a, side):
    n = a.shape[0]
    o = (n - side) // 2
    return np.ascontiguousarray(a[o:o + side, o:o + side])


def padded(a, side):
    n = a.shape[0]
    p = side - n
    return np.pad(a, ((0, p), (0, p)), mode="edge")


def variants(scene):
    from hcc.compression.zorder import expand_tile
    cell_ids, cell_data, utm_crop, meta, npz = build_tiles.load_tile(scene)
    utm = np.asarray(utm_crop, dtype=np.float32)
    hpx = np.asarray(expand_tile(np.asarray(cell_data, dtype=np.float32)),
                     dtype=np.float32)
    n = utm.shape[0]
    return {
        f"utm_{n}": utm,
        "utm_512": central(utm, 512),
        "utm_640": padded(utm, 640),
        f"hpx_{hpx.shape[0]}": hpx,
        f"hpx_{n}": central(hpx, n),
        "hpx_512": central(hpx, 512),
    }


def sweep(arr, scene, codecs):
    rows = []
    for run in codec_lib.iter_runs(arr, codecs=codecs, on_error="warn"):
        rows.append(dict(scene=scene, codec=run.codec,
                         config=json.dumps(run.config, sort_keys=True),
                         total_bytes=int(run.total_bytes),
                         rmse=codec_lib.rmse(arr, run.reconstruction)))
    return pd.DataFrame(rows)


def region_rows(scene, codecs, lossless, targets):
    rows = []
    for name, arr in variants(scene).items():
        df = sweep(arr, scene, list(codecs) + list(lossless))
        n_px = float(arr.size)
        ll = df[df.codec.isin(lossless)]
        row = dict(scene=scene, variant=name, grid=name.split("_")[0],
                   side=int(arr.shape[0]), n_pixels=int(n_px),
                   lossless_bytes_per_px=(float(ll.total_bytes.min()) / n_px
                                          if len(ll) else np.nan))
        for t in targets:
            s = st.best_of_codecs_size(df, scene, list(codecs), t)
            br = st.oracle_bracketed(df, scene, list(codecs), t, "rmse")
            row[f"bytes_per_px@{t:g}"] = (s / n_px if np.isfinite(s) else np.nan)
            row[f"bracketed@{t:g}"] = bool(br)
        rows.append(row)
    return rows


def digest_of(df, targets):
    cols = ["lossless_bytes_per_px"] + [f"bytes_per_px@{t:g}" for t in targets]
    piv = {c: df.pivot(index="scene", columns="variant", values=c) for c in cols}
    utm_full = [v for v in df.variant.unique() if v.startswith("utm_")
                and v not in ("utm_512", "utm_640")][0]
    hpx_full = [v for v in df.variant.unique() if v.startswith("hpx_")
                and v not in ("hpx_512",) and v != utm_full.replace("utm", "hpx")][0]
    hpx_odd = utm_full.replace("utm", "hpx")
    pairs = {
        "utm_odd_over_utm_512": (utm_full, "utm_512"),
        "utm_640_over_utm_odd": ("utm_640", utm_full),
        "hpx_odd_over_hpx_512": (hpx_odd, "hpx_512"),
        "hpx_full_over_hpx_512": (hpx_full, "hpx_512"),
        "hpx_full_over_hpx_odd": (hpx_full, hpx_odd),
    }
    out = {"n_regions": int(df.scene.nunique()), "variants": sorted(df.variant.unique()),
           "ratios_bytes_per_px": {}}
    for c in cols:
        entry = {}
        for key, (a, b) in pairs.items():
            r = (piv[c][a] / piv[c][b]).dropna()
            if len(r):
                entry[key] = dict(median=float(r.median()), lo=float(r.min()),
                                  hi=float(r.max()), n=int(r.size))
        out["ratios_bytes_per_px"][c] = entry
    # the bias estimate: what the odd side costs UTM beyond what the same
    # reduction in size costs the HEALPix raster
    try:
        e = out["ratios_bytes_per_px"]["lossless_bytes_per_px"]
        out["odd_side_excess_lossless"] = (e["utm_odd_over_utm_512"]["median"]
                                           / e["hpx_odd_over_hpx_512"]["median"])
    except KeyError:
        out["odd_side_excess_lossless"] = None
    return out


def summarize(d):
    out = ["# Block-size and odd-side bias", "",
           f"{d['n_regions']} regions; bytes per pixel, median over regions of "
           "the per-region ratio.", "",
           "| quantity | pair | median | range |", "|---|---|---|---|"]
    for c, entry in d["ratios_bytes_per_px"].items():
        for key, v in entry.items():
            out.append(f"| {c} | {key} | {v['median']:.4f} | "
                       f"[{v['lo']:.4f}, {v['hi']:.4f}] |")
    x = d.get("odd_side_excess_lossless")
    if x is not None:
        out += ["", f"Odd-side excess (UTM 637/512 divided by HPX 637/512, "
                    f"lossless): {x:.4f}, i.e. {100*(x-1):+.2f}% per pixel."]
    return "\n".join(out) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=5)
    ap.add_argument("--codecs", nargs="*", default=cfg.COMPRESSOR_CODECS)
    ap.add_argument("--lossless", nargs="*", default=cfg.LOSSLESS_CODECS)
    ap.add_argument("--targets", type=float, nargs="*", default=list(TARGETS))
    args = ap.parse_args(argv)

    scenes = args.regions or list(region_sites.load_sites()["scene"])
    scenes = [s for s in scenes if (cfg.TILE_DIR / f"{s}.npz").exists()]
    if args.limit:
        # one region per class first, so five regions are not five lakes
        by_cls = {}
        for s in scenes:
            by_cls.setdefault(s.split("__")[1], []).append(s)
        picked, i = [], 0
        while len(picked) < args.limit and any(by_cls.values()):
            for cls in sorted(by_cls):
                if by_cls[cls] and len(picked) < args.limit:
                    picked.append(by_cls[cls].pop(0))
            i += 1
        scenes = picked
    if not scenes:
        print(f"no tiles under {cfg.TILE_DIR}; run build_tiles.py first")
        return 1
    print(f"{len(scenes)} region(s): {scenes}")
    rows = []
    for scene in scenes:
        print(f"  {scene}")
        rows += region_rows(scene, args.codecs, args.lossless, args.targets)
    df = pd.DataFrame(rows)
    cfg.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(cfg.RESULTS_DIR / "block_size_bias.csv", index=False)
    d = digest_of(df, args.targets)
    (cfg.RESULTS_DIR / "block_size_bias.json").write_text(json.dumps(d, indent=1))
    md = summarize(d)
    (cfg.RESULTS_DIR / "block_size_bias.md").write_text(md, encoding="utf-8")
    print(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
