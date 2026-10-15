"""Tables of the ten whole Sentinel-2 L2A products (``s2_full_products.py``).

    python s2_ten_tiles_tables.py s2_ten_tiles.json --fetch     # copy the JSON results from the bucket
    python s2_ten_tiles_tables.py s2_ten_tiles.json             # write the tables

The JSON results of every finished tile are copied to
``results/s2_ten_tiles/<site>/`` (sizes, convert, quality, mirror, steps); the
tables go to ``tables/table_operational_*.tex`` of the manuscript and a Markdown
summary to ``results/s2_ten_tiles/summary.md``.

The manuscript's only source is ``tex/`` on the Makepaper branch of
jmdelouis/healpix-compress: ``$HCC_PAPER_TEX`` if set, else the ``tex/`` of the
repository this script sits in (Makepaper, where it lives in ``test/``), else a
Makepaper checkout next to this repository (healpix-compress-paper, where it
lives in ``benchmarks/``). The results sit in ``results/`` next to the script
(Makepaper) or one level up (healpix-compress-paper).
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = next((d for d in (HERE / "results" / "s2_ten_tiles", HERE.parent / "results" / "s2_ten_tiles")
                if d.is_dir()), HERE.parent / "results" / "s2_ten_tiles")


def _tex_dir() -> Path:
    if os.environ.get("HCC_PAPER_TEX"):
        return Path(os.environ["HCC_PAPER_TEX"]).expanduser()
    own = HERE.parent / "tex"
    return own if (own / "tgrs_healpix_compress.tex").exists() else HERE.parent.parent / "healpix-compress" / "tex"


TABLES = _tex_dir() / "tables"
NAMES = ("sizes", "convert", "quality", "mirror", "steps")
LANDSCAPE = {
    "paris_urban": "Urban, cropland", "alps_mountain_snow": "High mountain, snow",
    "andalusia_agriculture": "Irrigated cropland", "landes_forest": "Pine forest",
    "almeria_arid_greenhouses": "Arid, greenhouses", "brittany_coast_sea": "Coast, open sea",
    "rotterdam_port_delta": "Port, delta", "iceland_volcanic_glacier": "Volcanic, glacier",
    "camargue_wetlands": "Wetlands, lagoons", "finland_boreal_lakes": "Boreal forest, lakes",
}
RES = ("r10m", "r20m", "r60m")


def fetch(cfg: dict) -> list[str]:
    import fsspec

    fs, root = fsspec.core.url_to_fs(cfg["dest"].rstrip("/"), **cfg.get("storage_options", {}))
    done = []
    for site in cfg["tiles"]:
        if not fs.exists(f"{root}/{site}/DONE"):
            continue
        d = RESULTS / site
        d.mkdir(parents=True, exist_ok=True)
        for n in NAMES:
            try:
                (d / f"{n}.json").write_bytes(fs.cat_file(f"{root}/{site}/{n}.json"))
            except FileNotFoundError:
                pass
        done.append(site)
    return done


def load(cfg: dict) -> dict:
    out = {}
    for site in cfg["tiles"]:
        d = RESULTS / site
        if all((d / f"{n}.json").exists() for n in ("sizes", "convert", "quality")):
            out[site] = {n: json.loads((d / f"{n}.json").read_text()) for n in NAMES if (d / f"{n}.json").exists()}
    return out


def mb(b: float) -> str:
    return f"{b / 1e6:,.0f}".replace(",", "\\,")


def site_row(site: str, r: dict, tile: dict) -> dict:
    s = r["sizes"]
    groups = s["groups"]
    refl = [g for g in groups if g.startswith("measurements/reflectance/")]
    utm_refl = sum(groups[g]["utm_bytes"] for g in refl)
    hpx_refl = sum(groups[g]["healpix_bytes"] for g in refl)
    jp2 = sum(v["jpeg2000_lossless_bytes"] for v in s["utm_reflectance_lossless"].values())
    jxl = sum(v["jpegxl_lossless_bytes"] for v in s["utm_reflectance_lossless"].values())
    jp2_10 = sum(v["jpeg2000_lossless_bytes"] for k, v in s["utm_reflectance_lossless"].items() if "/r10m/" in k)
    hpx_10 = groups["measurements/reflectance/r10m"]["healpix_bytes"]
    other_utm = s["utm_total_bytes"] - utm_refl
    other_hpx = s["healpix_total_bytes"] - hpx_refl
    preset = sum(v.get("healpix_preset_default_bytes", v["healpix_bytes"]) for v in groups.values())
    return dict(site=site, tile=tile["tile"], date=tile["product"].split("_")[2][:8], landscape=LANDSCAPE.get(site, site),
                utm_total=s["utm_total_bytes"], hpx_total=s["healpix_total_bytes"], hpx_preset_total=preset,
                utm_refl=utm_refl, jp2_refl=jp2, jxl_refl=jxl, hpx_refl=hpx_refl, jp2_10=jp2_10, hpx_10=hpx_10,
                other_utm=other_utm, other_hpx=other_hpx,
                # wall time of the conversion step, as the progress log reports it
                convert_h=r.get("steps", {}).get("convert", r["convert"]["seconds"]) / 3600)


def table_sizes(rows: list[dict]) -> str:
    lines = [
        r"\begin{table*}[t]",
        r"\centering",
        r"\caption{Whole Sentinel-2 L2A products, EOPF Zarr as distributed against the HEALPix Zarr of the operational chain"
        r" (legacy-converters preset, reflectance under \texttt{hcc-sentinel2-pyramid} with $r=0.57$). Sizes in MB."
        r" JP2: the UTM reflectance under lossless JPEG~2000 per source chunk, the SAFE reference."
        r" ``Other'': masks, classification, atmosphere, probabilities, footprints, geometry, meteorology, metadata.}",
        r"\label{tab:operational_sizes}",
        r"\small",
        r"\begin{tabular}{llrrrrrrrr}",
        r"\toprule",
        r" & & \multicolumn{4}{c}{Reflectance (all bands, 10/20/60~m)} & \multicolumn{2}{c}{Other groups} & \multicolumn{2}{c}{Whole product} \\",
        r"\cmidrule(lr){3-6}\cmidrule(lr){7-8}\cmidrule(lr){9-10}",
        r"Landscape & Tile & EOPF & JP2 & HEALPix & $\times$JP2 & EOPF & HEALPix & EOPF & HEALPix \\",
        r"\midrule",
    ]
    for r in rows:
        lines.append(f"{r['landscape']} & {r['tile']} & {mb(r['utm_refl'])} & {mb(r['jp2_refl'])} & {mb(r['hpx_refl'])} & "
                     f"{r['hpx_refl'] / r['jp2_refl']:.2f} & {mb(r['other_utm'])} & {mb(r['other_hpx'])} & "
                     f"{mb(r['utm_total'])} & {mb(r['hpx_total'])} \\\\")
    tot = {k: sum(r[k] for r in rows) for k in ("utm_refl", "jp2_refl", "hpx_refl", "other_utm", "other_hpx", "utm_total", "hpx_total")}
    ratios = [r["hpx_refl"] / r["jp2_refl"] for r in rows]
    lines += [
        r"\midrule",
        f"All {len(rows)} & & {mb(tot['utm_refl'])} & {mb(tot['jp2_refl'])} & {mb(tot['hpx_refl'])} & "
        f"{tot['hpx_refl'] / tot['jp2_refl']:.2f} & {mb(tot['other_utm'])} & {mb(tot['other_hpx'])} & "
        f"{mb(tot['utm_total'])} & {mb(tot['hpx_total'])} \\\\",
        f"\\multicolumn{{5}}{{l}}{{Median (min--max) of $\\times$JP2 over tiles}} & {np.median(ratios):.2f} "
        f"({min(ratios):.2f}--{max(ratios):.2f}) & & & & \\\\",
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table*}",
    ]
    return "\n".join(lines) + "\n"


def band_quality(data: dict) -> dict:
    """Per resolution and band: pixel-weighted RMS over the tiles of the per-tile pooled RMSE."""
    out = {}
    for res in RES:
        g = f"measurements/reflectance/{res}"
        bands = sorted({b for r in data.values() for b in r["quality"].get(g, {}).get("per_band", {})})
        for b in bands:
            per = [r["quality"][g]["per_band"][b] for r in data.values() if b in r["quality"].get(g, {}).get("per_band", {})]
            if not per:
                continue

            def rms(key):
                return math.sqrt(np.mean([p[key] ** 2 for p in per]))

            ratio = [p["added_rmse_dn"] / p["transfer_rmse_dn"] for p in per if p["transfer_rmse_dn"] > 0]
            out[(res, b)] = dict(tiles=len(per), std=rms("native_std_dn"), transfer=rms("transfer_rmse_dn"),
                                 added=rms("added_rmse_dn"), total=rms("total_rmse_dn"),
                                 ratio_med=float(np.median(ratio)) if ratio else float("nan"),
                                 ratio_min=min(ratio) if ratio else float("nan"), ratio_max=max(ratio) if ratio else float("nan"))
    return out


def table_quality(q: dict, n_tiles: int, n_chunks: int) -> str:
    lines = [
        r"\begin{table}[t]",
        r"\centering",
        r"\caption{Error on the native UTM pixels of the ten products, in digital numbers (1~DN $=10^{-4}$ reflectance),"
        fr" on {n_chunks} random interior level-11 chunks per tile and resolution: standard deviation of the native pixels,"
        r" error of the PSF transfer alone, error added by \texttt{hcc-sentinel2-pyramid} (stored against uncompressed"
        r" transfer, both read back with the transfer's operator) and total. RMS over the tiles; the last column is the"
        r" median (min--max) over the tiles of added/transfer.}",
        r"\label{tab:operational_quality}",
        r"\small",
        r"\setlength{\tabcolsep}{3.5pt}",
        r"\begin{tabular}{llrrrrl}",
        r"\toprule",
        r"Res. & Band & $\sigma$ & Transfer & Added & Total & Added/transfer \\",
        r"\midrule",
    ]
    last = None
    for (res, b), v in q.items():
        if last is not None and res != last:
            lines.append(r"\midrule")
        last = res
        lines.append(f"{res[1:]} & {b.upper()} & {v['std']:.0f} & {v['transfer']:.1f} & {v['added']:.1f} & {v['total']:.1f} & "
                     f"{v['ratio_med']:.2f} ({v['ratio_min']:.2f}--{v['ratio_max']:.2f}) \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(lines) + "\n"


def markdown(rows: list[dict], q: dict) -> str:
    md = ["| Paysage | Tuile | Date | EOPF total (MB) | HEALPix total (MB) | HEALPix/EOPF | Réflectance JP2 (MB) | Réflectance HEALPix (MB) | ×JP2 | r10m ×JP2 | Autres EOPF→HEALPix (MB) | Conversion (h) |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['landscape']} | {r['tile']} | {r['date']} | {r['utm_total'] / 1e6:.0f} | {r['hpx_total'] / 1e6:.0f} | "
                  f"{r['hpx_total'] / r['utm_total']:.2f} | {r['jp2_refl'] / 1e6:.0f} | {r['hpx_refl'] / 1e6:.0f} | "
                  f"{r['hpx_refl'] / r['jp2_refl']:.2f} | {r['hpx_10'] / r['jp2_10']:.2f} | "
                  f"{r['other_utm'] / 1e6:.0f} → {r['other_hpx'] / 1e6:.0f} | {r['convert_h']:.1f} |")
    md += ["", "| Rés. | Bande | σ (DN) | Transfert (DN) | Ajouté codec (DN) | Total (DN) | Ajouté/transfert médiane (min–max) |",
           "|---|---|---|---|---|---|---|"]
    for (res, b), v in q.items():
        md.append(f"| {res} | {b} | {v['std']:.0f} | {v['transfer']:.1f} | {v['added']:.1f} | {v['total']:.1f} | "
                  f"{v['ratio_med']:.2f} ({v['ratio_min']:.2f}–{v['ratio_max']:.2f}) |")
    return "\n".join(md) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config", type=Path)
    ap.add_argument("--fetch", action="store_true")
    args = ap.parse_args(argv)
    cfg = json.loads(args.config.read_text())
    if args.fetch:
        print("fetched:", fetch(cfg))
    data = load(cfg)
    if not data:
        print("no finished tile yet")
        return 1
    rows = [site_row(s, r, cfg["tiles"][s]) for s, r in data.items()]
    q = band_quality(data)
    TABLES.mkdir(parents=True, exist_ok=True)
    (TABLES / "table_operational_sizes.tex").write_text(table_sizes(rows))
    (TABLES / "table_operational_quality.tex").write_text(table_quality(q, len(rows), cfg.get("quality_chunks", 12)))
    (RESULTS / "summary.md").write_text(markdown(rows, q))
    print(markdown(rows, q))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
