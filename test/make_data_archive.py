"""Build the data deposit that the paper cites (revision 6, item 1.7).

Four tables were taken out of the manuscript in revision 6 and replaced by a
sentence saying their per-region values are "in the data repository": the
region list, the reading-operator tests, the byte-and-type audit and the
window radiometry. That sentence is a commitment, so this script refuses to
build an archive that does not carry them.

What it collects
----------------
    regions.csv         one row per region (scene, class, coordinates, the
                        pinned Sentinel-2 product), from region_sites, which
                        reads healpix-resample's frozen data_manifest.csv
    results/            every CSV, JSON and Markdown digest the sweeps and
                        controls wrote, including bench_environment.json (the
                        codec library versions and the healpix-resample
                        checkout each run used)
    MANIFEST.csv        path, bytes, sha256 and the paper element each file
                        feeds, so a reader can go from a table to its data
    README.md           provenance, licences, and the command that rebuilds
                        each table from these files
    zenodo_data.json    deposit metadata (the creators are read from the
                        repository's .zenodo.json, so the two deposits cannot
                        drift apart)

    --tiles             also the stored level-20 tiles and UTM crops
                        (data/tiles/*.npz, about 8 MB per region)
    --windows           also the native windows (data/regions/*.zarr, large;
                        zipped one archive per region)

Usage
-----
    python make_data_archive.py --tag tgrs-submission
    python make_data_archive.py --tag tgrs-submission --tiles --zip

Writes dist/healpix-compress-data-<tag>/ and, with --zip, the same as a
single .zip beside it. Nothing is uploaded: the archive is what you attach
to the Zenodo deposit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd

import paper_config as cfg
import region_sites

# What each results file feeds, by glob. The value is the paper element, named
# by its LaTeX label so the mapping survives a renumbering of the tables.
FEEDS = {
    "regions.csv": "the region list (Section IV-A; table removed from the "
                   "manuscript in revision 6)",
    "operator_tests.csv": "the reading-operator characterisation "
                          "(Section V-J, Appendix B-B; table removed)",
    "codec_audit.csv": "the byte and type audit of the codec adapters "
                       "(Appendix B-C; table removed)",
    "window_radiometry.csv": "is each window still the product its tile came "
                             "from (Appendix B-E; table removed)",
    "hpx_codec_bench.csv": "the HEALPix-side rate-distortion sweep behind "
                           "every ratio (tab:synthesis, tab:error_at_size)",
    "utm_codec_bench.csv": "the native-UTM sweep behind every ratio",
    "hpx_codec_bench__*.csv": "the same sweep on another transfer or reader "
                              "(suffix = transfer, then reader)",
    "paper_numbers_self.json": "tab:synthesis, tab:error_at_size, "
                               "tab:fixed_codec, tab:latitude_effect, "
                               "tab:codec_ranking, fig:headline",
    "paper_numbers_common*.json": "tab:common_two_readers, Section V-I",
    "gamma_analysis_self.csv": "the critical overhead and f(e) "
                               "(tab:synthesis, fig:verdict_vs_gamma)",
    "headline.json": "fig:headline",
    "support_controls_*.csv": "the support, synthetic and smoothed-native "
                              "controls (tab:support_controls, Section V-G)",
    "common_reference_floor*.csv": "the transfer-and-read floors "
                                   "(tab:floors, Section V-I)",
    "reader_floor.csv": "the three reads of the same carrier (tab:floors, "
                        "Section V-J)",
    "reader_rescore.csv": "tab:common_two_readers",
    "level19_rd*.csv": "the level-19 carrier under compression "
                       "(tab:level19_rd, Section V-K)",
    "level19_rd*.json": "tab:level19_rd, Section V-K",
    "level19_rd_diagnose.md": "why the aggregate column is flat (Section V-K)",
    "level19_direct*.csv": "the two level-19 carriers (tab:floors, Section V-L)",
    "level19_direct*.json": "tab:floors, Section V-L",
    "level_choice.csv": "what level 20 buys over level 19 (Section V-K)",
    "censoring_*.csv": "the censored-point sensitivity (tab:censoring, "
                       "Appendix B-D)",
    "censoring_digest.json": "tab:censoring, Appendix B-D",
    "superres_cost_*.csv": "what the transfer costs to store (tab:transfer, "
                           "Section V-M)",
    "transfer_cost_digest.json": "tab:transfer, Section V-M",
    "compactability.csv": "tab:transfer, second block",
    "kernel_ablation.csv": "tab:transfer, third block",
    "resampling_error*.csv": "each transfer against the reference it answers "
                             "to (tab:resampling_error, Appendix C)",
    "cross_term.csv": "the cross term between transfer and compression error "
                      "(tab:cross_term, Appendix D)",
    "interscale*.csv": "the quadtree-native coder test (Appendix A-B)",
    "granule_overhead.json": "the measured tiling overhead (Section IV-C)",
    "granule_overhead.md": "the measured tiling overhead (Section IV-C)",
    "block_size_bias.json": "the odd-side divisibility control (Section IV-C)",
    "bench_environment.json": "the codec library versions and the "
                              "healpix-resample checkout each run used "
                              "(Section VI-C)",
}

# The manuscript promises these by name. Building an archive without them
# would publish a broken cross-reference.
REQUIRED = ("regions.csv", "operator_tests.csv", "codec_audit.csv",
            "window_radiometry.csv")

DATA_LICENCE = (
    "The derived results (every CSV, JSON and Markdown file under results/) "
    "are released under CC-BY-4.0. The Sentinel-2 L2A extracts they were "
    "computed from are Copernicus Sentinel data, free to use under the "
    "Copernicus licence; the product identifier of each region is in "
    "regions.csv.")


def feeds_for(name: str) -> str:
    """The paper element a results file feeds; exact match, then glob, then
    the Markdown digest of a described CSV or JSON (same stem)."""
    if name in FEEDS:
        return FEEDS[name]
    for pattern, what in FEEDS.items():
        if "*" in pattern and Path(name).match(pattern):
            return what
    if name.endswith(".md"):
        stem = name[:-3]
        for ext in (".csv", ".json"):
            what = feeds_for(stem + ext)
            if what:
                return f"{what} (Markdown digest)"
    return ""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def code_version() -> dict:
    """Which checkout of this repository wrote the archive."""
    out = {}
    for key, cmd in (("commit", ["git", "rev-parse", "HEAD"]),
                     ("describe", ["git", "describe", "--tags", "--dirty",
                                   "--always"])):
        try:
            out[key] = subprocess.run(cmd, cwd=cfg.REPO_ROOT, check=True,
                                      capture_output=True, text=True
                                      ).stdout.strip()
        except Exception:                                     # noqa: BLE001
            out[key] = "unknown"
    return out


def write_regions(dest: Path) -> Path:
    """One row per region, with the Sentinel-2 product each one is pinned to."""
    sites = region_sites.load_sites()
    keep = [c for c in ("scene", "scene_class", "region_id", "lat", "lon",
                        "hpx_regime", "catalog", "product_id")
            if c in sites.columns]
    out = dest / "regions.csv"
    sites[keep].to_csv(out, index=False)
    return out


def collect_results(dest: Path) -> list[Path]:
    """Every digest the sweeps and controls wrote, flat, under results/."""
    src = cfg.RESULTS_DIR
    if not src.exists():
        raise SystemExit(f"{src} does not exist; nothing to archive")
    target = dest / "results"
    target.mkdir(parents=True, exist_ok=True)
    copied = []
    for p in sorted(src.iterdir()):
        if p.is_dir() or p.suffix not in (".csv", ".json", ".md"):
            continue
        if ".pre_extension" in p.name or "DO_NOT_PUBLISH" in p.name:
            continue
        shutil.copy2(p, target / p.name)
        copied.append(target / p.name)
    return copied


def collect_tiles(dest: Path) -> list[Path]:
    target = dest / "tiles"
    target.mkdir(parents=True, exist_ok=True)
    copied = []
    for p in sorted(cfg.TILE_DIR.glob("*.npz")):
        shutil.copy2(p, target / p.name)
        copied.append(target / p.name)
    return copied


def collect_windows(dest: Path) -> list[Path]:
    """One zip per region: a .zarr store is a directory, not a file."""
    target = dest / "windows"
    target.mkdir(parents=True, exist_ok=True)
    copied = []
    for store in sorted(cfg.REGION_DIR.glob("*_data.zarr")):
        if not store.is_dir():
            continue
        out = target / f"{store.stem}.zip"
        shutil.make_archive(str(out.with_suffix("")), "zip", root_dir=store)
        copied.append(out)
    return copied


def readme(dest: Path, version: dict, counts: dict) -> str:
    lines = [
        "# Data and results for *When Does an Oversampled HEALPix Carrier "
        "Cost Less Than Native UTM?*", "",
        "This deposit holds the measurements behind every table, figure and "
        "quoted value of the manuscript, plus the per-region values of four "
        "tables the manuscript cites here rather than printing.", "",
        f"Written from healpix-compress at commit `{version['commit']}` "
        f"(`{version['describe']}`). The code is archived separately; its DOI "
        "is in the related identifiers of this deposit.", "",
        "## Contents", "",
        "| Path | Files | What it is |",
        "|---|---|---|",
        "| `regions.csv` | 1 | one row per evaluation region: class, "
        "coordinates, HEALPix cell-shape regime and the exact Sentinel-2 "
        "product the window was cut from |",
        f"| `results/` | {counts.get('results', 0)} | every CSV, JSON and "
        "Markdown digest the sweeps and controls wrote |",
    ]
    if counts.get("tiles"):
        lines.append(f"| `tiles/` | {counts['tiles']} | the stored level-20 "
                     "carrier tiles and their equal-area UTM crops, one "
                     "`.npz` per region |")
    if counts.get("windows"):
        lines.append(f"| `windows/` | {counts['windows']} | the native "
                     "Sentinel-2 windows each tile was built from, one zipped "
                     "Zarr store per region |")
    lines += [
        "| `MANIFEST.csv` | 1 | path, size, SHA-256 and the paper element "
        "each file feeds |", "",
        "## Reading a table back to its data", "",
        "`MANIFEST.csv` names, for each file, the table or section it feeds. "
        "The tables themselves are rebuilt from these files by the code "
        "deposit: `update_paper_numbers.py` (the per-target tables), "
        "`merged_tables.py` (synthesis, floors, transfer), "
        "`gamma_analysis.py` (critical overhead), `reader_rescore.py` (the "
        "common-reference readers), `level19_rd.py`, `alpha_budget.py`, "
        "`cross_term.py` and `make_figures.py` (the figures). "
        "`check_paper_numbers.py` recomputes every value quoted in the text "
        "from these same files and fails if one of them no longer supports "
        "it.", "",
        "## The four tables the manuscript cites here", "",
        "- the region list: `regions.csv`;",
        "- the reading operator on single-frequency fields: "
        "`results/operator_tests.csv`;",
        "- the byte and type audit of the codec adapters: "
        "`results/codec_audit.csv`;",
        "- the window radiometry audit: `results/window_radiometry.csv`.", "",
        "## Licence and provenance", "",
        DATA_LICENCE, "",
        "The regions, acquisitions and radiometric conventions are shared "
        "with the companion study of PSF-aware resampling onto HEALPix, whose "
        "own dataset is archived under its own DOI; see the related "
        "identifiers.", "",
    ]
    return "\n".join(lines)


def zenodo_metadata(version: dict) -> dict:
    """Deposit metadata, with the creators taken from the code deposit."""
    creators = []
    for name in (".zenodo.json", "zenodo_code.json"):
        src = cfg.REPO_ROOT / name
        if src.exists():
            creators = json.loads(src.read_text(encoding="utf-8")).get("creators", [])
            break
    if not creators:
        print("  no .zenodo.json at the repository root: zenodo_data.json "
              "will carry no creators, fill them in the Zenodo form")
    return {
        "upload_type": "dataset",
        "title": ("Data and results for \"When Does an Oversampled HEALPix "
                  "Carrier Cost Less Than Native UTM? A Rate-Distortion and "
                  "Critical-Overhead Study on Sentinel-2\""),
        "description": (
            "Measurements behind every table, figure and quoted value of the "
            "manuscript: the rate-distortion sweeps of five general-purpose "
            "codecs on HEALPix level-20 carrier tiles and on area-matched "
            "native Sentinel-2 UTM crops, the critical-overhead analysis, the "
            "support, transfer and reading-operator controls, and the "
            "per-region values of four tables the manuscript cites rather "
            "than printing (the region list, the reading-operator tests, the "
            "byte and type audit, and the window radiometry audit). "
            "MANIFEST.csv maps every file to the table or section it feeds. "
            f"Written from healpix-compress at commit {version['commit']}."),
        "access_right": "open",
        "license": "cc-by-4.0",
        "language": "eng",
        "keywords": ["HEALPix", "DGGS", "data compression", "rate-distortion",
                     "Sentinel-2", "resampling", "Earth observation"],
        "creators": creators,
        "notes_for_the_author": [
            "Add in the Zenodo form, once the DOIs exist: related identifier "
            "'isSupplementTo' the TGRS article, 'isCompiledBy' the code "
            "deposit, 'isDerivedFrom' the Copernicus Sentinel-2 products.",
            "Add the funding entry with the funder picker (ESA contract "
            "4000147951/25/I-NS).",
        ],
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", default="tgrs-submission",
                    help="suffix of the archive directory (default: "
                         "tgrs-submission)")
    ap.add_argument("--tiles", action="store_true",
                    help="also archive the stored tiles and UTM crops")
    ap.add_argument("--windows", action="store_true",
                    help="also archive the native windows (large)")
    ap.add_argument("--zip", action="store_true",
                    help="also write the archive as a single .zip")
    ap.add_argument("--out", default=None,
                    help="output directory (default: dist/ beside the repo)")
    args = ap.parse_args(argv)

    root = Path(args.out) if args.out else cfg.REPO_ROOT / "dist"
    dest = root / f"healpix-compress-data-{args.tag}"
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)

    version = code_version()
    print(f"healpix-compress at {version['describe']} ({version['commit'][:12]})")

    regions = write_regions(dest)
    results = collect_results(dest)
    counts = {"results": len(results)}
    files = [regions, *results]
    if args.tiles:
        t = collect_tiles(dest)
        counts["tiles"] = len(t)
        files += t
    if args.windows:
        w = collect_windows(dest)
        counts["windows"] = len(w)
        files += w

    missing = [n for n in REQUIRED
               if not any(p.name == n for p in files)]
    if missing:
        shutil.rmtree(dest)
        raise SystemExit(
            "the manuscript cites these by name and they are absent from "
            f"results/: {', '.join(missing)}.\nRun the scripts that write "
            "them (operator_tests.py, codec_audit.py, window_radiometry.py) "
            "before building the archive; an archive without them would "
            "publish a broken cross-reference.")

    rows = []
    for p in files:
        rows.append(dict(path=str(p.relative_to(dest)).replace("\\", "/"),
                         bytes=p.stat().st_size, sha256=sha256(p),
                         feeds=feeds_for(p.name)))
    man = pd.DataFrame(rows).sort_values("path")
    man.to_csv(dest / "MANIFEST.csv", index=False)

    (dest / "README.md").write_text(readme(dest, version, counts),
                                    encoding="utf-8")
    (dest / "zenodo_data.json").write_text(
        json.dumps(zenodo_metadata(version), indent=2), encoding="utf-8")

    total = sum(r["bytes"] for r in rows)
    undescribed = man[man.feeds == ""]
    print(f"{len(rows)} file(s), {total / 1e6:.1f} MB -> {dest}")
    if len(undescribed):
        print(f"  {len(undescribed)} file(s) with no entry in FEEDS: "
              f"{', '.join(undescribed.path.head(10))}")
    if args.zip:
        archive = shutil.make_archive(str(dest), "zip", root_dir=dest.parent,
                                      base_dir=dest.name)
        print(f"zip -> {archive} ({Path(archive).stat().st_size / 1e6:.1f} MB)")
    print("\nNext: attach the archive to a Zenodo deposit; the metadata to "
          "paste is in zenodo_data.json, and reviews/ZENODO.md has the order "
          "in which the two deposits and the paper cite each other.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
