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
    python make_data_archive.py --tag tgrs-submission-v2
    python make_data_archive.py --tag tgrs-submission-v2 --tiles --zip
    python make_data_archive.py --tag tgrs-submission-v2 --zip \
        --reuse "https://zenodo.org/records/23388572/files/healpix-compress-data-tgrs-submission.zip?download=1"

--reuse copies tiles/ and windows/ byte for byte from an earlier archive
(a local .zip or its Zenodo URL, read by HTTP range requests), each member
checked against that archive's own MANIFEST.csv.

--update-of builds an UPDATE archive instead: regions.csv, results/,
README.md and a MANIFEST.csv that also lists, with their SHA-256, the tiles
and windows of the earlier archive, which stays in the deposit unchanged.
Version 2 of the deposit was built that way, so that its 2 GB of tiles and
windows did not move:

    python make_data_archive.py --tag tgrs-submission-v2 --zip \
        --update-of "https://zenodo.org/records/23388572/files/healpix-compress-data-tgrs-submission.zip?download=1"

Writes dist/healpix-compress-data-<tag>/ and, with --zip, the same as a
single .zip beside it. Nothing is uploaded: the archive is what you attach
to the Zenodo deposit.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd

import paper_config as cfg
import region_sites

# What each results file feeds, by glob. The value is the paper element, named
WITHDRAWN = ("withdrawn: scored with a psf read whose kernel was half as wide "
             "as the transfer's (see README); no longer used by the paper")

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
    "hpx_codec_bench__psf__psf.csv": WITHDRAWN,
    "paper_numbers_common_psf_psf.json": WITHDRAWN,
    "paper_numbers_common_psf_psf.md": WITHDRAWN,
    "common_reference_floor__psf__psf.csv": WITHDRAWN,
    "hpx_codec_bench__*.csv": "the same sweep on another transfer or reader "
                              "(suffix = transfer, then reader)",
    "paper_numbers_self.json": "tab:synthesis, tab:error_at_size, "
                               "tab:fixed_codec, fig:headline, and the "
                               "latitude and codec-ranking figures quoted in "
                               "the text",
    "paper_numbers_common*.json": "tab:common_two_readers (common reference)",
    "psf_matched_floor.json": "the PSF-aware tiles read through their own "
                              "forward operator on four regions (Section on "
                              "the reading operator, trial appendix)",
    "gamma_analysis_self.csv": "the critical overhead and f(e) "
                               "(tab:synthesis, fig:verdict_vs_gamma)",
    "headline.json": "fig:headline",
    "support_controls_*.csv": "the support, synthetic and smoothed-native "
                              "controls (tab:support_controls, Section V-G)",
    "common_reference_floor*.csv": "the transfer-and-read floors "
                                   "(tab:floors, Section V-I)",
    "reader_floor.csv": "the reads of the same carrier (tab:floors); its "
                        "psf rows used the mis-scaled read and are unused",
    "reader_rescore.csv": "tab:common_two_readers",
    "level19_rd*.csv": "the level-19 carrier under compression "
                       "(tab:level19_rd, Section V-K)",
    "level19_rd*.json": "tab:level19_rd, and level 19 against level 20 "
                        "under the common reference (block common_reference)",
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
    "kernel_ablation.csv": "the comparison of five kernels, cited by the "
                           "text and kept here (former third block of "
                           "tab:transfer)",
    "resampling_error*.csv": "each transfer against the reference it answers "
                             "to (tab:resampling_error, Appendix C)",
    "cross_term.csv": "the cross term between transfer and compression error "
                      "(appendix on the error budget; table kept here)",
    "interscale*.csv": "the quadtree-native coder test (Appendix A-B)",
    "granule_overhead.json": "the measured tiling overhead (Section IV-C)",
    "granule_overhead.md": "the measured tiling overhead (Section IV-C)",
    "block_size_bias.json": "the odd-side divisibility control (Section IV-C)",
    "bench_environment.json": "the codec library versions and the "
                              "healpix-resample checkout each run used "
                              "(Section VI-C)",
}

# What changed since the previous version of the deposit, written into the
# README. Empty for a first version.
LAYOUT_REUSE = (
    "Version 2 goes with the revised manuscript (healpix-compress tag "
    "`tgrs-submission-v2`). `tiles/` and `windows/` are those of version 1, "
    "byte for byte (their SHA-256 are unchanged in `MANIFEST.csv`).")
LAYOUT_UPDATE = (
    "Version 2 goes with the revised manuscript (healpix-compress tag "
    "`tgrs-submission-v2`) and holds two files. "
    "`healpix-compress-data-tgrs-submission.zip` is the version-1 archive, "
    "unchanged: its `tiles/` and `windows/` are still current. This update "
    "archive, `healpix-compress-data-tgrs-submission-v2.zip`, carries "
    "`regions.csv`, `results/`, this README and a `MANIFEST.csv` that "
    "supersede those of the version-1 archive; the MANIFEST also lists the "
    "tiles and windows, with the SHA-256 they have in the version-1 archive.")
CHANGES = [
    LAYOUT_REUSE,
    "`results/level19_rd.json` gains a block `common_reference`: level 19 "
    "against level 20 when both are scored against the observation. No other "
    "value of the file changed.",
    "`results/reader_rescore.csv` loses the rows of the PSF-aware chain, and "
    "`tab:common_two_readers` its second block.",
    "New: `results/psf_matched_floor.json`, the PSF-aware tiles of four "
    "regions read through their own forward operator.",
    "Withdrawn: `hpx_codec_bench__psf__psf.csv`, "
    "`paper_numbers_common_psf_psf.*`, `common_reference_floor__psf__psf.csv` "
    "and the psf rows of `reader_floor.csv` were scored with a read whose "
    "Gaussian was half as wide as the transfer's (the reader passed the "
    "standard deviation where healpix-resample expects twice it). They stay "
    "in the deposit for the record and are marked so in `MANIFEST.csv`.",
]

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
                     ("describe", ["git", "describe", "--tags", "--match",
                                   "tgrs-*", "--dirty", "--always"])):
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
    sites[keep].to_csv(out, index=False, lineterminator="\n")
    return out


def _tracked_results() -> set[str] | None:
    """Names of the results files git tracks, or None outside a git checkout.

    The results directory also collects files no run of the paper wrote (a
    draft issue, a scratch digest); only what the repository records is
    archived.
    """
    try:
        out = subprocess.run(["git", "ls-files", "--", str(cfg.RESULTS_DIR)],
                             cwd=cfg.REPO_ROOT, check=True, capture_output=True,
                             text=True).stdout
    except Exception:                                         # noqa: BLE001
        return None
    return {Path(line).name for line in out.splitlines()
            if Path(line).parent.name == cfg.RESULTS_DIR.name}


def collect_results(dest: Path) -> list[Path]:
    """Every digest the sweeps and controls wrote, flat, under results/.

    Text files are written with LF line endings whatever the checkout's
    convention, so a file whose content did not change keeps the SHA-256 it
    has in an earlier version of the deposit.
    """
    src = cfg.RESULTS_DIR
    if not src.exists():
        raise SystemExit(f"{src} does not exist; nothing to archive")
    target = dest / "results"
    target.mkdir(parents=True, exist_ok=True)
    tracked = _tracked_results()
    copied = []
    for p in sorted(src.iterdir()):
        if p.is_dir() or p.suffix not in (".csv", ".json", ".md"):
            continue
        if ".pre_extension" in p.name or "DO_NOT_PUBLISH" in p.name:
            continue
        if tracked is not None and p.name not in tracked:
            print(f"  skipped {p.name} (not tracked by git)")
            continue
        (target / p.name).write_bytes(p.read_bytes().replace(b"\r\n", b"\n"))
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


class _HTTPRangeFile(io.RawIOBase):
    """A read-only, seekable file over HTTP range requests (no full download)."""

    def __init__(self, url):
        import requests
        super().__init__()
        self.s = requests.Session()
        h = self.s.head(url, allow_redirects=True, timeout=60)
        h.raise_for_status()
        self.url, self.size, self.pos = h.url, int(h.headers["Content-Length"]), 0

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]
        return self.pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        if n == 0 or self.pos >= self.size:
            return b""
        end = min(self.pos + n, self.size) - 1
        r = self.s.get(self.url, headers={"Range": f"bytes={self.pos}-{end}"},
                       timeout=600)
        r.raise_for_status()
        self.pos += len(r.content)
        return r.content

    def readinto(self, b):
        data = self.read(len(b))
        b[:len(data)] = data
        return len(data)


def _open_archive(source: str):
    import zipfile
    if source.startswith("http"):
        raw = io.BufferedReader(_HTTPRangeFile(source), buffer_size=8 << 20)
    else:
        raw = open(source, "rb")
    z = zipfile.ZipFile(raw)
    return z, z.namelist()[0].split("/")[0]


def previous_manifest_rows(source: str) -> list[dict]:
    """The tiles/ and windows/ rows of an earlier archive's MANIFEST.csv."""
    import csv
    z, root = _open_archive(source)
    with z.open(f"{root}/MANIFEST.csv") as m:
        rows = [r for r in csv.DictReader(io.TextIOWrapper(m, encoding="utf-8"))
                if r["path"].startswith(("tiles/", "windows/"))]
    name = source.split("/")[-1].split("?")[0]
    for r in rows:
        r["bytes"] = int(r["bytes"])
        r["feeds"] = (r.get("feeds") or "") + f" (in {name}, unchanged)"
    return rows


def reuse_from(source: str, dest: Path) -> list[Path]:
    """Copy tiles/ and windows/ from an earlier archive, checked by SHA-256.

    `source` is a local .zip or an URL; members are read one at a time and
    each one must match the SHA-256 the earlier archive's MANIFEST.csv gives.
    """
    import csv
    import io
    import zipfile
    if source.startswith("http"):
        raw = io.BufferedReader(_HTTPRangeFile(source), buffer_size=8 << 20)
    else:
        raw = open(source, "rb")
    z = zipfile.ZipFile(raw)
    root = z.namelist()[0].split("/")[0]
    with z.open(f"{root}/MANIFEST.csv") as m:
        sums = {r["path"]: r["sha256"]
                for r in csv.DictReader(io.TextIOWrapper(m, encoding="utf-8"))}
    copied = []
    for path in sorted(sums):
        if not (path.startswith("tiles/") or path.startswith("windows/")):
            continue
        out = dest / path
        out.parent.mkdir(parents=True, exist_ok=True)
        with z.open(f"{root}/{path}") as src, out.open("wb") as dst:
            shutil.copyfileobj(src, dst, 1 << 20)
        if sha256(out) != sums[path]:
            raise SystemExit(f"{path}: SHA-256 differs from {source}'s MANIFEST")
        copied.append(out)
        print(f"  reused {path}")
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
    ]
    if CHANGES:
        lines += ["", "## Changes from the previous version", ""]
        lines += [f"- {c}" for c in CHANGES]
    lines += ["",
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
        "related_identifiers": [
            {"identifier": "10.5281/zenodo.23812788", "relation": "isCompiledBy",
             "resource_type": "software"},
            {"identifier": "10.5281/zenodo.22210945", "relation": "isDerivedFrom",
             "resource_type": "dataset"},
        ],
        "notes_for_the_author": [
            "This is a NEW VERSION of record 10.5281/zenodo.23388572 (concept "
            "10.5281/zenodo.23388571): open that record and click 'New "
            "version', do not create a new upload.",
            "Keep the related identifiers of version 1 (isCompiledBy the code, "
            "isDerivedFrom the Copernicus Sentinel-2 products) and the ESA "
            "funding (contract 4000147951/25/I-NS); the code is cited by its "
            "concept DOI 10.5281/zenodo.23812788.",
        ],
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", default="tgrs-submission-v2",
                    help="suffix of the archive directory (default: "
                         "tgrs-submission-v2)")
    ap.add_argument("--reuse", default=None,
                    help="take tiles/ and windows/ from this earlier archive "
                         "(local .zip or URL), checked against its MANIFEST")
    ap.add_argument("--update-of", default=None,
                    help="build an update archive: list the tiles/ and "
                         "windows/ of this earlier archive (local .zip or "
                         "URL) in MANIFEST.csv without copying them")
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
    previous_rows = []
    if args.update_of:
        if args.tiles or args.windows or args.reuse:
            raise SystemExit("--update-of excludes --tiles, --windows and --reuse")
        previous_rows = previous_manifest_rows(args.update_of)
        counts["tiles"] = sum(1 for r in previous_rows if r["path"].startswith("tiles/"))
        counts["windows"] = sum(1 for r in previous_rows if r["path"].startswith("windows/"))
        CHANGES[0] = LAYOUT_UPDATE
    if args.reuse:
        if args.tiles or args.windows:
            raise SystemExit("--reuse replaces --tiles and --windows")
        r = reuse_from(args.reuse, dest)
        counts["tiles"] = sum(1 for x in r if x.parent.name == "tiles")
        counts["windows"] = sum(1 for x in r if x.parent.name == "windows")
        files += r

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
    rows += previous_rows
    man = pd.DataFrame(rows).sort_values("path")
    man.to_csv(dest / "MANIFEST.csv", index=False)

    (dest / "README.md").write_text(readme(dest, version, counts),
                                    encoding="utf-8")
    (dest / "zenodo_data.json").write_text(
        json.dumps(zenodo_metadata(version), indent=2), encoding="utf-8")

    total = sum(int(r["bytes"]) for r in rows if r not in previous_rows)
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
