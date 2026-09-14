"""Verify that the benchmark chain is completely and consistently installed.

The pipeline is split across several sibling modules in this directory, and a
partial sync between machines (laptop -> mounted folder -> JupyterHub) shows up
much later as a bare `ModuleNotFoundError` in the middle of a run. This checks
the whole set at once, and prints a size/sha256 per file so two machines can be
compared with a single command.

    python check_install.py            # files, imports, third-party deps
    python check_install.py --sha      # sha256 only, to diff across machines
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

MODULES = [
    ("paper_config", "geometry and protocol constants"),
    ("region_sites", "the 40 shared regions, read from healpix-resample"),
    ("window_geometry", "which tile an anchor selects, and what window it needs"),
    ("build_region_windows", "acquisition of the native windows"),
    ("build_tiles", "HEALPix tile + equal-area UTM crop"),
    ("codec_lib", "codec adapters and sweeps"),
    ("common_reference", "forward operator for the common UTM reference"),
    ("run_codec_bench", "rate-distortion sweep over both grids"),
    ("paper_stats", "interpolation, oracle codec, Wilcoxon, Holm"),
    ("update_paper_numbers", "tables and digest"),
    ("make_figures", "paper figures"),
    ("codec_audit", "byte/type audit of the codec adapters"),
    ("support_controls", "controls for the footprint support mismatch"),
    ("gamma_analysis", "verdict over a continuous overhead ratio"),
    ("superres_cost", "storage surcharge of a deconvolving transfer"),
    ("cross_term", "is the compression error orthogonal to the transfer error"),
    ("level_choice", "is level 20 the right carrier level"),
    ("interscale", "does the link between scales pay"),
    ("reader_floor", "how much of the floor is the read"),
    ("reader_rescore", "the common-reference verdict under two readers"),
    ("check_paper_numbers", "manuscript figures vs the run that produced them"),
    ("operator_tests", "what the HEALPix->UTM operator does beyond a constant"),
    ("kernel_ablation", "transfer-kernel ablation"),
    ("resampling_error", "reconstruction error of each transfer"),
    ("alpha_budget", "the verdict in units of the product's own error"),
]

THIRD_PARTY = [
    ("numpy", True), ("scipy", True), ("pandas", True), ("matplotlib", True),
    ("xarray", True), ("zarr", True), ("healpix_geo", True),
    ("healpix_resample", True), ("hcc", True),
    ("pyproj", True), ("rasterio", False), ("pystac_client", False),
    ("zfpy", False), ("pysz", False), ("imagecodecs", False),
    ("PIL", False), ("pywt", False), ("zstandard", False),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sha", action="store_true",
                    help="print only 'sha256  size  name' lines, for diffing")
    args = ap.parse_args(argv)

    if args.sha:
        for name, _ in MODULES:
            p = HERE / f"{name}.py"
            print(f"{sha256(p) if p.exists() else '-'*64}  "
                  f"{p.stat().st_size if p.exists() else 0:>7}  {name}.py")
        return 0

    print(f"directory: {HERE}\npython:    {sys.version.split()[0]}\n")
    missing, broken = [], []
    print("pipeline modules")
    for name, what in MODULES:
        path = HERE / f"{name}.py"
        if not path.exists():
            missing.append(name)
            print(f"  MISSING  {name}.py  -- {what}")
            continue
        try:
            importlib.import_module(name)
            status = "ok      "
        except Exception as exc:                          # noqa: BLE001
            broken.append((name, exc))
            status = "IMPORT !"
        print(f"  {status} {name}.py  {path.stat().st_size:>6} B  "
              f"{sha256(path)[:12]}  -- {what}")

    print("\nthird-party packages")
    hard = []
    for mod, required in THIRD_PARTY:
        try:
            importlib.import_module(mod)
            print(f"  ok       {mod}")
        except Exception:                                 # noqa: BLE001
            print(f"  {'MISSING ' if required else 'absent  '} {mod}"
                  f"{'' if required else '  (only needed for some steps)'}")
            if required:
                hard.append(mod)

    print()
    if missing:
        print(f"{len(missing)} pipeline file(s) absent from this checkout: "
              + ", ".join(f"{m}.py" for m in missing))
        print("  This is a sync problem, not a code problem. Commit the files "
              "on the machine that has them and pull here:")
        print("    git add test/*.py *.md && git commit -m 'paper pipeline' "
              "&& git push")
        print("  then compare both machines with:  python check_install.py --sha")
    for name, exc in broken:
        print(f"{name}.py present but fails to import: {type(exc).__name__}: {exc}")
    if hard:
        print(f"required packages missing: {', '.join(hard)}  "
              f"(pip install -e '.[paper]' from the repository root)")
    if not (missing or broken or hard):
        print("Everything the pipeline needs is present. "
              "Next: python run_codec_bench.py --limit 2 --grid both")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
