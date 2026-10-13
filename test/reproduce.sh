#!/usr/bin/env bash
# Rebuild everything the paper quotes, from the native windows on disk to the
# compiled manuscript. Stops at the first failing step (set -e), so a partial
# run leaves results/ in a state the next run can resume from.
#
# Prerequisites (REPRODUCE.md, steps 1 and 2.2): the environment is installed,
# HEALPIX_RESAMPLE_ROOT points to the healpix-resample checkout, and the native
# windows are in test/data/regions/ (build_region_windows.py, the only step
# that needs the network, is NOT run here).
#
# Usage, from test/:
#     bash reproduce.sh                       # the 40 regions, everything
#     REGIONS="paris amazon" bash reproduce.sh   # smoke test on two regions
#     SKIP_TEX=1 bash reproduce.sh            # stop before pdflatex
#
# Every step prints its name and its wall time; the full log is written to
# results/reproduce.log.

set -euo pipefail
cd "$(dirname "$0")"

REGIONS="${REGIONS:-}"          # empty = the 40 regions
SKIP_TEX="${SKIP_TEX:-0}"
LOG="results/reproduce.log"
mkdir -p results

# --regions ... is accepted by the per-region scripts; the digest and table
# scripts take no region argument and work on whatever is in results/.
R=()
if [[ -n "$REGIONS" ]]; then
    # shellcheck disable=SC2206
    R=(--regions $REGIONS)
fi

step() {                        # step "<label>" <command...>
    local label="$1"; shift
    local t0
    t0=$(date +%s)
    printf '\n==== %s\n     %s\n' "$label" "$*" | tee -a "$LOG"
    if "$@" 2>&1 | tee -a "$LOG"; then
        printf '     ok, %ss\n' "$(( $(date +%s) - t0 ))" | tee -a "$LOG"
    else
        printf '\n**** FAILED after %ss: %s\n' "$(( $(date +%s) - t0 ))" "$label" | tee -a "$LOG"
        exit 1
    fi
}

printf '\n######## reproduce.sh  %s  regions=%s\n' "$(date -Is)" "${REGIONS:-all}" | tee -a "$LOG"

# ---- 0. environment -------------------------------------------------------
step "check_install"           python check_install.py
step "paper_config"            python paper_config.py
step "region_sites"            python region_sites.py
step "windows on disk"         python build_region_windows.py --check "${R[@]}"

# ---- 2.3 tiles and crops --------------------------------------------------
step "tiles: clough_tocher"    python build_tiles.py "${R[@]}"
step "tiles: psf"              python build_tiles.py --resampler psf "${R[@]}"
step "tiles: linear"           python build_tiles.py --resampler linear "${R[@]}"
step "tiles: audit"            python build_tiles.py --audit "${R[@]}"

# ---- 2.4 dataset audits ---------------------------------------------------
step "window_radiometry"       python window_radiometry.py "${R[@]}"
step "tile_crop_consistency"   python tile_crop_consistency.py "${R[@]}"
step "codec_audit"             python codec_audit.py --grid both "${R[@]}"
step "operator_tests"          python operator_tests.py "${R[@]}"

# ---- 3. sweeps ------------------------------------------------------------
step "bench: both grids"       python run_codec_bench.py --grid both "${R[@]}"
step "bench: extend ladders"   python extend_sweep.py "${R[@]}"
step "bench: common ref, nearest"  python run_codec_bench.py --grid hpx --common-reference --force "${R[@]}"
step "bench: common ref, bilinear" python run_codec_bench.py --grid hpx --common-reference --cr-mode bilinear --force "${R[@]}"
step "bench: psf transfer"     python run_codec_bench.py --grid hpx --resampler psf "${R[@]}"
step "bench: linear transfer"  python run_codec_bench.py --grid hpx --resampler linear "${R[@]}"

step "resampling_error"        python resampling_error.py "${R[@]}"
step "resampling_error noise0" python resampling_error.py --noise 0 "${R[@]}"
step "superres_cost"           python superres_cost.py --kernel psf --control linear --tex
step "kernel_ablation"         python kernel_ablation.py --rd --limit 0 "${R[@]}"

step "level19_rd"              python level19_rd.py "${R[@]}"
step "level19_direct"          python level19_direct.py "${R[@]}"
step "level_choice"            python level_choice.py "${R[@]}"

step "controls: synthetic"     python support_controls.py --control synthetic "${R[@]}"
step "controls: synthetic via transfer" python support_controls.py --control synthetic --transfer clough_tocher "${R[@]}"
step "controls: smoothed 1/2"  python support_controls.py --control smoothed "${R[@]}"
step "controls: smoothed 1/4"  python support_controls.py --control smoothed --smooth-shift 0.25 0.25 "${R[@]}"
step "controls: common, shift" python support_controls.py --control common shift "${R[@]}"
step "block_size_bias"         python block_size_bias.py --limit 0 "${R[@]}"
step "granule_overhead"        python granule_overhead.py

step "reader_floor"            python reader_floor.py "${R[@]}"
step "alpha_budget"            python alpha_budget.py --no-publish
step "cross_term"              python cross_term.py "${R[@]}"
step "diag_bicubic"            python diag_bicubic.py --all
step "interscale"              python interscale.py "${R[@]}"

# ---- 4. statistics and digests -------------------------------------------
step "paper_numbers self"      python update_paper_numbers.py
step "paper_numbers common"    python update_paper_numbers.py --metric common
step "paper_numbers common/bilinear" python update_paper_numbers.py --metric common --cr-mode bilinear
step "censoring_sensitivity"   python censoring_sensitivity.py
step "gamma_analysis"          python gamma_analysis.py
step "headline"                python headline.py
step "reader_rescore"          python reader_rescore.py

# ---- 5. tables and figures ------------------------------------------------
step "tables: synthesis"       python merged_tables.py --only synthesis
step "tables: floors"          python merged_tables.py --only floors
step "tables: transfer"        python merged_tables.py --only transfer
step "tables: rebuild"         python rebuild_tables.py
step "fig_workflow"            python fig_workflow.py
step "make_figures"            python make_figures.py

# ---- 6. checker and manuscript -------------------------------------------
step "glossary terms"          python bold_glossary_terms.py ../tex/sections/*.tex
step "check_paper_numbers"     python check_paper_numbers.py

if [[ "$SKIP_TEX" != "1" ]]; then
    cd ../tex
    step "pdflatex 1"          pdflatex -interaction=nonstopmode tgrs_healpix_compress.tex
    step "bibtex"              bibtex tgrs_healpix_compress
    step "pdflatex 2"          pdflatex -interaction=nonstopmode tgrs_healpix_compress.tex
    step "pdflatex 3"          pdflatex -interaction=nonstopmode tgrs_healpix_compress.tex
    if grep -q "undefined" tgrs_healpix_compress.log; then
        echo "**** undefined references in tgrs_healpix_compress.log" | tee -a "../test/$LOG"
        exit 1
    fi
    grep -c "Overfull" tgrs_healpix_compress.log | sed 's/^/     overfull boxes: /' || true
    cd ../test
fi

printf '\n######## done  %s\n' "$(date -Is)" | tee -a "$LOG"
