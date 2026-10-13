#!/usr/bin/env bash
# List the files under tex/ that changed since the last upload to Overleaf.
#
# Overleaf is the official copy of the manuscript and this repository is its
# mirror; without Overleaf's git features the sync is a manual upload. This
# prints what to upload (status: M modified, A added, D deleted in Overleaf).
# After uploading, mark the point the two copies agree on:
#
#     bash tex/overleaf_changes.sh        # what to upload
#     git tag -f overleaf-synced          # once Overleaf matches this commit
#
# Build products (the PDF, .aux, .log, ...) are never uploaded.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
if git rev-parse -q --verify overleaf-synced >/dev/null; then
    ref=overleaf-synced
else
    echo "no 'overleaf-synced' tag yet: showing changes of the last commit; set the tag after the first upload" >&2
    ref=HEAD~1
fi
echo "changes in tex/ since $ref ($(git log -1 --format="%h %as" "$ref")):"
{ git diff --name-status "$ref" -- tex/; git status --short -- tex/ | sed 's/^ *//'; } \
  | grep -Ev '\.(pdf|aux|log|bbl|blg|out|synctex\.gz|fdb_latexmk|fls)$|tgrs_healpix_compress\.pdf' \
  | grep -Ev '^\?\? ' \
  | sort -u -k2 || true
echo "(figures/*.pdf are listed separately below: they are inputs, not build products)"
{ git diff --name-status "$ref" -- tex/figures/; git status --short -- tex/figures/ | sed 's/^ *//'; } | sort -u -k2 || true
