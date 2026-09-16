"""Remove the explanatory note from under every generated table, once.

Why
---
Each generated table carried a paragraph of prose under it, emitted by
`paper_config.tex_note` as a minipage. IEEE's style manual asks for a short
title above the table, centred, in small caps and set as an inverted pyramid;
it does not provide for a paragraph below one. Forty-two of them also made the
manuscript impossible to size at a glance.

`tex_note` now returns nothing, so any table regenerated from here on comes out
clean. This script applies the same change to the tables already on disk, so
that the manuscript can be compiled without re-running forty-two analyses to
get there. It is idempotent, and it is meant to be run once.

The note text is not thrown away: each one is written to
`<tables>/notes/<table>.txt`, because it is the raw material for the prose that
has to introduce the table instead, and because it records in one place what
each column of that table means.

Usage
-----
    python strip_table_notes.py --dry-run
    python strip_table_notes.py
"""

from __future__ import annotations

import argparse
import pathlib
import re

import paper_config as cfg

# What tex_note used to emit: a vspace, a minipage at column or text width,
# \footnotesize, the note on its own line, and the close.
NOTE = re.compile(
    r"\n?\\vspace\{3pt\}\n"
    r"\\begin\{minipage\}\{[^}]*\}\n"
    r"\\footnotesize\n"
    r"(?P<text>.*?)\n"
    r"\\end\{minipage\}\n",
    re.S)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--tables", default=None,
                    help=f"default: {cfg.TEX_TABLES}")
    args = ap.parse_args(argv)

    root = pathlib.Path(args.tables) if args.tables else cfg.TEX_TABLES
    if not root.is_dir():
        raise SystemExit(f"{root} is not a directory")
    notes = root / "notes"
    stripped = words = 0
    for p in sorted(root.glob("*.tex")):
        s = p.read_text(encoding="utf-8")
        m = NOTE.search(s)
        if not m:
            continue
        text = m.group("text")
        words += len(text.split())
        stripped += 1
        print(f"  {p.name}: {len(text.split())} words")
        if args.dry_run:
            continue
        notes.mkdir(parents=True, exist_ok=True)
        (notes / f"{p.stem}.txt").write_text(text, encoding="utf-8")
        p.write_text(s[:m.start()] + "\n" + s[m.end():], encoding="utf-8")
    verb = "would strip" if args.dry_run else "stripped"
    print(f"\n{verb} {stripped} note(s), {words} words")
    if stripped and not args.dry_run:
        print(f"note text kept under {notes}")
    if not stripped:
        print("nothing to do (already clean, or the tables predate tex_note)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
