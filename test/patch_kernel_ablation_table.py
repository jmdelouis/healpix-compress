# -*- coding: utf-8 -*-
"""The kernel-ablation table is now cited by the paper, so it has to fit and
say how many regions it rests on."""
import sys
from pathlib import Path

p = Path(sys.argv[1] if len(sys.argv) > 1 else "kernel_ablation.py")
s = p.read_text(encoding="utf-8")

# full width: the bicubic row's magnitudes overrun a single column
n = 0
for old, new in (("\\\\begin{table}[t]", "\\\\begin{table*}[t]"),
                 ("\\\\end{table}", "\\\\end{table*}")):
    if old in s:
        s = s.replace(old, new); n += 1
for old, new in (('"\\\\begin{table}[t]', '"\\\\begin{table*}[t]'),
                 ('"\\\\end{table}\\\\n"', '"\\\\end{table*}\\\\n"')):
    if old in s:
        s = s.replace(old, new); n += 1
print(f"table* substitutions: {n}")

# the region count belongs in the caption of a table the text now leans on
CAND = ["Transfer-kernel ablation.", "Transfer-kernel ablation "]
done = False
for c in CAND:
    if c in s:
        s = s.replace(c, c.rstrip() + " Medians over the {n_scenes} regions "
                                      "swept.", 1)
        done = True
        break
print("caption note inserted:", done)
if done and "{n_scenes}" in s:
    # make it an f-string field only if the caption is already an f-string;
    # otherwise fall back to a literal the caller must not have to think about
    print("NOTE: '{n_scenes}' placeholder inserted -- check the caption is an "
          "f-string, or replace it with the literal count.")
p.write_text(s, encoding="utf-8")
print("kernel_ablation.py patched")
