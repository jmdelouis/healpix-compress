"""Set the glossary terms of Table II in italics at every occurrence (revision 7).

Wraps each term in \term{...} (defined in the main file as \emph; tile and
overhead are left plain at the author's request) in the
plain-text runs of the section files only: never inside math, inside the
argument of a command (\ref, \cite, \emph, \textbf, \caption, ...), inside
a table environment or in a section title. Idempotent: an already wrapped
occurrence is left alone.
"""
import re, sys, pathlib

TERMS = [
    # (regex, replacement) applied in this order; longer phrases first
    (r"\b(native UTM squares?|equal-area UTM squares?|equal-area squares?|native squares?|UTM squares?|equal-area crops?|native crops?|UTM crops?)\b(?!-)", r"\\term{\1}"),
    (r"\b(carrier grids?)\b", r"\\term{\1}"),
    (r"\b(critical overheads?)\b", r"\\term{\1}"),
    (r"\b(size ratios?)\b", r"\\term{\1}"),
    (r"\b(target RMSE|RMSE targets?)\b", r"\\term{\1}"),
    (r"\b(matched(?:-| )errors?|matched(?:-| )budgets?|matched(?:-| )storage|matched(?:-| )RMSE)\b", r"\\term{\1}"),
    (r"\b(common(?:-| )references?|reconstruction references?)\b", r"\\term{\1}"),
    (r"\b(carriers?)\b(?!-)", r"\\term{\1}"),
    (r"\b(transfers?)\b(?!-and-read)(?!-induced)(?!-side)", r"\\term{\1}"),
    (r"\b(nearest|bilinear|interpolating|containing-cell|matched|same|the|The|a|A|each|Each|that|this|either|every|both|one|its|whose|which|two|per|mismatched) (reads?)\b(?!\s+(?:through|back|on|it|them|at|by|as|the|a|each|all|in|is|are|was|were|and|only|directly|identically|with|from|patterns?|pixel)\b)",
     r"\1 \\term{\2}"),
    (r"(?<!nearest-cell )(?<!nearest )(?<!bilinear )\b(reading operator)\b", r"\\term{read}ing operator"),
    (r"(?<!to )\b(targets?)\b(?!-)(?!\s+(?:a |the product|on the product))", r"\\term{\1}"),  # noqa
    (r"(?<!shared )(?<!cross-)(?<!spatial )(?<!coordinate )(?<!single )(?<!shared\n)\b(references?)\b(?!-)(?! frame)(?!, the raw size)", r"\\term{\1}"),
]

SKIP_ENV = ("table", "table*", "figure", "figure*", "equation", "align")
CMD_WITH_ARG = re.compile(r"\\[a-zA-Z@]+\*?\s*(?:\[[^\]]*\])?\s*\{")


def split_runs(text):
    """Yield (is_text, chunk) preserving order; text chunks get the regexes."""
    out = []
    i = 0
    n = len(text)
    env_depth = 0
    buf = []
    while i < n:
        c = text[i]
        if c == "%":                       # comment to end of line
            j = text.find("\n", i)
            j = n if j < 0 else j
            out.append((False, text[i:j])); i = j; continue
        if c == "$":                       # inline math
            j = text.find("$", i + 1)
            j = n if j < 0 else j + 1
            out.append((False, text[i:j])); i = j; continue
        if c == "\\":
            m = CMD_WITH_ARG.match(text, i)
            if m:
                name = re.match(r"\\([a-zA-Z@]+)", text[i:]).group(1)
                # find the matching closing brace of the first argument
                j = m.end(); depth = 1
                while j < n and depth:
                    if text[j] == "{": depth += 1
                    elif text[j] == "}": depth -= 1
                    j += 1
                if name in ("begin",):
                    env = text[m.end():j - 1]
                    if env in SKIP_ENV:
                        k = text.find("\\end{" + env + "}", j)
                        k = n if k < 0 else k + len("\\end{" + env + "}")
                        out.append((False, text[i:k])); i = k; continue
                if name in ("section", "subsection", "subsubsection", "paragraph",
                            "caption", "label", "ref", "cite", "eqref", "emph",
                            "textbf", "texttt", "term", "input", "SI", "num",
                            "todo", "url", "footnote", "thanks"):
                    out.append((False, text[i:j])); i = j; continue
                # other commands: keep the command token opaque, recurse into nothing
                out.append((False, text[i:m.end()])); i = m.end(); continue
            m2 = re.match(r"\\[a-zA-Z@]+\*?|\\.", text[i:])
            out.append((False, m2.group(0))); i += len(m2.group(0)); continue
        # plain text run
        j = i
        while j < n and text[j] not in "\\$%":
            j += 1
        out.append((True, text[i:j])); i = j
    return out


def bold(text):
    runs = split_runs(text)
    res = []
    for is_text, chunk in runs:
        if is_text:
            store = []
            def keep(m):
                store.append(m.group(0))
                return "\x00%d\x00" % (len(store) - 1)
            for pat, rep in TERMS:
                # wrap, then hide the wrapped text from the later patterns
                chunk = re.sub(pat, rep, chunk)
                chunk = re.sub(r"\\term\{[^{}]*\}(?:ing operator)?", keep, chunk)
            for word in ("Carrier", "Transfer", "Reference", "Target"):
                chunk = re.sub(r"\b(%s)\b(?!-)" % word, r"\\term{\1}", chunk)
                chunk = re.sub(r"\\term\{[^{}]*\}", keep, chunk)
            chunk = re.sub(r"\x00(\d+)\x00", lambda m: store[int(m.group(1))], chunk)
        res.append(chunk)
    return "".join(res)


if __name__ == "__main__":
    for f in sys.argv[1:]:
        p = pathlib.Path(f)
        s = p.read_text(encoding="utf-8")
        # idempotence: strip existing wrappers first
        s = re.sub(r"\\term\{([^{}]*)\}", r"\1", s)
        s = s.replace("\\term{read}ing operator", "reading operator")
        t = bold(s)
        p.write_text(t, encoding="utf-8")
        print(f, t.count("\\term{"))
