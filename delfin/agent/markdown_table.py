"""Recognising a GitHub-style pipe table, for every surface that draws one.

Input: the lines of an answer. Output: whether a table starts at a line,
its column alignment, and its rows split into cells -- no drawing. The
dashboard turns the result into HTML (``delfin.dashboard.chat_markdown``),
the terminal into aligned columns (``repl_render.MarkdownStream``); both
read the same recognition so a table is a table on both or on neither.

Stdlib only: the terminal imports it, and must not pay for the dashboard.
"""

from __future__ import annotations

import re

_DELIM_CELL = re.compile(r"^:?-+:?$")
_NUMERIC = re.compile(
    r"^[~≈±+\-−]?\s?\d[\d.,'  ]*(?:[eE][+\-−]?\d+)?\s?%?$")


def split_row(line: str) -> list[str]:
    """The cells of a pipe-table row, inline markup left in place.

    A pipe splits unless it is escaped (``\\|``) or inside a code span whose
    closing backticks are on the same line. Inside a code span an escaped
    pipe becomes a plain one, as GitHub renders it.
    """
    s = line.strip()
    if s.startswith("|"):
        s = s[1:]
    if s.endswith("|") and not s.endswith("\\|"):
        s = s[:-1]
    cells: list[str] = []
    buf: list[str] = []
    i, n, code_run = 0, len(s), 0
    while i < n:
        ch = s[i]
        if ch == "\\" and i + 1 < n and s[i + 1] == "|":
            buf.append("|" if code_run else "\\|")
            i += 2
            continue
        if ch == "`":
            j = i
            while j < n and s[j] == "`":
                j += 1
            run = j - i
            if code_run == 0:
                if s.find("`" * run, j) >= 0:
                    code_run = run
            elif run == code_run:
                code_run = 0
            buf.append(s[i:j])
            i = j
            continue
        if ch == "|" and code_run == 0:
            cells.append("".join(buf))
            buf = []
            i += 1
            continue
        buf.append(ch)
        i += 1
    cells.append("".join(buf))
    return [c.strip() for c in cells]


def has_pipe(line: str) -> bool:
    """Whether *line* holds an unescaped pipe."""
    return bool(re.search(r"(?<!\\)\|", line))


def alignment(header: str, delim: str) -> list[str] | None:
    """Column alignments ("l", "c", "r" or "") when *header* followed by
    *delim* opens a table, else None.

    The alignment row must have as many cells as the header: that is what
    tells ``a | b`` over ``---`` (a sentence and a rule) from a table.
    """
    if not has_pipe(header):
        return None
    if "-" not in delim or not re.fullmatch(r"[\s|:\-]+", delim):
        return None
    dcells = split_row(delim)
    if not dcells or not all(_DELIM_CELL.match(c) for c in dcells):
        return None
    if len(dcells) != len(split_row(header)):
        return None
    out = []
    for c in dcells:
        if c.startswith(":") and c.endswith(":"):
            out.append("c")
        elif c.endswith(":"):
            out.append("r")
        elif c.startswith(":"):
            out.append("l")
        else:
            out.append("")
    return out


def is_numeric(cell: str) -> bool:
    """A cell holding one number (sign, digits, separators, exponent, %)."""
    return bool(_NUMERIC.match(cell.replace("`", "").strip()))


def normalise(header: list[str], rows: list[list[str]], aligns: list[str]
              ) -> tuple[list[str], list[list[str]], list[str]]:
    """Every row as wide as the widest, and a default alignment.

    Short rows are padded with empty cells; extra cells are kept rather
    than dropped (GitHub drops them, which loses what the model wrote). A
    column the model left unaligned but that holds only numbers is
    right-aligned, so the decimal places stand under each other.
    """
    ncols = max([len(header)] + [len(r) for r in rows])
    header = header + [""] * (ncols - len(header))
    rows = [r + [""] * (ncols - len(r)) for r in rows]
    aligns = aligns + [""] * (ncols - len(aligns))
    for c in range(ncols):
        if not aligns[c]:
            body = [r[c] for r in rows if r[c]]
            if body and all(is_numeric(v) for v in body):
                aligns[c] = "r"
    return header, rows, aligns
