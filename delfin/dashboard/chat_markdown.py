"""Markdown to HTML for the agent chat.

Input: the text of a chat message as a model wrote it. Output: an HTML
fragment for an ``ipywidgets.HTML`` value, which the browser assigns with
``innerHTML`` -- nothing between this function and the page sanitises it.

Semantics, block level: GitHub-style pipe tables (alignment row, escaped
``\\|``, a pipe inside a code span stays in its cell, a table may follow a
sentence without a blank line), fenced code with ``` or ~~~ and any fence
length (a longer fence holds a shorter one; an unclosed fence runs to the
end), ATX headings 1-6, nested bullet and ordered lists (content of an item
is rendered recursively, so a table or a fence under an item stays one),
block quotes, horizontal rules, ``$$`` blocks shown as their source, and
paragraphs whose single newlines are kept as ``<br>`` -- models write line
breaks they mean.

Inline: code spans of any backtick length, backslash escapes, ``**strong**``,
``*em*``, ``~~del~~``, links and autolinks to ``http(s)``/``mailto`` only,
bare web addresses. ``$...$`` is kept verbatim so a product inside it is
not read as emphasis. Underscores are never emphasis: ``calc_a_b`` and
``__init__.py`` occur in answers far more often than ``_word_`` does.

Safety: every character of model text reaches the output either
HTML-escaped or not at all; the only markup is what this module emits,
with constant class names. An ``href`` is emitted only for an allow-listed
scheme and is attribute-escaped. NUL characters are removed first because
they delimit the internal placeholders. Fenced code is handed to the
caller's ``code_block(lang, code)``, which owns its own escaping.

Rejected: mistune (present only as a transitive dependency of nbconvert,
with a different API in 2.x and 3.x, and its escape/sanitise behaviour
would then be configuration rather than code read here); line-by-line
regex substitution, which is what failed -- a table, a list or a fence is
a block, and only a block parser can see where it ends.
"""

from __future__ import annotations

import html as _html
import re
from typing import Callable

CodeBlockFn = Callable[[str, str], str]

# -- block patterns ----------------------------------------------------------

_FENCE_OPEN = re.compile(r"^( {0,3})(`{3,}|~{3,})[ \t]*([^\s`]*)[^`]*$")
_HEADING = re.compile(r"^ {0,3}(#{1,6})(?:[ \t]+(.*?))?[ \t]*$")
_HR = re.compile(r"^ {0,3}(?:(?:-[ \t]*){3,}|(?:\*[ \t]*){3,}|(?:_[ \t]*){3,})$")
_QUOTE = re.compile(r"^ {0,3}> ?(.*)$")
_ITEM = re.compile(r"^( *)([-*+]|\d{1,9}[.)])(?:([ \t]+)(.*))?$")
_DELIM_CELL = re.compile(r"^:?-+:?$")
_LANG = re.compile(r"[\w+#.\-]{1,32}")
_NUMERIC = re.compile(
    r"^[~≈±+\-−]?\s?\d[\d.,'\u2009\u202f]*(?:[eE][+\-−]?\d+)?\s?%?$")

# -- inline patterns ---------------------------------------------------------

_CODE_SPAN = re.compile(r"(`+)(.+?)(?<!`)\1(?!`)")
_ESCAPE = re.compile(r"\\([!-/:-@\[-`{-~])")
_MATH = re.compile(r"\$\$[^$\n]+?\$\$|\$[^\s$](?:[^$\n]*?[^\s$])?\$")
_AUTOLINK = re.compile(r"<((?:https?://|mailto:)[^\s<>\"']+)>")
_LINK = re.compile(
    r"!?\[([^\]\n]+)\]\(\s*([^)\s]+)(?:\s+\"[^\"\n]*\")?\s*\)")
_BARE_URL = re.compile(
    r"(?<![\w/\"'=@])(https?://[^\s<>\"'`]*[^\s<>\"'`.,;:!?)\]])")
_STRONG_EM = re.compile(r"\*\*\*(?=\S)(.+?)(?<=\S)\*\*\*")
_STRONG = re.compile(r"\*\*(?=\S)(.+?)(?<=\S)\*\*")
_EM = re.compile(r"(?<![\w*])\*(?=[^\s*])(.+?)(?<=[^\s*])\*(?![\w*])")
_DEL = re.compile(r"~~(?=\S)(.+?)(?<=\S)~~")
_SLOT = re.compile(r"\x00(\d+)\x00")
_SAFE_SCHEME = re.compile(r"^(?:https?://|mailto:)", re.IGNORECASE)


def render(text: str, code_block: CodeBlockFn) -> str:
    """The HTML fragment for *text*; see the module docstring."""
    src = str(text or "").replace("\x00", "")
    src = src.replace("\r\n", "\n").replace("\r", "\n").expandtabs(4)
    # Leading and trailing blank lines would only pad the bubble.
    src = src.strip("\n \t")
    if not src:
        return ""
    return _blocks(src.split("\n"), code_block)


# -- inline ------------------------------------------------------------------

def _link(href: str, label_html: str) -> str:
    return (f'<a href="{_html.escape(href, quote=True)}" target="_blank" '
            f'rel="noopener noreferrer">{label_html}</a>')


def inline(text: str) -> str:
    """Inline markdown of one line, escaped."""
    slots: list[str] = []

    def put(fragment: str) -> str:
        slots.append(fragment)
        return f"\x00{len(slots) - 1}\x00"

    def code(m: re.Match) -> str:
        body = m.group(2)
        if len(body) > 1 and body[0] == " " and body[-1] == " ":
            body = body[1:-1]
        return put("<code>" + _html.escape(body, quote=False) + "</code>")

    def link(m: re.Match) -> str:
        label, href = m.group(1), m.group(2)
        if not _SAFE_SCHEME.match(href):
            return m.group(0)               # shown as written, escaped below
        return put(_link(href, inline(label)))

    s = text.replace("\x00", "")
    s = _CODE_SPAN.sub(code, s)
    s = _ESCAPE.sub(lambda m: put(_html.escape(m.group(1), quote=False)), s)
    s = _MATH.sub(lambda m: put(_html.escape(m.group(0), quote=False)), s)
    s = _AUTOLINK.sub(
        lambda m: put(_link(m.group(1), _html.escape(m.group(1),
                                                      quote=False))), s)
    s = _LINK.sub(link, s)
    s = _BARE_URL.sub(
        lambda m: put(_link(m.group(1), _html.escape(m.group(1),
                                                      quote=False))), s)
    s = _html.escape(s, quote=False)
    s = _STRONG_EM.sub(r"<strong><em>\1</em></strong>", s)
    s = _STRONG.sub(r"<strong>\1</strong>", s)
    s = _EM.sub(r"<em>\1</em>", s)
    s = _DEL.sub(r"<del>\1</del>", s)
    return _SLOT.sub(lambda m: slots[int(m.group(1))], s)


# -- tables ------------------------------------------------------------------

def split_row(line: str) -> list[str]:
    """The cells of a pipe-table row, unstripped of inline markup.

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


def _has_pipe(line: str) -> bool:
    return bool(re.search(r"(?<!\\)\|", line))


def _table_at(lines: list[str], i: int) -> list[str] | None:
    """The alignment of the table starting at line *i*, or None."""
    if i + 1 >= len(lines) or not _has_pipe(lines[i]):
        return None
    delim = lines[i + 1]
    if "-" not in delim or not re.fullmatch(r"[\s|:\-]+", delim):
        return None
    dcells = split_row(delim)
    if not dcells or not all(_DELIM_CELL.match(c) for c in dcells):
        return None
    if len(dcells) != len(split_row(lines[i])):
        return None
    aligns = []
    for c in dcells:
        if c.startswith(":") and c.endswith(":"):
            aligns.append("c")
        elif c.endswith(":"):
            aligns.append("r")
        elif c.startswith(":"):
            aligns.append("l")
        else:
            aligns.append("")
    return aligns


def _table_html(header: list[str], rows: list[list[str]],
                aligns: list[str]) -> str:
    width = len(header)
    rows = [(r + [""] * width)[:max(width, len(r))] for r in rows]
    ncols = max([width] + [len(r) for r in rows])
    header = header + [""] * (ncols - width)
    aligns = aligns + [""] * (ncols - len(aligns))

    def numeric(cell: str) -> bool:
        return bool(_NUMERIC.match(cell.replace("`", "")))

    # A column the model did not align but that holds only numbers reads
    # right-aligned, so the decimal places line up under each other.
    for c in range(ncols):
        if not aligns[c]:
            body = [r[c] for r in rows if c < len(r) and r[c]]
            if body and all(numeric(v) for v in body):
                aligns[c] = "r"

    def cell(tag: str, text: str, col: int) -> str:
        classes = []
        if aligns[col] in ("c", "r"):
            classes.append("al-" + aligns[col])
        if tag == "td" and text and numeric(text):
            classes.append("num")
        attr = f' class="{" ".join(classes)}"' if classes else ""
        return f"<{tag}{attr}>{inline(text)}</{tag}>"

    head = "".join(cell("th", h, c) for c, h in enumerate(header))
    body = "".join(
        "<tr>" + "".join(cell("td", v, c) for c, v in enumerate(r)) + "</tr>"
        for r in rows)
    return ('<div class="delfin-md-table-wrap"><table class="delfin-md-table">'
            f"<thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>")


# -- blocks ------------------------------------------------------------------

def _indent(line: str) -> int:
    return len(line) - len(line.lstrip(" "))


def _item(line: str):
    m = _ITEM.match(line)
    if not m:
        return None
    if m.group(3) is None and m.group(4) is None and line.strip() != m.group(2):
        return None
    marker = m.group(2)
    ordered = marker[-1] in ".)"
    spaces = len(m.group(3) or " ")
    if spaces > 4:
        spaces = 1
    content_col = len(m.group(1)) + len(marker) + spaces
    return {"indent": len(m.group(1)), "marker": marker, "ordered": ordered,
            "col": content_col, "text": m.group(4) or ""}


def _starts_block(lines: list[str], i: int, in_paragraph: bool) -> bool:
    line = lines[i]
    if _FENCE_OPEN.match(line) or _HEADING.match(line) or _HR.match(line):
        return True
    if _QUOTE.match(line) or _table_at(lines, i) is not None:
        return True
    if line.strip() == "$$":
        return True
    it = _item(line) if _indent(line) <= 3 else None
    if it and it["text"]:
        # An ordered item interrupts a sentence only when it starts at 1:
        # "in 2026. Then" on its own line is not a list.
        return not (in_paragraph and it["ordered"]
                    and it["marker"][:-1] != "1")
    return False


def _same_list(a: dict, b: dict) -> bool:
    """Whether item *b* continues the list *a* started: same kind and the
    same bullet character or ordinal delimiter, as CommonMark has it."""
    if a["ordered"] != b["ordered"]:
        return False
    return a["marker"][-1] == b["marker"][-1]


def _list(lines: list[str], i: int, cb: CodeBlockFn) -> tuple[str, int]:
    first = _item(lines[i])
    ordered = first["ordered"]
    items: list[list[str]] = []
    n = len(lines)
    loose = False
    while i < n:
        it = _item(lines[i])
        if (not it or not _same_list(first, it)
                or it["indent"] > first["indent"] + 3):
            break
        col = it["col"]
        body = [it["text"]]
        i += 1
        while i < n:
            line = lines[i]
            if not line.strip():
                # A blank line belongs to the item only if the item goes on.
                j = i
                while j < n and not lines[j].strip():
                    j += 1
                if j < n and _indent(lines[j]) >= min(col, first["indent"] + 2):
                    body.extend([""] * (j - i))
                    i = j
                    loose = True
                    continue
                break
            ind = _indent(line)
            sub = _item(line)
            if sub and ind < col:
                break                       # a sibling, or a new list
            if ind >= min(col, first["indent"] + 2):
                body.append(line[min(ind, col):])
                i += 1
                continue
            # Lazy continuation of the item's paragraph.
            if body[-1].strip() and not _starts_block(lines, i, True):
                body.append(line.strip())
                i += 1
                continue
            break
        items.append(body)
        if i < n and not lines[i].strip():
            # Blank lines between two items of the same list keep it one
            # list (a loose one); anything else ends it.
            j = i
            while j < n and not lines[j].strip():
                j += 1
            nxt = _item(lines[j]) if j < n else None
            if (nxt and _same_list(first, nxt)
                    and nxt["indent"] <= first["indent"] + 3):
                i = j
                loose = True
                continue
            break
    parts = []
    for body in items:
        inner = _blocks(body, cb)
        if not loose and inner.startswith("<p>"):
            end = inner.find("</p>")
            inner = inner[3:end] + inner[end + 4:]
        inner = re.sub(r"^\[ \] ", "\u2610 ", inner)
        inner = re.sub(r"^\[[xX]\] ", "\u2611 ", inner)
        parts.append(f"<li>{inner}</li>")
    if ordered:
        start = int(first["marker"][:-1])
        tag_open = f'<ol start="{start}">' if start != 1 else "<ol>"
        return tag_open + "".join(parts) + "</ol>", i
    return "<ul>" + "".join(parts) + "</ul>", i


def _blocks(lines: list[str], cb: CodeBlockFn) -> str:
    out: list[str] = []
    para: list[str] = []

    def flush() -> None:
        if para:
            out.append("<p>" + "<br>".join(inline(p.strip()) for p in para)
                       + "</p>")
            para.clear()

    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        if not line.strip():
            flush()
            i += 1
            continue

        m = _FENCE_OPEN.match(line)
        if m:
            flush()
            indent, fence = len(m.group(1)), m.group(2)
            lang = m.group(3).lower()
            if not _LANG.fullmatch(lang or "x"):
                lang = ""
            close = re.compile(
                r"^ {0,3}" + re.escape(fence[0]) + "{" + str(len(fence))
                + r",}[ \t]*$")
            body = []
            i += 1
            while i < n and not close.match(lines[i]):
                body.append(lines[i][min(indent, _indent(lines[i])):])
                i += 1
            i += 1                          # the closing fence, if any
            out.append(cb(lang, "\n".join(body).strip("\n")))
            continue

        if line.strip() == "$$" or (line.strip().startswith("$$")
                                    and line.strip().endswith("$$")
                                    and len(line.strip()) > 4):
            flush()
            s = line.strip()
            if s == "$$":
                body = []
                i += 1
                while i < n and lines[i].strip() != "$$":
                    body.append(lines[i])
                    i += 1
                i += 1
                src = "$$\n" + "\n".join(body) + "\n$$"
            else:
                src = s
                i += 1
            out.append('<div class="delfin-md-math">'
                       + _html.escape(src, quote=False) + "</div>")
            continue

        m = _HEADING.match(line)
        if m:
            flush()
            level = len(m.group(1))
            content = re.sub(r"[ \t]+#+[ \t]*$", "", m.group(2) or "")
            out.append(f"<h{level}>{inline(content)}</h{level}>")
            i += 1
            continue

        if _HR.match(line):
            flush()
            out.append("<hr>")
            i += 1
            continue

        aligns = _table_at(lines, i)
        if aligns is not None:
            flush()
            header = split_row(lines[i])
            rows = []
            i += 2
            while (i < n and lines[i].strip() and _has_pipe(lines[i])
                   and not _FENCE_OPEN.match(lines[i])):
                rows.append(split_row(lines[i]))
                i += 1
            out.append(_table_html(header, rows, aligns))
            continue

        if _QUOTE.match(line):
            flush()
            body = []
            while i < n and _QUOTE.match(lines[i]):
                body.append(_QUOTE.match(lines[i]).group(1))
                i += 1
            out.append("<blockquote>" + _blocks(body, cb) + "</blockquote>")
            continue

        it = _item(line)
        if it and (not para or _starts_block(lines, i, True)):
            flush()
            html_list, i = _list(lines, i, cb)
            out.append(html_list)
            continue

        para.append(line)
        i += 1
    flush()
    return "".join(out)
