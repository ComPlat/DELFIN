"""Rendered deliverables: durable reports instead of chat-only prose.

An agent answer scrolls away; a report file survives, can be opened in a
browser, attached to an ELN entry, or shared. ``publish_report`` writes a
markdown report under ``<workspace>/reports/`` and renders a standalone
HTML file next to it with a deterministic, dependency-free converter
(headers, fenced code, inline code, bold/italic, lists, links, tables,
rules). All content is HTML-escaped first — report text can never inject
markup into the rendered page.
"""

from __future__ import annotations

import html as _html
import re
import time
from pathlib import Path


_SLUG_RE = re.compile(r"[^a-z0-9]+")


def _slugify(title: str, max_len: int = 48) -> str:
    slug = _SLUG_RE.sub("-", (title or "report").lower()).strip("-")
    return (slug or "report")[:max_len].rstrip("-")


_CSS = (
    "body{font-family:system-ui,sans-serif;max-width:56rem;margin:2rem auto;"
    "padding:0 1rem;line-height:1.55;color:#1a1a1a;background:#fff}"
    "code{background:#f2f2f2;padding:.1em .3em;border-radius:3px;"
    "font-size:.92em}pre{background:#f6f6f6;padding:.8em;border-radius:6px;"
    "overflow-x:auto}pre code{background:none;padding:0}"
    "table{border-collapse:collapse}td,th{border:1px solid #ccc;"
    "padding:.3em .6em}h1,h2,h3{line-height:1.25}hr{border:none;"
    "border-top:1px solid #ddd;margin:1.5em 0}"
    "blockquote{border-left:3px solid #ccc;margin-left:0;padding-left:1em;"
    "color:#444}"
)

_INLINE_CODE = re.compile(r"`([^`\n]+)`")
_BOLD = re.compile(r"\*\*([^*\n]+)\*\*")
_ITALIC = re.compile(r"(?<!\*)\*([^*\n]+)\*(?!\*)")
_LINK = re.compile(r"\[([^\]\n]+)\]\((https?://[^)\s]+)\)")


def _inline(text: str) -> str:
    """Inline markdown on an ALREADY-ESCAPED line."""
    text = _INLINE_CODE.sub(r"<code>\1</code>", text)
    text = _BOLD.sub(r"<strong>\1</strong>", text)
    text = _ITALIC.sub(r"<em>\1</em>", text)
    text = _LINK.sub(r'<a href="\2">\1</a>', text)
    return text


def md_to_html(markdown: str, *, title: str = "") -> str:
    """Deterministic markdown -> standalone HTML page. Content is escaped
    before any markup is generated; unsupported constructs degrade to
    plain paragraphs rather than breaking the page."""
    out: list[str] = []
    in_code = False
    in_list = False
    in_table = False

    def _close_list() -> None:
        nonlocal in_list
        if in_list:
            out.append("</ul>")
            in_list = False

    def _close_table() -> None:
        nonlocal in_table
        if in_table:
            out.append("</table>")
            in_table = False

    for raw in (markdown or "").splitlines():
        line = raw.rstrip("\n")
        if line.strip().startswith("```"):
            _close_list(); _close_table()
            out.append("<pre><code>" if not in_code else "</code></pre>")
            in_code = not in_code
            continue
        if in_code:
            out.append(_html.escape(line))
            continue
        esc = _html.escape(line)
        stripped = esc.strip()
        if not stripped:
            _close_list(); _close_table()
            continue
        m = re.match(r"^(#{1,4})\s+(.*)$", stripped)
        if m:
            _close_list(); _close_table()
            level = len(m.group(1))
            out.append(f"<h{level}>{_inline(m.group(2))}</h{level}>")
            continue
        if re.match(r"^(-{3,}|\*{3,})$", stripped):
            _close_list(); _close_table()
            out.append("<hr>")
            continue
        if stripped.startswith("|") and stripped.endswith("|"):
            cells = [c.strip() for c in stripped.strip("|").split("|")]
            if all(re.match(r"^:?-{2,}:?$", c) for c in cells):
                continue  # separator row
            if not in_table:
                _close_list()
                out.append("<table>")
                in_table = True
            out.append(
                "<tr>" + "".join(f"<td>{_inline(c)}</td>" for c in cells)
                + "</tr>")
            continue
        _close_table()
        m = re.match(r"^[-*]\s+(.*)$", stripped)
        if m:
            if not in_list:
                out.append("<ul>")
                in_list = True
            out.append(f"<li>{_inline(m.group(1))}</li>")
            continue
        _close_list()
        if stripped.startswith("&gt;"):
            out.append(
                f"<blockquote>{_inline(stripped[4:].strip())}</blockquote>")
            continue
        out.append(f"<p>{_inline(stripped)}</p>")

    if in_code:
        out.append("</code></pre>")
    _close_list(); _close_table()
    t = _html.escape(title or "DELFIN report")
    return (
        "<!doctype html><html><head><meta charset=\"utf-8\">"
        f"<title>{t}</title><style>{_CSS}</style></head><body>\n"
        + "\n".join(out) + "\n</body></html>\n"
    )


def publish_report(
    workspace: Path | str,
    *,
    title: str,
    markdown: str,
    fmt: str = "both",
) -> dict:
    """Write a report under ``<workspace>/reports/`` and return the paths.

    ``fmt``: "md" | "html" | "both". Filenames are date-stamped and
    slug-derived; an existing file gets a numeric suffix instead of being
    overwritten (raw results are never clobbered).
    """
    ws = Path(workspace)
    reports = ws / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%d")
    base = f"{stamp}-{_slugify(title)}"
    def _free(suffix: str) -> Path:
        p = reports / f"{base}{suffix}"
        n = 2
        while p.exists():
            p = reports / f"{base}-{n}{suffix}"
            n += 1
        return p
    result: dict = {"title": title}
    body = (markdown or "").strip() + "\n"
    if fmt in ("md", "both"):
        mp = _free(".md")
        mp.write_text(f"# {title}\n\n{body}", encoding="utf-8")
        result["md_path"] = str(mp)
    if fmt in ("html", "both"):
        hp = _free(".html")
        hp.write_text(md_to_html(f"# {title}\n\n{body}", title=title),
                      encoding="utf-8")
        result["html_path"] = str(hp)
    return result
