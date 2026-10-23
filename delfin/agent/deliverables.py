"""Rendered deliverables: durable reports instead of chat-only prose.

An agent answer scrolls away; a report file survives, can be opened in a
browser, attached to an ELN entry, or shared. ``publish_report`` writes a
markdown report under ``<workspace>/reports/`` and renders a standalone
HTML file next to it. The body is converted by the chat's renderer
(``delfin.dashboard.chat_markdown.render``), so a report draws tables,
lists, fences, quotes and links exactly as the chat drew the same answer;
this module adds the page around it and the page's stylesheet. Every
character of report text reaches the page escaped or not at all, and a
link only to http(s) or mailto -- report text can never inject markup.
"""

from __future__ import annotations

import html as _html
import time
from pathlib import Path

from . import german as _german


def _slugify(title: str, max_len: int = 48) -> str:
    """The report slug. The framework's one slug rule — see german.py.

    This was a SECOND implementation, and it disagreed with the memory
    store's about the same title: ``Müller GmbH`` came out ``m-ller-gmbh``
    here and ``mller-gmbh`` there, and neither told ``Möller GmbH`` apart
    from it.
    """
    return _german.slugify(title, max_len=max_len, fallback="report")


_CSS = (
    "body{font-family:system-ui,sans-serif;max-width:56rem;margin:2rem auto;"
    "padding:0 1rem;line-height:1.55;color:#1a1a1a;background:#fff}"
    "code{background:#f2f2f2;padding:.1em .3em;border-radius:3px;"
    "font-size:.92em}pre{background:#f6f6f6;padding:.8em;border-radius:6px;"
    "overflow-x:auto}pre code{background:none;padding:0}"
    "h1,h2,h3{line-height:1.25}hr{border:none;"
    "border-top:1px solid #ddd;margin:1.5em 0}"
    "blockquote{border-left:3px solid #ccc;margin-left:0;padding-left:1em;"
    "color:#444}"
    # The classes chat_markdown emits for a table: the wrapper scrolls a
    # wide table inside itself, al-c/al-r carry the alignment row, num
    # lines figures up digit under digit.
    ".delfin-md-table-wrap{overflow-x:auto;max-width:100%;margin:.6em 0}"
    ".delfin-md-table{border-collapse:collapse}"
    ".delfin-md-table th,.delfin-md-table td{border:1px solid #ccc;"
    "padding:.3em .6em;text-align:left;vertical-align:top}"
    ".delfin-md-table th{background:#f2f2f2;font-weight:600}"
    ".delfin-md-table .al-r{text-align:right}"
    ".delfin-md-table .al-c{text-align:center}"
    ".delfin-md-table td.num{font-variant-numeric:tabular-nums;"
    "white-space:nowrap}"
    ".delfin-md-math{background:#f6f6f6;padding:.4em .8em;border-radius:4px;"
    "font-family:monospace;white-space:pre-wrap}"
)


def _code_block(lang: str, code: str) -> str:
    """A fenced block as the report shows it: escaped, no highlighting.

    *lang* is not emitted: the page has no highlighter that would read it.
    """
    return "<pre><code>" + _html.escape(code, quote=False) + "</code></pre>"


def md_to_html(markdown: str, *, title: str = "") -> str:
    """Report markdown -> a standalone HTML page.

    Input: the report text and its title. Output: a complete document
    (doctype, charset, escaped ``<title>``, this module's stylesheet)
    whose body is ``chat_markdown.render(markdown)`` with plain escaped
    code blocks.

    Rejected: the line-by-line converter this replaced. It drew a table
    as header-less ``<td>`` rows, dropped the alignment row without
    effect, split cells at ``\\|``, and showed a table under a sentence
    or a list item, and every ordered list, as paragraphs -- a second
    parser that disagreed with the chat about the same answer.
    """
    # Imported here: the dashboard package is not needed to import this
    # module, only to render a page.
    from delfin.dashboard.chat_markdown import render

    body = render(markdown or "", _code_block)
    t = _html.escape(title or "DELFIN report")
    return (
        "<!doctype html><html><head><meta charset=\"utf-8\">"
        f"<title>{t}</title><style>{_CSS}</style></head><body>\n"
        + body + "\n</body></html>\n"
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
