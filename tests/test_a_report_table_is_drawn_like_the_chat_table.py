"""A published report draws markdown the way the chat draws it.

Input: report markdown as ``publish_report`` receives it from the model.
Output: the standalone HTML page ``deliverables.md_to_html`` writes.

Semantics pinned here: the page body is ``chat_markdown.render`` of the
same text, so a table is recognised by the same rules
(``delfin.agent.markdown_table``) and emitted as the same element tree --
alignment row honoured, escaped ``\\|`` and a pipe inside a code span kept
in their cell, a table directly under a sentence or indented under a list
item still a table. The page keeps its own guarantees: a complete
document with its own stylesheet (now also styling the alignment and
number classes), the title escaped, every character of model text
escaped or dropped, links only to http(s) and mailto, no script, no event
attribute.

Reproduction (before this change): the page used a line-by-line converter
of its own. ``| a | b |`` lines became ``<tr><td>`` rows with no header,
the alignment row was dropped without effect, ``\\|`` split a cell, a
table under a sentence or under a list item came out as a paragraph of
pipes, and an ordered list as paragraphs.
"""

from __future__ import annotations

import html as _html
import re

import pytest

from delfin.agent.deliverables import md_to_html
from delfin.dashboard import chat_markdown
from tests.test_a_table_in_the_chat_is_a_table import PAYLOADS, _parse


def _plain_code(lang: str, code: str) -> str:
    return "<pre><code>" + _html.escape(code, quote=False) + "</code></pre>"


def _tables(html: str) -> list[str]:
    return re.findall(
        r'<div class="delfin-md-table-wrap">.*?</table></div>', html, re.S)


TABLES = {
    "aligned": ("| Rechnung | Methode | E / Eh | Status |\n"
                "|:--|:-:|--:|:-:|\n"
                "| run_a | PBE0 | -613.4171289 | ok |\n"
                "| run_b | B3LYP | -613.4012 | fehlgeschlagen |\n"),
    "escaped_pipe": ("| Ausdruck | Bedeutung |\n|---|---|\n"
                     "| `a | b` | oder |\n| x \\| y | literal |\n"),
    "after_sentence": ("Ergebnis:\n| a | b |\n|---|---|\n| 1 | 2 |\n"),
    "under_list_item": ("- Lauf A:\n\n  | E | dE |\n  |---|---:|\n"
                        "  | -1.5 | 0.0 |\n"),
    "short_rows": ("| a | b | c |\n|---|---|---|\n| 1 |\n| 1 | 2 | 3 | 4 |\n"),
    "numbers_unaligned": ("| T / K | k / s-1 |\n|---|---|\n"
                          "| 298.15 | 1.2e-3 |\n| 310 | −2515 |\n"),
}


@pytest.mark.parametrize("name", sorted(TABLES))
def test_a_report_table_is_the_chat_table(name):
    md = TABLES[name]
    chat = _tables(chat_markdown.render(md, _plain_code))
    page = _tables(md_to_html(md, title="t"))
    assert chat, f"chat renders no table for {name}"
    assert page == chat


def test_the_alignment_row_reaches_the_page():
    page = md_to_html(TABLES["aligned"], title="t")
    head = re.search(r"<thead>.*?</thead>", page, re.S).group(0)
    assert re.findall(r"<th(?: class=\"([^\"]*)\")?>", head) == [
        "", "al-c", "al-r", "al-c"]
    assert '<td class="al-r num">-613.4171289</td>' in page
    assert "<thead>" in page and "<tbody>" in page


def test_the_page_styles_what_the_table_carries():
    page = md_to_html(TABLES["aligned"], title="t")
    css = re.search(r"<style>(.*?)</style>", page, re.S).group(1)
    for rule in (".al-r", ".al-c", ".num", ".delfin-md-table"):
        assert rule in css, rule


def test_an_ordered_list_and_a_nested_list_are_lists():
    page = md_to_html("1. eins\n2. zwei\n   - a\n   - b\n", title="t")
    assert "<ol><li>eins</li><li>zwei<ul><li>a</li><li>b</li></ul></li></ol>" \
        in page


def test_the_page_is_a_complete_document_with_an_escaped_title():
    page = md_to_html("text", title="<script>alert(1)</script>")
    assert page.startswith("<!doctype html><html><head><meta charset=\"utf-8\">")
    assert "<title>&lt;script&gt;alert(1)&lt;/script&gt;</title>" in page
    assert page.rstrip().endswith("</body></html>")


_ALLOWED = {"html", "head", "meta", "title", "style", "body",
            "p", "br", "strong", "em", "del", "code", "pre", "div",
            "table", "thead", "tbody", "tr", "th", "td", "ul", "ol", "li",
            "blockquote", "hr", "a", "h1", "h2", "h3", "h4", "h5", "h6"}


@pytest.mark.parametrize("payload", PAYLOADS)
def test_model_text_never_becomes_markup_in_a_report(payload):
    page = md_to_html("vorher\n\n" + payload + "\n\nnachher", title=payload)
    parsed = _parse(page)
    names = {t for t, _ in parsed.tags}
    assert names <= _ALLOWED, f"unexpected elements: {names - _ALLOWED}"
    assert [t for t, _ in parsed.tags].count("style") == 1
    for tag, attrs in parsed.tags:
        for k, v in attrs.items():
            assert not k.startswith("on"), (tag, k, v)
            assert k in {"class", "href", "target", "rel", "start",
                         "charset"}, (tag, k)
        if tag == "a":
            assert re.match(r"^(https?://|mailto:)", attrs["href"]), attrs
    assert "nachher" in page
