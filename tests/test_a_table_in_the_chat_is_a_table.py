"""What the dashboard chat makes of a model's markdown.

Input: answer text as models write it (GitHub-style pipe tables, fenced
code, nested lists, headings, links). Output: the HTML the agent tab puts
into the page through ``_md_to_html``.

Reproduction this pins (real page, 2026-10-15): every table row came out as
one ``<code>`` line with the cells joined by `` | ``; the alignment row
left an empty line; a table indented under a list item was printed as raw
pipes; ``\\|`` split a cell; a four-backtick fence holding a three-backtick
one broke apart; ``~~~`` fences, links, block quotes, ``####`` headings and
nested lists were shown as their markup; ``2 * 3 * 4`` turned the 3 into
italics.

The sanitisation half is asserted on STRUCTURE (which elements and
attributes exist), not on substrings: escaped text may legitimately contain
``onerror=``.
"""

from __future__ import annotations

import re
from html.parser import HTMLParser

import pytest

from delfin.dashboard import tab_agent as T


def _html(md: str) -> str:
    return T._md_to_html(md)


class _Tags(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.tags: list[tuple[str, dict]] = []
        self.text: list[str] = []

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))

    def handle_data(self, data):
        self.text.append(data)


def _parse(html: str) -> _Tags:
    p = _Tags()
    p.feed(html)
    p.close()
    return p


def _cells(html: str, tag: str) -> list[str]:
    return [re.sub(r"<[^>]+>", "", c).strip()
            for c in re.findall(rf"<{tag}\b[^>]*>(.*?)</{tag}>", html, re.S)]


# -- tables ----------------------------------------------------------------

BASIC = (
    "Ergebnis:\n\n"
    "| Rechnung | Methode | Energie / Eh | Status |\n"
    "|:---------|:-------:|-------------:|--------|\n"
    "| `run_a` | PBE0/def2-SVP | -613.4171289 | ok |\n"
    "| `a | b` | **fett** | 1.5e-3 | a \\| b |\n"
    "\nDas Minimum ist `run_b`."
)


def test_a_pipe_table_becomes_a_table_element():
    out = _html(BASIC)
    assert "<table" in out and "</table>" in out
    assert _cells(out, "th") == ["Rechnung", "Methode", "Energie / Eh",
                                 "Status"]
    rows = re.findall(r"<tr>(.*?)</tr>", out.split("<tbody>", 1)[1], re.S)
    assert len(rows) == 2, "the alignment row is not a data row"


def test_the_alignment_row_sets_the_column_alignment():
    out = _html(BASIC)
    heads = re.findall(r"<th\b([^>]*)>", out)
    assert "al-c" in heads[1] and "al-r" in heads[2]
    assert "al-" not in heads[0] or "al-l" in heads[0]


def test_an_escaped_pipe_and_a_pipe_in_code_stay_inside_their_cell():
    out = _html(BASIC)
    body = out.split("<tbody>", 1)[1]
    second = re.findall(r"<tr>(.*?)</tr>", body, re.S)[1]
    tds = _cells(second, "td")
    assert len(tds) == 4, tds
    assert tds[0] == "a | b"
    assert tds[3] == "a | b"
    assert "<strong>fett</strong>" in second


def test_a_table_right_after_a_sentence_is_still_a_table():
    out = _html("Die Tabelle folgt direkt:\n| A | B |\n|---|---|\n| 1 | 2 |\n"
                "Und danach geht der Text weiter.")
    assert "<table" in out
    assert _cells(out, "td") == ["1", "2"]
    after = out.split("</table>", 1)[1]
    assert "Und danach geht der Text weiter." in after


def test_a_table_without_outer_pipes_is_a_table():
    out = _html("A | B\n--- | ---\n1 | 2")
    assert _cells(out, "td") == ["1", "2"]


def test_a_table_under_a_list_item_is_a_table():
    out = _html("1. Energien:\n\n   | Datei | E / Eh |\n   |---|---:|\n"
                "   | a.out | -76.42 |\n\n2. Weiter.")
    assert "<table" in out
    assert out.index("<table") > out.index("<li")
    assert "| Datei |" not in out


def test_lines_with_pipes_but_no_alignment_row_stay_text():
    """Not a table in any markdown dialect; shown as written."""
    out = _html("a | b\nc | d")
    assert "<table" not in out
    assert "a | b" in out


def test_a_wide_table_scrolls_inside_its_own_box():
    head = "| " + " | ".join(f"Spalte {i}" for i in range(10)) + " |"
    sep = "|" + "---|" * 10
    row = "| " + " | ".join("x" * 30 for _ in range(10)) + " |"
    out = _html("\n".join([head, sep, row]))
    assert re.search(r'<div class="delfin-md-table-wrap">\s*<table', out)


def test_a_number_cell_is_marked_for_tabular_figures():
    out = _html("| x | E |\n|---|---:|\n| a | -613.4171289 |\n| b | 1.5e-3 |")
    nums = re.findall(r'<td class="([^"]*)"[^>]*>(-613\.4171289|1\.5e-3)', out)
    assert len(nums) == 2 and all("num" in c for c, _ in nums)


def test_a_short_row_is_padded_to_the_header():
    out = _html("| a | b | c |\n|---|---|---|\n| 1 |")
    assert _cells(out, "td") == ["1", "", ""]


# -- code fences -----------------------------------------------------------

def test_a_longer_fence_can_hold_a_shorter_one():
    out = _html("````markdown\n```bash\norca x.inp\n```\n````\nnach")
    code = re.findall(r"<code>(.*?)</code>", out, re.S)
    assert any("```bash" in c and "orca x.inp" in c for c in code), out
    assert "nach" in out.split("</pre>")[-1]


def test_a_tilde_fence_is_a_code_block():
    out = _html("~~~bash\ngrep -n x out | sort\n~~~")
    assert "delfin-code-wrap" in out
    assert "grep" in out and "~~~" not in out


def test_an_unclosed_fence_runs_to_the_end():
    out = _html("Text\n\n```python\nx = 1")
    assert "delfin-code-wrap" in out and "```" not in out


def test_a_language_with_punctuation_is_kept_as_the_label():
    out = _html("```c++\nint x;\n```")
    assert "c++" in out and "int x;" in out
    assert "```" not in out


# -- blocks ----------------------------------------------------------------

def test_headings_one_to_four_are_headings():
    out = _html("# Eins\n## Zwei\n### Drei\n#### Vier")
    for n, word in enumerate(["Eins", "Zwei", "Drei", "Vier"], start=1):
        assert re.search(rf"<h{n}\b[^>]*>{word}</h{n}>", out), (n, out)
    assert "####" not in out


def test_a_nested_list_is_nested():
    out = _html("- eins\n- zwei\n  - a\n  - b\n    1. tief\n    2. tiefer\n- drei")
    tags = [t for t, _ in _parse(out).tags]
    assert tags.count("ul") == 2 and tags.count("ol") == 1
    assert tags.count("li") == 7
    assert re.search(r"<li>zwei\s*<ul>", out), out


def test_an_ordered_list_keeps_its_start_number():
    out = _html("3. drei\n4. vier")
    assert '<ol start="3">' in out


def test_a_block_quote_is_a_quote():
    out = _html("> Zitat\n> zweite Zeile")
    assert "<blockquote>" in out and "&gt;" not in out


def test_a_link_with_a_web_address_is_a_link():
    out = _html("Siehe [Doku](https://example.org/a_b?x=1&y=2).")
    links = [a for t, a in _parse(out).tags if t == "a"]
    assert len(links) == 1
    assert links[0]["href"] == "https://example.org/a_b?x=1&y=2"
    assert links[0].get("target") == "_blank"
    assert "noopener" in links[0].get("rel", "")


def test_paragraphs_are_separated_and_lines_kept():
    out = _html("Zeile eins\nZeile zwei\n\nNeuer Absatz")
    assert out.count("<p>") == 2
    assert "Zeile eins<br>Zeile zwei" in out


# -- inline ------------------------------------------------------------------

def test_a_product_is_not_italics():
    assert not re.search(r"<(em|i)>", _html("2 * 3 * 4 = 24"))


def test_underscores_inside_a_name_are_not_emphasis():
    out = _html("calc_a_b und a_b_c")
    assert "calc_a_b" in out and not re.search(r"<(em|i)>", out)


def test_inline_math_is_left_alone():
    out = _html("mit $a*b*c$ und $$x$$")
    assert not re.search(r"<(em|i)>", out) and "$a*b*c$" in out


def test_bold_italic_strike_and_code():
    out = _html("**fett** *kursiv* ~~weg~~ `x < y`")
    assert "<strong>fett</strong>" in out
    assert "<em>kursiv</em>" in out
    assert "<del>weg</del>" in out
    assert "<code>x &lt; y</code>" in out


def test_umlauts_survive():
    out = _html("Größe | Übergänge\n---|---\nÄpfel | ß")
    assert "Größe" in out and "Äpfel" in out and "ß" in out


# -- sanitisation ------------------------------------------------------------

# Every element the renderer itself emits. Anything else in the output came
# from the model's text, which must never happen.
_ALLOWED = {"p", "br", "strong", "em", "del", "code", "pre", "span", "div",
            "button", "table", "thead", "tbody", "tr", "th", "td", "ul",
            "ol", "li", "blockquote", "hr", "a", "h1", "h2", "h3", "h4",
            "h5", "h6"}

PAYLOADS = [
    "<script>alert(1)</script>",
    '<img src=x onerror="alert(1)">',
    '<iframe src="https://evil"></iframe>',
    "[klick](javascript:alert(1))",
    "[klick](JaVaScRiPt:alert(1))",
    "[d](data:text/html,<script>alert(2)</script>)",
    '[q](https://x.org/" onmouseover="alert(1))',
    "| <script>alert(3)</script> | <b>x</b> |\n|---|---|\n| <svg onload=alert(4)> | y |",
    "# <img src=x onerror=alert(5)>",
    "- <a href=javascript:alert(6)>x</a>\n  - <style>*{}</style>",
    "> <object data=x></object>",
    '```"><script>alert(7)</script>\ncode\n```',
    "`<script>alert(8)</script>`",
    "<https://x.org/\"onmouseover=alert(9)>",
    "\x00CB0\x00 \x00IL0\x00",
]


@pytest.mark.parametrize("payload", PAYLOADS)
def test_model_text_never_becomes_markup(payload):
    out = _html("vorher\n\n" + payload + "\n\nnachher")
    parsed = _parse(out)
    names = {t for t, _ in parsed.tags}
    assert names <= _ALLOWED, f"unexpected elements: {names - _ALLOWED}"
    for tag, attrs in parsed.tags:
        for k, v in attrs.items():
            if tag == "button" and k == "onclick":
                continue                   # the renderer's own copy button
            assert not k.startswith("on"), (tag, k, v)
            assert k in {"class", "href", "target", "rel", "start",
                         "data-codeidx", "data-mmd"}, (tag, k)
        if tag == "a":
            assert re.match(r"^(https?://|mailto:)", attrs["href"]), attrs
    assert "nachher" in out


def test_a_refused_link_is_shown_as_its_text():
    out = _html("[klick](javascript:alert(1))")
    assert "<a " not in out
    assert "klick" in out and "javascript:alert(1)" in out
