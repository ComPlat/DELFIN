"""A pipe table in the plain terminal answer is drawn as aligned columns.

Input: an answer streamed through ``MarkdownStream`` with colour on (a
terminal). Output: what is written to stdout.

Before: tables were explicitly not handled, so a model's table reached
the screen as its source -- ragged ``| a | bb |`` lines and a
``|---|---:|`` row -- while everything around it was styled.

The terminal cannot repaint, so a table is held back from its first row
until the first line that does not start with ``|``, then drawn at once.
Rules kept from the stream contract: colour off changes no byte; where a
table does not fit the width, the rows are written as the model wrote
them; how the input is chunked never changes the result.
"""

from __future__ import annotations

import re
import unicodedata

import pytest

from delfin.agent.repl_render import MarkdownStream, Theme

ON = Theme(enabled=True)
OFF = Theme(enabled=False)

TABLE = (
    "Ergebnis:\n"
    "| Run | Energie / Eh | Status |\n"
    "|:----|-------------:|--------|\n"
    "| `A` | -25.1844 | ✅ ok |\n"
    "| B | -125.2019 | **nicht** konvergiert |\n"
    "Danach Text.\n"
)


def _plain(styled: str) -> str:
    return re.sub(r"\x1b\[[0-9;]*m", "", styled)


def _cols(s: str) -> int:
    return sum(2 if unicodedata.east_asian_width(c) in "WF" else 1
               for c in s if not unicodedata.combining(c))


def _whole(text, theme=ON, width=100):
    md = MarkdownStream(theme, width=width)
    return md.feed(text) + md.flush()


def _by_char(text, theme=ON, width=100):
    md = MarkdownStream(theme, width=width)
    return "".join(md.feed(c) for c in text) + md.flush()


def _table_lines(out: str) -> list[str]:
    return [ln for ln in _plain(out).splitlines() if "│" in ln or "─" in ln]


def test_a_table_is_drawn_in_aligned_columns():
    out = _whole(TABLE)
    rows = _table_lines(out)
    assert len(rows) == 4, _plain(out)          # header, rule, two rows
    assert len({_cols(r) for r in rows}) == 1, rows
    assert "|---" not in _plain(out) and "-----:" not in _plain(out)
    # The column separators sit at the same columns in every row.
    seps = [[i for i, ch in enumerate(r) if ch in "│┼"] for r in rows
            if "✅" not in r]
    assert all(s == seps[0] for s in seps), rows


def test_markup_inside_a_cell_is_styled_not_printed():
    out = _whole(TABLE)
    plain = _plain(out)
    assert "`A`" not in plain and "**nicht**" not in plain
    assert "nicht konvergiert" in plain


def test_a_right_aligned_column_is_right_aligned():
    rows = _table_lines(_whole(TABLE))
    a = next(r for r in rows if "-25.1844" in r)
    b = next(r for r in rows if "-125.2019" in r)
    assert a.index("-25.1844") + len("-25.1844") == \
        b.index("-125.2019") + len("-125.2019")


def test_the_text_after_the_table_is_not_a_row():
    plain = _plain(_whole(TABLE))
    assert plain.splitlines()[0] == "Ergebnis:"
    assert plain.rstrip("\n").splitlines()[-1] == "Danach Text."


def test_the_chunking_does_not_change_the_table():
    ref = _whole(TABLE)
    assert _by_char(TABLE) == ref
    for i in range(1, len(TABLE), 7):
        md = MarkdownStream(ON, width=100)
        assert md.feed(TABLE[:i]) + md.feed(TABLE[i:]) + md.flush() == ref, i


def test_colour_off_changes_no_byte():
    assert _whole(TABLE, OFF) == TABLE
    assert _by_char(TABLE, OFF) == TABLE


def test_a_table_too_wide_for_the_terminal_is_written_as_is():
    head = "| " + " | ".join(f"Spalte {i}" for i in range(10)) + " |\n"
    sep = "|" + "---|" * 10 + "\n"
    row = "| " + " | ".join("wert" * 3 for _ in range(10)) + " |\n"
    src = head + sep + row
    assert _plain(_whole(src, width=60)) == src


def test_pipes_without_an_alignment_row_are_text():
    src = "| a | b |\n| c | d |\n"
    assert _plain(_whole(src)) == src


def test_a_table_inside_a_fence_is_code():
    src = "```\n| a | b |\n|---|---|\n```\n"
    plain = _plain(_whole(src))
    assert "| a | b |\n|---|---|\n" in plain and "│" not in plain


def test_a_table_at_the_very_end_is_drawn_on_flush():
    src = "| a | b |\n|---|---|\n| 1 | 2 |"
    rows = _table_lines(_whole(src))
    assert len(rows) == 3


@pytest.mark.parametrize("width", [40, 80, 120])
def test_no_drawn_row_is_wider_than_the_terminal(width):
    for r in _table_lines(_whole(TABLE, width=width)):
        assert _cols(r) <= width


def test_a_held_table_is_written_before_the_next_tool_line():
    """The table ends the model's text; a tool line that follows it in the
    same turn must not overtake it on the screen."""
    import io

    from delfin.agent import repl

    screen = io.StringIO()                  # stdout and stderr: one terminal
    t = repl.Transcript(screen, screen, theme=Theme(enabled=True),
                        width=100, color="always")
    t.render(repl.RenderItem("text", "Tabelle:\n| a | b |\n|---|---|\n"
                                     "| 1 | 2 |\n"))
    t.render(repl.RenderItem("tool_use", '{"command": "ls"}', name="bash"))
    t.finish()
    plain = _plain(screen.getvalue())
    assert "│" in plain
    assert plain.index("│") < plain.index("ls"), plain
