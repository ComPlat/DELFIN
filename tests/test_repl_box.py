"""The framed input box, as a table of cases.

Every row of the table is one layout question a framed input line has
to answer, checked against a hand-drawn expectation. The renderer is
pure — no terminal anywhere — so each case is just input in, rows out.
"""

import pytest

from delfin.agent.repl_box import (
    MIN_WIDTH, PROMPT, BoxView, char_width, render_box, string_width,
)

HINT = "esc interrupt · shift+tab approval mode · /help"


def box_rows(text, cursor=None, width=60, hint=""):
    view = render_box(text, len(text) if cursor is None else cursor,
                      width, hint)
    return view


# -- the table -----------------------------------------------------------

CASES = [
    pytest.param(
        # name, text, cursor, width, hint, expected content rows,
        # expected cursor (row, col)
        "empty", "", 0, 40, "",
        [PROMPT.strip()],
        (0, 2),
        id="empty",
    ),
    pytest.param(
        "one short line", "hello", 5, 40, "",
        ["> hello"],
        (0, 7),
        id="short",
    ),
    pytest.param(
        "cursor mid-word", "hello", 2, 40, "",
        ["> hello"],
        (0, 4),
        id="cursor-middle",
    ),
    pytest.param(
        "wraps at the inner width",
        "a" * 39, 39, 40, "",
        ["> " + "a" * 34, "a" * 5],
        (1, 5),
        id="wrap",
    ),
    pytest.param(
        "cursor on the first wrapped row",
        "a" * 39, 10, 40, "",
        ["> " + "a" * 34, "a" * 5],
        (0, 12),
        id="cursor-first-row",
    ),
    pytest.param(
        "CJK wide chars measure two columns",
        "日本語", 3, 40, "",
        ["> 日本語"],
        (0, 8),
        id="cjk",
    ),
    pytest.param(
        "seventeen wide chars exactly fill the row",
        "あ" * 17, 34, 40, "",
        # inner = 36; "> " + 17 wide chars = 36 columns exactly.
        ["> " + "あ" * 17],
        (0, 36),
        id="wide-exact-fill",
    ),
    pytest.param(
        "wide char straddling moves down whole",
        "a" * 34 + "あ", 35, 40, "",
        # inner=36, "> " + 34 a = 36 columns; the wide char cannot fit
        # in 0 remaining columns and moves down.
        ["> " + "a" * 34, "あ"],
        (1, 2),
        id="wide-straddle",
    ),
    pytest.param(
        "newline hard-wraps (a pasted block)",
        "one\ntwo", 7, 40, "",
        ["> one", "two"],
        (1, 3),
        id="newline",
    ),
    pytest.param(
        "cursor clamped past the end",
        "hi", 99, 40, "",
        ["> hi"],
        (0, 4),
        id="cursor-clamped",
    ),
    pytest.param(
        "empty text, hint present",
        "", 0, 40, HINT,
        [PROMPT.strip()],
        (0, 2),
        id="hint-empty-text",
    ),
]


@pytest.mark.parametrize(
    "name,text,cursor,width,hint,expected_content,expected_cursor", CASES)
def test_box_table(name, text, cursor, width, hint,
                   expected_content, expected_cursor):
    view = render_box(text, cursor, width, hint)
    content = view.rows[1:-1] if hint == "" else view.rows[1:-2]
    # Content rows carry the frame; compare the inner text.
    for i, row in enumerate(content):
        assert row.startswith("│ "), name
        assert row.endswith(" │"), name
        inner = row[2:-2]
        # Strip padding for the comparison, but the measured width must
        # be exact: the right border lines up because every inner row is
        # padded to the same width.
        assert string_width(inner.rstrip()) == string_width(
            expected_content[i].rstrip()), (name, inner)
        assert inner.rstrip() == expected_content[i].rstrip(), name
        assert string_width(inner) == width - 4, (name, row)
    assert view.cursor == expected_cursor, (name, view.cursor)


# -- the frame -----------------------------------------------------------

def test_frame_rows_and_alignment():
    view = render_box("hi", 2, 40, HINT)
    assert view.rows[0] == "╭" + "─" * 38 + "╮"
    assert view.rows[-2] == "╰" + "─" * 38 + "╯"
    # Every bordered row is exactly the terminal width, so nothing
    # wraps. (The hint row is deliberately shorter — it sits under the
    # frame like the example, not inside it.)
    for row in view.rows[:view.hint_row or len(view.rows)]:
        if row.startswith(("╭", "│", "╰")):
            assert string_width(row) == 40, row


def test_hint_is_last_row_and_never_wraps():
    view = render_box("", 0, 40, HINT)
    assert view.hint_row == len(view.rows) - 1
    assert string_width(view.rows[view.hint_row]) <= 40


def test_hint_that_does_not_fit_keeps_its_end():
    view = render_box("", 0, 24, "esc interrupt · shift+tab approval · /help")
    row = view.rows[view.hint_row]
    assert string_width(row) <= 24
    assert row.endswith("/help"), row
    assert row.lstrip().startswith("…"), row


def test_no_hint_leaves_no_hint_row():
    view = render_box("x", 1, 40)
    assert view.hint_row is None
    assert len(view.rows) == 3


# -- below the minimum width --------------------------------------------

def test_below_min_width_degenerates_to_one_row():
    view = render_box("hello world this is long", 24, 10, HINT)
    assert len(view.rows) == 1
    assert string_width(view.rows[0]) <= 10
    # END of the text survives.
    assert view.rows[0].endswith("long"), view.rows[0]


def test_narrow_row_cursor_is_visible():
    view = render_box("hello world this is long", 24, 10)
    assert view.cursor == (0, 9), view.cursor


def test_at_min_width_has_a_frame():
    view = render_box("", 0, MIN_WIDTH)
    assert view.rows[0].startswith("╭")
    assert len(view.rows) == 3


def test_width_zero_does_not_crash():
    view = render_box("anything", 3, 0, HINT)
    assert isinstance(view, BoxView)


# -- width measurement ---------------------------------------------------

@pytest.mark.parametrize("ch,expected", [
    ("a", 1), ("あ", 2), ("Ａ", 2),  # fullwidth latin
    ("́", 0),                          # combining acute
    (" ", 1),
])
def test_char_width(ch, expected):
    assert char_width(ch) == expected


def test_string_width_sums_columns_not_characters():
    assert string_width("aあb") == 4
