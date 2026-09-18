"""The input area between two rules, as a table of cases.

Every row of the table is one layout question the input line has to
answer, checked against a hand-drawn expectation. The renderer is
pure — no terminal anywhere — so each case is just input in, rows out.

There is no frame: a content row is the raw wrapped text, with no side
borders and no padding, and it wraps at ``width - 1`` columns — the
length of the rules above and below it.
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
        # The row IS the prompt, trailing space and all: there is no
        # padding to strip and the cursor sits in the column after it.
        [PROMPT],
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
        # inner = width - 1 = 39: "> " plus 37 of the 39 a's fill row 0.
        "a" * 39, 39, 40, "",
        ["> " + "a" * 37, "a" * 2],
        (1, 2),
        id="wrap",
    ),
    pytest.param(
        "cursor on the first wrapped row",
        "a" * 39, 10, 40, "",
        ["> " + "a" * 37, "a" * 2],
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
        "nineteen wide chars exactly fill the row",
        # Width 41, not 40: the case exists to check the EXACT fill, and
        # "> " plus wide characters is always an even number of columns,
        # so it can only land on an even inner width. inner = 40;
        # "> " + 19 wide chars = 40 columns exactly.
        "あ" * 19, 19, 41, "",
        ["> " + "あ" * 19],
        (0, 40),
        id="wide-exact-fill",
    ),
    pytest.param(
        "wide char straddling moves down whole",
        "a" * 36 + "あ", 37, 40, "",
        # inner=39, "> " + 36 a = 38 columns; the wide char needs two
        # and one column is left, so it moves down whole.
        ["> " + "a" * 36, "あ"],
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
        [PROMPT],
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
    assert len(content) == len(expected_content), (name, content)
    for i, row in enumerate(content):
        # A content row is the raw wrapped text — no left border, no
        # right border, no padding — so it must match the expectation
        # CHARACTER for character. (The old form compared the two
        # rstripped and then checked the padding separately, because the
        # row was padded out to the right border; there is no padding to
        # allow for now, and an exact match is what says so.)
        assert row == expected_content[i], (name, row)
        # What the padding used to buy — every row ending in the same
        # column — the wrap width buys instead: no row is wider than the
        # rules above and below it, so the terminal never wraps one.
        assert string_width(row) <= width - 1, (name, row)
    assert view.cursor == expected_cursor, (name, view.cursor)


# -- the rules -----------------------------------------------------------

def test_rule_rows_and_alignment():
    view = render_box("hi", 2, 40, HINT)
    rule = "─" * 39
    assert view.rows[0] == rule
    assert view.rows[-2] == rule
    # Alignment used to mean a right border standing in one column on
    # every row. With no right border the property that survives is the
    # one the rules define: the two are identical and one column short
    # of the terminal, so neither wraps, and no row between them runs
    # past them. (The hint row sits under the lower rule, as before.)
    for row in view.rows[:view.hint_row or len(view.rows)]:
        assert string_width(row) <= string_width(rule), row


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


def test_at_min_width_the_rules_are_drawn():
    # At exactly MIN_WIDTH the full form is still used, not the
    # degenerate single row: a rule above, the text, a rule below.
    view = render_box("", 0, MIN_WIDTH)
    assert view.rows[0] == "─" * (MIN_WIDTH - 1)
    assert view.rows[-1] == "─" * (MIN_WIDTH - 1)
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
