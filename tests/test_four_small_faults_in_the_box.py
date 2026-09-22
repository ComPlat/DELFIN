"""The four the drawing fix named and left standing.

None of them destroys the transcript — that was the other pass — but
each is wrong where somebody would notice it, and all four were found by
reading rather than by any test:

  viewport's cursor    a "…" marker row above the window is a row of the
                       picture, so every content row moves down by one;
                       the cursor did not, and sat on a full line of
                       text instead of at its end the moment the box
                       scrolled

  viewport's rule      with no hint row, `rows[1:]` swept the BOTTOM
                       rule in as content and then appended it again:
                       two closing rules, one above the other

  the narrow form      below MIN_WIDTH the box draws no rules, and the
                       caller added the border's two columns anyway, so
                       the cursor stood two columns past the text. The
                       same row measured its cut in characters, which is
                       half a glyph out for CJK

  Ctrl+G at idle       the decoder takes the line into the event and
                       clears its buffer; `_collect` had no branch, so
                       the typed line vanished and nothing happened.
                       Ctrl+O and Ctrl+T did the same
"""

from __future__ import annotations

import pytest

from delfin.agent import repl_box as rb


HINT = "esc interrupt · /help"


def _rule_rows(view, width: int = 40) -> list[int]:
    """Indices of the rule rows in *view*: the top and bottom of the box.

    The two rules are the same string, so a rule is found by matching it
    whole rather than by a corner character. Not by ``hint_row`` either:
    ``viewport`` copies that index from the view it windowed, so on a
    windowed view it points into a picture that is no longer there.
    """
    rule = rb._HORIZONTAL * (width - 1)
    return [i for i, r in enumerate(view.rows) if r == rule]


# -- the window's cursor ----------------------------------------------------

def test_a_scrolled_window_keeps_the_cursor_on_the_text():
    text = "y" * 350
    view = rb.viewport(rb.render_box(text, len(text), 40, HINT), 6)
    row, col = view.cursor
    # The cursor's row must be the LAST content row of the window: that
    # is where the end of the text is.
    # A content row has no border character to be recognised by any
    # more, so it is named by position instead: everything between the
    # two rules that is not a "…" marker.
    rules = _rule_rows(view)
    content = [i for i in range(rules[0] + 1, rules[-1])
               if view.rows[i] != "…"]
    assert content, view.rows
    last_content_index = content[-1]          # index into view.rows
    assert row + 1 == last_content_index, (
        f"cursor row {row} is not the last content row "
        f"{last_content_index - 1}:\n" + "\n".join(view.rows))


def test_an_unscrolled_window_is_unchanged():
    view = rb.render_box("short", 5, 40, HINT)
    assert rb.viewport(view, 6) == view


# -- the window's rules -----------------------------------------------------

# The two rules are the same string now, so "how many rows start with ╰"
# has no answer: the count that carries the same meaning is WHERE the
# rule rows are. Exactly two, the first and the last of the picture —
# the duplicated closing rule showed up as a third, between them.

def test_a_window_without_a_hint_has_one_closing_rule():
    text = "z" * 350
    view = rb.viewport(rb.render_box(text, len(text), 40, ""), 6)
    assert _rule_rows(view) == [0, len(view.rows) - 1], "\n".join(view.rows)


def test_a_window_with_a_hint_still_has_one():
    text = "z" * 350
    view = rb.viewport(rb.render_box(text, len(text), 40, HINT), 6)
    # With a hint the closing rule is second from the bottom; the hint
    # sits under it.
    assert _rule_rows(view) == [0, len(view.rows) - 2], "\n".join(view.rows)
    assert view.rows[-1].strip().startswith("esc")


# -- the narrow form --------------------------------------------------------

def test_the_narrow_form_says_it_has_no_border():
    view = rb.render_box("hello", 5, 6)
    assert view.border is False, "the caller must not offset for a frame"


def test_the_wide_form_reports_the_column_the_caller_can_use():
    """Replaces "a framed box says it has one".

    There is no framed form left: the wide form draws two rules and no
    side borders, so `border` is False for it too and the flag no longer
    tells the two forms apart. What it still has to do is the thing the
    flag was for — tell the caller NOT to add a border's two columns,
    because the reported column is measured from the start of the row
    itself. So that is what is checked: the flag, and the column it
    promises, against the row it indexes into.
    """
    view = rb.render_box("hello", 5, 40, HINT)
    assert view.border is False, "the caller must not offset for a frame"
    assert view.rows[1] == "> hello"
    assert view.cursor == (0, rb.string_width("> hello"))


def test_the_narrow_cursor_is_measured_in_columns():
    """A CJK character is one character and two columns wide."""
    text = "日本語テキスト"
    narrow = rb.render_box(text, len(text), 10)
    row = narrow.rows[0]
    assert narrow.cursor[1] <= rb.string_width(row), (
        f"cursor column {narrow.cursor[1]} past the row's "
        f"{rb.string_width(row)} columns: {row!r}")


# -- the keys that vanished -------------------------------------------------

def test_the_idle_loop_puts_the_line_back_for_turn_time_keys():
    """Ctrl+G, Ctrl+O and Ctrl+T belong to a running turn. At the idle
    prompt there is none — and the line the decoder took must come
    back, because a key that empties the box and does nothing is the
    one outcome a key must never have."""
    import ast
    import inspect
    from delfin.agent import repl as R

    tree = ast.parse(inspect.getsource(R))
    collect = next(n for n in ast.walk(tree)
                   if isinstance(n, ast.FunctionDef) and n.name == "_collect")
    src = ast.unparse(collect)
    assert "rk.STEER" in src and "rk.EXPAND" in src and "rk.TASKS" in src, (
        "a turn-time key still falls through and eats the line")
    assert "decoder.buffer = text" in src
