"""The bug control is folded away until it is wanted.

The agent tab kept a text field and a labelled button standing open in
the control row — the row where the session is DRIVEN: mode, provider,
model, effort, permissions, stop. Two controls' worth of width, on every
screen, for the one thing a user does when something has already gone
wrong.

The viewer had already solved this: one beetle that opens the line, the
send button, and a sentence saying what is held and where it goes.

These build the REAL tab and work the REAL toggle. An earlier draft of
this file asserted substrings of `create_tab`'s source, which is the
mistake this codebase has been bitten by twice: a test that pins the
text of a rule passes when the rule has stopped being wired to anything.
"""

from __future__ import annotations

import pytest

pytest.importorskip("ipywidgets")

import ipywidgets as widgets                                # noqa: E402

from delfin.dashboard import tab_agent                      # noqa: E402
from delfin.dashboard.context import DashboardContext       # noqa: E402


def _walk(widget):
    yield widget
    for child in getattr(widget, "children", ()) or ():
        yield from _walk(child)


@pytest.fixture
def bug_group(tmp_path):
    """The four controls, as the built tab actually holds them."""
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    tab, _refs = tab_agent.create_tab(ctx)

    for node in _walk(tab):
        kids = getattr(node, "children", ()) or ()
        if any(isinstance(k, widgets.ToggleButton) and k.description == "🐞"
               for k in kids):
            return node.children
    raise AssertionError("no beetle in the agent tab")


def test_the_beetle_stands_alone_until_it_is_pressed(bug_group):
    """The point of the change: nothing but the beetle takes room in the
    driving row until somebody wants to file something."""
    toggle, note, send, where = bug_group
    assert isinstance(toggle, widgets.ToggleButton)
    assert [w.layout.display for w in (note, send, where)] == ["none"] * 3


def test_pressing_it_opens_the_line_and_pressing_again_puts_it_away(bug_group):
    toggle, note, send, where = bug_group

    toggle.value = True
    assert [w.layout.display for w in (note, send, where)] == ["flex"] * 3

    toggle.value = False
    assert [w.layout.display for w in (note, send, where)] == ["none"] * 3


def test_opening_says_what_is_held_and_where_a_press_would_put_it(bug_group):
    """A button that could equally be a start switch has to say it is not
    one, and the amount already held is the smallest thing that shows the
    difference."""
    toggle, _note, _send, where = bug_group
    assert where.value == ""

    toggle.value = True
    text = where.value
    assert "Kept as you work" in text
    assert "message(s)" in text
    # Where it goes is named, one way or the other — a destination, or
    # the fact that there is no further one.
    assert "→ <code>" in text


def test_the_tooltip_says_which_way_round_the_button_works(bug_group):
    toggle = bug_group[0]
    assert "does not start" in (toggle.tooltip or "")


def test_the_note_and_the_send_button_are_still_the_ones_that_work(bug_group):
    """Folding is presentation. The control it folds has to remain the
    control that files the report — a beetle that opens an inert line
    would look identical and do nothing."""
    _toggle, note, send, _where = bug_group
    assert isinstance(note, widgets.Text)
    assert isinstance(send, widgets.Button)
    assert send.description == "Send"
    # `on_click` handlers live in the widget's callback registry; an
    # unwired button has none, which is exactly the failure above.
    assert send._click_handlers.callbacks
