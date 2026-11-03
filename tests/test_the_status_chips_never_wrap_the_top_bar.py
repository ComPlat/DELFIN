"""The status chips in the header never decide whether the top bar wraps.

The agent's working line is mirrored into the header (``agent_status_html``)
and rewritten every second while a turn runs, with text of changing length.
It sat in the right-hand group of the top bar, which wraps rather than
shrinks: where the bar was nearly full, a longer text moved that group to a
second row and the whole agent tab below the header with it, and a shorter
one moved them back.

Measured 2026-11-03 in a real page (headless Chromium 151, DPR 1.5, a
scripted endpoint, the branch label "HEAD (detached)"), origin/main at
1130x720 through one turn: the header 54 -> 94 px when the line appeared,
the tab 40 px down; 1 s later the tab-top measurement shortened the chat
viewport by 40 px; at the end of the turn the tab went up, down and up
again by 40 px within 0.95 s as the line changed, and the viewport got its
40 px back in between.

Input: :func:`delfin.dashboard.header_ui.top_bar` with stand-in widgets.
Output asserted: the chips sit in the spacer between the two sides, whose
flex basis is 0 and whose minimum width is fixed, so their text never
counts when the bar decides on a second row; the chips shrink and cut their
text with an ellipsis instead of growing; a chip's full text is its
tooltip; the dashboard builds its bar with this function.
"""

from __future__ import annotations

import ast
import pathlib
import re

import pytest

pytest.importorskip("ipywidgets")

from delfin.dashboard import header_ui as HU  # noqa: E402

_INIT = (pathlib.Path(__file__).resolve().parents[1]
         / "delfin" / "dashboard" / "__init__.py")


def _classes(widget) -> set:
    return set(getattr(widget, "_dom_classes", ()) or ())


def _bar():
    import ipywidgets as widgets
    title = widgets.HTML("DELFIN")
    branch = widgets.HBox([widgets.Button(description="main")])
    busy = widgets.HTML("")
    status = widgets.HTML("")
    usage = widgets.HTML("n/a")
    logout = widgets.HTML("Log out")
    right = widgets.HBox([usage, logout])
    bar = HU.top_bar([title, branch], [busy, status], right)
    return bar, title, branch, busy, status, right


def test_the_chips_sit_in_the_spacer_not_in_the_right_hand_group():
    bar, title, branch, busy, status, right = _bar()
    assert "delfin-topbar" in _classes(bar)
    kids = list(bar.children)
    assert kids[:2] == [title, branch]
    spacer = kids[2]
    assert kids[3] is right and len(kids) == 4
    assert "delfin-topbar-spacer" in _classes(spacer)
    # Same order as before, at the right end of the spacer, next to the
    # group they used to start: one shown copy per chip, which follows it.
    shown = list(spacer.children)
    assert len(shown) == 2
    assert all("delfin-topbar-chip" in _classes(s) for s in shown)
    busy.value = "<b>busy</b>"
    status.value = "<i>status</i>"
    assert "<b>busy</b>" in shown[0].value
    assert "<i>status</i>" in shown[1].value
    assert busy not in right.children and status not in right.children
    assert "delfin-group" in _classes(right)


def test_the_spacer_never_counts_its_content_in_the_wrap():
    """The wrap is decided on each item's base size clamped by its minimum.
    A basis of 0 and a fixed minimum make that size the same whatever the
    chips say. The minimum keeps the start of a chip, and with it its
    tooltip, on screen when the rest of the row is full."""
    bar, *_ = _bar()
    spacer = bar.children[2]
    assert spacer.layout.flex == "1 1 0"
    assert spacer.layout.min_width == "2em"
    assert spacer.layout.justify_content == "flex-end"


def test_a_cut_chip_keeps_its_full_text_in_the_tooltip():
    """The tooltip is in the shown markup: ipywidgets puts an HTML widget's
    ``tooltip`` on its description label, hidden while that is empty."""
    bar, _title, _branch, busy, status, _right = _bar()
    shown_busy, shown_status = bar.children[2].children
    line = ('<span style="font-family:monospace; color:#60a5fa;">'
            '⏳ waiting for kit.glm &amp; co … 12 s</span>')
    status.value = line
    assert line in shown_status.value
    assert ('title="⏳ waiting for kit.glm &amp; co … 12 s"'
            in shown_status.value)
    # What the writers wrote stays as they wrote it.
    assert status.value == line
    status.value = ""
    assert shown_status.value == ""
    # Nothing to read, nothing added: the spinner keeps its own title.
    spinner = "<span class='delfin-busy' title='Working'></span>"
    busy.value = spinner
    assert shown_busy.value == spinner


def test_the_chip_rules_cut_rather_than_grow():
    css = " ".join(HU.HEADER_CSS.split())
    rule = re.search(r"\.delfin-topbar-chip\s*\{([^}]*)\}", css)
    assert rule, "no rule for the chips"
    body = rule.group(1)
    for decl in ("min-width: 0", "overflow: hidden", "flex: 0 1 auto"):
        assert decl in body, decl
    content = re.search(
        r"\.delfin-topbar-chip\s*>\s*\.widget-html-content\s*\{([^}]*)\}", css)
    assert content, "no rule for the chip's text"
    for decl in ("white-space: nowrap", "text-overflow: ellipsis",
                 "overflow: hidden", "min-width: 0"):
        assert decl in content.group(1), decl


def test_the_dashboard_builds_its_bar_with_top_bar():
    tree = ast.parse(_INIT.read_text(encoding="utf-8"))
    calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
             and isinstance(n.func, ast.Attribute) and n.func.attr == "top_bar"]
    assert len(calls) == 1, "the dashboard must build its bar with top_bar"
    chips = ast.unparse(calls[0].args[1])
    assert "busy_indicator" in chips and "ctx.agent_status_html" in chips
    # And nothing else puts the agent's line into the right-hand group.
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "HBox" and node.args):
            assert "agent_status_html" not in ast.unparse(node.args[0])
