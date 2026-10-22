"""The agent tab's bars are rounded, and the session divider spans the column.

Input: the agent tab's and the session shell's stylesheets, and the
divider's init script. Before (measured in the real page, 2026-10-22):
- the composer's shell and Files buttons, the Details bar, and the
  working line, task line and Background box at the foot of the chat had
  square corners, while the bubbles, cards, chat frame, input and Send had
  12 px corners: the strips ran edge to edge under a rounded frame;
- the divider between the session list and the chat had a height written
  by a script, measured from the shell's top to the Send row's bottom and
  re-measured only when the Send row or the window changed size. Its top
  sat 10 px below the column's and its bottom 10 px below the Send row;
  anything that moved the Send row without resizing it left the length
  stale, so it varied from one moment to the next.
After: the buttons and the Details bar take the 12 px corner, the foot
strips are inset 8 px with 8 px corners, and the divider stretches over
the whole column in CSS, from the control row to the Details bar, with no
height written by a script. The column fits the window since the
agent tab's height script, so the full column no longer means the full
transcript, which is what the height measure was written against.
"""
from __future__ import annotations

import re


def _rule(css: str, selector: str) -> str:
    css = " ".join(re.sub(r"/\*.*?\*/", "", css, flags=re.S).split())
    found = ""
    for head, body in re.findall(r"([^{}]+)\{([^{}]*)\}", css):
        if selector in [s.strip() for s in head.split(",")]:
            found += body + ";"
    return found


def test_the_divider_stretches_over_the_column():
    from delfin.dashboard import agent_sessions as AS
    host = _rule(AS._SIDEBAR_CSS, ".delfin-splitter-host")
    assert "align-self: stretch" in host and "sticky" not in host
    bar = _rule(AS._SIDEBAR_CSS, ".delfin-session-splitter")
    assert "flex: 1 1 auto" in bar and "100vh" not in bar
    assert "splitter.style.height" not in AS._SPLITTER_INIT_JS


def test_the_composer_buttons_and_details_bar_are_rounded():
    from delfin.dashboard import tab_agent
    css = tab_agent._AGENT_CSS
    assert "border-radius: 12px" in _rule(
        css, ".delfin-agent-root .delfin-agent-send-row .jupyter-button")
    assert "border-radius: 12px" in _rule(
        css, ".delfin-agent-root .delfin-agent-below > .jupyter-widget-Collapse "
             "> .jupyter-widget-Collapse-header")


def test_the_foot_strips_are_inset_with_rounded_corners():
    from delfin.dashboard import tab_agent
    css = tab_agent._AGENT_CSS
    for sel in (".delfin-agent-chat-frame .delfin-agent-working",
                ".delfin-agent-chat-frame > .delfin-agent-task-strip",
                ".delfin-agent-chat-frame > .delfin-agent-background"):
        body = _rule(css, sel)
        assert "border-radius: 8px" in body, sel
        assert "margin: 0 8px 6px" in body, sel
