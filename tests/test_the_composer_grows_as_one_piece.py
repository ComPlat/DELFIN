"""The message box and its buttons grow as one piece, and scroll only past the cap.

Input: the agent tab's composer (message box, shell and Files buttons,
Send) and the script that grows the box. Measured in the real page with a
five-line message (2026-10-22): the box grew to 98 px while its content
needed 98 px and only 96 px were inside its 1 px borders, so a grey
native scrollbar with arrow buttons appeared for 2 px; the buttons beside
it stayed at 80 px; focus drew a hard blue border.

After:
- the box's height is its content plus its borders, so a message up to
  the cap (280 px) shows no scrollbar; past the cap the box stops growing
  and scrolls, with a thin rounded scrollbar and no arrow buttons;
- the shell, Files and Send buttons stretch to the box's height (80 px
  at least);
- focus is a soft ring around a light blue border.
"""
from __future__ import annotations

import re
from pathlib import Path


def _css() -> str:
    from delfin.dashboard import tab_agent
    return " ".join(re.sub(r"/\*.*?\*/", "", tab_agent._AGENT_CSS,
                           flags=re.S).split())


def _rule(selector: str) -> str:
    out = ""
    for head, body in re.findall(r"([^{}]+)\{([^{}]*)\}", _css()):
        if selector in [s.strip() for s in head.split(",")]:
            out += body + ";"
    return out


def test_the_box_height_counts_its_borders():
    from delfin.dashboard import tab_agent
    src = Path(tab_agent.__file__).read_text(encoding="utf-8")
    grow = src[src.index("// --- Auto-grow the message box"):]
    grow = grow[:grow.index("})();")]
    assert "ta.offsetHeight - ta.clientHeight" in grow
    assert "Math.min(ta.scrollHeight + edge, MAXH)" in grow


def test_past_the_cap_it_scrolls_with_a_thin_scrollbar():
    bar = _rule(".delfin-agent-root .delfin-agent-input "
                "textarea::-webkit-scrollbar")
    assert "width: 8px" in bar
    assert "border-radius" in _rule(
        ".delfin-agent-root .delfin-agent-input "
        "textarea::-webkit-scrollbar-thumb")
    assert "display: none" in _rule(
        ".delfin-agent-root .delfin-agent-input "
        "textarea::-webkit-scrollbar-button")


def test_the_buttons_stretch_with_the_box():
    assert "align-items: stretch !important" in _rule(
        ".delfin-agent-root .delfin-agent-input-area > .widget-hbox")
    assert "align-items: stretch !important" in _rule(
        ".delfin-agent-root .delfin-agent-send-row")
    buttons = _rule(".delfin-agent-root .delfin-agent-input-area "
                    ".jupyter-button")
    assert "height: auto !important" in buttons
    assert "min-height: 80px" in buttons


def test_focus_is_a_soft_ring():
    focus = _rule(".delfin-agent-root .delfin-agent-input textarea:focus")
    assert "border-color: #93c5fd" in focus
    assert "0 0 0 4px" in focus
