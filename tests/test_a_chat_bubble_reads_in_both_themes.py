"""A chat bubble's text is readable whatever theme the page runs in.

Input: the agent tab's stylesheet (``_AGENT_CSS``). Output: the rules the
bubbles are drawn with.

Reproduction (real page, ``delfin-voila --dark``, 2026-10-15): the agent
and user bubbles have fixed light backgrounds (#f3f4f6, #dbeafe) but no
text colour, so they inherit the dark theme's white -- measured computed
colour rgb(255, 255, 255) on rgb(243, 244, 246), a contrast of about 1.1:1.
And no rule sets a line height, so a bubble inherits the widget label's
28 px: every line of an answer, and every list item, sits a blank line's
height apart from the next.
"""

from __future__ import annotations

import re

from delfin.dashboard import tab_agent as T


def _rule(selector: str) -> str:
    """The declarations of the FIRST rule whose selector is exactly this."""
    for m in re.finditer(r"([^{}]+)\{([^}]*)\}", T._AGENT_CSS):
        sels = [s.strip() for s in m.group(1).split(",")]
        if selector in [s.splitlines()[-1].strip() for s in sels]:
            return m.group(2)
    raise AssertionError(f"no rule for {selector}")


def _channel(hex6: str) -> float:
    c = int(hex6, 16) / 255
    return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4


def _luminance(hex_color: str) -> float:
    h = hex_color.lstrip("#")
    r, g, b = (_channel(h[i:i + 2]) for i in (0, 2, 4))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _contrast(a: str, b: str) -> float:
    la, lb = sorted((_luminance(a), _luminance(b)), reverse=True)
    return (la + 0.05) / (lb + 0.05)


def _decl(body: str, prop: str) -> str | None:
    m = re.search(rf"(?<![-\w]){prop}\s*:\s*([^;]+)", body)
    return m.group(1).strip() if m else None


def test_the_answer_and_question_bubbles_set_their_own_text_colour():
    for sel in (".delfin-chat-agent", ".delfin-chat-user"):
        body = _rule(sel)
        bg, fg = _decl(body, "background"), _decl(body, "color")
        assert fg, f"{sel} inherits the theme's text colour over {bg}"
        assert _contrast(bg, fg) >= 7.0, (sel, bg, fg)


def test_a_bubble_has_its_own_line_height():
    lh = _decl(_rule(".delfin-chat-msg"), "line-height")
    assert lh, "the widget label's 28 px line height leaks into the bubble"
    assert float(re.match(r"[\d.]+", lh).group(0)) <= 1.7
