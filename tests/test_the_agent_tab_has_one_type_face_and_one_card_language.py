"""The chat speaks one type face and draws every card the same way.

Measured in a real page at 1920x1080 (2026-10-22): the transcript mixed the
Jupyter UI font, 13 px text and monospace system lines; the answer bubble
was a flat grey box without an edge, the same family of grey as the tool
lines around it; file, plot and molecule cards had 6 px corners, the
bubbles 8 px, the input a square box with a monospace placeholder.

Setting the new face on the whole tab cut the control labels ("Provider:"
to "Provid...", the shell button to ">_ ..."): the controls are sized to
the Jupyter face, so the face is scoped to the chat and the input.
"""

from __future__ import annotations

import re

from delfin.dashboard import tab_agent


def _rules(selector: str) -> list[str]:
    """Every declaration block whose selector list names ``selector``."""
    css = re.sub(r"/\*.*?\*/", "", tab_agent._AGENT_CSS, flags=re.S)
    css = " ".join(css.split())
    found = []
    for head, body in re.findall(r"([^{}]+)\{([^{}]*)\}", css):
        if selector in [s.strip() for s in head.split(",")]:
            found.append(body)
    return found


def _decl(selector: str, prop: str) -> str | None:
    """The last value ``selector`` gives ``prop`` (the cascade's winner)."""
    value = None
    for body in _rules(selector):
        for decl in body.split(";"):
            name, _, val = decl.partition(":")
            if name.strip() == prop:
                value = val.strip()
    return value


def test_the_chat_and_the_input_share_one_sans_serif_face():
    stack = _decl(".delfin-agent-root", "--delfin-ui-font")
    assert stack and stack.endswith("sans-serif")
    for selector in (".delfin-agent-root .delfin-agent-chat-host",
                     ".delfin-agent-root .delfin-agent-input textarea",
                     ".delfin-chat-system"):
        assert _decl(selector, "font-family") == "var(--delfin-ui-font)", \
            selector


def test_the_controls_keep_the_face_they_are_sized_to():
    css = " ".join(re.sub(r"/\*.*?\*/", "", tab_agent._AGENT_CSS,
                          flags=re.S).split())
    heads = [head.strip() for head, body
             in re.findall(r"([^{}]+)\{([^{}]*)\}", css)
             if "var(--delfin-ui-font)" in body]
    allowed = {".delfin-agent-root .delfin-agent-chat-host",
               ".delfin-agent-root .delfin-agent-input textarea",
               ".delfin-chat-system"}
    for head in heads:
        assert {s.strip() for s in head.split(",")} <= allowed, head


def test_an_answer_is_a_card_set_apart_from_the_question():
    assert _decl(".delfin-chat-agent", "border").startswith("1px solid")
    assert _decl(".delfin-chat-agent", "box-shadow")
    user = _decl(".delfin-chat-user", "background")
    agent = _decl(".delfin-chat-agent", "background")
    assert user and agent and user != agent
    assert _decl(".delfin-chat-msg", "border-radius") == "12px"


def test_every_card_in_the_transcript_has_the_same_corner():
    assert _decl(".delfin-agent-chat-host .delfin-artifact",
                 "border-radius") == "12px !important"
    assert _decl(".delfin-agent-root > .delfin-agent-chat-frame",
                 "border-radius") == "12px"
    assert _decl(".delfin-agent-root .delfin-agent-input textarea",
                 "border-radius") == "12px"


def test_the_input_shows_where_the_cursor_is():
    focus = ".delfin-agent-root .delfin-agent-input textarea:focus"
    assert _decl(focus, "box-shadow")
    assert _decl(focus, "border-color")
