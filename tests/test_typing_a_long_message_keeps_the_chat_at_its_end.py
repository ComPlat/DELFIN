"""Typing a message of several lines keeps the chat at its end.

The message box grows with its text (``grow`` in the agent tab's page
script): it writes ``height: auto``, reads ``scrollHeight`` and writes the
final height. The read forces a layout with the box at its 80 px floor, so
the chat frame above it is taller by the rest of the box for that layout,
and the browser clamps a transcript scrolled to its end to the taller box.
The final height does not give the offset back.

Measured 2026-11-03 in a real page (headless Chromium 151, 1920x1080, a
scripted endpoint), the reader at the end of the chat, typing a seven-line
draft (Shift+Enter) and then letters and spaces for 20 s:

  origin/main       painted gaps to the end of 20, 40, 60, 12 and 102 px;
                    the chat stayed 102 px short of its end until the box
                    was cleared, and following was switched off (the clamp
                    of 79 px, within 300 ms of a space, read as the reader
                    scrolling away; the "Newest" control appeared)
  with this fix     every painted frame at the end (gap 0), following on
                    throughout; what still moves is the transcript rising
                    by a line when the box grows by a line

The page script is JavaScript in a Python string and there is no JS runtime
among DELFIN's dependencies, so its order of operations is asserted: the
chat offsets are read before the box is collapsed and put back after its
final height is written, in the same function.
"""

from __future__ import annotations

import ast
import pathlib

_TAB_AGENT = (pathlib.Path(__file__).resolve().parents[1]
              / "delfin" / "dashboard" / "tab_agent.py")


def _grow_js() -> str:
    """The body of ``grow`` as the browser receives it."""
    src = _TAB_AGENT.read_text(encoding="utf-8")
    blocks = [
        node.value for node in ast.walk(ast.parse(src))
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
        and "window.__delfinInputGrow = true" in node.value
    ]
    assert len(blocks) == 1, "expected one script to grow the message box"
    js = blocks[0]
    start = js.index("function grow(ta)")
    end = js.index("function inputTA()", start)
    return js[start:end]


def test_the_chat_offsets_are_put_back_after_the_box_is_measured():
    body = _grow_js()
    kept = body.index("querySelectorAll('.delfin-agent-chat-host')")
    collapse = body.index("ta.style.height = 'auto'")
    final = body.index("ta.style.height = Math.min(")
    put_back = body.index("p[0].scrollTop = p[1]")
    assert kept < collapse < final < put_back, (
        "the chat offsets must be read before the box collapses and put "
        "back after its final height, or a chat at its end is left short "
        "of it by the collapse")


def test_only_an_offset_the_measurement_lowered_is_put_back():
    """A box that really grew leaves the chat shorter for good; the chat
    then follows its end through its resize observer. Putting back a
    higher offset would undo that."""
    body = _grow_js()
    assert "if (p[0].scrollTop < p[1]) p[0].scrollTop = p[1];" in body
    # Hidden sessions keep their own transcript in the page; a host
    # without a viewport has nothing to keep.
    assert "if (h.clientHeight) kept.push([h, h.scrollTop]);" in body


def test_the_box_still_grows_to_its_text_up_to_the_cap():
    body = _grow_js()
    assert "var edge = ta.offsetHeight - ta.clientHeight;" in body
    assert "Math.min(ta.scrollHeight + edge, MAXH)" in body
