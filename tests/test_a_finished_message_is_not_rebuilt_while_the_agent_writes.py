"""A finished message stays on the page while the agent writes the next one.

Input: the agent chat (delfin/dashboard/tab_agent.py). Python renders the
whole transcript as one HTML string and assigns it to ``chat_html.value``
about four times a second while the agent streams.

Before: ipywidgets puts every assignment in with ``innerHTML``, so every
message, card and plot on the page was destroyed and built again four
times a second. Measured in the real page (headless Chromium, butane
drive, three 3D cards, a reader scrolled up into them): the number of
canvases in the chat went 3 -> 0 -> 1 -> 2 -> 3 on every update
(timeline ``33301230123333123012333312301233``), so the molecules a reader
was looking at blinked; a text selection, a hover and a drag in progress
were lost with the node under them.

Now the value carries the transcript inside an inert ``<template>`` and a
small script after it. The script keeps the shown messages in a container
of its own beside the widget's content (which ipywidgets never touches),
and replaces only the messages from the first one whose HTML changed: a
streaming answer replaces itself, a new message is appended, and
everything above stays the same DOM nodes.

Rejected: one widget per message. It would change the DOM every
stylesheet rule and script of the chat is written against, for the same
effect.
"""
from __future__ import annotations

import ast
import pathlib
import re

import ipywidgets as widgets

from delfin.dashboard import tab_agent as T
from delfin.dashboard.context import DashboardContext

_SRC = (pathlib.Path(__file__).resolve().parents[1]
        / "delfin" / "dashboard" / "tab_agent.py")


def _ctx(tmp_path):
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           agent_dir=tmp_path / "agent_workspace",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    return ctx


def _walk(w):
    yield w
    for c in getattr(w, "children", ()) or ():
        yield from _walk(c)


def _host(widget):
    hosts = [w for w in _walk(widget) if isinstance(w, widgets.HTML)
             and "delfin-agent-chat-host" in (w._dom_classes or ())]
    assert len(hosts) == 1, hosts
    return hosts[0]


def _template_body(value: str) -> str:
    m = re.fullmatch(r'<template class="delfin-chat-src">(.*)</template>'
                     r'(<img [^>]*>)', value, flags=re.S)
    assert m, value[:200]
    assert "__delfinChatMorph" in m.group(2)
    return m.group(1)


def test_the_transcript_reaches_the_page_only_inside_an_inert_template():
    body = _template_body(T._chat_value('<div class="delfin-agent-chat">'
                                        '<div class="delfin-chat-msg">hi</div>'
                                        '</div>'))
    assert body.startswith('<div class="delfin-agent-chat">')
    assert "hi" in body


def test_text_that_spells_the_closing_tag_cannot_end_the_template():
    """A message that quotes ``</template>`` must not close the wrapper:
    what followed it would be parsed live and shown twice."""
    value = T._chat_value('<div class="delfin-agent-chat"><p>a </template> '
                          '<img src=x onerror=alert(1)> b</p></div>')
    body = _template_body(value)
    assert "</template" not in body.lower()
    assert value.lower().count("</template") == 1


def test_a_built_tab_shows_its_chat_through_the_template(tmp_path,
                                                         monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    ctx = _ctx(tmp_path)
    widget, refs = T.create_tab(ctx)
    try:
        host = _host(widget)
        assert "Start a conversation" in _template_body(host.value)
        refs["state"]["chat_messages"].extend([
            {"role": "user", "content": "first question"},
            {"role": "assistant", "content": "first answer"},
            {"role": "assistant", "content": "streaming"},
        ])
        refs["refresh_chat"](streaming=False)
        body = _template_body(host.value)
        assert "first question" in body and "first answer" in body
        refs["refresh_chat"](streaming=True)
        assert "streaming" in _template_body(host.value)
    finally:
        refs["shutdown"]()


def test_no_assignment_of_the_chat_bypasses_the_template():
    tree = ast.parse(_SRC.read_text(encoding="utf-8"))
    bad = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        for t in node.targets:
            if (isinstance(t, ast.Attribute) and t.attr == "value"
                    and isinstance(t.value, ast.Name)
                    and t.value.id == "chat_html"):
                v = node.value
                if not (isinstance(v, ast.Call)
                        and getattr(v.func, "id", "") == "_chat_value"):
                    bad.append(node.lineno)
    assert not bad, f"chat_html.value assigned without _chat_value: {bad}"


def test_the_shown_messages_are_laid_out_and_the_source_is_not():
    css = T._AGENT_CSS
    assert re.search(r"\.delfin-agent-chat-host\s*>\s*\.widget-html-content"
                     r"\s*\{[^}]*display:\s*none", css)
    assert re.search(r"\.delfin-agent-chat-host\s*>\s*\.delfin-agent-chat"
                     r"\s*\{[^}]*min-height:\s*100%", css)


def test_the_follow_logic_watches_the_shown_messages():
    """The ResizeObserver read the replaced content's size; that element
    is now empty, so late-loading images in the shown messages must be
    observed where they are."""
    blocks = [n.value for n in ast.walk(ast.parse(_SRC.read_text()))
              if isinstance(n, ast.Constant) and isinstance(n.value, str)
              and "__delfinChatScroll" in n.value]
    assert len(blocks) == 1
    js = blocks[0]
    assert "__delfinWatchChat" in js
    assert "> .delfin-agent-chat" in js


def _mol3d_js() -> str:
    blocks = [n.value for n in ast.walk(ast.parse(_SRC.read_text()))
              if isinstance(n, ast.Constant) and isinstance(n.value, str)
              and "__delfinMol3DInstalled" in n.value]
    assert len(blocks) == 1
    return blocks[0]


def test_the_wheel_over_a_molecule_scrolls_the_chat():
    """3Dmol takes every wheel event over its canvas for zoom and cancels
    it. Measured in the real page at 1920x1080 once the chat kept its cards
    (9d8f5cca): a 25-frame trajectory card fills the chat, and forty wheel
    steps down from inside it left the reader ~1000 px short of the end.
    Before that, every update rebuilt the card, so the canvas under the
    pointer was missing part of the time and the wheel fell through to the
    chat at random. A wheel
    over a molecule in the chat now scrolls the chat; it zooms with Ctrl
    (a trackpad pinch arrives as one) or after the molecule was clicked,
    until the pointer leaves it."""
    js = _mol3d_js()
    i = js.index("'wheel'")
    guard = js[max(0, i - 1500): i + 1500]
    assert "capture: true" in guard and "passive: true" in guard
    assert "ctrlKey" in guard
    assert "stopPropagation" in guard
    assert "preventDefault" not in guard
    assert ".delfin-agent-chat-host" in guard
    assert "'pointerdown'" in guard
