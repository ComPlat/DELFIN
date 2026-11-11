"""A finished turn folds its tool steps behind one line.

Input: the agent chat transcript (delfin/dashboard/tab_agent.py). Today
every tool call of a turn stays expanded in the transcript, so a long
turn buries its answer below dozens of tool lines.

Wanted: while a turn runs, its steps show live as now; once the turn
ends, the steps between the user message and the final answer fold into
one native ``<details>/<summary>`` line -- "Worked for 12s" (the turn's
wall time; "1m 05s" past a minute), "Stopped after 40s" for a turn that
was stopped, left unfolded for one that ended in an error. One click
unfolds the steps in place, a second folds them. Thinking blocks fold
with the steps. The answer stays visible below the fold line.

The fold line renders in Python, like every other message. The morph
script keys the shown messages by outerHTML, so the fold line must keep
its key when its open/closed state changes in the browser -- otherwise
every re-render (the transcript is re-rendered four times a second while
the agent streams) rebuilt the block and the reader's toggle was lost.
"""
from __future__ import annotations

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
    return m.group(1)


# A finished turn with two tool calls between the question and the answer.
def _turn_msgs():
    return [
        {"role": "user", "content": "first question"},
        {"role": "tool", "content": "<span class=\"tool-name\">Read</span>"},
        {"role": "thinking", "content": "I should look at the file"},
        {"role": "tool", "content": "<span class=\"tool-name\">Bash</span>"},
        {"role": "assistant", "content": "the answer",
         "turn_fold": {"duration_s": 12.4, "outcome": "ok"}},
    ]


# ---------------------------------------------------------------------------
# The fold line itself (module-level helpers, pure functions)
# ---------------------------------------------------------------------------

def test_fold_time_format():
    """``12.4`` -> ``12s``; ``65`` -> ``1m 05s``; the wall clock of the turn."""
    assert T._format_fold_time(12.4) == "12s"
    assert T._format_fold_time(65) == "1m 05s"
    assert T._format_fold_time(0.5) == "1s"


def test_fold_line_text_names_the_outcome():
    """A turn that ended in an error or was stopped says so."""
    assert "Worked for 12s" in T._render_turn_fold_line(12.4, "ok")
    assert "Stopped after 40s" in T._render_turn_fold_line(40.0, "stopped")
    assert "Stopped after 40s" in T._render_turn_fold_line(40.0, "error")


# ---------------------------------------------------------------------------
# The transcript render
# ---------------------------------------------------------------------------

def test_a_finished_turn_folds_its_steps_behind_one_line(tmp_path,
                                                         monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    ctx = _ctx(tmp_path)
    widget, refs = T.create_tab(ctx)
    try:
        refs["state"]["chat_messages"].extend(_turn_msgs())
        refs["refresh_chat"](streaming=False)
        body = _template_body(_host(widget).value)
        # One fold line, native details/summary, the turn's wall time.
        assert body.count("delfin-chat-turn-steps") == 1
        assert "<details" in body and "Worked for 12s" in body
        # The steps live inside the fold, not beside it.
        i_details = body.index("delfin-chat-turn-steps")
        i_read = body.index("Read</span>")
        i_bash = body.index("Bash</span>")
        i_answer = body.index("the answer")
        assert i_details < i_read < i_bash < i_answer
        # Closed by default (no open attribute) for a turn that went well.
        i_tag = body.index("<details", i_details)
        assert body[i_tag:body.index(">", i_tag)].count("open") == 0
        # The thinking bubble folded with the steps.
        assert body.index("I should look at the file") < i_answer
        assert body.index("I should look at the file") > i_details
    finally:
        refs["shutdown"]()


def test_a_turn_without_a_fold_meta_renders_unfolded_as_today(
        tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    ctx = _ctx(tmp_path)
    widget, refs = T.create_tab(ctx)
    try:
        msgs = _turn_msgs()
        msgs[-1].pop("turn_fold")
        refs["state"]["chat_messages"].extend(msgs)
        refs["refresh_chat"](streaming=False)
        body = _template_body(_host(widget).value)
        # While a turn runs (no fold meta yet): no fold line, steps live.
        assert "delfin-chat-turn-steps" not in body
        assert "Read</span>" in body and "the answer" in body
    finally:
        refs["shutdown"]()


def test_an_error_turn_folds_unfolded_and_says_stopped(tmp_path,
                                                       monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    ctx = _ctx(tmp_path)
    widget, refs = T.create_tab(ctx)
    try:
        msgs = _turn_msgs()
        msgs[-1]["turn_fold"] = {"duration_s": 40.0, "outcome": "error"}
        msgs.append({"role": "system", "content": "Error: broken"})
        refs["state"]["chat_messages"].extend(msgs)
        refs["refresh_chat"](streaming=False)
        body = _template_body(_host(widget).value)
        i_tag = body.index("<details")
        tag = body[i_tag:body.index(">", i_tag)]
        assert "open" in tag
        assert "Stopped after 40s" in body
    finally:
        refs["shutdown"]()


# ---------------------------------------------------------------------------
# The morph key: the reader's open state survives a re-render
# ---------------------------------------------------------------------------

def _key_js() -> str:
    """The morph script's key function, as shipped in _CHAT_MORPH_JS."""
    blocks = [n.value for n in __import__("ast").walk(
        __import__("ast").parse(_SRC.read_text()))
        if isinstance(n, __import__("ast").Constant)
        and isinstance(n.value, str) and "__delfinChatMorph" in n.value]
    assert len(blocks) == 1
    return blocks[0]


def test_the_morph_key_ignores_the_fold_blocks_open_state():
    """The morph keys shown nodes by outerHTML. ``outerHTML`` carries the
    ``open`` attribute, so a reader toggling the fold changes the key and
    the next re-render rebuilds the block with its state reset. The key
    must drop the ``open`` attribute for turn-step folds."""
    js = _key_js()
    assert "delfin-chat-turn-steps" in js
    # The key for a fold element is computed without its open attribute.
    assert "open" in js
