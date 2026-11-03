"""The end of the chat holds still, and Stop goes out with the turn.

Two defects, both measured per animation frame in a real page (2026-11-03,
headless Chromium, 1920x1080, a scripted endpoint):

1. Stop stayed enabled and red after every turn. The worker's cleanup
   updates the buttons from the turn's own thread, and Stop stays enabled
   while that thread is alive -- which it is, while it asks. Nothing asked
   again after the thread had exited, so the header said "working" for
   good: Stop enabled 180 s after the answer, with no turn running.

2. Every 1.5 s the chat host scrolled 12 px up from its end (scrollTop
   1511 -> 1499), idle or streaming, with the reader at the end only. The
   tab-top measurement writes the untrimmed top, reads scrollHeight (a
   forced layout with the tab 12 px taller, which clamps a host at its
   end), then writes the trimmed top. The chat's scroll listener put it
   back in the next frame unless the reader had used wheel or pointer in
   the last 300 ms; then 97 painted frames in 20 s idle carried the
   offset. The measurement now writes nothing while its input is
   unchanged, and puts hosts at their end back at the end when it does
   re-measure.

The first is driven through the real tab; the second is asserted on the
script's order of operations, which is the mechanism (there is no JS
runtime among DELFIN's dependencies, and a fake DOM has no layout to clamp).
"""

from __future__ import annotations

import json
import re
import time
from types import SimpleNamespace

import pytest


def _widgets(tab):
    import ipywidgets as widgets
    found, stack = {}, [tab]
    while stack:
        node = stack.pop()
        if isinstance(node, (widgets.Button, widgets.Textarea)):
            found.setdefault((type(node).__name__, node.description), node)
        stack.extend(getattr(node, "children", ()) or ())
    return found


def test_stop_goes_out_when_the_turn_has_answered(tmp_path, monkeypatch):
    pytest.importorskip("ipywidgets")
    from delfin import user_settings
    from delfin.agent import api_client as A
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "sk-test")
    path = user_settings.get_settings_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"agent": {
        "provider": "kit", "model": "kit.deepseek-v4-flash",
        "permission_profile": "all_free", "probe_endpoint": False}}))

    def _answers(self, system="", messages=None, **kw):
        yield SimpleNamespace(type="text_delta", text="an answer")
        yield SimpleNamespace(type="message_delta", input_tokens=10,
                              output_tokens=2, cost_usd=0.0,
                              stop_reason="end_turn", text="")

    monkeypatch.setattr(A.OpenAIClient, "stream_message", _answers)
    for name in ("calc", "archive", "office", "ws"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script, **kw: None
    ctx.agent_workspace = str(tmp_path / "ws")
    tab, refs = tab_agent.create_tab(ctx)
    state = refs["state"]
    w = _widgets(tab)
    stop = w[("Button", "Stop")]
    try:
        assert stop.disabled is True
        w[("Textarea", "")].value = "Hi"
        w[("Button", "Send")].click()
        deadline = time.monotonic() + 30
        while state.get("_turn_thread") is None and time.monotonic() < deadline:
            time.sleep(0.02)
        thread = state.get("_turn_thread")
        assert thread is not None, "the turn never started"
        thread.join(30)
        assert not thread.is_alive()
        assert any("an answer" in str(m.get("content"))
                   for m in state["chat_messages"]), "the turn did not answer"
        assert state["streaming"] is False
        # Nothing runs any more: Stop is greyed out and not red.
        assert stop.disabled is True, (
            "Stop is still enabled after the turn's thread has ended")
        assert stop.button_style == ""
    finally:
        refs["shutdown"](save=False)


def _measure_js():
    from delfin.dashboard.tab_agent import _MEASURE_AGENT_TOP_JS
    return _MEASURE_AGENT_TOP_JS


def _body_of_m(js):
    start = js.index("function m(){")
    end = js.index("window.addEventListener('resize',m)")
    return js[start:end]


def test_the_tab_top_is_sent_with_the_startup_scripts(tmp_path):
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    S._GLOBAL = S.Scheduler(path=tmp_path / "cron.json")
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script, **kw: None
    _tab, refs = tab_agent.create_tab(ctx)
    try:
        assert _measure_js() in ctx.init_js_parts
    finally:
        refs["shutdown"](save=False)


def test_an_unchanged_layout_writes_nothing():
    """The periodic tick returns before its first write when the tab's
    top, the window height and the overflow at the current value are what
    they were after the last measurement."""
    body = _body_of_m(_measure_js())
    guard = body.index("if(key===done)return;")
    writes = [m.start() for m in re.finditer(
        r"setProperty\('--delfin-agent-top'", body)]
    assert len(writes) == 2
    assert guard < writes[0]
    key = body[body.index("var key="):guard]
    assert "window.innerHeight" in key and "overflow(ss)" in key
    # The key is recorded from the layout AFTER the writes, or the next
    # tick would compare against the intermediate value and re-measure.
    assert body.index("done=") > writes[-1]
    assert "overflow(ss)" in body[body.index("done="):]


def test_a_remeasure_keeps_a_chat_at_its_end():
    """The intermediate (taller) layout clamps a host at its end; the hosts
    at their end are noted before the first write and put back after the
    last, in the same task, so no frame paints the clamped offset."""
    body = _body_of_m(_measure_js())
    writes = [m.start() for m in re.finditer(
        r"setProperty\('--delfin-agent-top'", body)]
    noted = body.index("querySelectorAll('.delfin-agent-chat-host')")
    put_back = body.index("ends[k].scrollTop=ends[k].scrollHeight")
    assert noted < writes[0] < writes[-1] < put_back
    # The interval is still there: a page that settles late is re-measured.
    assert "setInterval(m,1500)" in _measure_js()
