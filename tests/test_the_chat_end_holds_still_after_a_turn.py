"""Stop goes out with the turn.

Measured per animation frame in a real page (2026-11-03, headless
Chromium, 1920x1080, a scripted endpoint): Stop stayed enabled and red
after every turn. The worker's cleanup updates the buttons from the turn's
own thread, and Stop stays enabled while that thread is alive -- which it
is, while it asks. Nothing asked again after the thread had exited, so the
header said "working" for good: Stop enabled 180 s after the answer, with
no turn running.
"""

from __future__ import annotations

import json
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
