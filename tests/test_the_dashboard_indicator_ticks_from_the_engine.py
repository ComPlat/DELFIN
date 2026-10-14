"""The agent tab's working line ticks from the engine's state, not its output.

Input: the real agent tab (``tab_agent.create_tab``) with the KIT
client's ``stream_message`` replaced by a script: a first token five
seconds late, a tool call that takes three seconds, an approval the
gate asks for. Output: the text of the working line under the chat,
read off the widget tree at fixed moments during the turn.

Asserted: during the silent first five seconds the line says it waits
for the model and its seconds advance (the label read "Thinking..." until
20 s had passed, then changed every 10 s); while the tool runs it names
the tool, its argument and a clock (it read "[1] $ sleep 3..." with no
clock, and kept saying so through the next request); while the gate asks,
it says the turn waits for the user's approval; after the turn it is gone.
"""

from __future__ import annotations

import html
import json
import re
import threading
import time

import pytest

pytest.importorskip("ipywidgets")

_SID = "b1c2d3e4f5a6b7c8d9e0f1a2b3c4d5e6"
MODEL = "kit.deepseek-v4-flash"


def _all(tab):
    import ipywidgets as widgets
    found, html_widgets, stack = {}, [], [tab]
    while stack:
        node = stack.pop()
        if isinstance(node, (widgets.Dropdown, widgets.Button,
                             widgets.Textarea)):
            found.setdefault((type(node).__name__, node.description), node)
        if isinstance(node, widgets.HTML):
            html_widgets.append(node)
        stack.extend(getattr(node, "children", ()) or ())
    return found, html_widgets


def _working_text(html_widgets) -> str:
    for w in html_widgets:
        value = str(w.value or "")
        if "delfin-agent-working" in value:
            return html.unescape(re.sub(r"<[^>]+>", " ", value)).strip()
    return ""


@pytest.fixture
def tab(tmp_path, monkeypatch):
    from delfin import user_settings
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "sk-test")
    path = user_settings.get_settings_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"agent": {
        "provider": "kit", "model": MODEL,
        "permission_profile": "all_free", "probe_endpoint": False}}))
    for name in ("calc", "archive", "office", "ws"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    ctx.agent_workspace = str(tmp_path / "ws")
    ctx.initial_session_id = _SID
    ctx.sessions_managed = True
    ctx.presence_key = "k1"
    widget, refs = tab_agent.create_tab(ctx)
    found, html_widgets = _all(widget)
    try:
        yield refs, found, html_widgets
    finally:
        refs["shutdown"](save=False)
        from delfin.agent import session_store as ss
        ss.release_session_lock(_SID)


def _script(monkeypatch, body):
    from delfin.agent import api_client as A
    monkeypatch.setattr(A.OpenAIClient, "stream_message", body)


def _send(found, text="go"):
    found[("Textarea", "")].value = text
    found[("Button", "Send")].click()


def _sample(html_widgets, times, t0):
    out = []
    for at in times:
        time.sleep(max(0.0, t0 + at - time.monotonic()))
        out.append(_working_text(html_widgets))
    return out


def _secs(text, phrase):
    m = re.search(re.escape(phrase) + r"[^…]*… (\d+) s", text)
    return int(m.group(1)) if m else None


def _ended(refs, timeout=20):
    deadline = time.monotonic() + timeout
    while refs["state"].get("streaming") and time.monotonic() < deadline:
        time.sleep(0.1)
    time.sleep(1.5)


def test_a_late_first_token_is_counted_from_the_send(tab, monkeypatch):
    from types import SimpleNamespace as NS
    refs, found, html_widgets = tab

    def _late(self, system="", messages=None, **kw):
        time.sleep(5.0)
        yield NS(type="text_delta", text="hello")
        yield NS(type="message_delta", input_tokens=1, output_tokens=1,
                 cost_usd=0.0, stop_reason="end_turn", text="")

    _script(monkeypatch, _late)
    t0 = time.monotonic()
    _send(found)
    early, late = _sample(html_widgets, (1.8, 4.2), t0)
    a, b = _secs(early, f"waiting for {MODEL}"), _secs(late, f"waiting for {MODEL}")
    assert a is not None and b is not None, (early, late)
    assert b > a, (early, late)
    _ended(refs)
    assert _working_text(html_widgets) == "", "the line outlived the turn"


def test_a_running_tool_has_its_name_argument_and_clock(tab, monkeypatch):
    from types import SimpleNamespace as NS
    refs, found, html_widgets = tab

    def _tool(self, system="", messages=None, **kw):
        yield NS(type="tool_use", tool_name="mcp__kit-coding__bash",
                 tool_input='{"command": "sleep 3"}', text="")
        time.sleep(3.0)
        yield NS(type="tool_result", tool_name="mcp__kit-coding__bash",
                 tool_output='{"exit_code": 0}', text="")
        time.sleep(2.5)
        yield NS(type="text_delta", text="done")
        yield NS(type="message_delta", input_tokens=1, output_tokens=1,
                 cost_usd=0.0, stop_reason="end_turn", text="")

    _script(monkeypatch, _tool)
    t0 = time.monotonic()
    _send(found)
    rows = _sample(html_widgets, (1.2, 2.6, 4.6), t0)
    secs = [_secs(r, "running bash: sleep 3") for r in rows[:2]]
    assert None not in secs and secs[1] > secs[0], rows
    assert f"waiting for {MODEL}" in rows[2], (
        "after the result the line goes back to the model: " + rows[2])
    _ended(refs)


def test_an_approval_reads_as_waiting_for_you(tab, monkeypatch):
    from types import SimpleNamespace as NS
    refs, found, html_widgets = tab
    asked = threading.Event()

    def _gated(self, system="", messages=None, **kw):
        yield NS(type="tool_use", tool_name="mcp__kit-coding__bash",
                 tool_input='{"command": "rm x"}', text="")
        broker = refs["state"].get("_kit_confirm_broker")
        asked.set()
        broker.callback("bash", {"command": "rm x"}, "rm x")
        yield NS(type="tool_result", tool_name="mcp__kit-coding__bash",
                 tool_output='{"error": "refused"}', text="")
        yield NS(type="text_delta", text="ok")
        yield NS(type="message_delta", input_tokens=1, output_tokens=1,
                 cost_usd=0.0, stop_reason="end_turn", text="")

    _script(monkeypatch, _gated)
    _send(found)
    assert asked.wait(10)
    time.sleep(1.5)
    text = _working_text(html_widgets)
    assert "waiting for your approval" in text, text
    broker = refs["state"]["_kit_confirm_broker"]
    for req in broker.pending_requests():
        broker.resolve(req, False)
    _ended(refs)
    assert _working_text(html_widgets) == ""
