"""One card per delegated task, in the parent's chat.

The chat showed a delegation today as a grey tool line plus a side panel;
each delegate is now one card in the transcript at the point it started:
title, model, a state icon (running / done / failed / stopped), elapsed
time, and the result summary behind "more" when done. A running card
updates in place — no new card per event — and "Open" / "Stop" on the
card go to the existing drill-in view and the existing stop path.

The state carried by the card is the SAME dict the existing subagent
pane already fills (state["subagent_calls"]), so the Details-fold call
log and the side panel keep their current behaviour untouched.
"""

from __future__ import annotations

import json
import time
from types import SimpleNamespace


def _widgets(tab):
    import ipywidgets as widgets
    found, stack = {}, [tab]
    while stack:
        node = stack.pop()
        if isinstance(node, (widgets.Button, widgets.Textarea)):
            found.setdefault((type(node).__name__, node.description), node)
        stack.extend(getattr(node, "children", ()) or ())
    return found


def _make_tab(tmp_path, monkeypatch):
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
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script, **kw: None
    ctx.agent_workspace = str(tmp_path / "ws")
    tab, refs = tab_agent.create_tab(ctx)
    return tab, refs, refs["state"]


# ---------------------------------------------------------------------------
# The card renderer (module-level, pure)
# ---------------------------------------------------------------------------

def _call(**over):
    call = {
        "subagent_type": "explore",
        "description": "map the sandbox layout",
        "model": "kit.deepseek-v4-flash",
        "sa_id": "sa-1",
        "status": "running",
        "started_at": time.time() - 12.0,
        "elapsed_s": None,
        "output": "",
    }
    call.update(over)
    return call


def test_card_renders_title_type_state_and_elapsed():
    from delfin.dashboard.tab_agent import _render_delegate_card_html
    html = _render_delegate_card_html(_call())
    assert "map the sandbox layout" in html
    assert "explore" in html
    assert "kit.deepseek-v4-flash" in html
    assert "12s" in html                    # elapsed time, not a timestamp
    # state icon: running is marked as such
    assert "delfin-delegate-running" in html


def test_card_carries_open_and_stop_for_a_running_delegate():
    from delfin.dashboard.tab_agent import _render_delegate_card_html
    html = _render_delegate_card_html(_call())
    # The card speaks to the page script through data attributes, the way
    # the chat's jump button does; no inline JS beyond the one bridge.
    assert 'data-delfin-delegate-open="sa-1"' in html
    assert 'data-delfin-delegate-stop="sa-1"' in html


def test_done_card_shows_result_summary_first_lines_behind_more():
    from delfin.dashboard.tab_agent import _render_delegate_card_html
    report = "line one\nline two\nline three\nline four"
    html = _render_delegate_card_html(
        _call(status="done", elapsed_s=41.2, output=json.dumps(
            {"sa_id": "sa-1", "result": report, "error": ""})))
    assert "delfin-delegate-done" in html
    assert "line one" in html                       # first lines shown
    assert "41s" in html
    assert "more" in html                           # the rest behind more
    assert "<details" in html


def test_failed_and_stopped_cards_show_their_state():
    from delfin.dashboard.tab_agent import _render_delegate_card_html
    failed = _render_delegate_card_html(
        _call(status="failed", elapsed_s=3.0,
              output=json.dumps({"sa_id": "sa-1", "result": "",
                                 "error": "backend down"})))
    assert "delfin-delegate-failed" in failed
    assert "backend down" in failed
    stopped = _render_delegate_card_html(_call(status="stopped"))
    assert "delfin-delegate-stopped" in stopped


def test_card_stays_a_chat_message_bubble_of_the_chat_css():
    """The card reuses the chat's own bubble class and its 12 px corner;
    it adds no copied stylesheet of its own."""
    from delfin.dashboard import tab_agent as T
    html = T._render_delegate_card_html(_call())
    assert 'class="delfin-chat-msg delfin-delegate-card' in html
    assert ".delfin-delegate-card" in T._AGENT_CSS
    assert "border-radius: 12px" in T._AGENT_CSS


# ---------------------------------------------------------------------------
# Emit path: one card per delegate, updated in place
# ---------------------------------------------------------------------------

def test_one_card_in_the_transcript_per_delegate_and_updates_in_place(
        tmp_path, monkeypatch):
    tab, refs, state = _make_tab(tmp_path, monkeypatch)
    try:
        emit = refs["emit_delegate_card"]
        emit({"subagent_type": "explore",
              "description": "first task", "model": "m1", "sa_id": "sa-1",
              "status": "running", "started_at": time.time()})
        cards = [m for m in state["chat_messages"]
                 if m.get("role") == "delegate"]
        assert len(cards) == 1
        # A second event for the SAME delegate (a state refresh, a
        # heartbeat) must not append a second card.
        emit({"sa_id": "sa-1", "status": "running",
              "started_at": cards[0].get("started_at")})
        cards = [m for m in state["chat_messages"]
                 if m.get("role") == "delegate"]
        assert len(cards) == 1
        # A different delegate appends its own card, stacked after it.
        emit({"subagent_type": "plan",
              "description": "second task", "model": "m1", "sa_id": "sa-2",
              "status": "running", "started_at": time.time()})
        cards = [m for m in state["chat_messages"]
                 if m.get("role") == "delegate"]
        assert [c.get("sa_id") for c in cards] == ["sa-1", "sa-2"]
    finally:
        refs["shutdown"](save=False)


def test_tool_result_marks_the_right_delegate_done_by_id(
        tmp_path, monkeypatch):
    """The result carries sa_id; the card updated is that delegate's —
    not the first running one (two parallel delegates used to collide)."""
    tab, refs, state = _make_tab(tmp_path, monkeypatch)
    try:
        emit = refs["emit_delegate_card"]
        emit({"subagent_type": "explore", "description": "A",
              "model": "m1", "sa_id": "sa-A", "status": "running",
              "started_at": time.time()})
        emit({"subagent_type": "explore", "description": "B",
              "model": "m1", "sa_id": "sa-B", "status": "running",
              "started_at": time.time()})
        payload = json.dumps({"sa_id": "sa-B", "result": "B report",
                              "error": "", "elapsed_s": 9.5})
        on_result = refs["on_tool_result"]
        on_result("Agent", payload)
        cards = {c.get("sa_id"): c for c in state["chat_messages"]
                 if c.get("role") == "delegate"}
        assert cards["sa-B"]["status"] == "done"
        assert cards["sa-B"]["output"] == payload
        assert cards["sa-A"]["status"] == "running"
    finally:
        refs["shutdown"](save=False)


def test_failed_payload_marks_the_card_failed(tmp_path, monkeypatch):
    tab, refs, state = _make_tab(tmp_path, monkeypatch)
    try:
        emit = refs["emit_delegate_card"]
        emit({"subagent_type": "explore", "description": "A",
              "model": "m1", "sa_id": "sa-A", "status": "running",
              "started_at": time.time()})
        payload = json.dumps({"sa_id": "sa-A", "result": "",
                              "error": "backend down", "elapsed_s": 3.0})
        refs["on_tool_result"]("Agent", payload)
        cards = [c for c in state["chat_messages"]
                 if c.get("role") == "delegate"]
        assert cards[0]["status"] == "failed"
    finally:
        refs["shutdown"](save=False)


def test_the_agent_call_emits_a_running_card_in_the_chat(
        tmp_path, monkeypatch):
    """The call path of the chat (the 'Agent' tool event) emits the card —
    the transcript shows the card where the delegation was started."""
    tab, refs, state = _make_tab(tmp_path, monkeypatch)
    try:
        on_call = refs["on_tool_call"]
        on_call("Agent", json.dumps({
            "subagent_type": "explore",
            "description": "investigate the flake",
            "prompt": "find it"}))
        cards = [m for m in state["chat_messages"]
                 if m.get("role") == "delegate"]
        assert len(cards) == 1
        assert cards[0]["description"] == "investigate the flake"
        assert cards[0]["status"] == "running"
        # The Details-fold call log stays as it is: the pane state entry
        # still exists for the side panel.
        assert state["subagent_calls"], "side panel state lost"
    finally:
        refs["shutdown"](save=False)


def test_stopping_a_delegate_stops_the_running_one(tmp_path, monkeypatch):
    """Stop on the card goes to the same path as the panel's stop."""
    tab, refs, state = _make_tab(tmp_path, monkeypatch)
    calls = []
    monkeypatch.setattr(
        "delfin.agent.background_view.cancel",
        lambda ws, group, item_id: calls.append((group, item_id)) or "ok")
    try:
        emit = refs["emit_delegate_card"]
        emit({"subagent_type": "explore", "description": "A",
              "model": "m1", "sa_id": "sa-A", "status": "running",
              "started_at": time.time()})
        refs["stop_delegate_card"]("sa-A")
        assert calls == [("agents", "sa-A")], calls
        cards = [c for c in state["chat_messages"]
                 if c.get("role") == "delegate"]
        assert cards[0]["status"] == "stopped"
    finally:
        refs["shutdown"](save=False)
