"""A turn that never answered leaves the history as it found it.

Input: a turn that started fresh -- the context was over the budget, so the
engine archived the history and kept only the fused prompt -- and then
ended before any answer text: a stop, an empty backend turn, or an
exception. Output: ``engine.messages`` holds the pre-turn history again.

The defect this pins (field report, dashboard, kit): a restored session of
about 470 rows, a switch to a model whose window it exceeded, a one-word
message, Stop before the first token. The fresh start had cut the history
down to the fused prompt; the stop branch then took that prompt back out,
so the engine held nothing. The next save wrote ``engine_messages: []``
beside the full chat transcript, and every reload after it showed the whole
conversation while the model started with none of it. The stop notice said
"Nothing was added to the history" -- the history had been removed.

The dashboard test drives that sequence through the real tab: restore,
model switch, Send, Stop, a second tab loading the saved file.
"""

from __future__ import annotations

import json
import threading
import time
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from delfin.agent.engine import AgentEngine

from tests.test_t2_engine_fresh_start import _pack_tree


def _history(n=30):
    return [{"role": "user" if i % 2 == 0 else "assistant",
             "content": f"turn {i} " + "x" * 40} for i in range(n)]


class _Client:
    """One turn's transport. ``how``: stop | empty | raise | answer."""

    def __init__(self, how):
        self.how = how
        self.engine = None
        self.sent = None

    def stream_message(self, *, messages, system="", max_tokens=0, **kw):
        self.sent = list(messages or [])
        if self.how == "stop":
            self.engine.request_stop()
            return
        if self.how == "raise":
            raise RuntimeError("connection reset before the first token")
        if self.how == "answer":
            yield SimpleNamespace(type="text_delta", text="an answer")
        yield SimpleNamespace(type="message_delta", input_tokens=10,
                              output_tokens=0, cost_usd=0.0,
                              stop_reason="end_turn", text="")


def _engine(tmp_path, how):
    tmp_path = _pack_tree(tmp_path)
    client = _Client(how)
    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=tmp_path, backend="cli", mode="quick",
                             pack_dir=tmp_path)
    engine.client = client
    client.engine = engine
    return engine, client


@pytest.fixture(autouse=True)
def _no_archive(monkeypatch):
    import delfin.agent.session_store as store
    monkeypatch.setattr(store, "archive_pre_compaction_transcript",
                        lambda *a, **k: None)


@pytest.mark.parametrize("how", ["stop", "empty", "raise"])
def test_a_fresh_turn_without_an_answer_restores_the_cut_history(tmp_path, how):
    engine, client = _engine(tmp_path, how)
    before = _history()
    engine.messages = [dict(m) for m in before]
    engine.start_fresh = "[Fresh context — full history cut at the token budget]"
    try:
        engine.stream_response(user_message="[Fresh context]\n\nHi")
    except Exception:
        pass
    # The turn really started fresh: one row went to the transport.
    assert client.sent is not None and len(client.sent) == 1
    assert [m["content"] for m in engine.messages] == [
        m["content"] for m in before], (how, len(engine.messages))
    assert engine.start_fresh is False


def test_a_fresh_turn_that_answered_keeps_the_fresh_start(tmp_path):
    """The restore is for a turn that produced nothing. A turn that answered
    ran against the fresh block, and that is the history from then on."""
    engine, _client = _engine(tmp_path, "answer")
    engine.messages = _history()
    engine.start_fresh = "[Fresh context — full history cut at the token budget]"
    engine.stream_response(user_message="[Fresh context]\n\nHi")
    assert [m["role"] for m in engine.messages] == ["user", "assistant"]


def test_a_later_unanswered_turn_does_not_bring_back_an_older_cut(tmp_path):
    """The cut belongs to the turn that made it: once a fresh turn has
    answered, a later stop takes back only its own message."""
    engine, client = _engine(tmp_path, "answer")
    engine.messages = _history()
    engine.start_fresh = "[Fresh context — full history cut at the token budget]"
    engine.stream_response(user_message="[Fresh context]\n\nHi")
    client.how = "stop"
    engine.stream_response(user_message="and then?")
    assert [m["role"] for m in engine.messages] == ["user", "assistant"]


# ---------------------------------------------------------------------------
# The dashboard: the field sequence through the real tab
# ---------------------------------------------------------------------------

_SID = "a0b1c2d3e4f5a6b7c8d9e0f1a2b3c4d5"


def _widgets(tab):
    import ipywidgets as widgets
    found, stack = {}, [tab]
    while stack:
        node = stack.pop()
        if isinstance(node, (widgets.Dropdown, widgets.Button,
                             widgets.Textarea)):
            found.setdefault((type(node).__name__, node.description), node)
        stack.extend(getattr(node, "children", ()) or ())
    return found


def _open(tmp_path, key):
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext
    for name in ("calc", "archive", "office", "ws"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    ctx.agent_workspace = str(tmp_path / "ws")
    ctx.initial_session_id = _SID
    ctx.sessions_managed = True
    ctx.presence_key = key
    tab, refs = tab_agent.create_tab(ctx)
    return refs, _widgets(tab)


def test_a_stop_after_a_model_switch_does_not_empty_the_saved_session(
        tmp_path, monkeypatch):
    pytest.importorskip("ipywidgets")
    from delfin.agent import api_client as A
    from delfin.agent import session_store as ss
    from delfin.agent import scheduler as S
    from delfin import user_settings

    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "sk-test")
    path = user_settings.get_settings_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"agent": {
        "provider": "kit", "model": "kit.deepseek-v4-flash",
        "permission_profile": "all_free", "probe_endpoint": False}}))

    # About 180k tokens at four characters each: inside a 524k window,
    # over the 131k of kit.glm-5.3.
    history = []
    for i in range(30):
        history += [{"role": "user", "content": f"task {i} " + "x" * 24000},
                    {"role": "assistant", "content": f"answer {i}"}]
    chat = [{"role": m["role"], "content": m["content"][:40], "role_label": ""}
            for m in history]
    ss.save_session(_SID, mode="solo", chat_messages=chat,
                    engine_messages=history, provider="kit",
                    model="kit.glm-5.3-flash", perm_profile="all_free",
                    workspace=str(tmp_path / "ws"))
    ss.release_session_lock(_SID)

    release = threading.Event()
    sent = []

    def _blocked_until_stopped(self, system="", messages=None, **kw):
        sent.append(len(messages or []))
        release.wait(20)
        return
        yield

    monkeypatch.setattr(A.OpenAIClient, "stream_message", _blocked_until_stopped)
    # The subject is the fresh start: the forced compaction tried before it
    # fails here, so the turn takes the fresh-start path (the compaction
    # path is pinned in test_a_compaction_before_a_fresh_start_saves_and_stops).
    _original_compact = AgentEngine._compact_history

    def _forced_compaction_fails(self, *, force=False):
        if force:
            raise RuntimeError("summariser unavailable")
        return _original_compact(self, force=force)
    monkeypatch.setattr(AgentEngine, "_compact_history",
                        _forced_compaction_fails)

    refs, w = _open(tmp_path, "k1")
    state = refs["state"]
    try:
        assert len(state["engine"].messages) == len(history)
        w[("Dropdown", "Model:")].value = "kit.glm-5.3"
        w[("Textarea", "")].value = "Hi"
        w[("Button", "Send")].click()
        deadline = time.monotonic() + 20
        while not sent and time.monotonic() < deadline:
            time.sleep(0.05)
        assert sent == [1], "the turn did not start fresh; the test lost its premise"
        w[("Button", "Stop")].click()
        release.set()
        deadline = time.monotonic() + 20
        while time.monotonic() < deadline:
            saved = ss.load_session(_SID) or {}
            if any("Stopped before any answer" in str(m.get("content"))
                   for m in saved.get("chat_messages") or []):
                break
            time.sleep(0.05)
        saved = ss.load_session(_SID) or {}
        assert len(saved.get("engine_messages") or []) == len(history), (
            "the saved session lost the conversation the chat still shows")
    finally:
        release.set()
        refs["shutdown"](save=False)
    ss.release_session_lock(_SID)

    refs2, _w2 = _open(tmp_path, "k2")
    try:
        assert len(refs2["state"]["engine"].messages) == len(history)
    finally:
        refs2["shutdown"](save=False)
