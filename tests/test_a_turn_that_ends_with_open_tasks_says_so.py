"""A turn that stops with work still on the list must say so.

THE INCIDENT. The auto-continue guard evaluated the open-task state one
line above the exit and threw the answer away: when the "did tools since
the last continue" flag was False, or the 12-continue cap was spent,
control fell through to a bare terminal event — a token count and a
cost. The state that PROVES the work is unfinished had just been
computed. The only other consumer was a prompt block that reaches the
model on the NEXT turn, if the user sends one, so from the user's side a
turn that stopped one task in and a turn that finished everything looked
identical.

``_has_pending_tasks`` also failed CLOSED: any exception — a corrupt
task store, an unreadable file — returned False, so a broken ledger and
finished work were the same answer. The tri-state exists for that: the
third value is "I could not tell".

The notice is emitted as stream text, which is also what covers the
headless path: it prints the stream and then tokens and cost, and never
reads the task store itself.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import mcp_client as M
from delfin.agent.agent_tasks import (
    format_open_tasks_notice, get_store, open_task_summary,
)


# --- Fake OpenAI streaming primitives (mirrors test_loop_hardening) ---------

class _FnDelta:
    def __init__(self, name=None, arguments=None):
        self.name = name
        self.arguments = arguments


class _TCDelta:
    def __init__(self, index, id=None, name=None, arguments=None):
        self.index = index
        self.id = id
        self.function = _FnDelta(name, arguments)


class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content = content
        self.tool_calls = tool_calls


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta = delta
        self.finish_reason = finish


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices = choices
        self.usage = usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


class _EmptyReg:
    def discover_all(self): return []
    def discover_resources(self): return []
    def discover_prompts(self): return []


def _final(text="I have created the plan."):
    return _Stream([
        _Chunk([_Choice(_Delta(content=text), finish="stop")],
               usage=_Usage()),
    ])


@pytest.fixture
def workspace(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    from delfin.agent import agent_tasks
    agent_tasks._STORES.clear()
    return ws


@pytest.fixture
def client(monkeypatch, workspace):
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    from delfin.agent.api_client import KitToolPermissions, create_client
    cl = create_client(backend="api", provider="ollama",
                       model="qwen2.5-coder:7b", cwd=str(workspace))
    cl.set_permissions(KitToolPermissions(
        workspace=workspace, mode="default", task_session_id="sess-1"))
    return cl


def _run(client, streams):
    calls = {"n": 0}

    def _create(**kwargs):
        s = streams[min(calls["n"], len(streams) - 1)]
        calls["n"] += 1
        return s

    client.client.chat.completions.create = _create
    return list(client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=100))


def _text(events) -> str:
    """What reached the user, model answer and harness notice alike. The
    open-tasks sentence is the framework speaking, so it rides the notice
    type and is deliberately not part of the assistant's message."""
    return "".join(getattr(e, "text", "") or "" for e in events
                   if getattr(e, "type", "") in ("text_delta", "notice"))


def _stop(events) -> str:
    md = [e for e in events if getattr(e, "type", "") == "message_delta"]
    return md[-1].stop_reason if md else ""


# ---------------------------------------------------------------------------
# The silent exit
# ---------------------------------------------------------------------------

def test_a_turn_ending_with_an_open_task_names_it(client, workspace):
    store = get_store(workspace)
    store.create("Wire the parser into ops_server", session_id="sess-1")
    events = _run(client, [_final()])
    text = _text(events)
    assert "open task" in text
    assert "Wire the parser into ops_server" in text


def test_the_terminal_reason_is_distinct(client, workspace):
    get_store(workspace).create("Add the regression test", session_id="sess-1")
    assert _stop(_run(client, [_final()])) == "end_turn_open_tasks"


def test_a_turn_that_really_finished_says_nothing_extra(client, workspace):
    """The notice must not fire on a clean turn, or it is noise within a
    day and unread within two."""
    store = get_store(workspace)
    t = store.create("Ship it", session_id="sess-1")
    store.update(t["id"], status="in_progress")
    store.update(t["id"], status="completed")
    events = _run(client, [_final()])
    assert "open task" not in _text(events)
    assert _stop(events) == "end_turn"


def test_another_sessions_tasks_do_not_raise_the_notice(client, workspace):
    get_store(workspace).create("someone else's work", session_id="other")
    events = _run(client, [_final()])
    assert "open task" not in _text(events)


def test_a_blocked_task_is_reported_with_what_it_waits_on(client, workspace):
    store = get_store(workspace)
    t = store.create("Deploy to the cluster", session_id="sess-1")
    store.update(t["id"], status="blocked",
                 blocked_reason="the KIT API key is missing")
    text = _text(_run(client, [_final()]))
    assert "Deploy to the cluster" in text
    assert "KIT API key" in text


# ---------------------------------------------------------------------------
# A ledger nobody can read is not "nothing to do"
# ---------------------------------------------------------------------------

def test_an_unreadable_task_store_is_unknown_not_empty(workspace):
    path = workspace / ".delfin" / "session_tasks.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("{not json at all", encoding="utf-8")
    summary = open_task_summary(workspace, "sess-1")
    assert summary["state"] == "unknown"
    assert summary["error"]


def test_the_notice_says_the_ledger_could_not_be_read(client, workspace):
    path = workspace / ".delfin" / "session_tasks.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("{not json at all", encoding="utf-8")
    events = _run(client, [_final()])
    assert "could not be read" in _text(events)
    assert _stop(events) == "end_turn_tasks_unknown"


def test_auto_continue_still_declines_on_an_unknown_ledger(client, workspace):
    """Unknown must not become a reason to keep driving the model — it is a
    reason to tell the user."""
    path = workspace / ".delfin" / "session_tasks.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("{not json at all", encoding="utf-8")
    assert client._has_pending_tasks() is False
    assert client._open_task_state()["state"] == "unknown"


def test_a_blocked_task_alone_does_not_drive_auto_continue(client, workspace):
    """Continuing into work that waits on a user answer is a loop, not
    progress — but it is still named at turn end."""
    store = get_store(workspace)
    t = store.create("Get the credential", session_id="sess-1")
    store.update(t["id"], status="blocked", blocked_reason="user answer")
    assert client._has_pending_tasks() is False
    assert client._open_task_state()["state"] == "open"
    assert "Get the credential" in _text(_run(client, [_final()]))


# ---------------------------------------------------------------------------
# The renderer
# ---------------------------------------------------------------------------

def test_the_notice_is_empty_only_for_a_readable_empty_list():
    assert format_open_tasks_notice({"state": "none"}) == ""
    assert format_open_tasks_notice({"state": "unknown", "error": "boom"})
    assert format_open_tasks_notice({}) != ""      # no state == not "none"


def test_the_notice_counts_everything_and_lists_a_bounded_number(workspace):
    store = get_store(workspace)
    for i in range(12):
        store.create(f"step {i}", session_id="s")
    notice = format_open_tasks_notice(open_task_summary(workspace, "s"))
    assert "12 open task(s)" in notice
    assert "more" in notice                       # the tail is summarised
    assert len(notice.splitlines()) <= 10
