"""Both completion drains lived inside the tool-calling branch only.

``_drain_background_events`` (finished bash/cluster jobs) and
``_drain_turn_steering`` (the engine's drain-backed blocks: finished
delegates, a late answer to a parked question) were called after a round
that executed tools, and nowhere else. The turn-end path drained the
USER steer inbox and then broke out.

So a turn that ended without tool calls surfaced nothing -- and that is
precisely the "I started the job, I'll wait" ending, after which the
agent never takes another turn: there is no timer, no completion
callback and no autonomous turn anywhere in the codebase. Nothing would
ever wake it, and the next drain that did run would consume the event
without anyone reading it.

The fix drains in the no-tool-calls path too, before the loop breaks,
and keeps the turn going so the model can act on what it was handed --
the drains are exactly-once, so an event delivered into a turn that ends
immediately afterwards is an event destroyed.
"""

from __future__ import annotations

import pytest

from delfin.agent import mcp_client as M


# --- Fake OpenAI streaming primitives (mirrors test_loop_hardening) ---------

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


def _final(text="I started the job; I'll wait for it."):
    return _Stream([
        _Chunk([_Choice(_Delta(content=text), finish="stop")], usage=_Usage()),
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


def _run(client, rounds=6):
    """Drive the loop; every round returns a final answer, no tool calls.

    The messages of every round are kept: what the MODEL was handed is
    the half of this that matters, and it is invisible in the event
    stream.
    """
    calls = {"n": 0, "messages": []}

    def _create(**kwargs):
        calls["n"] += 1
        calls["messages"].append(list(kwargs.get("messages") or []))
        return _final()

    client.client.chat.completions.create = _create
    events = list(client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=100))
    return events, calls


def _text(events) -> str:
    """What the model itself said. A harness sentence must never appear
    here — that is the separation #146 installed, and asserting a
    completion notice in this channel would ask for it back."""
    return "".join(getattr(e, "text", "") or "" for e in events
                   if getattr(e, "type", "") == "text_delta")


def _shown(events) -> str:
    """Everything the user is shown, on whichever channel."""
    return "".join(getattr(e, "text", "") or "" for e in events
                   if getattr(e, "type", "") in ("text_delta", "notice"))


def _handed_to_model(calls) -> str:
    return "\n".join(str(m.get("content", ""))
                     for round_msgs in calls["messages"]
                     for m in round_msgs)


def _one_job_event(monkeypatch, jid="5c294b67"):
    """A single finished job, drained exactly once (as the real store is)."""
    seen = {"drained": False}

    def _drain(ws):
        if seen["drained"]:
            return []
        seen["drained"] = True
        return [{"job_id": jid, "exit_code": 0, "description": "orca opt"}]

    monkeypatch.setattr("delfin.agent.bash_jobs.drain_finished_events", _drain)
    return seen


# ---------------------------------------------------------------------------
# The event reaches a turn that called no tools
# ---------------------------------------------------------------------------

def test_a_job_that_finished_reaches_a_turn_without_tool_calls(
        client, monkeypatch):
    _one_job_event(monkeypatch)
    events, calls = _run(client)
    assert "5c294b67" in _shown(events), (
        "the job finished and the turn ended without ever mentioning it")
    assert "5c294b67" in _handed_to_model(calls), (
        "the user saw it and the model did not — the model is the one that "
        "has to act on it")
    assert calls["n"] >= 2, "the turn ended instead of letting the model react"


def test_the_completion_is_not_the_models_own_answer(client, monkeypatch):
    """It reaches the user, and not through the channel a benchmark or a
    verifier reads as what the model said."""
    _one_job_event(monkeypatch)
    events, _calls = _run(client)
    assert "5c294b67" not in _text(events)
    assert any(getattr(e, "type", "") == "notice"
               and "5c294b67" in (getattr(e, "text", "") or "")
               for e in events)


def test_the_completion_is_not_repeated(client, monkeypatch):
    """Exactly-once: the drain is asked again and answers nothing."""
    _one_job_event(monkeypatch)
    events, _calls = _run(client)
    assert _shown(events).count("5c294b67") == 1


def test_a_turn_with_nothing_pending_ends_at_once(client, monkeypatch):
    monkeypatch.setattr("delfin.agent.bash_jobs.drain_finished_events",
                        lambda ws: [])
    _events, calls = _run(client)
    assert calls["n"] == 1, "an empty drain still cost the turn an extra round"


def test_a_delegate_report_reaches_the_same_path(client, monkeypatch):
    """The engine's drain-backed blocks ride the same turn-end drain."""
    monkeypatch.setattr("delfin.agent.bash_jobs.drain_finished_events",
                        lambda ws: [])
    # Hand the block over exactly once, the way the engine's drain does.
    handed = {"n": 0}

    def _provider():
        if handed["n"]:
            return []
        handed["n"] += 1
        return ["- sub-agent explore [bg7] found the parser in editblock.py"]

    client.steering_provider = _provider
    _events, calls = _run(client)
    assert handed["n"] == 1
    assert calls["n"] >= 2, "a finished delegate did not keep the turn alive"


def test_the_extra_rounds_are_capped(client, monkeypatch):
    """A drain that keeps producing must not be able to run forever."""
    monkeypatch.setattr("delfin.agent.bash_jobs.drain_finished_events",
                        lambda ws: [{"job_id": "endless", "exit_code": 0}])
    _events, calls = _run(client)
    assert calls["n"] <= 5


def test_a_broken_drain_does_not_end_the_turn(client, monkeypatch):
    monkeypatch.setattr(
        "delfin.agent.bash_jobs.drain_finished_events",
        lambda ws: (_ for _ in ()).throw(OSError("registry down")))
    events, calls = _run(client)
    assert calls["n"] == 1
    assert any(getattr(e, "type", "") == "message_delta" for e in events)


# ---------------------------------------------------------------------------
# Both drains are on the exit path
# ---------------------------------------------------------------------------

def test_both_drains_run_before_the_loop_breaks(client, monkeypatch):
    """A finished job and a run note are two different inboxes, and the
    turn-end path has to empty both.

    This was a source-text assertion: it read api_client.py and required
    one drain call to appear after another. It broke the moment the code
    around it moved, while the behaviour was unchanged — the failure said
    the mechanism was gone when only its line number was.
    """
    _one_job_event(monkeypatch)
    client.push_run_note("⚙️ permission mode → acceptEdits")
    _events, calls = _run(client)
    handed = _handed_to_model(calls)
    assert "5c294b67" in handed
    assert "acceptEdits" in handed


def test_the_steer_inbox_is_emptied_on_the_same_path(client, monkeypatch):
    monkeypatch.setattr("delfin.agent.bash_jobs.drain_finished_events",
                        lambda ws: [])
    client.push_steer("stop after this one")
    _events, calls = _run(client)
    assert "stop after this one" in _handed_to_model(calls)
