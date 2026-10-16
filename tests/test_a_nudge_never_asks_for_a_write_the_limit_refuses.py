"""The memory nudge and the memory-write limit agree; a peer session is
found by the name the sender used or by the key it had before a restore.

Field data (33 sessions, 2026-09-16..10-14): 41 nudges, 17 answered by an
immediate remember; of 60 remember/forget calls 9 were refused by the
per-turn limit ("this turn has already recorded 2 fact(s)") and 4 by plan
mode. A nudge that asks for a write the next check refuses costs the round
it provoked. session_message failed 10 times in 140 calls: by a title the
session list does not hold verbatim ("Session A", "Session A -- denials
report" for a session titled "Hallo Session A") and by a key the peer had
before it was re-opened.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import memory_nudge as mn


# ---------------------------------------------------------------------------
# maybe_nudge: silent when the write it asks for would be refused
# ---------------------------------------------------------------------------

def _worked(state):
    return mn.maybe_nudge(state, tool_calls=mn.NUDGE_TOOL_CALLS,
                          chars=mn.NUDGE_CHARS, settings={})


def test_no_nudge_when_no_memory_write_is_left_this_turn():
    st = mn.NudgeState()
    assert mn.maybe_nudge(st, tool_calls=mn.NUDGE_TOOL_CALLS,
                          chars=mn.NUDGE_CHARS, settings={},
                          writes_left=0) == ""
    # The work is not thrown away: it nudges once a write is possible.
    assert mn.maybe_nudge(st, settings={}, writes_left=1)


def test_no_nudge_in_plan_mode():
    st = mn.NudgeState()
    assert mn.maybe_nudge(st, tool_calls=mn.NUDGE_TOOL_CALLS,
                          chars=mn.NUDGE_CHARS, settings={},
                          plan_mode=True) == ""


def test_the_default_call_is_unchanged():
    assert _worked(mn.NudgeState())


# ---------------------------------------------------------------------------
# Driven through the real client loop (stub backend, real executor)
# ---------------------------------------------------------------------------

from tests.test_memory_context_nudge_integration import (  # noqa: E402
    _count_new_nudges, _read_args, _round_files, _run_turn, _stub_client,
)


@pytest.fixture
def workspace(tmp_path):
    for i in range(24):
        (tmp_path / f"big{i}.txt").write_text(f"finding {i} " * 600)
    return tmp_path


def test_driven_a_turn_at_its_write_limit_is_not_nudged(workspace):
    from delfin.agent.model_profiles import get_profile

    cap = int(get_profile("kit.glm-5.3").max_memory_writes_per_turn or 0)
    assert cap > 0
    remember = [("remember", json.dumps({"fact": f"note {i}"}))
                for i in range(cap + 1)]          # one over: refused
    # 18 reads: the allowance earns one write per 10 other calls, so
    # after the over-limit round no write is left until the 20th read.
    script = ([remember]
              + [[("read_file", _read_args(f))]
                 for f in _round_files(workspace, 18)]
              + ["done"])
    client, state = _stub_client(script, workspace_root=workspace)
    _run_turn(client)
    assert _count_new_nudges(state["requests"]) == 0


def test_driven_a_plan_mode_turn_is_not_nudged(workspace):
    from delfin.agent import api_client as ac

    script = ([[("read_file", _read_args(f))]
               for f in _round_files(workspace, 20)] + ["done"])
    client, state = _stub_client(script, workspace_root=workspace)
    client._permissions = ac.KitToolPermissions(
        workspace=workspace, mode="plan")
    _run_turn(client)
    assert _count_new_nudges(state["requests"]) == 0


# ---------------------------------------------------------------------------
# session_message: the name the sender used, and a key from before a restore
# ---------------------------------------------------------------------------

class _Perms:
    def __init__(self, key):
        self.presence_key = key


@pytest.fixture
def executor():
    from delfin.agent.api_client import _DocToolExecutor

    return object.__new__(_DocToolExecutor)


@pytest.fixture
def peers(monkeypatch):
    import delfin.agent.session_messages as msgs
    import delfin.agent.session_presence as pres

    rows = [
        {"key": "2a15fb71", "title": "Hallo Session A", "session_id": "s-a"},
        {"key": "f2100819", "title": "Hallo Session B", "session_id": "s-b"},
    ]
    sent = []
    monkeypatch.setattr(pres, "open_sessions",
                        lambda **kw: [r for r in rows
                                      if r["key"] != kw.get("exclude_key")])
    monkeypatch.setattr(msgs, "deliverable", lambda key: False)
    monkeypatch.setattr(msgs, "send",
                        lambda to, text, **kw: sent.append((to, text)))
    return sent


def _send(executor, to, me="3eb6c29f"):
    return json.loads(executor._execute_session_message(
        {"to": to, "message": "hi"}, _Perms(me)))


@pytest.mark.parametrize("to, key", [
    ("Session A", "2a15fb71"),
    ("session b", "f2100819"),
    ("Session A — denials report", "2a15fb71"),
    ("Session B - cost report", "f2100819"),
])
def test_a_name_that_identifies_one_session_is_delivered(
        executor, peers, to, key):
    out = _send(executor, to)
    assert out.get("status") == "sent", out
    assert peers == [(key, "hi")]


@pytest.mark.parametrize("to", ["Session", "Hallo", "se"])
def test_a_name_shared_by_several_sessions_is_not_guessed(
        executor, peers, to):
    out = _send(executor, to)
    assert "error" in out and peers == []
    assert {r["key"] for r in out["sessions"]} == {"2a15fb71", "f2100819"}


def test_a_key_from_before_a_restore_reaches_the_session(
        tmp_path, monkeypatch, executor):
    """Real presence records: a session first announced under a fresh key
    (no conversation id yet), closed, and re-opened under its stable key."""
    import delfin.agent.session_messages as msgs
    import delfin.agent.session_presence as pres

    monkeypatch.setattr(pres, "_DIR", tmp_path / "session_presence")
    monkeypatch.setattr(pres, "_last_written", {})
    monkeypatch.setattr(pres, "repository_of",
                        lambda ws: {"root": "", "common_dir": "",
                                    "branch": ""})
    sent = []
    monkeypatch.setattr(msgs, "send",
                        lambda to, text, **kw: sent.append((to, text)))
    sid = "7cdf74d4fdda4a539b1fa3e8e60a7067"
    pres.announce("e1e2e3e4", session_id=sid, title="Hallo Session A")
    pres.withdraw("e1e2e3e4")
    pres.announce(sid[:8], session_id=sid, title="Hallo Session A")
    out = _send(executor, "e1e2e3e4")
    assert out.get("status") == "sent", out
    assert sent == [(sid[:8], "hi")]
