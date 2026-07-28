"""Crash-safety tests for agent session persistence.

Covers the three protections against unclean process death and
concurrent writers:

1. Atomic saves — every session_store write goes through temp file +
   ``os.replace``, so a simulated failure mid-write never leaves a torn
   session file (the old complete version survives).
2. Single-writer guard — a per-session ``<sid>.lock`` (pid + timestamp)
   makes a second live writer's save raise ``SessionLockedError``
   instead of silently clobbering; stale locks (dead pid or >1h old)
   are broken silently.
3. Mid-turn checkpoints — ``<sid>.turn.json`` written during long tool
   loops and cleared at turn end; a survivor is proof of an unclean
   death and injects a "[recovered]" note on the next restore.
"""

from __future__ import annotations

import json
import os
import subprocess
import time
from pathlib import Path

import pytest

from delfin.agent import session_store as ss


@pytest.fixture
def sessions_dir(monkeypatch, tmp_path):
    d = tmp_path / "agent_sessions"
    d.mkdir()
    monkeypatch.setattr(ss, "_SESSIONS_DIR", d)
    return d


def _dead_pid() -> int:
    """Spawn-and-reap a child so we hold a pid that is guaranteed dead."""
    p = subprocess.Popen(["true"])
    p.wait()
    return p.pid


# ---------------------------------------------------------------------------
# 1. Atomic saves
# ---------------------------------------------------------------------------


def _boom_replace(src, dst):
    raise OSError("simulated crash during replace")


def test_atomic_write_failure_preserves_old_file(sessions_dir):
    """A failure at the replace step must leave the previous complete
    file intact and no *.tmp litter behind."""
    target = sessions_dir / "s.json"
    ss._atomic_write_text(target, json.dumps({"v": 1}))

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(os, "replace", _boom_replace)
        with pytest.raises(OSError):
            ss._atomic_write_text(target, json.dumps({"v": 2}))

    # Old content survives, still valid JSON
    assert json.loads(target.read_text())["v"] == 1
    # No temp-file litter
    assert not list(sessions_dir.glob("*.tmp"))


def test_save_session_failure_leaves_previous_save_parseable(sessions_dir):
    ss.save_session("sid-a", mode="solo",
                    chat_messages=[{"role": "user", "content": "one"}])

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(os, "replace", _boom_replace)
        with pytest.raises(OSError):
            ss.save_session("sid-a", mode="solo",
                            chat_messages=[{"role": "user", "content": "two"}])

    data = ss.load_session("sid-a")
    assert data is not None
    assert data["chat_messages"][0]["content"] == "one"
    assert not list(sessions_dir.glob("*.tmp"))


# ---------------------------------------------------------------------------
# 2. Single-writer guard
# ---------------------------------------------------------------------------


def test_second_writer_with_live_pid_is_blocked(sessions_dir):
    """A fresh lock held by a DIFFERENT live pid (pid 1 is always alive)
    must make save_session raise SessionLockedError and skip the write."""
    (sessions_dir / "sid-b.lock").write_text(
        json.dumps({"pid": 1, "ts": time.time()}))
    with pytest.raises(ss.SessionLockedError) as exc:
        ss.save_session("sid-b", mode="solo",
                        chat_messages=[{"role": "user", "content": "x"}])
    assert exc.value.holder_pid == 1
    assert exc.value.session_id == "sid-b"
    # The save really was skipped
    assert not (sessions_dir / "sid-b.json").exists()


def test_dead_pid_lock_is_broken_silently(sessions_dir):
    (sessions_dir / "sid-c.lock").write_text(
        json.dumps({"pid": _dead_pid(), "ts": time.time()}))
    p = ss.save_session("sid-c", mode="solo",
                        chat_messages=[{"role": "user", "content": "x"}])
    assert p.exists()
    # Lock now belongs to us
    lock = json.loads((sessions_dir / "sid-c.lock").read_text())
    assert lock["pid"] == os.getpid()


def test_stale_lock_broken_by_age_even_with_live_pid(sessions_dir):
    (sessions_dir / "sid-d.lock").write_text(
        json.dumps({"pid": 1, "ts": time.time() - 7200}))
    p = ss.save_session("sid-d", mode="solo",
                        chat_messages=[{"role": "user", "content": "x"}])
    assert p.exists()


def test_own_pid_refreshes_lock_and_saves(sessions_dir):
    ss.save_session("sid-e", mode="solo",
                    chat_messages=[{"role": "user", "content": "x"}])
    t1 = json.loads((sessions_dir / "sid-e.lock").read_text())["ts"]
    time.sleep(0.02)
    ss.save_session("sid-e", mode="solo",
                    chat_messages=[{"role": "user", "content": "y"}])
    t2 = json.loads((sessions_dir / "sid-e.lock").read_text())["ts"]
    assert t2 > t1
    assert ss.load_session("sid-e")["chat_messages"][0]["content"] == "y"


def test_corrupt_lock_file_is_broken(sessions_dir):
    (sessions_dir / "sid-f.lock").write_text("{not json")
    assert ss.save_session(
        "sid-f", mode="solo",
        chat_messages=[{"role": "user", "content": "x"}],
    ).exists()


def test_release_session_lock_only_removes_own(sessions_dir):
    ss.acquire_session_lock("sid-g")
    ss.release_session_lock("sid-g")
    assert not (sessions_dir / "sid-g.lock").exists()
    # Someone else's lock stays put
    (sessions_dir / "sid-h.lock").write_text(
        json.dumps({"pid": 1, "ts": time.time()}))
    ss.release_session_lock("sid-h")
    assert (sessions_dir / "sid-h.lock").exists()


# ---------------------------------------------------------------------------
# 3. Mid-turn checkpoints
# ---------------------------------------------------------------------------


def test_checkpoint_write_load_clear(sessions_dir):
    p = ss.save_turn_checkpoint("sid-k", {
        "user_message": "long refactor",
        "partial_response": "did steps 1-3",
        "tool_calls": 12,
    })
    assert p is not None and p.name == "sid-k.turn.json"
    ck = ss.load_turn_checkpoint("sid-k")
    assert ck["user_message"] == "long refactor"
    assert ck["tool_calls"] == 12
    assert ck["ts"] > 0          # timestamp defaulted
    assert ss.clear_turn_checkpoint("sid-k") is True
    assert ss.load_turn_checkpoint("sid-k") is None
    assert ss.clear_turn_checkpoint("sid-k") is False


def test_checkpoint_empty_sid_is_noop(sessions_dir):
    assert ss.save_turn_checkpoint("", {"user_message": "x"}) is None
    assert ss.load_turn_checkpoint("") is None
    assert ss.clear_turn_checkpoint("") is False


def test_checkpoint_files_never_appear_as_sessions(sessions_dir):
    ss.save_session("sid-l", mode="solo",
                    chat_messages=[{"role": "user", "content": "x"}])
    ss.save_turn_checkpoint("sid-l", {"user_message": "x", "tool_calls": 1})
    ids = [s["session_id"] for s in ss.list_sessions()]
    assert ids == ["sid-l"]


def test_fresh_checkpoint_yields_recovery_note_once(sessions_dir):
    ss.save_session("sid-m", mode="solo",
                    chat_messages=[{"role": "user", "content": "goal"}])
    time.sleep(0.02)   # checkpoint strictly newer than the save
    ss.save_turn_checkpoint("sid-m", {
        "user_message": "harden the parser",
        "partial_response": "patched lexer, tests pending",
        "tool_calls": 47,
    })
    note = ss.consume_crash_recovery_note("sid-m")
    assert note is not None
    assert note.startswith("[recovered]")
    assert "harden the parser" in note
    assert "47 tool calls" in note
    assert "patched lexer, tests pending" in note
    assert "verify state before continuing" in note
    # One-shot: checkpoint consumed
    assert ss.load_turn_checkpoint("sid-m") is None
    assert ss.consume_crash_recovery_note("sid-m") is None


def test_stale_checkpoint_older_than_save_is_dropped(sessions_dir):
    ss.save_turn_checkpoint("sid-n", {
        "user_message": "old work", "tool_calls": 3,
        "ts": time.time() - 100,
    })
    ss.save_session("sid-n", mode="solo",
                    chat_messages=[{"role": "user", "content": "x"}])
    assert ss.consume_crash_recovery_note("sid-n") is None
    # Stale checkpoint was still cleaned up
    assert ss.load_turn_checkpoint("sid-n") is None


def test_checkpoint_without_session_file_still_recovers(sessions_dir):
    """Process died before the FIRST turn-boundary save: no session file
    exists, but the checkpoint alone must still produce the note."""
    ss.save_turn_checkpoint("sid-o", {
        "user_message": "first task", "partial_response": "step 1 done",
        "tool_calls": 5,
    })
    note = ss.consume_crash_recovery_note("sid-o")
    assert note is not None and "first task" in note


def test_delete_session_removes_lock_and_checkpoint(sessions_dir):
    ss.save_session("sid-p", mode="solo",
                    chat_messages=[{"role": "user", "content": "x"}])
    ss.save_turn_checkpoint("sid-p", {"user_message": "x", "tool_calls": 1})
    assert ss.delete_session("sid-p") is True
    assert not (sessions_dir / "sid-p.lock").exists()
    assert not (sessions_dir / "sid-p.turn.json").exists()


# ---------------------------------------------------------------------------
# Engine restore path — recovery note injection
# ---------------------------------------------------------------------------


def _bare_engine(mode: str = "solo"):
    from delfin.agent.engine import AgentEngine
    eng = AgentEngine.__new__(AgentEngine)
    eng.mode = mode
    eng.route = ["solo_agent"]
    eng.current_role_index = 0
    eng.role_outputs = {}
    eng.compaction_summaries = {}
    eng.messages = []
    eng.token_usage = {"input": 0, "output": 0}
    eng.cost_usd = 0.0
    eng.session_id = ""
    eng._project_dir = ""
    eng._last_input_tokens = 0
    return eng


def test_restore_state_injects_recovery_note_pair(sessions_dir):
    ss.save_session("sid-q", mode="solo",
                    chat_messages=[{"role": "user", "content": "goal"}])
    time.sleep(0.02)
    ss.save_turn_checkpoint("sid-q", {
        "user_message": "migrate the config loader",
        "partial_response": "moved 3 of 7 call sites",
        "tool_calls": 23,
    })
    eng = _bare_engine()
    eng.restore_state({
        "mode": "solo",
        "engine_messages": [
            {"role": "user", "content": "goal"},
            {"role": "assistant", "content": "done with step 0"},
        ],
        "session_id": "sid-q",
    })
    # Note pair appended: user note + assistant ack (house pattern, so
    # message alternation survives _sanitize_messages on the next turn)
    assert len(eng.messages) == 4
    assert eng.messages[2]["role"] == "user"
    assert eng.messages[2]["content"].startswith("[recovered]")
    assert "migrate the config loader" in eng.messages[2]["content"]
    assert "23 tool calls" in eng.messages[2]["content"]
    assert "moved 3 of 7 call sites" in eng.messages[2]["content"]
    assert eng.messages[3]["role"] == "assistant"
    # Checkpoint consumed — restoring again injects nothing
    eng2 = _bare_engine()
    eng2.restore_state({
        "mode": "solo",
        "engine_messages": [{"role": "user", "content": "goal"}],
        "session_id": "sid-q",
    })
    assert len(eng2.messages) == 1


def test_restore_state_without_checkpoint_adds_nothing(sessions_dir):
    eng = _bare_engine()
    eng.restore_state({
        "mode": "solo",
        "engine_messages": [{"role": "user", "content": "hi"}],
        "session_id": "sid-r",
    })
    assert len(eng.messages) == 1


# ---------------------------------------------------------------------------
# Engine stream path — checkpoint written mid-loop, cleared at turn end
# ---------------------------------------------------------------------------


@pytest.fixture
def agent_tree(tmp_path):
    """Minimal agent pack tree (same shape as test_agent_engine.py)."""
    import textwrap
    agent_dir = tmp_path / "pack"
    shared = agent_dir / "shared"
    shared.mkdir(parents=True)
    agents = agent_dir / "agents"
    agents.mkdir()
    (shared / "delfin_context.md").write_text("# Context")
    (shared / "work_cycle_rules.md").write_text("# Rules")
    (shared / "goal_decomposition_rules.md").write_text("# Goal Decomposition")
    (shared / "universal_input_template.md").write_text("")
    (shared / "minimal_final_verdict.md").write_text("")
    (agents / "session_manager.md").write_text("# Session Manager")
    (agents / "builder_agent.md").write_text("# Builder Agent")
    (agents / "test_agent.md").write_text("# Test Agent")
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "quick.md").write_text("# quick mode")
    (lite_dir / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: quick
            file: modes/quick.md
            route:
              - session_manager
              - builder_agent
              - test_agent
    """))
    return tmp_path


def test_stream_writes_throttled_checkpoint_and_clears_at_turn_end(
    sessions_dir, agent_tree, monkeypatch,
):
    """11 tool rounds in one turn: the 10th tool_result crosses the
    throttle so exactly one mid-turn checkpoint is written; the normal
    turn end then clears it (only a SIGKILL would leave it behind)."""
    from unittest.mock import MagicMock, patch
    from delfin.agent.api_client import StreamEvent

    def fake_stream(system, messages, max_tokens=4096, session_id="",
                    thinking_budget=0):
        yield StreamEvent(type="session_init", text="crash-sid-1")
        yield StreamEvent(type="message_start", input_tokens=100)
        yield StreamEvent(type="text_delta", text="working... ")
        for i in range(11):
            yield StreamEvent(type="tool_use", tool_name="Bash",
                              tool_input='{"command": "ls"}')
            yield StreamEvent(type="tool_result", tool_name="Bash",
                              tool_output=f"round {i}")
        yield StreamEvent(type="text_delta", text="done")
        yield StreamEvent(type="message_delta", output_tokens=10, cost_usd=0.01)

    client = MagicMock()
    client.stream_message = MagicMock(side_effect=fake_stream)

    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=agent_tree, backend="cli",
                             mode="quick", pack_dir=agent_tree)

    saved: list[tuple[str, dict]] = []
    real_save = ss.save_turn_checkpoint

    def _recording_save(sid, payload):
        saved.append((sid, dict(payload)))
        return real_save(sid, payload)

    monkeypatch.setattr(ss, "save_turn_checkpoint", _recording_save)

    engine.stream_response("do a long refactor")

    assert len(saved) == 1                      # 11 results → one throttled write
    sid, payload = saved[0]
    assert sid == "crash-sid-1"
    assert payload["user_message"] == "do a long refactor"
    assert payload["tool_calls"] == 10          # count at the 10th tool_result
    assert "working" in payload["partial_response"]
    # Normal turn end cleared the checkpoint
    assert ss.load_turn_checkpoint("crash-sid-1") is None
    assert not (sessions_dir / "crash-sid-1.turn.json").exists()


# ---------------------------------------------------------------------------
# Restore completeness — project_dir + last_input_tokens round-trip
# ---------------------------------------------------------------------------


def test_export_restore_roundtrips_project_dir_and_token_floor(sessions_dir):
    eng = _bare_engine()
    eng.session_id = "sid-s"
    eng._project_dir = "/work/proj/calc_07"
    eng._last_input_tokens = 48_213
    state = eng.export_state()
    assert state["project_dir"] == "/work/proj/calc_07"
    assert state["last_input_tokens"] == 48_213

    eng2 = _bare_engine()
    eng2.restore_state(state)
    assert eng2._project_dir == "/work/proj/calc_07"
    assert eng2._last_input_tokens == 48_213


def test_restore_defaults_new_fields_for_legacy_sessions(sessions_dir):
    """Sessions saved before these fields existed must restore clean."""
    eng = _bare_engine()
    eng._project_dir = "leftover"      # must be overwritten, not kept
    eng._last_input_tokens = 999
    eng.restore_state({
        "mode": "solo",
        "engine_messages": [],
        "session_id": "sid-t",
    })
    assert eng._project_dir == ""
    assert eng._last_input_tokens == 0
