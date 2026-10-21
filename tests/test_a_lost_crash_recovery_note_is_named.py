"""A silent except in restore_state hides a lost crash-recovery note.

`AgentEngine.restore_state` reads the surviving mid-turn checkpoint
(session_store.consume_crash_recovery_note) and appends the "[recovered]"
note pair to the restored conversation. That whole block was wrapped in
``except Exception: pass`` — if the import or the read blew up, the note
was silently dropped and the resumed session continued believing the
previous turn had finished cleanly. Nothing in the restored report or the
stream said so.

These tests pin the new contract: when the recovery-note read fails inside
restore_state, the loss is NAMED in the RestoreReport (which the CLI prints
for an incomplete restore) instead of being swallowed, and the checkpoint
file is left in place — a failed read must not delete the crash evidence
it failed to read. A successful read stays exactly as silent as before.
"""

from __future__ import annotations

import json
import time
from unittest.mock import patch

import pytest

from delfin.agent import session_store as ss


@pytest.fixture
def sessions_dir(monkeypatch, tmp_path):
    d = tmp_path / "agent_sessions"
    d.mkdir()
    monkeypatch.setattr(ss, "_SESSIONS_DIR", d)
    return d


def _bare_engine():
    from delfin.agent.engine import AgentEngine
    eng = AgentEngine.__new__(AgentEngine)
    eng.mode = "solo"
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


def _saved_checkpoint(sid: str, sessions_dir) -> None:
    """A session plus a checkpoint newer than the save: crash evidence."""
    ss.save_session(sid, mode="solo",
                    chat_messages=[{"role": "user", "content": "goal"}])
    time.sleep(0.02)
    ss.save_turn_checkpoint(sid, {
        "user_message": "migrate the config loader",
        "partial_response": "moved 3 of 7 call sites",
        "tool_calls": 23,
    })
    assert ss.load_turn_checkpoint(sid) is not None


def test_import_failure_of_the_recovery_note_read_is_said(
        sessions_dir, capsys):
    """If the store module cannot even be imported, restore_state raises
    (ImportError from the early restore read — this is the upstream half
    of the read that lives before the note block). The test pins that the
    early import is the one that decides, and that a failing checkpoint
    READ inside the note block is covered by the second test below."""
    _saved_checkpoint("sid-imp", sessions_dir)
    eng = _bare_engine()

    def _boom_import(name, *args, **kwargs):
        if name == "delfin.agent.session_store" or name == ".session_store":
            raise RuntimeError("simulated import failure")
        real = type(__import__)  # the real builtins.__import__ builtin
        with patch("builtins.__import__", real):
            return real(name, *args, **kwargs)

    notices: list[str] = []
    with patch("builtins.__import__", side_effect=_boom_import):
        with pytest.raises(ImportError):
            eng.restore_state({
                "mode": "solo",
                "engine_messages": [{"role": "user", "content": "goal"}],
                "session_id": "sid-imp",
            })
    # The checkpoint evidence is NOT deleted by the failed restore read:
    # a later, fixed restore can still recover from it.
    assert ss.load_turn_checkpoint("sid-imp") is not None


def test_checkpoint_read_failure_is_said_in_the_report(
        sessions_dir, monkeypatch):
    """A crash checkpoint exists but reading it blows up: the note is
    lost, and the loss must reach the RestoreReport (which the CLI prints
    on an incomplete restore)."""
    _saved_checkpoint("sid-read", sessions_dir)
    eng = _bare_engine()

    def _boom_consume(session_id):
        raise OSError("simulated checkpoint read failure")

    monkeypatch.setattr(
        "delfin.agent.session_store.consume_crash_recovery_note",
        _boom_consume)
    eng.restore_state({
        "mode": "solo",
        "engine_messages": [{"role": "user", "content": "goal"}],
        "session_id": "sid-read",
    })
    report = eng.last_restore_report
    assert report is not None
    assert any("recovery note" in str(d) for d in report.failed), \
        "a lost crash-recovery note must be named in the restore report"
    # The checkpoint itself is untouched — a failed read must not delete
    # the crash evidence it failed to read.
    assert ss.load_turn_checkpoint("sid-read") is not None


def test_working_read_stays_silent(sessions_dir):
    """The normal path gains no new machinery speech: a successful
    restore with the note pair appended reports nothing extra."""
    _saved_checkpoint("sid-ok", sessions_dir)
    eng = _bare_engine()
    eng.restore_state({
        "mode": "solo",
        "engine_messages": [{"role": "user", "content": "goal"}],
        "session_id": "sid-ok",
    })
    report = eng.last_restore_report
    assert report is not None
    assert not [d for d in report.failed if "recovery note" in str(d)]
    # And the note pair is still injected, as before: the goal, the
    # recovery note, the ack — note pair appended after the restored
    # engine_messages (3 restored rows already carry the goal).
    assert any(m["content"].startswith("[recovered]") for m in eng.messages)
    assert eng.messages[-1]["content"].startswith("Understood. I will verify")
