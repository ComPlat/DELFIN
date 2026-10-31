"""The step journal stays bounded, knows its turns, and leaves with its session.

Input: delfin/agent/step_journal.py, the DR1 hook in
api_client._DocToolExecutor._execute_pinned, the engine's turn start
(AgentEngine.stream_response) and session_store.delete_session.

Before:
- every tool call appended two lines to ``<session>.steps.jsonl`` and
  nothing ever removed one, so a long unattended session grew the file
  without a bound (readers look at the last 256 KB only);
- no caller passed ``turn``, so every step of every turn was turn 1. A
  crash note then listed an open step left by an EARLIER crash, already
  reported and dealt with, as if it belonged to the turn that just died;
- delete_session removed the lock, the checkpoint and the elided store
  but left the journal behind.
"""
from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest

from delfin.agent import api_client, session_store, step_journal


@pytest.fixture()
def sessions(tmp_path, monkeypatch):
    d = tmp_path / "sess"
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", d)
    return d


# ---------------------------------------------------------------------------
# Bound
# ---------------------------------------------------------------------------

def test_the_journal_on_disk_stays_bounded(sessions, monkeypatch):
    monkeypatch.setattr(step_journal, "_ROTATE_BYTES", 64 * 1024)
    big = {"command": "x" * 150}
    for _ in range(2000):
        rec = step_journal.auto_start("s1", tool="bash", arguments=big)
        step_journal.auto_finish("s1", rec, ok=True, result="done")
    files = sorted(sessions.glob("s1.steps.jsonl*"))
    total = sum(f.stat().st_size for f in files)
    assert total <= 2 * 64 * 1024 + 4096, (files, total)
    assert all(f.stat().st_mode & 0o077 == 0 for f in files)


def test_a_rotation_keeps_what_a_read_sees(sessions, monkeypatch):
    """Right after a rotation the live file holds two lines; the open step
    written just before it, now in the rotated file, must still be found."""
    for _ in range(10):
        rec = step_journal.auto_start("s1", tool="read_file",
                                      arguments={"path": "a"})
        step_journal.auto_finish("s1", rec, ok=True, result="")
    step_journal.auto_start("s1", tool="bash",
                            arguments={"command": "sbatch job.sh"})
    live = step_journal.journal_path("s1")
    # The next append finds the file at its bound and rotates first.
    monkeypatch.setattr(step_journal, "_ROTATE_BYTES", live.stat().st_size)
    rec = step_journal.auto_start("s1", tool="read_file",
                                  arguments={"path": "c"})
    step_journal.auto_finish("s1", rec, ok=True, result="")
    assert live.with_name(live.name + ".1").is_file()
    assert len(live.read_text().splitlines()) == 2
    open_out, _ = step_journal.open_steps("s1")
    assert [e["tool"] for e in open_out] == ["bash"]
    assert rec["step"] == 12


def test_a_read_that_ends_at_the_file_boundary_keeps_the_first_line(
        sessions, monkeypatch):
    """The live file exactly as long as the read window: its first line
    is whole, and is not dropped as a cut one."""
    for i in range(3):
        rec = step_journal.auto_start("s1", tool="bash",
                                      arguments={"command": f"c{i}"})
        step_journal.auto_finish("s1", rec, ok=True)
    live = step_journal.journal_path("s1")
    live.with_name(live.name + ".1").write_text('{"kind": "start"}\n')
    monkeypatch.setattr(step_journal, "_TAIL_BYTES", live.stat().st_size)
    steps, _ = step_journal.read_steps("s1")
    assert [(e["kind"], e["step"]) for e in steps][:2] == [
        ("start", 1), ("finish", 1)]


def test_the_step_numbers_continue_across_a_rotation(sessions, monkeypatch):
    monkeypatch.setattr(step_journal, "_ROTATE_BYTES", 4 * 1024)
    seen = []
    for _ in range(60):
        rec = step_journal.auto_start("s1", tool="bash",
                                      arguments={"command": "true"},
                                      turn_key="T")
        step_journal.auto_finish("s1", rec, ok=True, result="")
        seen.append(rec["step"])
    assert seen == list(range(1, 61))


# ---------------------------------------------------------------------------
# Turns
# ---------------------------------------------------------------------------

def test_a_new_turn_key_opens_the_next_turn(sessions):
    a1 = step_journal.auto_start("s1", tool="bash", turn_key="A")
    step_journal.auto_finish("s1", a1, ok=True)
    a2 = step_journal.auto_start("s1", tool="bash", turn_key="A")
    step_journal.auto_finish("s1", a2, ok=True)
    b1 = step_journal.auto_start("s1", tool="bash", turn_key="B")
    assert (a1["turn"], a1["step"]) == (1, 1)
    assert (a2["turn"], a2["step"]) == (1, 2)
    assert (b1["turn"], b1["step"]) == (2, 1)


def test_the_note_names_only_the_open_steps_of_the_last_turn(sessions):
    """A step left open by an earlier crash was reported then; the note of
    the next crash is about the turn that died."""
    step_journal.auto_start("s1", tool="bash",
                            arguments={"command": "old submit"},
                            turn_key="A")                  # crash 1
    rec = step_journal.auto_start("s1", tool="read_file",
                                  arguments={"path": "x"}, turn_key="B")
    step_journal.auto_finish("s1", rec, ok=True)
    step_journal.auto_start("s1", tool="bash",
                            arguments={"command": "new submit"},
                            turn_key="C")                  # crash 2
    note = step_journal.render_resume_note("s1")
    assert note and "new submit" in note
    assert "old submit" not in note


def test_a_crash_before_the_first_call_of_a_turn_still_names_the_open_step(
        sessions):
    """A restarted session that dies again before any tool call has not
    resolved the open step: the newest journalled turn is still the one."""
    step_journal.auto_start("s1", tool="bash",
                            arguments={"command": "sbatch job.sh"},
                            turn_key="A")
    note = step_journal.render_resume_note("s1")
    assert note and "sbatch job.sh" in note


def test_each_engine_turn_hands_the_executor_a_new_turn_key(tmp_path):
    from delfin.agent import engine as E
    from delfin.agent.api_client import KitToolPermissions, StreamEvent

    keys = []

    class _Client:
        def __init__(self):
            self._permissions = KitToolPermissions(workspace=str(tmp_path))

        def stream_message(self, *_a, **_k):
            keys.append(self._permissions.step_turn_key)
            yield StreamEvent(type="text_delta", text="ok")
            yield StreamEvent(type="message_delta", stop_reason="end_turn",
                              input_tokens=1, output_tokens=1)

    client = _Client()
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = E.AgentEngine(repo_dir=str(tmp_path), backend="api",
                            provider="kit", model="kit.qwen3.5-397b-A17b",
                            mode="solo")
    eng.stream_response("one")
    eng.stream_response("two")
    assert len(keys) >= 2 and all(keys), keys
    assert keys[0] != keys[-1], keys


def test_the_executor_journals_under_the_turn_it_is_given(sessions, tmp_path,
                                                          monkeypatch):
    exc = api_client._DocToolExecutor.__new__(api_client._DocToolExecutor)

    class _Perms:
        task_session_id = "probe"
        step_turn_key = "K1"
        mode = "default"
        pre_tool_hook = None
        post_tool_hook = None
        workspace = str(Path(tmp_path) / "ws")
        read_tracker = {}
        session_allowed_tools = None
        session_denied_tools = None

    perms = _Perms()
    monkeypatch.setattr(
        exc, "_dispatch_under_rules",
        lambda name, arguments, permissions, rule_hit: json.dumps({"ok": 1}),
        raising=False)
    exc._execute_pinned("probe_tool", {}, perms)
    perms.step_turn_key = "K2"
    exc._execute_pinned("probe_tool", {}, perms)
    starts = [e for e in step_journal.read_steps("probe")[0]
              if e["kind"] == "start"]
    assert [(e["turn"], e["step"]) for e in starts] == [(1, 1), (2, 1)]


# ---------------------------------------------------------------------------
# Deletion
# ---------------------------------------------------------------------------

def test_a_deleted_session_takes_its_journal_with_it(sessions, monkeypatch):
    monkeypatch.setattr(step_journal, "_ROTATE_BYTES", 2 * 1024)
    sessions.mkdir(parents=True, exist_ok=True)
    (sessions / "s1.json").write_text("{}")
    for _ in range(40):
        rec = step_journal.auto_start("s1", tool="bash",
                                      arguments={"command": "true"})
        step_journal.auto_finish("s1", rec, ok=True)
    assert list(sessions.glob("s1.steps.jsonl*"))
    assert session_store.delete_session("s1")
    assert not list(sessions.glob("s1.steps.jsonl*"))
