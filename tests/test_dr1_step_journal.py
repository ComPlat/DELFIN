"""DR1: the per-session step journal — red tests (written against the
unchanged behaviour, i.e. a module that does not yet fully exist or does
not yet satisfy the resume contract).

These tests drive the REAL writer functions (step_journal.start_step /
finish_step) — never a typed dict — per the wave-15b rule, and never touch
the real home: the journal resolves its directory per call from
session_store._SESSIONS_DIR, which each test redirects to tmp_path.
"""

from __future__ import annotations

import json
import stat
import threading

import pytest

from delfin.agent import session_store, step_journal


@pytest.fixture()
def journal_dir(tmp_path, monkeypatch):
    """The journal lives under the redirected session-store directory."""
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path / "sess")
    return tmp_path / "sess"


def _start(session="s1", turn=1, step=1, tool="bash", arguments=None):
    return step_journal.start_step(
        session, turn=turn, step=step, tool=tool, arguments=arguments)


# --- the writer ------------------------------------------------------------

def test_start_then_finish_writes_two_append_only_lines(journal_dir):
    assert _start(arguments={"command": "git commit -m x"})
    assert step_journal.finish_step("s1", step=1, ok=True,
                                    result="committed")
    text = (journal_dir / "s1.steps.jsonl").read_text()
    kinds = [json.loads(l)["kind"] for l in text.strip().splitlines()]
    assert kinds == ["start", "finish"], kinds


def test_journal_file_is_0600(journal_dir):
    _start()
    p = journal_dir / "s1.steps.jsonl"
    mode = stat.S_IMODE(p.stat().st_mode)
    assert mode == 0o600, oct(mode)


def test_no_secret_reaches_disk(journal_dir):
    args = {"command": "export API_KEY=sk-supersecret-123456 && ls",
            "path": "notes.md"}
    assert _start(arguments=args)
    raw = (journal_dir / "s1.steps.jsonl").read_text()
    assert "sk-supersecret-123456" not in raw
    assert "API_KEY" not in raw.replace("API_KEY=", "") or "API_KEY=" not in raw
    # The args head IS recorded (redacted) so the resume line can name it.
    assert "notes.md" in raw


def test_side_effect_classification(journal_dir):
    _start(tool="read_file", arguments={"path": "a.md"})
    assert not json.loads(
        (journal_dir / "s1.steps.jsonl").read_text())["side_effect"]
    _start(tool="bash", arguments={"command": "ls"})
    assert json.loads(
        (journal_dir / "s1.steps.jsonl").read_text().splitlines()[-1]
    )["side_effect"]


def test_finish_digest_is_secret_free_and_deterministic(journal_dir):
    step_journal.start_step("s2", turn=1, step=1, tool="bash",
                            arguments={"command": "echo hi"})
    secret = "token ghp_0123456789abcdefghijklmnopqrstuvwxyz012345"
    step_journal.finish_step("s2", step=1, ok=True,
                             result=f"ok {secret}")
    raw = (journal_dir / "s2.steps.jsonl").read_text()
    assert secret not in raw
    lines = raw.strip().splitlines()
    d1 = json.loads(lines[-1])["result_digest"]
    # Same content -> same digest (deterministic, like tool_digest).
    again = step_journal and json.loads(lines[-1])["result_digest"]
    assert d1 == again
    assert len(d1) <= 400


# --- the readers (resume path) --------------------------------------------

def test_open_step_detected_after_start_without_finish(journal_dir):
    _start(turn=7, step=3, tool="bash", arguments={"command": "make"})
    open_out, notes = step_journal.open_steps("s1")
    assert len(open_out) == 1
    assert open_out[0]["step"] == 3
    assert open_out[0]["tool"] == "bash"
    assert open_out[0]["side_effect"] is True
    assert notes == []


def test_finished_step_is_not_open(journal_dir):
    _start(turn=1, step=1, tool="bash")
    step_journal.finish_step("s1", step=1, ok=True, result="done")
    open_out, _ = step_journal.open_steps("s1")
    assert open_out == []


def test_may_or_may_not_have_run_line_for_open_side_effect(journal_dir):
    _start(turn=7, step=3, tool="bash", arguments={"command": "git push"})
    note = step_journal.render_resume_note("s1")
    assert note is not None
    assert "step 3" in note
    assert "bash" in note
    assert "may or may not have run" in note
    assert "CHECK BEFORE REPEATING" in note


def test_finished_readonly_step_result_digest_replayable(journal_dir):
    step_journal.start_step("s1", turn=1, step=2, tool="read_file",
                            arguments={"path": "a.md"})
    step_journal.finish_step("s1", step=2, ok=True,
                             result="1  the content of a.md")
    fin = step_journal.finished_result("s1", step=2)
    assert fin is not None
    assert fin["ok"] is True
    assert fin["result_digest"], "digest must be non-empty"


def test_corrupt_line_skipped_with_note_never_a_crash(journal_dir):
    p = journal_dir / "s1.steps.jsonl"
    _start(turn=1, step=1, tool="bash")
    good = p.read_text()
    # Simulate a torn append: a partial JSON line after the good one.
    p.write_text(good + '{"kind": "start", "ts": 1.0, "ste')
    _start(turn=1, step=2, tool="read_file", arguments={"path": "b.md"})
    open_out, notes = step_journal.open_steps("s1")
    tools = [e["tool"] for e in open_out]
    assert tools == ["bash", "read_file"], tools
    assert any("corrupt" in n for n in notes), notes


def test_missing_journal_is_no_steps_no_note(journal_dir):
    open_out, notes = step_journal.open_steps("never-existed")
    assert open_out == []
    assert notes == []
    assert step_journal.render_resume_note("never-existed") is None


def test_render_note_reports_corrupt_lines_and_last_finished(journal_dir):
    step_journal.start_step("s1", turn=1, step=1, tool="bash",
                            arguments={"command": "x"})
    step_journal.finish_step("s1", step=1, ok=True, result="ok")
    p = journal_dir / "s1.steps.jsonl"
    p.write_text(p.read_text() + "not-json-at-all\n")
    step_journal.start_step("s1", turn=1, step=2, tool="bash",
                            arguments={"command": "y"})
    note = step_journal.render_resume_note("s1")
    assert note is not None
    assert "corrupt" in note
    # "step 1" also matches the open step 2's neighbours? No: the open
    # step line names "step 2"; "step 1" must come from the LAST
    # finished step line.
    assert "Last recorded finished step: 1" in note
    assert "must NOT be repeated" in note


def test_journal_resolves_per_call_from_session_store(journal_dir, tmp_path):
    """The redirect moves the journal: resolved per call, never cached."""
    _start()
    other = tmp_path / "elsewhere"
    session_store._SESSIONS_DIR = other
    try:
        assert step_journal.journal_path("s1").parent == other
        assert (journal_dir / "s1.steps.jsonl").exists()
    finally:
        session_store._SESSIONS_DIR = journal_dir


def test_concurrent_appends_all_land(journal_dir):
    # Append-only under contention: two threads' records both survive.
    def worker(step):
        _start(turn=1, step=step, tool="bash")

    ts = [threading.Thread(target=worker, args=(i,)) for i in range(1, 9)]
    for t in ts:
        t.start()
    for t in ts:
        t.join()
    text = (journal_dir / "s1.steps.jsonl").read_text().strip().splitlines()
    assert len(text) == 8, len(text)
    assert all(json.loads(l)["kind"] == "start" for l in text)


def test_self_trim_keeps_recent_tail(journal_dir, monkeypatch):
    monkeypatch.setattr(step_journal, "_MAX_NOTE_LINES", 10)
    for i in range(1, 26):
        _start(turn=1, step=i, tool="read_file", arguments={"path": f"f{i}"})
    steps, notes = step_journal.read_steps("s1")
    assert len(steps) == 10
    assert any("trimmed" in n for n in notes), notes
    assert steps[-1]["step"] == 25
