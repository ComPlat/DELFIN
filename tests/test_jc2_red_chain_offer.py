"""JC2 target 3: the chain offer in the crash-recovery note — red first.

When an agent session restarts after a crash, the recovery note must
list that session's chains that are NOT finished: name, work dir, one
status line per unfinished step, and the exact command
``delfin chain resume <dir>``.  It must say the chain is NOT resumed
automatically and that the agent must ask the user first.  Facts only:
a chain whose state file is gone or fully finished is not offered; a
session without chains gets no chain text at all; an error in the
chain part never breaks the resume (best-effort, like the step-journal
part).

Control: red on c4591afa4 — consume_crash_recovery_note knows nothing
about chains.  Nothing resumes on its own (target 4): the restart path
only produces TEXT; pinned by the same assertions (no state change)
and test_jc2_adversarial_no_auto_resume.py.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest


@pytest.fixture()
def home(tmp_path, monkeypatch):
    monkeypatch.setattr("pathlib.Path.home", lambda: tmp_path / "home")
    import delfin.agent.session_store as ss
    monkeypatch.setattr(
        ss, "_SESSIONS_DIR", tmp_path / "home" / ".delfin" / "agent_sessions")
    return tmp_path / "home"


def _fake_run_step(step_name, geometry=None, work_dir=None, **kwargs):
    from delfin.tools.pipeline import StepResult, StepStatus

    out = Path(work_dir) / f"{step_name}.out"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("fake run output\n")
    return StepResult(step_name=step_name, status=StepStatus.SUCCESS,
                      work_dir=Path(work_dir))


def _make_session_file(home, session_id):
    """A minimal saved session so load_session finds it and the
    checkpoint is newer than its save (the note's precondition)."""
    import delfin.agent.session_store as ss
    ss.save_session(session_id, title="t", messages=[])
    ck_dir = home / ".delfin" / "agent_sessions"
    ck_dir.mkdir(parents=True, exist_ok=True)
    ck_path = ck_dir / f"{session_id}.turn.json"
    ck_path.write_text(json.dumps({
        "ts": 9e9,  # far newer than any save
        "user_message": "run the chain",
        "tool_calls": 2,
        "partial_response": "",
    }), encoding="utf-8")


def _run_chain(wd, monkeypatch, session_id, *, finish=False):
    """A real chain under *session_id*; the first step SUCCEEDS (the
    fake writes output) and, when finish=False (the crash scenario), a
    second step is left UNFINISHED — as a chain killed mid-way looks."""
    monkeypatch.setattr("delfin.tools._runner.run_step", _fake_run_step)
    from delfin.tools.pipeline import Pipeline

    geom = Path(wd) / "start.xyz"
    geom.parent.mkdir(parents=True, exist_ok=True)
    geom.write_text("2\n\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8")
    pipe = Pipeline(name="demo")
    pipe.add("xtb_opt", charge=0)
    if not finish:
        # A step whose run crashes: the process died inside it —
        # persisted as submitted/running, never finished.  Only the
        # SECOND step crashes (a chain killed midway), the first
        # finishes cleanly.
        def selective_run_step(step_name, geometry=None, work_dir=None,
                               **kwargs):
            if step_name == "orca_sp":
                raise RuntimeError("simulated node death mid-step")
            return _fake_run_step(step_name, geometry=geometry,
                                  work_dir=work_dir, **kwargs)
        monkeypatch.setattr("delfin.tools._runner.run_step",
                            selective_run_step)
    pipe.add("orca_sp", charge=0, method="PBE0", basis="def2-SVP")
    if not finish:
        import pytest
        with pytest.raises(RuntimeError):
            pipe.run(work_dir=wd, geometry=geom,
                     agent_session_id=session_id)
    else:
        pipe.run(work_dir=wd, geometry=geom, agent_session_id=session_id)
    return wd


def test_offer_lists_unfinished_chains_with_resume_command(
        home, monkeypatch):
    import delfin.agent.session_store as ss

    session_id = "sess-crash"
    wd = _run_chain(home / "runs" / "a", monkeypatch, session_id)
    _make_session_file(home, session_id)
    note = ss.consume_crash_recovery_note(session_id)
    assert note is not None
    assert "[recovered]" in note
    # The offer:
    assert "unfinished chain" in note
    assert "demo" in note
    assert str(wd) in note
    assert f"delfin chain resume {wd}" in note
    # Not automatic:
    assert "NOT resumed automatically" in note
    assert "ask the user first" in note
    # One status line per unfinished step:
    assert "orca_sp" in note


def test_finished_chain_is_not_offered(home, monkeypatch):
    import delfin.agent.session_store as ss

    session_id = "sess-done"
    wd = _run_chain(home / "runs" / "ok", monkeypatch, session_id,
                    finish=True)
    _make_session_file(home, session_id)
    note = ss.consume_crash_recovery_note(session_id)
    assert note is None or "unfinished chain" not in note


def test_gone_chain_is_not_offered(home, monkeypatch):
    import delfin.agent.session_store as ss
    from delfin.tools import chain_index

    session_id = "sess-gone"
    chain_index.record("gone", home / "runs" / "vanished",
                       agent_session_id=session_id)
    _make_session_file(home, session_id)
    note = ss.consume_crash_recovery_note(session_id)
    assert note is None or "unfinished chain" not in note


def test_session_without_chains_gets_no_chain_text(home, monkeypatch):
    import delfin.agent.session_store as ss

    session_id = "sess-plain"
    _make_session_file(home, session_id)
    note = ss.consume_crash_recovery_note(session_id)
    assert note is not None
    assert "unfinished chain" not in note.lower()
    assert "[unfinished chains]" not in note


def test_other_sessions_chains_are_not_offered(home, monkeypatch):
    import delfin.agent.session_store as ss

    _run_chain(home / "runs" / "a", monkeypatch, session_id="sess-1")
    _make_session_file(home, "sess-2")
    note = ss.consume_crash_recovery_note("sess-2")
    assert note is None or "unfinished chain" not in note


def test_chain_error_never_breaks_the_note(home, monkeypatch):
    import delfin.agent.session_store as ss

    session_id = "sess-err"
    _make_session_file(home, session_id)
    # A poisoned index: record() entries unreadable -> the chain part
    # must fail silently while the note itself survives.
    p = home / ".delfin" / "chain_index.jsonl"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("{torn\n", encoding="utf-8")
    note = ss.consume_crash_recovery_note(session_id)
    assert note is not None
    assert "[recovered]" in note


def test_nothing_resumes_on_its_own(home, monkeypatch):
    """Target 4 pin through the restart path: the note is TEXT only —
    after consuming it, every chain's state file is byte-identical."""
    import delfin.agent.session_store as ss

    session_id = "sess-pin"
    wd = _run_chain(home / "runs" / "a", monkeypatch, session_id)
    state_file = wd / "chain_state.json"
    before = state_file.read_text(encoding="utf-8")
    _make_session_file(home, session_id)
    ss.consume_crash_recovery_note(session_id)
    assert state_file.read_text(encoding="utf-8") == before, (
        "the restart path must only produce text — never touch a chain")
