"""F1 crash-resume: a killed worker comes back to its in-flight turn.

Package F1 target 1. When a worker (or the whole supervisor) dies in the
middle of a turn, a restart must bring the session back to the turn it was
in -- from the mid-turn checkpoint and the step journal, not from the start.
An OPEN side-effecting step is never run again automatically: the durable
layer resolves it by CHECKING whether the side effect happened and else
surfaces it as MUST_ASK. A finished step is never repeated.

The records are produced by the REAL writers (step_journal.start_step /
finish_step over the real JSONL, session_store.save_turn_checkpoint over a
redirected store root), and the "crash" is a REAL SIGKILL of a child process
-- not a mock of the kill. The decision function under test lives in
delfin/agent/durable.py (new in this package).
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

from delfin.agent import session_store, step_journal
from delfin.agent import durable

from conftest import child_env


# A side-effecting bash command whose only observable is a file written
# OUTSIDE the scratch tree. The resumed session must check the file, never
# re-run the command automatically.
COMMIT_CMD = 'git commit -am "land target"'
FILE_CMD = "touch /tmp/f1-side-effect.sentinel"


@pytest.fixture
def store_root(tmp_path, monkeypatch):
    """Redirect the session store and journal under one tmp root."""
    root = tmp_path / "store"
    root.mkdir()
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", root)
    return root


def _record_finished_commit(store_root, session):
    """Write a real journal: one finished side-effecting commit step."""
    rec = step_journal.auto_start(session, tool="bash", arguments=COMMIT_CMD, turn=1, turn_key="t1")
    assert rec["recorded"]
    step_journal.auto_finish(session, rec, ok=True, result="done")
    return rec


def _parked_child(store_root, session, env):
    """A child that writes a real checkpoint+journal then parks for SIGKILL.

    ``env`` comes from the suite's ``child_env`` fixture: the child is a
    real process and must resolve every ``~/.delfin`` sink of its own
    outside the real home (CI guard test_a_child_process_writes_into_
    no_real_home.py).

    Returns (proc, checkpoint_payload). The child writes the journal via the
    real step_journal, a checkpoint via the real session_store writer, then
    blocks in a sleep so the parent can SIGKILL it mid-"turn".
    """
    payload = {"task": "the in-flight turn", "turn": 1}
    code = (
        "import sys, time\n"
        "from pathlib import Path\n"
        "import sys; sys.path.insert(0, %r)\n"
        "from delfin.agent import session_store, step_journal\n"
        "store_root = Path(%r)\n"
        "session_store._SESSIONS_DIR = store_root\n"
        "step_journal.start_step(%r, turn=1, step=1, tool='bash', arguments='touch /tmp/f1-side-effect.sentinel', turn_key='t1')\n"
        "session_store.save_turn_checkpoint(%r, %r)\n"
        "sys.stdout.write('parked\\n'); sys.stdout.flush()\n"
        "time.sleep(600)\n"
    ) % (str(Path(__file__).resolve().parents[1]), str(store_root), session, session, {"task": "the in-flight turn", "turn": 1})
    proc = subprocess.Popen(
        [sys.executable, "-c", code],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        cwd=str(Path(__file__).resolve().parents[1]),
        env=env,
    )
    line = proc.stdout.readline().strip()  # type: ignore[union-attr]
    assert line == "parked", f"child did not park, got {line!r}: {proc.stderr.read()}"
    return proc


def test_resume_from_checkpoint_turn_not_start(store_root):
    """A restart resumes the in-flight turn from the checkpoint, not from start."""
    session = "sess-a"
    # a checkpoint says the session was in turn 2
    path = session_store.save_turn_checkpoint(session, {"task": "x", "turn": 2})
    assert path is not None

    plan = durable.decide_recovery(session, root=store_root)
    assert plan.resume_turn == 2
    assert plan.in_flight is True


def test_open_side_effecting_step_never_auto_rerun_file_absent(store_root, tmp_path):
    """An OPEN side-effecting step whose effect is absent is MUST_ASK, never auto-rerun."""
    session = "sess-b"
    out_file = tmp_path / "f1-side-effect.sentinel"
    # real journal: an open bash step that would create out_file, never finished
    step_journal.start_step(session, turn=3, step=1, tool="bash",
                            arguments=f"touch {out_file}", turn_key="t3")
    session_store.save_turn_checkpoint(session, {"task": "x", "turn": 3})
    assert not out_file.exists()

    plan = durable.decide_recovery(session, root=store_root)
    decision = plan.decision_for_open_step(step=1, turn=3)
    assert decision is not None
    assert decision.verdict is durable.RecoveryVerdict.MUST_ASK
    assert decision.auto_rerun is False
    assert plan.would_block_rerun(session, tool="bash", step=1, turn=3) is True


def test_open_side_effecting_step_effect_present_never_repeated(store_root, tmp_path):
    """An OPEN step whose side effect IS present is VERIFIED_DONE: never re-run."""
    session = "sess-c"
    out_file = tmp_path / "f1-side-effect.sentinel"
    out_file.write_text("made it")
    step_journal.start_step(session, turn=4, step=1, tool="bash",
                            arguments=f"touch {out_file}", turn_key="t4")
    session_store.save_turn_checkpoint(session, {"task": "x", "turn": 4})

    def checker(step):
        # a real check: does the file exist?
        return "present" if out_file.exists() else "absent"

    plan = durable.decide_recovery(session, root=store_root, checkers={"bash": checker})
    decision = plan.decision_for_open_step(step=1, turn=4)
    assert decision is not None
    assert decision.verdict is durable.RecoveryVerdict.VERIFIED_DONE
    assert decision.auto_rerun is False
    assert plan.would_block_rerun(session, tool="bash", step=1, turn=4) is True


def test_finished_step_never_repeated(store_root, tmp_path):
    """A step that finished is never repeated, whatever its tool."""
    session = "sess-d"
    rec = _record_finished_commit(store_root, session)

    plan = durable.decide_recovery(session, root=store_root)
    # the finished step is not an OPEN step: no decision, and no block listing it
    assert plan.decision_for_open_step(step=rec["step"], turn=rec["turn"]) is None
    # a finished step appears only as already-done, never as something to rerun
    for d in plan.decisions:
        if d.step == rec["step"] and d.turn == rec["turn"]:
            assert d.verdict is durable.RecoveryVerdict.FINISHED
    # a finished step is NEVER re-dispatched: would_block_rerun is True even
    # though decision_for_open_step is None (it is not open, so no recovery
    # decision; but the guard must still refuse to run it again). Review
    # finding D: the phase-1 guard answered False here.
    assert plan.would_block_rerun(session, tool="bash",
                                  step=rec["step"], turn=rec["turn"]) is True


def test_real_sigkill_mid_turn_recovers_via_checkpoint(store_root, tmp_path):
    """A REAL SIGKILL of a child mid-turn: the durable plan resumes, asks."""
    session = "sess-kill"
    proc = _parked_child(store_root, session, child_env(tmp_path))
    # real kill: SIGKILL the child process while it is parked (simulating the
    # moment the worker dies inside an open side-effecting step)
    proc.send_signal(signal.SIGKILL)
    proc.wait()
    assert proc.returncode == -signal.SIGKILL

    # the journal has an OPEN step (started, never finished) and a checkpoint
    open_steps, _ = step_journal.open_steps(session, root=store_root)
    assert len(open_steps) == 1
    cp = session_store.load_turn_checkpoint(session)
    assert cp is not None and cp.get("turn") == 1

    plan = durable.decide_recovery(session, root=store_root)
    assert plan.in_flight is True
    # the sole open step is side-effecting and unverifiable -> MUST_ASK
    decision = plan.decision_for_open_step(step=1, turn=1)
    assert decision is not None
    assert decision.verdict is durable.RecoveryVerdict.MUST_ASK
    assert decision.auto_rerun is False


def test_read_only_open_step_may_rerun(store_root, tmp_path):
    """A read-only OPEN step (no side effect) is NOT auto-blocked."""
    session = "sess-ro"
    step_journal.start_step(session, turn=5, step=1, tool="read_file",
                            arguments={"path": str(tmp_path / "x.txt")}, turn_key="t5")
    session_store.save_turn_checkpoint(session, {"task": "x", "turn": 5})

    plan = durable.decide_recovery(session, root=store_root)
    decision = plan.decision_for_open_step(step=1, turn=5)
    assert decision is not None and decision.side_effecting is False
    assert plan.would_block_rerun(session, tool="read_file", step=1, turn=5) is False


# --- Phase 2 hardening (F1, agent/s1-f1b) -----------------------------------
# Each test was written to be RED on the phase-1 code (6d1789304) and is only
# green with the phase-2 hardening in delfin/agent/durable.py.

def test_stale_checkpoint_never_drags_resume_turn_backwards(store_root):
    """A stale checkpoint (turn 1) cannot pull the resume point below the
    journal's highest completed turn (3)."""
    session = "sess-stale-cp"
    step_journal.start_step(session, turn=3, step=1, tool="read_file",
                            arguments={"path": "x"}, turn_key="t3")
    step_journal.finish_step(session, turn=3, step=1, tool="read_file", ok=True)
    # a stale checkpoint left behind from an earlier crash says turn 1
    session_store.save_turn_checkpoint(session, {"task": "old", "turn": 1})

    plan = durable.decide_recovery(session, root=store_root)
    assert plan.resume_turn == 3  # not 1


def test_corrupted_side_effect_flag_never_unlocks_rerun(store_root, tmp_path):
    """A journal record lying about side_effect=False must not unlock an
    automatic re-run of a tool the classifier calls side-effecting."""
    session = "sess-lie"
    out_file = tmp_path / "f1-lied-sentinel"
    # Fabricate a CORRUPTED start record: the STORED flag says side_effect
    # False, but the tool is "bash", which the side-effect classifier calls
    # side-effecting. The stored flag is not the authority -- the classifier
    # is -- so this must still read as side-effecting and be blocked.
    rec = {
        "kind": "start",
        "ts": 0.0,
        "turn": 2,
        "step": 1,
        "tool": "bash",
        "side_effect": False,  # the lie
        "args_digest": f"touch {out_file}",
    }
    journal = step_journal.journal_path(session, root=store_root)
    journal.parent.mkdir(parents=True, exist_ok=True)
    with journal.open("a") as fh:
        fh.write(json.dumps(rec, sort_keys=True) + "\n")
    session_store.save_turn_checkpoint(session, {"task": "x", "turn": 2})

    plan = durable.decide_recovery(session, root=store_root)
    decision = plan.decision_for_open_step(step=1, turn=2)
    assert decision is not None
    assert decision.side_effecting is True  # classifier wins over the lie
    assert plan.would_block_rerun(session, tool="bash", step=1, turn=2) is True


def test_broken_checker_degrades_to_unknown_not_present(store_root):
    """A side-effect checker that raises degrades to 'unknown' (MUST_ASK),
    never to 'absent' (which would unlock a re-run)."""
    session = "sess-broken"
    step_journal.start_step(session, turn=4, step=1, tool="bash",
                            arguments="git commit -am x", turn_key="t4")
    session_store.save_turn_checkpoint(session, {"task": "x", "turn": 4})

    def broken(step):
        raise RuntimeError("corrupt side-effect index")

    plan = durable.decide_recovery(session, root=store_root,
                                   checkers={"bash": broken})
    decision = plan.decision_for_open_step(step=1, turn=4)
    assert decision is not None
    assert decision.verdict is durable.RecoveryVerdict.MUST_ASK
    assert decision.auto_rerun is False
