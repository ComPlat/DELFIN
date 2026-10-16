"""Every way out of ``delfin-agent chat`` is time-limited and keeps its save.

Driven through the real front door on a pty (``agent_terminal_probes.
pty_chat``), with faults injected into the child by
``agent_terminal_probes/leave_shim/sitecustomize.py``.

* The terminal gone (SIGHUP, SIGTERM, a hung-up line): the whole leave,
  session-end stages included, ends within the leave deadline
  (``repl._LEAVE_DEADLINE_S``, 10 s).
* An attended leave (Ctrl+D at the prompt, /exit, a second idle Ctrl+C):
  the save keeps that short bound; the session-end stages after it (the
  report, skill learning -- a model call --, indexing, the cycle outcome,
  the episode, the eval loop) get their own budget
  (``agent.leave_stage_budget_s``, 300 s by default), with a line saying
  what is being finished and that Ctrl+C skips it. A stage that takes
  15 s completes; Ctrl+C during a stage skips the rest and leaves at
  once; a stage that hangs is cut at its budget.
* A SIGHUP, SIGTERM or SIGINT that arrives while the leave is saving the
  session does not cut the save. The handlers that turned a signal into
  a leave used to be restored to the defaults before the save, so a
  repeated signal ended the process mid-save.
"""

from __future__ import annotations

import json
import os
import signal
import time

import pytest

pytest.importorskip("pexpect")
pytest.importorskip("termios")

from agent_terminal_probes.pty_chat import REPO, Chat, close_all  # noqa: E402

SHIM = REPO / "tests" / "agent_terminal_probes" / "leave_shim"

#: The session's own deadline plus a margin for a loaded runner.
BOUND_S = 25.0


@pytest.fixture
def shimmed(tmp_path):
    opened: list[Chat] = []

    def _open(mode: str, agent_settings: dict | None = None) -> Chat:
        marks = tmp_path / "marks"
        marks.mkdir(exist_ok=True)
        root = tmp_path / f"s{len(opened)}"
        if agent_settings:
            (root / "home").mkdir(parents=True, exist_ok=True)
            (root / "home" / ".delfin_settings.json").write_text(
                json.dumps({"agent": dict(agent_settings)}))
        session = Chat(root, args=("--permission-mode", "default"),
                       env={"PYTHONPATH": f"{SHIM}{os.pathsep}{REPO}",
                            "LEAVE_SHIM": mode,
                            "LEAVE_SHIM_DIR": str(marks)})
        session.marks = marks
        opened.append(session)
        assert session.wait_for("idle", 60) == "idle", (
            "no prompt within 60 s:\n" + "\n".join(session.visible()))
        for mark in marks.iterdir():
            mark.unlink()
        return session

    yield _open
    close_all(opened)


def _reached_the_stage(session: Chat) -> None:
    assert session.wait_until(
        lambda s: (s.marks / "stage-start").exists(), 30), (
        "the leave never reached the session-end stages")


# ---------------------------------------------------------------------------
# An attended leave finishes its session-end stages, under their budget
# ---------------------------------------------------------------------------

def test_a_slow_session_end_stage_completes_on_ctrl_d(shimmed):
    session = shimmed("slow-stage")
    session.child.send("\x04")
    _reached_the_stage(session)
    assert session.ended_within(60), "still running 60 s after Ctrl+D"
    assert (session.marks / "stage-done").exists(), (
        "a 15 s session-end stage was cut on an attended leave:\n"
        + "\n".join(session.visible()[-6:]))
    assert session.child.exitstatus == 0


def test_ctrl_c_during_the_session_end_skips_it_and_leaves(shimmed):
    session = shimmed("hang-after-save")
    session.child.send("\x04")
    _reached_the_stage(session)
    assert session.wait_until(lambda s: "Ctrl+C skips" in s.text(), 15), (
        "no line says what is being finished:\n"
        + "\n".join(session.visible()[-6:]))
    session.child.send("\x03")
    assert session.ended_within(8), "still running 8 s after Ctrl+C"
    assert session.child.exitstatus == 0


def test_a_hung_session_end_stage_is_cut_at_its_budget(shimmed):
    session = shimmed("hang-after-save",
                      agent_settings={"leave_stage_budget_s": 5})
    session.child.send("\x04")
    t0 = time.monotonic()
    _reached_the_stage(session)
    assert session.ended_within(15 - (time.monotonic() - t0)), (
        "still running 15 s after Ctrl+D with a 5 s stage budget")
    assert "did not finish within 5 s" in session.text(), (
        "\n".join(session.visible()[-6:]))
    assert session.child.exitstatus == 0


# ---------------------------------------------------------------------------
# The terminal gone keeps the short bound
# ---------------------------------------------------------------------------

def test_a_hangup_with_a_hung_stage_ends_within_the_deadline(shimmed):
    session = shimmed("hang-after-save")
    os.kill(session.child.pid, signal.SIGHUP)
    _reached_the_stage(session)
    assert session.ended_within(BOUND_S), (
        f"still running {BOUND_S:g} s after SIGHUP with a hung stage")
    assert session.child.exitstatus == 129


# ---------------------------------------------------------------------------
# A repeated signal does not cut the save
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("repeat", ["SIGHUP", "SIGTERM", "SIGINT"])
def test_a_signal_during_the_leave_does_not_cut_the_save(shimmed, repeat):
    session = shimmed("slow-save")
    session.child.send("\x04")
    assert session.wait_until(
        lambda s: (s.marks / "save-start").exists(), 30), (
        "the leave never started saving the session")
    os.kill(session.child.pid, getattr(signal, repeat))
    assert session.ended_within(BOUND_S)
    assert (session.marks / "save-done").exists(), (
        f"{repeat} during the save ended the process before the save "
        f"finished (exit {session.child.exitstatus}, "
        f"signal {session.child.signalstatus})")
    assert session.child.signalstatus is None
