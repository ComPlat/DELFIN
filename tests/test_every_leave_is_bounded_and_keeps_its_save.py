"""Every way out of ``delfin-agent chat`` is time-limited and keeps its save.

Driven through the real front door on a pty (``agent_terminal_probes.
pty_chat``), with faults injected into the child by
``agent_terminal_probes/leave_shim/sitecustomize.py``.

* Ctrl+D at the prompt leaves within the leave deadline
  (``repl._LEAVE_DEADLINE_S``, 10 s) even when a session-end stage after
  the save hangs. A hangup or SIGTERM was bounded this way; Ctrl+D was
  not, so a hung stage kept the process alive indefinitely.
* A SIGHUP, SIGTERM or SIGINT that arrives while the leave is saving the
  session does not cut the save. The handlers that turned a signal into
  a leave were restored to the defaults before the save, so a repeated
  signal ended the process mid-save, by the signal's default action.
"""

from __future__ import annotations

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

    def _open(mode: str) -> Chat:
        marks = tmp_path / "marks"
        marks.mkdir(exist_ok=True)
        session = Chat(tmp_path / f"s{len(opened)}",
                       args=("--permission-mode", "default"),
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


def test_ctrl_d_at_the_prompt_leaves_within_the_deadline(shimmed):
    session = shimmed("hang-after-save")
    session.child.send("\x04")
    t0 = time.monotonic()
    assert session.wait_until(
        lambda s: (s.marks / "stage-start").exists(), 30), (
        "the leave never reached the session-end stages")
    assert session.ended_within(BOUND_S - (time.monotonic() - t0)), (
        f"still running {BOUND_S:g} s after Ctrl+D with a hung "
        "session-end stage")
    assert session.child.exitstatus == 0


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
