"""A negative exit code is not a failure of the work.

It is POSIX shorthand for "killed by signal |n|" — something outside the
command ended it. Left as a bare ``exit_code: -9`` it reads as an
ordinary failure, and the session has to know the convention to get any
further.

Measured on 2026-09-18: two full test suites on a shared login node were
killed at 24 min and at 205 s, both far inside their own timeouts. The
session reasoned its way to "something on this host kills long runs" and
still had to file it as a supposition, because nothing in the status said
so. The evidence was there; the status simply did not carry it.

Now it does: the signal, its name, and what such a kill usually means —
out of memory, a watchdog on a shared machine, or a kill from elsewhere —
with the explicit warning that this says nothing about the work.
"""

from __future__ import annotations

import pytest

from delfin.agent import bash_jobs as BJ


class _Proc:
    def __init__(self, rc, pid=424242):
        self.returncode = rc
        self.pid = pid

    def poll(self):
        return self.returncode


def _job(rc, command="python -m pytest tests/ -q"):
    job = BJ.BashJob.__new__(BJ.BashJob)
    job.job_id = "01e5b151"
    job.command = command
    job.description = "the full suite"
    job.cwd = "/tmp"
    job.stdout_path = "/tmp/out.log"
    job.stderr_path = "/tmp/err.log"
    job.started_at = 0.0
    job.finished_at = 1440.0
    job.proc = _Proc(rc)
    job.exit_code = rc
    return job


@pytest.fixture(autouse=True)
def _no_children(monkeypatch):
    monkeypatch.setattr(BJ, "_group_children_alive", lambda *a, **k: False)


def test_a_signalled_job_names_the_signal():
    st = _job(-9).status_dict()
    assert st["exit_code"] == -9, "the raw code is still there"
    assert st["killed_by_signal"] == 9
    assert st["signal_name"] == "SIGKILL"
    assert "not the command's own exit status" in st["note"]


def test_the_note_names_the_usual_causes():
    note = _job(-9).status_dict()["note"]
    assert "memory" in note and "watchdog" in note, (
        "a session that has to guess files its conclusion as a guess")


def test_an_ordinary_failure_is_left_alone():
    st = _job(1).status_dict()
    assert st["exit_code"] == 1
    assert "killed_by_signal" not in st
    assert "note" not in st, "exit 1 IS the command's own answer"


def test_success_is_left_alone():
    st = _job(0).status_dict()
    assert st["exit_code"] == 0 and "killed_by_signal" not in st


def test_a_running_job_says_nothing_about_signals():
    st = _job(None).status_dict()
    assert st["running"] is True and "killed_by_signal" not in st


def test_an_unknown_signal_number_still_answers():
    """Never raises: this decorates an error path, and a platform without
    a name for the number must not turn a report into a traceback."""
    assert BJ._signal_name(9) == "SIGKILL"
    assert "129" in BJ._signal_name(129)
    assert BJ._signal_name("nonsense")
