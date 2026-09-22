"""Looking inside a background shell was the agent's privilege alone.

``bash_output`` is a tool: the model can read what a long run has
written, and the person who started it could not. ``/bash`` listed jobs
and killed them, and that was all. A run that takes half an hour is
exactly the one somebody wants to watch — and on 2026-09-18 four
sessions spent 154 minutes inside calls longer than five minutes, with
nothing to look at while they did.

``/bash <job_id>`` answers the four questions worth asking: what it is,
whether it still runs, how long it has been going, and what it has
written.

The second half of this file is the end-of-turn line about work nobody
committed. Four sessions that day produced four commits between them,
all at the very end, and two full suites were killed mid-session by a
node watchdog: what is committed survives that. It is a line, not a
demand, and it says nothing outside a git checkout.
"""

from __future__ import annotations

import subprocess

import pytest

from delfin.agent import repl_commands as RC
from delfin.agent import api_client as A


class _Job:
    def __init__(self, job_id, running=True, code=None, killed=None):
        self.job_id = job_id
        self.command = "python -m pytest tests/ -q"
        self.started_at = 0.0
        self._st = {
            "job_id": job_id,
            "running": running,
            "exit_code": code,
            "elapsed_s": 419.0,
            "command": self.command,
            "description": "the full suite",
            "cwd": "/work",
        }
        if killed:
            self._st.update(killed_by_signal=killed, signal_name="SIGKILL",
                            note="killed by SIGKILL (9)")

    def status_dict(self):
        return dict(self._st)


class _Registry:
    def __init__(self, jobs):
        self._jobs = jobs

    def list_jobs(self, include_finished=False):
        return list(self._jobs)


@pytest.fixture()
def bash(monkeypatch):
    def _run(args, jobs, output=None):
        from delfin.agent import bash_jobs as bj
        monkeypatch.setattr(bj, "get_registry", lambda: _Registry(jobs))
        if output is not None:
            monkeypatch.setattr(bj, "read_output", lambda *a, **k: output)
        return RC._bash(None, args).output
    return _run


# -- looking inside ---------------------------------------------------------

def test_a_job_id_shows_what_it_is_doing(bash):
    out = bash("01e5b151", [_Job("01e5b151")],
               output={"stdout": "collected 900 items\n...", "stderr": "",
                       "stdout_total_lines": 2})
    assert "Status:" in out and "running" in out
    assert "Runtime:" in out and "6m" in out, out
    assert "pytest" in out
    assert "collected 900 items" in out, "the output is the point"


def test_a_finished_job_names_its_exit(bash):
    out = bash("01e5b151", [_Job("01e5b151", running=False, code=1)],
               output={"stdout": "1 failed", "stderr": ""})
    assert "finished (exit 1)" in out


def test_a_killed_job_says_it_was_killed(bash):
    out = bash("01e5b151",
               [_Job("01e5b151", running=False, code=-9, killed=9)],
               output={"stdout": "", "stderr": ""})
    assert "killed by SIGKILL" in out


def test_the_peek_stays_small_and_says_how_to_see_more(bash):
    """A running suite writes thousands of lines. Pasting them into the
    chat buries the four facts above them, and the count says more about
    the run than the lines do."""
    out = bash("01e5b151", [_Job("01e5b151")],
               output={"stdout": "\n".join(f"line {i}" for i in range(900)),
                       "stderr": "", "stdout_total_lines": 900})
    body = [ln for ln in out.splitlines() if ln.startswith("    line ")]
    assert len(body) == 5, f"{len(body)} lines pasted into the chat"
    assert "last 5 of 900 lines" in out
    assert "/bash 01e5b151 200   for more of it" in out


def test_asking_for_more_gives_more(bash):
    out = bash("01e5b151 50", [_Job("01e5b151")],
               output={"stdout": "\n".join(f"line {i}" for i in range(900)),
                       "stderr": "", "stdout_total_lines": 900})
    assert len([ln for ln in out.splitlines()
                if ln.startswith("    line ")]) == 50


def test_a_short_output_needs_no_invitation(bash):
    out = bash("01e5b151", [_Job("01e5b151")],
               output={"stdout": "done", "stderr": "", "stdout_total_lines": 1})
    assert "for more of it" not in out


def test_a_job_that_wrote_nothing_says_so(bash):
    out = bash("01e5b151", [_Job("01e5b151")],
               output={"stdout": "", "stderr": ""})
    assert "nothing written yet" in out


def test_an_unknown_id_says_how_to_find_the_right_one(bash):
    out = bash("nope", [_Job("01e5b151")], output={"stdout": "", "stderr": ""})
    assert "no background job" in out and "/bash lists them" in out


def test_the_list_still_lists_and_offers_the_look(bash):
    out = bash("", [_Job("01e5b151")], output={"stdout": "", "stderr": ""})
    assert "01e5b151" in out
    assert "/bash <job_id>" in out, "the new way has to be discoverable"


def test_kill_is_untouched(bash):
    out = bash("kill", [_Job("01e5b151")], output={"stdout": "", "stderr": ""})
    assert "usage: /bash kill <job_id>" in out


# -- the end-of-turn line ---------------------------------------------------

@pytest.fixture()
def repo(tmp_path):
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run(["git", "config", "user.email", "t@example.com"],
                   cwd=tmp_path, check=True)
    subprocess.run(["git", "config", "user.name", "t"], cwd=tmp_path,
                   check=True)
    (tmp_path / "a.py").write_text("x = 1\n", encoding="utf-8")
    subprocess.run(["git", "add", "a.py"], cwd=tmp_path, check=True)
    subprocess.run(["git", "commit", "-qm", "first"], cwd=tmp_path, check=True)
    return tmp_path


def test_a_clean_checkout_says_nothing(repo):
    assert A._uncommitted_note(repo) == ""


def test_changed_and_uncommitted_is_named(repo):
    (repo / "a.py").write_text("x = 2\n", encoding="utf-8")
    note = A._uncommitted_note(repo)
    assert "1 tracked file changed" in note
    assert "control run that fails" in note, "it says which commit is worth most"


def test_an_untracked_file_alone_is_not_a_reminder(repo):
    """A scratch file is not work somebody forgot to commit."""
    (repo / "scratch.log").write_text("noise\n", encoding="utf-8")
    assert A._uncommitted_note(repo) == ""


def test_outside_a_checkout_it_says_nothing(tmp_path):
    assert A._uncommitted_note(tmp_path / "nowhere") == ""


def test_it_never_raises():
    assert A._uncommitted_note(None) == ""
    assert A._uncommitted_note("/proc/nonexistent/x") == ""
