"""Code a running background job uses is not edited under it.

Input: a background job started through the real registry and a write
through the real executor. Output under test: the write gate's answer
and the run note.

Pinned here:

* a file the job's command names (the script, a ``path::node`` test id,
  the tests directory, the module of ``-m``) is refused while the job
  runs, with the job's id and the three ways on;
* the same write goes through once the job has ended;
* a file elsewhere in the job's checkout is written and produces one
  note per job;
* a job of another session does not lock this session's files;
* a file outside every running job is untouched by any of this.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time

import pytest

from delfin.agent import bash_jobs as BJ
from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def repo(tmp_path):
    root = tmp_path / "repo"
    (root / "tests").mkdir(parents=True)
    (root / "pkg").mkdir()
    (root / "pkg" / "__init__.py").write_text("")
    (root / "pkg" / "tool.py").write_text("X = 1\n")
    (root / "run.py").write_text("import time\ntime.sleep(60)\n")
    (root / "tests" / "test_a.py").write_text("def test_a():\n    pass\n")
    (root / "notes.py").write_text("Y = 2\n")
    subprocess.run(["git", "init", "-q", str(root)], check=True)
    BJ._CHECKOUT_TOPS.clear()
    BJ._IN_USE_NOTED.clear()
    BJ.take_in_use_notes()
    return root


def _start(root, command, session_id="s-lock"):
    # The job's process is a sleep; its COMMAND is what names the files.
    job = BJ.get_registry().start(
        command=command, cwd=str(root), workspace=str(root),
        session_id=session_id,
        argv=[sys.executable, "-c", "import time; time.sleep(60)"])
    return job


def _stop(job):
    BJ.get_registry().kill(job.job_id)
    deadline = time.time() + 10
    while job.poll() is None and time.time() < deadline:
        time.sleep(0.05)


def _perms(root, sid="s-lock"):
    p = KitToolPermissions(workspace=str(root), mode="acceptEdits")
    p.task_session_id = sid
    return p


def _write(root, rel, perms, text="Z = 3\n"):
    _doc_executor.execute("read_file", {"path": str(root / rel)}, perms)
    return _doc_executor.execute(
        "write_file", {"path": str(root / rel), "content": text}, perms)


@pytest.mark.parametrize("command,rel", [
    ("python run.py", "run.py"),
    ("python -m pytest tests/test_a.py::test_a -q", "tests/test_a.py"),
    ("pytest tests -q", "tests/test_a.py"),
    ("python -m pkg.tool --fast", "pkg/tool.py"),
])
def test_a_file_the_job_names_is_refused_while_it_runs(repo, command, rel):
    job = _start(repo, command)
    try:
        out = _write(repo, rel, _perms(repo))
        assert f"is used by background job {job.job_id}" in out, out
        assert "bash_status" in out and "bash_kill" in out
        assert "worktree" in out
        assert (repo / rel).read_text() != "Z = 3\n"
    finally:
        _stop(job)
    out = _write(repo, rel, _perms(repo))
    assert "background job" not in out, out
    assert (repo / rel).read_text() == "Z = 3\n"


def test_a_file_elsewhere_in_the_checkout_gets_one_note(repo):
    job = _start(repo, "python run.py")
    try:
        perms = _perms(repo)
        out = _write(repo, "notes.py", perms)
        assert "background job" not in out
        assert (repo / "notes.py").read_text() == "Z = 3\n"
        _write(repo, "pkg/__init__.py", perms, "Q = 1\n")
        notes = BJ.take_in_use_notes()
        assert len(notes) == 1
        assert notes[0].startswith("[running job]")
        assert job.job_id in notes[0]
    finally:
        _stop(job)


def test_a_job_of_another_session_does_not_lock(repo):
    job = _start(repo, "python run.py", session_id="someone-else")
    try:
        out = _write(repo, "run.py", _perms(repo))
        assert "background job" not in out
    finally:
        _stop(job)


def test_job_using_answers_none_outside_every_job(repo, tmp_path):
    other = tmp_path / "elsewhere.py"
    other.write_text("x")
    job = _start(repo, "python run.py")
    try:
        assert BJ.job_using(other, session_id="s-lock") is None
        use = BJ.job_using(repo / "run.py", session_id="s-lock")
        assert use and use["kind"] == "named"
        assert json.dumps(use)
    finally:
        _stop(job)
