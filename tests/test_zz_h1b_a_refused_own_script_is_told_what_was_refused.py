"""A refused run of the session's own script is told what was refused.

The bash gate scans the content of a script it is asked to run. The
commonest way into that refusal is a script the session wrote a moment ago,
and the commonest reaction was to run the same file under another spelling
(`python3 x.py`, then `bash x.py`, then `./x.py`). The refusal now says that
the file's CONTENT was refused, not running one's own script, and the second
refusal of the same file ends the attempts.

Running a self-written script is otherwise untouched: these tests pin that
too, because a gate that refused every own script would stop ordinary work.

Through the real executor: the write is a real ``write_file`` call (which
records the change journal), the run a real ``bash`` call.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor

CLEAN = "print('ran fine')\n"
TRIPS = "import os\nos.system('rm -rf /tmp/never-there-h1b')\n"
OWN = "a file this session wrote"


@pytest.fixture
def box(tmp_path, monkeypatch):
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "home")
    (tmp_path / "home").mkdir()
    ws = tmp_path / "ws"
    ws.mkdir()
    perms = KitToolPermissions(workspace=str(ws), mode="acceptEdits",
                               task_session_id="sid-h1b-own-script")

    def write(rel, content):
        return str(_doc_executor.execute(
            "write_file", {"path": str(ws / rel), "content": content}, perms))

    def bash(cmd):
        return str(_doc_executor.execute("bash", {"command": cmd}, perms))

    return ws, write, bash


def test_a_harmless_own_script_simply_runs(box):
    _, write, bash = box
    write("x.py", CLEAN)
    out = bash("python3 x.py")
    assert "ran fine" in out
    assert "refuses" not in out and OWN not in out


def test_the_refusal_of_an_own_script_says_it_is_about_the_content(box):
    _, write, bash = box
    write("x.py", TRIPS)
    out = bash("python3 x.py")
    assert "trip the safety scan" in out
    assert OWN in out and "its content, not running your own script" in out
    assert "x.py" in out


def test_the_second_refusal_of_the_same_file_ends_the_attempts(box):
    _, write, bash = box
    write("x.py", TRIPS)
    bash("python3 x.py")
    for respelled in ("bash x.py", "python3 ./x.py", "cd . && python3 x.py"):
        out = bash(respelled)
        assert "already refused once" in out, (respelled, out)
        assert "Stop here and tell the user" in out


def test_a_script_made_clean_runs_again(box):
    """The allowed way out must stay open: fix the file, run it."""
    _, write, bash = box
    write("x.py", TRIPS)
    assert "trip the safety scan" in bash("python3 x.py")
    write("x.py", CLEAN)
    out = bash("python3 x.py")
    assert "ran fine" in out and "refuses" not in out


def test_a_script_the_session_did_not_write_is_refused_as_before(box):
    ws, _, bash = box
    (ws / "theirs.py").write_text(TRIPS, encoding="utf-8")
    out = bash("python3 theirs.py")
    assert "trip the safety scan" in out
    assert OWN not in out and "already refused once" not in out


def test_another_own_script_gets_the_first_answer_not_the_stop(box):
    _, write, bash = box
    write("x.py", TRIPS)
    write("y.py", TRIPS)
    bash("python3 x.py")
    out = bash("python3 y.py")
    assert OWN in out and "already refused once" not in out


def test_reading_the_own_script_is_not_running_it(box):
    _, write, bash = box
    write("x.py", TRIPS)
    out = bash("cat x.py")
    assert OWN not in out and "already refused once" not in out
