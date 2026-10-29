"""A session scratch directory holds what one command of a session wrote so
the NEXT command of the SAME session and the read tool can use it, and no
other session can see it; /tmp itself stays private and fresh.

Input: the scratch root built per session (``_session_scratch_dir``), the
cage argv builders (bwrap in ``delfin/agent/sandbox.py`` and the api_client
bash cage, Landlock, Seatbelt) and the file-tool root checks.

Output measured by running the real cage on this host (bwrap and Landlock
are functional here; skipped where not):

* a file one caged command writes to ``$TMPDIR`` survives to the next
  command of the same session, and ``read_file`` can read it;
* ``TMPDIR``/``TMP``/``TEMP`` inside every cage point at the scratch dir;
* the scratch dir is 0700 and named after the session's own private temp
  dir (never shared between two processes);
* a write tool call into the scratch dir is admitted by the write gate;
* /tmp remains a fresh tmpfs (a file in /tmp does NOT survive).
"""
from __future__ import annotations

import json
import os
import stat
import subprocess
from pathlib import Path

import pytest

from conftest import requires_a_working_sandbox
from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def landlock_host():
    if not A._landlock_functional():
        pytest.skip("no Landlock ABI on this kernel")
    yield


def _perms(ws: Path, **kw) -> KitToolPermissions:
    return KitToolPermissions(mode="default", workspace=str(ws), **kw)


def _run_caged(cmd: str, ws: Path) -> dict:
    """Run a bash command through the executor's bash gate (real cage)."""
    out = json.loads(_doc_executor.execute(
        "bash", {"command": cmd}, _perms(ws)))
    return out


def _read_via_tool(ws: Path, rel: str) -> str:
    out = json.loads(_doc_executor.execute(
        "read_file", {"path": rel}, _perms(ws)))
    assert "error" not in out, out
    return out.get("content", out.get("data", ""))


def test_scratch_dir_exists_private_and_per_process():
    d = Path(A._session_scratch_dir())
    assert d.is_dir()
    mode = stat.S_IMODE(os.stat(d).st_mode)
    assert mode == 0o700, mode
    # Under the process's private temp dir: another process's scratch is a
    # different tree, so two sessions never share one.
    assert A._private_tmp_dir() in str(d.resolve())


def test_tmpdir_points_at_scratch_in_the_cage():
    ws = Path(A._session_scratch_dir()).parent.parent  # any existing dir
    out = _run_caged("echo \"$TMPDIR|$TMP|$TEMP\"", ws)
    assert out.get("returncode", 0) == 0, out
    line = out.get("stdout", "")
    parts = line.strip().split("|")
    scratch = A._session_scratch_dir()
    for p in parts:
        assert p == scratch, (parts, scratch)


@requires_a_working_sandbox
def test_file_written_by_one_command_survives_to_the_next():
    ws = Path(A._session_scratch_dir()).parent
    first = _run_caged("printf 'kept' > \"$TMPDIR/note.txt\"", ws)
    assert first.get("returncode", 0) == 0, first
    second = _run_caged("cat \"$TMPDIR/note.txt\"", ws)
    assert "kept" in second.get("stdout", ""), second
    # And the read tool can see it too, via a relative path into the scratch
    # is not possible (outside the workspace), so use the absolute one:
    content = _read_via_tool(ws, A._session_scratch_dir() + "/note.txt")
    assert "kept" in content


def test_another_process_cannot_see_the_same_scratch():
    """A fresh process (a different session) gets a different scratch tree."""
    mine = Path(A._session_scratch_dir()).resolve()
    # A child process of this one is a separate session process: its scratch
    # is its own private temp dir, not this one.
    code = (
        "from delfin.agent import api_client as A;"
        "print(A._session_scratch_dir())"
    )
    env = dict(os.environ)
    proc = subprocess.run(
        ["python", "-I", "-c", code], capture_output=True, text=True,
        env=env, cwd=os.getcwd())
    theirs = proc.stdout.strip()
    assert theirs and Path(theirs).resolve() != mine, (mine, theirs)


@requires_a_working_sandbox
def test_write_gate_admits_the_scratch_dir(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    scratch = A._session_scratch_dir()
    out = json.loads(_doc_executor.execute(
        "write_file",
        {"path": os.path.join(scratch, "agent.txt"),
         "content": "hello"}, _perms(ws)))
    assert "error" not in out, out
    assert (Path(scratch) / "agent.txt").read_text() == "hello"


@requires_a_working_sandbox
def test_landlock_cage_grants_the_scratch_dir():
    plain = ["bash", "-c", "echo ok"]
    argv = A._landlock_argv(plain, _perms(Path("/tmp")), extra_write=())
    scratch = A._session_scratch_dir()
    assert any("--tmpdir" == a and scratch == argv[i + 1]
               for i, a in enumerate(argv) if a == "--tmpdir"), argv
    write_args = [argv[i + 1] for i, a in enumerate(argv) if a == "--write"]
    assert scratch in write_args, (scratch, write_args)


@requires_a_working_sandbox
def test_tmp_is_still_fresh_a_file_there_does_not_survive():
    ws = Path(A._session_scratch_dir()).parent
    first = _run_caged("printf 'ephemeral' > /tmp/gone.txt", ws)
    assert first.get("returncode", 0) == 0, first
    second = _run_caged("cat /tmp/gone.txt", ws)
    stdout = second.get("stdout", "") + second.get("stderr", "")
    assert "ephemeral" not in stdout, second
