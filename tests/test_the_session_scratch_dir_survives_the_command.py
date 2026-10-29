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
import sys
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


def _cage_runs() -> bool:
    """Whether a real bwrap cage can actually run ON THIS HOST.

    ``requires_a_working_sandbox`` only checks that bwrap *could* confine
    (the minimal probe); a host inside a nested container declares bwrap
    functional and then refuses every real wrap at run time — the full DELFIN
    argv wraps the command in an inner Landlock/socket-guard helper, and the
    host's seccomp rejects it with `[Errno 16] seccomp filter: Device or
    resource busy` (the command exits 126 with a DELFIN refusal). The minimal
    probe cannot see that, so gate on the REAL argv: build it and try to run
    it. On such a host the executed-command property cannot be measured, so
    the tests assert the argv instead.
    """
    try:
        from delfin.agent.api_client import _bwrap_functional
        if not _bwrap_functional():
            return False
        ws = Path(A._session_scratch_dir()).parent
        argv = _bash_argv("true", ws)
        if not argv or argv[0] != "bwrap":
            return False
        r = subprocess.run(argv, capture_output=True, timeout=30)
        return r.returncode == 0
    except Exception:
        return False


def _bash_argv(cmd: str, ws: Path) -> list:
    """The argv the executor would run for a bash command (real cage)."""
    return A._bash_isolation_argv(cmd, str(ws), _perms(ws))


def _argv_binds_scratch(argv: list, scratch: str) -> bool:
    """True when argv binds scratch read-write and sets TMPDIR/TMP/TEMP at it.

    ``--setenv`` alternates KEY VALUE and appears many times (TMPDIR/TMP/TEMP
    and each package-cache dir); accumulating into a dict keyed by the literal
    ``--setenv`` (as a ``{argv[i]: argv[i+1]}`` comprehension would) keeps only
    the last pair. Walk the argv instead and read the KEY and VALUE in step.
    """
    ks = list(argv)
    writes = [ks[i + 1] for i, a in enumerate(ks) if a in ("--bind-try", "--bind")
              and i + 1 < len(ks)]
    if not any(w == scratch for w in writes):
        return False
    env: dict[str, str] = {}
    for i, a in enumerate(ks):
        if a == "--setenv" and i + 2 < len(ks):
            env[ks[i + 1]] = ks[i + 2]
    return env.get("TMPDIR") == scratch and env.get("TMP") == scratch \
        and env.get("TEMP") == scratch


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
    scratch = A._session_scratch_dir()
    argv = _bash_argv('echo "$TMPDIR"', ws)
    argv_carries_scratch(ws, argv, scratch)
    if _cage_runs():
        out = _run_caged("echo \"$TMPDIR|$TMP|$TEMP\"", ws)
        assert out.get("returncode", 0) == 0, out
        parts = out.get("stdout", "").strip().split("|")
        for p in parts:
            assert p == scratch, (parts, scratch)


def argv_carries_scratch(ws, argv, scratch):
    """Assert the executor's real argv binds scratch read-write with TMPDIR at it."""
    assert _argv_binds_scratch(argv, scratch), (argv, scratch)


@requires_a_working_sandbox
def test_file_written_by_one_command_survives_to_the_next():
    ws = Path(A._session_scratch_dir()).parent
    scratch = A._session_scratch_dir()
    # The argv DELFIN would run binds the scratch read-write and points
    # TMPDIR at it (the security property this change controls).
    argv_carries_scratch(ws, _bash_argv("true", ws), scratch)
    # Where the cage can actually run, prove the persistence end to end.
    if _cage_runs():
        first = _run_caged("printf 'kept' > \"$TMPDIR/note.txt\"", ws)
        assert first.get("returncode", 0) == 0, first
        second = _run_caged("cat \"$TMPDIR/note.txt\"", ws)
        assert "kept" in second.get("stdout", ""), second
        content = _read_via_tool(ws, scratch + "/note.txt")
        assert "kept" in content


def test_another_process_cannot_see_the_same_scratch():
    """A fresh process (a different session) gets a different scratch tree."""
    mine = Path(A._session_scratch_dir()).resolve()
    # A child process of this one is a separate session process: its scratch
    # is its own private temp dir, not this one. Run it with the SAME
    # interpreter and PYTHONPATH so it imports this tree's delfin ('-I' would
    # drop PYTHONPATH, so it is not used).
    code = (
        "from delfin.agent import api_client as A;"
        "print(A._session_scratch_dir())"
    )
    env = dict(os.environ)
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1])
    proc = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True,
        env=env, cwd=os.getcwd())
    assert proc.returncode == 0, proc.stderr
    theirs = proc.stdout.strip()
    assert theirs and Path(theirs).resolve() != mine, (mine, theirs, proc.stderr)


@requires_a_working_sandbox
def test_write_gate_admits_the_scratch_dir(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    scratch = A._session_scratch_dir()
    # write_file returns the human-readable diff, not JSON. A path outside
    # the workspace roots is refused before the write: admission means the
    # write happened (the scratch dir is on all_workspace_roots).
    out = _doc_executor.execute(
        "write_file",
        {"path": os.path.join(scratch, "agent.txt"),
         "content": "hello"}, _perms(ws))
    assert "error" not in out, out
    content = (Path(scratch) / "agent.txt").read_text()
    assert content == "hello", content


@requires_a_working_sandbox
def test_landlock_cage_grants_the_scratch_dir():
    """landlock_exec.main keeps the LAST --tmpdir (tmpdir = value, then
    sets TMPDIR/TMP/TEMP and adds it to the write roots). The argv must
    therefore have exactly one --tmpdir and it must be the scratch dir."""
    plain = ["bash", "-c", "echo ok"]
    argv = A._landlock_argv(plain, _perms(Path("/tmp")), extra_write=())
    scratch = A._session_scratch_dir()
    tmpdirs = [argv[i + 1] for i, a in enumerate(argv) if a == "--tmpdir"]
    assert tmpdirs == [scratch], (tmpdirs, scratch)
    write_args = [argv[i + 1] for i, a in enumerate(argv) if a == "--write"]
    assert scratch in write_args, (scratch, write_args)


@requires_a_working_sandbox
def test_tmp_is_still_fresh_a_file_there_does_not_survive():
    ws = Path(A._session_scratch_dir()).parent
    # /tmp stays a fresh tmpfs: the bwrap argv makes /tmp a tmpfs and never
    # binds the scratch there, so a file written to /tmp does not persist.
    argv = _bash_argv("true", ws)
    assert "--tmpfs" in argv, argv
    # Where the cage actually runs, prove /tmp stays fresh end to end.
    if _cage_runs():
        first = _run_caged("printf 'ephemeral' > /tmp/gone.txt", ws)
        assert first.get("returncode", 0) == 0, first
        second = _run_caged("cat /tmp/gone.txt", ws)
        stdout = second.get("stdout", "") + second.get("stderr", "")
        assert "ephemeral" not in stdout, second
