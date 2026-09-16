"""A command the agent runs ends with everything it started, on every host.

Without the bubblewrap process cage -- a CI runner, an HPC node with user
namespaces off, macOS -- subprocess.run killed only bash at a timeout, so
an ``&`` job, a pipeline member or a nohup survived the command and the
stop, and the command shared the agent's terminal (audit 2026-09-16).
These tests run the host-independent floor with the cage switched off.
"""
import json
import os
import subprocess
import sys
import time

import pytest

from delfin.agent import api_client as A
from delfin.agent import contained_run as C


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    try:
        with open(f"/proc/{pid}/stat") as fh:
            return fh.read().split(")")[-1].split()[0] != "Z"
    except OSError:
        return True


def _wait_dead(pid, s=5.0):
    end = time.monotonic() + s
    while time.monotonic() < end and _alive(pid):
        time.sleep(0.05)
    return not _alive(pid)


def test_what_a_command_left_running_ends_with_it(tmp_path):
    pidfile = tmp_path / "pid"
    out = C.run(["/bin/bash", "-c", f"sleep 60 & echo $! > {pidfile}; echo done"], timeout=30)
    assert out.returncode == 0 and "done" in out.stdout
    assert _wait_dead(int(pidfile.read_text()))


def test_a_timeout_ends_the_whole_group(tmp_path):
    pidfile = tmp_path / "pid"
    with pytest.raises(subprocess.TimeoutExpired) as exc:
        C.run(["/bin/bash", "-c", f"sleep 60 & echo $! > {pidfile}; echo started; sleep 60"],
              timeout=1)
    assert "started" in (exc.value.stdout or "")
    assert _wait_dead(int(pidfile.read_text()))


def test_the_command_has_no_terminal_and_a_session_of_its_own():
    # Asked in Python: `ps -o sid=` is Linux procps, not macOS ps.
    out = C.run([sys.executable, "-c",
                 "import os, sys; print('tty' if sys.stdin.isatty() else 'notty', os.getsid(0))"],
                timeout=30)
    first, sid = out.stdout.split()
    assert first == "notty"
    assert int(sid) != os.getsid(0)


def test_output_exit_code_and_input_come_through():
    out = C.run("read x; echo got:$x; echo err >&2; exit 3", shell=True,
                input="hello\n", timeout=30)
    assert (out.returncode, out.stdout.strip(), out.stderr.strip()) == (3, "got:hello", "err")


def test_the_agents_shell_uses_it_where_there_is_no_cage(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_bash_isolation_argv", lambda cmd, cwd, perms: ["/bin/bash", "-c", cmd])
    perms = A.KitToolPermissions(workspace=str(tmp_path), mode="bypassPermissions")
    pidfile = tmp_path / "pid"
    res = json.loads(A._doc_executor._execute_bash(
        {"command": f"sleep 60 & echo $! > {pidfile}; echo ok"}, perms))
    assert "ok" in res.get("stdout", ""), res
    assert _wait_dead(int(pidfile.read_text()))


def test_the_terminal_agent_ends_its_background_shells_on_the_way_out():
    import inspect
    from delfin.agent import cli
    src = inspect.getsource(cli.main)
    assert "_atexit.register(_stop_own_background_shells)" in src
    assert "_stop_own_background_shells()" in src and "_lifeline.watch(_end_with_everything)" in src
