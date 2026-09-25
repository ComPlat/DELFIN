"""A background job outlives the tool call that started it -- not the process.

Seen while measuring on 2026-09-21: an agent started `sleep 311` as a
background job with a 24-hour timeout, and 4.4 seconds later the job was
gone, killed by SIGKILL. Nothing in DELFIN had sent it: the job's own
timeout path sends SIGTERM first.

The cage runs a job under bubblewrap's --die-with-parent, which is
PR_SET_PDEATHSIG -- and Linux ties that to the THREAD that forked, not to
the process. Tool calls run on threads that end when the call returns, so
every caged background job was killed as soon as its call came back.
Reproduced without a model: the same cage and registry, started from a
thread that then ended, dead after 0.1 s; started from the main thread,
alive.

The cage is right to end a job with the agent -- that is the rule this
suite also holds (a session ends with everything it started). What was
wrong is WHICH thread counts as the parent. A child that must end with the
process is now forked from a thread that lives as long as the process.
The same goes for an MCP server under isolation: it is started on first
use, from a tool call's thread, and died with it.
"""

from __future__ import annotations

import os
import subprocess
import sys
import textwrap
import threading
import time
from pathlib import Path

import pytest

from conftest import child_env

from delfin.agent import api_client as A
from delfin.agent import bash_jobs as BJ
from delfin.agent import mcp_client as M

pytestmark = pytest.mark.skipif(
    not sys.platform.startswith("linux"),
    reason="PR_SET_PDEATHSIG is Linux; elsewhere nothing ties a child to a thread")


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


def _dies_with_its_parent_thread(marker: Path, *, then: str) -> list[str]:
    """What --die-with-parent does, without needing bubblewrap: the child
    asks for SIGKILL when the thread that forked it ends, says so, then
    does *then*."""
    code = textwrap.dedent(f"""
        import ctypes, signal, sys, time
        ctypes.CDLL(None, use_errno=True).prctl(1, signal.SIGKILL)
        open({str(marker)!r}, "w").write("armed")
        {then}
    """)
    return [sys.executable, "-c", code]


def _from_a_thread_that_ends(start, marker: Path):
    """Run *start* on a thread that waits until the child is armed, then ends."""
    box: dict = {}

    def _call():
        box["out"] = start()
        end = time.monotonic() + 10
        while time.monotonic() < end and not marker.exists():
            time.sleep(0.02)

    t = threading.Thread(target=_call, name="a-tool-call")
    t.start()
    t.join()
    assert marker.exists(), "the child never armed its death signal"
    return box["out"]


@pytest.fixture
def jobs():
    started: list = []
    yield started
    for job in started:
        try:
            BJ.get_registry().kill(job.job_id)
        except Exception:
            pass


def test_a_job_started_from_a_call_is_still_running_after_it(tmp_path, jobs):
    marker = tmp_path / "armed"
    argv = _dies_with_its_parent_thread(marker, then="time.sleep(60)")
    job = _from_a_thread_that_ends(
        lambda: BJ.get_registry().start(
            command="sleep 60", cwd=str(tmp_path), workspace=str(tmp_path),
            timeout_s=120, argv=argv),
        marker)
    jobs.append(job)
    time.sleep(1.0)
    assert job.proc.poll() is None, (
        f"the job died with the call that started it (exit {job.proc.poll()})")


def test_a_caged_job_is_still_running_after_its_call(tmp_path, jobs):
    if not A._bwrap_functional():
        pytest.skip("no working bubblewrap here")
    perms = A.KitToolPermissions(workspace=str(tmp_path))
    argv = A._bash_isolation_argv("sleep 60", tmp_path, perms)
    box: dict = {}

    def _call():
        box["job"] = BJ.get_registry().start(
            command="sleep 60", cwd=str(tmp_path), workspace=str(tmp_path),
            timeout_s=120, argv=argv)
        time.sleep(0.5)          # the cage is up before the call returns

    t = threading.Thread(target=_call, name="a-tool-call")
    t.start()
    t.join()
    job = box["job"]
    jobs.append(job)
    time.sleep(1.0)
    assert job.proc.poll() is None, (
        f"the caged job died with its call (exit {job.proc.poll()})")


def test_the_job_still_ends_with_its_process(tmp_path):
    """The other half of the rule: the job must not outlive the agent."""
    marker = tmp_path / "armed"
    pidfile = tmp_path / "pid"
    child = _dies_with_its_parent_thread(marker, then="time.sleep(60)")
    script = textwrap.dedent(f"""
        import os, sys, threading, time
        from delfin.agent import bash_jobs as BJ
        box = {{}}
        def call():
            box["job"] = BJ.get_registry().start(
                command="sleep 60", cwd={str(tmp_path)!r},
                workspace={str(tmp_path)!r}, timeout_s=120, argv={child!r})
            while not os.path.exists({str(marker)!r}):
                time.sleep(0.02)
        t = threading.Thread(target=call); t.start(); t.join()
        open({str(pidfile)!r}, "w").write(str(box["job"].proc.pid))
        os._exit(0)
    """)
    subprocess.run([sys.executable, "-c", script], timeout=60, check=True,
                   env={**child_env(tmp_path), "PYTHONPATH": os.getcwd()})
    pid = int(pidfile.read_text())
    end = time.monotonic() + 5
    while time.monotonic() < end and _alive(pid):
        time.sleep(0.05)
    assert not _alive(pid), "the job outlived the process that started it"


def test_an_mcp_server_started_on_first_use_outlives_that_call(tmp_path):
    marker = tmp_path / "armed"
    argv = _dies_with_its_parent_thread(
        marker, then="sys.stdin.read()")
    server = M.MCPServer(name="probe", command=argv[0], args=argv[1:])
    try:
        _from_a_thread_that_ends(server.start, marker)
        time.sleep(1.0)
        assert server.proc is not None and server.proc.poll() is None, (
            "the server died with the call that started it")
    finally:
        try:
            server.stop()
        except Exception:
            pass
