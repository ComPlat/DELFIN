"""An agent process cannot be read by the commands it runs, on any host.

The model's key lives in the agent's memory and often its environment.
Without the bubblewrap PID namespace, a command could read
/proc/<agent>/environ, and where ptrace_scope is 0 (RHEL, Rocky, many HPC
nodes) attach a debugger to the agent (audit 2026-09-16). A protected
process must also stay findable by the emergency stop, which read its
environment to recognise it.
"""
import json
import os
import subprocess
import sys
import time

import pytest

pytestmark = pytest.mark.skipif(not sys.platform.startswith("linux"),
                                reason="PR_SET_DUMPABLE is Linux")

_PROBE = """\
import ctypes, sys
pid = int(sys.argv[1])
try:
    open(f"/proc/{pid}/environ", "rb").read()
    env = "READ"
except PermissionError:
    env = "denied"
libc = ctypes.CDLL(None, use_errno=True)
rc = libc.ptrace(16, pid, None, None)          # PTRACE_ATTACH
attach = "ATTACHED" if rc == 0 else "denied"
if rc == 0:
    libc.ptrace(17, pid, None, None)           # PTRACE_DETACH
own = "own-ok" if open("/proc/self/environ", "rb").read() else "own-empty"
print(env, attach, own)
"""

_AGENT = """\
import os, subprocess, sys, time
from delfin.agent import process_guard
if os.environ.get("AGENT_UNPROTECTED") != "1":
    assert process_guard.protect("test agent")
probe = subprocess.run([sys.executable, sys.argv[1], str(os.getpid())],
                       capture_output=True, text=True)
print(probe.stdout.strip() or probe.stderr.strip(), flush=True)
if len(sys.argv) > 2:
    time.sleep(float(sys.argv[2]))
"""


def _scripts(tmp_path):
    probe, agent = tmp_path / "probe.py", tmp_path / "agent.py"
    probe.write_text(_PROBE)
    agent.write_text(_AGENT)
    return str(agent), str(probe)


def _env(tmp_path):
    env = dict(os.environ)
    env.pop("DELFIN_PROCESS_GUARD", None)
    env["HOME"] = str(tmp_path)
    env["DELFIN_TEST_SECRET"] = "in-the-environment"
    return env


def test_a_command_cannot_read_the_agents_environment_or_attach(tmp_path):
    agent, probe = _scripts(tmp_path)
    out = subprocess.run([sys.executable, agent, probe], env=_env(tmp_path),
                         capture_output=True, text=True, timeout=60)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip().splitlines()[-1] == "denied denied own-ok"


def test_without_the_guard_a_command_reads_it(tmp_path):
    """The control on this very host: the same probe against an agent that
    did not protect itself reads its environment."""
    agent, probe = _scripts(tmp_path)
    env = _env(tmp_path)
    env["AGENT_UNPROTECTED"] = "1"
    out = subprocess.run([sys.executable, agent, probe], env=env,
                         capture_output=True, text=True, timeout=60)
    assert out.stdout.strip().splitlines()[-1].startswith("READ ")


def test_the_emergency_stop_still_finds_a_protected_process(tmp_path, monkeypatch):
    agent, probe = _scripts(tmp_path)
    child = subprocess.Popen([sys.executable, agent, probe, "30"], env=_env(tmp_path),
                             stdout=subprocess.PIPE, text=True, start_new_session=True)
    try:
        child.stdout.readline()
        from delfin.agent import process_guard, stop_all
        monkeypatch.setattr(process_guard, "_DIR", tmp_path / ".delfin" / "agent_processes")
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline and not process_guard.registered_here():
            time.sleep(0.1)
        assert [r["pid"] for r in process_guard.registered_here()] == [child.pid]
        assert child.pid in stop_all._own_agent_processes(time.time() + 60)
    finally:
        child.kill()
        child.wait()
    assert process_guard.registered_here() == []      # gone process, record dropped


def test_off_leaves_a_process_debuggable(monkeypatch):
    from delfin.agent import process_guard
    monkeypatch.setenv("DELFIN_PROCESS_GUARD", "off")
    assert process_guard.protect("anything") is False
    assert process_guard.is_protected() is False
