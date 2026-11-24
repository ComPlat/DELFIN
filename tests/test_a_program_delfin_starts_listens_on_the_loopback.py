"""A program DELFIN starts listens on the loopback only.

Measured 2026-11-24: a two-process ORCA single point listened on
0.0.0.0:1024-1025 and [::]:1024-1025 (OpenMPI's TCP transport) and mpirun
on four more wildcard ports, reachable from the network of a host without
a firewall. An agent's command is held by its cage; an ORCA that DELFIN
itself starts -- a workflow, the ops server's run tools, the dashboard's
input check and saddle search -- was held by nothing.
``delfin.inbound_guard.guarded`` starts such a program under the socket
guard's inbound rule. These tests pin when it does and when it leaves the
argv alone (every machine), and, where the kernel has seccomp user
notification, that the program then binds only the loopback, still stops
when the caller stops it, and ends with its supervisor.
"""
from __future__ import annotations

import os
import signal
import subprocess
import sys
import time

import pytest

from conftest import child_env

from delfin import inbound_guard as G
from delfin.agent import socket_guard as SG



def _runs_as_before_without_the_guard():
    """On a kernel without seccomp user notification the guard does not run:
    the probe must say so, and a program then starts exactly as it would
    have before -- checked instead of skipping."""
    saved = G._usable
    G._usable = None
    try:
        assert G.usable() is False
        assert G.guarded(["orca", "job.inp"], {}) == ["orca", "job.inp"]
    finally:
        G._usable = saved


@pytest.fixture
def probe_said(monkeypatch):
    """The probe's answer, set by the test instead of run."""
    def set_to(value):
        monkeypatch.setattr(G, "_usable", value)
    return set_to


# -- when the argv is wrapped, on every machine ---------------------------------

def test_a_local_program_runs_under_the_loopback_rule(probe_said):
    probe_said(True)
    argv = G.guarded(["orca", "job.inp"], {})
    assert argv[:2] == [sys.executable, "-I"] and argv[2].endswith("landlock_exec.py")
    tail = argv[argv.index("--inbound"):]
    assert tail == ["--inbound", "loopback", "--", "orca", "job.inp"]
    assert ["--net", "open"] == argv[argv.index("--net"):argv.index("--net") + 2]
    assert ["--fs", "0"] == argv[argv.index("--fs"):argv.index("--fs") + 2]


def test_a_multi_node_allocation_keeps_its_binds(probe_said):
    probe_said(True)
    env = {"SLURM_JOB_ID": "11", "SLURM_JOB_NUM_NODES": "3"}
    assert G.guarded(["orca", "job.inp"], env)[-4:-2] == ["any", "--"]
    one = {"SLURM_JOB_ID": "11", "SLURM_JOB_NUM_NODES": "1"}
    assert G.guarded(["orca", "job.inp"], one)[-4:-2] == ["loopback", "--"]


@pytest.mark.parametrize("value", ["held", "off", "OFF", "0", "no"])
def test_a_held_or_switched_off_program_is_left_as_it_is(probe_said, value):
    probe_said(True)
    assert G.guarded(["orca", "job.inp"], {"DELFIN_INBOUND_GUARD": value}) == ["orca", "job.inp"]


def test_where_the_guard_cannot_run_the_program_runs_as_before(probe_said):
    probe_said(False)
    assert G.guarded(["orca", "job.inp"], {}) == ["orca", "job.inp"]


def test_a_probe_that_could_not_run_is_asked_again(monkeypatch):
    """An exception is no answer: it costs that one start its guard, and the
    next start probes again instead of running unguarded for good."""
    monkeypatch.setattr(G, "_usable", None)
    monkeypatch.setattr(G.sys, "platform", "linux")

    def broken(*a, **k):
        raise RuntimeError("a fake Popen")

    monkeypatch.setattr(G.subprocess, "run", broken)
    assert G.usable() is False and G._usable is None

    class Done:
        returncode = 0

    monkeypatch.setattr(G.subprocess, "run", lambda *a, **k: Done())
    assert G.usable() is True and G._usable is True


def test_a_probe_without_the_marker_is_no_guard(monkeypatch):
    """The helper's fallback (no user notification) runs the probe without
    the marker and refuses Unix sockets outright: that is not usable."""
    monkeypatch.setattr(G, "_usable", None)
    monkeypatch.setattr(G.sys, "platform", "linux")

    class Fell:
        returncode = 3

    monkeypatch.setattr(G.subprocess, "run", lambda *a, **k: Fell())
    assert G.usable() is False


def test_the_allocation_rule_is_one_rule():
    from delfin.agent import api_client as A
    for env in ({}, {"SLURM_JOB_ID": "1", "SLURM_NNODES": "2"},
                {"SLURM_JOB_ID": "1", "SLURM_JOB_NUM_NODES": "x"}):
        assert A._in_multi_node_allocation(env) is G.multi_node_allocation(env)


# -- live, under the real guard ----------------------------------------------------

_SERVE = r'''
import os, signal, socket, sys, time
s = socket.socket(); s.bind(("0.0.0.0", 0))
if s.getsockname()[0] != "127.0.0.1":
    print("WILDCARD", flush=True); sys.exit(1)
s.listen(1)
print(s.getsockname()[0], os.environ.get("DELFIN_INBOUND_GUARD"), flush=True)
def bye(*_):
    open(sys.argv[1], "w").write("terminated"); sys.exit(0)
signal.signal(signal.SIGTERM, bye)
time.sleep(60)
'''


def _start(tmp_path, argv):
    return subprocess.Popen(argv, cwd=str(tmp_path), env=child_env(tmp_path),
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            start_new_session=True)


def test_the_program_binds_the_loopback_and_stops_when_its_caller_stops_it(tmp_path):
    if not SG.available():
        return _runs_as_before_without_the_guard()
    monkey_usable = G._usable
    G._usable = None
    try:
        assert G.usable() is True
    finally:
        G._usable = monkey_usable
    note = tmp_path / "note"
    proc = _start(tmp_path, G._wrap([sys.executable, "-c", _SERVE, str(note)], "loopback"))
    try:
        first = proc.stdout.readline().split()
        assert first == ["127.0.0.1", "held"], (first, proc.stderr.read() if proc.poll() else "")
        proc.terminate()                     # the supervisor, as a caller would
        proc.wait(timeout=30)
        assert note.read_text() == "terminated"
    finally:
        if proc.poll() is None:
            os.killpg(proc.pid, signal.SIGKILL)


def test_the_program_ends_with_its_supervisor(tmp_path):
    if not SG.available():
        return _runs_as_before_without_the_guard()
    proc = _start(tmp_path, G._wrap([sys.executable, "-c",
                                     "import os, time; print(os.getpid(), flush=True); time.sleep(60)"],
                                    "loopback"))
    try:
        child = int(proc.stdout.readline())
        os.kill(proc.pid, signal.SIGKILL)    # the supervisor dies hard
        proc.wait(timeout=30)
        deadline = time.time() + 10
        while time.time() < deadline and os.path.exists(f"/proc/{child}"):
            try:
                with open(f"/proc/{child}/stat") as fh:
                    if fh.read().rsplit(")", 1)[1].split()[0] == "Z":
                        break
            except OSError:
                break
            time.sleep(0.1)
        alive = os.path.exists(f"/proc/{child}")
        if alive:
            with open(f"/proc/{child}/stat") as fh:
                alive = fh.read().rsplit(")", 1)[1].split()[0] != "Z"
        assert not alive, "the program outlived its supervisor"
    finally:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except OSError:
            pass


def test_a_program_already_held_is_not_wrapped_again(tmp_path):
    """A DELFIN under the guard (an agent's command) starts its programs
    without a second guard: they are held by the first."""
    if not SG.available():
        return _runs_as_before_without_the_guard()
    inner = ("import json, sys; from delfin import inbound_guard as G; "
             "print(json.dumps(G.guarded(['orca', 'job.inp'])))")
    env = dict(child_env(tmp_path), PYTHONPATH=os.pathsep.join(
        p for p in (str(G._HELPER.parents[2]), os.environ.get("PYTHONPATH", "")) if p))
    out = subprocess.run(G._wrap([sys.executable, "-c", inner], "loopback"), cwd=str(tmp_path),
                         env=env, capture_output=True, text=True, timeout=120)
    assert out.stdout.strip() == '["orca", "job.inp"]', out.stderr[-1500:]


_LAUNCHER = r'''#!PYTHON
"""A stand-in MPI launcher: starts -np N copies of its command, as mpirun
starts ranks, and fails if any rank fails."""
import subprocess, sys
args = sys.argv[1:]
n = int(args[args.index("-np") + 1]); cmd = args[args.index("-np") + 2:]
ranks = [subprocess.Popen(cmd) for _ in range(n)]
sys.exit(max(r.wait() for r in ranks))
'''


def test_ranks_a_launcher_starts_listen_on_the_loopback(tmp_path):
    """Two ranks of an MPI-style launcher (a stand-in on PATH; the real
    OpenMPI was measured by hand, see the module notes), each binding the
    wildcard as a TCP transport does: under the rule both land on the
    loopback."""
    if not SG.available():
        return _runs_as_before_without_the_guard()
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    launcher = bin_dir / "mpirun"
    launcher.write_text(_LAUNCHER.replace("PYTHON", sys.executable))
    launcher.chmod(0o755)
    rank = ("import socket; s = socket.socket(); s.bind(('0.0.0.0', 0)); "
            "print(s.getsockname()[0], flush=True); s.close()")
    env = dict(child_env(tmp_path))
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    out = subprocess.run(G._wrap(["mpirun", "-np", "2", sys.executable, "-c", rank], "loopback"),
                         cwd=str(tmp_path), env=env, capture_output=True, text=True,
                         timeout=120)
    assert out.returncode == 0, out.stderr[-1500:]
    assert out.stdout.split() == ["127.0.0.1", "127.0.0.1"], out.stdout
