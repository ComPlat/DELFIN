"""On a login node, a command nobody was asked about runs on two CPUs.

Input: the argv the bash executors start (foreground and background), the
host kind, and whether a person approved the command in a dialog.

calc_policy decides a calculation by reading its command. A command it
cannot read -- `python3 run.py` that starts `orca` with PAL16, `make -j`,
a shell script -- passed the auto-allow list on its first word and took
whatever cores it asked for. On a shared login node that is the overload
the user ruled out (decision 2026-10-22: locally at most two cores).

After: on a login node, a command that runs without a dialog starts under
a CPU affinity of LOGIN_LOCAL_CORES CPUs (taskset where installed, else a
small launcher calling sched_setaffinity), inherited by everything it
starts. A command a person approved in the dialog runs as before. A
workstation and a Slurm allocation (whose cgroup already binds it) are
unchanged, and a host without an affinity call starts the command as
before.
"""
from __future__ import annotations

import json
import os
import sys

import pytest

from delfin.agent import api_client as A
from delfin.agent import calc_policy as C
from delfin.agent import cpu_share


# -- the parts -------------------------------------------------------------

def test_a_login_node_gives_an_unasked_run_two_cpus(monkeypatch):
    monkeypatch.setattr(C.os, "sched_getaffinity",
                        lambda pid: set(range(16)), raising=False)
    cpus = C.unasked_cpu_set("login_node")
    assert len(cpus) == C.LOGIN_LOCAL_CORES
    assert set(cpus) <= set(range(16))


@pytest.mark.parametrize("kind", ["workstation", "slurm_job"])
def test_other_hosts_are_not_bound(monkeypatch, kind):
    monkeypatch.setattr(C.os, "sched_getaffinity",
                        lambda pid: set(range(16)), raising=False)
    assert C.unasked_cpu_set(kind) is None


def test_a_host_without_the_affinity_call_is_not_bound(monkeypatch):
    monkeypatch.delattr(C.os, "sched_getaffinity", raising=False)
    assert C.unasked_cpu_set("login_node") is None


def test_a_process_already_on_two_cpus_is_left_alone(monkeypatch):
    monkeypatch.setattr(C.os, "sched_getaffinity", lambda pid: {3, 7},
                        raising=False)
    assert C.unasked_cpu_set("login_node") is None


def test_taskset_is_used_where_installed(monkeypatch):
    monkeypatch.setattr(C.shutil, "which",
                        lambda name: "/usr/bin/taskset"
                        if name == "taskset" else None)
    assert C.bound_argv(["bash", "-c", "x"], [2, 5]) == [
        "/usr/bin/taskset", "-c", "2,5", "bash", "-c", "x"]


def test_the_launcher_stands_in_for_taskset(monkeypatch):
    monkeypatch.setattr(C.shutil, "which", lambda name: None)
    assert C.bound_argv(["bash", "-c", "x"], [2, 5]) == [
        sys.executable, "-m", "delfin.agent.cpu_share", "2,5",
        "bash", "-c", "x"]


def test_the_launcher_binds_itself_then_runs_the_command(monkeypatch):
    seen = {}
    monkeypatch.setattr(cpu_share.os, "sched_setaffinity",
                        lambda pid, cpus: seen.update(cpus=set(cpus)),
                        raising=False)
    monkeypatch.setattr(cpu_share.os, "execvp",
                        lambda prog, argv: seen.update(argv=argv))
    cpu_share.main(["2,5", "bash", "-c", "x"])
    assert seen == {"cpus": {2, 5}, "argv": ["bash", "-c", "x"]}


def test_the_launcher_runs_the_command_where_it_cannot_bind(monkeypatch):
    seen = {}
    monkeypatch.delattr(cpu_share.os, "sched_setaffinity", raising=False)
    monkeypatch.setattr(cpu_share.os, "execvp",
                        lambda prog, argv: seen.update(argv=argv))
    cpu_share.main(["2,5", "bash", "-c", "x"])
    assert seen == {"argv": ["bash", "-c", "x"]}


# -- through the executor --------------------------------------------------

@pytest.fixture
def started(tmp_path, monkeypatch):
    """The argv the foreground executor hands to contained_run."""
    from delfin.agent import contained_run
    monkeypatch.setattr(A, "_BASH_ISOLATION_OVERRIDE", "off")
    monkeypatch.setattr(C, "unasked_cpu_set", lambda kind=None: [0, 1])
    monkeypatch.setattr(C.shutil, "which", lambda name: None)
    seen: list = []

    def fake_run(argv, **kw):
        seen.append(list(argv))
        import subprocess
        return subprocess.CompletedProcess(argv, 0, "", "")
    monkeypatch.setattr(contained_run, "run", fake_run)
    return seen


def _perms(ws, confirm=None):
    perms = A.KitToolPermissions(workspace=ws, mode="default",
                                 agent_role="solo_agent")
    perms.confirm_callback = confirm
    return perms


def test_an_auto_allowed_command_starts_bound(tmp_path, started):
    json.loads(A._doc_executor.execute("bash", {"command": "ls"},
                                       _perms(tmp_path)))
    assert started and started[-1][:4] == [
        sys.executable, "-m", "delfin.agent.cpu_share", "0,1"]


def test_an_approved_command_starts_unbound(tmp_path, started):
    json.loads(A._doc_executor.execute(
        "bash", {"command": "gcc --version"},
        _perms(tmp_path, confirm=lambda *a, **k: True)))
    assert started and started[-1][:3] != [
        sys.executable, "-m", "delfin.agent.cpu_share"]


def test_a_background_job_starts_bound(tmp_path, monkeypatch, started):
    from delfin.agent import bash_jobs
    seen = {}

    def fake_start(self, **kw):
        seen.update(kw)
        raise ValueError("stopped by the test")
    monkeypatch.setattr(bash_jobs._Registry, "start", fake_start)
    A._doc_executor.execute("bash_background", {"command": "ls"},
                            _perms(tmp_path))
    assert seen["argv"][:4] == [
        sys.executable, "-m", "delfin.agent.cpu_share", "0,1"]
