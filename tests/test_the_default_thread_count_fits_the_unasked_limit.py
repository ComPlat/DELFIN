"""A command that names no thread count gets one the gate admits.

Input: calc_policy.default_threads / local_core_limit, the agent's shell
environment (api_client._scrubbed_bash_env) and the bash gate through the
real executor, on a host kind and CPU count set per test.

Before: on a workstation the default was ``min(4, os.cpu_count())`` and the
unasked limit ``available CPUs // 2``. On a 4-core machine every plain
``xtb mol.xyz --opt`` therefore asked for 4 cores against a limit of 2 and
was refused (field test on a 4-core workstation), while the shell
environment set OMP_NUM_THREADS=4 for the same command.

Pinned: on every host kind and CPU count the default is at most the
unasked local limit; on a 4-core workstation a plain xtb run is admitted
with 2 threads and its environment says 2; on a login node 2; inside a
Slurm job the allocation. Controls: an explicit ``-P`` above the limit and
the user's own OMP_NUM_THREADS above it are still refused.
"""
from __future__ import annotations

import json
import os

import pytest

from delfin.agent import api_client as A
from delfin.agent import calc_policy as C

H2 = "2\n\nH 0 0 0\nH 0 0 0.74\n"
_THREAD_VARS = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "SLURM_JOB_ID", "SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE")


@pytest.fixture(autouse=True)
def _no_thread_env(monkeypatch):
    for name in _THREAD_VARS:
        monkeypatch.delenv(name, raising=False)


def _host(monkeypatch, kind, cpus):
    monkeypatch.setattr(C, "host_kind", lambda env=None: kind)
    monkeypatch.setattr(C, "_available_cpus", lambda: cpus)
    # The machine's own count must not decide it either.
    monkeypatch.setattr(os, "cpu_count", lambda: cpus)


def _perms(ws):
    perms = A.KitToolPermissions(workspace=ws, mode="default",
                                 agent_role="solo_agent")
    perms.confirm_callback = None
    return perms


def _gate(ws, cmd, notes=None):
    return A._doc_executor._gate_bash_call("bash", {"command": cmd},
                                           _perms(ws), [], notes)


@pytest.mark.parametrize("cpus", [1, 2, 3, 4, 6, 8, 16, 64, 384])
@pytest.mark.parametrize("kind", ["workstation", "login_node", "slurm_job"])
def test_the_default_never_exceeds_the_unasked_limit(monkeypatch, kind, cpus):
    _host(monkeypatch, kind, cpus)
    assert 1 <= C.default_threads(kind, {}) <= C.local_core_limit(kind)


@pytest.mark.parametrize("cpus", [1, 2, 4, 8, 64])
@pytest.mark.parametrize("kind", ["workstation", "login_node"])
def test_a_plain_xtb_run_is_admitted_on_every_host(tmp_path, monkeypatch,
                                                   kind, cpus):
    _host(monkeypatch, kind, cpus)
    (tmp_path / "mol.xyz").write_text(H2)
    notes: list = []
    assert _gate(tmp_path, "xtb mol.xyz --opt", notes) is None
    assert notes and f"{C.default_threads(kind, {})} core(s)" in notes[0]


def test_on_a_four_core_workstation_xtb_gets_two_threads(tmp_path,
                                                         monkeypatch):
    _host(monkeypatch, "workstation", 4)
    (tmp_path / "mol.xyz").write_text(H2)
    notes: list = []
    assert _gate(tmp_path, "xtb mol.xyz --opt > out 2>&1", notes) is None
    assert "xtb mol.xyz, 2 core(s) on this machine" in notes[0]
    assert A._scrubbed_bash_env()["OMP_NUM_THREADS"] == "2"


def test_the_environment_the_command_gets_says_two(tmp_path, monkeypatch):
    """Through the real executor: the shell the agent's command runs in."""
    _host(monkeypatch, "workstation", 4)
    out = json.loads(A._doc_executor.execute(
        "bash", {"command": "printenv OMP_NUM_THREADS"}, _perms(tmp_path)))
    assert str(out.get("stdout", "")).strip() == "2", out


def test_on_a_login_node_it_stays_two(tmp_path, monkeypatch):
    _host(monkeypatch, "login_node", 64)
    (tmp_path / "mol.xyz").write_text(H2)
    notes: list = []
    assert _gate(tmp_path, "xtb mol.xyz --opt", notes) is None
    assert "2 core(s) on this login node" in notes[0]
    assert A._scrubbed_bash_env()["OMP_NUM_THREADS"] == "2"


@pytest.mark.parametrize("var", ["SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"])
def test_in_a_slurm_job_it_gets_the_allocation(tmp_path, monkeypatch, var):
    """The host kind and the CPU count from the job's own variables, as the
    gate and the shell read them (no patched host kind here)."""
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setenv(var, "6")
    monkeypatch.setattr(os, "cpu_count", lambda: 64)
    (tmp_path / "mol.xyz").write_text(H2)
    assert C.host_kind() == "slurm_job"
    notes: list = []
    assert _gate(tmp_path, "xtb mol.xyz --opt", notes) is None
    assert "6 core(s) in this Slurm allocation" in notes[0]
    assert A._scrubbed_bash_env()["OMP_NUM_THREADS"] == "6"


# -- controls: the limit itself still holds --------------------------------

def test_an_explicit_width_above_the_limit_is_refused(tmp_path, monkeypatch):
    _host(monkeypatch, "workstation", 4)
    (tmp_path / "mol.xyz").write_text(H2)
    out = _gate(tmp_path, "xtb mol.xyz --opt -P 4")
    assert out is not None and "asks for 4 cores" in out


def test_the_users_own_thread_setting_is_kept_and_judged(tmp_path,
                                                         monkeypatch):
    _host(monkeypatch, "workstation", 4)
    monkeypatch.setenv("OMP_NUM_THREADS", "8")
    (tmp_path / "mol.xyz").write_text(H2)
    assert A._scrubbed_bash_env()["OMP_NUM_THREADS"] == "8"
    out = _gate(tmp_path, "xtb mol.xyz --opt")
    assert out is not None and "asks for 8 cores" in out
