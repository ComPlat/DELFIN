"""A command the agent runs starts with a thread cap unless one is set.

Input: the environment every agent child is built from
(``_scrubbed_bash_env``: shell tool, background job, test run).
Before: no OpenMP/BLAS thread variable was set, so an OpenMP program such
as xtb took every core the host has. Measured on a 384-core host: a
37-point GFN2-xTB dihedral scan of butane ran at 29 264 % CPU (44 052
CPU-s in 151 s wall). On a cluster login node that is every user's cores.
After: OMP_NUM_THREADS, MKL_NUM_THREADS and OPENBLAS_NUM_THREADS default
to the Slurm allocation's cores when the agent runs inside one
(SLURM_CPUS_PER_TASK), else to min(4, the unasked local limit) --
calc_policy.default_threads, so the gate judges the number the shell
sets (half the available CPUs on a workstation, 2 on a login node); a
value the user's environment already sets is kept. A Slurm job inherits
the variable at submission unless its script sets its own. The host kind
and CPU count are set per test, so the expectation holds on any machine.
"""
from __future__ import annotations

import pytest

from delfin.agent import calc_policy as C
from delfin.agent.api_client import _scrubbed_bash_env

CAPPED = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")


def _clear(monkeypatch):
    for name in CAPPED + ("SLURM_CPUS_PER_TASK",):
        monkeypatch.delenv(name, raising=False)


def _host(monkeypatch, kind, cpus):
    monkeypatch.setattr(C, "host_kind", lambda env=None: kind)
    monkeypatch.setattr(C, "_available_cpus", lambda: cpus)


def _a_few(kind, cpus):
    """The rule: at most four, never above the unasked local limit."""
    limit = C.LOGIN_LOCAL_CORES if kind == "login_node" else max(1, cpus // 2)
    return str(max(1, min(4, limit)))


@pytest.mark.parametrize("cpus", [1, 2, 4, 8, 384])
@pytest.mark.parametrize("kind", ["workstation", "login_node"])
def test_an_unset_cap_defaults_to_a_few_cores(monkeypatch, kind, cpus):
    _clear(monkeypatch)
    _host(monkeypatch, kind, cpus)
    env = _scrubbed_bash_env()
    want = _a_few(kind, cpus)
    assert all(env.get(name) == want for name in CAPPED), \
        {n: env.get(n) for n in CAPPED}


def test_inside_an_allocation_the_cap_is_the_allocation(monkeypatch):
    _clear(monkeypatch)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
    env = _scrubbed_bash_env()
    assert all(env.get(name) == "16" for name in CAPPED)


def test_a_cap_the_user_set_is_kept(monkeypatch):
    _clear(monkeypatch)
    monkeypatch.setenv("OMP_NUM_THREADS", "32")
    env = _scrubbed_bash_env()
    assert env["OMP_NUM_THREADS"] == "32"


@pytest.mark.parametrize("cpus", [1, 4, 8, 64])
def test_a_malformed_allocation_falls_back_to_the_default(monkeypatch, cpus):
    _clear(monkeypatch)
    _host(monkeypatch, "workstation", cpus)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "all")
    env = _scrubbed_bash_env()
    assert env["OMP_NUM_THREADS"] == _a_few("workstation", cpus)
