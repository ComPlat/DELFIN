"""A command the agent runs starts with a thread cap unless one is set.

Input: the environment every agent child is built from
(``_scrubbed_bash_env``: shell tool, background job, test run).
Before: no OpenMP/BLAS thread variable was set, so an OpenMP program such
as xtb took every core the host has. Measured on a 384-core host: a
37-point GFN2-xTB dihedral scan of butane ran at 29 264 % CPU (44 052
CPU-s in 151 s wall). On a cluster login node that is every user's cores.
After: OMP_NUM_THREADS, MKL_NUM_THREADS and OPENBLAS_NUM_THREADS default
to the Slurm allocation's cores when the agent runs inside one
(SLURM_CPUS_PER_TASK), else to min(4, cores); a value the user's
environment already sets is kept. A Slurm job inherits the variable at
submission unless its script sets its own.
"""
from __future__ import annotations

import os

from delfin.agent.api_client import _scrubbed_bash_env

CAPPED = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")


def _clear(monkeypatch):
    for name in CAPPED + ("SLURM_CPUS_PER_TASK",):
        monkeypatch.delenv(name, raising=False)


def test_an_unset_cap_defaults_to_a_few_cores(monkeypatch):
    _clear(monkeypatch)
    env = _scrubbed_bash_env()
    want = str(max(1, min(4, os.cpu_count() or 1)))
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


def test_a_malformed_allocation_falls_back_to_the_default(monkeypatch):
    _clear(monkeypatch)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "all")
    env = _scrubbed_bash_env()
    assert env["OMP_NUM_THREADS"] == str(max(1, min(4, os.cpu_count() or 1)))
