"""SL1 phase 2: auto-PAL honours the local-run cap on a non-cluster host.

Reviewer point 4 / inventory item 4, now a measured defect, not a
hypothesis: auto_configure_resources on a machine that is neither in a
scheduler job nor detected as a cluster sets PAL = os.cpu_count() with
no cap (cluster_utils.py:194-197). On a workstation with 192 CPUs the
probe returned PAL=192 while DELFIN's own calc_policy documents
LOGIN_LOCAL_CORES = 2 as the cap for a local run (calc_policy.py:61).
The consumer path is real: cli.py:1817 -> ctx.PAL (cli.py:1956) feeds
local execution through the scheduling manager.

The fix: when nothing indicates a scheduler allocation, the auto-set
PAL is capped at calc_policy.LOGIN_LOCAL_CORES. An explicit PAL in the
config is never touched (only `not config.get('PAL')` enters), and a
SLURM/PBS/LSF job keeps its allocation's full size.
"""
from __future__ import annotations

import os

import pytest

from delfin.cluster_utils import _detect_system_resources, auto_configure_resources


def _clear_scheduler_env(monkeypatch):
    for var in ("SLURM_JOB_ID", "SLURM_CPUS_PER_TASK", "SLURM_NPROCS",
                "PBS_JOBID", "LSB_JOBID", "PBS_NUM_PPN", "NCPUS",
                "LSB_DJOB_NUMPROC"):
        monkeypatch.delenv(var, raising=False)


def _big_machine(monkeypatch):
    _clear_scheduler_env(monkeypatch)
    monkeypatch.setattr("delfin.cluster_utils.psutil", None)
    monkeypatch.setattr(os, "cpu_count", lambda: 192, raising=False)
    import os as _os
    monkeypatch.setattr(_os, "cpu_count", lambda: 192)


def test_auto_pal_on_a_big_local_host_is_capped(monkeypatch):
    _big_machine(monkeypatch)
    info = _detect_system_resources()
    assert info["cpus_available"] == 192
    cfg = auto_configure_resources({"maxcore": 1000})
    from delfin.agent.calc_policy import LOGIN_LOCAL_CORES
    assert cfg["PAL"] == LOGIN_LOCAL_CORES


def test_explicit_pal_is_never_touched(monkeypatch):
    _big_machine(monkeypatch)
    cfg = auto_configure_resources({"PAL": 4})
    assert cfg["PAL"] == 4


def test_slurm_job_allocation_is_not_capped(monkeypatch):
    _big_machine(monkeypatch)
    monkeypatch.setenv("SLURM_JOB_ID", "12345")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "64")
    cfg = auto_configure_resources({})
    assert cfg["PAL"] == 64
