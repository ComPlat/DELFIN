"""The background-job cap must come from the machine, not from a memory.

The cap was the constant 4. It was not measured on the machine the field
reports come from: every timing and capacity number in this project was
taken on the workstation, while the agent runs on a shared cluster node.
Four is not wrong there so much as unrelated to it -- on a two-core VM it
oversubscribes, and inside a 64-core allocation it wastes most of what
the user was granted.

The number that matters is what THIS process may use, and the system
already states it three ways: a SLURM allocation says how many CPUs the
job was given, a cgroup quota bounds a container, and the scheduler
affinity mask reflects both plus any `taskset`. Reading them costs
nothing and is correct on every host, which a re-measurement of one host
would not be.

An explicit setting still wins -- the user knows things the machine does
not, such as sharing a node with a colleague's calculation.
"""

from __future__ import annotations

import pytest

from delfin.agent import bash_jobs as B


@pytest.fixture(autouse=True)
def no_user_setting(monkeypatch):
    """Isolate the derivation from whatever this installation configured."""
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda *a, **kw: {})
    for var in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"):
        monkeypatch.delenv(var, raising=False)


# ---------------------------------------------------------------------------
# Where the number comes from
# ---------------------------------------------------------------------------

def test_a_slurm_allocation_decides(monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
    monkeypatch.setattr(B, "_affinity_cpus", lambda: 128)
    assert B._available_cpus() == 16


def test_the_node_wide_variable_is_used_when_the_task_one_is_absent(monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "8")
    assert B._available_cpus() == 8


def test_without_slurm_the_affinity_mask_decides(monkeypatch):
    """It reflects cgroup cpusets and taskset, which nproc does not."""
    monkeypatch.setattr(B, "_affinity_cpus", lambda: 6)
    assert B._available_cpus() == 6


def test_nonsense_in_the_environment_is_ignored(monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "eight")
    monkeypatch.setattr(B, "_affinity_cpus", lambda: 4)
    assert B._available_cpus() == 4


# ---------------------------------------------------------------------------
# What is made of it
# ---------------------------------------------------------------------------

def test_half_the_cpus_with_a_floor_of_one(monkeypatch):
    monkeypatch.setattr(B, "_available_cpus", lambda: 1)
    assert B._max_background_jobs() == 1
    monkeypatch.setattr(B, "_available_cpus", lambda: 2)
    assert B._max_background_jobs() == 1


def test_a_big_allocation_does_not_become_a_fork_bomb(monkeypatch):
    """Background jobs are a convenience, not the way to use a node --
    that is what the batch scheduler is for."""
    monkeypatch.setattr(B, "_available_cpus", lambda: 256)
    assert B._max_background_jobs() <= B._MAX_BG_JOBS_CEILING


def test_a_typical_workstation_keeps_the_number_it_had(monkeypatch):
    """Continuity: the old constant was reasonable for the machine it was
    chosen on, and the formula must not quietly change that machine."""
    monkeypatch.setattr(B, "_available_cpus", lambda: 8)
    assert B._max_background_jobs() == 4


def test_an_explicit_setting_still_wins(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda *a, **kw: {"agent": {"max_background_jobs": 9}})
    monkeypatch.setattr(B, "_available_cpus", lambda: 2)
    assert B._max_background_jobs() == 9


def test_a_typo_in_the_setting_cannot_disable_the_cap(monkeypatch):
    for junk in ("", "many", None, 0, -3):
        monkeypatch.setattr("delfin.user_settings.load_settings",
                            lambda *a, **kw: {"agent": {"max_background_jobs": junk}})
        assert B._max_background_jobs() >= 1


def test_the_cap_survives_a_machine_that_answers_nothing(monkeypatch):
    def boom():
        raise OSError("no /proc here")
    monkeypatch.setattr(B, "_affinity_cpus", boom)
    assert B._max_background_jobs() >= 1
