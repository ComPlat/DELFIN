"""The agent sized its work from a remembered number, not a measured one.

Machine facts were written down once, on one host, and then recalled
everywhere. They were wrong the moment the session moved -- and they were
recorded on a host that was not the one the work runs on, which is how
the question came up.

The cap on background jobs already derives from the real grant
(``sched_getaffinity``, or the SLURM allocation when there is one), so
the MECHANISM was already honest. What was missing is that the model
never saw it, so anything the model sizes itself -- ``pytest -n``, a PAL
value, an sbatch header -- was a guess against the machine's total.

Stated once per session, not per turn: a grant does not change while a
session runs, and the per-turn steering budget is paid for by blocks
that do change. The line is short for the same reason.
"""

from __future__ import annotations

import pytest

from delfin.agent.engine import AgentEngine


def _block(monkeypatch, *, cpus: int, slurm: bool, first_turn: bool = True):
    from delfin.agent import bash_jobs as bj
    monkeypatch.setattr(bj, "_available_cpus", lambda: cpus)
    monkeypatch.setattr(bj, "_affinity_cpus", lambda: cpus if not slurm else 999)
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = [{"role": "user", "content": "go"}] if first_turn else [
        {"role": "user", "content": "a"}, {"role": "assistant", "content": "b"},
        {"role": "user", "content": "c"}]
    return eng._build_machine_grant_block()


def test_the_measured_grant_is_stated(monkeypatch):
    assert "48" in _block(monkeypatch, cpus=48, slurm=False)


def test_it_says_what_the_number_is_for(monkeypatch):
    text = _block(monkeypatch, cpus=48, slurm=False).lower()
    assert "cpu" in text or "core" in text


def test_an_allocation_is_named_as_such(monkeypatch):
    """A SLURM grant and a whole idle machine are not the same situation."""
    assert "slurm" in _block(monkeypatch, cpus=8, slurm=True).lower()


def test_a_plain_machine_does_not_claim_an_allocation(monkeypatch):
    assert "slurm" not in _block(monkeypatch, cpus=48, slurm=False).lower()


def test_it_is_said_once_not_every_turn(monkeypatch):
    """The grant does not change while a session runs; the steering budget
    is paid for by blocks that do."""
    assert _block(monkeypatch, cpus=48, slurm=False, first_turn=False) == ""


def test_it_stays_short(monkeypatch):
    """It rides on every request of the first turn."""
    assert len(_block(monkeypatch, cpus=48, slurm=False)) < 240


def test_a_failure_to_measure_says_nothing(monkeypatch):
    from delfin.agent import bash_jobs as bj

    def _boom():
        raise OSError("no /proc")
    monkeypatch.setattr(bj, "_available_cpus", _boom)
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = [{"role": "user", "content": "go"}]
    assert eng._build_machine_grant_block() == ""


def test_a_single_cpu_is_still_worth_saying(monkeypatch):
    """One core is exactly when a parallel flag does damage."""
    assert "1" in _block(monkeypatch, cpus=1, slurm=False)


def test_the_block_is_part_of_the_steering(monkeypatch):
    import pathlib
    src = (pathlib.Path(__file__).resolve().parents[1] / "delfin" / "agent"
           / "engine.py").read_text(encoding="utf-8")
    i = src.index("def _steering_blocks(")
    j = src.index("def _drain_turn_steering(")
    assert "_build_machine_grant_block" in src[i:j], (
        "the block exists but is never delivered")
