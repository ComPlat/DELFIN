"""A contained command that multiplies past its budget is ended, and says so.

On 2026-09-25 one gate run fanned out into ~75 nested login shells on a
shared login node; the timeout could not help (each shell was quick) and
only the operator's outside watcher noticed. `contained_run.run(budget=
True)` -- how the agent's bash tool and the sandbox run every command --
now watches the call's tree against `process_budget` and ends the group.
"""

from __future__ import annotations

import time

from delfin.agent import contained_run, process_budget

_FAN_OUT = "for i in 1 2 3 4 5 6 7 8; do sleep 30 & done; wait"


def _small(monkeypatch, n=4):
    monkeypatch.setattr(
        process_budget, "_default_limits",
        lambda: process_budget.BudgetLimits(max_processes=n))


def test_a_fan_out_past_the_budget_is_ended_and_named(monkeypatch):
    _small(monkeypatch)
    t0 = time.monotonic()
    done = contained_run.run(["bash", "-c", _FAN_OUT], timeout=60,
                             budget=True)
    assert time.monotonic() - t0 < 20, "the fan-out ran on past its budget"
    assert "[process budget]" in done.stderr
    assert "max_processes" in done.stderr


def test_within_the_budget_nothing_changes(monkeypatch):
    _small(monkeypatch, n=16)
    done = contained_run.run(["bash", "-c", "echo hi; sleep 0.2"],
                             timeout=30, budget=True)
    assert done.returncode == 0
    assert done.stdout.strip() == "hi"
    assert "[process budget]" not in done.stderr


def test_without_budget_the_fan_out_is_not_watched(monkeypatch):
    _small(monkeypatch)
    done = contained_run.run(
        ["bash", "-c", "for i in 1 2 3 4 5 6; do sleep 1 & done; wait"],
        timeout=30)
    assert done.returncode == 0
    assert "[process budget]" not in done.stderr
