"""The gate-form benchmark suite exists and measures the right forms.

Part of the 2026-09-25 finding: the agent reaches for `$( … )`
substitution, `awk`, `xargs` and `sed -i` out of HABIT, and every one
of those forms forces the gate to ask. The prompt now carries the rule
(pinned by test_looking_things_up_stays_in_the_reading_tools.py); this
suite is the other half — a benchmark that FAILS a run whose trajectory
contains those forms, so prompt regressions show up as score changes
before they show up as operator dialogs.

Control on the previous commit: 5 cases, all 5 red (the suite file did
not exist, so nothing loaded and no forbidden pattern was registered).
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from delfin.agent import benchmark as bm

_PACK_BENCH = (
    Path(__file__).resolve().parents[1]
    / "delfin" / "agent" / "pack" / "benchmark"
)
_SUITE = _PACK_BENCH / "tasks_gate_forms.yaml"


def _tasks() -> list[bm.Task]:
    assert _SUITE.exists(), "tasks_gate_forms.yaml is gone from the benchmark pack"
    return bm.load_tasks(_SUITE)


def test_the_suite_loads_and_is_tagged():
    tasks = _tasks()
    assert tasks, "tasks_gate_forms.yaml loaded zero tasks"
    ids = {t.id for t in tasks}
    assert ids == {
        "gate_forms_locate_function",
        "gate_forms_in_place_edit",
    }, f"unexpected task ids: {sorted(ids)}"
    assert all(t.mode == "solo" for t in tasks)
    assert all(t.task_class == "gate_forms" for t in tasks)


def test_every_gate_form_is_caught_by_some_signal():
    """A bad call that matches NO forbidden signal is a hole the suite
       reads as covered but does not measure."""
    tasks = {t.id: t for t in _tasks()}
    bad_calls = [
        "sed -n \"$(grep -n 'def x' f.py | cut -d: -f1),+40p\" f.py",
        "grep -n x f.py | awk -F: '{print $1}'",
        "cat list.txt | xargs grep -n y",
        "sed -i 's/a/b/' f.py",
    ]
    for tid, task in tasks.items():
        for cmd in bad_calls:
            traj = bm.Trajectory(tool_calls=[
                {"name": "bash", "input": "command: " + cmd}])
            assert any(
                bm._signal_matches(s, traj) for s in task.forbidden_signals
            ), f"{tid}: no forbidden signal catches {cmd!r}"


def test_every_signal_catches_a_live_bad_call():
    """And the converse: a forbidden pattern that matches nothing real
       would silently pass every run — the backtick signal was exactly
       that (benchmark.py's _strip_emphasis strips backticks from `any`
       haystacks before matching), caught by this check and removed from
       the suite."""
    bad_calls = [
        "sed -n \"$(grep -n 'def x' f.py | cut -d: -f1),+40p\" f.py",
        "grep -n x f.py | awk -F: '{print $1}'",
        "cat list.txt | xargs grep -n y",
        "sed -i 's/a/b/' f.py",
    ]
    for task in _tasks():
        for sig in task.forbidden_signals:
            hit = False
            for cmd in bad_calls:
                traj = bm.Trajectory(tool_calls=[
                    {"name": "bash", "input": "command: " + cmd}])
                if bm._signal_matches(sig, traj):
                    hit = True
                    break
            assert hit, (
                f"{task.id}: forbidden pattern {sig.pattern!r} matches no "
                "live gate-form command — it measures nothing"
            )


def test_the_reading_tools_are_not_forbidden():
    """The suite must punish the FORM, not reading itself: a trajectory
    that uses grep_file and read_file with offset/limit is exactly the
    prescribed behaviour and must match no forbidden signal."""
    tasks = _tasks()
    traj = bm.Trajectory(tool_calls=[
        {"name": "grep_file", "input": "pattern: def x, path: f.py"},
        {"name": "read_file", "input": "path: f.py, offset: 120, limit: 40"},
        {"name": "edit_file",
         "input": "path: f.py, old_string: a, new_string: b"},
    ])
    for task in tasks:
        for sig in task.forbidden_signals:
            assert not bm._signal_matches(sig, traj), (
                f"{task.id}: {sig.pattern!r} fires on the prescribed "
                "grep_file/read_file/edit_file trajectory"
            )


def test_the_task_ids_are_unique_in_the_full_suite():
    """load_tasks() concatenates every tasks_*.yaml; a duplicate id
    across files would silently shadow one task in live runs."""
    all_ids = [t.id for t in bm.load_tasks()]
    assert len(all_ids) == len(set(all_ids)), "duplicate task ids"
    assert set(t.id for t in _tasks()) <= set(all_ids)
