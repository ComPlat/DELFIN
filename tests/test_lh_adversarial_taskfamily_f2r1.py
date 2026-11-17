"""F2 task-family contract adversarial tests (reviewer nacht-s4).

These encode the contract that the F2 long-horizon task family
(tasks_w28_lh.yaml, target item 1) must satisfy so the yardstick F3 must
beat cannot silently be a partial package. They run against the READ-ONLY
benchmark harness (delfin/agent/benchmark.py, delfin/agent/benchmark_runner.py)
with tmp_path YAML fixtures built at runtime -- no fixture in a shared
knowledge dir, no live model, no ORCA data.

Findings these pin (each would be a real defect if the builder's yaml broke
the invariant):
  EXACT  mechanism: ExpectedValue.matches uses a RELATIVE margin
         (benchmark.py:272-273, margin = max(|value|*|tolerance|,
         |absolute|)). Default tolerance is 0.005 (benchmark.py:682), so a
         fact planted without tolerance:0/absolute:0 is 'matched' by a
         near-but-wrong figure. F2's spec says "exact numbers count only if
         exact" -- the task-family yaml must declare planted exact facts
         with tolerance:0 (or they false-positive a wrong late-quoted
         number). This test pins the mechanism it must use.
  LOAD   _load_one (benchmark.py:752-761) does safe_load + a list-comp over
         _coerce_task with NO per-task guard, and _coerce_task injests
         float(raw["value"])/raw["id"] (benchmark.py:680,704). A malformed
         F2 yaml among the tasks_*.yaml siblings therefore RAISES and breaks
         the WHOLE suite's load_tasks() -- the F2 yaml must be schema-valid,
         and the harness must fail loudly (not silently return empty).
  ONEARG run_acceptance (benchmark_runner.py:1332) invokes the accept script
         as [sys.executable, script, str(workspace)] -- exactly one argument.
         An F2 accept script built for a different signature breaks the
         harness contract.
"""

from __future__ import annotations

import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from delfin.agent.benchmark import ExpectedValue, load_tasks, _load_one
from delfin.agent.benchmark_runner import run_acceptance


# --- EXACT: near-but-wrong numbers must NOT count for exact facts ----------


def test_exact_fact_tolerance_zero_rejects_near_but_wrong():
    # A planted exact fact (tolerance:0, absolute:0) must match ONLY the
    # exact figure. Near-but-wrong must not count.
    exact = ExpectedValue(value=28.5, tolerance=0.0, absolute=0.0)
    assert not exact.matches(28.6)
    assert not exact.matches(28.4)
    assert not exact.matches(28.4999)
    assert exact.matches(28.5)


def test_default_tolerance_false_positives_near_but_wrong():
    # The DEFAULT tolerance (0.005) is a HALF-PERCENT band. A task-family
    # yaml that declares a planted exact fact WITHOUT tolerance:0/absolute:0
    # false-positives a near-but-wrong late-quoted number -- which is exactly
    # what the F2 spec forbids ("exact numbers count only if exact").
    # This test pins the failure mode so the builder must opt into exactness.
    lax = ExpectedValue(value=28.5)  # default tolerance 0.005, absolute 0.0
    # 28.6 is within 0.005*28.5 = 0.1425 of 28.5 -> falsely "matched"
    assert lax.matches(28.6)
    # Show the user what a half-percent band accepts.
    assert not lax.matches(28.9)  # just outside the band


def test_exact_fact_judge_rejects_wrong_text_absent_on_empty():
    # judge() must say "wrong" when the transcript states only a near-but
    # wrong figure for an exact fact, and "absent" when no figure at all.
    exact = ExpectedValue(value=1000.0, tolerance=0.0, absolute=0.0)
    assert exact.judge("the yield was 999.9") == "wrong"
    assert exact.judge("the yield was 1000") == "matched"
    assert exact.judge("no number here") == "absent"


def test_malformed_task_yaml_fails_loudly_not_silent(monkeypatch, tmp_path):
    # A malformed F2 yaml among the tasks_*.yaml siblings must RAISE in
    # _load_one, so the suite FAILS LOUDLY -- never silently return a
    # partial/empty task list (which would let a broken task family pass
    # as "green" when it simply never loaded its tasks).
    bad = tmp_path / "bad.yaml"
    # A task dict with a non-numeric "value" on an expected figure.
    bad.write_text(
        "tasks:\n"
        "  - id: lh_late_exact\n"
        "    task_class: long_horizon\n"
        "    mode: solo\n"
        "    prompt: plant fact\n"
        "    expected_values:\n"
        "      - value: not-a-number\n"
        "        tolerance: 0\n"
        "    max_tool_calls: 200\n",
        encoding="utf-8",
    )
    with pytest.raises((ValueError, TypeError, KeyError)):
        _load_one(bad)


def test_accept_script_name_is_contained_inside_the_pack(tmp_path):
    # run_acceptance resolves a task-named accept script inside the pack dir
    # (benchmark_runner.py:1199: <root>/delfin/agent/pack/benchmark/accept/),
    # and the child process is invoked as [sys.executable, script, ws] --
    # exactly one workspace argument (benchmark_runner.py:1350). A task yaml
    # may name a script, so a name that escapes the pack is data-driven path
    # traversal; acceptance_path must reject it. Pin containment.
    from delfin.agent.benchmark_runner import acceptance_path

    pack = tmp_path / "_root"
    pack.mkdir()
    base = acceptance_path("f2_sample_verify", root=pack)
    # acceptance_path returns <base>/<name> WITHOUT appending ".py"
    # (benchmark_runner.py:1200). The contract is containment: the
    # resolved path must stay inside the pack's accept dir.
    expect_dir = (pack / "delfin" / "agent" / "pack" / "benchmark" / "accept")
    assert str(base).startswith(str(expect_dir.resolve()))
    # a traversal name must be refused, not resolved outside the pack
    for bad in ("../escape", "..", "/../../etc/passwd"):
        with pytest.raises(ValueError):
            acceptance_path(bad, root=pack)
