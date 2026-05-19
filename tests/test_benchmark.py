"""Tests for the canned-task benchmark scorer + persistence + comparison.

The scorer is the heart of iterative agent optimisation: it must be
deterministic, must flag forbidden patterns, must give a stable 0-100
quality score, and must round-trip cleanly through JSONL.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import benchmark as bm


# ---------------------------------------------------------------------------
# Task loading
# ---------------------------------------------------------------------------


def test_load_tasks_packaged_suite_has_expected_tasks():
    tasks = bm.load_tasks()
    assert len(tasks) >= 8, f"packaged suite too small: {len(tasks)}"
    ids = {t.id for t in tasks}
    # Sanity-check a few well-known tasks
    assert "dash_nav_calc_typo" in ids
    assert "solo_research_setworking" in ids
    assert "debug_action_parse" in ids


def test_load_tasks_returns_empty_on_missing_path(tmp_path):
    assert bm.load_tasks(tmp_path / "nonexistent.yaml") == []


def test_load_tasks_coerces_signal_shorthand(tmp_path):
    """A bare-string signal under expected_signals must be coerced
    to a Signal(pattern=...) rather than crashing."""
    p = tmp_path / "tasks.yaml"
    p.write_text(
        "tasks:\n"
        "  - id: t1\n"
        "    task_class: misc\n"
        "    mode: solo\n"
        "    prompt: hi\n"
        "    expected_signals:\n"
        "      - hello\n",
        encoding="utf-8",
    )
    tasks = bm.load_tasks(p)
    assert len(tasks) == 1
    assert tasks[0].expected_signals[0].pattern == "hello"
    assert tasks[0].expected_signals[0].against == "any"


# ---------------------------------------------------------------------------
# Signal matching
# ---------------------------------------------------------------------------


def _make_task(**kw):
    defaults = dict(
        id="t", task_class="misc", mode="solo", prompt="p",
        expected_signals=(), forbidden_signals=(),
        max_duration_s=60.0, max_cost_usd=0.10, max_tool_calls=5,
    )
    defaults.update(kw)
    return bm.Task(**defaults)


def test_signal_against_action_only_scans_actions():
    task = _make_task(expected_signals=(
        bm.Signal(pattern=r"/tab\s+calc", against="action"),
    ))
    # Action present
    traj_hit = bm.Trajectory(text="some chatter", actions=["/tab calc"])
    r = bm.score_outcome(task, traj_hit)
    assert r.success is True
    # Same pattern in text but NOT in actions — must not satisfy
    traj_miss = bm.Trajectory(text="here is /tab calc", actions=[])
    r2 = bm.score_outcome(task, traj_miss)
    assert r2.success is False


def test_signal_against_tool_name_matches_tool_only():
    task = _make_task(expected_signals=(
        bm.Signal(pattern="Read", against="tool_name"),
    ))
    traj = bm.Trajectory(
        text="i'll read it", tool_calls=[{"name": "Read", "input": {}}],
    )
    assert bm.score_outcome(task, traj).success
    traj2 = bm.Trajectory(text="i'll Read it", tool_calls=[])
    assert not bm.score_outcome(task, traj2).success


def test_optional_signal_does_not_break_success():
    task = _make_task(expected_signals=(
        bm.Signal(pattern="must-have", against="text"),
        bm.Signal(pattern="nice-bonus", against="text", optional=True),
    ))
    traj = bm.Trajectory(text="must-have only")
    r = bm.score_outcome(task, traj)
    assert r.success is True
    # The optional miss is reported but doesn't kill the run
    assert any(":optional" in m for m in r.missing_signals)


def test_forbidden_signal_flips_success_to_false():
    task = _make_task(
        expected_signals=(bm.Signal(pattern="ok", against="text"),),
        forbidden_signals=(bm.Signal(
            pattern="(?i)cannot", against="text",
        ),),
    )
    traj = bm.Trajectory(text="ok but i cannot help")
    r = bm.score_outcome(task, traj)
    assert r.success is False
    assert r.violated_signals


def test_error_field_kills_success_even_if_all_signals_match():
    task = _make_task(expected_signals=(
        bm.Signal(pattern="ok", against="text"),
    ))
    traj = bm.Trajectory(text="ok", error="provider 500")
    r = bm.score_outcome(task, traj)
    assert r.success is False


# ---------------------------------------------------------------------------
# Quality 0-100 components
# ---------------------------------------------------------------------------


def test_quality_score_high_on_clean_success():
    task = _make_task(expected_signals=(
        bm.Signal(pattern="ok", against="text"),
    ), max_duration_s=10.0, max_cost_usd=0.10, max_tool_calls=3)
    traj = bm.Trajectory(
        text="ok done", duration_s=1.0, cost_usd=0.005,
        tool_calls=[{"name": "Read", "input": {}}],
    )
    r = bm.score_outcome(task, traj)
    assert r.quality_0_100 >= 85, (
        f"expected ≥85, got {r.quality_0_100} components={r.components}"
    )
    assert r.components["success_pts"] == 40


def test_quality_score_low_on_failure():
    task = _make_task(expected_signals=(
        bm.Signal(pattern="never-appears", against="text"),
    ))
    traj = bm.Trajectory(text="something else", duration_s=30.0)
    r = bm.score_outcome(task, traj)
    assert r.quality_0_100 <= 40
    assert r.components["success_pts"] == 0


def test_quality_score_penalises_duration_overage():
    task = _make_task(
        expected_signals=(bm.Signal(pattern="ok", against="text"),),
        max_duration_s=10.0,
    )
    fast = bm.score_outcome(task, bm.Trajectory(text="ok", duration_s=1.0))
    slow = bm.score_outcome(task, bm.Trajectory(text="ok", duration_s=9.5))
    assert slow.components["speed_pts"] < fast.components["speed_pts"]


def test_budget_violations_reported_but_not_fatal():
    """Going over duration/cost budget shouldn't flip success — the
    rubric reports it as a separate dimension."""
    task = _make_task(
        expected_signals=(bm.Signal(pattern="ok", against="text"),),
        max_duration_s=5.0, max_cost_usd=0.05, max_tool_calls=1,
    )
    traj = bm.Trajectory(text="ok", duration_s=20.0, cost_usd=0.20,
                         tool_calls=[{"name": "Read"}, {"name": "Bash"}])
    r = bm.score_outcome(task, traj)
    assert r.success is True
    assert len(r.budget_violations) == 3      # duration + cost + tool_calls


# ---------------------------------------------------------------------------
# Persistence round-trip
# ---------------------------------------------------------------------------


def test_write_and_read_run_roundtrip(tmp_path):
    task = _make_task(expected_signals=(
        bm.Signal(pattern="ok", against="text"),
    ))
    results = [
        bm.score_outcome(task, bm.Trajectory(text="ok", duration_s=1.0),
                         model="opus", profile_name="strong"),
        bm.score_outcome(task, bm.Trajectory(text="ko", duration_s=2.0),
                         model="opus", profile_name="strong"),
    ]
    path = bm.write_run(results, model="opus", runs_dir=tmp_path)
    assert path.exists()
    loaded = bm.read_run(path)
    assert len(loaded) == 2
    assert loaded[0]["model"] == "opus"
    assert loaded[0]["task_id"] == "t"
    assert loaded[0]["success"] is True
    assert loaded[1]["success"] is False


def test_read_run_handles_missing_file(tmp_path):
    assert bm.read_run(tmp_path / "nonexistent.jsonl") == []


# ---------------------------------------------------------------------------
# Run-vs-run comparison
# ---------------------------------------------------------------------------


def _result_dict(task_id, *, quality, success, cost=0.01, duration=2.0):
    return {
        "task_id": task_id, "model": "m", "profile_name": "",
        "mode": "solo", "ts": 0.0,
        "success": success, "quality_0_100": quality,
        "components": {}, "duration_s": duration, "cost_usd": cost,
        "input_tokens": 0, "output_tokens": 0, "tool_calls": 0,
        "matched_signals": [], "violated_signals": [],
        "missing_signals": [], "budget_violations": [], "error": "",
    }


def test_compare_runs_verdict_better_on_clear_quality_lift():
    baseline = [
        _result_dict("a", quality=60, success=True),
        _result_dict("b", quality=40, success=False),
        _result_dict("c", quality=70, success=True),
        _result_dict("d", quality=50, success=True),
    ]
    candidate = [
        _result_dict("a", quality=80, success=True),
        _result_dict("b", quality=70, success=True),
        _result_dict("c", quality=85, success=True),
        _result_dict("d", quality=75, success=True),
    ]
    cmp = bm.compare_runs(baseline, candidate)
    assert cmp["verdict"] == "better"
    assert cmp["summary"]["n_better"] == 4
    assert cmp["summary"]["n_worse"] == 0


def test_compare_runs_verdict_worse_on_regression():
    baseline = [_result_dict(f"t{i}", quality=80, success=True) for i in range(4)]
    candidate = [_result_dict(f"t{i}", quality=50, success=False) for i in range(4)]
    cmp = bm.compare_runs(baseline, candidate)
    assert cmp["verdict"] == "worse"
    assert cmp["summary"]["n_worse"] == 4


def test_compare_runs_verdict_thin_on_few_overlap():
    baseline = [_result_dict("a", quality=60, success=True)]
    candidate = [_result_dict("a", quality=80, success=True)]
    cmp = bm.compare_runs(baseline, candidate)
    assert cmp["verdict"] == "thin"


def test_compare_runs_credits_cost_drop_at_flat_quality():
    """If quality is unchanged but cost halved meaningfully, that's a win."""
    baseline = [
        _result_dict(f"t{i}", quality=70, success=True, cost=0.10)
        for i in range(4)
    ]
    candidate = [
        _result_dict(f"t{i}", quality=70, success=True, cost=0.02)
        for i in range(4)
    ]
    cmp = bm.compare_runs(baseline, candidate)
    assert cmp["verdict"] == "better"
    for row in cmp["per_task"]:
        assert row["class"] == "better"


def test_summarise_run_computes_basic_aggregates():
    rows = [
        _result_dict("a", quality=80, success=True,  cost=0.01, duration=1.0),
        _result_dict("b", quality=60, success=False, cost=0.02, duration=2.0),
        _result_dict("c", quality=90, success=True,  cost=0.03, duration=3.0),
    ]
    s = bm.summarise_run(rows)
    assert s["n_tasks"] == 3
    assert s["n_pass"] == 2
    assert s["pass_rate"] == pytest.approx(2 / 3)
    assert s["avg_quality"] == pytest.approx((80 + 60 + 90) / 3)
    assert s["total_cost_usd"] == pytest.approx(0.06)
    assert s["total_duration_s"] == pytest.approx(6.0)
