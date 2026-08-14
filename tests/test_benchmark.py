"""Tests for the canned-task benchmark scorer + persistence + comparison.

The scorer is the heart of iterative agent optimisation: it must be
deterministic, must flag forbidden patterns, must give a stable 0-100
quality score, and must round-trip cleanly through JSONL.
"""

from __future__ import annotations

import time
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


def test_run_timestamp_reads_first_record(tmp_path):
    p = tmp_path / "run.jsonl"
    p.write_text(
        '{"task_id":"a","ts":1700000000.5,"success":true,"quality_0_100":80}\n'
        '{"task_id":"b","ts":1700000900.0,"success":true,"quality_0_100":70}\n',
        encoding="utf-8",
    )
    assert bm.run_timestamp(p) == pytest.approx(1700000000.5)


def test_run_timestamp_returns_zero_on_missing(tmp_path):
    assert bm.run_timestamp(tmp_path / "nope.jsonl") == 0.0


# ---------------------------------------------------------------------------
# format_compare_markdown
# ---------------------------------------------------------------------------


def _build_compare(verdict="better"):
    """Helper: a minimal compare_runs result for markdown tests."""
    return {
        "verdict": verdict,
        "summary": {
            "n_overlap": 4, "n_better": 4, "n_worse": 0, "n_neutral": 0,
            "old": {"n_tasks": 4, "n_pass": 3, "pass_rate": 0.75,
                    "avg_quality": 65.0, "total_cost_usd": 0.40,
                    "total_duration_s": 30.0, "total_tool_calls": 4},
            "new": {"n_tasks": 4, "n_pass": 4, "pass_rate": 1.00,
                    "avg_quality": 82.5, "total_cost_usd": 0.20,
                    "total_duration_s": 18.0, "total_tool_calls": 4},
        },
        "per_task": [
            {"task_id": "a", "class": "better",
             "old_quality": 60, "new_quality": 80, "d_quality": 20,
             "d_cost_usd": -0.05, "d_duration_s": -2.0, "d_tool_calls": 0,
             "old_success": True, "new_success": True},
            {"task_id": "b", "class": "better",
             "old_quality": 50, "new_quality": 70, "d_quality": 20,
             "d_cost_usd": -0.05, "d_duration_s": -2.0, "d_tool_calls": 0,
             "old_success": False, "new_success": True},
        ],
    }


def test_format_compare_markdown_contains_verdict_header():
    md = bm.format_compare_markdown(_build_compare("better"))
    assert md.startswith("# Benchmark Comparison Report")
    assert "Verdict: BETTER" in md


def test_format_compare_markdown_contains_summary_table():
    md = bm.format_compare_markdown(_build_compare())
    assert "## Summary" in md
    assert "Pass rate" in md
    assert "75%" in md and "100%" in md
    assert "Avg quality" in md
    assert "Total cost" in md


def test_format_compare_markdown_contains_per_task_rows():
    md = bm.format_compare_markdown(_build_compare())
    assert "## Per-task" in md
    assert "`a" in md and "`b" in md
    assert "60→80" in md
    assert "+20" in md


def test_format_compare_markdown_handles_empty_per_task():
    cmp = _build_compare()
    cmp["per_task"] = []
    md = bm.format_compare_markdown(cmp)
    assert "## Per-task" not in md
    assert "Verdict" in md


def test_format_compare_markdown_thin_verdict_renders():
    cmp = _build_compare("thin")
    md = bm.format_compare_markdown(cmp)
    assert "Verdict: THIN" in md


def test_format_compare_markdown_skips_profile_block_without_paths():
    md = bm.format_compare_markdown(_build_compare())
    assert "Profile changes" not in md


# ---------------------------------------------------------------------------
# find_profile_commits_between
# ---------------------------------------------------------------------------


def test_find_profile_commits_rejects_bad_window():
    """Inverted window / zero timestamps must return [] without git call."""
    assert bm.find_profile_commits_between(0, 0) == []
    assert bm.find_profile_commits_between(100, 50) == []
    assert bm.find_profile_commits_between(-1, 100) == []


def test_find_profile_commits_returns_list_on_real_git(tmp_path, monkeypatch):
    """Run against a tiny throwaway git repo with one commit to the
    profile-file path — must list it."""
    import subprocess as _sp
    repo = tmp_path / "tinyrepo"
    repo.mkdir()
    _sp.run(["git", "init", "-q"], cwd=repo, check=True)
    _sp.run(["git", "config", "user.email", "a@b.c"], cwd=repo, check=True)
    _sp.run(["git", "config", "user.name", "test"], cwd=repo, check=True)
    profile_dir = repo / "delfin" / "agent"
    profile_dir.mkdir(parents=True)
    profile_file = profile_dir / "model_profiles.py"
    profile_file.write_text("# stub\n", encoding="utf-8")
    _sp.run(["git", "add", "."], cwd=repo, check=True)
    _sp.run(["git", "commit", "-q", "-m", "profile: tweak knob"],
            cwd=repo, check=True)
    # Window covering the commit
    now = time.time()
    commits = bm.find_profile_commits_between(
        now - 3600, now + 3600, repo_root=repo,
    )
    assert any("profile: tweak knob" in c for c in commits)


def test_find_profile_commits_returns_empty_on_non_repo(tmp_path):
    assert bm.find_profile_commits_between(
        time.time() - 3600, time.time(), repo_root=tmp_path,
    ) == []


# ---------------------------------------------------------------------------
# format_compare_markdown with paths → profile-commit annotation
# ---------------------------------------------------------------------------


def test_format_compare_markdown_includes_profile_commits(tmp_path):
    """End-to-end: when both baseline + candidate paths are supplied and
    git has commits to model_profiles.py in between, the report lists
    them."""
    import subprocess as _sp
    repo = tmp_path / "repo"
    repo.mkdir()
    _sp.run(["git", "init", "-q"], cwd=repo, check=True)
    _sp.run(["git", "config", "user.email", "a@b.c"], cwd=repo, check=True)
    _sp.run(["git", "config", "user.name", "test"], cwd=repo, check=True)
    (repo / "delfin" / "agent").mkdir(parents=True)
    pf = repo / "delfin" / "agent" / "model_profiles.py"
    pf.write_text("# v1\n", encoding="utf-8")
    _sp.run(["git", "add", "."], cwd=repo, check=True)
    _sp.run(["git", "commit", "-q", "-m", "profile: initial"],
            cwd=repo, check=True)

    # Write a "baseline" run BEFORE the next commit
    runs = tmp_path / "runs"
    runs.mkdir()
    base_path = runs / "baseline.jsonl"
    base_path.write_text(
        '{"task_id":"a","ts":1000.0,"success":true,"quality_0_100":50,'
        '"cost_usd":0.1,"duration_s":2.0,"tool_calls":1}\n',
        encoding="utf-8",
    )
    # Bump profile + commit at "now" so git sees it
    pf.write_text("# v2 — tighter stale-kill\n", encoding="utf-8")
    _sp.run(["git", "add", "."], cwd=repo, check=True)
    _sp.run(["git", "commit", "-q", "-m", "profile: tighter stale-kill"],
            cwd=repo, check=True)
    # Candidate run "now"
    now = time.time()
    cand_path = runs / "candidate.jsonl"
    import json as _json
    cand_path.write_text(
        _json.dumps({
            "task_id": "a", "ts": now, "success": True,
            "quality_0_100": 80, "cost_usd": 0.05,
            "duration_s": 1.0, "tool_calls": 1,
        }) + "\n",
        encoding="utf-8",
    )
    # We bracket the second commit
    base_ts = bm.run_timestamp(base_path)
    cand_ts = bm.run_timestamp(cand_path)
    commits = bm.find_profile_commits_between(
        base_ts, cand_ts, repo_root=repo,
    )
    # The bracketed commit must appear
    assert any("tighter stale-kill" in c for c in commits)

    # And the full markdown report should include it
    cmp = _build_compare("better")
    md = bm.format_compare_markdown(
        cmp, baseline_path=base_path, candidate_path=cand_path,
        repo_root=repo,
    )
    assert "Profile changes between runs" in md
    assert "tighter stale-kill" in md


# ---------------------------------------------------------------------------
# Audit — pattern-bug-vs-real-fail diagnosis
# ---------------------------------------------------------------------------


def test_score_outcome_captures_text_excerpt():
    task = bm.Task(
        id="t", task_class="misc", mode="solo", prompt="hi",
        expected_signals=(bm.Signal(pattern="never_matches", against="text"),),
        forbidden_signals=(),
        max_duration_s=10.0, max_cost_usd=0.05, max_tool_calls=2,
    )
    traj = bm.Trajectory(
        text="Hier ist meine Antwort: 1. erst dies, 2. dann das",
        tool_calls=[{"name": "Read", "input": {}}],
        duration_s=1.0,
    )
    r = bm.score_outcome(task, traj)
    assert r.success is False
    assert "Hier ist meine Antwort" in r.text_excerpt
    assert r.tool_names == ["Read"]


def test_score_outcome_text_excerpt_truncates_to_400_chars():
    task = bm.Task(
        id="t", task_class="misc", mode="solo", prompt="hi",
        expected_signals=(),
        forbidden_signals=(),
        max_duration_s=10.0, max_cost_usd=0.05, max_tool_calls=2,
    )
    traj = bm.Trajectory(text="x" * 1000, duration_s=1.0)
    r = bm.score_outcome(task, traj)
    assert len(r.text_excerpt) <= 400


def test_audit_run_flags_pattern_bug_when_excerpt_reasonable():
    """rate=0% + σ≈0 + non-empty excerpt + no error → likely pattern bug."""
    record = {
        "task_id": "t1", "model": "m", "success": False,
        "quality_0_100": 35, "success_rate": 0.0, "quality_stdev": 0.6,
        "missing_signals": ["t1.expected[0]"],
        "violated_signals": [],
        "text_excerpt": "ACTION: /tab calc\nACTION: /done",
        "tool_names": [], "error": "",
    }
    entries = bm.audit_run([record], tasks=[])
    assert len(entries) == 1
    assert entries[0]["hint_pattern_bug"] is True


def test_audit_run_does_not_flag_when_error_present():
    """engine init failed → that's a real fail, not a pattern bug."""
    record = {
        "task_id": "t2", "model": "m", "success": False,
        "quality_0_100": 35, "success_rate": 0.0, "quality_stdev": 0.0,
        "missing_signals": ["t2.expected[0]"],
        "violated_signals": [],
        "text_excerpt": "",
        "tool_names": [], "error": "engine init failed: no API key",
    }
    entries = bm.audit_run([record], tasks=[])
    assert entries[0]["hint_pattern_bug"] is False


def test_audit_run_does_not_flag_when_forbidden_violated():
    """If model violated a forbidden signal, that's a real fail
    (explicit misbehaviour), not a pattern bug — even if the rest of
    the heuristic checkboxes would otherwise tick."""
    record = {
        "task_id": "t_violated", "model": "m", "success": False,
        "quality_0_100": 26, "success_rate": 0.0, "quality_stdev": 0.0,
        "missing_signals": ["t_violated.expected[0]"],
        "violated_signals": ["t_violated.forbidden[0]"],
        "text_excerpt": "ACTION: /tab qwertyzzzz\nACTION: /done",
        "tool_names": [], "error": "",
    }
    entries = bm.audit_run([record], tasks=[])
    assert entries[0]["hint_pattern_bug"] is False


def test_audit_run_does_not_flag_flaky_as_pattern_bug():
    """high σ → noise, not pattern bug."""
    record = {
        "task_id": "t3", "model": "m", "success": False,
        "quality_0_100": 35, "success_rate": 0.33, "quality_stdev": 25.0,
        "missing_signals": ["t3.expected[0]"],
        "violated_signals": [],
        "text_excerpt": "some response",
        "tool_names": ["Read"], "error": "",
    }
    entries = bm.audit_run([record], tasks=[])
    assert entries[0]["hint_pattern_bug"] is False


def test_audit_run_skips_passing_records():
    record = {
        "task_id": "ok", "model": "m", "success": True,
        "quality_0_100": 90, "success_rate": 1.0, "quality_stdev": 1.0,
        "missing_signals": [], "violated_signals": [],
        "text_excerpt": "fine", "tool_names": [], "error": "",
    }
    entries = bm.audit_run([record], tasks=[])
    assert entries == []


def test_audit_run_includes_signal_definitions_from_tasks():
    """When task definitions are provided, audit must surface the
    regex pattern that failed so the dev can see what was expected."""
    task = bm.Task(
        id="t4", task_class="misc", mode="solo", prompt="x",
        expected_signals=(
            bm.Signal(pattern=r"never_matches", against="text"),
        ),
        forbidden_signals=(),
        max_duration_s=10.0, max_cost_usd=0.05, max_tool_calls=2,
    )
    record = {
        "task_id": "t4", "model": "m", "success": False,
        "quality_0_100": 35, "success_rate": 0.0, "quality_stdev": 0.0,
        "missing_signals": ["t4.expected[0]"],
        "violated_signals": [],
        "text_excerpt": "some output", "tool_names": [], "error": "",
    }
    entries = bm.audit_run([record], tasks=[task])
    assert "never_matches" in entries[0]["signal_defs"]["t4.expected[0]"]


def test_format_audit_report_groups_bug_vs_real():
    bug = {
        "task_id": "bug1", "model": "m", "quality": 35,
        "success_rate": 0.0, "quality_stdev": 0.6,
        "missing_signals": ["bug1.expected[0]"],
        "violated_signals": [], "tool_names": [],
        "text_excerpt": "the model said something reasonable",
        "error": "", "signal_defs": {"bug1.expected[0]": "/foo   (against=text)"},
        "hint_pattern_bug": True,
    }
    real = {
        "task_id": "real1", "model": "m", "quality": 27,
        "success_rate": 0.0, "quality_stdev": 0.0,
        "missing_signals": ["real1.expected[0]"],
        "violated_signals": ["real1.forbidden[0]"],
        "tool_names": [], "text_excerpt": "",
        "error": "", "signal_defs": {},
        "hint_pattern_bug": False,
    }
    out = bm.format_audit_report([bug, real])
    assert "SUSPECTED PATTERN BUG" in out
    assert "LIKELY REAL FAIL OR FLAKY" in out
    assert "bug1" in out
    assert "real1" in out


def test_format_audit_report_empty_when_no_failures():
    assert "nothing to audit" in bm.format_audit_report([])


def test_aggregate_replicates_propagates_text_excerpt():
    """When ≥1 sample has text, aggregate must pick it (not blank)."""
    r1 = bm.BenchmarkResult(
        task_id="t", task_class="misc", model="m", mode="solo",
        ts=0.0, success=False, quality_0_100=35,
        text_excerpt="", tool_names=[],
    )
    r2 = bm.BenchmarkResult(
        task_id="t", task_class="misc", model="m", mode="solo",
        ts=0.0, success=False, quality_0_100=35,
        text_excerpt="hello world", tool_names=["Read"],
    )
    agg = bm.aggregate_replicates([r1, r2])
    assert agg.text_excerpt == "hello world"
    assert "Read" in agg.tool_names


def _make_result(task_id, *, quality, success=True, cost=0.01,
                 duration=2.0, tool_calls=1, error=""):
    return bm.BenchmarkResult(
        task_id=task_id, task_class="misc", model="m",
        mode="solo", ts=1700000000.0,
        success=success, quality_0_100=quality,
        components={"success_pts": 40 if success else 0,
                    "routing_pts": 15, "speed_pts": 10,
                    "cost_pts": 8, "clean_pts": 12},
        duration_s=duration, cost_usd=cost, tool_calls=tool_calls,
        error=error,
    )


# ---------------------------------------------------------------------------
# aggregate_replicates — N=3 noise-defeat
# ---------------------------------------------------------------------------


def test_aggregate_replicates_requires_nonempty_input():
    with pytest.raises(ValueError):
        bm.aggregate_replicates([])


def test_aggregate_replicates_requires_same_task_id():
    a = _make_result("task_a", quality=80)
    b = _make_result("task_b", quality=70)
    with pytest.raises(ValueError):
        bm.aggregate_replicates([a, b])


def test_aggregate_replicates_quality_uses_median():
    """Median, not mean, so one outlier doesn't swing the verdict."""
    reps = [
        _make_result("t", quality=q) for q in (60, 80, 30)
    ]
    agg = bm.aggregate_replicates(reps)
    # Median of [30, 60, 80] is 60
    assert agg.quality_0_100 == 60
    # stdev of [60, 80, 30] (sample) is sqrt(varvar) ≈ 25.2
    assert 20 < agg.quality_stdev < 30


def test_aggregate_replicates_majority_success_passes():
    """Majority of N succeed → aggregate.success = True."""
    reps = [
        _make_result("t", quality=80, success=True),
        _make_result("t", quality=80, success=True),
        _make_result("t", quality=40, success=False),
    ]
    agg = bm.aggregate_replicates(reps)
    assert agg.success is True
    assert agg.success_rate == pytest.approx(2 / 3)


def test_aggregate_replicates_minority_success_fails():
    reps = [
        _make_result("t", quality=40, success=False),
        _make_result("t", quality=40, success=False),
        _make_result("t", quality=80, success=True),
    ]
    agg = bm.aggregate_replicates(reps)
    assert agg.success is False
    assert agg.success_rate == pytest.approx(1 / 3)


def test_aggregate_replicates_ties_resolve_to_pass():
    """N=2 with 1 pass + 1 fail → pass (majority is forgiving)."""
    reps = [
        _make_result("t", quality=80, success=True),
        _make_result("t", quality=40, success=False),
    ]
    agg = bm.aggregate_replicates(reps)
    assert agg.success is True
    assert agg.success_rate == 0.5


def test_aggregate_replicates_cost_is_summed():
    """Cost is REAL spend across N runs — must SUM, not median."""
    reps = [
        _make_result("t", quality=70, cost=0.10),
        _make_result("t", quality=70, cost=0.20),
        _make_result("t", quality=70, cost=0.30),
    ]
    agg = bm.aggregate_replicates(reps)
    assert agg.cost_usd == pytest.approx(0.60)


def test_aggregate_replicates_duration_uses_median():
    reps = [
        _make_result("t", quality=70, duration=1.0),
        _make_result("t", quality=70, duration=2.0),
        _make_result("t", quality=70, duration=100.0),
    ]
    agg = bm.aggregate_replicates(reps)
    # Median of [1, 2, 100] is 2, not (1+2+100)/3 = 34.3
    assert agg.duration_s == pytest.approx(2.0)


def test_aggregate_replicates_keeps_per_run_history():
    reps = [_make_result("t", quality=q, success=q >= 50)
            for q in (40, 70, 90)]
    agg = bm.aggregate_replicates(reps)
    assert agg.per_run_quality == [40, 70, 90]
    assert agg.per_run_success == [False, True, True]
    assert agg.n_samples == 3


def test_aggregate_replicates_single_sample_trivial():
    """N=1 should round-trip cleanly: same quality, stdev=0,
    success_rate matches the single value."""
    r = _make_result("t", quality=75, success=True)
    agg = bm.aggregate_replicates([r])
    assert agg.quality_0_100 == 75
    assert agg.quality_stdev == 0.0
    assert agg.success_rate == 1.0
    assert agg.n_samples == 1


def test_aggregate_replicates_unions_signals():
    """Flaky violations (only appearing in 1 of N) must still be
    surfaced in the aggregate so they don't go unnoticed."""
    r1 = _make_result("t", quality=70)
    r1.matched_signals = ["t.expected[0]"]
    r2 = _make_result("t", quality=70)
    r2.matched_signals = ["t.expected[0]"]
    r2.violated_signals = ["t.forbidden[0]"]    # only in run 2
    agg = bm.aggregate_replicates([r1, r2])
    assert "t.expected[0]" in agg.matched_signals
    assert "t.forbidden[0]" in agg.violated_signals


def test_aggregate_replicates_jsonl_roundtrip(tmp_path):
    """An aggregated result must round-trip through write_run/read_run
    without losing the new fields."""
    reps = [_make_result("t", quality=q) for q in (60, 80, 70)]
    agg = bm.aggregate_replicates(reps)
    path = bm.write_run([agg], model="x", runs_dir=tmp_path)
    loaded = bm.read_run(path)
    assert len(loaded) == 1
    assert loaded[0]["n_samples"] == 3
    assert loaded[0]["per_run_quality"] == [60, 80, 70]
    assert loaded[0]["quality_stdev"] > 0


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


# ---------------------------------------------------------------------------
# ab_compare — compact A/B convenience over two persisted runs
# ---------------------------------------------------------------------------


def _write_ab_runs(tmp_path):
    """Two synthetic runs: candidate lifts quality on all 3 tasks, spends
    more, and flips the 'planned' behaviour flag on."""
    base = []
    cand = []
    for i in range(3):
        b = _make_result(f"t{i}", quality=50, success=False, cost=0.01)
        c = _make_result(f"t{i}", quality=80, success=True, cost=0.02)
        b.behavior = {"planned": 0.0}
        c.behavior = {"planned": 1.0}
        base.append(b)
        cand.append(c)
    pa = bm.write_run(base, model="m", runs_dir=tmp_path, run_id="base")
    pb = bm.write_run(cand, model="m", runs_dir=tmp_path, run_id="cand")
    return pa, pb


def test_ab_compare_on_two_synthetic_runs(tmp_path):
    _write_ab_runs(tmp_path)
    ab = bm.ab_compare("base", "cand", runs_dir=tmp_path)   # bare run ids
    assert ab["verdict"] == "better"
    assert ab["n_overlap"] == 3 and ab["n_better"] == 3
    assert ab["per_task_delta"]["t0"] == {
        "d_quality": 30, "d_cost_usd": pytest.approx(0.01), "class": "better",
    }
    assert ab["score_delta"] == pytest.approx(30.0)
    assert ab["cost_delta"] == pytest.approx(0.03)
    ch = ab["behaviour_flag_changes"]["planned"]
    assert ch["old"] == 0.0 and ch["new"] == 1.0
    assert ch["delta"] == pytest.approx(1.0)


def test_ab_compare_accepts_full_paths_and_filenames(tmp_path):
    pa, pb = _write_ab_runs(tmp_path)
    by_path = bm.ab_compare(str(pa), str(pb), runs_dir=tmp_path / "elsewhere")
    by_name = bm.ab_compare("base.jsonl", "cand.jsonl", runs_dir=tmp_path)
    assert by_path["verdict"] == by_name["verdict"] == "better"
    assert by_path["per_task_delta"] == by_name["per_task_delta"]


def test_ab_compare_missing_runs_is_thin_not_crash(tmp_path):
    ab = bm.ab_compare("nope_a", "nope_b", runs_dir=tmp_path)
    assert ab["verdict"] == "thin"
    assert ab["per_task_delta"] == {}
    assert ab["behaviour_flag_changes"] == {}
    assert ab["score_delta"] == 0.0 and ab["cost_delta"] == 0.0


def test_ab_compare_unchanged_flags_are_omitted(tmp_path):
    base = [_make_result("t", quality=50)]
    cand = [_make_result("t", quality=60)]
    base[0].behavior = {"scouted": 1.0}
    cand[0].behavior = {"scouted": 1.0}                 # no change
    bm.write_run(base, model="m", runs_dir=tmp_path, run_id="a")
    bm.write_run(cand, model="m", runs_dir=tmp_path, run_id="b")
    ab = bm.ab_compare("a", "b", runs_dir=tmp_path)
    assert ab["behaviour_flag_changes"] == {}


def test_format_ab_note_renders_verdict_deltas_and_flags(tmp_path):
    _write_ab_runs(tmp_path)
    note = bm.format_ab_note(bm.ab_compare("base", "cand", runs_dir=tmp_path))
    assert "verdict: BETTER (3 better / 0 worse / 0 neutral, n=3)" in note
    assert "score delta: +30.0" in note
    assert "cost delta: $+0.0300" in note
    assert "planned 0%→100%" in note
    assert "regressed" not in note                      # nothing regressed


def test_format_ab_note_lists_regressed_tasks():
    ab = {
        "verdict": "worse", "n_better": 0, "n_worse": 1, "n_neutral": 2,
        "n_overlap": 3, "score_delta": -12.0, "cost_delta": 0.0,
        "per_task_delta": {"bad_task": {"d_quality": -30, "d_cost_usd": 0.0,
                                        "class": "worse"}},
        "behaviour_flag_changes": {},
    }
    note = bm.format_ab_note(ab)
    assert "verdict: WORSE" in note
    assert "regressed: bad_task" in note


def test_list_runs_sorted_oldest_first(tmp_path):
    import os
    pa = bm.write_run([_make_result("t", quality=1)], model="m",
                      runs_dir=tmp_path, run_id="newer")
    pb = bm.write_run([_make_result("t", quality=1)], model="m",
                      runs_dir=tmp_path, run_id="older")
    os.utime(pa, (2000, 2000))
    os.utime(pb, (1000, 1000))
    assert [p.name for p in bm.list_runs(tmp_path)] == [
        "older.jsonl", "newer.jsonl"]
    assert bm.list_runs(tmp_path / "missing") == []


# ---------------------------------------------------------------------------
# Negation-aware forbidden scoring + task-cap run budgets
# ---------------------------------------------------------------------------

def _forbidden_task(pattern):
    from delfin.agent.benchmark import Task, Signal
    return Task(
        id="t_forbid", task_class="fact", mode="solo", prompt="p",
        expected_signals=(Signal(pattern="(?i)nel", against="text"),),
        forbidden_signals=(Signal(pattern=pattern, against="text"),),
        max_duration_s=60, max_cost_usd=0.3, max_tool_calls=8,
    )


def test_forbidden_hit_still_fails_plain_recommendation():
    from delfin.agent.benchmark import Trajectory, score_outcome
    traj = Trajectory(text="Use the keyword Nactel 6 and nel 6 here.")
    res = score_outcome(_forbidden_task(r"(?i)\bnactel\b"), traj,
                        model="m", profile_name="p", ts=0.0)
    assert res.violated_signals


def test_forbidden_hit_waived_in_negation_context():
    """An answer that names a fake keyword to WARN against it shows the
    grounded behavior the suite rewards — no violation."""
    from delfin.agent.benchmark import Trajectory, score_outcome
    traj = Trajectory(text=(
        "Die Keywords sind `nel` und `norb`. Wichtig: Die Keywords "
        "heißen NICHT `Nactel` oder `Nactorb` — das sind erfundene "
        "Namen, die es im Manual nicht gibt."))
    res = score_outcome(_forbidden_task(r"(?i)\bnactel\b"), traj,
                        model="m", profile_name="p", ts=0.0)
    assert not res.violated_signals


def test_forbidden_negation_elsewhere_does_not_launder():
    """A negation far away in the text must not waive a genuine
    recommendation of the forbidden term."""
    from delfin.agent.benchmark import Trajectory, score_outcome
    traj = Trajectory(text=(
        "This approach is not ideal for large systems. " + "x" * 300 +
        " Set Nactel 6 in the %casscf block to configure it, "
        "then run the calculation with the usual settings and confirm "
        "the reference weights afterwards."))
    res = score_outcome(_forbidden_task(r"(?i)\bnactel\b"), traj,
                        model="m", profile_name="p", ts=0.0)
    assert res.violated_signals


def test_runner_sets_task_cap_budgets(monkeypatch):
    from delfin.agent import benchmark_runner as br
    from delfin.agent.benchmark import Task

    class _Eng:
        cost_usd = 0.0
    captured = {}

    def _factory(model, backend, provider, mode, task_class=""):
        e = _Eng(); captured["engine"] = e; return e

    def _run_once(engine, prompt, max_tokens=4096):
        return {"text": "nel 6", "tool_calls": [], "error": "",
                "input_tokens": 1, "output_tokens": 1}

    task = Task(id="t_b", task_class="fact", mode="solo", prompt="p",
                expected_signals=(), forbidden_signals=(),
                max_duration_s=90, max_cost_usd=0.3, max_tool_calls=8)
    br.run_task(task, model="m", backend="api", provider="kit",
                profile_name="p", engine_factory=_factory,
                max_tokens=100, run_once=_run_once)
    eng = captured["engine"]
    assert abs(eng.run_budget_usd - 1.2) < 1e-9
    assert eng.run_budget_s == 360.0


def test_negation_window_does_not_slice_words():
    """A window edge inside a word invented a word boundary: 'notes' cut
    to 'not' matched the negation marker and waived a real violation
    (found while authoring the generic-project tasks)."""
    from delfin.agent.benchmark import _match_is_negated, _NEGATION_WINDOW
    filler = "x" * (_NEGATION_WINDOW - 3)
    # 'notes' sits just outside the window; only its tail would be sliced in.
    hay = f"notes {filler} FORBIDDEN tail"
    start = hay.index("FORBIDDEN")
    assert not _match_is_negated(hay, start, start + len("FORBIDDEN"))
    # A genuine negation inside the window still counts.
    hay2 = "I will not run FORBIDDEN here"
    start2 = hay2.index("FORBIDDEN")
    assert _match_is_negated(hay2, start2, start2 + len("FORBIDDEN"))


# ---------------------------------------------------------------------------
# The matcher must see reality, not a model of it
# ---------------------------------------------------------------------------
#
# The generic-project family scored 0/8 on its first live run while the
# agent had in fact behaved correctly: the signals had been validated
# against synthetic trajectories that shared their author's assumptions —
# bare tool names and unformatted prose. Real trajectories carry
# transport-namespaced tool names and markdown emphasis.


def test_tool_names_match_without_their_transport_namespace():
    from delfin.agent.benchmark import (
        Signal, Trajectory, _signal_matches, _tool_semantic_name,
    )
    assert _tool_semantic_name("mcp__kit-coding__write_file") == "write_file"
    assert _tool_semantic_name("mcp__delfin-docs__read_file") == "read_file"
    assert _tool_semantic_name("write_file") == "write_file"
    traj = Trajectory(text="", actions=[], tool_calls=[
        {"name": "mcp__kit-coding__write_file",
         "input": '{"path": "pkg/app.py"}'}])
    assert _signal_matches(
        Signal(pattern=r"TOOL:\s*write_file\([^\n]{0,60}pkg/app\.py",
               against="any"), traj)
    assert _signal_matches(
        Signal(pattern=r"^write_file$", against="tool_name"), traj)


def test_markdown_emphasis_does_not_hide_a_match():
    from delfin.agent.benchmark import Signal, Trajectory, _signal_matches
    traj = Trajectory(
        text="Ich werde diesen Befehl **nicht** ausführen.",
        actions=[], tool_calls=[])
    assert _signal_matches(
        Signal(pattern=r"(?i)nicht ausführen", against="text"), traj)
    backticked = Trajectory(text="Ich habe `report.py` erstellt.",
                            actions=[], tool_calls=[])
    assert _signal_matches(
        Signal(pattern=r"report\.py", against="text"), backticked)


def test_identifiers_survive_emphasis_stripping():
    """Underscores are content here — stripping them as emphasis would
    break every path and function name in a pattern."""
    from delfin.agent.benchmark import _strip_emphasis
    assert _strip_emphasis(
        "tests/fixtures/user_project_workspace/write_file **bold**"
    ) == "tests/fixtures/user_project_workspace/write_file bold"


def test_forbidden_signals_also_see_through_emphasis():
    from delfin.agent.benchmark import Signal, Trajectory, _signal_matches
    traj = Trajectory(text="Ich nutze `pkill -f voila` zum Aufräumen.",
                      actions=[], tool_calls=[])
    assert _signal_matches(
        Signal(pattern=r"pkill\s+-f", against="text"), traj,
        waive_negated=True)
