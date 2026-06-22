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
