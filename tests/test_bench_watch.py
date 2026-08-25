"""Tests for the unattended benchmark-watch pipeline (bench_watch).

Everything runs OFFLINE: the live suite run, the suspect rerun and the
attention emission are injected/monkeypatched, and the real run-file
shape is covered by two committed fixture runs (complete 48-task runs
of the same model taken on the same day — the pair that exhibits the
measured pass/fail flip noise the thresholds are designed around).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import bench_watch as bw
from delfin.agent import cli as agent_cli
from delfin.agent import scheduler_daemon as sd
from delfin.agent.scheduler import ScheduleEntry


FIXTURE_DIR = Path(__file__).parent / "fixtures" / "bench_runs"
REAL_MODEL = "kit.qwen3.5-397b-A17b"
REAL_OLD = FIXTURE_DIR / f"1785305351_{REAL_MODEL}.jsonl"
REAL_NEW = FIXTURE_DIR / f"1785306861_{REAL_MODEL}.jsonl"

MODEL = "test.model-x"


# ---------------------------------------------------------------------------
# Helpers: synthetic run records / files with the real JSONL shape
# ---------------------------------------------------------------------------


def _record(task_id, *, success=True, quality=75, cost=0.10,
            n_samples=1, per_run_quality=None, per_run_success=None,
            error=""):
    return {
        "task_id": task_id,
        "task_class": "synthetic",
        "model": MODEL,
        "mode": "solo",
        "ts": 1785300000.0,
        "success": bool(success),
        "quality_0_100": int(quality),
        "duration_s": 2.0,
        "cost_usd": float(cost),
        "n_samples": int(n_samples),
        "quality_stdev": 0.0,
        "success_rate": 1.0 if success else 0.0,
        "per_run_quality": per_run_quality or [],
        "per_run_success": per_run_success or [],
        "matched_signals": [],
        "violated_signals": [],
        "missing_signals": [],
        "error": error,
    }


def _write_run(dirpath: Path, ts: int, records, model: str = MODEL) -> Path:
    dirpath.mkdir(parents=True, exist_ok=True)
    path = dirpath / f"{ts}_{model}.jsonl"
    with path.open("w", encoding="utf-8") as fh:
        for rec in records:
            fh.write(json.dumps(rec) + "\n")
    return path


def _baseline(median=80.0, spread=0.0, majority=True, n=5, cost=0.10):
    return {
        "n_files": n, "n_samples": n,
        "majority_success": majority,
        "pass_fraction": 1.0 if majority else 0.0,
        "median_quality": float(median),
        "quality_spread": float(spread),
        "median_cost_usd": float(cost),
    }


# ---------------------------------------------------------------------------
# load_history — rolling baseline math
# ---------------------------------------------------------------------------


def test_list_model_runs_orders_and_filters(tmp_path):
    _write_run(tmp_path, 100, [_record("a")])
    _write_run(tmp_path, 200, [_record("a")])
    _write_run(tmp_path, 150, [_record("a")], model="other.model")
    # Custom-run_id files (non-numeric prefix) must not enter the window.
    (tmp_path / f"behavior_baseline_300_{MODEL}.jsonl").write_text(
        json.dumps(_record("a")) + "\n")
    files = bw.list_model_runs(MODEL, runs_dir=tmp_path)
    assert [f.name for f in files] == [
        f"200_{MODEL}.jsonl", f"100_{MODEL}.jsonl"]


def test_load_history_majority_median_spread(tmp_path):
    qualities = [80, 90, 70, 85, 75]
    passes = [True, True, False, True, True]
    for i, (q, s) in enumerate(zip(qualities, passes)):
        _write_run(tmp_path, 100 + i, [_record("t1", quality=q, success=s)])
    hist = bw.load_history(MODEL, runs_dir=tmp_path, last_k=5)
    base = hist["tasks"]["t1"]
    assert base["n_samples"] == 5
    assert base["majority_success"] is True
    assert base["pass_fraction"] == pytest.approx(0.8)
    assert base["median_quality"] == 80.0
    assert base["quality_spread"] == pytest.approx(7.9057, abs=1e-3)


def test_load_history_window_is_per_task_last_k(tmp_path):
    # 7 files for t1: the two OLDEST have quality 0 and must fall outside
    # a last_k=5 window.
    for i in range(7):
        q = 0 if i < 2 else 80
        _write_run(tmp_path, 100 + i, [_record("t1", quality=q)])
    hist = bw.load_history(MODEL, runs_dir=tmp_path, last_k=5)
    base = hist["tasks"]["t1"]
    assert base["n_files"] == 5
    assert base["median_quality"] == 80.0
    assert base["quality_spread"] == 0.0


def test_load_history_tolerates_torn_lines_and_missing_tasks(tmp_path):
    _write_run(tmp_path, 100, [_record("t1", quality=70), _record("t2")])
    p = _write_run(tmp_path, 101, [_record("t1", quality=90)])
    with p.open("a", encoding="utf-8") as fh:
        fh.write('{"task_id": "t2", "quality_0_100": TORN\n')
    hist = bw.load_history(MODEL, runs_dir=tmp_path, last_k=5)
    assert hist["tasks"]["t1"]["n_samples"] == 2
    assert hist["tasks"]["t1"]["median_quality"] == 80.0
    # t2 only parsed from the intact file
    assert hist["tasks"]["t2"]["n_samples"] == 1


def test_load_history_expands_aggregated_replicates(tmp_path):
    _write_run(tmp_path, 100, [_record(
        "t1", quality=70, cost=0.30, n_samples=3,
        per_run_quality=[60, 70, 80],
        per_run_success=[True, True, False])])
    base = bw.load_history(MODEL, runs_dir=tmp_path)["tasks"]["t1"]
    assert base["n_samples"] == 3
    assert base["median_quality"] == 70.0
    assert base["pass_fraction"] == pytest.approx(2 / 3)
    # Aggregate cost is the SUM over samples — split back per sample.
    assert base["median_cost_usd"] == pytest.approx(0.10)


def test_load_history_exclude_drops_one_file(tmp_path):
    _write_run(tmp_path, 100, [_record("t1", quality=50)])
    newest = _write_run(tmp_path, 200, [_record("t1", quality=90)])
    hist = bw.load_history(MODEL, runs_dir=tmp_path, exclude=newest)
    assert hist["tasks"]["t1"]["median_quality"] == 50.0
    assert str(newest) not in hist["files"]


def test_load_history_missing_dir_is_empty(tmp_path):
    hist = bw.load_history(MODEL, runs_dir=tmp_path / "nope")
    assert hist["tasks"] == {}


# ---------------------------------------------------------------------------
# suspect_threshold + compare_run classification
# ---------------------------------------------------------------------------


def test_suspect_threshold_floor_and_scaling():
    assert bw.suspect_threshold(0.0) == 8.0
    assert bw.suspect_threshold(4.0) == 8.0     # 1.5*4=6 < floor
    assert bw.suspect_threshold(20.0) == 30.0
    assert bw.suspect_threshold(-3.0) == 8.0    # garbage clamps to floor


def test_compare_run_classifies_all_statuses(tmp_path):
    history = {"tasks": {
        "base_pass": _baseline(median=80, spread=0, majority=True),
        "base_fail": _baseline(median=30, spread=0, majority=False),
        "noisy": _baseline(median=60, spread=20, majority=True),
    }}
    run = _write_run(tmp_path, 300, [
        _record("base_pass", success=False, quality=40),   # fail flip
        _record("base_fail", success=True, quality=50),    # pass flip
        _record("noisy", success=True, quality=35),        # within 30-band
        _record("brand_new", success=True, quality=90),    # no baseline
    ])
    cmp = bw.compare_run(run, history=history)
    by_id = {e["task_id"]: e for e in cmp["per_task"]}
    assert by_id["base_pass"]["status"] == "suspect_regression"
    assert by_id["base_pass"]["reason"] == "fail_flip"
    assert by_id["base_fail"]["status"] == "improved"
    assert by_id["base_fail"]["reason"] == "pass_flip"
    assert by_id["noisy"]["status"] == "stable"
    assert by_id["noisy"]["threshold"] == 30.0
    assert by_id["brand_new"]["status"] == "new_task"
    assert cmp["counts"]["suspect_regression"] == 1
    assert cmp["counts"]["new_task"] == 1


def test_compare_run_quality_drop_beyond_threshold(tmp_path):
    history = {"tasks": {"t1": _baseline(median=80, spread=0)}}
    run = _write_run(tmp_path, 300, [
        _record("t1", success=True, quality=68)])   # -12 <= -8
    cmp = bw.compare_run(run, history=history)
    assert cmp["per_task"][0]["status"] == "suspect_regression"
    assert cmp["per_task"][0]["reason"] == "quality_drop"


def test_compare_run_small_jitter_is_stable(tmp_path):
    """The +-7 point scoring jitter between identical runs must never
    become a suspect — the whole point of the 8-point floor."""
    history = {"tasks": {"t1": _baseline(median=80, spread=0)}}
    run = _write_run(tmp_path, 300, [
        _record("t1", success=True, quality=73)])   # -7 > -8
    assert bw.compare_run(run, history=history)["per_task"][0][
        "status"] == "stable"


def test_compare_run_missing_file_is_empty_not_crash(tmp_path):
    cmp = bw.compare_run(tmp_path / "absent.jsonl", history={"tasks": {}})
    assert cmp["n_tasks"] == 0
    assert cmp["suspects"] == []


# ---------------------------------------------------------------------------
# The measured-noise scenario: 5 pass/fail flips between identical runs
# ---------------------------------------------------------------------------


def _five_flip_setup(tmp_path):
    """5 baseline runs of 48 passing tasks; current run flips 5 to FAIL."""
    task_ids = [f"task_{i:02d}" for i in range(48)]
    for k in range(5):
        _write_run(tmp_path, 100 + k,
                   [_record(t, success=True, quality=75 + (k % 3))
                    for t in task_ids])
    flipped = task_ids[:5]
    current = [
        _record(t, success=(t not in flipped),
                quality=15 if t in flipped else 76)
        for t in task_ids
    ]
    run = _write_run(tmp_path, 200, current)
    return task_ids, flipped, run


def test_five_flip_scenario_flags_exactly_the_flips(tmp_path):
    _, flipped, run = _five_flip_setup(tmp_path)
    hist = bw.load_history(MODEL, runs_dir=tmp_path, exclude=run)
    cmp = bw.compare_run(run, history=hist)
    assert sorted(e["task_id"] for e in cmp["suspects"]) == sorted(flipped)
    assert all(e["reason"] == "fail_flip" for e in cmp["suspects"])
    assert cmp["counts"]["stable"] == 43


def test_five_flip_scenario_recheck_downgrades_to_noise(tmp_path):
    """A passing recheck must classify all 5 flips as noise — the exact
    cry-wolf case observed between two identical real runs."""
    _, flipped, run = _five_flip_setup(tmp_path)
    hist = bw.load_history(MODEL, runs_dir=tmp_path, exclude=run)
    suspects = bw.compare_run(run, history=hist)["suspects"]
    rerun_calls = []

    def _rerun(task_id):
        rerun_calls.append(task_id)
        return _record(task_id, success=True, quality=75, cost=0.1,
                       n_samples=3)

    out = bw.recheck_suspects(
        suspects, model=MODEL, repeats=3, max_tasks=6, max_cost_usd=3.0,
        rerun=_rerun)
    assert sorted(rerun_calls) == sorted(flipped)
    assert out["confirmed"] == []
    assert len(out["noise"]) == 5
    assert out["skipped"] == []


def test_recheck_confirms_persistent_failure(tmp_path):
    _, flipped, run = _five_flip_setup(tmp_path)
    hist = bw.load_history(MODEL, runs_dir=tmp_path, exclude=run)
    suspects = bw.compare_run(run, history=hist)["suspects"]

    def _rerun(task_id):
        # First flipped task keeps failing; the rest recover.
        broken = task_id == flipped[0]
        return _record(task_id, success=not broken,
                       quality=10 if broken else 75, n_samples=3)

    out = bw.recheck_suspects(suspects, model=MODEL, rerun=_rerun)
    assert [c["task_id"] for c in out["confirmed"]] == [flipped[0]]
    assert out["confirmed"][0]["verdict"] == "confirmed_regression"
    assert len(out["noise"]) == 4


# ---------------------------------------------------------------------------
# recheck_suspects — caps, ordering, error isolation
# ---------------------------------------------------------------------------


def _suspect(task_id, *, reason="fail_flip", delta=-40.0, median=80.0,
             cost=0.10, threshold=8.0):
    return {
        "task_id": task_id, "status": "suspect_regression",
        "reason": reason, "success": False, "quality": median + delta,
        "cost_usd": cost, "threshold": threshold, "delta_quality": delta,
        "baseline": _baseline(median=median, cost=cost),
    }


def test_recheck_task_cap_skips_least_severe(tmp_path):
    suspects = (
        [_suspect(f"drop_{i}", reason="quality_drop", delta=-10 - i)
         for i in range(4)]
        + [_suspect(f"flip_{i}") for i in range(4)]
    )
    seen = []

    def _rerun(task_id):
        seen.append(task_id)
        return _record(task_id, success=True, quality=80)

    out = bw.recheck_suspects(
        suspects, model=MODEL, max_tasks=5, max_cost_usd=100.0,
        rerun=_rerun)
    # All 4 fail flips first, then the single largest drop.
    assert set(seen[:4]) == {f"flip_{i}" for i in range(4)}
    assert seen[4] == "drop_3"
    assert len(out["skipped"]) == 3
    assert all("task cap" in s["reason"] for s in out["skipped"])
    assert {s["task_id"] for s in out["skipped"]} == {
        "drop_0", "drop_1", "drop_2"}


def test_recheck_cost_cap_pre_estimate(tmp_path):
    # Each task estimates 3 * $1.00 — only one fits a $3 budget.
    suspects = [_suspect("a", cost=1.0), _suspect("b", cost=1.0)]
    out = bw.recheck_suspects(
        suspects, model=MODEL, repeats=3, max_cost_usd=3.0,
        rerun=lambda tid: _record(tid, success=True, quality=80, cost=3.0))
    assert len(out["checked"]) == 1
    assert len(out["skipped"]) == 1
    assert "cost cap" in out["skipped"][0]["reason"]


def test_recheck_cost_cap_on_actual_spend(tmp_path):
    # Estimates are tiny but actual reruns cost $1.60 each: after two
    # tasks the live spend (3.20) exceeds the $3 budget — task 3 skipped.
    suspects = [_suspect(t, cost=0.01) for t in ("a", "b", "c")]
    out = bw.recheck_suspects(
        suspects, model=MODEL, repeats=3, max_cost_usd=3.0,
        rerun=lambda tid: _record(tid, success=True, quality=80, cost=1.6))
    assert len(out["checked"]) == 2
    assert out["spent_usd"] == pytest.approx(3.2)
    assert len(out["skipped"]) == 1
    assert "mid-recheck" in out["skipped"][0]["reason"]


def test_recheck_rerun_error_is_isolated(tmp_path):
    suspects = [_suspect("boom"), _suspect("fine")]

    def _rerun(task_id):
        if task_id == "boom":
            raise RuntimeError("engine exploded")
        return _record(task_id, success=True, quality=80)

    out = bw.recheck_suspects(suspects, model=MODEL, rerun=_rerun)
    verdicts = {c["task_id"]: c["verdict"] for c in out["checked"]}
    assert verdicts["boom"] == "recheck_error"
    assert verdicts["fine"] == "noise"
    assert out["errors"]


def test_recheck_empty_suspects_never_touches_runner():
    out = bw.recheck_suspects(
        [], model=MODEL,
        rerun=lambda tid: pytest.fail("must not be called"))
    assert out["checked"] == []
    assert out["skipped"] == []


# ---------------------------------------------------------------------------
# nightly — full offline cycle
# ---------------------------------------------------------------------------


def _offline_nightly(tmp_path, *, rerun_success, recheck=True, emit=None,
                     suite_exc=None):
    """Run nightly with all live stages injected. Baseline: 5 clean runs;
    the 'new' run flips task_00 to FAIL."""
    for k in range(5):
        _write_run(tmp_path, 100 + k,
                   [_record(f"task_{i:02d}", success=True, quality=75)
                    for i in range(6)])

    def _suite(model, provider, backend, task_ids, repeats, runs_dir):
        if suite_exc is not None:
            raise suite_exc
        recs = [_record(f"task_{i:02d}", success=(i != 0),
                        quality=12 if i == 0 else 76)
                for i in range(6)]
        return _write_run(Path(runs_dir), 999, recs)

    def _rerun(task_id):
        return _record(task_id, success=rerun_success,
                       quality=75 if rerun_success else 10, n_samples=3)

    return bw.nightly(
        MODEL, "kit", "api",
        recheck=recheck, runs_dir=tmp_path,
        run_suite_fn=_suite, rerun=_rerun, emit=emit,
    )


def test_nightly_confirmed_regression_emits_one_attention(tmp_path):
    events = []

    def _emit(kind, **kwargs):
        events.append((kind, kwargs))
        return "att-test-1"

    summary = _offline_nightly(tmp_path, rerun_success=False, emit=_emit)
    assert summary["ok"] is True
    assert [c["task_id"] for c in summary["confirmed"]] == ["task_00"]
    assert summary["attention_id"] == "att-test-1"
    assert len(events) == 1
    kind, kwargs = events[0]
    assert kind == "run_failed"
    assert "task_00" in kwargs["detail"]
    assert MODEL in kwargs["title"]
    # Report written next to the runs.
    report = Path(summary["report_path"])
    assert report.parent == tmp_path / "reports"
    text = report.read_text(encoding="utf-8")
    assert "task_00" in text
    assert "CONFIRMED" in text


def test_nightly_noise_emits_nothing(tmp_path):
    events = []
    summary = _offline_nightly(
        tmp_path, rerun_success=True,
        emit=lambda kind, **kw: events.append(kind) or "x")
    assert summary["confirmed"] == []
    assert events == []
    assert summary["attention_id"] == ""


def test_nightly_no_recheck_leaves_suspects_unconfirmed(tmp_path):
    events = []
    summary = _offline_nightly(
        tmp_path, rerun_success=False, recheck=False,
        emit=lambda kind, **kw: events.append(kind) or "x")
    assert summary["recheck"] is None
    assert summary["confirmed"] == []
    assert events == []
    assert (summary["comparison"]["counts"]["suspect_regression"] == 1)
    assert "UNCONFIRMED" in Path(summary["report_path"]).read_text(
        encoding="utf-8")


def test_nightly_never_raises_on_suite_failure(tmp_path):
    summary = _offline_nightly(
        tmp_path, rerun_success=True,
        suite_exc=RuntimeError("no API key"),
        emit=lambda *a, **kw: pytest.fail("no attention on failed suite"))
    assert summary["ok"] is False
    assert any("suite run failed" in e for e in summary["errors"])
    # Report still written, with the error inside.
    assert "no API key" in Path(summary["report_path"]).read_text(
        encoding="utf-8")


def test_nightly_uses_attention_module_when_emit_not_injected(
        tmp_path, monkeypatch):
    from delfin.agent import attention as att
    calls = []
    monkeypatch.setattr(
        att, "emit_attention",
        lambda kind, **kw: calls.append((kind, kw)) or "att-mp-1")
    summary = _offline_nightly(tmp_path, rerun_success=False, emit=None)
    assert summary["attention_id"] == "att-mp-1"
    assert calls and calls[0][0] == "run_failed"


def test_format_report_sections(tmp_path):
    summary = _offline_nightly(tmp_path, rerun_success=False)
    text = bw.format_report(summary)
    assert text.startswith(f"# Benchmark watch — {MODEL}")
    assert "## Run vs rolling baseline" in text
    assert "### Suspects" in text
    assert "## Recheck" in text
    assert "## Verdict" in text
    assert "confirmed_regression" in text


def test_format_report_tolerates_empty_summary():
    assert "Benchmark watch" in bw.format_report({})


# ---------------------------------------------------------------------------
# Real fixture runs — the measured noise, end to end
# ---------------------------------------------------------------------------


def test_real_fixture_files_are_complete_runs():
    assert REAL_OLD.exists() and REAL_NEW.exists()
    for p in (REAL_OLD, REAL_NEW):
        lines = [ln for ln in p.read_text().splitlines() if ln.strip()]
        assert len(lines) == 48


def test_real_fixture_history_shape():
    hist = bw.load_history(REAL_MODEL, runs_dir=FIXTURE_DIR, last_k=5)
    assert len(hist["tasks"]) == 48
    assert all(b["n_samples"] == 2 for b in hist["tasks"].values())


def test_real_fixture_compare_known_outcome():
    """Baseline = the older run only; classify the newer identical-config
    run. Known measured outcome: 5 suspects (2 fail flips + 3 quality
    drops), 5 improved, 38 stable — and no crash anywhere."""
    hist = bw.load_history(
        REAL_MODEL, runs_dir=FIXTURE_DIR, exclude=REAL_NEW)
    cmp = bw.compare_run(REAL_NEW, history=hist)
    assert cmp["n_tasks"] == 48
    assert sorted(e["task_id"] for e in cmp["suspects"]) == [
        "dash_nav_submit",
        "dash_orca_set_functional",
        "fact_orca_casscf_keywords",
        "plan_chemistry_workflow",
        "workflow_parallel_checks",
    ]
    reasons = {e["task_id"]: e["reason"] for e in cmp["suspects"]}
    assert reasons["dash_nav_submit"] == "fail_flip"
    assert reasons["dash_orca_set_functional"] == "fail_flip"
    assert cmp["counts"]["improved"] == 5
    assert cmp["counts"]["stable"] == 38
    assert cmp["counts"]["new_task"] == 0


def test_real_fixture_suspects_downgrade_with_passing_recheck():
    hist = bw.load_history(
        REAL_MODEL, runs_dir=FIXTURE_DIR, exclude=REAL_NEW)
    suspects = bw.compare_run(REAL_NEW, history=hist)["suspects"]
    out = bw.recheck_suspects(
        suspects, model=REAL_MODEL,
        rerun=lambda tid: _record(tid, success=True, quality=80,
                                  cost=0.15, n_samples=3))
    assert out["confirmed"] == []
    assert len(out["noise"]) == 5


# ---------------------------------------------------------------------------
# CLI — compare-only path, nightly wiring, opt-in scheduling hint
# ---------------------------------------------------------------------------


def test_cli_bench_compare_history_mode_is_offline(capsys):
    rc = agent_cli.main([
        "bench", "compare", "--model", REAL_MODEL,
        "--runs-dir", str(FIXTURE_DIR),
    ])
    out = capsys.readouterr().out
    assert rc == 0
    assert "SUSPECT dash_nav_submit" in out
    assert "no API cost" in out


def test_cli_bench_compare_two_file_mode_unchanged(capsys):
    rc = agent_cli.main([
        "bench", "compare", str(REAL_OLD), str(REAL_NEW)])
    assert rc == 0
    assert "Verdict:" in capsys.readouterr().out


def test_cli_bench_compare_without_args_errors(capsys):
    rc = agent_cli.main(["bench", "compare"])
    assert rc == 2
    assert "--model" in capsys.readouterr().err


def test_cli_bench_nightly_prints_optin_scheduling_hint(
        capsys, monkeypatch):
    summary = {
        "ok": True, "model": REAL_MODEL, "comparison": {
            "n_tasks": 48, "counts": {"stable": 48, "improved": 0,
                                      "suspect_regression": 0,
                                      "new_task": 0}},
        "recheck": None, "confirmed": [], "attention_id": "",
        "errors": [], "report_path": "/tmp/report.md",
    }
    seen = {}

    def _fake_nightly(model, provider, backend, **kwargs):
        seen.update(model=model, provider=provider, backend=backend,
                    **kwargs)
        return summary

    monkeypatch.setattr(bw, "nightly", _fake_nightly)
    rc = agent_cli.main([
        "bench", "nightly", "--model", REAL_MODEL, "--provider", "kit",
        "--repeats", "2", "--last-k", "3", "--no-recheck",
    ])
    out = capsys.readouterr().out
    assert rc == 0
    assert seen["model"] == REAL_MODEL
    assert seen["repeats"] == 2
    assert seen["recheck"] is False
    assert seen["last_k"] == 3
    # The exact opt-in scheduling instruction + cost estimate.
    assert ("scheduler add-bench --model " + REAL_MODEL) in out
    assert "delfin-agent scheduler start" in out
    assert "~$8 per nightly run" in out
    assert "opt-in" in out


def test_cli_scheduler_add_bench_is_explicit_optin(capsys, monkeypatch):
    from delfin.agent import scheduler as sched_mod
    created = {}

    class _FakeScheduler:
        def schedule_interval(self, **kwargs):
            created.update(kwargs)
            return ScheduleEntry(
                id="fake123", kind="interval",
                prompt=kwargs["prompt"],
                every_seconds=kwargs["every_seconds"],
                workspace=kwargs.get("workspace", ""))

    monkeypatch.setattr(sched_mod, "Scheduler", _FakeScheduler)
    rc = agent_cli.main([
        "scheduler", "add-bench", "--model", REAL_MODEL,
        "--provider", "kit", "--every", "24h",
    ])
    out = capsys.readouterr().out
    assert rc == 0
    assert created["every_seconds"] == 86400
    assert created["prompt"] == (
        f"[bench-nightly] model={REAL_MODEL} provider=kit")
    assert "~$8 per run" in out


def test_cli_scheduler_add_bench_rejects_bad_interval(capsys):
    rc = agent_cli.main([
        "scheduler", "add-bench", "--model", "m", "--every", "soon"])
    assert rc == 2
    assert "interval" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# scheduler_daemon bench hook
# ---------------------------------------------------------------------------


def test_bench_entry_prompt_roundtrip():
    prompt = sd.format_bench_entry_prompt(
        model=REAL_MODEL, provider="kit", backend="api", repeats=2)
    cfg = sd.parse_bench_entry(prompt)
    assert cfg == {"model": REAL_MODEL, "provider": "kit",
                   "backend": "api", "repeats": 2}
    assert sd.parse_bench_entry("summarise the nightly logs") is None
    assert sd.parse_bench_entry("") is None


def test_fire_callback_routes_bench_entry(tmp_path, monkeypatch):
    calls = {}

    def _fake_nightly(model, provider, backend, **kwargs):
        import os
        calls.update(model=model, provider=provider, backend=backend,
                     cwd=os.getcwd())
        return {"ok": True, "confirmed": [], "errors": [],
                "report_path": str(tmp_path / "r.md")}

    monkeypatch.setattr(bw, "nightly", _fake_nightly)
    from delfin.agent import attention as att
    events = []
    monkeypatch.setattr(
        att, "emit_attention",
        lambda kind, **kw: events.append((kind, kw)) or "att-x")
    # run_entry (the LLM path) must never be touched for a bench entry.
    monkeypatch.setattr(
        sd, "run_entry",
        lambda *a, **kw: pytest.fail("bench entry took the LLM path"))

    ws = tmp_path / "workspace"
    ws.mkdir()
    entry = ScheduleEntry(
        id="e1", kind="interval",
        prompt=sd.format_bench_entry_prompt(model=REAL_MODEL,
                                            provider="kit"),
        every_seconds=86400, workspace=str(ws))
    fire = sd.make_fire_callback(log=lambda s: None)
    fire(entry)
    assert calls["model"] == REAL_MODEL
    assert Path(calls["cwd"]).resolve() == ws.resolve()
    assert [k for k, _ in events] == ["run_finished"]


def test_run_bench_entry_degraded_cycle_raises(tmp_path, monkeypatch):
    monkeypatch.setattr(
        bw, "nightly",
        lambda *a, **kw: {"ok": False, "errors": ["suite run failed: x"],
                          "confirmed": []})
    entry = ScheduleEntry(
        id="e2", kind="interval",
        prompt=sd.format_bench_entry_prompt(model=REAL_MODEL),
        every_seconds=86400, workspace=str(tmp_path))
    with pytest.raises(RuntimeError):
        sd.run_bench_entry(entry)


def test_run_bench_entry_without_model_disables():
    from delfin.agent.scheduler import DisableEntry
    entry = ScheduleEntry(
        id="e3", kind="interval", prompt="[bench-nightly]",
        every_seconds=86400, workspace="/tmp")
    with pytest.raises(DisableEntry):
        sd.run_bench_entry(entry)
