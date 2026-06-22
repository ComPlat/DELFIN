"""Tests for the live-engine benchmark runner.

The runner glues the scorer to an ``AgentEngine``.  We don't spin up a
real engine in unit tests — we inject a fake engine_factory and a fake
``_run_once`` so we can verify routing, scoring, and error handling
deterministically.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import benchmark as bm
from delfin.agent import benchmark_runner as br


# ---------------------------------------------------------------------------
# extract_actions
# ---------------------------------------------------------------------------


def test_extract_actions_canonical_form():
    text = "Sure — switching now.\nACTION: /tab calc\n"
    assert br.extract_actions(text) == ["/tab calc"]


def test_extract_actions_no_space_after_colon():
    """ACTION:/tab calc — the tolerant form the dashboard accepts."""
    assert br.extract_actions("ACTION:/tab calc") == ["/tab calc"]


def test_extract_actions_lowercase_action_prefix():
    assert br.extract_actions("Action: /tab submit") == ["/tab submit"]


def test_extract_actions_bare_slash_with_known_prefix():
    """A bare slash line is accepted only for whitelisted dashboard prefixes."""
    text = "I'll switch tabs.\n/tab calc\n"
    assert "/tab calc" in br.extract_actions(text)


def test_extract_actions_rejects_bare_slash_unknown_prefix():
    """A bare ``/home/user/file.py`` line must NOT become an action."""
    text = "Open this: /home/user/file.py for reference"
    assert "/home/user/file.py" not in br.extract_actions(text)


def test_extract_actions_dedup():
    text = "ACTION: /tab calc\nand again /tab calc\n"
    out = br.extract_actions(text)
    assert out.count("/tab calc") == 1


def test_extract_actions_multiple_distinct():
    """ACTION must be at line-start (canonical), not mid-line — both
    extracted in order."""
    text = "Sure.\nACTION: /tab calc\nACTION: /control key value\n"
    out = br.extract_actions(text)
    assert "/tab calc" in out
    assert "/control key value" in out


def test_extract_actions_handles_empty():
    assert br.extract_actions("") == []
    assert br.extract_actions(None) == []  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# trajectory_from_run
# ---------------------------------------------------------------------------


def test_trajectory_from_run_extracts_actions_from_text():
    raw = {
        "text": "ACTION: /tab calc",
        "tool_calls": [], "input_tokens": 10, "output_tokens": 5, "error": "",
    }
    traj = br.trajectory_from_run(raw, duration_s=2.5)
    assert traj.actions == ["/tab calc"]
    assert traj.duration_s == 2.5
    assert traj.input_tokens == 10


def test_trajectory_from_run_preserves_tool_calls():
    raw = {
        "text": "looking",
        "tool_calls": [{"name": "Read", "input": {"file_path": "/x"}}],
        "input_tokens": 0, "output_tokens": 0, "error": "",
    }
    traj = br.trajectory_from_run(raw, duration_s=1.0, cost_usd=0.01)
    assert len(traj.tool_calls) == 1
    assert traj.tool_calls[0]["name"] == "Read"
    assert traj.cost_usd == pytest.approx(0.01)


# ---------------------------------------------------------------------------
# run_task — happy path + error path
# ---------------------------------------------------------------------------


class _FakeEngine:
    def __init__(self, cost_usd: float = 0.0):
        self.cost_usd = cost_usd


def test_run_task_happy_path_records_success():
    task = bm.Task(
        id="t1", task_class="dashboard_nav", mode="dashboard",
        prompt="öffne Calc",
        expected_signals=(
            bm.Signal(pattern=r"/tab\s+calc", against="action"),
        ),
        forbidden_signals=(),
        max_duration_s=10.0, max_cost_usd=0.05, max_tool_calls=2,
    )

    def fake_factory(model, backend, provider, mode):
        return _FakeEngine(cost_usd=0.0)

    def fake_run_once(engine, prompt, *, max_tokens=4096):
        engine.cost_usd = 0.01
        return {
            "text": "ACTION: /tab calc",
            "tool_calls": [],
            "input_tokens": 100, "output_tokens": 20, "error": "",
        }

    t = iter([100.0, 101.5, 101.6])

    result = br.run_task(
        task, model="kit.qwen3.5-397b-A17b",
        engine_factory=fake_factory, run_once=fake_run_once,
        clock=lambda: next(t),
    )
    assert result.success is True
    assert result.quality_0_100 >= 80
    assert result.model == "kit.qwen3.5-397b-A17b"
    assert result.duration_s == pytest.approx(1.5)
    assert result.cost_usd == pytest.approx(0.01)


def test_run_task_failure_when_required_signal_missing():
    task = bm.Task(
        id="t2", task_class="dashboard_nav", mode="dashboard",
        prompt="öffne Calc",
        expected_signals=(
            bm.Signal(pattern=r"/tab\s+calc", against="action"),
        ),
        forbidden_signals=(),
        max_duration_s=10.0, max_cost_usd=0.05, max_tool_calls=2,
    )

    def fake_factory(model, backend, provider, mode):
        return _FakeEngine()

    def fake_run_once(engine, prompt, *, max_tokens=4096):
        return {
            "text": "I'm not sure which tab you want.",
            "tool_calls": [], "input_tokens": 50, "output_tokens": 10,
            "error": "",
        }

    result = br.run_task(
        task, model="weak", engine_factory=fake_factory,
        run_once=fake_run_once,
    )
    assert result.success is False
    assert any(":expected" in m or "expected[" in m
               for m in result.missing_signals)


def test_run_task_engine_init_failure_yields_error_result():
    task = bm.Task(
        id="t3", task_class="misc", mode="solo", prompt="hi",
        expected_signals=(bm.Signal(pattern="ok"),),
        forbidden_signals=(),
        max_duration_s=10.0, max_cost_usd=0.05, max_tool_calls=2,
    )

    def bad_factory(*a, **kw):
        raise RuntimeError("no api key")

    def never_called_run_once(*a, **kw):
        raise AssertionError("must not be called when factory fails")

    result = br.run_task(
        task, model="x", engine_factory=bad_factory,
        run_once=never_called_run_once,
    )
    assert result.success is False
    assert "engine init failed" in result.error


def test_run_task_run_once_exception_captured():
    task = bm.Task(
        id="t4", task_class="misc", mode="solo", prompt="hi",
        expected_signals=(bm.Signal(pattern="ok"),),
        forbidden_signals=(),
        max_duration_s=10.0, max_cost_usd=0.05, max_tool_calls=2,
    )

    def fake_factory(*a, **kw):
        return _FakeEngine()

    def bad_run_once(engine, prompt, *, max_tokens=4096):
        raise RuntimeError("rate-limited")

    result = br.run_task(
        task, model="x", engine_factory=fake_factory,
        run_once=bad_run_once,
    )
    assert result.success is False
    assert "rate-limited" in result.error


# ---------------------------------------------------------------------------
# run_suite — iteration + progress callback
# ---------------------------------------------------------------------------


def test_run_suite_invokes_progress_for_each_task():
    tasks = [
        bm.Task(id=f"t{i}", task_class="misc", mode="solo",
                prompt="hi",
                expected_signals=(bm.Signal(pattern="ok"),),
                forbidden_signals=(),
                max_duration_s=10.0, max_cost_usd=0.05,
                max_tool_calls=2)
        for i in range(3)
    ]

    def fake_factory(*a, **kw):
        return _FakeEngine()

    def fake_run_once(engine, prompt, *, max_tokens=4096):
        return {"text": "ok", "tool_calls": [],
                "input_tokens": 1, "output_tokens": 1, "error": ""}

    seen = []

    def progress(task, result):
        seen.append((task.id, result.success))

    results = br.run_suite(
        tasks, model="x", engine_factory=fake_factory,
        run_once=fake_run_once, progress=progress,
    )
    assert len(results) == 3
    assert len(seen) == 3
    assert [s[0] for s in seen] == ["t0", "t1", "t2"]
    assert all(s[1] for s in seen)


def test_run_suite_progress_exception_does_not_abort():
    """A flaky progress callback must not stop the suite."""
    tasks = [
        bm.Task(id="t0", task_class="misc", mode="solo", prompt="hi",
                expected_signals=(bm.Signal(pattern="ok"),),
                forbidden_signals=(),
                max_duration_s=10.0, max_cost_usd=0.05,
                max_tool_calls=2),
    ]

    def fake_factory(*a, **kw):
        return _FakeEngine()

    def fake_run_once(engine, prompt, *, max_tokens=4096):
        return {"text": "ok", "tool_calls": [],
                "input_tokens": 0, "output_tokens": 0, "error": ""}

    def boom_progress(task, result):
        raise RuntimeError("printer broken")

    results = br.run_suite(
        tasks, model="x", engine_factory=fake_factory,
        run_once=fake_run_once, progress=boom_progress,
    )
    assert len(results) == 1
    assert results[0].success is True


# ---------------------------------------------------------------------------
# run_task — N=3 retry / replicates
# ---------------------------------------------------------------------------


def test_run_task_repeats_runs_engine_N_times():
    """``repeats=3`` must invoke the factory + run_once exactly 3×."""
    task = bm.Task(
        id="rep1", task_class="misc", mode="solo", prompt="hi",
        expected_signals=(bm.Signal(pattern="ok"),),
        forbidden_signals=(),
        max_duration_s=10.0, max_cost_usd=0.05, max_tool_calls=2,
    )
    factory_calls = []
    run_once_calls = []

    def factory(model, backend, provider, mode):
        factory_calls.append(model)
        return _FakeEngine()

    def fake_run_once(engine, prompt, *, max_tokens=4096):
        run_once_calls.append(prompt)
        return {"text": "ok", "tool_calls": [],
                "input_tokens": 0, "output_tokens": 0, "error": ""}

    result = br.run_task(
        task, model="X", engine_factory=factory,
        run_once=fake_run_once, repeats=3,
    )
    assert len(factory_calls) == 3
    assert len(run_once_calls) == 3
    assert result.n_samples == 3
    assert result.success is True
    assert result.success_rate == 1.0


def test_run_task_repeats_majority_success_when_one_fails():
    task = bm.Task(
        id="rep2", task_class="misc", mode="solo", prompt="hi",
        expected_signals=(bm.Signal(pattern="ok"),),
        forbidden_signals=(),
        max_duration_s=10.0, max_cost_usd=0.05, max_tool_calls=2,
    )

    def factory(*a, **kw):
        return _FakeEngine()

    # 2/3 produce "ok" (PASS), 1/3 produces "nope" (FAIL)
    counter = [0]

    def fake_run_once(engine, prompt, *, max_tokens=4096):
        i = counter[0]
        counter[0] += 1
        text = "ok" if i != 1 else "nope"
        return {"text": text, "tool_calls": [],
                "input_tokens": 0, "output_tokens": 0, "error": ""}

    result = br.run_task(
        task, model="X", engine_factory=factory,
        run_once=fake_run_once, repeats=3,
    )
    assert result.success is True       # 2-of-3 majority
    assert result.success_rate == pytest.approx(2 / 3)
    assert result.n_samples == 3


def test_run_task_repeats_on_replicate_callback_fires_per_sample():
    task = bm.Task(
        id="rep3", task_class="misc", mode="solo", prompt="hi",
        expected_signals=(bm.Signal(pattern="ok"),),
        forbidden_signals=(),
        max_duration_s=10.0, max_cost_usd=0.05, max_tool_calls=2,
    )

    def factory(*a, **kw):
        return _FakeEngine()

    def fake_run_once(*a, **kw):
        return {"text": "ok", "tool_calls": [],
                "input_tokens": 0, "output_tokens": 0, "error": ""}

    seen = []

    def cb(idx, result):
        seen.append((idx, result.success))

    br.run_task(
        task, model="X", engine_factory=factory,
        run_once=fake_run_once, repeats=3, on_replicate=cb,
    )
    assert [s[0] for s in seen] == [0, 1, 2]
    assert all(s[1] for s in seen)


def test_run_task_repeats_one_is_passthrough():
    """repeats=1 should not engage aggregation (no n_samples>1 fields)."""
    task = bm.Task(
        id="rep4", task_class="misc", mode="solo", prompt="hi",
        expected_signals=(bm.Signal(pattern="ok"),),
        forbidden_signals=(),
        max_duration_s=10.0, max_cost_usd=0.05, max_tool_calls=2,
    )

    def factory(*a, **kw):
        return _FakeEngine()

    def fake_run_once(*a, **kw):
        return {"text": "ok", "tool_calls": [],
                "input_tokens": 0, "output_tokens": 0, "error": ""}

    result = br.run_task(
        task, model="X", engine_factory=factory,
        run_once=fake_run_once, repeats=1,
    )
    assert result.n_samples == 1
    assert result.quality_stdev == 0.0


# ---------------------------------------------------------------------------
# resolve_profile_name
# ---------------------------------------------------------------------------


def test_resolve_profile_name_for_known_model():
    """Primary KIT model has a profile registered — must return its
    notes (≠ 'default')."""
    name = br.resolve_profile_name("kit.qwen3.5-397b-A17b")
    assert name and name != "default"


def test_resolve_profile_name_for_unknown_model_falls_back():
    """Truly unknown model — falls back to the WEAK_DEFAULT profile,
    whose notes string is non-empty.  Either way must not crash."""
    name = br.resolve_profile_name("completely-made-up-model-xyz-2099")
    assert isinstance(name, str)


# ---------------------------------------------------------------------------
# CLI smoke: argparse wires bench subcommand
# ---------------------------------------------------------------------------


def test_cli_bench_list_runs_without_engine(capsys):
    from delfin.agent import cli as agent_cli
    rc = agent_cli.main(["bench", "list"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "dash_nav_calc_typo" in out
    assert "solo_research_setworking" in out


def test_cli_bench_run_requires_model(capsys):
    from delfin.agent import cli as agent_cli
    # No --model → argparse rejects at parse time (rc=2)
    with pytest.raises(SystemExit) as exc:
        agent_cli.main(["bench", "run"])
    assert exc.value.code == 2


def test_cli_bench_compare_missing_files(capsys, tmp_path):
    from delfin.agent import cli as agent_cli
    missing = tmp_path / "nope.jsonl"
    rc = agent_cli.main(["bench", "compare", str(missing), str(missing)])
    assert rc == 2


def test_cli_bench_compare_happy_path(capsys, tmp_path):
    from delfin.agent import cli as agent_cli
    base = tmp_path / "base.jsonl"
    cand = tmp_path / "cand.jsonl"
    base.write_text(
        '{"task_id": "a", "model": "m", "success": true, "quality_0_100": 60, "cost_usd": 0.10, "duration_s": 5.0, "tool_calls": 1}\n'
        '{"task_id": "b", "model": "m", "success": true, "quality_0_100": 50, "cost_usd": 0.20, "duration_s": 6.0, "tool_calls": 1}\n'
        '{"task_id": "c", "model": "m", "success": true, "quality_0_100": 70, "cost_usd": 0.30, "duration_s": 7.0, "tool_calls": 1}\n'
        '{"task_id": "d", "model": "m", "success": true, "quality_0_100": 80, "cost_usd": 0.40, "duration_s": 8.0, "tool_calls": 1}\n',
        encoding="utf-8",
    )
    cand.write_text(
        '{"task_id": "a", "model": "m", "success": true, "quality_0_100": 80, "cost_usd": 0.05, "duration_s": 3.0, "tool_calls": 1}\n'
        '{"task_id": "b", "model": "m", "success": true, "quality_0_100": 70, "cost_usd": 0.10, "duration_s": 4.0, "tool_calls": 1}\n'
        '{"task_id": "c", "model": "m", "success": true, "quality_0_100": 85, "cost_usd": 0.15, "duration_s": 5.0, "tool_calls": 1}\n'
        '{"task_id": "d", "model": "m", "success": true, "quality_0_100": 90, "cost_usd": 0.20, "duration_s": 6.0, "tool_calls": 1}\n',
        encoding="utf-8",
    )
    rc = agent_cli.main(["bench", "compare", str(base), str(cand)])
    out = capsys.readouterr().out
    assert rc == 0
    assert "BETTER" in out.upper()
