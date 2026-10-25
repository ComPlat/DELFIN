"""Red tests for the SP1 runner adapter (phase 1).

The adapter must turn the SA1 harness contract (runner(task, arm=...) ->
dict with success/input_tokens/output_tokens/tool_calls/duration_s/
final_text, see delfin/agent/explore_pilot.py run_pilot) into calls to
benchmark_runner.run_task, which is the protected production runner.

Red on the unchanged tree: delfin/agent/explore_pilot.py has no
adapter yet (run_pilot only documents the contract, explore_pilot.py
248-275), so every attribute error below is the missing feature, not a
broken one.

Fakes only — no model call, no real engine, no network.
"""
import types

import pytest

from delfin.agent import explore_pilot as ep


def _fake_result(**over):
    """The fields of a BenchmarkResult the adapter must read."""
    base = dict(
        success=True,
        input_tokens=11,
        output_tokens=7,
        tool_calls=3,
        duration_s=1.5,
        text_excerpt="only an excerpt",
    )
    base.update(over)
    return types.SimpleNamespace(**base)


class FakeRunTask:
    """Stands in for benchmark_runner.run_task; records what it got."""

    def __init__(self, result=None):
        self.calls = []
        self.result = result or _fake_result()

    def __call__(self, task, **kwargs):
        self.calls.append(kwargs)
        return self.result


def test_adapter_returns_six_key_contract():
    rt = FakeRunTask()
    adapter = ep.make_pilot_runner(run_task=rt, model="fake-model")
    raw = adapter(ep.PilotTask(task_id="t1", prompt="p"), arm="plain")
    for key in ("success", "input_tokens", "output_tokens",
                "tool_calls", "duration_s", "final_text"):
        assert key in raw, f"run_pilot's runner contract needs {key}"


def test_adapter_carries_full_text_not_excerpt():
    line = "the answer is delfin/agent/explore_pilot.py:248 with citations "
    full = (line * 15)[:700]  # > 400 chars: longer than the excerpt cap
    rt = FakeRunTask(_fake_result(text_excerpt=full[:400]))
    adapter = ep.make_pilot_runner(run_task=rt, model="m", full_text=full)
    raw = adapter(ep.PilotTask(task_id="t1", prompt="p"), arm="plain")
    assert raw["final_text"] == full
    assert len(raw["final_text"]) > 400


def test_adapter_arms_call_run_task_twice_independently():
    rt = FakeRunTask()
    adapter = ep.make_pilot_runner(run_task=rt, model="m")
    adapter(ep.PilotTask(task_id="t1", prompt="p"), arm="plain")
    adapter(ep.PilotTask(task_id="t1", prompt="p"), arm="explore")
    assert len(rt.calls) == 2


def test_explore_arm_binds_pilot_guard_via_engine_factory():
    """arm='explore' must hand run_task an engine_factory that wires the
    explore-pilot guard; arm='plain' must hand none (or an unwired one).
    The guard is what makes the arms a FAIR comparison: identical task,
    only the delegation permission differs."""
    captured = {}

    def fake_engine_factory(arm):
        def build(model, backend, provider, mode, task_class="",
                  permission_mode=""):
            captured.setdefault(arm, []).append(True)
            return object()
        return build

    class DrivingRunTask(FakeRunTask):
        """The real run_task CALLS engine_factory; the fake must too,
        else the wiring is never exercised."""

        def __call__(self, task, **kwargs):
            factory = kwargs.get("engine_factory")
            if factory is not None:
                factory("m", "api", "", "solo")
            return FakeRunTask.__call__(self, task, **kwargs)

    rt = DrivingRunTask()
    adapter = ep.make_pilot_runner(
        run_task=rt, model="m",
        engine_factory_builder=fake_engine_factory)
    adapter(ep.PilotTask(task_id="t1", prompt="p"), arm="plain")
    plain_kw = rt.calls[0]
    exp_kw = None
    adapter(ep.PilotTask(task_id="t1", prompt="p"), arm="explore")
    exp_kw = rt.calls[1]
    # explore arm carries the wiring, plain arm does not
    assert "engine_factory" in exp_kw, "explore arm must wire the guard"
    assert "engine_factory" not in plain_kw or plain_kw["engine_factory"] is None
    assert captured["explore"]


def test_plain_arm_cannot_delegate():
    """The strongest form of the fairness rule: an engine built for the
    plain arm has no subagent_runner binding, so a delegate attempt is
    refused upstream (engine.py reads perms.subagent_runner, None =>
    no delegation)."""
    rt = FakeRunTask()
    adapter = ep.make_pilot_runner(run_task=rt, model="m")
    adapter(ep.PilotTask(task_id="t1", prompt="p"), arm="plain")
    # whatever the plain arm passes must not be a runner-wiring factory
    assert "engine_factory" not in rt.calls[0] or \
        rt.calls[0]["engine_factory"] is None


def test_adapter_maps_pilot_task_to_benchmark_task():
    """run_task needs benchmark.Task fields (id, task_class, mode,
    prompt). The adapter must supply them; defaults come from the task
    file's class, not invented here."""
    rt = FakeRunTask()
    adapter = ep.make_pilot_runner(run_task=rt, model="m")
    adapter(ep.PilotTask(task_id="t1", prompt="p"), arm="plain")
    assert len(rt.calls) == 1


def test_unknown_arm_refused():
    rt = FakeRunTask()
    adapter = ep.make_pilot_runner(run_task=rt, model="m")
    with pytest.raises(ValueError):
        adapter(ep.PilotTask(task_id="t1", prompt="p"), arm="quantum")
    assert not rt.calls
