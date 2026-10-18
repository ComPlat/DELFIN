"""Red tests for the operator's three findings on the KG4 steer hook —
now against the REAL delfin.agent.benchmark / benchmark_runner, which
carry the built hook (operator build fb4157c0). Same assertions the px
copies were proven with; the px dependency is gone so the file passes in
any checkout and in CI.
"""
from __future__ import annotations

import json
from pathlib import Path


def _fake_engine_factory():
    class _FakeEngine:
        token_usage = {"input": 1, "output": 1}
        cost_usd = 0.0
        messages = []

    return lambda *a, **k: _FakeEngine()


_RUN_ONCE_OK = lambda engine, prompt, **k: {  # noqa: E731
    "text": "ok",
    "tool_calls": [{"name": "write_file"}],
    "input_tokens": 1, "output_tokens": 1, "error": ""}


def test_no_mid_run_messages_leaves_no_facts_file(tmp_path, monkeypatch):
    """Operator finding 2: a task WITHOUT mid_run_messages must leave no
    .longrun_facts.json in its workspace — the first draft wrote one into
    every task's workspace, which flips accepts that check nothing else
    was written."""
    from delfin.agent import benchmark_runner as br
    from delfin.agent.benchmark import Task
    task = Task(id="plain_task", task_class="office_ops", mode="solo",
                prompt="plain")
    monkeypatch.chdir(tmp_path)
    result = br._run_task_once(
        task, model="test-model", backend="test", provider="test",
        profile_name="test", engine_factory=_fake_engine_factory(),
        max_tokens=128, run_once=_RUN_ONCE_OK, clock=lambda: 0.0)
    found = list(tmp_path.rglob(".longrun_facts.json"))
    assert found == [], found


def test_intervention_count_survives_the_late_overwrite(tmp_path,
                                                        monkeypatch):
    """Operator finding 1: the existing raw['user_interventions'] =
    len(_denied) line used to erase the per-steer count. The built hook
    adds the steers AFTER that line; the count must survive."""
    from delfin.agent import benchmark_runner as br
    from delfin.agent.benchmark import Task
    task = Task(id="lr_steered_goal", task_class="office",
                mode="solo", prompt="old goal",
                mid_run_messages=(("after",
                                   "PERSON (mid-task steering, new goal): "
                                   "follow the new goal"),))
    monkeypatch.chdir(tmp_path)
    result = br._run_task_once(
        task, model="test-model", backend="test", provider="test",
        profile_name="test", engine_factory=_fake_engine_factory(),
        max_tokens=128, run_once=_RUN_ONCE_OK, clock=lambda: 0.0)
    assert result.user_interventions == 1, result.user_interventions


def test_facts_file_minimal_shape(tmp_path, monkeypatch):
    """Operator finding 3 + what the steer-accept judges: the file carries
    the steer text and tool call NAMES — no transcript text (the file
    sits in a workspace a later agent can read). Observed at the accept
    call site via a spy on run_acceptance, which the real
    _run_task_once calls INSIDE the pristine-workspace guard."""
    from delfin.agent import benchmark_runner as br
    from delfin.agent.benchmark import Task
    task = Task(id="lr_steered_goal", task_class="office",
                mode="solo", prompt="old goal", verify="never_called.py",
                mid_run_messages=(("after",
                                   "PERSON (mid-task steering, new goal): "
                                   "follow the new goal"),))
    monkeypatch.chdir(tmp_path)
    (tmp_path / "tests" / "fixtures" / "office_workspace").mkdir(
        parents=True, exist_ok=True)
    captured: dict = {}

    def spy_acceptance(name, workspace, **kw):
        p = Path(workspace) / ".longrun_facts.json"
        captured["exists"] = p.is_file()
        if p.is_file():
            captured["facts"] = json.loads(p.read_text(encoding="utf-8"))
        return None, "spy: not a real accept"

    monkeypatch.setattr(br, "run_acceptance", spy_acceptance)
    result = br._run_task_once(
        task, model="test-model", backend="test", provider="test",
        profile_name="test", engine_factory=_fake_engine_factory(),
        max_tokens=128, run_once=_RUN_ONCE_OK, clock=lambda: 0.0)
    assert captured.get("exists") is True, captured
    facts = captured["facts"]
    assert facts["steer_delivered"][0]["text"].startswith("PERSON")
    assert facts["actions"][0]["tool_calls"] == ["write_file"]
    assert "turns" not in facts  # operator finding 3: minimal content
    assert result.verify_ok is None  # the spy returned (None, ...) —
    # proof the spy, not the original, ran at the real call site.
