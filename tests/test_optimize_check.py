"""Tests for the safe-optimisation pre-flight check."""

from __future__ import annotations

from types import SimpleNamespace

from delfin.agent import optimize_check as oc


def _sig(pattern, against="text"):
    return SimpleNamespace(pattern=pattern, against=against, optional=False)


def _task(tid, prompt="do x", mode="dashboard", expected=None, forbidden=None):
    return SimpleNamespace(
        id=tid, prompt=prompt, mode=mode, task_class="x",
        expected_signals=expected if expected is not None else [_sig("ok")],
        forbidden_signals=forbidden or [],
    )


def test_real_repo_is_safe():
    # The committed suite + ground-truth + prompts must pass with no errors.
    issues = oc.run_checks()
    errors = [i for i in issues if i.severity == "error"]
    assert errors == [], f"unexpected errors: {[str(e) for e in errors]}"


def test_detects_bad_regex(monkeypatch):
    monkeypatch.setattr(
        "delfin.agent.benchmark.load_tasks",
        lambda: [_task("t1", expected=[_sig("(unclosed")])],
    )
    issues = oc.check_benchmark_tasks()
    assert any("bad regex" in i.message for i in issues if i.severity == "error")


def test_detects_duplicate_ids(monkeypatch):
    monkeypatch.setattr(
        "delfin.agent.benchmark.load_tasks",
        lambda: [_task("dup"), _task("dup")],
    )
    issues = oc.check_benchmark_tasks()
    assert any("duplicate task id" in i.message for i in issues)


def test_detects_empty_prompt_and_bad_against(monkeypatch):
    monkeypatch.setattr(
        "delfin.agent.benchmark.load_tasks",
        lambda: [_task("t", prompt="  ", expected=[_sig("x", against="bogus")])],
    )
    issues = oc.check_benchmark_tasks()
    assert any("empty prompt" in i.message for i in issues)
    assert any("invalid against" in i.message for i in issues)


def test_todo_placeholder_is_warned(monkeypatch):
    monkeypatch.setattr(
        "delfin.agent.benchmark.load_tasks",
        lambda: [_task("t", expected=[_sig("TODO-fill-me")])],
    )
    issues = oc.check_benchmark_tasks()
    assert any("TODO" in i.message and i.severity == "warn" for i in issues)


def test_prompts_present():
    issues = oc.check_agent_prompts()
    assert [i for i in issues if i.severity == "error"] == []


def test_ground_truth_intact():
    issues = oc.check_ground_truth()
    assert [i for i in issues if i.severity == "error"] == []
