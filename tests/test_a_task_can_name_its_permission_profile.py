"""Plan is a permission profile, and the suite could not express one.

``Task.mode`` is the agent mode — solo, dashboard, office, research — and
``load_mode`` raises on anything else. Its comment listed "plan" among the
values, so a task written against it failed at construction rather than
testing anything. Plan is a *permission* profile: read freely, write
nothing until the user approves. It is a safety boundary, and until Task
carried a profile nothing in the suite measured whether a model respects
it.

Measured once it could be: DeepSeek 3/3 at average quality 88.3, GLM 3/3
at 91.3. The behaviour was already right; what was missing was the ability
to see it.
"""

import pytest

from delfin.agent.benchmark import Task, _coerce_task


def test_a_task_carries_a_permission_profile():
    t = _coerce_task({"id": "x", "mode": "solo", "prompt": "p",
                      "permission_mode": "plan"})
    assert t.permission_mode == "plan"


def test_an_older_task_names_none_and_keeps_the_engine_default():
    """Every task written before this field existed must run as it did."""
    t = _coerce_task({"id": "y", "mode": "solo", "prompt": "p"})
    assert t.permission_mode == ""


def test_plan_is_not_an_agent_mode():
    """The comment that said it was sent tasks down a path that raises."""
    from delfin.agent.prompt_loader import PromptLoader

    with pytest.raises(ValueError):
        PromptLoader().load_mode("plan")


def test_the_factory_takes_the_profile_and_only_when_given():
    import inspect

    from delfin.agent import benchmark_runner as br

    sig = inspect.signature(br._default_engine_factory)
    assert "permission_mode" in sig.parameters
    src = inspect.getsource(br._default_engine_factory)
    # Empty means "leave the engine's own default alone".
    assert "if permission_mode:" in src


def test_the_runner_passes_it_and_tolerates_an_older_factory():
    """A test double written before this field must keep working."""
    import inspect

    from delfin.agent import benchmark_runner as br

    src = inspect.getsource(br)
    assert "permission_mode=task.permission_mode" in src
    assert "except TypeError:" in src


def test_the_shipped_plan_tasks_ask_for_the_profile():
    import glob
    from pathlib import Path

    import yaml

    found = []
    for f in glob.glob("delfin/agent/pack/benchmark/*.yaml"):
        rows = yaml.safe_load(Path(f).read_text())
        rows = rows if isinstance(rows, list) else rows.get("tasks", [])
        found += [t for t in rows if t.get("permission_mode") == "plan"]
    assert len(found) >= 3, [t["id"] for t in found]
    for t in found:
        assert t.get("mode") in ("solo", "dashboard", "office", "research")


def test_the_plan_tasks_judge_the_claim_and_not_the_attempt():
    """The permission layer exists so an attempt is harmless and records a
    denial. Two earlier versions of these tasks flagged denied attempts and
    failed an agent whose own words were that it was locked in plan mode
    and waiting for approval."""
    import glob
    from pathlib import Path

    import yaml

    for f in glob.glob("delfin/agent/pack/benchmark/*.yaml"):
        rows = yaml.safe_load(Path(f).read_text())
        rows = rows if isinstance(rows, list) else rows.get("tasks", [])
        for t in rows:
            if t.get("permission_mode") != "plan":
                continue
            for sig in (t.get("forbidden_signals") or []):
                pat = sig["pattern"]
                assert "write_file" not in pat, (t["id"], pat)
                assert "\\brm\\b" not in pat, (t["id"], pat)
