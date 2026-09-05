"""The benchmark demanded a CONTROL key that DELFIN does not have.

``workflow_plan_before_act`` expected the agent to write

    /control key <something containing "dispersion"> D3BJ

There is no such key. CONTROL.txt calls it ``disp_corr``; ``dispersion``
is the name of a DASHBOARD dropdown field, which is a different surface
reached by a different command. So the pattern rewarded a command that
would not work and marked the working one wrong.

It stayed hidden because the model used to write the broken form. When
it started writing the real key the task went from 100% to 0%, and the
run looked like a regression caused by that day's changes -- five
failures out of five, low variance, exactly the shape of a real defect.
It took reading the actual CONTROL keys to see that the agent had got
BETTER and the benchmark was punishing it for it.

A benchmark that scores against names the product does not use will
quietly drive the agent away from correct behaviour, and every wrong step
it causes looks like evidence. So the invariant is asserted directly:
every control-key pattern has to be satisfiable by a key that exists.
"""

from __future__ import annotations

import pathlib
import re

import pytest
import yaml


_TASKS = (pathlib.Path(__file__).resolve().parents[1] / "delfin" / "agent"
          / "pack" / "benchmark" / "tasks.yaml")
_AGENT_DOC = (pathlib.Path(__file__).resolve().parents[1] / "delfin" / "agent"
              / "pack" / "agents" / "dashboard_agent.md")


def _real_control_keys() -> set[str]:
    """The CONTROL.txt keys, read from the reference the agent is given."""
    text = _AGENT_DOC.read_text(encoding="utf-8")
    i = text.index("CONTROL.txt — quick reference")
    block = text[i:i + 500]
    return {k.lower() for k in re.findall(r"`([A-Za-z_][A-Za-z0-9_]*)`", block)}


def _control_key_fragments() -> list[tuple[str, str]]:
    """(task_id, fragment) for every ``\\w*NAME\\w*`` in a control-key
    pattern -- the part that has to match a real key name."""
    data = yaml.safe_load(_TASKS.read_text(encoding="utf-8"))
    tasks = data if isinstance(data, list) else data.get("tasks", [])
    out: list[tuple[str, str]] = []
    for task in tasks:
        for sig in task.get("expected_signals") or []:
            pattern = str(sig.get("pattern", ""))
            if "control" not in pattern or "key" not in pattern:
                continue
            for frag in re.findall(r"\\w\*(\w+)\\w\*", pattern):
                out.append((task["id"], frag))
    return out


def test_the_reference_lists_the_keys():
    keys = _real_control_keys()
    assert "disp_corr" in keys
    assert "functional" in keys
    assert "dispersion" not in keys, (
        "if CONTROL really gained this key, the pattern was right after all")


def test_the_benchmark_has_control_key_patterns_at_all():
    """A guard on the guard: if the extraction stops finding anything,
    the test below would pass by doing nothing."""
    assert _control_key_fragments()


@pytest.mark.parametrize("task_id,fragment", _control_key_fragments())
def test_every_expected_control_key_exists(task_id, fragment):
    keys = _real_control_keys()
    assert any(fragment.lower() in key for key in keys), (
        f"{task_id} expects a control key containing '{fragment}', and no "
        f"real key does — the task scores correct behaviour as wrong")


def test_the_real_dispersion_key_satisfies_its_pattern():
    """The concrete command the agent should write must pass the check."""
    data = yaml.safe_load(_TASKS.read_text(encoding="utf-8"))
    tasks = data if isinstance(data, list) else data.get("tasks", [])
    task = next(t for t in tasks if t["id"] == "workflow_plan_before_act")
    patterns = [s["pattern"] for s in task["expected_signals"]
                if "disp" in s["pattern"]]
    assert patterns, "the dispersion expectation disappeared"
    assert re.search(patterns[0], "ACTION: /control key disp_corr D3BJ")


def test_the_broken_form_is_no_longer_required():
    """Accepting the real key must not mean ONLY the real key: a model
    that writes the older field name is still doing the right thing by the
    task's intent, and the point here is to stop scoring on a typo."""
    data = yaml.safe_load(_TASKS.read_text(encoding="utf-8"))
    tasks = data if isinstance(data, list) else data.get("tasks", [])
    task = next(t for t in tasks if t["id"] == "workflow_plan_before_act")
    patterns = [s["pattern"] for s in task["expected_signals"]
                if "disp" in s["pattern"]]
    assert re.search(patterns[0], "ACTION: /control key dispersion D3BJ")
