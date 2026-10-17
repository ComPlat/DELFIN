"""Green proof for the KG4 steer-injection + transport patch, run against
the .gate/px/ patched copies. The patched module source is exec'd with the
real package context (delfin.agent), so its relative imports resolve, and
the patched classes/functions take effect. The operator applies the patch
to the real tree; this proves the patched CODE is green, not the tree.
"""
from __future__ import annotations

import json
import sys
import types
from pathlib import Path

PX = Path(__file__).resolve().parents[1] / ".gate" / "px"


def _load_patched(pkg_name: str, file_name: str, real_name: str):
    """Exec the patched file's source inside the real package context.

    ``sys.modules[real_name]`` is swapped BEFORE the source runs, so the
    patched file's own ``from .benchmark import …`` resolves to the
    patched sibling, not the real one.
    """
    import delfin.agent as agent_pkg
    real = sys.modules[real_name]
    patched = types.ModuleType(real_name)
    patched.__package__ = real.__package__
    patched.__file__ = str(PX / file_name)
    sys.modules[real_name] = patched
    setattr(agent_pkg, pkg_name, patched)
    src = (PX / file_name).read_text(encoding="utf-8")
    exec(compile(src, patched.__file__, "exec"), patched.__dict__)
    return patched


def test_patched_task_and_runner_deliver_the_steer(tmp_path):
    _load_patched("benchmark", "new_benchmark.py",
                  "delfin.agent.benchmark")
    _load_patched("benchmark_runner", "new_benchmark_runner.py",
                  "delfin.agent.benchmark_runner")
    from delfin.agent.benchmark import Task
    from delfin.agent import benchmark_runner as br

    task = Task(
        id="lr_steered_goal", task_class="office_ops", mode="solo",
        prompt="old goal",
        mid_run_messages=(("after", "PERSON (mid-task steering, new "
                                   "goal): follow the new goal"),),
    )
    calls = []

    def fake_run_once(engine, prompt, **kwargs):
        calls.append(prompt)
        return {"text": "ok", "tool_calls": [], "input_tokens": 1,
                "output_tokens": 1, "error": ""}

    class _FakeEngine:
        token_usage = {"input": 1, "output": 1}
        cost_usd = 0.0
        messages = [{"role": "assistant", "content":
                     "wrote answer.txt with the new goal"}]

    result = br._run_task_once(
        task, model="test-model", backend="test", provider="test",
        profile_name="test",
        engine_factory=lambda *a, **k: _FakeEngine(),
        max_tokens=128, run_once=fake_run_once, clock=lambda: 0.0,
    )
    assert len(calls) == 2, calls
    assert "new goal" in calls[1]


def test_patched_loader_parses_mid_run_messages_from_yaml(tmp_path):
    _load_patched("benchmark", "new_benchmark.py",
                  "delfin.agent.benchmark")
    from delfin.agent import benchmark as bench

    yaml_text = """\
- id: lr_steered_goal
  task_class: office_ops
  mode: solo
  prompt: old goal
  mid_run_messages:
    - when: after
      text: >-
        PERSON (mid-task steering, new goal): follow the new goal
"""
    p = tmp_path / "tasks.yaml"
    p.write_text(yaml_text, encoding="utf-8")
    tasks = bench.load_tasks(p)
    assert len(tasks) == 1
    t = tasks[0]
    assert len(t.mid_run_messages) == 1
    when, text = t.mid_run_messages[0]
    assert when == "after"
    assert "new goal" in text
    # A task WITHOUT the declaration loads empty — no regression.
    p2 = tmp_path / "tasks2.yaml"
    p2.write_text("- id: plain\n  task_class: office_ops\n  mode: solo\n"
                  "  prompt: nothing\n", encoding="utf-8")
    plain = bench.load_tasks(p2)[0]
    assert plain.mid_run_messages == ()
