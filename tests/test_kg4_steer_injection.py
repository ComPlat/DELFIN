"""Red control for the steer-injection hook (KG4 phase 3 condition 1).

lr_steered_goal needs the harness to deliver a scripted person message
mid-run: after the first turn's run_once, a second turn whose user
message is the steering text, recorded as an intervention. Today
benchmark_runner._run_task_once runs exactly one turn (line 1049) and
there is no field to declare one — the red test proves both ends.
"""
from __future__ import annotations

from pathlib import Path

BENCH = Path(__file__).resolve().parents[1] / "delfin" / "agent"


def _fake_task():
    """A minimal stand-in with the steering declaration the yaml sets."""
    from delfin.agent.benchmark import Task

    return Task(
        id="lr_steered_goal",
        task_class="generalist",
        mode="solo",
        prompt="do the old goal",
        mid_run_messages=(("after", "PERSON (mid-task steering, new goal): "
                           "answer only the new goal"),),
    )


def test_task_accepts_mid_run_messages_field():
    from delfin.agent.benchmark import Task

    task = _fake_task()
    assert len(task.mid_run_messages) == 1
    when, text = task.mid_run_messages[0]
    assert when == "after"
    assert "new goal" in text


def test_run_task_once_delivers_the_steer_as_a_second_turn():
    """After turn 1, the declared steering message runs as its own turn
    and the fact lands on the raw run (user_interventions / steer text),
    so the accept script can require steer_followed from the records."""
    from delfin.agent import benchmark_runner as br

    task = _fake_task()
    calls = []

    def fake_run_once(engine, prompt, **kwargs):
        calls.append(prompt)
        return {"text": "ok", "tool_calls": [], "input_tokens": 1,
                "output_tokens": 1, "error": ""}

    class _FakeEngine:
        token_usage = {"input": 1, "output": 1}
        cost_usd = 0.0
        messages = []

    result = br._run_task_once(
        task,
        model="test-model", backend="test", provider="test",
        profile_name="test",
        engine_factory=lambda *a, **k: _FakeEngine(),
        max_tokens=128, run_once=fake_run_once, clock=lambda: 0.0,
    )
    # The steering text was delivered as a SECOND turn, after turn 1.
    assert len(calls) == 2, calls
    assert "new goal" in calls[1]
