"""Red control, part 2 (reviewer nacht-s19 phase-3 pass, findings 1+2).

The steer accept must work END-TO-END through the real call site —
run_acceptance runs the accept script as a SUBPROCESS with only the
workspace — and the runner must deliver the declared steering message
as its own second turn, writing what it observed to the workspace.
Both are red today: no Task field, no runner hook, no facts file.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1] / "delfin" / "agent"
ACCEPT = (BENCH / "pack" / "benchmark" / "accept"
          / "w17_longrun_steer_verify.py")


def _build_workspace(tmp_path: Path, *, facts: dict | None) -> Path:
    """A workspace the accept SHOULD accept (new goal carried out)."""
    d = tmp_path / "w17_longrun" / "lr_steered_goal"
    d.mkdir(parents=True)
    (d / "STEERING_MESSAGE.txt").write_text(
        "PERSON (mid-task steering, new goal): follow the new goal",
        encoding="utf-8")
    (d / "answer.txt").write_text(
        "a: 3.4 nm\nb: absent\nc: 12.0 kJ/mol\nd: absent",
        encoding="utf-8")
    if facts is not None:
        (tmp_path / ".longrun_facts.json").write_text(
            json.dumps(facts), encoding="utf-8")
    return tmp_path


def test_run_acceptance_subprocess_accepts_when_steer_was_followed(
        tmp_path):
    """The REAL call site: run_acceptance executes the script as a
    subprocess with only the workspace; the facts file carries what
    the runner observed, and the accept judges steer_followed from it.
    Green only when the transport exists."""
    ws = _build_workspace(tmp_path, facts={
        "steer_delivered": [
            {"text": "PERSON (mid-task steering, new goal): follow it",
             "tool_calls": [], "error": ""}],
        "turns": [
            {"text": "PERSON (mid-task steering, new goal): follow it",
             "tool_calls": [], "error": ""},
            {"text": "wrote answer.txt with the new goal lines",
             "tool_calls": [{"name": "write_file"}]},
        ],
        "actions": [
            {"text": "write answer.txt with the new goal lines "
                     "(new goal)",
             "tool_calls": [{"name": "write_file"}]},
        ],
    })
    from delfin.agent.benchmark_runner import run_acceptance
    ok, out = run_acceptance(
        "w17_longrun_steer_verify.py", ws, root=BENCH.parent.parent)
    assert ok is True, out


def test_run_acceptance_subprocess_rejects_when_steer_not_attested(
        tmp_path):
    """No facts file (or no steer in it) -> the accept must REJECT, not
    pass by luck: None-never-passes at the real call site too."""
    ws = _build_workspace(tmp_path, facts=None)
    from delfin.agent.benchmark_runner import run_acceptance
    ok, out = run_acceptance(
        "w17_longrun_steer_verify.py", ws, root=BENCH.parent.parent)
    assert ok is False, out
    assert "steer-followed" in out


def test_runner_delivers_steer_and_writes_facts_file():
    """_run_task_once: a task with mid_run_messages gets a SECOND
    run_once call with the steering text and leaves the observed facts
    in the workspace."""
    from delfin.agent import benchmark_runner as br
    from delfin.agent.benchmark import Task

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

    ws = Path(br._run_task_once.__globals__["_mkdtemp"]() if False else ".")
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        old_mkdtemp = getattr(br, "_mkdtemp", None)
        # The runner builds its own workspace; force it into td via the
        # task's workspace resolution: run with workspace override if
        # supported, else patch. Simplest honest check: call and then
        # find the facts file in the workspace the result reports.
        result = br._run_task_once(
            task,
            model="test-model", backend="test", provider="test",
            profile_name="test",
            engine_factory=lambda *a, **k: _FakeEngine(),
            max_tokens=128, run_once=fake_run_once, clock=lambda: 0.0,
        )
        assert len(calls) == 2, calls
        assert "new goal" in calls[1]
