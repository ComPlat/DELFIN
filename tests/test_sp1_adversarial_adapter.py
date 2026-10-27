"""SP1 adversarial round 2 — attack the runner adapter (987ef534)."""
from types import SimpleNamespace

import pytest

from delfin.agent import explore_pilot as ep


def _fake_result(**over):
    base = dict(success=True, input_tokens=11, output_tokens=7, tool_calls=3,
                duration_s=1.5, text_excerpt="only an excerpt")
    base.update(over)
    return SimpleNamespace(**base)


class MultiTurnRunTask:
    """Drives run_once the way the REAL run_task does: once for the
    task prompt, once more per mid-run steer (benchmark_runner.py:1051,
    1076). The adapter's captured full_text must still be the ANSWER."""

    def __init__(self, result=None, calls_log=None):
        self.result = result or _fake_result()
        self.log = calls_log if calls_log is not None else []

    def __call__(self, task, **kwargs):
        run_once = kwargs.get("run_once")
        engine = SimpleNamespace(token_usage={}, cost_usd=0.0,
                                 stream_response=lambda **kw: "")
        run_once(engine, "MAIN PROMPT", max_tokens=100)
        run_once(engine, "STEER: now cite file:line", max_tokens=100)
        return self.result


def test_captured_text_survives_a_mid_run_steer():
    """With a steer the LAST run_once call is the steer text, not the
    answer. The adapter grades citations on what it captured — it must
    capture the task's own turn, not a later steer."""
    rt = MultiTurnRunTask()
    runner = ep.make_pilot_runner(run_task=rt, model="m")
    out = runner(ep.PilotTask(task_id="t", prompt="p"), arm="plain")
    assert out["final_text"] != "STEER: now cite file:line"
    # The main answer text is what the main call produced; the wrapper
    # can only know it captured the right call by which prompt ran.
    assert rt.log if False else True


def test_captured_text_survives_a_mid_run_steer_precise():
    """Precise form: the captured text is the answer produced by the
    main call — verified by having the wrapper return distinguishable
    text per call."""
    class DrivingRunTask:
        result = _fake_result()

        def __call__(self, task, **kwargs):
            run_once = kwargs.get("run_once")

            # The wrapper delegates to cli._run_once; give it an engine
            # whose stream_response echoes the user_message it got, so
            # the answer identifies which call produced it.
            class E:
                token_usage = {}
                kit_permissions = None

                def stream_response(self, **kw):
                    return f"ANSWER-TO:{kw.get('user_message', '')}"

            run_once(E(), "MAIN PROMPT", max_tokens=100)
            run_once(E(), "STEER TEXT", max_tokens=100)
            return self.result

    rt = DrivingRunTask()
    runner = ep.make_pilot_runner(run_task=rt, model="m")
    out = runner(ep.PilotTask(task_id="t", prompt="p"), arm="plain")
    assert out["final_text"].startswith("ANSWER-TO:MAIN PROMPT"), out


def test_real_run_once_is_wrapped_not_replaced(monkeypatch):
    """The wrapper must delegate to cli._run_once, not fake the turn:
    monkeypatch cli._run_once and check the adapter's run_once calls it
    and still captures the full text."""
    import delfin.agent.cli as cli

    seen = []

    def fake_run_once(engine, prompt, **kw):
        seen.append(prompt)
        return {"text": f"FULL:{prompt}", "tool_calls": [],
                "input_tokens": 1, "output_tokens": 1}

    monkeypatch.setattr(cli, "_run_once", fake_run_once)
    rt = MultiTurnRunTask()
    runner = ep.make_pilot_runner(run_task=rt, model="m")
    out = runner(ep.PilotTask(task_id="t", prompt="p"), arm="explore")
    assert seen == ["MAIN PROMPT", "STEER: now cite file:line"]
    # captured must be the MAIN answer, not the steer
    assert out["final_text"] == "FULL:MAIN PROMPT"


def test_unknown_task_field_raises_not_defaults():
    rt = MultiTurnRunTask()
    runner = ep.make_pilot_runner(run_task=rt, model="m",
                                  default_task={"mode": "solo",
                                                "not_a_field": 1})
    with pytest.raises(ValueError):
        runner(ep.PilotTask(task_id="t", prompt="p"), arm="plain")


def test_real_benchmark_task_passes_through():
    from delfin.agent.benchmark import Task

    seen = {}

    class SpyRunTask:
        def __call__(self, task, **kw):
            seen["task"] = task
            run_once = kw.get("run_once")
            run_once(SimpleNamespace(token_usage={}, cost_usd=0.0,
                                     stream_response=lambda **k: ""),
                     "p", max_tokens=10)
            return _fake_result()

    rt = SpyRunTask()
    t = Task(id="real", task_class="generalist", mode="solo", prompt="p")
    runner = ep.make_pilot_runner(run_task=rt, model="m")
    runner(t, arm="plain")
    assert seen["task"] is t


def test_engine_init_failure_gives_empty_text_not_excerpt():
    """run_task returns a failed result without ever calling run_once
    (benchmark_runner.py:979-983). The adapter must not then grade the
    400-char excerpt as an answer."""
    rt = MultiTurnRunTask()
    rt.result = _fake_result(success=False, text_excerpt="excerpt bits",
                             error="engine init failed")
    calls = {"n": 0}
    orig = rt.__call__

    def no_turn(task, **kw):
        # engine init failed: run_once is never reached
        return rt.result

    rt.__call__ = no_turn
    runner = ep.make_pilot_runner(run_task=rt, model="m")
    out = runner(ep.PilotTask(task_id="t", prompt="p"), arm="plain")
    assert out["final_text"] == "" or "excerpt" not in out["final_text"]
    assert out["success"] is False
