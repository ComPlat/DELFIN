"""One activity state, written by the engine, read by both surfaces.

Input: ``delfin.agent.activity.Activity`` driven by a scripted clock, by
the real ``AgentEngine.stream_response`` over a stub client with scripted
delays, and by the terminal's confirm broker. Output: the snapshot and
the two renderers (``line``, ``slow_hint``).

Asserted: the phase follows what the turn is doing (request out, first
token, tool dispatched, tool returned, person asked, turn ended); each
phase has its own clock; the clock advances with no model output at all;
a stop and an error end the phase as stopped/failed so no surface keeps
an indicator up; the slow hint names the stop key of the surface asking.
"""

from __future__ import annotations

import threading
import time
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from delfin.agent import activity as act


class _Clock:
    def __init__(self):
        self.t = 100.0

    def __call__(self):
        return self.t


# ---------------------------------------------------------------------------
# The state machine, on a scripted clock
# ---------------------------------------------------------------------------

def test_the_phases_and_their_own_clocks():
    clock = _Clock()
    a = act.Activity(now=clock)
    assert act.line(a.snapshot()) == "" and not a.snapshot().active
    a.begin_turn("kit.glm-5.3")
    clock.t += 45
    s = a.snapshot()
    assert s.phase == act.MODEL_WAIT and s.first
    assert act.line(s) == "waiting for kit.glm-5.3 … 45 s"
    a.output()
    clock.t += 12
    assert act.line(a.snapshot()) == "kit.glm-5.3 is answering … 12 s"
    a.tool_start("mcp__kit-coding__bash", '{"command": "pytest -q tests"}')
    clock.t += 75
    assert act.line(a.snapshot()) == \
        "running bash: pytest -q tests … 1 min 15 s"
    a.tool_end("bash")
    clock.t += 3
    s = a.snapshot()
    assert s.phase == act.MODEL_WAIT and s.first, (
        "the request carrying the result waits for its own first token")
    assert act.line(s) == "waiting for kit.glm-5.3 … 3 s"
    a.output()
    a.end_turn()
    s = a.snapshot()
    assert not s.active and act.line(s) == "done in 2 min 15 s"


def test_an_approval_inside_a_tool_gives_the_tool_its_clock_back():
    clock = _Clock()
    a = act.Activity(now=clock)
    a.begin_turn("m")
    a.tool_start("write_file", {"path": "/x/y/z/report.md"})
    clock.t += 2
    with a.awaiting(act.APPROVAL):
        clock.t += 30
        assert act.line(a.snapshot()) == "waiting for your approval … 30 s"
        assert act.slow_hint(a.snapshot()) == "", (
            "waiting on the person is not slowness to report")
    clock.t += 1
    assert act.line(a.snapshot()) == "running write_file: z/report.md … 33 s"
    with a.awaiting(act.ANSWER):
        assert act.line(a.snapshot()).startswith("waiting for your answer")


@pytest.mark.parametrize("how,word", [("error", "failed after"),
                                      ("stop", "stopped after")])
def test_a_turn_that_fails_or_is_stopped_ends_the_indicator(how, word):
    clock = _Clock()
    a = act.Activity(now=clock)
    a.begin_turn("m")
    a.tool_start("bash", {"command": "sleep 100"})
    clock.t += 4
    a.end_turn(error="boom" if how == "error" else "",
               stopped=how == "stop")
    s = a.snapshot()
    assert not s.active and act.line(s) == f"{word} 4 s"


def test_the_slow_hint_names_what_is_slow_and_the_surfaces_stop():
    clock = _Clock()
    a = act.Activity(now=clock)
    a.begin_turn("kit.glm-5.3")
    clock.t += 59
    assert act.slow_hint(a.snapshot()) == ""
    clock.t += 2
    term = act.slow_hint(a.snapshot(), "terminal")
    dash = act.slow_hint(a.snapshot(), "dashboard")
    assert "kit.glm-5.3" in term and "no token" in term and "Esc" in term
    assert "Stop" in dash and "Esc" not in dash
    a.output()
    a.tool_start("bash", {"command": "pytest"})
    clock.t += 119
    assert act.slow_hint(a.snapshot()) == "", "a shell gets two minutes"
    clock.t += 2
    hint = act.slow_hint(a.snapshot(), "terminal")
    assert "bash" in hint and "once it returns" in hint, (
        "a stop cannot cut a running call, and the hint must not claim it")
    a.tool_end()
    a.tool_start("read_file", {"path": "a.py"})
    clock.t += 31
    assert "read_file" in act.slow_hint(a.snapshot(), "dashboard")


# ---------------------------------------------------------------------------
# The engine writes it, with no output needed for the clock to move
# ---------------------------------------------------------------------------

class _Client:
    model = "stub-model"

    def __init__(self, script):
        self.script = script
        self.engine = None

    def stream_message(self, *, messages, system="", max_tokens=0, **kw):
        for step in self.script:
            if step[0] == "sleep":
                time.sleep(step[1])
            elif step[0] == "raise":
                raise RuntimeError(step[1])
            elif step[0] == "stop":
                self.engine.request_stop()
            elif step[0] == "text":
                yield SimpleNamespace(type="text_delta", text=step[1])
            elif step[0] == "tool_use":
                yield SimpleNamespace(type="tool_use", tool_name=step[1],
                                      tool_input=step[2], text="")
            elif step[0] == "tool_result":
                yield SimpleNamespace(type="tool_result", tool_name=step[1],
                                      tool_output="ok", text="")
        yield SimpleNamespace(type="message_delta", input_tokens=10,
                              output_tokens=1, cost_usd=0.0,
                              stop_reason="end_turn", text="")


def _engine(tmp_path, script):
    from delfin.agent.engine import AgentEngine
    from tests.test_t2_engine_fresh_start import _pack_tree
    root = _pack_tree(tmp_path)
    client = _Client(script)
    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=root, backend="cli", mode="quick",
                             pack_dir=root)
    engine.client = client
    client.engine = engine
    return engine


@pytest.fixture(autouse=True)
def _no_archive(monkeypatch):
    import delfin.agent.session_store as store
    monkeypatch.setattr(store, "archive_pre_compaction_transcript",
                        lambda *a, **k: None)


def _sampled(engine, run, period=0.25):
    seen: list[tuple[float, str, str]] = []
    done = threading.Event()

    def _sample():
        t0 = time.monotonic()
        while not done.is_set():
            s = engine.activity.snapshot()
            seen.append((time.monotonic() - t0, s.phase, act.line(s)))
            done.wait(period)

    sampler = threading.Thread(target=_sample, daemon=True)
    sampler.start()
    try:
        run()
    except Exception:
        pass
    finally:
        done.set()
        sampler.join(2)
    return seen


def test_a_five_second_first_token_is_counted_while_nothing_arrives(tmp_path):
    engine = _engine(tmp_path, [("sleep", 5.0), ("text", "hello")])
    seen = _sampled(engine, lambda: engine.stream_response(
        user_message="hi"))
    waiting = [ln for _t, ph, ln in seen if ph == act.MODEL_WAIT]
    secs = sorted({int(ln.rsplit("… ", 1)[1].split()[0]) for ln in waiting})
    assert waiting and waiting[0].startswith("waiting for stub-model …")
    assert secs[0] <= 1 and secs[-1] >= 4, (
        f"the clock did not run through the silence: {secs}")
    s = engine.activity.snapshot()
    assert not s.active and act.line(s).startswith("done in")


def test_a_tool_run_reads_as_running_then_back_to_waiting(tmp_path):
    engine = _engine(tmp_path, [
        ("tool_use", "mcp__kit-coding__bash", '{"command": "sleep 3"}'),
        ("sleep", 1.5),
        ("tool_result", "mcp__kit-coding__bash"),
        ("sleep", 1.0),
        ("text", "done")])
    seen = _sampled(engine, lambda: engine.stream_response(
        user_message="run it"))
    lines = [ln for _t, _ph, ln in seen]
    assert any(ln.startswith("running bash: sleep 3 …") for ln in lines), lines
    phases = [ph for _t, ph, _ln in seen]
    order = [ph for i, ph in enumerate(phases) if i == 0 or ph != phases[i - 1]]
    assert act.TOOL in order and order[order.index(act.TOOL) + 1:][:1] == [
        act.MODEL_WAIT], (
        "after the result the next request waits for its first token: "
        f"{order}")


@pytest.mark.parametrize("script,word", [
    ([("sleep", 0.3), ("raise", "connection reset")], "failed after"),
    ([("sleep", 0.3), ("stop", None)], "stopped after"),
])
def test_an_error_or_stop_clears_the_indicator(tmp_path, script, word):
    engine = _engine(tmp_path, script)
    try:
        engine.stream_response(user_message="hi")
    except Exception:
        pass
    s = engine.activity.snapshot()
    assert not s.active, s
    assert act.line(s).startswith(word), act.line(s)


# ---------------------------------------------------------------------------
# The surfaces that ask the person mark the wait
# ---------------------------------------------------------------------------

def test_the_terminal_broker_marks_an_approval_and_a_question(tmp_path,
                                                             monkeypatch):
    from delfin.agent import terminal_confirm as tc
    monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "pending")
    a = act.Activity()
    a.begin_turn("m")
    a.tool_start("bash", {"command": "ls"})
    broker = tc.TerminalConfirmBroker(session_key="t")
    broker.activity = a
    for call, want in ((lambda: broker.callback("bash", {"command": "ls"},
                                                "ls"), "approval"),
                       (lambda: broker.ask_user({"question": "which?"}),
                        "answer")):
        worker = threading.Thread(target=call, daemon=True)
        worker.start()
        deadline = time.monotonic() + 5
        req = None
        while req is None and time.monotonic() < deadline:
            req = broker.take()
            time.sleep(0.02)
        assert req is not None
        assert act.line(a.snapshot()).startswith(f"waiting for your {want}")
        broker.resolve(req, False if want == "approval" else {"answers": []})
        worker.join(5)
        assert a.snapshot().phase == act.TOOL


def test_the_dashboard_broker_marks_an_approval():
    from delfin.agent.kit_confirm import KitConfirmBroker
    a = act.Activity()
    a.begin_turn("m")
    broker = KitConfirmBroker(default_timeout_s=10)
    broker.follow_activity(lambda: a)
    worker = threading.Thread(
        target=lambda: broker.callback("bash", {"command": "ls"}, "ls"),
        daemon=True)
    worker.start()
    deadline = time.monotonic() + 5
    while not broker.pending_requests() and time.monotonic() < deadline:
        time.sleep(0.02)
    assert act.line(a.snapshot()).startswith("waiting for your approval")
    broker.resolve(broker.pending_requests()[0], False)
    worker.join(5)
    assert a.snapshot().phase == act.MODEL_WAIT
