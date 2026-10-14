"""A TUI session ends when it must, and a cancel ends only what it names.

The second review of the worker reproduced these; each test is one of its
findings, made deterministic where a race could be staged directly:

* a shutdown or a cancel landing as a turn starts must not be undone by
  the REPL lifting the abort at turn start (the question would wait with
  nobody to answer it);
* the window going away ends a session that waits on a long tool call or
  a ``!command``, under the REPL's leave deadline;
* a cancel refuses the questions of what it cancelled, not those of the
  next command; it skips the line it named, not one sent after it; and the
  engine's stop request does not outlive it;
* prompts keep the order they were sent in; a steer that comes too late is
  sent as the next message; a queued line is stripped like a typed one, and
  drops an automatic follow-up that was about to start;
* nothing a question or an answer carries leaks a credential.
"""

from __future__ import annotations

import io
import json
import time

import pytest

from delfin.agent.service import protocol
from tui_worker_client import WorkerClient, WorkerError


def _in_process(tmp_path, monkeypatch):
    from delfin.agent.service.fake_engine import FakeEngine
    from delfin.agent.service.worker import Worker

    monkeypatch.setenv("DELFIN_AGENT_SERVICE_DIR", str(tmp_path / "service"))
    out = io.BytesIO()
    worker = Worker(FakeEngine(str(tmp_path)), writer=protocol.LineWriter(out),
                    reader=io.BytesIO(), workspace=tmp_path, banner=[],
                    notices=[], fake=True, save=False)
    return worker, out


def _wait_idle(w: WorkerClient, after: int, timeout: float = 20) -> int:
    _, at = w.wait_for(lambda e: e["type"] == "session.state"
                       and e["state"] == "idle", after=after, timeout=timeout)
    return at


def _notice(w: WorkerClient, after: int, needle: str) -> dict:
    found, _ = w.wait_for(lambda e: e["type"] == "notice" and needle in e["text"],
                          after=after)
    return found


# -- a shutdown or a cancel as a turn starts ------------------------------------

def test_a_turn_starting_during_shutdown_cannot_reopen_the_gate(tmp_path, monkeypatch):
    worker, _ = _in_process(tmp_path, monkeypatch)
    worker.shutdown("shutdown")
    # What the REPL's turn() does first, had it got past the closing check:
    worker.agent._begin_turn_resets_the_abort()
    worker.broker.reset_abort()
    answer = worker.broker.callback("bash", {"command": "make deploy"}, "$ make deploy")
    assert answer is False                     # refused at once, never waits
    assert worker.broker.pending() == []


def test_a_turn_starting_during_a_cancel_keeps_its_questions_refused(tmp_path, monkeypatch):
    worker, _ = _in_process(tmp_path, monkeypatch)
    agent = worker.agent
    agent.busy = True                          # the turn has begun ...
    worker._m_session_cancel({})               # ... and Esc lands now,
    agent._begin_turn_resets_the_abort()       # just before the REPL's reset
    answer = worker.broker.callback("bash", {"command": "rm -r build"}, "$ rm -r build")
    assert answer is False
    agent.busy = False
    agent.idle = False
    agent._line_ready()                        # back at the prompt: lifted
    assert worker.broker.aborted is False


def test_a_cancel_racing_the_reset_at_turn_start_still_holds(tmp_path, monkeypatch):
    # The cancel arrives while the REPL is lifting the abort for the new
    # turn: whichever comes first, the cancel's refusal must be in force
    # afterwards.
    import threading

    from delfin.agent import repl

    worker, _ = _in_process(tmp_path, monkeypatch)
    agent = worker.agent
    original = repl.TerminalAgent._begin_turn_resets_the_abort
    racer: list[threading.Thread] = []

    def lifting_while_a_cancel_arrives(self):
        racer.append(threading.Thread(target=worker._m_session_cancel, args=({},)))
        racer[0].start()
        racer[0].join(0.3)
        original(self)

    monkeypatch.setattr(repl.TerminalAgent, "_begin_turn_resets_the_abort",
                        lifting_while_a_cancel_arrives)
    agent.busy = True
    agent._begin_turn_resets_the_abort()
    racer[0].join(5)
    assert worker.broker.aborted is True


# -- the window going away ---------------------------------------------------------

def test_the_window_going_away_during_a_long_tool_call_ends_promptly(tmp_path):
    w = WorkerClient(tmp_path)
    start = len(w.events)
    w.request("session.prompt", {"text": "longtool please"})
    w.wait_for("tool.started", after=start)
    began = time.monotonic()
    w.proc.stdin.close()
    assert w.ended(timeout=30) == 129
    assert time.monotonic() - began < 12        # not the tool's 20 s
    assert w.events[-1]["type"] == "session.ended"
    assert w.events[-1]["reason"] == "eof"


def test_the_window_going_away_during_a_long_command_ends_promptly(tmp_path):
    w = WorkerClient(tmp_path)
    start = len(w.events)
    w.request("session.prompt", {"text": "!sleep 20"})
    w.wait_for(lambda e: e["type"] == "session.state" and e["state"] == "running",
               after=start)
    time.sleep(0.5)
    began = time.monotonic()
    w.proc.stdin.close()
    assert w.ended(timeout=30) == 129
    assert time.monotonic() - began < 12


# -- a cancel ends what it names, and nothing after it ---------------------------

def test_a_cancel_does_not_refuse_the_next_commands_question(tmp_path):
    with WorkerClient(tmp_path, "--permission-mode", "default", delay=0.05) as w:
        start = len(w.events)
        w.request("session.prompt", {"text": "slow please"})
        w.wait_for("text.started", after=start)
        w.request("session.cancel")
        _, end = w.wait_for("turn.ended", after=start)
        _wait_idle(w, end)
        at = len(w.events)
        w.request("session.prompt", {"text": "!make deploy"})
        asked, _ = w.wait_for("permission.asked", after=at, timeout=10)
        w.request("permission.reply", {"requestId": asked["requestId"], "choice": "y"})
        _notice(w, at, "(fake) make deploy")


def test_a_cancel_during_a_command_does_not_skip_a_message_sent_after_it(tmp_path):
    with WorkerClient(tmp_path) as w:
        start = len(w.events)
        w.request("session.prompt", {"text": "!sleep 2"})
        w.wait_for(lambda e: e["type"] == "session.state" and e["state"] == "running",
                   after=start)
        w.request("session.cancel")
        w.request("session.prompt", {"text": "hello after the cancel"})
        _, end = w.wait_for("turn.ended", after=start, timeout=20)
        assert w.events[end].get("error") is None


def test_a_cancel_leaves_no_stop_behind_for_the_next_command(tmp_path):
    # A stop request left armed ends the next !command with 130.
    with WorkerClient(tmp_path, delay=0.05) as w:
        for when in ("during the turn", "before it started"):
            start = len(w.events)
            if when == "during the turn":
                w.request("session.prompt", {"text": "slow please"})
                w.wait_for("text.started", after=start)
            else:
                w.send("session.prompt", {"text": "slow please"})
            w.request("session.cancel")
            _, end = w.wait_for("turn.ended", after=start)
            _wait_idle(w, end)
            at = len(w.events)
            w.request("session.prompt", {"text": "!sleep 1"})
            done, _ = w.wait_for("notice", after=at)
            assert "(fake) sleep 1" in done["text"], (when, done["text"])


def test_a_message_queued_after_a_cancel_runs_its_own_turn(tmp_path):
    with WorkerClient(tmp_path, delay=0.05) as w:
        start = len(w.events)
        w.request("session.prompt", {"text": "slow please"})
        w.wait_for("text.started", after=start)
        w.request("session.cancel")
        w.request("session.prompt", {"text": "hello again"})
        _, first = w.wait_for("turn.ended", after=start)
        _, second = w.wait_for("turn.ended", after=first + 1)
        assert w.events[second].get("error") is None


# -- prompts in order; late steers; stripping; follow-ups -------------------------

def _send_together(w: WorkerClient, *texts: str) -> None:
    lines = []
    for text in texts:
        w._next_id += 1
        lines.append(json.dumps({"jsonrpc": "2.0", "id": w._next_id,
                                 "method": "session.prompt",
                                 "params": {"text": text}}) + "\n")
    w.proc.stdin.write("".join(lines).encode())
    w.proc.stdin.flush()


def test_prompts_sent_together_are_shown_and_run_in_the_order_sent(tmp_path):
    with WorkerClient(tmp_path) as w:
        for n in range(10):
            start = len(w.events)
            sent = [f"first {n}", f"second {n}", f"third {n}"]
            _send_together(w, *sent)
            at = start
            for _ in sent:
                _, at = w.wait_for("turn.ended", after=at)
                at += 1
            shown = [e["text"] for e in w.events[start:at]
                     if e["type"] == "user.message"]
            assert shown == sent
            # Each turn answers its own prompt (the fake says which).
            ran = [e["text"].rsplit("(re: ", 1)[-1].rstrip(")\n")
                   for e in w.events[start:at] if e["type"] == "text.ended"]
            assert ran == sent
            _wait_idle(w, at - 1)


def test_a_steer_that_comes_after_the_answer_goes_out_next(tmp_path):
    with WorkerClient(tmp_path) as w:
        for n in range(5):
            start = len(w.events)
            w.request("session.prompt", {"text": f"hello {n}"})
            w.wait_for("text.ended", after=start)
            got = w.request("session.prompt", {"text": f"late steer {n}",
                                               "how": "steer"})
            if got.get("steered"):
                # Accepted while the engine's turn was still live: it must
                # then reach that turn or go out next, never a later one.
                _, end = w.wait_for("turn.ended", after=start)
                _, nxt = w.wait_for("turn.ended", after=end + 1, timeout=20)
            else:
                # Queued behind the ending turn, or -- when the turn had
                # already ended -- an ordinary message, handled next. Either
                # way it is the next turn's.
                assert got.get("queued") or got == {"accepted": True}
                _, end = w.wait_for("turn.ended", after=start)
                _, nxt = w.wait_for("turn.ended", after=end + 1, timeout=20)
                said = [e["text"] for e in w.events[end:nxt] if e["type"] == "user.message"]
                assert f"late steer {n}" in said or got.get("queued")
            _wait_idle(w, nxt)
            # Nothing of it is left to surface in an unrelated turn.
            s, e = w.turn("steer test")
            assert "No message arrived." in w.text_of(s, e)
            _wait_idle(w, e)


def test_a_queued_line_is_stripped_like_a_typed_one(tmp_path):
    with WorkerClient(tmp_path, delay=0.05) as w:
        start = len(w.events)
        w.request("session.prompt", {"text": "slow please"})
        w.wait_for("turn.started", after=start)
        w.request("session.prompt", {"text": "   !echo hi   "})
        _notice(w, start, "(fake) echo hi")             # a shell escape


def test_a_line_sent_as_a_follow_up_starts_drops_the_follow_up(tmp_path, monkeypatch):
    # The REPL looks for queued lines before it shows and starts a
    # follow-up; a line queued in between must still go first.
    from delfin.agent.repl import HarnessText

    worker, out = _in_process(tmp_path, monkeypatch)
    agent = worker.agent
    agent.idle = False
    agent.queued.append("hello instead")
    agent._show_user_input("continue", queued=True)
    result = agent.turn(HarnessText("continue"))
    assert not result.error and not result.text
    sent = [json.loads(line)["params"] for line in out.getvalue().splitlines()]
    assert not any(e["type"] in ("turn.started", "user.message") for e in sent)
    assert agent.queued == ["hello instead"]


def test_a_wake_up_is_not_dropped_for_a_queued_line(tmp_path, monkeypatch):
    # A finished job or a message from another session is not a follow-up:
    # its turn runs, and the queued line after it.
    from delfin.agent.repl import HarnessText

    worker, out = _in_process(tmp_path, monkeypatch)
    agent = worker.agent
    agent.idle = False
    agent._wake_turn = True
    agent.queued.append("hello after")
    agent.turn(HarnessText("[job finished] the ORCA run ended"))
    sent = [json.loads(line)["params"] for line in out.getvalue().splitlines()]
    assert [e["source"] for e in sent if e["type"] == "turn.started"] == ["harness"]


# -- credentials -----------------------------------------------------------------

def test_nothing_a_question_or_an_answer_carries_leaks_a_held_credential(tmp_path):
    fake = "fakeKitKeyForTheTuiTests0123456789"           # not a real key
    w = WorkerClient(tmp_path, "--permission-mode", "default",
                     env={"KIT_TOOLBOX_API_KEY": fake})
    try:
        start = len(w.events)
        w.request("session.prompt", {"text": f"!DEPLOY_KEY={fake} make deploy"})
        asked, _ = w.wait_for("permission.asked", after=start)
        again = w.request("initialize", {"protocolVersion": 1})
        with pytest.raises(WorkerError) as helped:
            w.request("permission.reply", {"requestId": asked["requestId"],
                                           "choice": "?"})
        kept = w.request("permission.reply", {"requestId": asked["requestId"],
                                              "choice": "k", "persist": True})
    finally:
        w.close()
    assert fake not in json.dumps(asked)
    assert fake not in json.dumps(again)
    assert fake not in helped.value.message
    assert fake not in json.dumps(kept)
    assert not any(fake.encode() in line for line in w.raw_lines)
    assert fake not in w.log_path().read_text()


# -- fake sessions; the schema ----------------------------------------------------

def test_a_fake_session_saves_nothing_unless_told(tmp_path):
    with WorkerClient(tmp_path) as w:
        w.turn("hello")
        session = w.init_result["session"]["sessionId"]
    saved = [p for p in (tmp_path / "child_scratch").rglob("*.json")
             if session in p.read_text(errors="replace")]
    assert saved == []


def test_every_event_and_result_fits_the_schema(tmp_path):
    jsonschema = pytest.importorskip("jsonschema")
    schema = protocol.json_schema()
    bad: list[str] = []

    def check(kind: str, name: str, value) -> None:
        try:
            jsonschema.validate(value, schema[kind][name])
        except (KeyError, jsonschema.ValidationError) as exc:
            bad.append(f"{kind} {name}: {str(exc)[:120]}")

    with WorkerClient(tmp_path, "--permission-mode", "default") as w:
        check("results", "initialize", w.init_result)
        for text in ("hello", "secret please", "latin please", "fail please"):
            w.turn(text)
        start = len(w.events)
        check("results", "session.prompt", w.request("session.prompt", {"text": "edit it"}))
        asked, _ = w.wait_for("permission.asked", after=start)
        check("results", "permission.reply",
              w.request("permission.reply", {"requestId": asked["requestId"],
                                             "choice": "y"}))
        w.wait_for("turn.ended", after=start)
        for method, params in (("session.status", {}), ("session.cancel", {}),
                               ("session.setMode", {"mode": "plan"})):
            check("results", method, w.request(method, params))
        replay = w.request("session.replay", {"turns": 2})
        check("results", "session.replay", replay)
        for e in replay["events"]:
            check("events", e["type"], e)
    for e in w.events:
        check("events", e["type"], e)
    assert bad == []
