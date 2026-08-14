"""A Stop must survive the next message the user types.

WHAT THE USER SAW. A long tool call was running, they pressed Stop, the
spinner went out, and they typed the next instruction. The reply was "A
turn is already running on this session" — and then, half a minute later,
the agent finished the work they had just abandoned and billed for it.
Pressing Stop again did nothing at all.

WHAT ACTUALLY HAPPENED. The stop was one boolean on the engine, owned by
nobody. Stop set it and deliberately did not wait for the worker, so the
widget's ``streaming`` flag was already False when the next message
arrived; the send therefore started a NEW turn instead of steering the
running one. That turn cleared the flag on its way in — the reset was the
first statement of ``stream_response``, long before the one-turn gate
refused it. So a turn that never ran, and never owned anything, deleted
the only brake the running turn had; the running turn polls exactly that
flag between rounds, found it clear, and carried on. Every later Stop was
erased the same way, by the same route.

WHY THE OBVIOUS FIX IS NOT THE FIX. "Do not clear the flag in the
dashboard" is not enough: the engine cleared it too, and any other caller
would keep the hazard alive. Nor is "make Stop join the worker" — Stop
must stay instant, and joining a turn that is inside a ten-minute command
would freeze the UI. What was missing is an OWNER: the stop belongs to the
turn it interrupted, so only a turn that has actually acquired the gate
may clear it, and only once the turn that owns it has finished.

THE SECOND HALF, same region. The stop check sat at the top of the event
loop, before dispatch. What a client sends the instant it sees the flag is
its closing notice and then the message_delta with the turn's tokens and
cost — so the break discarded the accounting for every stopped turn. A
turn could burn several rounds, be stopped, and report $0.00; the
run-budget gate never saw a cent of it.
"""

from __future__ import annotations

import threading
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import engine as E
from delfin.agent.api_client import StreamEvent


@pytest.fixture
def eng(tmp_path):
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        return E.AgentEngine(
            repo_dir=tmp_path, backend="api", provider="kit",
            model="kit.qwen3.5-397b-A17b", mode="solo")


# ---------------------------------------------------------------------------
# The incident, end to end
# ---------------------------------------------------------------------------

def test_the_next_message_cannot_erase_a_running_turns_stop(eng):
    """Stop mid-tool, send the next message, and the stopped turn stays
    stopped instead of resuming the abandoned plan."""
    in_the_tool = threading.Event()
    release = threading.Event()
    resumed = threading.Event()

    def long_tool_call(**kw):
        yield StreamEvent(type="text_delta", text="running the long job")
        in_the_tool.set()
        release.wait(timeout=5)
        # Exactly what the real tool loop does between rounds.
        if eng.client.should_stop():
            yield StreamEvent(type="text_delta", text="\nstopped\n")
            yield StreamEvent(type="message_delta", input_tokens=10,
                              output_tokens=5, cost_usd=0.10,
                              stop_reason="stopped")
            return
        resumed.set()
        yield StreamEvent(type="text_delta", text=" ...and finished it")
        yield StreamEvent(type="message_delta", input_tokens=10,
                          output_tokens=900, cost_usd=2.00)

    eng.client.stream_message = MagicMock(side_effect=long_tool_call)
    first: list[str] = []
    worker = threading.Thread(
        target=lambda: first.append(eng.stream_response("do the big thing")))
    worker.start()
    assert in_the_tool.wait(timeout=5)

    eng.request_stop()                     # the user presses Stop
    eng.clear_stop()                       # the next send's worker
    second = eng.stream_response("actually, do this instead")

    assert eng._stop_requested is True, (
        "the refused turn erased the running turn's stop")
    assert "already running" in second.lower()

    release.set()
    worker.join(timeout=10)
    assert not resumed.is_set(), (
        "the stopped turn resumed and executed the abandoned plan")
    assert first and "finished it" not in first[0]


def test_a_refused_turn_changes_no_state_it_never_owned(eng):
    """The gate refuses before anything is touched — no reset, no message."""
    started = threading.Event()
    release = threading.Event()

    def slow(**kw):
        started.set()
        release.wait(timeout=5)
        yield StreamEvent(type="text_delta", text="done")

    eng.client.stream_message = MagicMock(side_effect=slow)
    t = threading.Thread(target=lambda: eng.stream_response("first"))
    t.start()
    assert started.wait(timeout=5)

    eng.request_stop()
    eng.stream_response("second")          # refused by the gate

    assert eng._stop_requested is True
    release.set()
    t.join(timeout=10)
    # The refused turn appended nothing. (The stopped turn's own message is
    # taken back out by the stop branch — it was never answered.)
    assert not any("second" in str(m.get("content")) for m in eng.messages)


def test_the_stop_is_released_once_its_own_turn_has_ended(eng):
    """The refusal is narrow: a finished turn's stop must not outlive it,
    or the first fix would simply re-create the original mute."""
    def fine(**kw):
        yield StreamEvent(type="text_delta", text="ok")

    eng.client.stream_message = MagicMock(side_effect=fine)
    eng.stream_response("hello")
    eng.request_stop()
    assert eng._stop_requested is True

    eng.clear_stop()                       # nothing is in flight any more
    assert eng._stop_requested is False
    assert eng.stream_response("again") == "ok"


def test_a_new_turn_clears_a_stop_left_over_from_an_older_one(eng):
    """Even without clear_stop: the turn that acquires the gate owns the
    flag from that moment on."""
    def fine(**kw):
        yield StreamEvent(type="text_delta", text="ok")

    eng.client.stream_message = MagicMock(side_effect=fine)
    eng.stream_response("first")
    eng.request_stop()
    assert eng.stream_response("second") == "ok"
    assert eng._stop_requested is False


def test_clear_stop_still_works_before_any_turn_has_run(eng):
    eng.request_stop()
    eng.clear_stop()
    assert eng._stop_requested is False


def test_the_reset_sits_behind_the_gate(eng):
    """Source-level: the reset above the gate is the defect itself."""
    import pathlib
    src = pathlib.Path(E.__file__).read_text(encoding="utf-8")
    body = src[src.index("def stream_response"):]
    body = body[:body.index("_turn_t0 = _time.monotonic()")]
    gate_at = body.index("with self._turn_gate:")
    reset_at = body.index("self._stop_requested = False")
    assert reset_at > gate_at, (
        "the stop reset runs before the turn gate again — a turn that is "
        "about to be refused can wipe a running turn's stop")


# ---------------------------------------------------------------------------
# A stopped turn still costs money
# ---------------------------------------------------------------------------

def test_a_stopped_turn_still_reports_what_it_spent(eng):
    """The accounting arrives after the flag is set, which is why the
    break at the top of the loop lost every stopped turn's spend."""
    notice: list[str] = []

    def stopped_mid_turn(**kw):
        yield StreamEvent(type="text_delta", text="round one")
        # The user presses Stop while the tool loop is between rounds.
        eng.request_stop()
        yield StreamEvent(type="text_delta", text="\nstopped, rounds kept\n")
        yield StreamEvent(type="message_delta", input_tokens=4000,
                          output_tokens=250, cost_usd=0.42,
                          stop_reason="stopped")

    eng.client.stream_message = MagicMock(side_effect=stopped_mid_turn)
    eng.stream_response("expensive task", on_token=notice.append)

    assert eng.cost_usd == pytest.approx(0.42), (
        "a stopped turn billed real money and reported $0.00 — the "
        "run-budget gate cannot see it")
    assert eng.token_usage["output"] == 250
    assert any("stopped" in n for n in notice), (
        "the client's own stop notice never reached the user")


def test_the_stop_notice_is_not_taken_for_the_models_answer(eng):
    """It goes to the caller, not into the conversation: the history must
    not gain an assistant message the model never wrote."""
    def stopped_mid_turn(**kw):
        eng.request_stop()
        yield StreamEvent(type="text_delta", text="stopped\n")
        yield StreamEvent(type="message_delta", output_tokens=1,
                          cost_usd=0.01, stop_reason="stopped")

    eng.client.stream_message = MagicMock(side_effect=stopped_mid_turn)
    eng.stream_response("hi")
    assert eng.messages == [], (
        "the unanswered message was left behind, or machinery text was "
        "stored as an answer")


def test_a_client_that_ignores_the_stop_is_still_cut_off(eng):
    """The drain is a budget, not a door: text that keeps coming after a
    stop must not keep being delivered."""
    delivered: list[str] = []

    def never_stops(**kw):
        eng.request_stop()
        for i in range(50):
            yield StreamEvent(type="text_delta", text=f"chunk{i} ")

    eng.client.stream_message = MagicMock(side_effect=never_stops)
    eng.stream_response("hi", on_token=delivered.append)
    assert len(delivered) <= 4, delivered


def test_a_normal_turn_is_unaffected(eng):
    def fine(**kw):
        yield StreamEvent(type="text_delta", text="the answer")
        yield StreamEvent(type="message_delta", input_tokens=100,
                          output_tokens=20, cost_usd=0.05)

    eng.client.stream_message = MagicMock(side_effect=fine)
    assert eng.stream_response("q") == "the answer"
    assert eng.cost_usd == pytest.approx(0.05)


# ---------------------------------------------------------------------------
# A stop that produced nothing must say that it produced nothing
# ---------------------------------------------------------------------------

def _stopped_before_any_text(engine):
    """Stop during thinking: the model streams reasoning and no answer."""
    def thinking_then_stop(*a, **kw):
        engine.request_stop()
        yield StreamEvent(type="thinking_delta", text="denke nach ...")
        yield StreamEvent(type="message_delta", output_tokens=0, cost_usd=0.0)
    engine.client.stream_message = MagicMock(side_effect=thinking_then_stop)


def test_a_stop_before_any_text_says_the_question_was_taken_back_out(eng):
    """Measured before the fix: the turn returned '', on_token was never
    called, and the history came back EMPTY -- the user's question removed
    without a word.

    Removing it is right; two consecutive user messages break the API on
    the next turn. Saying nothing is not: the stop notice the user does
    get reads "the rounds completed so far are kept; send a message to
    continue from here", which is true of a stop inside the tool loop and
    false here. Nothing was kept and there is no question left to continue
    from. The empty-turn branch next to this one already names what it
    took out; this one now does too."""
    seen: list[str] = []
    _stopped_before_any_text(eng)

    out = eng.stream_response("Berechne bitte die Bindungslängen",
                              on_token=seen.append)

    assert "send it again" in out
    assert "not in it" in out
    # ... and the user is told, not just the return value's caller.
    assert any("send it again" in s for s in seen)
    # The question really is out of the history -- the sentence is not a
    # consolation for something that did not happen.
    assert not any("Bindungsl" in str(m.get("content", ""))
                   for m in eng.messages)


def test_a_stop_after_real_text_keeps_the_answer_and_says_nothing_extra(eng):
    """The other half. A stop that interrupted a turn which HAD already
    answered keeps both sides of the exchange, and must not have the
    "send it again" sentence bolted onto a real answer."""
    def text_then_stop(*a, **kw):
        yield StreamEvent(type="text_delta", text="Teilantwort")
        eng.request_stop()
        yield StreamEvent(type="message_delta", output_tokens=3, cost_usd=0.0)
    eng.client.stream_message = MagicMock(side_effect=text_then_stop)

    out = eng.stream_response("frage zwei")

    assert out == "Teilantwort"
    assert "send it again" not in out
    assert [m["role"] for m in eng.messages] == ["user", "assistant"]


# ---------------------------------------------------------------------------
# The engine wins the race, so the engine owes the hand-back
# ---------------------------------------------------------------------------

def test_a_stop_hands_back_the_typed_message_the_client_never_saw(eng):
    """TRACED live, 2026-08-14, against kit.qwen3.5-397b-A17b:

        0.42s  turn starts
        4.16s  _stop_was_requested() -> False      (round 1, top)
       13.13s  TOOL DISPATCHED + stop requested
       13.14s  turn returned

    Ten milliseconds. The engine checks the stop flag BETWEEN STREAM
    EVENTS and breaks out there; the client's tool loop is a generator, and
    a consumer that stops iterating simply abandons it. So the client's own
    stop branch — the one that drains the queue and hands the typed message
    back — never executes. Three runs out of three left the message sitting
    in the queue, from where it would surface in a later, unrelated turn,
    out of order with whatever the user typed after pressing Stop.

    Whoever ends the turn owes the hand-back."""
    queue = ["und bitte auch die kleinste Zahl"]
    eng.client._drain_steer = lambda: [queue.pop(0)] if queue else []

    def thinking_then_stop(*a, **kw):
        eng.request_stop()
        yield StreamEvent(type="thinking_delta", text="denke nach ...")
        yield StreamEvent(type="message_delta", output_tokens=0, cost_usd=0.0)
    eng.client.stream_message = MagicMock(side_effect=thinking_then_stop)

    out = eng.stream_response("die ursprüngliche Frage")

    assert "kleinste Zahl" in out
    assert "send it again" in out
    # ... and it is out of the queue, so no later turn replays it.
    assert queue == []


def test_a_run_note_is_left_where_it_is(eng):
    """The other queue is NOT handed back: a run note is a fact about the
    session that is still true, the client outlives the turn, and the next
    turn drains it at its first round. Handing it to the user as their own
    words is the defect that was fixed earlier today."""
    eng.client._drain_steer = lambda: []
    notes = ["[permissions] switched to acceptEdits"]
    eng.client._drain_run_notes = lambda: list(notes)

    def thinking_then_stop(*a, **kw):
        eng.request_stop()
        yield StreamEvent(type="thinking_delta", text="...")
        yield StreamEvent(type="message_delta", output_tokens=0, cost_usd=0.0)
    eng.client.stream_message = MagicMock(side_effect=thinking_then_stop)

    out = eng.stream_response("frage")

    assert "[permissions]" not in out
    assert notes == ["[permissions] switched to acceptEdits"]


def test_nothing_pending_adds_no_sentence(eng):
    """The stop message must not grow a dangling clause when there is
    nothing to hand back."""
    eng.client._drain_steer = lambda: []

    def thinking_then_stop(*a, **kw):
        eng.request_stop()
        yield StreamEvent(type="thinking_delta", text="...")
        yield StreamEvent(type="message_delta", output_tokens=0, cost_usd=0.0)
    eng.client.stream_message = MagicMock(side_effect=thinking_then_stop)

    out = eng.stream_response("frage")
    assert "Not delivered" not in out


def test_a_drain_that_returns_something_odd_is_not_rendered(eng):
    """A client whose drain does not return a list of strings — an older
    backend, a stub, a test double — must not have its repr pasted into the
    user's answer. That has happened once already: MagicMock text reached a
    real transcript.

    Asserted with a bare MagicMock, which is exactly the shape that slips
    through a truthiness check: it is truthy, it slices, and it formats."""
    eng.client._drain_steer = MagicMock(return_value=MagicMock())

    def thinking_then_stop(*a, **kw):
        eng.request_stop()
        yield StreamEvent(type="thinking_delta", text="...")
        yield StreamEvent(type="message_delta", output_tokens=0, cost_usd=0.0)
    eng.client.stream_message = MagicMock(side_effect=thinking_then_stop)

    out = eng.stream_response("frage")

    assert "Mock" not in out
    assert "Not delivered" not in out
