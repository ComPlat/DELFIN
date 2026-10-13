"""A session that ends three turns in a row on the same dead call stops.

The model sometimes writes a tool call as text. The engine asks again once
per turn; then the harness sends a follow-up turn on its own, the model
writes the call as text again, and so on -- sessions sat in that loop for
hours. Nothing inside one model call can see it: each call is one round.
The engine counts how consecutive turns ENDED and, at three identical dead
endings with no tool run and no message from a person in between, says so
and sends no further follow-up.

Driven through ``AgentEngine.stream_response`` with a client that leaks.
"""
from __future__ import annotations

from unittest.mock import patch

import pytest

from delfin.agent import engine as E
from delfin.agent.api_client import StreamEvent


def _leak(tool="write_file", prose="Let me write the file."):
    return (f'{prose}\n<invoke name="{tool}">\n<parameter name="path">a.txt'
            f'</parameter>\n<parameter name="content">x</parameter>\n</invoke>')


class _Client:
    """Answers every call with the next scripted text (the last one forever)."""

    def __init__(self, *answers):
        self.answers = list(answers)
        self.calls = 0

    def stream_message(self, *_a, **_k):
        self.calls += 1
        text = self.answers[min(self.calls - 1, len(self.answers) - 1)]
        yield StreamEvent(type="text_delta", text=text)
        yield StreamEvent(type="message_delta", stop_reason="end_turn",
                          input_tokens=1, output_tokens=1)


@pytest.fixture
def make(tmp_path):
    def _make(*answers):
        client = _Client(*answers)
        with patch("delfin.agent.engine.create_client", return_value=client):
            eng = E.AgentEngine(repo_dir=str(tmp_path), backend="api",
                                provider="kit", model="kit.qwen3.5-397b-A17b",
                                mode="solo")
        return eng, client
    return _make


def test_three_dead_turns_in_a_row_stop_with_a_sentence_naming_the_tool(make):
    eng, _ = make(_leak())
    eng.stream_response("go")
    eng.stream_response("harness follow-up 1")
    out = eng.stream_response("harness follow-up 2")
    assert eng._stall_guard_stopped
    assert "write_file" in eng._stall_guard_stopped
    assert "written as text" in out and "3 turns in a row" in out


def test_two_dead_turns_do_not_stop(make):
    eng, _ = make(_leak())
    eng.stream_response("go")
    eng.stream_response("harness follow-up 1")
    assert eng._stall_guard_stopped == ""


def test_different_words_around_the_same_dead_call_are_the_same_ending(make):
    """The model never repeats itself word for word; the call that does not
    run is what repeats."""
    eng, _ = make(_leak(prose="Let me write the file."),
                  _leak(prose="Now I will create a.txt."),
                  _leak(prose="Writing it now, for real this time."),
                  _leak(prose="One more try."),
                  _leak(prose="Again."), _leak(prose="And again."))
    for i in range(3):
        eng.stream_response(f"turn {i}")
    assert eng._stall_guard_stopped


def test_the_third_dead_turn_is_not_asked_again(make):
    """Each of the first two turns re-requests once (two model calls); the
    stopping turn must not spend another call on the same request."""
    eng, client = make(_leak())
    eng.stream_response("go")
    eng.stream_response("f1")
    before = client.calls
    eng.stream_response("f2")
    assert client.calls == before + 1


def test_a_stopped_session_produces_no_follow_up_note(make):
    eng, _ = make(_leak())
    eng._open_task_count = lambda: 1
    eng._note_turn_continuation("Let me run the tests now.")
    assert eng.pending_turn_continuation, "control: the note exists before a stop"
    eng.clear_turn_continuation()
    for i in range(3):
        eng.stream_response(f"turn {i}")
    assert eng._stall_guard_stopped
    eng.clear_turn_continuation()
    eng._note_turn_continuation("Let me run the tests now.")
    assert eng.pending_turn_continuation == ""


def test_a_message_from_a_person_starts_the_count_again(make):
    eng, _ = make(_leak())
    eng.stream_response("go")
    eng.stream_response("f1")
    eng.note_typed("please continue")
    eng.stream_response("please continue")
    assert eng._stall_guard_stopped == ""


def test_a_message_from_a_person_lifts_a_stop(make):
    eng, _ = make(_leak())
    for i in range(3):
        eng.stream_response(f"turn {i}")
    assert eng._stall_guard_stopped
    eng.note_typed("go on")
    assert eng._stall_guard_stopped == ""
    eng._open_task_count = lambda: 1
    eng._note_turn_continuation("Let me run the tests now.")
    assert eng.pending_turn_continuation


def test_a_follow_up_the_harness_sent_does_not_count_as_a_person(make):
    """A follow-up turn enters through the same door as a typed one. If it
    reset the count, three in a row could never be reached."""
    eng, _ = make(_leak())
    for i in range(3):
        eng.stream_response("System note: this turn ended by announcing an "
                            "action and did not take it. Take it now.")
    assert eng._stall_guard_stopped


def test_an_ordinary_answer_between_dead_turns_starts_the_count_again(make):
    eng, _ = make(_leak(), _leak(), _leak(), _leak(),
                  "Here is the plain answer you asked for.",
                  _leak(), _leak(), _leak(), _leak())
    eng.stream_response("t0")          # leak + re-request leak
    eng.stream_response("t1")          # leak + re-request leak
    eng.stream_response("t2")          # a real answer
    assert eng._stall_guard_stopped == ""
    eng.stream_response("t3")
    eng.stream_response("t4")
    assert eng._stall_guard_stopped == "", (
        "two dead turns after a real answer are two, not four")


def test_every_engine_counts_for_itself(make):
    one, _ = make(_leak())
    two, _ = make(_leak())
    one.stream_response("a"); two.stream_response("a")
    one.stream_response("b"); two.stream_response("b")
    one.stream_response("c")
    assert one._stall_guard_stopped and two._stall_guard_stopped == ""


def test_a_turn_whose_second_request_is_answered_is_not_a_dead_turn(tmp_path):
    """The engine asks once more inside the turn. When that request is
    answered properly, the turn ended well: five such turns in a row are
    five good turns, each with its one extra request."""

    class _LeaksThenAnswers:
        calls = 0

        def stream_message(self, system=None, messages=None, *_a, **_k):
            type(self).calls += 1
            last_user = ""
            for m in messages or []:
                if m.get("role") == "user":
                    last_user = str(m.get("content", ""))
            text = ("Done, the file is written."
                    if "expressed as literal markup" in last_user else _leak())
            yield StreamEvent(type="text_delta", text=text)
            yield StreamEvent(type="message_delta", stop_reason="end_turn",
                              input_tokens=1, output_tokens=1)

    client = _LeaksThenAnswers()
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = E.AgentEngine(repo_dir=str(tmp_path), backend="api",
                            provider="kit", model="kit.qwen3.5-397b-A17b",
                            mode="solo")
    for i in range(5):
        out = eng.stream_response(f"turn {i}")
        assert "Done, the file is written." in out
    assert eng._stall_guard_stopped == ""
    assert _LeaksThenAnswers.calls == 10


def test_an_engine_that_never_ran_its_constructor_still_notes_a_follow_up():
    eng = object.__new__(E.AgentEngine)
    eng._continuation_fired = False
    eng._open_task_count = lambda: 1
    eng.pending_turn_continuation = ""
    eng._note_turn_continuation("Let me run the tests now.")   # must not raise


def test_the_stop_reaches_the_screen_and_the_history(make):
    """A surface prints what arrives through on_token (a terminal cannot
    print the return value afterwards), and the next turn's model reads
    the history: the stop sentence must be in both."""
    eng, _ = make(_leak())
    eng.stream_response("go")
    eng.stream_response("f1")
    streamed: list[str] = []
    eng.stream_response("f2", on_token=streamed.append)
    assert "stops asking" in "".join(streamed)
    last = [m for m in eng.messages if m.get("role") == "assistant"][-1]
    assert "stops asking" in str(last.get("content"))


def test_the_turn_that_stops_leaves_no_follow_up_behind(make):
    """The follow-up note is recorded before the dead-call check runs; the
    stopping turn must not leave it pending for the surface to send."""
    eng, _ = make(_leak())
    eng._open_task_count = lambda: 1
    eng.stream_response("go")
    assert eng.pending_turn_continuation, "control: a dead turn leaves a note"
    eng.clear_turn_continuation()
    eng.stream_response("f1")
    eng.clear_turn_continuation()
    eng.stream_response("f2")
    assert eng._stall_guard_stopped
    assert eng.pending_turn_continuation == ""


def test_a_new_conversation_or_a_loaded_session_starts_without_a_stop(make):
    """The dashboard's new cycle and session load reuse the engine object;
    the count belongs to the conversation that produced it."""
    eng, client = make(_leak())
    for i in range(3):
        eng.stream_response(f"t{i}")
    saved = eng.export_state()
    eng.reset_cycle()
    assert eng._stall_guard_stopped == "" and eng._stall_guard.streak == 0
    for i in range(3):
        eng.stream_response(f"t{i}")
    eng.restore_state(saved)
    assert eng._stall_guard_stopped == "" and eng._stall_guard.streak == 0


def test_a_line_typed_in_the_terminal_lifts_a_stop(make):
    """Driven through the terminal's own loop: the typed line starts a
    dead loop of follow-ups that stops on the third turn; no fourth
    follow-up goes out, and the next typed line gets its re-request."""
    import io

    from delfin.agent import repl, repl_render as rr

    class _Tty(io.StringIO):
        def isatty(self):
            return True

    eng, client = make(_leak())
    out = io.StringIO()
    agent = repl.TerminalAgent(eng, out=out, err=_Tty())
    agent.transcript.theme = rr.Theme(enabled=False)
    reads = iter(["write a.txt", "call write_file as a tool", "/exit"])
    agent._read_line = lambda _p: next(reads)
    calls_per_turn: list[int] = []
    real_turn = agent.turn

    def turn(prompt):
        before = client.calls
        result = real_turn(prompt)
        calls_per_turn.append(client.calls - before)
        return result

    agent.turn = turn
    agent.run()
    assert calls_per_turn == [2, 2, 1, 2], calls_per_turn
    assert "stops asking" in out.getvalue()
