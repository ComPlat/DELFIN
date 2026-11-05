"""A question without an answer says why; a typed phone line answers it.

Field report (Pi): ask_user_question returned ``{"answers": []}`` to the
model, which then guessed. Two causes, measured:

* the broker's refusal of a question (expiry, abort, a spent code) was a
  bare empty list, the same value as a person dismissing it;
* a plain line typed on the phone while a question waited became a new
  prompt for the blocked turn instead of the answer, and the question
  expired.

Output, each with its control:

* expiry: ``timed_out`` plus ``unanswered`` naming the seconds waited;
  abort: ``unanswered`` naming the stop; the tool result carries the
  reason and a note that the empty list is not a choice;
* a dismissal (an answer with no labels) is named as dismissed;
* a plain phone line in the session chat while a question waits is its
  answer (an option label or free text); with an approval waiting instead,
  the line is still a prompt;
* the same for a session outside the service, through its room.
"""
from __future__ import annotations

import asyncio
import json

import pytest

import test_one_answer_closes_a_request_everywhere as one
from delfin.agent import terminal_confirm as tc
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


QUESTION = one.QUESTION


@pytest.fixture
def acct(tmp_path, monkeypatch):
    return one.make_account(tmp_path, monkeypatch)


def test_an_expired_question_says_it_timed_out():
    broker = tc.TerminalConfirmBroker(timeout_s=0.05, poll_s=0.01)
    result = broker.ask_user(QUESTION)
    assert result["answers"] == [] and result.get("timed_out") is True
    assert "timed out after" in result.get("unanswered", "")


def test_a_stopped_question_says_it_was_stopped():
    broker = tc.TerminalConfirmBroker(timeout_s=5, poll_s=0.01)
    th, box = one._ask_in_thread(broker.ask_user, QUESTION)
    assert one._wait_for(lambda: broker.take() is not None)
    broker.abort_all()
    th.join(5)
    assert box["result"]["answers"] == []
    assert "stopped" in box["result"].get("unanswered", "")
    assert not box["result"].get("timed_out")


def _tool(tmp_path, result):
    perms = KitToolPermissions(workspace=tmp_path, mode="default")
    perms.ask_user_callback = lambda payload: result
    return json.loads(_DocToolExecutor()._execute_ask_user_question({
        "question": "Which solvent?",
        "options": [{"label": "water"}, {"label": "toluene"}]}, perms))


def test_the_tool_result_names_the_reason_and_says_it_is_no_choice(tmp_path):
    out = _tool(tmp_path, {"answers": [], "unanswered": "no answer (the question was stopped)"})
    assert out["answers"] == [] and "stopped" in out["unanswered"]
    assert "not a choice" in out["note"].lower()
    timed = _tool(tmp_path, {"answers": [], "timed_out": True})
    assert "timed out" in timed["unanswered"]


def test_a_dismissal_is_named_as_dismissed(tmp_path):
    out = _tool(tmp_path, {"answers": []})
    assert "dismissed" in out["unanswered"]
    control = _tool(tmp_path, {"answers": ["water"]})
    assert control["answers"] == ["water"] and "unanswered" not in control


def test_a_plain_phone_line_answers_the_waiting_question(acct):
    broker = acct.service.broker
    th, box = one._ask_in_thread(broker.ask_user, QUESTION)
    assert one._wait_for(lambda: acct.service.pending)
    acct.say("21", "#a1", "toluene")
    th.join(5)
    assert box["result"]["answers"] == ["toluene"]
    assert not [p for m, p in acct.service.asked if m == "session.prompt"]


def test_free_text_on_the_phone_is_a_free_answer(acct):
    broker = acct.service.broker
    th, box = one._ask_in_thread(broker.ask_user, QUESTION)
    assert one._wait_for(lambda: acct.service.pending)
    acct.say("21", "#a1", "use ethanol instead")
    th.join(5)
    assert box["result"]["answers"] == ["use ethanol instead"]


def test_with_an_approval_waiting_a_plain_line_is_still_a_prompt(acct):
    broker = acct.service.broker
    th, box = one._ask_in_thread(broker.callback, "write_file", {"path": "n.txt"}, "n.txt: +1")
    assert one._wait_for(lambda: acct.service.pending)
    acct.say("21", "#a1", "also add a title")
    assert acct.service.pending, "a plain line never approves"
    assert [p for m, p in acct.service.asked if m == "session.prompt"]
    broker.abort_all()
    th.join(5)
    assert box["result"] is False


def test_a_plain_line_answers_a_question_of_a_session_outside_the_service(acct):
    from delfin.messenger import activation
    activation.configure("s-repl", directory=acct.activation_dir("anna"))
    repl = tc.TerminalConfirmBroker(session_id="s-repl", poll_s=0.02)
    th, box = one._ask_in_thread(repl.ask_user, QUESTION)
    assert one._wait_for(lambda: list(tc._PENDING_DIR.glob("*.request.json")))
    asyncio.run(one._refresh(acct))
    acct.say("21", "@21", "water")
    th.join(5)
    assert box["result"]["answers"] == ["water"]
