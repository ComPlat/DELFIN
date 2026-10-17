"""A wait registered in one turn suppresses nudges and wakes on its event.

Phase 2, package KG1. What exists (verified): turn_continuation.should_continue
(turn_continuation.py:310) nudges on any announced action; the engine starts at
most one follow-up turn (engine.py:7880). The gap: a turn that WAITS — for a
message from a session, a job id, an approval id, a time — must not be nudged
("Let me know when X arrives" would be hand-back text today; "Let me check X
every minute" would be nudged forever). And the wait must WAKE the session when
its event arrives, with a one-line note naming the event.

The registry is pure + a small store: waits persist per session key in a
JSON file under the redirected user state (never a real ~/.delfin), so a
wake event arriving between turns can be matched against a wait registered
in an earlier turn.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import state_paths, turn_continuation
from delfin.agent.wait_registry import (
    WAIT_KINDS,
    SessionWaits,
    due,
    register,
    wait_note,
)


@pytest.fixture
def home(tmp_path):
    with state_paths.RedirectedUserState(tmp_path):
        yield tmp_path / ".delfin"


# -- kind vocabulary ---------------------------------------------------------

def test_kinds_are_the_briefs_kinds(home):
    assert set(WAIT_KINDS) == {"message_from", "job_id", "approval_id", "time"}


# -- registration ------------------------------------------------------------

def test_register_returns_a_wait_with_a_kind_and_a_value(home):
    w = register("message_from", "nacht-s25")
    assert w["kind"] == "message_from"
    assert w["value"] == "nacht-s25"


def test_register_refuses_an_unknown_kind(home):
    with pytest.raises(ValueError):
        register("weather", "rain")


def test_register_refuses_an_empty_value(home):
    with pytest.raises(ValueError):
        register("message_from", "  ")


def test_waits_persist_between_turns(home):
    w = register("message_from", "nacht-s25")
    again = SessionWaits().list()
    assert [x["id"] for x in again] == [w["id"]]


def test_a_wait_can_be_cleared_by_the_person(home):
    w = register("job_id", "4711")
    assert SessionWaits().clear(w["id"]) is True
    assert SessionWaits().list() == []
    assert SessionWaits().clear("no-such-id") is False


# -- due(now, events) --------------------------------------------------------

def test_due_matches_a_message_from_event(home):
    register("message_from", "nacht-s25")
    events = [("message_from", "nacht-s25")]
    hits = due(time.time(), events)
    assert [(h["kind"], h["value"]) for h in hits] == [("message_from", "nacht-s25")]


def test_due_matches_a_time_wait_that_has_passed(home):
    register("time", "0")  # epoch: in the past from the moment it is set
    hits = due(time.time(), [])
    assert [(h["kind"], h["value"]) for h in hits] == [("time", "0")]


def test_due_matches_a_job_id_event(home):
    register("job_id", "4711")
    assert [(h["kind"], h["value"]) for h in due(time.time(), [("job_id", "4711")])] == \
        [("job_id", "4711")]


def test_due_ignores_unrelated_events(home):
    register("message_from", "nacht-s25")
    hits = due(time.time(), [("message_from", "nacht-s13"),
                             ("job_id", "4711")])
    assert hits == []


def test_due_does_not_match_another_message_sender(home):
    register("message_from", "nacht-s25")
    assert due(time.time(), [("message_from", "NACHT-S25")]) == [], \
        "sender matching is exact, not case-folded: a wait names one sender"


def test_due_consumes_what_it_fires(home):
    register("message_from", "nacht-s25")
    first = due(time.time(), [("message_from", "nacht-s25")])
    assert first
    assert due(time.time(), [("message_from", "nacht-s25")]) == [], \
        "a fired wait is consumed: the wake happens once, not every round"


def test_due_fires_an_event_wait_only_once(home):
    register("job_id", "4711")
    assert due(time.time(), [("job_id", "4711")])
    assert due(time.time(), [("job_id", "4711")]) == []
    # but the same id arriving AGAIN after the wait was re-registered fires again
    register("job_id", "4711")
    assert due(time.time(), [("job_id", "after-first")]) == []
    assert due(time.time(), [("job_id", "4711")])


def test_due_fires_a_time_wait_only_once(home):
    register("time", "0")
    assert due(time.time(), [])
    assert due(time.time(), []) == []


def test_due_orders_time_waits_before_event_waits(home):
    register("job_id", "4711")
    register("time", "0")
    hits = due(time.time(), [("job_id", "4711")])
    assert [h["kind"] for h in hits] == ["time", "job_id"]


# -- suppression of nudges while a wait is registered ------------------------

def test_a_wait_suppresses_the_continuation_nudge(home):
    cont, _ = turn_continuation.should_continue(
        "Let me write the tests now.", 0, session_key="kg1-s20")
    assert cont is True
    register("message_from", "nacht-s25", key="kg1-s20")
    cont, note = turn_continuation.should_continue(
        "Let me write the tests now.", 0, session_key="kg1-s20")
    assert cont is False
    assert "wait" in note.lower()


def test_suppression_ends_when_the_wait_is_consumed(home):
    register("message_from", "nacht-s25", key="kg1-s20")
    due(time.time(), [("message_from", "nacht-s25")], key="kg1-s20")
    cont, note = turn_continuation.should_continue(
        "Let me write the tests now.", 0, session_key="kg1-s20")
    assert cont is True
    assert note  # the ordinary continuation note returns


def test_an_empty_registry_suppresses_nothing(home):
    cont, note = turn_continuation.should_continue(
        "Let me write the tests now.", 0, session_key="kg1-s20")
    assert cont is True


def test_a_caller_without_a_key_suppresses_nothing(home):
    # session_key=None (the engine's current call shape): no suppression,
    # so the wiring change cannot turn a nudge off by accident.
    register("message_from", "nacht-s25", key="kg1-s20")
    cont, _ = turn_continuation.should_continue(
        "Let me write the tests now.", 0)
    assert cont is True


# -- one session's wait never suppresses another's (finding #4) --------------

def test_a_wait_is_scoped_to_its_session(home):
    register("message_from", "nacht-s25", key="kg1-s20")
    cont, _ = turn_continuation.should_continue(
        "Let me write the tests now.", 0, session_key="kg2-s13")
    assert cont is True, "another session's wait must not suppress this one"
    assert due(time.time(), [("message_from", "nacht-s25")], key="kg2-s13") == []
    assert due(time.time(), [("message_from", "nacht-s25")], key="kg1-s20")


# -- the wake note names the event -------------------------------------------

def test_wait_note_names_kind_and_value(home):
    w = register("message_from", "nacht-s25")
    note = wait_note(w)
    assert "message" in note
    assert "nacht-s25" in note


def test_wait_note_for_a_time_wait(home):
    w = register("time", "0")
    note = wait_note(w)
    assert "wait" in note
