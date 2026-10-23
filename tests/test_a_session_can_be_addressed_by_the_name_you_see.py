"""Writing to a peer by the name on screen cost two rounds, every time.

``session_message`` resolved ``to`` against the session KEY only --- an
opaque eight-hex handle that appears nowhere in a conversation. The roster
the tool hands back carries the key and the title side by side, and the
title is what a sender has actually been told, so a session asked to
write to "Session A" addressed it by that name, was refused, and reissued
the identical message with the key.

Observed across three sessions in one afternoon (2026-10-07): two refused
sends, then two successful ones carrying the same text. The refusal was
already helpful --- it returns the roster rather than making the sender ask
for it --- so the cost was a round, not a lost message. A round each way,
on every first contact.

An exact, case-folded title is now accepted when exactly ONE open session
carries it. Not more than one: titles are the first line of a conversation
and two sessions can easily share one, and a message in the wrong inbox is
worse than the extra round. An ambiguous name is refused by name, with the
candidates.
"""

from __future__ import annotations

import json

import pytest


class _Perms:
    def __init__(self, key):
        self.presence_key = key


@pytest.fixture
def executor():
    from delfin.agent.api_client import _DocToolExecutor

    return object.__new__(_DocToolExecutor)


@pytest.fixture
def peers(monkeypatch):
    """Two open peers, and a mailbox that accepts and records."""
    import delfin.agent.session_messages as msgs
    import delfin.agent.session_presence as pres

    rows = [
        {"key": "2a15fb71", "title": "Hallo Session A",
         "workspace": "/w/a", "branch": "session/a"},
        {"key": "f2100819", "title": "Hallo Session B",
         "workspace": "/w/b", "branch": "session/b"},
    ]
    sent = []

    def _open(**kw):
        if kw.get("exclude_key"):
            return [r for r in rows if r["key"] != kw["exclude_key"]]
        return rows

    monkeypatch.setattr(pres, "open_sessions", _open)
    monkeypatch.setattr(msgs, "deliverable", lambda key: False)
    monkeypatch.setattr(
        msgs, "send",
        lambda to, text, **kw: sent.append((to, text)))
    return sent


def _send(executor, to, text="hi", me="3eb6c29f"):
    return json.loads(executor._execute_session_message(
        {"to": to, "message": text}, _Perms(me)))


# ---------------------------------------------------------------------------
# The name on screen works
# ---------------------------------------------------------------------------

def test_a_unique_title_is_delivered(executor, peers):
    out = _send(executor, "Hallo Session A")
    assert out.get("status") == "sent"
    assert out["to"] == "2a15fb71"
    assert peers == [("2a15fb71", "hi")]


def test_the_key_still_works(executor, peers):
    out = _send(executor, "f2100819")
    assert out.get("status") == "sent"
    assert peers == [("f2100819", "hi")]


def test_the_title_is_matched_case_insensitively(executor, peers):
    out = _send(executor, "hallo session a")
    assert out.get("status") == "sent"
    assert out["to"] == "2a15fb71"


def test_surrounding_space_is_not_a_different_name(executor, peers):
    out = _send(executor, "  Hallo Session A  ")
    assert out.get("status") == "sent"
    assert out["to"] == "2a15fb71"


def test_the_reply_names_the_title_so_the_sender_can_check_it(executor, peers):
    out = _send(executor, "Hallo Session B")
    assert out["title"] == "Hallo Session B"


# ---------------------------------------------------------------------------
# An ambiguous name is refused, not guessed
# ---------------------------------------------------------------------------

def test_two_sessions_with_one_name_are_not_guessed_between(
        executor, monkeypatch):
    import delfin.agent.session_messages as msgs
    import delfin.agent.session_presence as pres

    rows = [
        {"key": "aaaa1111", "title": "continue", "workspace": "/w/1"},
        {"key": "bbbb2222", "title": "continue", "workspace": "/w/2"},
    ]
    sent = []
    monkeypatch.setattr(pres, "open_sessions",
                        lambda **kw: [r for r in rows
                                      if r["key"] != kw.get("exclude_key")])
    monkeypatch.setattr(msgs, "deliverable", lambda key: False)
    monkeypatch.setattr(msgs, "send",
                        lambda to, text, **kw: sent.append((to, text)))

    out = _send(executor, "continue")
    assert "error" in out
    assert sent == [], "a message was delivered to a guess"
    assert "2 open sessions" in out["error"]
    keys = {row["key"] for row in out["sessions"]}
    assert keys == {"aaaa1111", "bbbb2222"}, "the candidates are not named"


def test_an_unknown_name_still_hands_back_the_roster(executor, peers):
    out = _send(executor, "Session Z")
    assert "error" in out
    assert {r["key"] for r in out["sessions"]} == {"2a15fb71", "f2100819"}
    assert peers == []


def test_the_refusal_no_longer_says_titles_do_not_work(executor, peers):
    out = _send(executor, "Session Z")
    assert "not by its title" not in out["error"]
    assert "title" in out["error"]


# ---------------------------------------------------------------------------
# Nothing else about the tool moved
# ---------------------------------------------------------------------------

def test_no_recipient_still_lists(executor, peers):
    out = json.loads(executor._execute_session_message(
        {}, _Perms("3eb6c29f")))
    assert {r["key"] for r in out["sessions"]} == {"2a15fb71", "f2100819"}
    assert peers == []


def test_to_all_still_reaches_everyone(executor, peers):
    out = json.loads(executor._execute_session_message(
        {"to": "all", "message": "x"}, _Perms("3eb6c29f")))
    assert out["status"] == "sent"
    assert {k for k, _ in peers} == {"2a15fb71", "f2100819"}


def test_a_message_is_still_required(executor, peers):
    out = _send(executor, "Hallo Session A", text="")
    assert "error" in out
    assert peers == []


def test_a_headless_turn_is_still_refused(executor, peers):
    out = json.loads(executor._execute_session_message(
        {"to": "Hallo Session A", "message": "x"}, _Perms("")))
    assert "only inside an open dashboard session" in out["error"]
    assert peers == []


def test_the_tool_description_says_a_title_works():
    """Discoverable without a refusal first --- the schema is what the model
    reads before it chooses an argument."""
    from delfin.agent.api_client import _DOC_TOOLS_OPENAI

    desc = next(t["function"]["description"] for t in _DOC_TOOLS_OPENAI
                if t["function"]["name"] == "session_message")
    assert "title" in desc
    assert "key" in desc
