"""Telling three peers the same thing takes one call, not three.

Of 541 tool calls across twelve sessions in one afternoon, 71 were
session_message — and a good part of that was the same sentence sent to
one peer and then to the next, because a message could only be addressed
to a single session. A session that addressed a peer by its TITLE was
told the key was unknown and had to make a second call just to be handed
the roster.

  to="all"                 one call reaches every other open session
  the roster comes back     with the refusal, not after another call
  a key still works         the ordinary case is unchanged
  nothing to say, nothing   no peers, or no text, is an error not a send
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.api_client import _doc_executor


class _Presence:
    def __init__(self, peers):
        self._peers = peers

    def open_sessions(self, exclude_key=None):
        return [p for p in self._peers if p.get("key") != exclude_key]


class _Messages:
    def __init__(self):
        self.sent = []

    def send(self, key, text, from_key="", from_title=""):
        self.sent.append((key, text))


@pytest.fixture
def wired(monkeypatch):
    peers = [
        {"key": "me", "title": "Session A", "session_id": "a1"},
        {"key": "b2", "title": "Session B — rendering", "session_id": "s-b"},
        {"key": "c3", "title": "Session C — command", "session_id": "s-c"},
    ]
    presence = _Presence(peers)
    messages = _Messages()
    import delfin.agent.session_presence as sp
    import delfin.agent.session_messages as sm
    monkeypatch.setattr(sp, "open_sessions", presence.open_sessions)
    monkeypatch.setattr(sm, "send", messages.send)
    return messages


def _call(args, key="me"):
    class _Perms:
        presence_key = key

    return json.loads(_doc_executor._execute_session_message(args, _Perms()))


def test_one_call_reaches_every_other_session(wired):
    out = _call({"to": "all", "message": "the contract is collect/format/main"})
    assert out["status"] == "sent"
    assert sorted(out["to"]) == ["b2", "c3"]
    assert [k for k, _t in wired.sent] == ["b2", "c3"]


def test_a_single_key_still_works(wired):
    out = _call({"to": "b2", "message": "yours is final?"})
    assert out["status"] == "sent" and out["to"] == "b2"
    assert wired.sent == [("b2", "yours is final?")]


def test_an_unknown_name_is_answered_with_the_roster(wired):
    out = _call({"to": "Session B — rendering", "message": "hi"})
    assert "error" in out
    assert [s["key"] for s in out["sessions"]] == ["b2", "c3"]
    assert wired.sent == []


def test_a_broadcast_with_nothing_to_say_sends_nothing(wired):
    assert "error" in _call({"to": "all", "message": "   "})
    assert wired.sent == []
