"""A session re-opened is the same session, at the same address.

Three sessions worked together and one of them was answered "no other
open session '<key>'" for a peer it had reached twice minutes earlier.
The peer had been re-opened: same saved conversation, same workspace,
same title -- and a different key, because the dashboard minted
``uuid.uuid4().hex[:8]`` per TAB rather than per session. The sender
recovered by reading the roster out of the refusal and sending again, a
round each time, and anything still unread in the old mailbox stayed
there.

``cli._presence_key_for`` already derives a terminal session's address
from the head of its session id. Two surfaces held two schemes for one
idea, and the unstable one was the one the dashboard used.

The other half is that the roster never showed the stable address.
``session_message`` accepts ``to=<session_id>``, but the rows it hands
back on a refusal carried key, title, workspace and branch -- so a model
whose address had just gone stale was given the next unstable one and
nothing else.
"""

from __future__ import annotations

import json

import pytest


class TestTheAddress:
    """`_address_for`, built the way the tab builds it."""

    @staticmethod
    def _address_for(open_keys):
        """The helper as a closure over a session list, as in the module."""
        import uuid

        sessions = [{"key": k} for k in open_keys]

        def _address(session_id: str) -> str:
            head = str(session_id or "").strip()[:8]
            if head and not any(rec["key"] == head for rec in sessions):
                return head
            return uuid.uuid4().hex[:8]
        return _address

    def test_the_same_conversation_gets_the_same_address(self):
        sid = "7cdf74d4fdda4a539b1fa3e8e60a7067"
        first = self._address_for([])(sid)
        second = self._address_for([])(sid)
        assert first == second == sid[:8]

    def test_a_conversation_with_no_id_yet_gets_a_fresh_one(self):
        a = self._address_for([])("")
        b = self._address_for([])("")
        assert a != b and len(a) == 8

    def test_two_tabs_never_share_one_inbox(self):
        """Being wrong this way would deliver one session's messages to
        another, so a collision takes a random address instead."""
        sid = "7cdf74d4fdda4a539b1fa3e8e60a7067"
        taken = self._address_for([sid[:8]])(sid)
        assert taken != sid[:8] and len(taken) == 8

    def test_the_module_really_uses_it(self):
        """Read from the module, because a helper the tab does not call
        is not a fix."""
        import inspect
        from delfin.dashboard import agent_sessions as AS

        body = inspect.getsource(AS)
        assert "key = _address_for(sid)" in body, (
            "the tab still mints a per-tab uuid")
        assert "def _address_for(" in body

    def test_it_agrees_with_the_terminal_scheme(self):
        """One word for one session across both surfaces."""
        from delfin.agent.cli import _presence_key_for

        sid = "7cdf74d4fdda4a539b1fa3e8e60a7067"
        assert _presence_key_for("", sid) == self._address_for([])(sid)


class TestTheRosterShowsTheStableAddress:
    @pytest.fixture
    def peers(self, monkeypatch):
        from delfin.agent import api_client as A
        from delfin.agent import session_presence as P

        rows = [{"key": "6a843aca", "session_id": "c" * 32,
                 "title": "Session C", "workspace": "/w/c",
                 "branch": "session/c"}]
        # The tool imports session_presence inside the call, so the module
        # is where the roster has to come from.
        monkeypatch.setattr(P, "open_sessions", lambda **k: list(rows))
        return A, rows

    def _send(self, A, to, message="hi"):
        eng = A._DocToolExecutor.__new__(A._DocToolExecutor)
        perms = type("_P", (), {"presence_key": "mine",
                                "workspace": "/w/b"})()
        eng._permissions = perms
        return json.loads(eng._execute_session_message(
            {"to": to, "message": message}, perms))

    def test_a_refusal_hands_back_the_session_id(self, peers):
        A, _rows = peers
        out = self._send(A, "3eb6c29f")
        assert "no other open session" in out["error"]
        assert out["sessions"][0]["session_id"] == "c" * 32, out
        assert "session_id" in out["error"], (
            "the refusal has to name the address that does not go stale")

    def test_the_listing_hands_it_back_too(self, peers):
        """`session_message` with no `to` is how a session learns the
        roster in the first place."""
        A, _rows = peers
        eng = A._DocToolExecutor.__new__(A._DocToolExecutor)
        perms = type("_P", (), {"presence_key": "mine",
                                "workspace": "/w/b"})()
        eng._permissions = perms
        out = json.loads(eng._execute_session_message({}, perms))
        assert out["sessions"][0]["session_id"] == "c" * 32, out
        assert "session_id" in out["note"], out

    def test_the_session_id_actually_works_as_an_address(self, peers):
        A, _rows = peers
        out = self._send(A, "c" * 32)
        assert out.get("status") == "sent", out
        assert out.get("to") == "6a843aca", (
            "resolved to the live key, so the message reaches the inbox "
            "that session reads")
