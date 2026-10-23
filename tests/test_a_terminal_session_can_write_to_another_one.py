"""A terminal session can write to another one, not only be written to.

Half of this landed already: a terminal session announces its presence
and takes what is left in its inbox. The other half was still refused.
``session_message`` reads ``perms.presence_key``, the dashboard sets it
and the CLI never did, so the tool answered

    session_message is available only inside an open dashboard session.

Seen for real on 2026-09-20: a session noticed it was editing files a
parallel session had left in the tree, went to say so, and was told it
does not exist. It carried on and committed around the other's work,
which was the right call -- but it made it without being able to ask.

The key is the one the operator's inbox already uses, so a session is
addressed the same way by a person and by another session.
"""

from __future__ import annotations

import pytest


class _Perms:
    def __init__(self, key=""):
        self.presence_key = key


def _key_for(name, session_id):
    from delfin.agent.cli import _presence_key_for
    return _presence_key_for(name, session_id)


class TestTheKey:
    def test_the_name_wins(self):
        assert _key_for("runde2-s1", "abcdef0123") == "runde2-s1"

    def test_without_a_name_the_session_id_serves(self):
        assert _key_for("", "abcdef0123") == "abcdef01"

    def test_nothing_at_all_is_empty_not_a_crash(self):
        assert _key_for("", "") == ""

    def test_whitespace_is_not_a_name(self):
        assert _key_for("   ", "abcdef0123") == "abcdef01"


class TestTheToolNoLongerRefuses:
    @pytest.fixture
    def executor(self):
        from delfin.agent.api_client import _DocToolExecutor
        return object.__new__(_DocToolExecutor)

    def test_a_session_with_a_key_is_allowed_to_look(self, executor, monkeypatch):
        import delfin.agent.session_presence as pres
        monkeypatch.setattr(pres, "open_sessions", lambda **k: [])
        out = executor._execute_session_message({}, _Perms("runde2-s1"))
        assert "only inside an open dashboard session" not in out

    def test_without_a_key_it_still_refuses(self, executor):
        out = executor._execute_session_message({}, _Perms(""))
        assert "only inside an open dashboard session" in out

    def test_it_lists_the_others(self, executor, monkeypatch):
        import delfin.agent.session_presence as pres
        monkeypatch.setattr(pres, "open_sessions",
                            lambda **k: [{"key": "runde2-s2", "title": "t",
                                          "workspace": "/w", "branch": "b"}])
        out = executor._execute_session_message({}, _Perms("runde2-s1"))
        assert "runde2-s2" in out
