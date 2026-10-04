"""T3 phase 2 — red control for the protected api_client change.

The api_client's _execute_session_message refuses any recipient that is not
among the open sessions ("no other open session"). This sharpens that: a
recipient that _msgs.deliverable() accepts (the reserved `operator` mailbox,
or a known but closed session) must be queued, not refused. This test is red
on the current api_client and goes green once .gate/t3_deliverable.patch is
built by the operator.
"""

from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from delfin.agent import session_messages as M
from delfin.agent import session_presence as P
from delfin.agent.api_client import _DocToolExecutor


@pytest.fixture(autouse=True)
def _dirs(tmp_path, monkeypatch):
    monkeypatch.setattr(M, "_DIR", tmp_path / "inbox")
    monkeypatch.setattr(P, "_DIR", tmp_path / "presence")
    P._last_written.clear()
    P._git_cache.clear()


def _call(perms, **arguments):
    return json.loads(_DocToolExecutor()._execute_session_message(arguments, perms))


def test_a_message_to_the_operator_is_queued_not_refused(tmp_path):
    """The operator needs no presence: a message to it is queued for its
    reserved mailbox instead of being refused with 'no other open session'."""
    perms = SimpleNamespace(presence_key="nacht-s16")
    # No session is open (no presence announces the operator).
    assert P.open_sessions() == []
    sent = _call(perms, to="operator", message="T3 phase 2 queued for you")
    assert sent["status"] == "sent"
    (message,) = M.take("operator")
    assert message["text"] == "T3 phase 2 queued for you"


def test_a_message_to_a_closed_session_is_queued_not_refused(tmp_path):
    """A session that was open but has gone stale is queued too, not refused."""
    P.announce("nacht-s12", title="builder", workspace=str(tmp_path))
    record_path = P._path("nacht-s12")
    old = json.loads(record_path.read_text(encoding="utf-8"))
    old["updated_at"] = 0.0
    record_path.write_text(json.dumps(old), encoding="utf-8")
    perms = SimpleNamespace(presence_key="nacht-s16")
    sent = _call(perms, to="nacht-s12", message="wake up on your next start")
    assert sent["status"] == "sent"
    (message,) = M.take("nacht-s12")
    assert message["text"] == "wake up on your next start"


def test_gibberish_is_still_refused(tmp_path):
    """An address that is neither reserved nor known is still refused, so a
    typo is not queued into a silent black hole."""
    perms = SimpleNamespace(presence_key="nacht-s16")
    reply = _call(perms, to="no-such-session-xyz", message="x")
    assert "error" in reply
