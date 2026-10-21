"""A session sees receipts only for mail it sent, and never a message body.

The session_message tool's status/ls and the `delfin-agent messages` CLI
answer from the asking session's own key. Another session's mail -- the
operator's included -- reads as unknown / nothing, and no row carries text.
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


def _tool(me, **arguments):
    perms = SimpleNamespace(presence_key=me)
    return json.loads(_DocToolExecutor()._execute_session_message(arguments, perms))


def test_another_sessions_mail_is_invisible_and_bodiless():
    sent = _tool("sess-a", to="operator", message="secret plan for the operator")
    assert sent["status"] == "sent"
    rows_b = _tool("sess-b", ls=True)["messages"]
    assert rows_b == []
    rows_a = _tool("sess-a", ls=True)["messages"]
    assert len(rows_a) == 1
    assert "secret plan" not in json.dumps(rows_a)
    mid = rows_a[0]["id"]
    assert _tool("sess-b", status=mid)["status"] == "unknown"
    assert _tool("sess-a", status=mid)["status"] in ("queued", "delivered", "read")


def test_a_session_cannot_queue_mail_to_itself():
    out = _tool("sess-a", to="sess-a", message="hello me")
    assert out.get("status") != "sent"
