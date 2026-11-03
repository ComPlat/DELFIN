"""Pending approvals of any session reach the paired messenger chat.

Input: the approval broker's request description
(delfin/agent/service/approvals.py), the bridge (delfin/agent/service/
bridge.py), the gateway's request texts (delfin/messenger/gateway.py)
and the schedule runner's session listener
(delfin/agent/service/schedule_runner.py).

Before, a scheduled (unattended) session's pending approval waited for
its timeout: the session is not /sync-activated, so the gateway's
refresh skipped it and the bridge refused answers.

- pending approvals of every session of this service -- interactive,
  scheduled, unattended -- reach the phone: one short message per
  request with its number, the tool, the command/path shortened, and
  the risk words the approval screen shows; several requests arriving
  together come in one batched message;
- answered with ``ok <n>`` / ``nein <n>`` (and ``zeige <n>`` for the
  full request); the answer carries the HMAC check the listing showed
  (trust_store), so only the paired device can answer; an answer from
  another chat, with a wrong or expired number is refused with a clear
  reply;
- a request the screen marks always-explicit (protected paths) is
  shown but approving it from the phone is refused: it stays explicit;
- the decision goes through the real broker (``permission.reply``),
  bound to the same live, unchanged request as every phone decision.
"""
from __future__ import annotations

import asyncio

import pytest

from delfin.agent.service.bridge import Bridge
from delfin.agent.service.redact import Redactor

CONFIRM = {
    "requestId": "p123-1", "kind": "confirm", "tool": "bash",
    "preview": "$ rm -r results", "purpose": "Räumt den Ordner auf",
    "flags": {"shell": True, "write": False, "protected": False,
              "outsideRead": False},
    "options": [{"key": "y", "label": "Once", "detail": ""},
                {"key": "n", "label": "No", "detail": ""}],
}
PROTECTED = {
    "requestId": "p123-2", "kind": "confirm", "tool": "write_file",
    "preview": "delfin/agent/sandbox.py: +2 lines",
    "flags": {"shell": False, "write": True, "protected": True,
              "outsideRead": False},
    "options": [{"key": "y", "label": "Once", "detail": ""},
                {"key": "n", "label": "No", "detail": ""}],
}


@pytest.fixture
def scheduled_bridge(tmp_path, monkeypatch):
    """A bridge whose service has a scheduled session with two pending
    approvals, and no messenger activation for it (the gap this fixes)."""
    from delfin.messenger import activation
    monkeypatch.setattr(activation, "_DIR", tmp_path / "activation")
    activation.configure("s1")           # an interactive session
    redactor = Redactor()
    redactor._known = ["SECRET-CREDENTIAL"]
    redactor._loaded = float("inf")
    bridge = Bridge(allowed_uid=1, redactor=redactor)
    requests = [CONFIRM, PROTECTED]
    calls = []

    async def request(method, params=None):
        calls.append((method, params))
        if method == "sessions.list":
            # s1 interactive (activated), s2 scheduled (NOT activated):
            return {"sessions": [{"sessionId": "s1", "title": "chat"},
                                 {"sessionId": "s2", "title": "morning run"}]}
        if method == "session.pending":
            return {"pending": list(requests) if params.get("sessionId") == "s2"
                    else []}
        if method == "permission.reply":
            requests.clear()
            return {"applied": True}
        return {"ok": True}

    bridge._request = request
    return bridge, requests, calls


def test_pending_approvals_of_a_scheduled_session_are_listed(scheduled_bridge):
    bridge, requests, _ = scheduled_bridge

    async def run():
        return await bridge._run_command({"cmd": "approvals"})

    reply = asyncio.run(run())
    items = reply["items"]
    assert [i["session"] for i in items] == ["morning run", "morning run"]
    # one short entry per request: number, tool, shortened action, risk word
    first = items[0]
    assert first["n"] >= 1
    assert first["tool"] == "bash"
    assert "rm -r results" in first["text"]
    assert first["approvable"] is True


def test_protected_request_is_shown_but_not_approvable(scheduled_bridge):
    bridge, _, _ = scheduled_bridge
    reply = asyncio.run(bridge._run_command({"cmd": "approvals"}))
    items = reply["items"]
    assert any("sandbox.py" in i["text"] for i in items)
    protected = next(i for i in items if "sandbox.py" in i["text"])
    assert protected["approvable"] is False


def test_ok_n_resolves_through_the_broker(scheduled_bridge):
    bridge, _, calls = scheduled_bridge

    async def run():
        listing = await bridge._run_command({"cmd": "approvals"})
        first = listing["items"][0]
        return await bridge._run_command(
            {"cmd": "answer", "n": first["n"], "ok": True,
             "check": listing["check"]})

    reply = asyncio.run(run())
    assert reply.get("ok"), reply
    assert any(m == "permission.reply" for m, _ in calls)


def test_wrong_or_missing_check_is_refused(scheduled_bridge):
    bridge, _, calls = scheduled_bridge

    async def run():
        listing = await bridge._run_command({"cmd": "approvals"})
        first = listing["items"][0]
        bad = await bridge._run_command(
            {"cmd": "answer", "n": first["n"], "ok": True, "check": "wrong"})
        missing = await bridge._run_command(
            {"cmd": "answer", "n": first["n"], "ok": True})
        stale = await bridge._run_command(
            {"cmd": "answer", "n": 9999, "ok": True, "check": listing["check"]})
        return bad, missing, stale

    bad, missing, stale = asyncio.run(run())
    for reply in (bad, missing, stale):
        assert "error" in reply
    assert not any(m == "permission.reply" for m, _ in calls)


def test_listing_numbers_are_stable_while_the_request_waits(scheduled_bridge):
    bridge, _, _ = scheduled_bridge

    async def run():
        first = await bridge._run_command({"cmd": "approvals"})
        second = await bridge._run_command({"cmd": "approvals"})
        return first, second

    a, b = asyncio.run(run())
    assert [i["n"] for i in a["items"]] == [i["n"] for i in b["items"]]


def test_an_expired_or_reused_number_cannot_approve(scheduled_bridge):
    bridge, _, calls = scheduled_bridge

    async def run():
        listing = await bridge._run_command({"cmd": "approvals"})
        n = listing["items"][0]["n"]
        await bridge._run_command(
            {"cmd": "answer", "n": n, "ok": True, "check": listing["check"]})
        # the request is gone now; the same number must not answer again
        return await bridge._run_command(
            {"cmd": "answer", "n": n, "ok": True, "check": listing["check"]})

    reply = asyncio.run(run())
    assert "error" in reply
    assert len([c for c in calls if c[0] == "permission.reply"]) == 1


def test_numbers_do_not_leak_credentials(scheduled_bridge):
    bridge, requests, _ = scheduled_bridge
    requests.append({**CONFIRM, "requestId": "p123-3",
                     "preview": "command SECRET-CREDENTIAL"})
    reply = asyncio.run(bridge._run_command({"cmd": "approvals"}))
    assert "SECRET-CREDENTIAL" not in str(reply)
