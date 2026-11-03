"""An old phone answer never frees a new action; the first answer wins.

Input: the messenger bridge (delfin/agent/service/bridge.py) against a real
agent service (tests/tui_service.ServiceProcess), a session whose turn asks
for a shell command, and the window answering locally.

Two races a phone answer can lose and must lose cleanly:
- the window answers first: the phone's code, issued for that request, is
  refused afterwards, and the decision that resolved the request is the
  window's;
- the same command is asked again in the next turn (an identical preview):
  the code issued for the first request does not approve the second, which
  stays open for the window.
"""
from __future__ import annotations

import asyncio

from delfin.agent.service.bridge import Bridge


def _asked(client, sid, after=None):
    event, _ = client.wait_for(lambda e: e.get("type") == "permission.asked"
                               and e.get("session") == sid
                               and e.get("requestId") != after)
    return event


def _setup(tmp_path, service):
    client = service.client()
    sid = service.session(client, mode="default")
    from delfin.messenger.activation import configure
    configure(sid)
    return client, sid, Bridge(allowed_uid=1, service_path=service.socket)


def _code(bridge, sid):
    async def run():
        reply = await bridge._run_command(
            {"cmd": "session", "session": sid, "action": "pending"})
        assert reply.get("items"), reply
        return reply["items"][0]["id"]
    return asyncio.run(run())


def _phone(bridge, sid, action, code):
    return asyncio.run(bridge._run_command(
        {"cmd": "session", "session": sid, "action": action, "id": code}))


def test_a_local_answer_first_leaves_the_phone_code_dead(tmp_path):
    from tui_service import ServiceProcess
    with ServiceProcess(tmp_path) as service:
        client, sid, bridge = _setup(tmp_path, service)
        try:
            client.request("session.prompt", {"sessionId": sid, "text": "bash"})
            asked = _asked(client, sid)
            code = _code(bridge, sid)
            local = client.request("permission.reply", {
                "sessionId": sid, "requestId": asked["requestId"],
                "choice": "n"})
            assert local.get("applied"), local
            event, _ = client.wait_for(
                lambda e: e.get("type") == "permission.resolved"
                and e.get("session") == sid)
            assert event["decision"] == "denied"
            late = _phone(bridge, sid, "approve", code)
            assert late.get("error") and not late.get("ok"), late
        finally:
            if bridge._service is not None:
                bridge._service.close()


def test_a_code_for_one_request_does_not_approve_its_twin(tmp_path):
    from tui_service import ServiceProcess
    with ServiceProcess(tmp_path) as service:
        client, sid, bridge = _setup(tmp_path, service)
        try:
            client.request("session.prompt", {"sessionId": sid, "text": "bash"})
            first = _asked(client, sid)
            old_code = _code(bridge, sid)
            client.request("permission.reply", {
                "sessionId": sid, "requestId": first["requestId"],
                "choice": "n"})
            client.wait_for(lambda e: e.get("type") == "turn.ended"
                            and e.get("session") == sid)
            client.request("session.prompt", {"sessionId": sid, "text": "bash"})
            second = _asked(client, sid, after=first["requestId"])
            assert second["requestId"] != first["requestId"]
            stale = _phone(bridge, sid, "approve", old_code)
            assert stale.get("error") and not stale.get("ok"), stale
            pending = client.request("session.pending", {"sessionId": sid})
            assert [p["requestId"] for p in pending.get("pending", [])] == [
                second["requestId"]], pending
        finally:
            if bridge._service is not None:
                bridge._service.close()
