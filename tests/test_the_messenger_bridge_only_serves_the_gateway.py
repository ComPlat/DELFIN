"""The messenger bridge: only the gateway's uid connects, only questions
are answered, only redacted and cut output leaves.

* a connection from another uid (or an unknown one) is closed without a
  word; the configured one is served;
* ``status``, ``sessions``, ``pending`` answer from the service, cut and
  redacted;
* a secret value in a session title or a question's text never reaches
  the line out;
* ``answer`` refuses a confirm (an approval) and an unknown id; a real
  question's answer reaches the service's attention.answer;
* an oversized line is refused; a line that is not JSON gets one error;
* ``stop`` shuts the service down; ``prompt`` without a named session
  reuses the one "phone" session, creating it once.
"""

from __future__ import annotations

import asyncio
import json
import os
import socket
import threading
from pathlib import Path

import pytest

from delfin.agent.service import bridge as bridge_mod
from delfin.agent.service.redact import Redactor


class FakeService:
    """Stands in for the agent service: answers what the bridge asks,
    records what it was asked."""

    def __init__(self):
        self.asked: list[tuple[str, dict]] = []
        self.events: list[dict] = []
        self.created: list[str] = []
        self.sessions: list[dict] = []
        self._cond = threading.Condition()
        self.closed = False
        self.hello = {"protocolVersion": 1}

    def request(self, method: str, params: dict | None = None,
                timeout: float = 30.0) -> dict:
        self.asked.append((method, params or {}))
        return self._answer(method, params or {})

    def _answer(self, method: str, params: dict) -> dict:
        if method == "service.status":
            return {"pid": 1, "host": "testhost",
                    "sessions": [{"id": "s1", "title": "calc 1",
                                  "state": "running", "pending": []}]}
        if method == "sessions.list":
            return {"sessions": [
                {"id": "s1", "title": "calc 1", "state": "running"},
                *self.sessions],
                "elsewhere": []}
        if method == "attention.list":
            return {"items": [
                {"id": "q1", "kind": "question", "sessionTitle": "calc 1",
                 "detail": "welches Lösemittel?"},
                {"id": "c1", "kind": "confirm", "sessionTitle": "calc 1",
                 "detail": "may run orca"}]}
        if method == "attention.answer":
            if params.get("id") == "c1":
                return {"ok": False, "error": "not a question"}
            return {"ok": True}
        if method == "session.create":
            self.created.append(str(params.get("title") or ""))
            self.sessions.append({"id": f"phone-{len(self.created)}",
                                  "title": str(params.get("title") or ""),
                                  "state": "running"})
            return {"sessionId": self.sessions[-1]["id"]}
        if method == "session.prompt":
            return {"accepted": True}
        if method == "service.shutdown":
            return {"stopping": True, "sessions": 1}
        raise AssertionError(f"unexpected method {method}")

    def close(self) -> None:
        self.closed = True


class _Rig:
    """One bridge serving one socket, with a fake service behind it."""

    def __init__(self, tmp_path: Path, allowed_uid: int):
        self.tmp = tmp_path
        self.service = FakeService()
        self.path = tmp_path / "messenger.sock"
        self.redactor = Redactor()
        self.redactor._known = ["TOPSECRET-KIT-KEY-1234"]
        self.redactor._loaded = float("inf")
        self.bridge = bridge_mod.Bridge(
            allowed_uid=allowed_uid, service_path=Path("/unused"),
            redactor=self.redactor)
        self.bridge._service = self.service
        self.server = None
        self.task = None

    async def start(self) -> None:
        self.server = await asyncio.start_unix_server(self.bridge._client_task,
                                                      path=str(self.path))
        os.chmod(self.path, 0o660)

    async def stop(self) -> None:
        if self.server is not None:
            self.server.close()
            await self.server.wait_closed()

    def conn(self) -> "BridgeConn":
        return BridgeConn(self.path)


class BridgeConn:
    """A client at the bridge socket."""

    def __init__(self, path: Path):
        self.sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        self.sock.settimeout(10.0)
        self.sock.connect(str(path))
        self.file = self.sock.makefile("rwb")

    def send(self, obj: dict, raw: bytes | None = None) -> None:
        data = raw if raw is not None else (json.dumps(obj) + "\n").encode()
        self.file.write(data)
        self.file.flush()

    def line(self) -> dict:
        data = self.file.readline()
        assert data, "connection closed without a reply"
        return json.loads(data)

    def ask(self, obj: dict) -> dict:
        self.send(obj)
        return self.line()

    def close(self) -> None:
        self.file.close()
        self.sock.close()


@pytest.fixture()
def rig(tmp_path):
    allowed = os.getuid()
    r = _Rig(tmp_path, allowed)
    loop = asyncio.new_event_loop()
    runner = threading.Thread(target=loop.run_forever, daemon=True)
    runner.start()
    fut = asyncio.run_coroutine_threadsafe(r.start(), loop)
    fut.result(5)
    yield r
    asyncio.run_coroutine_threadsafe(r.stop(), loop).result(5)
    loop.call_soon_threadsafe(loop.stop)


def test_a_connection_from_another_uid_is_closed(rig, tmp_path):
    # Directly the gate: an unknown peer uid (including None, where the
    # platform cannot say) refuses, the configured one passes.
    assert not bridge_mod.peer_uid_allowed(123, os.getuid())
    assert not bridge_mod.peer_uid_allowed(None, os.getuid())
    assert not bridge_mod.peer_uid_allowed(os.getuid(), None)
    assert bridge_mod.peer_uid_allowed(os.getuid(), os.getuid())

    # End to end: a bridge that allows nobody closes a peer silently --
    # no reply line, the connection just ends.
    refused = _Rig(tmp_path / "r2", allowed_uid=-1)
    refused.path = tmp_path / "refused.sock"

    async def run():
        server = await asyncio.start_unix_server(refused.bridge._client_task,
                                                 path=str(refused.path))
        async with server:
            try:
                reader, writer = await asyncio.open_unix_connection(str(refused.path))
            except (ConnectionRefusedError, ConnectionResetError):
                return "refused at connect"
            try:
                writer.write(b'{"cmd": "status"}\n')
                await writer.drain()
                await asyncio.wait_for(reader.readline(), timeout=2.0)
                return "answered"
            except (ConnectionResetError, asyncio.TimeoutError, asyncio.IncompleteReadError):
                return "closed"

    assert asyncio.run(run()) == "closed"
    # The refused peer never reached the service.
    assert not refused.service.asked


def test_status_and_sessions_are_cut_and_redacted(rig):
    real_answer = rig.service._answer
    rig.service._answer = lambda method, params: {
        "sessions": [{"id": "s1",
                      "title": "TOPSECRET-KIT-KEY-1234 in the title " + "x" * 900,
                      "state": "running", "pending": []}],
        "host": "h", "pid": 1}
    c = rig.conn()
    try:
        reply = c.ask({"cmd": "status"})
        assert reply.get("ok") is True
        title = reply["sessions"][0]["title"]
        assert "TOPSECRET" not in title
        assert "[redacted:credential]" in title
        assert len(title) <= 401
        # The request that reached the service was the service's own.
        assert rig.service.asked[0][0] == "service.status"
    finally:
        rig.service._answer = real_answer
        c.close()


def test_answer_refuses_an_approval_and_an_unknown_id(rig):
    c = rig.conn()
    try:
        reply = c.ask({"cmd": "answer", "id": "c1", "text": "yes"})
        assert reply.get("ok") is not True
        assert "approval" in reply["error"] or "never" in reply["error"]
        # The service's attention.answer was never asked for the confirm.
        assert not any(m == "attention.answer" and p.get("id") == "c1"
                       for m, p in rig.service.asked)
        reply = c.ask({"cmd": "answer", "id": "nope", "text": "hi"})
        assert reply.get("ok") is not True
        assert "no pending item" in reply["error"]
        # A question goes through.
        reply = c.ask({"cmd": "answer", "id": "q1", "text": "DMF"})
        assert reply.get("ok") is True
        assert ("attention.answer", {"id": "q1", "text": "DMF"}) in rig.service.asked
    finally:
        c.close()


def test_prompt_reuses_the_phone_session(rig):
    c = rig.conn()
    try:
        reply = c.ask({"cmd": "prompt", "text": "start the benzene run"})
        assert reply.get("ok") is True
        assert reply["session"] == "phone-1"
        first = [p for m, p in rig.service.asked if m == "session.create"]
        reply = c.ask({"cmd": "prompt", "text": "again"})
        second = [p for m, p in rig.service.asked if m == "session.create"]
        assert len(first) == len(second) == 1, "the phone session is created once"
        # The second prompt went to the same session.
        prompts = [p for m, p in rig.service.asked if m == "session.prompt"]
        assert all(p["sessionId"] == "phone-1" for p in prompts)
        # An unknown named session is refused.
        reply = c.ask({"cmd": "prompt", "session": "ghost", "text": "hi"})
        assert reply.get("ok") is not True
    finally:
        c.close()


def test_stop_shuts_the_service(rig):
    c = rig.conn()
    try:
        reply = c.ask({"cmd": "stop"})
        assert reply.get("ok") is True
        assert reply["stopping"] is True
        assert rig.service.asked[-1][0] == "service.shutdown"
    finally:
        c.close()


def test_oversized_and_broken_lines_get_one_error(rig):
    c = rig.conn()
    try:
        c.send(None, raw=b"x" * (bridge_mod.MAX_LINE_BYTES + 10) + b"\n")
        reply = c.line()
        assert "error" in reply
        c.send(None, raw=b"not json\n")
        assert "error" in c.line()
        c.send(None, raw=b'{"nonsense": 1}\n')
        assert "error" in c.line()
        # Still alive afterwards.
        assert c.ask({"cmd": "status"})["ok"] is True
    finally:
        c.close()


def test_unknown_command_lists_the_commands(rig):
    c = rig.conn()
    try:
        reply = c.ask({"cmd": "hacker"})
        assert reply.get("ok") is not True
        assert "commands" in reply
    finally:
        c.close()


def test_pending_names_approvals_without_granting_them(rig):
    c = rig.conn()
    try:
        reply = c.ask({"cmd": "pending"})
        items = {i["id"]: i for i in reply["items"]}
        assert items["c1"]["approvable"] is True
        assert items["q1"]["approvable"] is False
    finally:
        c.close()


def test_idle_gateway_receives_service_events(tmp_path, monkeypatch):
    """Exercise serve's real background pump, without a gateway command."""
    service = FakeService()
    bridge = bridge_mod.Bridge(allowed_uid=os.getuid())
    def connect():
        bridge._service = service
        return service
    monkeypatch.setattr(bridge, '_connect_service', connect)

    async def run():
        path = tmp_path / 'events.sock'
        task = asyncio.create_task(bridge.serve(path))
        try:
            for _ in range(100):
                if path.exists():
                    break
                await asyncio.sleep(.01)
            reader, writer = await asyncio.open_unix_connection(str(path))
            try:
                for _ in range(100):
                    if bridge._subscribers:
                        break
                    await asyncio.sleep(.01)
                with service._cond:
                    service.events.append({'type': 'session.summary', 'session': 's1',
                                           'summary': {'state': 'idle', 'pending': 0}})
                    service._cond.notify_all()
                event = json.loads(await asyncio.wait_for(reader.readline(), 2))
                assert event['event'] == 'session'
                assert event['state'] == 'idle'
            finally:
                writer.close()
                await writer.wait_closed()
        finally:
            task.cancel()
            await asyncio.gather(task, return_exceptions=True)
    asyncio.run(run())


def test_named_prompt_resolves_title_to_session_id(rig):
    c = rig.conn()
    try:
        assert c.ask({'cmd': 'prompt', 'session': 'calc 1', 'text': 'run'})['ok']
        assert rig.service.asked[-1][1]['sessionId'] == 's1'
    finally:
        c.close()


def test_status_accepts_real_service_summary(rig):
    rig.service._answer = lambda *args: {'sessions': [
        {'sessionId': 's1', 'state': 'idle', 'pending': 2}]}
    c = rig.conn()
    try:
        reply = c.ask({'cmd': 'status'})
        assert reply['sessions'][0]['pending'] == 2
        assert reply['sessions'][0]['id'] == 's1'
    finally:
        c.close()
