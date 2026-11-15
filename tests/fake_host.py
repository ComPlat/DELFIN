"""A dashboard kernel's host link, spoken by hand, for the service tests.

``hello``, then ``session.host``; from its answer on the connection is the
session's pipe: the service's requests arrive here and are answered by
``handlers`` (method -> callable(params) -> result), and ``emit`` sends
the session's events. Every request received is kept in ``asked``.
Independent of ``delfin.agent.service.host`` on purpose: what the service
does with a host is tested against the wire, not against the kernel's
own client.
"""

from __future__ import annotations

import json
import socket
import threading
import time
from pathlib import Path

from delfin.agent.service import protocol
from delfin.agent.service.client import ServiceError


class FakeHost:
    def __init__(self, socket_path: Path, sid: str, cwd: str, **params) -> None:
        self.sid = sid
        self.sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        self.sock.settimeout(30)
        self.sock.connect(str(socket_path))
        self.file = self.sock.makefile("rb")
        self.asked: list[dict] = []
        self.handlers: dict = {
            "session.status": lambda p: {"state": "idle", "approval": "default", "queued": 0},
            "session.prompt": lambda p: {"accepted": True},
            "session.cancel": lambda p: {"stopping": False, "refused": 0},
            "permission.reply": lambda p: {"applied": True, "message": ""},
            "session.setTitle": lambda p: {"title": str(p.get("title") or "")},
            "session.shutdown": lambda p: {"closing": True},
        }
        self.lock = threading.Lock()
        self.cond = threading.Condition()
        self.closed = False
        self.seq = int(params.pop("seq", 0) or 0)
        self._send({"jsonrpc": "2.0", "id": 1, "method": "hello",
                    "params": {"protocolVersion": protocol.PROTOCOL_VERSION,
                               "client": "delfin-dashboard-host"}})
        self.hello = self._answer(1)
        body = {"sessionId": sid, "cwd": cwd, **params}
        self._send({"jsonrpc": "2.0", "id": 2, "method": "session.host", "params": body})
        self.answer = self._answer(2)
        self.sock.settimeout(None)
        self.reader = threading.Thread(target=self._read, daemon=True, name="fake-host")
        self.reader.start()

    def _send(self, message: dict) -> None:
        with self.lock:
            self.sock.sendall(protocol.encode(message))

    def _answer(self, msg_id: int):
        while True:
            line = self.file.readline()
            if not line:
                raise ConnectionError("the service closed the host link")
            message = json.loads(line)
            if message.get("id") == msg_id and "method" not in message:
                if "error" in message:
                    err = message["error"]
                    raise ServiceError(int(err.get("code", -1)), str(err.get("message")))
                return message.get("result")

    def _read(self) -> None:
        try:
            for line in self.file:
                message = json.loads(line)
                method = message.get("method")
                if not method or "id" not in message:
                    continue
                params = message.get("params") or {}
                with self.cond:
                    self.asked.append({"method": method, **params})
                    self.cond.notify_all()
                handler = self.handlers.get(method)
                if handler is None:
                    self._send({"jsonrpc": "2.0", "id": message["id"],
                                "error": {"code": -32601, "message": f"no {method}"}})
                else:
                    self._send({"jsonrpc": "2.0", "id": message["id"],
                                "result": handler(params)})
        except (OSError, ValueError):
            pass
        finally:
            with self.cond:
                self.closed = True
                self.cond.notify_all()

    def emit(self, type_: str, **fields) -> None:
        event = {"type": type_, "session": self.sid, "ts": round(time.time(), 3), **fields}
        if type_ in protocol.DURABLE_EVENTS:
            self.seq += 1
            event["seq"] = self.seq
        self._send({"jsonrpc": "2.0", "method": "event", "params": event})

    def asked_for(self, method: str) -> list[dict]:
        with self.cond:
            return [r for r in self.asked if r["method"] == method]

    def wait_asked(self, method: str, timeout: float = 15.0) -> dict:
        deadline = time.monotonic() + timeout
        with self.cond:
            while True:
                for req in self.asked:
                    if req["method"] == method:
                        return req
                left = deadline - time.monotonic()
                if left <= 0 or self.closed:
                    raise TimeoutError(f"no {method} reached the host; asked {self.asked}")
                self.cond.wait(min(left, 0.2))

    def wait_closed(self, timeout: float = 15.0) -> bool:
        deadline = time.monotonic() + timeout
        with self.cond:
            while not self.closed:
                left = deadline - time.monotonic()
                if left <= 0:
                    return False
                self.cond.wait(min(left, 0.2))
        return True

    def close(self) -> None:
        try:
            self.sock.shutdown(socket.SHUT_RDWR)
        except OSError:
            pass
        self.sock.close()
