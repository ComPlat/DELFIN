"""A session hosted in the agent service from another process.

A dashboard kernel runs its agent tab's sessions itself -- the tab's own
turn code, slash commands, uploads and dashboard mode -- and hosts each of
them in the person's agent service, so that every window lists, joins and
answers them (``supervisor.HostedSession``):

    link = HostLink(sid, cwd, handle=on_request, title=..., keep=...)
    link.emit("user.message", text="...", source="user")   # durable: logged
    link.emit_live("dashboard.chat", v=3, base=2, keep=120, tail="...")
    link.set_keep(False)
    link.release(end=False)                                 # x in its list

Input: the session id and folder, a ``handle(method, params) -> result``
for the requests the service forwards (session.prompt, session.cancel,
permission.reply, ...), and the session's events. Output: one connection
to the service; requests are answered on a thread each, so a slow one
(a slash command) never holds an Esc behind it.

Events are masked here as a worker masks its own (the worker's
``Redactor``); the service numbers the durable ones and writes them to
the session's event log, so a window that joins later replays them and a
worker that takes the session over later continues the same log. This
process never writes into the service's files.

Raises :class:`HostRefused` when the service will not host the session
(it runs there already; this process may not act), ``OSError`` /
``ServiceError`` when no service can be reached.
"""

from __future__ import annotations

import json
import socket
import threading
import time
from pathlib import Path
from typing import Any, Callable

from . import protocol, rundir
from .client import ServiceError, check_peer, check_socket_dir
from .redact import Redactor

#: What a host link says it is in its hello.
CLIENT_NAME = "delfin-dashboard-host"

#: Fields that identify things; never rewritten by the redactor (the
#: worker's own list).
_IDENTIFIERS = frozenset({"type", "session", "seq", "ts", "turnId", "partId",
                          "callId", "requestId", "sessionId", "source",
                          "kind", "state", "decision", "by", "flags"})

#: Live-only events not masked: the model's stream (shown as it arrives
#: in the hosting tab too) and numbers.
_UNREDACTED = frozenset({"text.delta", "reasoning.delta", "usage", "session.state"})


class HostRefused(Exception):
    """The service will not host the session; the text says why."""

    def __init__(self, code: int, message: str):
        super().__init__(message)
        self.code = code
        self.message = message


class HostLink:
    def __init__(self, session_id: str, cwd: str, *,
                 handle: Callable[[str, dict], Any],
                 path: Path | None = None, title: str = "", keep: bool = True,
                 info: dict | None = None, fake: bool = False,
                 on_close: Callable[[], None] | None = None,
                 timeout: float = 15.0) -> None:
        self.session_id = str(session_id)
        self.handle = handle
        self.on_close = on_close
        self.closed = False
        self._lock = threading.Lock()
        self._released = False
        self.redactor = Redactor()
        self.path = Path(path) if path is not None else rundir.socket_path()
        check_socket_dir(self.path)
        self.sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        self.sock.settimeout(timeout)
        try:
            self.sock.connect(str(self.path))
            check_peer(self.sock, self.path)
            self._file = self.sock.makefile("rb")
            self._send(protocol.request(1, "hello", {
                "protocolVersion": protocol.PROTOCOL_VERSION, "client": CLIENT_NAME}))
            self.hello = self._answer(1)
            params: dict[str, Any] = {"sessionId": self.session_id, "cwd": str(cwd),
                                      "keep": bool(keep), "info": dict(info or {})}
            if title:
                params["title"] = str(title)
            if fake:
                params["fake"] = True
            self._send(protocol.request(2, "session.host", params))
            self.answer = self._answer(2)
        except BaseException:
            self._close_socket()
            raise
        #: Nothing durable in the session's log yet (the service keeps it):
        #: the session never ran there -- a resume from before, a new one
        #: -- so its history goes in once.
        self.log_was_empty = not self.answer.get("history")
        self.sock.settimeout(None)
        self._reader = threading.Thread(target=self._read, daemon=True,
                                        name=f"host-{self.session_id[:8]}")
        self._reader.start()

    # -- the wire ---------------------------------------------------------------------
    def _send(self, message: dict) -> None:
        data = protocol.encode(message)
        with self._lock:
            self.sock.sendall(data)

    def _answer(self, msg_id: int) -> dict:
        while True:
            line = self._file.readline(protocol.MAX_LINE_BYTES + 1)
            if not line:
                raise ConnectionError("the agent service closed the connection")
            try:
                message = json.loads(line)
            except ValueError:
                continue
            if not isinstance(message, dict) or "method" in message:
                continue
            if message.get("id") != msg_id:
                continue
            if "error" in message:
                err = message.get("error") or {}
                code = int(err.get("code", protocol.INTERNAL_ERROR))
                text = str(err.get("message") or "refused")
                if msg_id == 2:
                    raise HostRefused(code, text)
                raise ServiceError(code, text)
            return message.get("result") or {}

    def _read(self) -> None:
        try:
            for message in protocol.read_messages(self._file):
                method = message.get("method")
                if not method or "id" not in message:
                    continue
                threading.Thread(target=self._serve, args=(message,), daemon=True,
                                 name=f"host-{self.session_id[:8]}-{method}").start()
        except (OSError, ValueError):
            pass
        finally:
            self._gone()

    def _serve(self, message: dict) -> None:
        msg_id = message.get("id")
        method = str(message.get("method") or "")
        params = message.get("params") or {}
        try:
            result = self.handle(method, params if isinstance(params, dict) else {})
            reply = protocol.result(msg_id, result)
        except HostRefused as exc:
            reply = protocol.error(msg_id, exc.code, exc.message)
        except Exception as exc:                  # noqa: BLE001 -- said to the window
            reply = protocol.error(msg_id, protocol.INTERNAL_ERROR,
                                   f"{type(exc).__name__}: {exc}")
        try:
            self._send(reply)
        except OSError:
            pass

    def _gone(self) -> None:
        with self._lock:
            if self.closed:
                return
            self.closed = True
        self._close_socket()
        if self.on_close is not None:
            try:
                self.on_close()
            except Exception:
                pass

    def _close_socket(self) -> None:
        try:
            self.sock.shutdown(socket.SHUT_RDWR)
        except OSError:
            pass
        try:
            self.sock.close()
        except OSError:
            pass

    # -- what the session says -----------------------------------------------------------
    def emit(self, type_: str, **fields: Any) -> bool:
        """One event of the session: stamped, masked, logged when durable,
        sent. False when the link is gone. Never raises."""
        if self.closed:
            return False
        try:
            event: dict[str, Any] = {"type": type_, "session": self.session_id,
                                     "ts": round(time.time(), 3)}
            event.update({k: v for k, v in fields.items() if v is not None})
            if type_ not in _UNREDACTED:
                event = self.redactor.event(event, _IDENTIFIERS)
            self._send(protocol.event(event))
            return True
        except Exception:                            # noqa: BLE001
            return False

    def emit_live(self, type_: str, **fields: Any) -> bool:
        """An event nothing keeps and nothing masks (the chat as the
        hosting tab shows it, to the person's own dashboards)."""
        if self.closed:
            return False
        try:
            event = {"type": type_, "session": self.session_id,
                     "ts": round(time.time(), 3), **fields}
            self._send(protocol.event(event))
            return True
        except Exception:                            # noqa: BLE001
            return False

    def set_keep(self, keep: bool) -> bool:
        return self.emit_live("host.keep", keep=bool(keep))

    def release(self, *, end: bool = False) -> None:
        """Let go of the session on purpose (the tab saved it and released
        its lock): it runs on in the service (Keep session on, or another
        window holds it), or ends -- ``end``: for every window."""
        if self._released:
            return
        self._released = True
        self.emit_live("host.released", end=bool(end))
        self.close()

    def close(self) -> None:
        self._close_socket()
        if threading.current_thread() is not self._reader:
            self._reader.join(timeout=5)
        self._gone()


def hosting_available(hello: dict) -> bool:
    """Whether a service that said ``hello`` hosts sessions (version 8)."""
    try:
        return int((hello or {}).get("serviceVersion") or 0) >= 8
    except (TypeError, ValueError):
        return False
