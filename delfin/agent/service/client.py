"""A client of the agent service, for Python: tests, `delfin-agent service`,
and the launcher's check that the service is up.

    with ServiceClient() as svc:
        svc.request("sessions.list")

The TUI (TypeScript) speaks the same protocol over the same socket.
"""

from __future__ import annotations

import socket
import threading
import time
from pathlib import Path
from typing import Any, Callable

from . import protocol, rundir


class ServiceError(Exception):
    def __init__(self, code: int, message: str):
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message


class ServiceClient:
    """One connection: requests and their answers, and the events sent."""

    def __init__(self, path: Path | None = None, *, timeout: float = 10.0,
                 client: str = "python", hello: bool = True):
        self.path = Path(path) if path is not None else rundir.socket_path()
        self.sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        self.sock.settimeout(timeout)
        self.sock.connect(str(self.path))
        self.sock.settimeout(None)
        self._file = self.sock.makefile("rb")
        self._lock = threading.Lock()
        self._cond = threading.Condition()
        self._responses: dict[Any, dict] = {}
        self._next_id = 0
        self.events: list[dict] = []
        self.closed = False
        self._reader = threading.Thread(target=self._read, daemon=True,
                                        name="service-client")
        self._reader.start()
        self.hello: dict = {}
        if hello:
            self.hello = self.request("hello", {"protocolVersion":
                                                protocol.PROTOCOL_VERSION,
                                                "client": client}, timeout=timeout)

    def _read(self) -> None:
        try:
            for message in protocol.read_messages(self._file):
                with self._cond:
                    if message.get("method") == "event":
                        self.events.append(message.get("params") or {})
                    else:
                        self._responses[message.get("id")] = message
                    self._cond.notify_all()
        except (OSError, ValueError):
            pass
        finally:
            with self._cond:
                self.closed = True
                self._cond.notify_all()

    def send(self, method: str, params: dict | None = None) -> int:
        with self._lock:
            self._next_id += 1
            msg_id = self._next_id
            self.sock.sendall(protocol.encode(protocol.request(msg_id, method,
                                                               params or {})))
        return msg_id

    def result(self, msg_id: int, timeout: float = 30.0) -> Any:
        deadline = time.monotonic() + timeout
        with self._cond:
            while msg_id not in self._responses:
                left = deadline - time.monotonic()
                if self.closed:
                    raise ServiceError(-1, "the service closed the connection")
                if left <= 0:
                    raise TimeoutError(f"no answer to request {msg_id}")
                self._cond.wait(left)
            message = self._responses.pop(msg_id)
        if "error" in message:
            err = message["error"]
            raise ServiceError(int(err.get("code", -1)), str(err.get("message")))
        return message.get("result")

    def request(self, method: str, params: dict | None = None,
                timeout: float = 30.0) -> Any:
        return self.result(self.send(method, params), timeout)

    def wait_for(self, match: str | Callable[[dict], bool], *, after: int = 0,
                 timeout: float = 30.0) -> tuple[dict, int]:
        pred = (lambda e: e.get("type") == match) if isinstance(match, str) else match
        deadline = time.monotonic() + timeout
        with self._cond:
            while True:
                for i in range(after, len(self.events)):
                    if pred(self.events[i]):
                        return self.events[i], i
                left = deadline - time.monotonic()
                if left <= 0 or self.closed:
                    seen = [e.get("type") for e in self.events[after:]][-12:]
                    raise TimeoutError(f"no {match!r} event; saw {seen}")
                self._cond.wait(min(left, 0.2))

    def close(self) -> None:
        try:
            self.sock.shutdown(socket.SHUT_RDWR)
        except OSError:
            pass
        try:
            self.sock.close()
        except OSError:
            pass
        self._reader.join(timeout=5)

    def __enter__(self) -> "ServiceClient":
        return self

    def __exit__(self, *exc) -> None:
        self.close()


def reachable(path: Path | None = None, timeout: float = 2.0) -> dict | None:
    """The service's hello when one answers on this host, else None."""
    return probe(path, timeout)[0]


def probe(path: Path | None = None, timeout: float = 2.0) -> tuple[dict | None, str]:
    """``(hello, "")`` when the service answers; ``(None, why)`` when not.

    ``why`` tells nothing listening (``"none"``) from a service that is
    there and does not say hello: one of another protocol version (an
    update with the old service still running), or one that does not
    answer at all.
    """
    try:
        with ServiceClient(path, timeout=timeout, client="probe") as svc:
            return dict(svc.hello), ""
    except (FileNotFoundError, ConnectionRefusedError):
        return None, "none"
    except ServiceError as exc:
        if exc.code == protocol.VERSION_MISMATCH:
            return None, f"it speaks another protocol version ({exc.message})"
        return None, f"it refused the hello ({exc.message})"
    except (TimeoutError, socket.timeout):
        return None, f"it did not answer within {timeout:.0f} s"
    except (OSError, ValueError) as exc:
        return None, f"it could not be reached ({exc})"
