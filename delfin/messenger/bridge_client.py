"""The gateway's side of the bridge: JSON lines over a Unix socket.

The bridge is the agent service's door (delfin.agent.service.bridge):
one JSON object per line both ways. A request is ``{"cmd": ...}``, its
answer the next line with either an ``error`` or the command's result;
event lines carry ``{"event": ...}`` and arrive between answers.
"""

from __future__ import annotations

import asyncio
import json
import os
import socket
from pathlib import Path

MAX_LINE_BYTES = 8_192


class BridgeDown(RuntimeError):
    """The bridge socket is not there or does not answer."""


class BridgeClient:
    def __init__(self, path: Path):
        self.path = Path(path)

    def connect(self) -> socket.socket:
        sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        sock.settimeout(15.0)
        try:
            sock.connect(str(self.path))
        except OSError as exc:
            sock.close()
            raise BridgeDown(f"{self.path}: {exc}") from None
        return sock

    def ask(self, command: dict) -> dict:
        """One command, the answer line. Events that arrive first are
        handed to ``on_event`` (called with the event dict), if given."""
        on_event = getattr(self, "on_event", None)
        sock = self.connect()
        try:
            line = (json.dumps(command, ensure_ascii=False) + "\n").encode()
            if len(line) > MAX_LINE_BYTES:
                raise ValueError("command too long")
            sock.sendall(line)
            with sock.makefile("rb") as stream:
                while True:
                    data = stream.readline(MAX_LINE_BYTES + 1)
                    if not data:
                        raise BridgeDown("the bridge closed the connection")
                    if len(data) > MAX_LINE_BYTES or not data.endswith(b"\n"):
                        raise BridgeDown("line too long")
                    message = json.loads(data)
                    if "event" in message:
                        if on_event is not None:
                            on_event(message)
                        continue
                    return message
        finally:
            sock.close()

    async def events(self):
        reader, writer = await asyncio.open_unix_connection(str(self.path), limit=MAX_LINE_BYTES)
        try:
            while line := await reader.readline():
                event = json.loads(line)
                if "event" in event:
                    yield event
        finally:
            writer.close()
            await writer.wait_closed()


def default_socket_path() -> Path:
    """Where the bridge serves on this host: the agent service's runtime
    directory, which the gateway may reach only as the configured group
    member (the provisioning chgrps it)."""
    explicit = os.environ.get("DELFIN_MESSENGER_SOCKET")
    if explicit:
        return Path(explicit)
    base = os.environ.get("DELFIN_RUNTIME_DIR", "").strip() or "/tmp"
    return Path(base) / f"delfin-{_service_uid()}" / "messenger.sock"


def _service_uid() -> int:
    """The uid the agent service runs under. The gateway runs as another
    user (delfin-gw), so it names the service's uid: DELFIN_MESSENGER_
    SERVICE_UID when set, else the same uid as here (same-user testing)."""
    raw = os.environ.get("DELFIN_MESSENGER_SERVICE_UID", "").strip()
    if raw:
        try:
            return int(raw)
        except ValueError:
            pass
    return os.getuid()
