"""The dashboard uses only a service of this release's capabilities.

A dashboard view relies on what SERVICE_VERSION 4 added (summary names,
the empty rename, permission.resolved.client). An older service still
running (an update, its sessions busy) answers hello with its own
serviceVersion: the dashboard goes through launcher.ensure_service (which
moves an idle older service over to this release) and, if what answers is
still older, falls back to its kernel and says why -- never a view that
silently lacks names and who answered.

Control: the real service of this release is connected to.
"""

from __future__ import annotations

import json
import os
import socket
import tempfile
import threading
from pathlib import Path

import pytest

from delfin.agent.service import protocol
from tui_service import ServiceProcess


class _OldService:
    """A socket that answers hello as a service of serviceVersion 3."""

    def __init__(self, directory: Path):
        self.path = directory / "agent.sock"
        self.sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        self.sock.bind(str(self.path))
        self.sock.listen(4)
        threading.Thread(target=self._serve, daemon=True).start()

    def _serve(self):
        while True:
            try:
                conn, _ = self.sock.accept()
            except OSError:
                return
            threading.Thread(target=self._talk, args=(conn,), daemon=True).start()

    def _talk(self, conn):
        f = conn.makefile("rb")
        for line in f:
            msg = json.loads(line)
            hello = {"protocolVersion": protocol.PROTOCOL_VERSION, "serviceVersion": 3,
                     "pid": 1, "host": "h", "mayAct": True, "sessions": []}
            conn.sendall(protocol.encode(protocol.result(msg.get("id"), hello)))

    def close(self):
        self.sock.close()


@pytest.fixture
def old_dir():
    d = Path(tempfile.mkdtemp(prefix="dov", dir="/tmp"))
    rt = d / f"delfin-{os.getuid()}"
    rt.mkdir(mode=0o700)
    yield d, rt
    import shutil
    shutil.rmtree(d, ignore_errors=True)


def test_an_older_service_is_not_used_and_the_reason_is_said(old_dir, monkeypatch, tmp_path):
    pytest.importorskip("ipywidgets")
    from delfin.agent.service import launcher
    from delfin.dashboard import live_sessions as LS

    base, rt = old_dir
    old = _OldService(rt)
    monkeypatch.setenv("DELFIN_RUNTIME_DIR", str(base))
    asked = []

    def ensure(**kw):
        # The move-over found the older one busy and left it running (its
        # hello is what ensure_service answers then).
        asked.append(kw)
        return {"serviceVersion": 3}
    monkeypatch.setattr(launcher, "ensure_service", ensure)
    try:
        with pytest.raises(LS.ServiceUnavailable) as no:
            LS.connect(start=True)
        assert "older" in str(no.value) and "3" in str(no.value)
        assert asked, "the launcher's move-over was not asked"
        with pytest.raises(LS.ServiceUnavailable):
            LS.connect(start=False)
    finally:
        old.close()
    # Control: this release's service is used.
    svc = ServiceProcess(tmp_path)
    try:
        client = LS.connect(start=False, path=svc.socket)
        assert client.hello["serviceVersion"] >= 4
        client.close()
    finally:
        svc.stop()
