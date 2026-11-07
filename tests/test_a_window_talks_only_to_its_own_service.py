"""A window of the agent service talks only to the person's own service.

``ServiceClient`` (the dashboard's LiveClient, the bridge, the launcher's
probe and `delfin-agent service` all go through it) refuses, before
anything is sent:

* a socket whose directory is a symlink, not a directory, another uid's,
  or open to others -- a planted socket;
* a peer that runs as another uid (SO_PEERCRED).

The dashboard then falls back to its kernel and says why; nothing is
started beside a planted socket. Each refusal has its control: the same
listener reached the proper way is connected to.
"""

from __future__ import annotations

import os
import socket
import tempfile
import threading
from pathlib import Path

import pytest

from delfin.agent.service import client as C
from delfin.agent.service import rundir
from tui_service import ServiceProcess


class _Planted:
    """A listening Unix socket that counts who connects (and answers no hello)."""

    def __init__(self, directory: Path):
        self.path = directory / "agent.sock"
        self.sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        self.sock.bind(str(self.path))
        self.sock.listen(4)
        self.accepted = 0
        self._thread = threading.Thread(target=self._serve, daemon=True)
        self._thread.start()

    def _serve(self):
        while True:
            try:
                conn, _ = self.sock.accept()
            except OSError:
                return
            self.accepted += 1
            conn.close()

    def close(self):
        self.sock.close()


@pytest.fixture
def short_tmp():
    # A Unix socket path must fit in 108 bytes: pytest's tmp_path does not.
    d = Path(tempfile.mkdtemp(prefix="dwt", dir="/tmp"))
    yield d
    import shutil
    shutil.rmtree(d, ignore_errors=True)


def test_a_socket_behind_a_symlinked_directory_is_refused(short_tmp):
    real = short_tmp / "real"
    real.mkdir(mode=0o700)
    planted = _Planted(real)
    link = short_tmp / "rt"
    link.symlink_to(real)
    try:
        with pytest.raises(C.UnsafeSocket, match="not a plain directory"):
            C.ServiceClient(link / "agent.sock", timeout=1)
        assert planted.accepted == 0
        # Control: the same listener, through its own directory, is reached.
        with pytest.raises((TimeoutError, OSError, C.ServiceError)):
            C.ServiceClient(real / "agent.sock", timeout=1)
        assert planted.accepted == 1
    finally:
        planted.close()


def test_a_directory_others_may_enter_is_refused(short_tmp):
    d = short_tmp / "open"
    d.mkdir(mode=0o700)
    planted = _Planted(d)
    try:
        os.chmod(d, 0o755)
        with pytest.raises(C.UnsafeSocket, match="open to others"):
            C.ServiceClient(d / "agent.sock", timeout=1)
        assert planted.accepted == 0
    finally:
        planted.close()


def test_a_peer_of_another_uid_is_refused(tmp_path, monkeypatch):
    svc = ServiceProcess(tmp_path)
    try:
        # Control: the real service, as it is, says hello.
        with C.ServiceClient(svc.socket, client="test") as ok:
            assert ok.hello["protocolVersion"]
        # The kernel says another uid listens there (simulated: the test
        # cannot run a process as someone else).
        monkeypatch.setattr(rundir, "peer_uid", lambda _sock: os.getuid() + 1)
        with pytest.raises(C.UnsafeSocket, match="runs as uid"):
            C.ServiceClient(svc.socket, client="test")
        hello, why = C.probe(svc.socket)
        assert hello is None and why.startswith("unsafe:")
    finally:
        svc.stop()


def test_the_dashboard_refuses_a_planted_socket_and_runs_in_its_kernel(short_tmp, monkeypatch):
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import live_sessions as LS

    real = short_tmp / "real"
    real.mkdir(mode=0o700)
    planted = _Planted(real)
    runtime = short_tmp / "base"
    runtime.mkdir()
    (runtime / f"delfin-{os.getuid()}").symlink_to(real)
    monkeypatch.setenv("DELFIN_RUNTIME_DIR", str(runtime))
    started = []
    from delfin.agent.service import launcher
    monkeypatch.setattr(launcher, "start_service", lambda **kw: started.append(kw))
    try:
        with pytest.raises(LS.ServiceUnavailable) as no:
            LS.connect(start=True)
        assert "not safely yours" in str(no.value) or "not a plain directory" in str(no.value)
        assert planted.accepted == 0 and started == []
    finally:
        planted.close()
