"""A dashboard's kernels talk over Unix sockets in the user's own directory.

Over TCP every kernel channel is a port on 127.0.0.1 that any user of the
machine can connect to, and the IOPub channel publishes to any subscriber:
the connection file's key signs messages, it does not encrypt them.
Measured 2026-10-20 (ipykernel 7.3, jupyter_client 8.10): a client without
the key read a kernel's printed output over TCP. ipykernel 7 prints
"Kernel is running over TCP without encryption" for it.

The test starts a real kernel through Jupyter's ServerApp with the
arguments the launcher passes and DELFIN's own kernel manager, and checks
where its sockets are.
"""

from __future__ import annotations

import asyncio
import os
import stat
import sys

import pytest

from delfin import cli_voila as V


@pytest.fixture
def short_runtime(monkeypatch):
    """A runtime directory short enough for a socket path on any machine
    (pytest's tmp_path alone can pass 107 bytes)."""
    import shutil
    import tempfile
    d = tempfile.mkdtemp(prefix="dk", dir="/tmp" if os.path.isdir("/tmp") else None)
    os.chmod(d, 0o700)
    monkeypatch.setenv("JUPYTER_RUNTIME_DIR", d)
    yield d
    shutil.rmtree(d, ignore_errors=True)


def test_the_server_is_told_to_use_ipc(monkeypatch, short_runtime):
    monkeypatch.delenv("DELFIN_KERNEL_TRANSPORT", raising=False)
    monkeypatch.setattr(V.os, "name", "posix")
    assert V._kernel_transport_args() == ["--KernelManager.transport=ipc"]


def test_a_runtime_too_deep_for_a_socket_keeps_tcp_and_says_so(
        monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(V.os, "name", "posix")
    monkeypatch.delenv("DELFIN_KERNEL_TRANSPORT", raising=False)
    deep = tmp_path / ("d" * 80)
    deep.mkdir()
    monkeypatch.setenv("JUPYTER_RUNTIME_DIR", str(deep))
    assert V._kernel_transport_args() == []
    assert "JUPYTER_RUNTIME_DIR" in capsys.readouterr().err


def test_tcp_can_be_chosen_and_windows_keeps_it(monkeypatch):
    monkeypatch.setattr(V.os, "name", "posix")
    monkeypatch.setenv("DELFIN_KERNEL_TRANSPORT", "tcp")
    assert V._kernel_transport_args() == []
    monkeypatch.delenv("DELFIN_KERNEL_TRANSPORT")
    monkeypatch.setattr(V.os, "name", "nt")
    assert V._kernel_transport_args() == []


@pytest.mark.skipif(sys.platform == "win32", reason="no ipc transport on Windows")
def test_a_kernel_of_the_server_listens_on_no_port(tmp_path, monkeypatch,
                                                   short_runtime):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.delenv("DELFIN_KERNEL_TRANSPORT", raising=False)
    from jupyter_server.serverapp import ServerApp

    app = ServerApp()
    app.initialize(argv=[
        *V._kernel_transport_args(),
        "--ServerApp.kernel_manager_class="
        "delfin.dashboard.resume_server.ResumeAwareMappingKernelManager",
        "--ServerApp.open_browser=False",
        f"--ServerApp.root_dir={tmp_path}",
    ], new_httpserver=False)
    km = app.kernel_manager

    async def start_and_read():
        kid = await km.start_kernel()
        try:
            info = km.get_kernel(kid).get_connection_info()
            client = km.get_kernel(kid).client()
            client.start_channels()
            await client.wait_for_ready(timeout=60)
            client.stop_channels()
            return info
        finally:
            await km.shutdown_kernel(kid, now=True)

    info = asyncio.run(start_and_read())
    assert info["transport"] == "ipc"
    sockets = os.path.dirname(info["ip"])
    assert os.path.realpath(sockets) == os.path.realpath(short_runtime)
    mode = stat.S_IMODE(os.stat(sockets).st_mode)
    assert mode & 0o077 == 0, oct(mode)
