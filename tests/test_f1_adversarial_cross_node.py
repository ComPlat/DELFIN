"""Adversarial contract tests for Target 3 (cross-node reachability) —
reviewer s2.

The brief: a supervisor on a compute node is unreachable through a unix
socket from the login node; approvals ls/approve, messages, service
status and pause must work for a session whose host is a compute node
via the SHARED filesystem that already exists (session store, pause
flags, approval rooms). Never a plain TCP port without authentication.

On main 04bf9703b the building blocks exist (file_confirm.py records
carry host at :182; session_pause.py flag files at :29) but no
cross-node service channel — the module-existence case is red; the rest
arms for the builder hash.
"""

from __future__ import annotations

import json
import socket
import time
from pathlib import Path

import pytest


def _module():
    try:
        from delfin.agent.service import durable_host
        return durable_host
    except ImportError:
        return None


class TestCrossNodeReachability:
    def test_cross_node_channel_module_exists(self):
        assert _module() is not None, (
            "Target 3: no durable_host module yet — the cross-node "
            "channel cannot exist")

    def _pause_path(self, key: str) -> Path:
        from delfin.agent import session_pause
        return session_pause._path(key)

    def test_pause_flag_is_visible_from_another_process_namespace(
            self, tmp_path, monkeypatch):
        """The pause channel must be a plain shared file: pause in one
        process, is_paused in another (fresh import) reads it."""
        from delfin.agent import session_pause
        monkeypatch.setattr(session_pause, "_DIR", tmp_path / "pause")
        session_pause.pause("node-session", reason="steer while away")
        # A second "process": re-resolve through the same per-call path
        # (no caching), the way a login-node CLI would.
        assert session_pause.is_paused("node-session") is True
        assert session_pause.resume("node-session") is True
        assert session_pause.is_paused("node-session") is False

    def test_pause_path_cannot_escape_its_dir(self, tmp_path, monkeypatch):
        from delfin.agent import session_pause
        monkeypatch.setattr(session_pause, "_DIR", tmp_path / "pause")
        p = session_pause._path("../../escape")
        assert tmp_path / "pause" in p.parents or p.parent == tmp_path / "pause"
        assert ".." not in [x.name for x in p.parents]

    def test_approval_record_names_its_host(self, tmp_path, monkeypatch):
        """A cross-node question must carry the asking host, so the
        answering side knows where it came from (file_confirm._write_request
        already records host at file_confirm.py:182)."""
        from delfin.agent import file_confirm
        room = tmp_path / "approvals"
        monkeypatch.setattr(file_confirm, "requests_dir", lambda: room)
        broker = file_confirm.FileConfirmBroker(session_id="sid-x",
                                                require_code=False)
        request_id, path = broker._write_request("bash", {"command": "ls"},
                                                 "preview")
        record = json.loads(path.read_text(encoding="utf-8"))
        assert record["host"], "the request must name the host that asked"

    def test_no_tcp_listener_contract(self):
        """The channel must be filesystem-only: no module in the durable
        surface may bind a TCP socket. Verified against the module source
        once it exists — for now the absence contract is trivially true,
        so this pins the CHECK, to run against the hash."""
        mod = _module()
        if mod is None:
            pytest.skip("no durable_host module on this branch (red state)")
        source = Path(mod.__file__).read_text(encoding="utf-8")
        for banned in ("socket.bind", "socketserver.TCPServer",
                       "http.server", "socket.listen"):
            assert banned not in source, (
                f"durable_host must not open a TCP channel: {banned} found")
