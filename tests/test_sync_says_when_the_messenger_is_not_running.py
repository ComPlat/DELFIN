"""/sync tells the truth about the messenger on this machine.

Measured on the Pi on 2026-11-12: the messenger bridge crash-looped (its
socket could not be created) and /sync still answered "SimpleX is on for this
session. The configured messenger invites you to its chat."; no invitation
came. The answer now says when no bridge serves on this machine; the switch
itself is still recorded, so the session reaches the phone once it runs.
"""
from __future__ import annotations

import socket

from delfin.messenger import activation


def _no_bridge(monkeypatch, tmp_path):
    monkeypatch.setenv("DELFIN_MESSENGER_SOCKET", str(tmp_path / "none.sock"))
    monkeypatch.setattr(activation, "PROVISIONED_SOCKET", str(tmp_path / "gone.sock"))
    monkeypatch.setenv("DELFIN_RUNTIME_DIR", str(tmp_path / "rt"))


def test_switching_on_without_a_bridge_says_so_and_records(monkeypatch, tmp_path):
    _no_bridge(monkeypatch, tmp_path)
    out = activation.configure("s1", "on", directory=tmp_path / "rec")
    assert "not running on this machine" in out and "invites" not in out
    assert activation.enabled("s1", directory=tmp_path / "rec")
    status = activation.configure("s1", "status", directory=tmp_path / "rec")
    assert status.startswith("SimpleX for this session: on.") and "not running" in status


def test_with_a_bridge_socket_the_invitation_is_announced(monkeypatch, tmp_path):
    _no_bridge(monkeypatch, tmp_path)
    path = tmp_path / "b.sock"
    srv = socket.socket(socket.AF_UNIX)
    srv.bind(str(path))
    try:
        monkeypatch.setenv("DELFIN_MESSENGER_SOCKET", str(path))
        out = activation.configure("s1", "on", directory=tmp_path / "rec")
        assert "invites you" in out and "not running" not in out
    finally:
        srv.close()


def test_a_switch_from_the_phone_never_warns(monkeypatch, tmp_path):
    _no_bridge(monkeypatch, tmp_path)
    out = activation.configure("s1", "on", directory=tmp_path / "rec", by="phone", phone="iPhone")
    assert "not running" not in out
