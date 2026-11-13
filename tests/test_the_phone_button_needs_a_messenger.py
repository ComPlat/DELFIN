"""The Phone button and /pair only where a phone messenger can be reached.

Input: the host's messenger -- a bridge socket (the root install's
``/run/delfin-messenger-bridge/messenger.sock``, DELFIN_MESSENGER_SOCKET,
the per-user default, or the rootless setup's
``$XDG_DATA_HOME/delfin-messenger/run/bridge.sock``) and the traces of an
install that is not running now (the root install's bridge unit, the
rootless setup's state directory or user units).

Output: ``activation.messenger_here()`` -- whether the agent tab shows the
Phone button -- and ``activation.pairing_unavailable()`` -- why /pair cannot
work here, or "". On a host with no messenger (an HPC login node) /pair
answered "Asking the messenger gateway ..." and, 30 s later, "No answer from
the messenger gateway"; it now says at once why, and writes no request.

Rejected: probing the bridge with a connection. The gateway's bridge admits
only the gateway's process; a dashboard connecting to it would be refused
and logged as an intruder.
"""
from __future__ import annotations

import shutil
import socket
import tempfile
from pathlib import Path

import pytest

from delfin.messenger import activation, pairing


@pytest.fixture
def host(tmp_path, monkeypatch):
    """A host with no messenger at all: every place it could be is empty."""
    short = Path(tempfile.mkdtemp(prefix="dphone-", dir="/tmp"))
    monkeypatch.delenv("DELFIN_MESSENGER_SOCKET", raising=False)
    monkeypatch.delenv("DELFIN_MESSENGER_SERVICE_UID", raising=False)
    monkeypatch.setenv("DELFIN_RUNTIME_DIR", str(short / "run"))
    monkeypatch.setenv("XDG_DATA_HOME", str(short / "data"))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    monkeypatch.setattr(activation, "PROVISIONED_SOCKET", str(tmp_path / "run-bridge" / "messenger.sock"))
    monkeypatch.setattr(activation, "PROVISIONED_UNIT", str(tmp_path / "etc" / "delfin-messenger-bridge.service"))
    monkeypatch.setattr(pairing, "_DIR", tmp_path / "pairing")
    yield short
    shutil.rmtree(short, ignore_errors=True)


def _listen(path: Path) -> socket.socket:
    path.parent.mkdir(parents=True, exist_ok=True)
    sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    sock.bind(str(path))
    return sock


def test_no_messenger_no_button_and_pair_says_why(host, tmp_path):
    assert activation.messenger_running() is False
    assert activation.messenger_here() is False
    why = activation.pairing_unavailable()
    assert why == activation.UNAVAILABLE
    assert "not available on this machine" in why
    out = pairing.session_command("", requester="session:s1", show=lambda _t: None,
                                  spawn=lambda fn: fn(), timeout=0)
    assert out == activation.UNAVAILABLE
    assert not (tmp_path / "pairing").exists()          # no request was written
    said = []
    assert pairing.cli_main(["pair"], out=said.append, timeout=0) == 2
    assert said == [activation.UNAVAILABLE]
    assert not (tmp_path / "pairing").exists()


def test_the_rootless_bridge_counts_as_running(host):
    sock = _listen(host / "data" / "delfin-messenger" / "run" / "bridge.sock")
    try:
        assert activation.messenger_running() is True
        assert activation.messenger_here() is True
        assert activation.pairing_unavailable() == ""
    finally:
        sock.close()


def test_the_root_install_bridge_counts_as_running(host, monkeypatch):
    path = host / "rb" / "messenger.sock"
    monkeypatch.setattr(activation, "PROVISIONED_SOCKET", str(path))
    sock = _listen(path)
    try:
        assert activation.messenger_here() is True
        assert activation.pairing_unavailable() == ""
    finally:
        sock.close()


@pytest.mark.parametrize("trace", ["unit", "rootless-state", "user-unit"])
def test_an_install_that_is_stopped_keeps_the_button_and_says_so(host, tmp_path, trace):
    if trace == "unit":
        unit = Path(activation.PROVISIONED_UNIT)
        unit.parent.mkdir(parents=True)
        unit.write_text("[Unit]\n")
    elif trace == "rootless-state":
        (host / "data" / "delfin-messenger" / "gateway").mkdir(parents=True)
    else:
        units = tmp_path / "config" / "systemd" / "user"
        units.mkdir(parents=True)
        (units / "delfin-messenger-bridge.service").write_text("[Unit]\n")
    assert activation.messenger_running() is False
    assert activation.messenger_here() is True
    assert activation.pairing_unavailable() == activation.STOPPED
    out = pairing.session_command("", requester="session:s2", show=lambda _t: None,
                                  spawn=lambda fn: fn(), timeout=0)
    assert out == activation.STOPPED
    assert not (tmp_path / "pairing").exists()


def test_cancel_and_usage_need_no_messenger(host):
    assert pairing.session_command("cancel", requester="session:s3", show=lambda _t: None) \
        == "No pairing open in this session."
    assert pairing.session_command("bogus", requester="session:s3", show=lambda _t: None) \
        .startswith("/pair ")


def test_the_agent_tab_adds_the_phone_button_only_where_a_messenger_is():
    source = (Path(__file__).resolve().parents[1] / "delfin" / "dashboard"
              / "tab_agent.py").read_text(encoding="utf-8")
    add = "git_group.children = tuple(git_group.children) + (messenger_phone_btn,)"
    assert source.count(add) == 1
    head = source[:source.index(add)].rstrip().splitlines()[-1].strip()
    assert head == "if _messenger_activation.messenger_here():"
