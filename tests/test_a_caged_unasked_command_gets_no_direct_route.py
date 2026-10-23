"""The network of a command the gate admitted without a question.

Input: ``_sandbox_network_args(perms, allow_slurm)`` and the bubblewrap argv
``_bash_isolation_argv(..., caged_unasked=True)``, with the egress proxy,
the SLURM controller list and the direct TCP endpoints (SSH git remotes,
license servers) replaced by fakes, for every network mode and with and
without the socket guard.

Output asserted: a caged_unasked command gets no ``--allow-tcp`` at all
(a job runs outside the cage; an SSH remote or a license server reaches
another host's processes) and never ``--net open``. The control: the same
build for an approved command carries both endpoints, so the assertion
sees them when they are there.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent import egress_proxy as EP
from delfin.agent import socket_guard as SG
from delfin.agent.api_client import KitToolPermissions

SLURM = "10.0.0.5:6817"
DIRECT = "192.0.2.7:22"


class _Proxy:
    socket_path = "/nonexistent-delfin-proxy/proxy.sock"

    def token_file_for(self, perms):
        return "/nonexistent-delfin-proxy/token"


@pytest.fixture
def fakes(monkeypatch):
    monkeypatch.setattr(A, "_SLURM_CONTROLLERS", None)
    monkeypatch.setattr(EP, "slurm_controllers", lambda: [SLURM])
    monkeypatch.setattr(EP, "direct_tcp_endpoints", lambda ws=None: [DIRECT])
    monkeypatch.setattr(EP, "get_proxy", lambda: _Proxy())

    def _mode(m):
        monkeypatch.setattr(EP, "network_mode", lambda: m)
    return _mode


def _net(argv: list[str]) -> list[str]:
    """The value of every --net in argv."""
    return [argv[i + 1] for i, a in enumerate(argv[:-1]) if a == "--net"]


def _tcp(argv: list[str]) -> list[str]:
    return [argv[i + 1] for i, a in enumerate(argv[:-1]) if a == "--allow-tcp"]


@pytest.mark.parametrize("mode", ["proxy", "open", "none"])
def test_the_args_hold_no_direct_endpoint_and_never_open(mode, fakes, tmp_path):
    fakes(mode)
    perms = KitToolPermissions(mode="default", workspace=str(tmp_path))
    args, _ = A._sandbox_network_args(perms, allow_slurm=False)
    assert _tcp(args) == []
    assert _net(args) and "open" not in _net(args)


def test_control_an_approved_command_keeps_both_endpoints(fakes, tmp_path):
    fakes("proxy")
    perms = KitToolPermissions(mode="default", workspace=str(tmp_path))
    args, _ = A._sandbox_network_args(perms)
    assert set(_tcp(args)) == {SLURM, DIRECT}


@pytest.mark.parametrize("guard", [True, False])
@pytest.mark.parametrize("mode", ["proxy", "open", "none"])
def test_the_bwrap_argv_of_a_caged_unasked_command(mode, guard, fakes,
                                                   tmp_path, monkeypatch):
    fakes(mode)
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: ("bwrap", "auto"))
    monkeypatch.setattr(A, "_landlock_functional", lambda: True)
    monkeypatch.setattr(SG, "available", lambda: guard)
    ws = tmp_path / "ws"
    ws.mkdir()
    perms = KitToolPermissions(mode="default", workspace=str(ws))
    argv = A._bash_isolation_argv("make", ws, perms, caged_unasked=True)
    assert argv[0] == "bwrap", argv
    assert _tcp(argv) == []
    assert _net(argv) and "open" not in _net(argv)
    if not guard:
        # Without the guard a network namespace of its own holds it.
        assert "--unshare-net" in argv
        assert "--netns" in argv


@pytest.mark.parametrize("guard", [True, False])
def test_control_the_bwrap_argv_of_an_approved_command(guard, fakes, tmp_path,
                                                      monkeypatch):
    fakes("proxy")
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: ("bwrap", "auto"))
    monkeypatch.setattr(SG, "available", lambda: guard)
    ws = tmp_path / "ws"
    ws.mkdir()
    perms = KitToolPermissions(mode="default", workspace=str(ws))
    argv = A._bash_isolation_argv("make", ws, perms)
    assert set(_tcp(argv)) == {SLURM, DIRECT}
