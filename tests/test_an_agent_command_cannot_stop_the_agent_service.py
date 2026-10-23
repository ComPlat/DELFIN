"""`delfin-agent service stop|restart` from an agent session's own process
is refused, on both routes. Over the socket the service already answers
such a peer watch-only; when the socket does not answer (a caged command
on a Landlock + socket-guard host sees connect refused) the command used
to fall back to SIGTERM/SIGKILL by pid and so ended the service and every
session in it. The permission gate refuses the command before it runs,
under any program name, as it refuses `approvals approve`.
"""

from __future__ import annotations

import argparse

import pytest

from delfin.agent import approval_guard
from delfin.agent.service import cli_service, client, rundir


@pytest.mark.parametrize("action", ["stop", "restart"])
def test_the_pid_fallback_refuses_inside_an_agent_session(monkeypatch, capsys, action):
    killed: list = []
    monkeypatch.setenv(approval_guard.ENV_MARKER, "4242")
    monkeypatch.setattr(client, "probe", lambda: (None, "refused"))
    monkeypatch.setattr(rundir, "service_pid", lambda: 999_999)
    monkeypatch.setattr(cli_service.os, "kill", lambda pid, sig: killed.append((pid, sig)))
    monkeypatch.setattr(cli_service, "_start_again", lambda: killed.append("start") or 0)
    code = cli_service.cmd_service(argparse.Namespace(action=action))
    assert code != 0
    assert killed == []
    assert "agent session" in capsys.readouterr().err


def test_control_outside_an_agent_session_the_fallback_still_signals(monkeypatch):
    killed: list = []
    monkeypatch.delenv(approval_guard.ENV_MARKER, raising=False)
    monkeypatch.setattr(client, "probe", lambda: (None, "refused"))
    monkeypatch.setattr(rundir, "service_pid", lambda: 999_999)
    monkeypatch.setattr(cli_service.os, "kill", lambda pid, sig: killed.append((pid, sig)))
    monkeypatch.setattr(cli_service, "_gone", lambda pid, seconds: True)
    assert cli_service.cmd_service(argparse.Namespace(action="stop")) == 0
    assert killed and killed[0][0] == 999_999


@pytest.mark.parametrize("cmd", [
    "delfin-agent service stop",
    "delfin-agent service restart",
    "python -m delfin.agent.cli service stop",
    "del'fin'-agent \"service\" re\\start",
    "cd /tmp && delfin-agent service --quiet stop",
])
def test_the_gate_refuses_stopping_the_service(cmd):
    assert approval_guard.answer_command(cmd)


@pytest.mark.parametrize("cmd", [
    "delfin-agent service status",
    "delfin-agent service logs",
    "systemctl status sshd",
    "grep -n 'def _stop' delfin/agent/service/cli_service.py",
])
def test_reading_the_service_stays_allowed(cmd):
    assert approval_guard.answer_command(cmd) is None
