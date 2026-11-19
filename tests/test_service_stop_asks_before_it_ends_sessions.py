"""`delfin-agent service stop|restart` asks before it ends the sessions.

Bug hunt 11.11.2026, #15: stop and restart ended every session of the
service -- running turns stopped at their next safe point -- without a
question, while `service end` and `stop-all` ask. Now: on a terminal they
ask y/N (``--yes`` skips it); without a terminal they refuse only while a
session works or waits for you, as idle sessions are saved and offered
again (what scripts and the service's own tests rely on).
"""

from __future__ import annotations

import argparse

import pytest

from delfin.agent import approval_guard
from delfin.agent.service import cli_service, client


def _row(sid, state, saves=True):
    return {"id": sid, "title": f"t-{sid}", "workspace": "/w", "state": state,
            "saves": saves}


@pytest.fixture()
def service(monkeypatch):
    """A service that answers, with the rows a test sets; records stops."""
    monkeypatch.delenv(approval_guard.ENV_MARKER, raising=False)
    state = {"rows": [], "stopped": [], "started": 0}
    monkeypatch.setattr(client, "probe", lambda: ({"pid": 1}, ""))

    class _Svc:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    monkeypatch.setattr(client, "ServiceClient", lambda **kw: _Svc())
    monkeypatch.setattr(cli_service, "_session_rows",
                        lambda svc: (list(state["rows"]), {}))
    monkeypatch.setattr(cli_service, "_stop",
                        lambda restart=False: state["stopped"].append(restart) or 0)

    def _again():
        state["started"] += 1
        return 0

    monkeypatch.setattr(cli_service, "_start_again", _again)
    return state


def _tty(monkeypatch, on: bool):
    monkeypatch.setattr(cli_service, "_interactive", lambda: on)


def _run(action, **kw):
    return cli_service.cmd_service(argparse.Namespace(action=action, **kw))


@pytest.mark.parametrize("action", ["stop", "restart"])
def test_a_terminal_is_asked_and_no_ends_nothing(service, monkeypatch, capsys, action):
    _tty(monkeypatch, True)
    service["rows"] = [_row("aaaa1111", "working"), _row("bbbb2222", "idle")]
    asked = []
    monkeypatch.setattr("builtins.input", lambda prompt: asked.append(prompt) or "n")
    assert _run(action) == 1
    assert service["stopped"] == [] and service["started"] == 0
    out = capsys.readouterr().out
    assert "[y/N]" in asked[0] and "2 sessions" in asked[0]
    assert "its turn is stopped at its next safe point" in out
    assert "nothing stopped" in out


@pytest.mark.parametrize("action", ["stop", "restart"])
def test_yes_on_the_terminal_stops(service, monkeypatch, action):
    _tty(monkeypatch, True)
    service["rows"] = [_row("aaaa1111", "idle")]
    monkeypatch.setattr("builtins.input", lambda prompt: "y")
    assert _run(action) == 0
    assert service["stopped"] == [action == "restart"]


@pytest.mark.parametrize("action", ["stop", "restart"])
def test_the_yes_flag_skips_the_question(service, monkeypatch, action):
    _tty(monkeypatch, True)
    service["rows"] = [_row("aaaa1111", "working")]
    monkeypatch.setattr("builtins.input",
                        lambda prompt: pytest.fail("asked despite --yes"))
    assert _run(action, yes=True) == 0
    assert service["stopped"] == [action == "restart"]


def test_no_terminal_with_a_working_session_refuses(service, monkeypatch, capsys):
    _tty(monkeypatch, False)
    service["rows"] = [_row("aaaa1111", "waiting for you"), _row("b", "idle")]
    assert _run("stop") == 2
    assert service["stopped"] == []
    err = capsys.readouterr().err
    assert "not asked" in err and "--yes" in err


def test_no_terminal_with_only_idle_sessions_stops_as_before(service, monkeypatch):
    _tty(monkeypatch, False)
    service["rows"] = [_row("aaaa1111", "idle"), _row("bbbb2222", "sleeping")]
    assert _run("restart") == 0
    assert service["stopped"] == [True] and service["started"] == 1


def test_nothing_to_end_is_not_asked(service, monkeypatch):
    _tty(monkeypatch, True)
    service["rows"] = [_row("aaaa1111", "sleeping")]
    monkeypatch.setattr("builtins.input",
                        lambda prompt: pytest.fail("asked with nothing to end"))
    assert _run("stop") == 0


def test_an_unanswering_service_asks_on_a_terminal(monkeypatch, capsys):
    from delfin.agent.service import rundir
    monkeypatch.delenv(approval_guard.ENV_MARKER, raising=False)
    monkeypatch.setattr(client, "probe", lambda: (None, "refused"))
    monkeypatch.setattr(rundir, "service_pid", lambda: 999_999)
    killed = []
    monkeypatch.setattr(cli_service.os, "kill", lambda pid, sig: killed.append(pid))
    _tty(monkeypatch, True)
    monkeypatch.setattr("builtins.input", lambda prompt: "")
    assert _run("stop") == 1
    assert killed == []
    assert "nothing stopped" in capsys.readouterr().out
