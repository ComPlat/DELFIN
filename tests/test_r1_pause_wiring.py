"""`delfin-agent pause <session>` reaches a running session.

Signals do not stop a supervised agent (SIGSTOP had no effect on 2026-10-03);
the pause flag does: the engine checks it at every tool boundary and takes
the ordinary stop path, the terminal does not wake a paused session, and
resume lifts it.
"""
from __future__ import annotations

import argparse

import pytest

from delfin.agent import cli, session_pause
from delfin.agent.engine import AgentEngine


@pytest.fixture(autouse=True)
def _own_home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setattr(session_pause, "_DIR", tmp_path / "pause", raising=False)


def _engine(key):
    eng = object.__new__(AgentEngine)
    eng.pause_key = key
    eng.stops = 0

    def _stop():
        eng.stops += 1
    eng.request_stop = _stop
    return eng


def test_a_paused_session_stops_at_the_next_tool_call():
    cli.cmd_pause(argparse.Namespace(session="nacht-s99", reason=""))
    eng = _engine("nacht-s99")
    assert eng._honour_pause() is True
    assert eng.stops == 1


def test_resume_lets_the_tools_run_again():
    cli.cmd_pause(argparse.Namespace(session="nacht-s99", reason=""))
    cli.cmd_resume(argparse.Namespace(session="nacht-s99"))
    eng = _engine("nacht-s99")
    assert eng._honour_pause() is False
    assert eng.stops == 0


def test_another_sessions_pause_does_not_stop_me():
    cli.cmd_pause(argparse.Namespace(session="nacht-s98", reason=""))
    eng = _engine("nacht-s99")
    assert eng._honour_pause() is False


def test_an_unnamed_engine_is_never_paused_by_name():
    cli.cmd_pause(argparse.Namespace(session="nacht-s99", reason=""))
    eng = _engine("")
    assert eng._honour_pause() is False
