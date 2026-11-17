"""`delfin-agent service end` says the state `service status` says.

With no service at all it said "it is not listening: it is starting or
stopping" -- a reason for a service that runs and does not answer. An
error names the real state (gh, kubectl).
"""
from __future__ import annotations

import pytest

from delfin.agent.service import cli_service, client, rundir


@pytest.fixture
def no_answer(monkeypatch):
    monkeypatch.setattr(cli_service, "_refused_inside_an_agent", lambda action: False)
    monkeypatch.setattr(client, "probe", lambda *a, **k: (None, "none"))
    monkeypatch.setattr(rundir, "hostname", lambda: "node-a")


def test_no_service_running_is_said_as_such(no_answer, monkeypatch, capsys):
    monkeypatch.setattr(rundir, "service_pid", lambda: None)
    assert cli_service._end("zzzz", all_idle=False, yes=True) == 1
    err = capsys.readouterr().err
    assert "no agent service runs on node-a" in err
    assert "starting or stopping" not in err


def test_a_service_that_does_not_answer_is_named_with_its_pid(no_answer, monkeypatch, capsys):
    monkeypatch.setattr(rundir, "service_pid", lambda: 4242)
    assert cli_service._end("zzzz", all_idle=False, yes=True) == 1
    err = capsys.readouterr().err
    assert "(pid 4242) does not answer" in err
