"""A schedule is made where the user is: a session, the command line, and
-- proposed only -- by the agent.

* ``/schedule add <when> <task>`` in a session makes a confirmed schedule
  in the session's folder; ``/schedule`` lists it with its next run.
* ``delfin-agent schedule ...`` does the same from a shell (and refuses
  from an agent's shell: test_a_user_schedule_runs_only_as_confirmed).
* The agent asked for "every day at 7:30" calls ``cron_create`` with
  ``at``: a proposal, which never runs before the user's ``ok`` and is
  capped in number; ``/schedule ok N`` shows it and confirms on yes.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import schedules


@pytest.fixture
def ws(tmp_path, monkeypatch):
    from delfin.agent import schedule_commands
    monkeypatch.setattr(schedules, "_DIR", tmp_path / "user_schedules")
    # Never the user's own agent service.
    monkeypatch.setattr(schedule_commands, "_service_line", lambda: "")
    folder = tmp_path / "ws"
    folder.mkdir()
    return folder


def _dispatch(line, ws, ask=None):
    from delfin.agent import repl_commands as rc
    return rc.dispatch(line, rc.ReplContext(workspace=ws, ask=ask)).output


def test_a_session_makes_and_lists_a_schedule(ws):
    out = _dispatch("/schedule add mo-fr 18:00 Prüfe die Rechnungen", ws)
    assert "schedule #1 made" in out
    entry = schedules.get(1)
    assert entry["confirmed"] and entry["created_by"] == "session"
    assert entry["workspace"] == str(ws.resolve())
    assert entry["task"] == "Prüfe die Rechnungen"
    listing = _dispatch("/schedule", ws)
    assert "#1" in listing and "mo-fr 18:00" in listing and "active" in listing


def test_a_session_pauses_runs_and_deletes(ws):
    _dispatch("/schedule add alle 2h Schau nach", ws)
    assert "paused" in _dispatch("/schedule pause 1", ws)
    assert schedules.get(1)["paused"]
    assert "resumed" in _dispatch("/schedule resume 1", ws)
    assert "within a minute" in _dispatch("/schedule run 1", ws)
    assert schedules.get(1)["run_requested"]
    assert "deleted" in _dispatch("/schedule delete 1", ws)
    assert schedules.entries() == []


def test_a_wrong_time_says_how_to_write_it(ws):
    out = _dispatch("/schedule add irgendwann mal x", ws)
    assert out.startswith("schedule:") and "täglich 07:30" in out


def test_the_agent_only_proposes(ws):
    from delfin.agent import api_client
    reply = json.loads(api_client._propose_user_schedule(
        {"at": "täglich 07:30", "prompt": "Fasse die Jobs zusammen"}, str(ws)))
    assert reply["status"] == "proposed" and "ok" in reply["note"]
    entry = schedules.get(reply["id"])
    assert not entry["confirmed"] and entry["created_by"] == "agent"
    assert schedules.take_due(entry["host"], now=time.time() + 3 * 86400) == []
    for _ in range(4):
        api_client._propose_user_schedule({"at": "alle 1h", "prompt": "x"}, str(ws))
    refused = json.loads(api_client._propose_user_schedule({"at": "alle 1h", "prompt": "x"}, str(ws)))
    assert "error" in refused and "ok" in refused["error"]
    bad = json.loads(api_client._propose_user_schedule({"at": "sometime", "prompt": "x"}, str(ws)))
    assert "error" in bad


def test_the_user_confirms_a_proposal_after_seeing_it(ws):
    from delfin.agent import api_client
    reply = json.loads(api_client._propose_user_schedule(
        {"at": "täglich 07:30", "prompt": "Fasse die Jobs zusammen"}, str(ws)))
    n = reply["id"]
    shown = []
    out = _dispatch(f"/schedule ok {n}", ws, ask=lambda text: shown.append(text) or "n")
    assert "unconfirmed" in out and "Fasse die Jobs zusammen" in shown[0]
    assert not schedules.get(n)["confirmed"]
    out = _dispatch(f"/schedule ok {n}", ws, ask=lambda text: "y")
    assert "confirmed" in out and schedules.get(n)["confirmed"]


def test_a_surface_that_cannot_ask_confirms_with_yes(ws):
    from delfin.agent import api_client
    n = json.loads(api_client._propose_user_schedule(
        {"at": "alle 1h", "prompt": "x"}, str(ws)))["id"]
    out = _dispatch(f"/schedule ok {n}", ws)
    assert f"/schedule ok {n} yes" in out and not schedules.get(n)["confirmed"]
    _dispatch(f"/schedule ok {n} yes", ws)
    assert schedules.get(n)["confirmed"]


def test_the_cron_tool_offers_the_time_of_day():
    from delfin.agent import api_client
    tool = next(t for t in api_client._DOC_TOOLS_OPENAI
                if t["function"]["name"] == "cron_create")
    params = tool["function"]["parameters"]
    assert "at" in params["properties"] and params["required"] == ["prompt"]


def test_the_command_line_has_the_command(ws, capsys):
    from delfin.agent import cli
    args = cli.build_parser().parse_args(["schedule", "add", "alle", "1h", "x",
                                          "--cwd", str(ws)])
    code = args.func(args)
    assert code == 0, capsys.readouterr()
    assert schedules.get(1)["created_by"] == "cli"
