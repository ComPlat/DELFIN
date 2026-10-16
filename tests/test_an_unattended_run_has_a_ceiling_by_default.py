"""An unattended run stops at a ceiling nobody had to remember to set.

Input: an engine built for a run nobody watches (``--unattended``, a
headless run with no terminal, a scheduler entry), the settings
(``agent.run_budget_*``, ``agent.max_tool_rounds``,
``agent.unattended_budget``) and the flags. Output: the ceilings on the
engine and its client, the notice naming them, and the turn's stop when
one is reached in the middle of a turn.

Semantics:

* unattended, every run ceiling the flags and settings leave at 0 gets
  its default (``run_limits.UNATTENDED_DEFAULTS``: USD, wall clock,
  turns), and an uncapped tool-round setting gets a round cap; a flag or
  a setting always wins, and 0 in ``agent.unattended_budget`` lifts one;
* attended sessions are unchanged;
* the wall-clock and USD ceilings are checked between tool rounds, not
  only before a turn: a one-turn run stops in that turn, with the
  setting named as its stop reason, and the session is saved by the
  ordinary end of the turn.
"""
from __future__ import annotations

import contextlib
import io
import time
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import api_client as A
from delfin.agent import cli as agent_cli
from delfin.agent import engine as E
from delfin.agent import run_limits as RL


@contextlib.contextmanager
def _settings(agent):
    with patch("delfin.user_settings.load_settings",
               lambda *a, **k: {"agent": dict(agent)}):
        yield


class _Client:
    model = "kit.qwen3.5-397b-A17b"


def _engine():
    eng = E.AgentEngine.__new__(E.AgentEngine)     # no client, no network
    eng.client = _Client()
    return eng


def test_an_unattended_run_gets_every_ceiling_it_was_not_given():
    eng = _engine()
    with _settings({}):
        notes = RL.apply_unattended(eng)
        usd, secs = eng._run_budget()
        turns = eng._run_budget_turns()
    d = RL.UNATTENDED_DEFAULTS
    assert (usd, secs, turns) == (d["usd"], d["seconds"], d["turns"])
    assert d["usd"] > 0 and d["seconds"] > 0 and d["turns"] > 0
    text = " ".join(notes)
    assert "agent.unattended_budget" in text and "unattended" in text


def test_a_flag_and_a_setting_win_over_the_default():
    eng = _engine()
    eng.run_budget_s = 60.0                       # --max-run-seconds 60
    with _settings({"run_budget_usd": 2.0}):
        RL.apply_unattended(eng)
        usd, secs = eng._run_budget()
    assert secs == 60.0 and usd == 2.0


def test_zero_in_the_unattended_budget_lifts_one_ceiling():
    eng = _engine()
    with _settings({"unattended_budget": {"seconds": 0, "usd": 3}}):
        RL.apply_unattended(eng)
        usd, secs = eng._run_budget()
    assert secs == 0.0 and usd == 3.0


def test_an_uncapped_round_setting_gets_the_unattended_cap():
    eng = _engine()
    with _settings({"max_tool_rounds": 0}):
        RL.apply_unattended(eng)
    assert eng.client.max_tool_rounds == RL.UNATTENDED_DEFAULTS["tool_rounds"]
    eng2 = _engine()
    with _settings({"max_tool_rounds": 1200}):
        RL.apply_unattended(eng2)
    assert not hasattr(eng2.client, "max_tool_rounds")
    eng3 = _engine()
    eng3.client.max_tool_rounds = 7               # --max-turns 7
    with _settings({"max_tool_rounds": 0}):
        RL.apply_unattended(eng3)
    assert eng3.client.max_tool_rounds == 7


def test_the_wall_clock_is_checked_inside_a_turn():
    eng = _engine()
    eng.run_budget_s = 10.0
    eng._run_started_at = time.time() - 20
    eng.cost_usd = 0.0
    with _settings({}):
        assert eng._run_ceiling_hit(0.0) == "run_budget_s"
        eng._run_started_at = time.time()
        assert eng._run_ceiling_hit(0.0) == ""


def test_the_spend_of_the_running_turn_counts():
    eng = _engine()
    eng.run_budget_usd = 1.0
    eng.cost_usd = 0.5
    eng._run_started_at = time.time()
    with _settings({}):
        assert eng._run_ceiling_hit(0.4) == ""
        assert eng._run_ceiling_hit(0.6) == "run_budget_usd"


def test_the_engine_installs_the_probe(tmp_path):
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = E.AgentEngine(repo_dir=tmp_path, backend="api", provider="kit",
                            model="kit.qwen3.5-397b-A17b", mode="solo")

    def fine(**kw):
        yield A.StreamEvent(type="text_delta", text="ok")

    eng.client.stream_message = MagicMock(side_effect=fine)
    eng.stream_response("hi")
    probe = getattr(eng.client, "run_ceiling", None)
    assert callable(probe) and probe(0.0) == ""


def test_the_loop_asks_and_a_broken_probe_does_not_end_the_turn():
    client = A.OpenAIClient.__new__(A.OpenAIClient)
    assert client._run_ceiling_reached(0.0) == ""
    client.run_ceiling = lambda spent: "run_budget_s"
    assert client._run_ceiling_reached(0.0) == "run_budget_s"

    def boom(spent):
        raise RuntimeError("engine gone")

    client.run_ceiling = boom
    assert client._run_ceiling_reached(0.0) == ""
    import pathlib
    src = pathlib.Path(A.__file__).read_text(encoding="utf-8")
    code = "\n".join(line for line in src.splitlines()
                     if not line.lstrip().startswith("#"))
    start = code.index("for _round in range(_MAX_TOOL_ROUNDS + 1):")
    assert "_run_ceiling_reached" in code[start:start + 4000]


class _StandIn:
    session_id = "s-1"
    token_usage = {"input": 1, "output": 1}

    def __init__(self, ws):
        self.kit_permissions = A.KitToolPermissions(workspace=ws, mode="default")
        self.messages = []
        self.client = _Client()

    def set_kit_permission_mode(self, mode):
        self.kit_permissions.mode = mode
        return True

    def record_cycle_outcome(self, *a, **k):
        return None

    def export_state(self):
        return {"engine_state": {}}


@pytest.mark.parametrize("unattended", [True, False])
def test_the_headless_door_applies_them(monkeypatch, tmp_path, capsys, unattended):
    ws = tmp_path / "ws"
    ws.mkdir()
    holder = {}

    def _build(args):
        holder["e"] = _StandIn(ws)
        return holder["e"]

    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    monkeypatch.setattr(agent_cli, "_build_engine", _build)
    monkeypatch.setattr(agent_cli, "_run_once", lambda *a, **k: {
        "text": "done", "tool_calls": [], "input_tokens": 1,
        "output_tokens": 1, "error": "", "session_id": "s-1"})
    monkeypatch.setattr(agent_cli, "_save_session", lambda *a, **k: "s-1")
    monkeypatch.setattr(agent_cli, "_resume_or_create", lambda *a, **k: None)
    monkeypatch.setattr("sys.stdin", io.TextIOWrapper(io.BytesIO(b"")))
    monkeypatch.chdir(Path.cwd())
    argv = ["-p", "do it", "--cwd", str(ws)]
    if unattended:
        argv.insert(2, "--unattended")
    with _settings({}):
        rc = agent_cli.main(argv)
    assert rc == 0
    # Without a terminal nobody watches either way: both get the ceilings.
    assert getattr(holder["e"], "run_budget_s", 0) == RL.UNATTENDED_DEFAULTS["seconds"]
    assert "agent.unattended_budget" in capsys.readouterr().err


def test_an_attended_engine_is_left_alone():
    eng = _engine()
    with _settings({}):
        notes = RL.apply_unattended(eng, unattended=False)
        assert eng._run_budget() == (0.0, 0.0)
    assert notes == []
