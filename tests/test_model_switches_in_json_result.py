"""``_run_once`` reports the model switches a turn made (PR #448 records
them on the client as ``_turn_model_switches``; each entry is
``{before, after, reason}``). A headless ``chat -p --output-format json``
result carries them under ``model_switches``.

The turn that switches models is the client the engine owns, so the value
is read from ``engine.client`` -- and a stub engine without a client (the
scheduler and benchmark test doubles) must keep working and report ``[]``.
"""
import json

from delfin.agent import cli as agent_cli
from delfin.agent import repl


class _StubClient:
    def __init__(self, switches):
        self._turn_model_switches = list(switches)


class _EngineWithClient:
    """Mirrors AgentEngine: token_usage + stream_response + a .client
    whose switches the headless result should report."""
    token_usage = {"input": 10, "output": 5}

    def __init__(self, switches):
        self.client = _StubClient(switches)

    def stream_response(self, *, user_message, on_token, on_tool_use,
                        on_permission_denied=None, max_tokens=4096,
                        **kwargs):
        on_token("one switch this turn")
        return "one switch this turn"


class _EngineWithoutClient:
    """A stub engine that never had a client (a scheduler double): the
    switch read must not reach for an attribute that is absent."""
    token_usage = {"input": 0, "output": 0}

    def stream_response(self, *, user_message, on_token, on_tool_use,
                        max_tokens=4096, **kwargs):
        on_token("no switch this turn")
        return "no switch this turn"


def test_run_once_reports_one_switch_from_the_client():
    from delfin.agent.cli import _run_once

    switch = {"before": "old", "after": "new", "reason": "model unavailable"}
    engine = _EngineWithClient([switch])
    out = _run_once(engine, "probe it")
    assert out["model_switches"] == [switch]


def test_run_once_reports_no_switch_as_an_empty_list():
    """A turn without a switch reports [], never null or an absent key:
    the JSON result must be a list for code that reads it."""
    from delfin.agent.cli import _run_once

    engine = _EngineWithClient([])
    out = _run_once(engine, "probe it")
    assert out["model_switches"] == []


def test_run_once_without_a_client_still_reports_an_empty_list():
    """A stub engine without a client must keep working and report []."""
    from delfin.agent.cli import _run_once

    engine = _EngineWithoutClient()
    out = _run_once(engine, "probe it")
    assert out["model_switches"] == []


def test_model_switches_reaches_the_json_one_shot_result(
        monkeypatch, capsys, tmp_path):
    """End to end through --output-format json: the switch is in the
    printed payload, keyed by repl.TURN_KEYS."""
    from pathlib import Path

    switch = {"before": "old", "after": "new", "reason": "model not found"}
    engine = _EngineWithClient([switch])
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    monkeypatch.setattr(agent_cli, "_build_engine", lambda args: engine)
    monkeypatch.setattr(agent_cli, "_save_session",
                        lambda e, root, **kw: "sid-1")
    _no_tty(monkeypatch)

    agent_cli.main(["-p", "do it", "--output-format", "json"])
    payload = json.loads(capsys.readouterr().out)
    assert payload["model_switches"] == [switch]
    assert payload["model_switches"] is not None
    assert set(payload) == set(repl.TURN_KEYS) | {"session_id"}


def test_run_once_return_keys_match_turn_keys():
    """The contract that _run_once's return is exactly repl.TURN_KEYS; a
    key outside it breaks the pinned JSON tests already in the suite."""
    from delfin.agent.cli import _run_once

    engine = _EngineWithClient([])
    out = _run_once(engine, "probe it")
    assert set(out) == set(repl.TURN_KEYS)


def _no_tty(monkeypatch):
    import io

    stream = io.StringIO("")
    stream.isatty = lambda: False
    monkeypatch.setattr("sys.stdin", stream)
