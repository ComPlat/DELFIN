"""Session-end stages that fail are silent in the headless CLI.

``cmd_run`` runs six best-effort stages after the answer — session
report, skill learning, session indexing, cycle-outcome recording,
episodic memory and the eval loop. Each one is wrapped in ``except
Exception: pass``. Best-effort is right (a report write must never
break the answer's exit code), but SILENT is wrong: every scheduled
unattended run goes through here with nobody watching stderr, and a
stage that fails on every run — bad disk, broken import, a store that
moved — disappears forever. The run still prints its answer, so the
failure is indistinguishable from a stage that did its job.

These tests pin the new contract: a failing stage leaves ONE visible
stderr line naming the stage and the exception, and a failing stage
does not change the exit code or break the remaining stages.
"""

from __future__ import annotations

import io
from pathlib import Path

import pytest

from delfin.agent import cli as agent_cli


class _StubEngine:
    session_id = "sess-cli-1"
    token_usage = {"input": 1, "output": 1}

    def __init__(self):
        self.messages = [{"role": "user", "content": "do it"}]

    def stream_response(self, user_message="", **kwargs):
        self.messages.append({"role": "assistant", "content": "ANSWER"})
        return "ANSWER"

    def record_cycle_outcome(self, *a, **k):
        return None

    def export_state(self):
        return {"engine_state": {}}


_RESULT = {"text": "ANSWER", "tool_calls": [], "input_tokens": 1,
           "output_tokens": 1, "error": "", "session_id": "sess-cli-1"}


def _run_cmd_run(monkeypatch, tmp_path, failing: str):
    """Run cmd_run with the named session-end stage raising."""
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    monkeypatch.setattr(agent_cli, "_build_engine",
                        lambda args: _StubEngine())
    monkeypatch.setattr(agent_cli, "_run_once",
                        lambda engine, prompt, **kw: dict(_RESULT))
    monkeypatch.setattr(agent_cli, "_save_session",
                        lambda *a, **k: "sess-cli-1")

    stages = {
        "report": "delfin.agent.session_report.write_session_report",
        "learn": "delfin.agent.session_end.learn_at_session_end",
        "index": "delfin.agent.session_end.index_at_session_end",
        "outcome": "delfin.agent.cli.<engine>record_cycle_outcome",
        "episode": "delfin.agent.episodes.save_episode",
        "eval": "delfin.agent.eval_loop.maybe_run_scheduled",
    }
    target = stages[failing]

    def _boom(*a, **k):
        raise OSError(f"simulated {failing} stage failure")

    if failing == "outcome":
        engine = _StubEngine()
        engine.record_cycle_outcome = _boom
        monkeypatch.setattr(agent_cli, "_build_engine",
                            lambda args: engine)
    else:
        monkeypatch.setattr(target, _boom)

    monkeypatch.setattr("sys.stdin", io.TextIOWrapper(io.BytesIO(b"")))
    rc = agent_cli.main(["-p", "do it"])
    return rc


@pytest.mark.parametrize("failing", [
    "report", "learn", "index", "outcome", "episode", "eval",
])
def test_a_failing_session_end_stage_is_said_on_stderr(
        monkeypatch, tmp_path, capsys, failing):
    rc = _run_cmd_run(monkeypatch, tmp_path, failing)
    captured = capsys.readouterr()
    assert rc == 0, "a session-end stage failure must not change the exit code"
    combined = captured.err + captured.out
    assert "session-end" in combined, \
        f"stage {failing!r} failed silently: no stderr line at all"
    # The stage's human name is what the warning carries (learn = "skill
    # learning", outcome = "cycle outcome", episode = "episodic memory");
    # the simulated exception text also echoes the test's key.
    assert f"simulated {failing} stage failure" in combined


def test_remaining_stages_still_run_after_one_fails(
        monkeypatch, tmp_path, capsys):
    """One stage failing must not short-circuit the ones after it."""
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    monkeypatch.setattr(agent_cli, "_build_engine",
                        lambda args: _StubEngine())
    monkeypatch.setattr(agent_cli, "_run_once",
                        lambda engine, prompt, **kw: dict(_RESULT))
    monkeypatch.setattr(agent_cli, "_save_session",
                        lambda *a, **k: "sess-cli-1")

    def _boom(*a, **k):
        raise OSError("simulated report stage failure")

    monkeypatch.setattr(
        "delfin.agent.session_report.write_session_report", _boom)

    calls: list[str] = []

    def fake_learn(msgs, *, session_id="", **kw):
        calls.append("learn")
        return None

    monkeypatch.setattr(
        "delfin.agent.session_end.learn_at_session_end", fake_learn)
    monkeypatch.setattr("sys.stdin", io.TextIOWrapper(io.BytesIO(b"")))
    rc = agent_cli.main(["-p", "do it"])
    captured = capsys.readouterr()
    assert rc == 0
    assert calls == ["learn"], "the next stage still runs"
    assert "session-end" in captured.err


def test_healthy_run_stays_quiet(monkeypatch, tmp_path, capsys):
    """The normal path gains no new machinery speech: a run where every
    session-end stage works prints nothing about session-end."""
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    monkeypatch.setattr(agent_cli, "_build_engine",
                        lambda args: _StubEngine())
    monkeypatch.setattr(agent_cli, "_run_once",
                        lambda engine, prompt, **kw: dict(_RESULT))
    monkeypatch.setattr(agent_cli, "_save_session",
                        lambda *a, **k: "sess-cli-1")
    monkeypatch.setattr("sys.stdin", io.TextIOWrapper(io.BytesIO(b"")))
    rc = agent_cli.main(["-p", "do it"])
    captured = capsys.readouterr()
    assert rc == 0
    assert "session-end" not in captured.err + captured.out
