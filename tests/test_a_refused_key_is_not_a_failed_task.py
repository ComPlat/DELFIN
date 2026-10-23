"""A turn the provider refused before any work is not a failed task.

Input: a turn that ended on a provider error before the model produced
anything. Field report 2026-10-22 (KIT cluster): the key had expired,
four turns ended on 401 AuthenticationError, and each was recorded as a
FAILED task. The next session's briefing then told the model "History:
5 outcomes, 0% pass rate · Common errors: AuthenticationError (4x) ·
Recent failure (error=AuthenticationError): Which conformation of
butane ... · WARNING: Low success rate for this task type" -- a lesson
from a login problem, quoting the user's task as a failure.

After: AuthenticationError, PermissionDeniedError, APIConnectionError,
NotFoundError and RateLimitError are not task outcomes. Such a turn is
not recorded and does not update the provider profile, and records of
them already in the history are skipped by every reader.
"""
from __future__ import annotations

import json
import textwrap
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import outcome_tracker as OT


@pytest.fixture
def history(tmp_path, monkeypatch):
    path = tmp_path / "outcome_history.jsonl"
    monkeypatch.setattr(OT, "_DEFAULT_PATH", path)
    return path


def _row(verdict, error, task="Which conformation of butane is more stable?"):
    return json.dumps({"task": task, "provider": "kit", "model": "m",
                       "mode": "solo", "verdict": verdict,
                       "cost_usd": 0.0, "duration_s": 1.0, "retries": 0,
                       "denied_commands": [], "error_type": error,
                       "task_class": "coding", "timestamp": "t"})


def test_old_records_of_a_refused_key_are_not_read(history):
    history.write_text("\n".join(
        [_row("FAIL", "AuthenticationError")] * 4
        + [_row("PASS", None, "hi")]) + "\n")
    loaded = OT.load_outcomes()
    assert [o.verdict for o in loaded] == ["PASS"]


def test_the_briefing_does_not_teach_a_login_problem(history):
    from delfin.agent.briefing import generate_briefing
    history.write_text("\n".join(
        [_row("FAIL", "AuthenticationError")] * 4
        + [_row("PASS", None, "fix the parser")] * 3) + "\n")
    text = generate_briefing("kit", "fix the parser", outcome_path=history)
    assert "AuthenticationError" not in text
    assert "Low success rate" not in text
    assert "butane" not in text


def test_a_refused_turn_is_not_recorded(history, tmp_path, monkeypatch):
    from delfin.agent import provider_profile
    from delfin.agent.engine import AgentEngine
    # The provider profile is the user's own state: never written here.
    updated = []
    monkeypatch.setattr(provider_profile, "update_from_outcome",
                        lambda provider, outcome: updated.append(outcome) or {})
    modes = tmp_path / "pack_lite" / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# quick mode")
    (tmp_path / "pack_lite" / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    client = MagicMock()
    client.model = "kit.deepseek-v4-flash"     # written to the record as text
    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=tmp_path, backend="cli", mode="quick",
                             pack_dir=tmp_path)
    for error in ("AuthenticationError", "PermissionDeniedError",
                  "APIConnectionError", "NotFoundError", "RateLimitError"):
        assert engine.record_cycle_outcome("FAIL", "a task", error) == {}
    assert not history.exists() or history.read_text().strip() == ""
    assert updated == []
    engine.record_cycle_outcome("FAIL", "a task", "crash")
    assert "crash" in history.read_text()
    assert len(updated) == 1
