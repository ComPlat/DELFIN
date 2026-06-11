"""Shared test fixtures.

Redirect the subagent state files (live registry, telemetry, finished
sessions) away from the real ``~/.delfin`` so test runs never leave
artifacts that show up in the user's dashboard (live panel, ``/agents``
listing, stats).
"""

from __future__ import annotations

import pytest


@pytest.fixture(autouse=True)
def _isolate_subagent_state(tmp_path, monkeypatch):
    from delfin.agent import subagents as sa
    monkeypatch.setattr(sa, "_RUNNING_PATH",
                        tmp_path / "subagent_running.json")
    monkeypatch.setattr(sa, "_TELEMETRY_PATH",
                        tmp_path / "subagent_telemetry.jsonl")
    monkeypatch.setattr(sa, "_SESSIONS_DIR",
                        tmp_path / "subagent_sessions")
    yield
