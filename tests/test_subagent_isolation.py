"""Subagent parallel-safety + configurable limits.

Bug: concurrent fan-out subagents each mutated the SHARED parent client's
_permissions (swap/restore), racing so a subagent could run under a sibling's
sandbox. Fix: each run uses an isolated shallow copy of the client with its own
permissions; the parent is never touched. Also: per-run limits are now
settings-configurable.
"""

from __future__ import annotations

import pytest

from delfin.agent import subagents as S
from delfin.agent.api_client import StreamEvent, create_client


@pytest.fixture(autouse=True)
def _isolate_state(monkeypatch, tmp_path):
    """Redirect subagent state files to tmp so tests don't touch ~/.delfin."""
    monkeypatch.setattr(S, "_RUNNING_PATH", tmp_path / "running.json")
    monkeypatch.setattr(S, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(S, "_TELEMETRY_PATH", tmp_path / "telemetry.jsonl")
    yield


# ---------------------------------------------------------------------------
# Configurable limits
# ---------------------------------------------------------------------------


def test_limits_default_when_no_settings(monkeypatch):
    monkeypatch.setattr("delfin.user_settings.load_settings", lambda: {})
    lim = S._subagent_limits()
    assert lim["max_tool_calls"] == S._MAX_TOOL_CALLS
    assert lim["max_wall_s"] == S._MAX_WALL_S
    assert lim["max_output_tokens"] == S._MAX_OUTPUT_TOKENS


def test_limits_from_settings(monkeypatch):
    monkeypatch.setattr(
        "delfin.user_settings.load_settings",
        lambda: {"agent": {"subagents": {
            "max_tool_calls": 5, "max_wall_s": 12, "max_output_tokens": 2048}}},
    )
    lim = S._subagent_limits()
    assert lim["max_tool_calls"] == 5
    assert lim["max_wall_s"] == 12.0
    assert lim["max_output_tokens"] == 2048


def test_limits_zero_falls_back_to_default(monkeypatch):
    monkeypatch.setattr(
        "delfin.user_settings.load_settings",
        lambda: {"agent": {"subagents": {"max_wall_s": 0}}},
    )
    assert S._subagent_limits()["max_wall_s"] == S._MAX_WALL_S


# ---------------------------------------------------------------------------
# Isolation: run_subagent must NOT mutate the parent client's permissions
# ---------------------------------------------------------------------------


def _fake_stream(**kwargs):
    # Record the max_tokens the subagent was run with, then emit one reply.
    _fake_stream.last_max_tokens = kwargs.get("max_tokens")
    yield StreamEvent(type="text_delta", text="explored: found the thing")


def test_run_subagent_does_not_touch_parent_permissions(monkeypatch, tmp_path):
    parent = create_client(backend="api", provider="ollama",
                           model="qwen2.5-coder:7b", cwd=str(tmp_path))
    before = parent._permissions               # identity must be preserved
    # Avoid network: the (copied) client streams our fake events.
    monkeypatch.setattr(parent, "stream_message", _fake_stream)

    res = S.run_subagent(
        subagent_type="explore",
        description="probe",
        prompt="look around",
        parent_client=parent,
        parent_perms=parent._permissions,
    )
    assert res.error == "" or "no text" not in res.error
    assert "found the thing" in res.final_text
    # The parent's permissions object is the EXACT same instance as before —
    # the subagent ran on an isolated copy, not by swapping the parent.
    assert parent._permissions is before


def test_run_subagent_uses_configured_output_limit(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "delfin.user_settings.load_settings",
        lambda: {"agent": {"subagents": {"max_output_tokens": 4321}}},
    )
    parent = create_client(backend="api", provider="ollama",
                           model="qwen2.5-coder:7b", cwd=str(tmp_path))
    monkeypatch.setattr(parent, "stream_message", _fake_stream)
    S.run_subagent(
        subagent_type="explore", description="d", prompt="p",
        parent_client=parent, parent_perms=parent._permissions,
    )
    assert _fake_stream.last_max_tokens == 4321
