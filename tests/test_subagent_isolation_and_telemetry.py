"""Tests for the two subagent improvements:

- ``isolation="worktree"`` parameter on ``run_subagent`` — propagated
  through the tool schema, the dispatcher, and into a fresh git
  worktree when the parent workspace is a git repo.
- Telemetry JSONL log written by ``_write_telemetry`` + the
  ``read_telemetry`` / ``summarize_telemetry`` helpers consumed by
  ``/agents stats``.
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pytest

from delfin.agent import subagents as sa


# ---------------------------------------------------------------------------
# Telemetry helpers
# ---------------------------------------------------------------------------

@pytest.fixture
def telemetry_path(monkeypatch, tmp_path):
    """Redirect telemetry writes into a tmp dir so tests can't pollute
    the real ``~/.delfin/subagent_telemetry.jsonl``."""
    p = tmp_path / "subagent_telemetry.jsonl"
    monkeypatch.setattr(sa, "_TELEMETRY_PATH", p)
    return p


def test_write_and_read_telemetry_roundtrip(telemetry_path):
    sa._write_telemetry({
        "ts": 1.0, "subagent_type": "explore", "elapsed_s": 2.5,
        "input_tokens": 100, "output_tokens": 50,
    })
    records = sa.read_telemetry(telemetry_path)
    assert len(records) == 1
    assert records[0]["subagent_type"] == "explore"
    assert records[0]["input_tokens"] == 100


def test_telemetry_silent_on_unwritable_path(monkeypatch, tmp_path):
    """The writer must NEVER raise — telemetry is best-effort."""
    bad = tmp_path / "doesnotexist" / "nope" / "deny.jsonl"
    monkeypatch.setattr(sa, "_TELEMETRY_PATH", bad)
    # Make the parent fail to mkdir
    monkeypatch.setattr(
        Path, "mkdir",
        lambda self, **kw: (_ for _ in ()).throw(OSError("denied"))
    )
    sa._write_telemetry({"ts": 1.0, "subagent_type": "x"})
    assert not bad.exists()


def test_summarize_telemetry_groups_by_type(telemetry_path):
    sa._write_telemetry({"subagent_type": "explore", "elapsed_s": 1.0,
                         "input_tokens": 10, "output_tokens": 5})
    sa._write_telemetry({"subagent_type": "explore", "elapsed_s": 3.0,
                         "input_tokens": 20, "output_tokens": 10})
    sa._write_telemetry({"subagent_type": "plan", "elapsed_s": 5.0,
                         "input_tokens": 100, "output_tokens": 50,
                         "error": "wall-clock exhausted", "truncated": True})
    records = sa.read_telemetry(telemetry_path)
    summary = sa.summarize_telemetry(records)
    assert summary["explore"]["count"] == 2
    assert summary["explore"]["elapsed_s_total"] == 4.0
    assert summary["explore"]["input_tokens_total"] == 30
    assert summary["plan"]["errors"] == 1
    assert summary["plan"]["truncated"] == 1


def test_read_telemetry_skips_corrupt_lines(telemetry_path, tmp_path):
    telemetry_path.parent.mkdir(parents=True, exist_ok=True)
    telemetry_path.write_text(
        '{"subagent_type": "ok", "elapsed_s": 1}\n'
        "this is not json\n"
        '{"subagent_type": "ok2", "elapsed_s": 2}\n',
        encoding="utf-8",
    )
    records = sa.read_telemetry(telemetry_path)
    assert len(records) == 2
    assert records[0]["subagent_type"] == "ok"
    assert records[1]["subagent_type"] == "ok2"


def test_read_telemetry_last_n_filter(telemetry_path):
    for i in range(10):
        sa._write_telemetry({"subagent_type": "x", "i": i})
    last3 = sa.read_telemetry(telemetry_path, last_n=3)
    assert len(last3) == 3
    assert [r["i"] for r in last3] == [7, 8, 9]


# ---------------------------------------------------------------------------
# Isolation parameter — schema + runner gating (no actual worktree spawned)
# ---------------------------------------------------------------------------

def test_subagent_schema_advertises_isolation():
    from delfin.agent.api_client import _DOC_TOOLS_OPENAI
    spec = next(
        t["function"] for t in _DOC_TOOLS_OPENAI
        if t.get("function", {}).get("name") == "subagent"
    )
    iso = spec["parameters"]["properties"].get("isolation")
    assert iso is not None
    assert set(iso["enum"]) == {"", "worktree"}


def test_run_subagent_accepts_isolation_kwarg_no_op_path():
    """Passing ``isolation=""`` must keep the legacy code path — no
    worktree spawned, no error, telemetry record still appended."""

    @dataclass
    class _FakeEvent:
        type: str = ""
        text: str = ""
        tool_name: str = ""
        tool_input: dict | None = None
        input_tokens: int = 0
        output_tokens: int = 0

    class _FakeClient:
        _permissions = None

        def set_permissions(self, p): pass

        def stream_message(self, messages, system, max_tokens):
            yield _FakeEvent(type="text_delta", text="done")
            yield _FakeEvent(type="message_delta",
                             input_tokens=12, output_tokens=4)

    res = sa.run_subagent(
        subagent_type="explore",
        description="x",
        prompt="self-contained briefing of more than 20 chars",
        parent_client=_FakeClient(),
        parent_perms=None,
        isolation="",
    )
    assert res.error == ""
    assert res.final_text == "done"
    assert res.input_tokens == 12
    # Worktree summary should be empty when no isolation requested
    assert res.worktree == {}


def test_run_subagent_unknown_subagent_type_returns_error():
    """Unknown subagent_type still returns an error result and never
    writes telemetry for a bogus run (only good runs go to disk)."""
    res = sa.run_subagent(
        subagent_type="no-such-agent",
        description="x",
        prompt="self-contained briefing of more than 20 chars",
        parent_client=None,
        parent_perms=None,
    )
    assert res.error
    assert "unknown" in res.error.lower()


def test_run_subagent_worktree_isolation_gracefully_handles_no_git(
    monkeypatch, tmp_path
):
    """If isolation='worktree' is requested but the parent workspace is
    not a git repo, the runner must not crash — it just runs without
    isolation and surfaces a soft warning in the worktree summary."""

    @dataclass
    class _Perms:
        mode: str = "default"
        workspace: Path = tmp_path  # not a git repo

    @dataclass
    class _Ev:
        type: str = ""
        text: str = ""
        tool_name: str = ""
        tool_input: dict | None = None
        input_tokens: int = 0
        output_tokens: int = 0

    class _Client:
        _permissions = None

        def set_permissions(self, p): pass

        def stream_message(self, messages, system, max_tokens):
            yield _Ev(type="text_delta", text="ok")

    res = sa.run_subagent(
        subagent_type="explore",
        description="x",
        prompt="self-contained briefing of more than 20 chars",
        parent_client=_Client(),
        parent_perms=_Perms(),
        isolation="worktree",
    )
    # Either fell back to no-isolation OR returned a worktree dict
    # containing the warning. No exception either way.
    assert res.error == ""
    assert isinstance(res.worktree, dict)
