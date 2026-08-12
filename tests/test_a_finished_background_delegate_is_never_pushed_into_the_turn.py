"""A background delegate finished and nobody told the parent.

Background bash jobs get a completion event pushed into the turn at the
round boundary. Background SUB-AGENTS got nothing: no drain existed
anywhere in the engine or the tool loop. The only route back was the
parent spending a whole tool round on ``subagent_result(sa_id)`` --
against per-model round budgets of 10 to 50, and against a tool result
whose own note ends "Continue other work meanwhile". Worst case the
report was never relayed and the user never learned it existed.

The mechanism here is the one the jobs block already uses: a store that
hands each finished report to exactly ONE caller, and a per-turn block
built from what it hands over. These tests hold it to that -- a
completion appears once, an explicit collection counts as the delivery,
a run still in flight says nothing, and a delegate that died without a
report says THAT rather than staying silent.
"""

from __future__ import annotations

import json
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import engine as E
from delfin.agent import subagents as sa


@pytest.fixture(autouse=True)
def _iso(monkeypatch, tmp_path):
    monkeypatch.setattr(sa, "_RUNNING_DIR", tmp_path / "running")
    monkeypatch.setattr(sa, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(sa, "_PENDING_DIR", tmp_path / "pending")
    monkeypatch.setattr(sa, "_TELEMETRY_PATH", tmp_path / "telemetry.jsonl")


def _finish(sa_id: str, report: str, *, interactions=None) -> None:
    """End a background run the way ``run_subagent`` does: the live entry
    goes away, the report is stored."""
    sa._running_update(sa_id, None)
    sa._save_subagent_session(
        sa_id, subagent_type="explore", description="read the docs",
        messages=[{"role": "user", "content": "go"},
                  {"role": "assistant", "content": report}],
        interactions=interactions or [],
    )


@pytest.fixture
def engine(tmp_path):
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = E.AgentEngine(
            repo_dir=tmp_path, backend="api", provider="kit",
            model="kit.qwen3.5-397b-A17b", mode="solo")
    return eng


# ---------------------------------------------------------------------------
# The report reaches the turn without anyone polling for it
# ---------------------------------------------------------------------------

def test_a_finished_delegate_reaches_the_next_turn(engine):
    sa.reserve_running("bg1", subagent_type="explore", description="read docs")
    _finish("bg1", "The parser lives in delfin/agent/editblock.py.")

    block = engine._build_finished_subagents_block()
    assert "bg1" in block
    assert "editblock.py" in block, "the report itself never reached the turn"


def test_a_running_delegate_says_nothing(engine):
    sa.reserve_running("bg2", subagent_type="explore", description="read docs")
    assert engine._build_finished_subagents_block() == ""


def test_a_delegate_that_died_says_so(engine):
    sa.reserve_running("bg3", subagent_type="explore", description="read docs")
    sa.mark_running_died("bg3", "the background sub-agent thread ended with: boom")

    block = engine._build_finished_subagents_block()
    assert "bg3" in block
    assert "boom" in block
    assert "start it again" in block.lower()


def test_nothing_finished_costs_nothing(engine):
    assert engine._build_finished_subagents_block() == ""


# ---------------------------------------------------------------------------
# Exactly once -- the hard part
# ---------------------------------------------------------------------------

def test_a_completion_is_announced_once(engine):
    sa.reserve_running("bg4", subagent_type="explore", description="d")
    _finish("bg4", "Done reading; the answer is 42.")

    first = engine._build_finished_subagents_block()
    assert "bg4" in first
    assert engine._build_finished_subagents_block() == "", (
        "the same completion was announced twice")


def test_collecting_it_explicitly_counts_as_the_delivery(engine):
    sa.reserve_running("bg5", subagent_type="explore", description="d")
    _finish("bg5", "the answer is 42")

    assert sa.get_subagent_result("bg5")["status"] == "finished"
    assert engine._build_finished_subagents_block() == "", (
        "a report the parent had already collected was pushed again")


def test_the_claim_survives_a_new_process(engine, tmp_path):
    """The marker is on disk, so a restart cannot replay the report."""
    sa.reserve_running("bg6", subagent_type="explore", description="d")
    _finish("bg6", "the answer is 42")
    assert "bg6" in engine._build_finished_subagents_block()

    # A fresh engine reads the same store; the marker is gone from disk.
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        other = E.AgentEngine(
            repo_dir=tmp_path, backend="api", provider="kit",
            model="kit.qwen3.5-397b-A17b", mode="solo")
    assert other._build_finished_subagents_block() == ""


def test_one_drain_takes_a_report_and_the_other_gets_none():
    """Two callers racing must not both announce the same completion."""
    sa.reserve_running("bg7", subagent_type="explore", description="d")
    _finish("bg7", "the answer is 42")

    first = sa.drain_finished_subagents()
    second = sa.drain_finished_subagents()
    assert [r["sa_id"] for r in first] == ["bg7"]
    assert second == []


def test_more_finished_than_one_drain_carries_are_not_dropped():
    """What a drain does not hand over stays pending, it is not destroyed."""
    for i in range(sa._PENDING_DRAIN_LIMIT + 2):
        sa.reserve_running(f"m{i}", subagent_type="explore", description="d")
        _finish(f"m{i}", f"report {i}")

    seen: list[str] = []
    for _ in range(5):
        seen.extend(r["sa_id"] for r in sa.drain_finished_subagents())
    assert sorted(seen) == sorted(
        f"m{i}" for i in range(sa._PENDING_DRAIN_LIMIT + 2))
    assert len(seen) == len(set(seen)), "a report was handed over twice"


def test_a_foreground_delegate_is_not_pushed(engine):
    """Its report was already the tool result; pushing it would repeat it."""
    _finish("fg1", "I read the file and here is what it says.")
    assert engine._build_finished_subagents_block() == ""


# ---------------------------------------------------------------------------
# Wiring: the block is a drain-backed steering block like the jobs one
# ---------------------------------------------------------------------------

def test_the_block_is_registered_as_drain_backed():
    assert "finished_subagents" in E._DRAINED_STEERING_KEYS
    assert "finished_subagents" in E._MID_TURN_STEERING_KEYS


def test_it_rides_the_mid_turn_refresh(engine):
    engine._build_current_system_prompt("", task_text="los")
    sa.reserve_running("bg8", subagent_type="explore", description="d")
    _finish("bg8", "the answer is 42")

    blocks = engine._drain_turn_steering()
    assert any("bg8" in b for b in blocks)
    assert not any("bg8" in b for b in engine._drain_turn_steering())


def test_a_broken_store_does_not_end_the_turn(engine, monkeypatch):
    def boom(*a, **k):
        raise RuntimeError("registry unreadable")
    monkeypatch.setattr(sa, "drain_finished_subagents", boom)
    assert engine._build_finished_subagents_block() == ""


def test_a_long_report_is_cut_with_a_pointer_at_the_rest(engine):
    sa.reserve_running("bg9", subagent_type="explore", description="d")
    _finish("bg9", "x" * (engine._SUBAGENT_REPORT_CHARS + 5000))
    block = engine._build_finished_subagents_block()
    assert "subagent_result(sa_id='bg9')" in block
    assert len(block) < engine._SUBAGENT_REPORT_CHARS + 2000


def test_the_kept_worktree_warning_travels_with_the_report(engine):
    """A background run returns an id, so its worktree account has to be
    stored -- otherwise the one thing the parent must act on is lost."""
    sa.reserve_running("bg10", subagent_type="general-purpose",
                       description="build it")
    sa._running_update("bg10", None)
    sa._save_subagent_session(
        "bg10", subagent_type="general-purpose", description="build it",
        messages=[{"role": "assistant", "content": "built it"}],
        interactions=[],
        worktree={"warning": "the isolated worktree was NOT removed: 1 "
                             "background job(s) are still running in it (ab12)",
                  "running_jobs": ["ab12"]},
    )
    block = engine._build_finished_subagents_block()
    assert "NOT removed" in block
    assert "ab12" in block


def test_the_stored_report_is_still_json(tmp_path):
    """The extra field must not break the session store's readers."""
    sa.reserve_running("bg11", subagent_type="explore", description="d")
    _finish("bg11", "fine")
    raw = (sa._SESSIONS_DIR / "bg11.json").read_text(encoding="utf-8")
    assert json.loads(raw)["worktree"] == {}
