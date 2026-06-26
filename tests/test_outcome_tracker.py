"""Tests for delfin.agent.outcome_tracker — append + load + B3 update."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent.outcome_tracker import (
    CycleOutcome,
    append_outcome,
    load_outcomes,
    update_last_outcome,
)


@pytest.fixture
def tmp_history(tmp_path):
    return tmp_path / "outcomes.jsonl"


# ---------------------------------------------------------------------------
# B3 — update_last_outcome (in-place mutation of the most recent entry)
# ---------------------------------------------------------------------------

def test_update_last_outcome_bumps_retries(tmp_history):
    """A retry attempt updates the existing row's retries counter."""
    append_outcome(
        CycleOutcome(task="t1", verdict="PASS", cost_usd=0.1,
                     cost_usd_delta=0.1, timestamp="2026-04-28T08:00:00"),
        path=tmp_history,
    )
    ok = update_last_outcome(retries=1, path=tmp_history)
    assert ok is True
    outcomes = load_outcomes(path=tmp_history)
    assert len(outcomes) == 1
    assert outcomes[0].retries == 1


def test_update_last_outcome_changes_verdict(tmp_history):
    """A retry that fails can change PASS → FAIL on the same row."""
    append_outcome(
        CycleOutcome(task="t1", verdict="PASS", cost_usd=0.1,
                     cost_usd_delta=0.1, timestamp="2026-04-28T08:00:00"),
        path=tmp_history,
    )
    update_last_outcome(retries=1, verdict="FAIL", path=tmp_history)
    outcomes = load_outcomes(path=tmp_history)
    assert outcomes[0].verdict == "FAIL"
    assert outcomes[0].retries == 1


def test_update_last_outcome_only_touches_provided_fields(tmp_history):
    """Other fields stay as they were."""
    append_outcome(
        CycleOutcome(task="t1", verdict="PASS", cost_usd=0.1,
                     cost_usd_delta=0.1, duration_s=42.0,
                     mode="solo", provider="claude",
                     timestamp="2026-04-28T08:00:00"),
        path=tmp_history,
    )
    update_last_outcome(retries=2, path=tmp_history)
    outcomes = load_outcomes(path=tmp_history)
    o = outcomes[0]
    assert o.task == "t1"
    assert o.cost_usd == 0.1
    assert o.duration_s == 42.0
    assert o.mode == "solo"
    assert o.provider == "claude"


def test_update_last_outcome_only_modifies_last_row(tmp_history):
    """Earlier rows are untouched when we bump the latest one."""
    append_outcome(
        CycleOutcome(task="first", verdict="PASS", cost_usd=0.1,
                     cost_usd_delta=0.1, timestamp="2026-04-28T08:00:00"),
        path=tmp_history,
    )
    append_outcome(
        CycleOutcome(task="second", verdict="PASS", cost_usd=0.3,
                     cost_usd_delta=0.2, timestamp="2026-04-28T08:05:00"),
        path=tmp_history,
    )
    update_last_outcome(retries=1, verdict="FAIL", path=tmp_history)
    outcomes = load_outcomes(path=tmp_history)
    assert outcomes[0].task == "first"
    assert outcomes[0].retries == 0
    assert outcomes[0].verdict == "PASS"
    assert outcomes[1].task == "second"
    assert outcomes[1].retries == 1
    assert outcomes[1].verdict == "FAIL"


def test_update_last_outcome_returns_false_for_missing_file(tmp_path):
    """No file → no-op, False return."""
    missing = tmp_path / "nope.jsonl"
    assert update_last_outcome(retries=1, path=missing) is False


def test_update_last_outcome_returns_false_for_empty_file(tmp_history):
    tmp_history.write_text("", encoding="utf-8")
    assert update_last_outcome(retries=1, path=tmp_history) is False


def test_update_last_outcome_handles_trailing_blank_line(tmp_history):
    """append_outcome leaves a trailing newline — last_idx must skip it."""
    append_outcome(
        CycleOutcome(task="real", verdict="PASS", cost_usd=0.1,
                     cost_usd_delta=0.1, timestamp="2026-04-28T08:00:00"),
        path=tmp_history,
    )
    # Add an extra trailing blank line manually
    tmp_history.write_text(
        tmp_history.read_text(encoding="utf-8") + "\n",
        encoding="utf-8",
    )
    update_last_outcome(retries=3, path=tmp_history)
    outcomes = load_outcomes(path=tmp_history)
    assert len(outcomes) == 1
    assert outcomes[0].retries == 3


def test_update_last_outcome_no_args_is_noop(tmp_history):
    """Calling with no kwargs returns False (nothing to update)."""
    append_outcome(
        CycleOutcome(task="t", verdict="PASS", cost_usd_delta=0.1,
                     timestamp="2026-04-28T08:00:00"),
        path=tmp_history,
    )
    assert update_last_outcome(path=tmp_history) is False


def test_update_last_outcome_recovers_from_corrupt_last_line(tmp_history):
    """A non-JSON last line returns False instead of raising."""
    tmp_history.write_text("this is not valid json\n", encoding="utf-8")
    assert update_last_outcome(retries=1, path=tmp_history) is False


def test_update_last_outcome_preserves_legacy_fields(tmp_history):
    """An entry without cost_usd_delta (legacy) keeps deserialising fine
    after an in-place update of unrelated fields."""
    # Hand-craft a legacy line lacking cost_usd_delta
    legacy = {
        "task": "old", "provider": "claude", "model": "opus",
        "mode": "solo", "verdict": "PASS", "cost_usd": 0.5,
        "duration_s": 30.0, "retries": 0, "denied_commands": [],
        "error_type": None, "task_class": "coding",
        "timestamp": "2026-04-15T08:00:00",
    }
    tmp_history.write_text(json.dumps(legacy) + "\n", encoding="utf-8")
    update_last_outcome(retries=2, path=tmp_history)
    outcomes = load_outcomes(path=tmp_history)
    assert outcomes[0].retries == 2
    # cost_usd_delta defaults to 0.0 because the field was absent
    assert outcomes[0].cost_usd_delta == 0.0
    # Other legacy fields remain
    assert outcomes[0].task == "old"
    assert outcomes[0].cost_usd == 0.5
