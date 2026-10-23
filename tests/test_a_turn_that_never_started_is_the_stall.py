"""The stall this module exists to expose was the one it reported as healthy.

``turn_metrics`` was built to make backend stalls visible after the fact: a
turn whose time went into waiting for the first token, not into generating or
into tool rounds. The engine stamps the first-token time only when a token
actually arrives, so a turn where the backend answered with silence records
``ttft_ms=None``.

``is_stall`` began with ``ttft is not None``. The worst possible outcome — no
first token at all, the exact failure the dashboard watchdog kills on its
first-token budget — therefore scored as "not a stall". Everything downstream
inherited that: ``format_summary`` printed the turn without its warning, and
``aggregate_turn_stats`` left it out of the ttft sample while still counting it
as a turn, so the eval report's "turns: 214, stalls: 0 / ttft: avg 1.2s, p90
2.4s" line read *cleanest* while the backend was at its worst — every turn that
never answered pulled the reported percentiles down and the stall count
nowhere. The one number the report offers as evidence for "is the backend
stalling?" improved as the backend stalled harder.

These tests pin the corrected rule: a turn that produced nothing at all counts
as a stall, provided the record can prove nothing arrived — a record that
predates the ``ttft_ms`` field must never be reinterpreted as a stall, and
neither must one where tokens demonstrably did arrive but the stamp was missed.
"""

from __future__ import annotations

import time
from pathlib import Path

import pytest

from delfin.agent import turn_metrics as tm


@pytest.fixture
def home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    monkeypatch.setattr(tm, "_DIR", tmp_path / ".delfin" / "turn_metrics")
    return tmp_path


def _never_started(total_ms: int = 90_000) -> dict:
    """The shape the engine writes when no token ever arrived."""
    return {"total_ms": total_ms, "ttft_ms": None, "output_chars": 0,
            "tool_calls": 0, "output_tokens": 0, "stopped": False,
            "error": ""}


# ---------------------------------------------------------------------------
# is_stall
# ---------------------------------------------------------------------------

def test_a_turn_that_never_produced_a_token_is_the_worst_stall_of_all():
    assert tm.is_stall(_never_started()) is True


def test_a_slow_turn_that_did_answer_is_still_a_stall_by_the_old_rule():
    assert tm.is_stall({"total_ms": 92_700, "ttft_ms": 92_000,
                        "output_chars": 6, "tool_calls": 0}) is True


def test_a_fast_turn_is_not_a_stall():
    assert tm.is_stall({"total_ms": 4_200, "ttft_ms": 800,
                        "output_chars": 6, "tool_calls": 0}) is False


def test_a_record_written_before_the_ttft_field_existed_is_never_a_stall():
    # No ttft key at all: the recorder that wrote this had no first-token
    # stamp to omit, so its silence about ttft is not evidence of silence
    # from the backend. Old logs must not be rewritten into stalls.
    legacy = {"total_ms": 90_000, "output_chars": 0, "tool_calls": 0}
    assert "ttft_ms" not in legacy
    assert tm.is_stall(legacy) is False


def test_a_turn_that_ended_at_once_without_a_token_did_not_wait_for_one():
    # Nothing arrived, but nothing was waited for either — an immediate
    # failure is a different defect and must not inflate the stall count.
    assert tm.is_stall(_never_started(total_ms=12)) is False


def test_a_turn_with_text_but_no_ttft_is_a_recording_gap_not_a_stall():
    gap = dict(_never_started(), output_chars=5_000)
    assert tm.is_stall(gap) is False


def test_a_turn_that_ran_tools_without_a_ttft_is_a_recording_gap_not_a_stall():
    gap = dict(_never_started(), tool_calls=4)
    assert tm.is_stall(gap) is False


def test_a_turn_that_streamed_only_thinking_tokens_was_not_silent():
    # Thinking deltas never stamp ttft, so a long reasoning turn looks
    # token-less in the text fields; the counted output tokens are what
    # distinguish "the backend said nothing" from "it said nothing aloud".
    thinking_only = dict(_never_started(), output_tokens=900)
    assert tm.is_stall(thinking_only) is False


def test_the_user_stopping_a_silent_turn_does_not_excuse_the_backend():
    # Whoever ended the turn, twenty seconds of silence already happened.
    assert tm.is_stall(dict(_never_started(), stopped=True)) is True


# ---------------------------------------------------------------------------
# format_summary
# ---------------------------------------------------------------------------

def test_the_per_turn_summary_warns_about_a_turn_that_never_started(home):
    tm.record("s", model="m", total_ms=90_000, ttft_ms=None, output_chars=0)
    out = tm.format_summary(tm.read("s"))
    assert "backend-stall" in out
    assert "ttft=—" in out                  # nothing to show, still flagged


# ---------------------------------------------------------------------------
# aggregate_turn_stats
# ---------------------------------------------------------------------------

def test_the_roll_up_counts_a_turn_that_never_started(home):
    tm.record("s", model="m", total_ms=90_000, ttft_ms=None, output_chars=0)
    assert tm.aggregate_turn_stats()["stalls"] == 1


def test_the_report_stops_looking_better_as_never_started_turns_pile_up(home):
    now = time.time()
    for _ in range(4):                      # a healthy baseline
        tm.record("s", model="m", total_ms=2_000, ttft_ms=500,
                  output_chars=200, tool_calls=0)
    before = tm.aggregate_turn_stats(7, now=now)
    assert before["stalls"] == 0

    for _ in range(6):                      # the backend goes silent
        tm.record("s", model="m", total_ms=90_000, ttft_ms=None,
                  output_chars=0, tool_calls=0)
    after = tm.aggregate_turn_stats(7, now=now)

    # The ttft percentiles are computed over answered turns only, so they
    # cannot get worse here — the counter that carries the damage must.
    assert after["p90_ttft_ms"] <= before["p90_ttft_ms"]
    assert after["stalls"] == 6
    assert after["stalls"] > before["stalls"]
    # Rendered side by side, the failure rate is now what moves.
    assert (after["stalls"] / after["turns"]) > (
        before["stalls"] / max(1, before["turns"]))


def test_legacy_records_do_not_turn_into_stalls_when_rolled_up(home):
    import json
    tm.record("s", model="m", total_ms=2_000, ttft_ms=500, output_chars=20)
    with tm.metrics_path("s").open("a", encoding="utf-8") as f:
        f.write(json.dumps({"ts": time.time(), "model": "m",
                            "total_ms": 90_000, "output_chars": 0,
                            "tool_calls": 0}) + "\n")
    s = tm.aggregate_turn_stats()
    assert s["turns"] == 2
    assert s["stalls"] == 0


def test_a_malformed_record_is_judged_not_a_stall_rather_than_raising():
    assert tm.is_stall({"ttft_ms": None, "total_ms": "not-a-number"}) is False
    assert tm.is_stall({}) is False
