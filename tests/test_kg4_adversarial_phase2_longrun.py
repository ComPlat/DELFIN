"""KG4 adversarial, phase 2: longrun_metric against hostile inputs.

Adversarial pass over delfin/agent/longrun_metric.py (builder 1ebf8a8e),
pre-registered angles from .gate/HANDOFF.md:
  * empty / duplicate / missing-ts / non-numeric-ts records,
  * future timestamps (the reducer must not depend on wall clock),
  * text-as-tool-call disguises (quoted, unicode look-alikes),
  * steer-followed judged on the NEXT turn's action, not the steer itself,
  * unmeasurable -> None (never a silent 0),
  * duplicate compaction records must not double-count redone work.
All inputs are hand-built records in the documented writer shapes; the
module is a pure reducer, so no I/O is touched.
"""

from __future__ import annotations

import pytest

from delfin.agent.longrun_metric import LongrunInputs, summarize


def _turn(ts, *, out_chars=10, tool_calls=1, out_tokens=5, crashed=False):
    e = {"ts": ts, "output_chars": out_chars, "tool_calls": tool_calls,
         "output_tokens": out_tokens}
    if crashed:
        e["error"] = "crash"
    return e


def _cmd(ts, text):
    return {"ts": ts, "command": text}


def _steer(ts):
    return {"ts": ts}


# ---------------------------------------------------------------- empty /
# unmeasurable -> None, never a guessed zero

def test_no_records_reports_zeros_and_nones():
    r = summarize(LongrunInputs())
    assert r.turns == 0
    assert r.turns_without_progress == 0
    # Idle is unmeasurable without a next-turn observation -> None.
    assert r.idle_wait_s is None
    assert r.drift_fraction is None
    assert r.steers_followed == 0


def test_only_turns_no_next_turn_ts_idle_stays_none():
    r = summarize(LongrunInputs(turns=[_turn(100.0)]))
    assert r.idle_wait_s is None


def test_non_numeric_and_missing_ts_are_dropped_not_crashed():
    records = [{"command": "ls", "ts": "not-a-number"},
               {"command": "ls"},  # no ts at all
               _cmd(10.0, "ls -l")]
    r = summarize(LongrunInputs(commands=records,
                                compactions=[{"ts": 20.0}]))
    assert r.redone_after_compaction == 0
    assert r.compactions_seen == 1


# ------------------------------------------------------- idle & nudges

def test_idle_measured_only_when_both_ends_observed():
    r = summarize(LongrunInputs(turns=[_turn(100.0)],
                                next_turn_ts=160.0, mail_waiting=True))
    assert r.idle_wait_s == pytest.approx(60.0)
    assert r.idle_with_mail_s == pytest.approx(60.0)
    assert r.idle_with_due_event_s == 0.0


def test_nudges_inside_window_counted_outside_not():
    r = summarize(LongrunInputs(
        turns=[_turn(100.0)], next_turn_ts=160.0,
        nudges=[{"ts": 120.0}, {"ts": 159.9}, {"ts": 99.0}, {"ts": 200.0}]))
    assert r.nudges_while_waiting == 2


def test_nudge_with_unreadable_ts_ignored():
    r = summarize(LongrunInputs(
        turns=[_turn(100.0)], next_turn_ts=160.0,
        nudges=[{"ts": "soon"}, {"ts": 120.0}]))
    assert r.nudges_while_waiting == 1


# --------------------------------------------- redone after compaction

def test_redone_matches_whitespace_normalized_text():
    r = summarize(LongrunInputs(
        commands=[_cmd(10.0, "grep  -n   pattern  file"),
                  _cmd(30.0, "grep -n pattern file")],
        compactions=[{"ts": 20.0}]))
    assert r.redone_after_compaction == 1


def test_unicode_and_quoting_look_alikes_not_counted_as_redone():
    # Deliberately narrow: a look-alike/quoted variant is a DIFFERENT
    # command text, and the reducer must not claim work was redone.
    r = summarize(LongrunInputs(
        commands=[_cmd(10.0, "grep -n pattern file"),
                  _cmd(30.0, "grеp -n pattern file"),   # cyrillic е
                  _cmd(40.0, "grep -n 'pattern' file")],
        compactions=[{"ts": 20.0}]))
    assert r.redone_after_compaction == 0


def test_command_run_both_sides_of_one_compaction_counts_once():
    # One command run before AND after a compaction: the before-map has
    # one key, the after-scan runs once -> exactly 1, not 2.
    r = summarize(LongrunInputs(
        commands=[_cmd(10.0, "make test"), _cmd(30.0, "make test")],
        compactions=[{"ts": 20.0}]))
    assert r.redone_after_compaction == 1


def test_two_compactions_use_only_the_first_boundary():
    r = summarize(LongrunInputs(
        commands=[_cmd(5.0, "ls"), _cmd(40.0, "ls")],
        compactions=[{"ts": 20.0}, {"ts": 30.0}]))
    assert r.compactions_seen == 2
    assert r.redone_after_compaction == 1


# --------------------------------------------------- drift and steering

def test_drift_fraction_only_over_judged_actions():
    r = summarize(LongrunInputs(actions=[
        {"ts": 1.0, "in_goal": True},
        {"ts": 2.0, "in_goal": False},
        {"ts": 3.0},                       # unjudged -> excluded
        {"ts": 4.0, "in_goal": None},      # unjudged -> excluded
    ]))
    assert r.drift_actions == 1
    assert r.drift_fraction == pytest.approx(1 / 2)


def test_steer_followed_requires_an_in_goal_action_after_it():
    r = summarize(LongrunInputs(
        steers=[{"ts": 10.0}],
        actions=[{"ts": 5.0, "in_goal": False},    # before the steer
                 {"ts": 15.0, "in_goal": True}]))  # after -> followed
    assert r.steers == 1
    assert r.steers_followed == 1


def test_steer_with_judged_actions_only_before_is_not_followed():
    r = summarize(LongrunInputs(
        steers=[{"ts": 10.0}],
        actions=[{"ts": 5.0, "in_goal": True}]))
    assert r.steers_followed == 0


def test_no_wall_clock_or_machine_speed_dependence():
    # The same inputs, re-summarized, give the same numbers; and an all-
    # future-timestamps input is handled like any other (no "now" used).
    inputs = LongrunInputs(
        turns=[_turn(1e12), _turn(1e12 + 60)],
        commands=[_cmd(1e12, "ls"), _cmd(1e12 + 10, "ls")],
        compactions=[{"ts": 1e12 + 5}],
        steers=[{"ts": 1e12 + 20}],
        actions=[{"ts": 1e12 + 30, "in_goal": True}])
    a, b = summarize(inputs), summarize(inputs)
    assert (a.redone_after_compaction, a.steers_followed, a.compactions_seen) \
        == (b.redone_after_compaction, b.steers_followed, b.compactions_seen)
    assert a.redone_after_compaction == 1
