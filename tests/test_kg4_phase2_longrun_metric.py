"""KG4 phase 2: longrun_metric — a pure reducer over a session's records.

Contract (from the autonomy_metric it builds on, autonomy_metric.py:207-230):
pure, no I/O, real records only, anything the records cannot attest is
reported as None/UNCLASSIFIED — never guessed. Every input field names its
real writer; a test feeds records produced by the REAL writers
(turn_metrics.record, compaction_log.record_compaction), never hand-typed
shadows.
"""
from __future__ import annotations

import json
import time

import pytest


def _real_turn_entries(tmp_path, session, entries):
    """Feed entries through the REAL writer turn_metrics.record and read
    them back — the reducer must survive real shapes, not typed dicts."""
    from delfin.agent import turn_metrics as tm
    for e in entries:
        tm.record(session, **e)
    return tm.read(session)


def _real_compaction_records(tmp_path, session, events):
    """Feed events through the REAL writer compaction_log.record_compaction
    and read them back via read_compactions."""
    from delfin.agent import compaction_log as cl
    out = []
    for ev in events:
        line = cl.record_compaction(session, **ev)
        assert line, "record_compaction returned None — nothing recorded"
    out = cl.read_compactions(session)
    assert out, "read_compactions returned nothing after real writes"
    return out


@pytest.fixture
def redirected_home(tmp_path, monkeypatch):
    """Redirect user-state into the test's scratch home with the real
    mechanism (state_paths.RedirectedUserState), so turn_metrics and
    compaction_log write where the reducer reads them under test."""
    from delfin.agent import state_paths
    home = tmp_path / "home"
    home.mkdir(exist_ok=True)
    monkeypatch.setenv("HOME", str(home))
    redirect = state_paths.RedirectedUserState(home)
    redirect.__enter__()
    try:
        yield home
    finally:
        redirect.__exit__(None, None, None)


def test_idle_time_with_mail_waiting_is_measured(redirected_home):
    """Idle window with mail waiting: the idle time is measured from the
    records (last turn end -> next turn start), and counted as waiting
    because mail was deliverable during it."""
    from delfin.agent import longrun_metric as lm

    turns = _real_turn_entries(redirected_home, "s-idle", [
        {"provider": "kit", "model": "m", "role": "s16",
         "total_ms": 500, "output_chars": 20, "output_tokens": 5},
    ])
    last_ts = turns[-1]["ts"]
    idle = lm.LongrunInputs(
        turns=turns,
        next_turn_ts=last_ts + 900.0,
        mail_waiting=True,          # writer: session_messages deliverable()
        due_event_waiting=False,
    )
    res = lm.summarize(idle)
    assert res.idle_wait_s == pytest.approx(900.0)
    assert res.idle_with_mail_s == pytest.approx(900.0)


def test_idle_without_mail_or_event_is_not_counted_as_waiting(
        redirected_home):
    from delfin.agent import longrun_metric as lm
    turns = _real_turn_entries(redirected_home, "s-idle2", [
        {"provider": "kit", "model": "m", "role": "s16",
         "total_ms": 100, "output_chars": 3, "output_tokens": 1},
    ])
    res = lm.summarize(lm.LongrunInputs(
        turns=turns, next_turn_ts=turns[-1]["ts"] + 600.0,
        mail_waiting=False, due_event_waiting=False))
    assert res.idle_wait_s == pytest.approx(600.0)
    assert res.idle_with_mail_s == 0.0
    assert res.idle_with_due_event_s == 0.0


def test_turns_without_progress_and_text_calls_counted(redirected_home):
    """A turn that delivered nothing (real writer's own definition via
    turn_metrics helpers) is a no-progress turn; a 'tool call written as
    text' is one where the model emitted the call as prose, recorded by
    the caller's probe."""
    from delfin.agent import longrun_metric as lm
    turns = _real_turn_entries(redirected_home, "s-prog", [
        {"provider": "kit", "model": "m", "role": "s16",
         "total_ms": 100, "output_chars": 50, "output_tokens": 10,
         "tool_calls": 1},
        {"provider": "kit", "model": "m", "role": "s16",
         "total_ms": 100, "output_chars": 0, "output_tokens": 0},
        {"provider": "kit", "model": "m", "role": "s16",
         "total_ms": 100, "output_chars": 0, "output_tokens": 0},
    ])
    res = lm.summarize(lm.LongrunInputs(turns=turns, text_calls=[{"ts": 1.0}]))
    assert res.turns == 3
    assert res.turns_without_progress == 2
    assert res.text_calls == 1


def test_nudges_while_waiting_are_counted(redirected_home):
    from delfin.agent import longrun_metric as lm
    turns = _real_turn_entries(redirected_home, "s-nudge", [
        {"provider": "kit", "model": "m", "role": "s16",
         "total_ms": 100, "output_chars": 5, "output_tokens": 1},
    ])
    res = lm.summarize(lm.LongrunInputs(
        turns=turns, next_turn_ts=turns[-1]["ts"] + 300.0,
        mail_waiting=True,
        nudges=[{"ts": turns[-1]["ts"] + 100.0, "kind": "wait"},
                {"ts": turns[-1]["ts"] + 200.0, "kind": "wait"}]))
    assert res.idle_with_mail_s == pytest.approx(300.0)
    assert res.nudges_while_waiting == 2


def test_work_redone_after_compaction_is_counted(redirected_home):
    """Work redone after a compaction or restart: the same command text
    appearing again AFTER a compaction event is counted as redone work."""
    from delfin.agent import longrun_metric as lm
    compactions = _real_compaction_records(redirected_home, "s-redo", [
        {"kind": "auto", "tokens_before": 500000, "tokens_after": 100000,
         "messages_compacted": 10},
    ])
    t_comp = compactions[0]["ts"]
    turns = _real_turn_entries(redirected_home, "s-redo", [
        {"provider": "kit", "model": "m", "role": "s16", "total_ms": 100,
         "output_chars": 5, "output_tokens": 1, "tool_calls": 1},
    ])
    commands = [
        {"ts": t_comp - 10.0, "command": "pytest tests/test_x.py -q"},
        {"ts": t_comp + 10.0, "command": "pytest tests/test_x.py -q"},
        {"ts": t_comp + 20.0, "command": "pytest tests/test_x.py -q"},
        {"ts": t_comp + 30.0, "command": "ruff check delfin/"},
    ]
    res = lm.summarize(lm.LongrunInputs(
        turns=turns, compactions=compactions, commands=commands))
    # commands 2 and 3 repeat the pre-compaction command after the
    # compaction boundary -> redone work counted twice.
    assert res.redone_after_compaction == 2
    assert res.compactions_seen == 1


def test_drift_and_steer_followed(redirected_home):
    """Drift = actions outside the stated goal (the caller supplies the
    per-command in_goal verdict from its own goal statement — the reducer
    never guesses it). A steer is 'followed' when the next turn's action
    is back in the goal after the steer message."""
    from delfin.agent import longrun_metric as lm
    turns = _real_turn_entries(redirected_home, "s-steer", [
        {"provider": "kit", "model": "m", "role": "s16", "total_ms": 100,
         "output_chars": 5, "output_tokens": 1},
        {"provider": "kit", "model": "m", "role": "s16", "total_ms": 100,
         "output_chars": 5, "output_tokens": 1},
        {"provider": "kit", "model": "m", "role": "s16", "total_ms": 100,
         "output_chars": 5, "output_tokens": 1},
    ])
    actions = [
        {"ts": 1.0, "in_goal": True},
        {"ts": 2.0, "in_goal": False},
        {"ts": 3.0, "in_goal": True},
    ]
    res = lm.summarize(lm.LongrunInputs(
        turns=turns, actions=actions,
        steers=[{"ts": 2.5, "text": "focus back on the package"}]))
    assert res.drift_actions == 1
    assert res.drift_fraction == pytest.approx(1.0 / 3.0)
    # the steer at ts 2.5 was followed: the next in-goal action (ts 3.0)
    # came after it
    assert res.steers == 1
    assert res.steers_followed == 1


def test_pure_no_io_and_unclassifiable_is_none_not_guessed():
    """The reducer is pure (no file access) and never invents figures:
    an input it cannot attest (no turns at all) reports None, not 0."""
    import inspect
    from delfin.agent import longrun_metric as lm
    src = inspect.getsource(lm)
    for banned in ("open(", "Path(", "os.remove", "shutil"):
        assert banned not in src, f"longrun_metric must stay pure, found {banned}"
    res = lm.summarize(lm.LongrunInputs(turns=[], next_turn_ts=None))
    assert res.idle_wait_s is None
    assert res.turns == 0
