"""The wait side of ops health is on the product path.

Before: ``session_health.waits_summary`` and
``session_health.watch_cycle_events`` existed, were tested directly, and
nothing in the product called them — ``delfin-agent ops health`` showed
the per-session verdicts only, so a session's registered waits and a
mutual wait between two sessions were invisible to the one command whose
job is to show what is stuck. Both helpers are reached here through the
real entry point, ``cli.main(["ops", "health", ...])``.
"""
from __future__ import annotations

import json

import pytest

from delfin.agent import cli, session_health
from delfin.agent.state_paths import RedirectedUserState
from delfin.agent.wait_registry import register


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    with RedirectedUserState(tmp_path):
        yield tmp_path


@pytest.fixture
def one_pass(monkeypatch):
    """Run the --wait loop for exactly one pass: collect is emptied (its
    records are not what these tests look at) and pause_watch raises the
    KeyboardInterrupt the loop documents as its stop."""
    monkeypatch.setattr(session_health, "collect", lambda *a, **k: [])
    monkeypatch.setattr(session_health, "pause_watch",
                        lambda seconds: (_ for _ in ()).throw(
                            KeyboardInterrupt()))


def test_ops_health_prints_the_waits_table(home, capsys):
    """A registered wait reaches the `ops health` output (target 1)."""
    register("message_from", "other", key="s16",
             on_whom="session", for_what="the answer")
    code = cli.main(["ops", "health", "--json"])
    assert code == 0
    out = capsys.readouterr().out
    json.loads(out.splitlines()[-1])       # the rows stay valid JSON
    assert any("s16" in ln and "the answer" in ln for ln in
               out.splitlines()), out


def test_a_session_with_two_waits_shows_both_of_them(home):
    """One row per WAIT, not one per session.

    The shape changed with the new scan: ``waits_summary`` used to return
    one row per session key with every wait joined by "; ", and now
    returns one row per registered wait, each naming its own session. Both
    waits must be visible either way -- the reason the display exists is
    to say who waits on whom -- so the new shape is pinned rather than
    left to a joined string nobody asserted.
    """
    register("message_from", "first", key="s17",
             on_whom="session", for_what="the plan")
    register("message_from", "second", key="s17",
             on_whom="session", for_what="the number")
    rows = session_health.waits_summary()
    mine = [r for r in rows if r["session"] == "s17"]
    assert len(mine) == 2, rows
    assert {"the plan", "the number"} <= {
        r["waits"].split(": ", 1)[-1] for r in mine}, mine


def test_ops_health_wait_loop_emits_one_event_per_new_cycle(
        home, capsys, monkeypatch, one_pass):
    """A mutual wait emits exactly one attention event through the real
    emit route, and the same ring again (a second pass) is deduped."""
    register("message_from", "builder", key="qa",
             on_whom="session", for_what="review")
    register("message_from", "qa", key="builder",
             on_whom="session", for_what="hash")
    passes = {"n": 0}

    def stop_on_second(seconds):
        passes["n"] += 1
        if passes["n"] >= 2:
            raise KeyboardInterrupt()
        # first pass falls through to the real (stubbed no-sleep) pause
        return None

    monkeypatch.setattr(session_health, "pause_watch", stop_on_second)
    code = cli.main(["ops", "health", "--wait"])
    assert code == 130          # the loop's documented Ctrl-C stop
    out = capsys.readouterr().out
    lines = [ln for ln in out.splitlines() if "mutual wait" in ln]
    assert len(lines) == 1, out
    assert "qa" in lines[0] and "builder" in lines[0]


def test_ops_health_wait_loop_survives_a_broken_registry(home, capsys,
                                                         monkeypatch,
                                                         one_pass):
    """The wait side is best-effort: a raising wait_edges must not take
    the verdict loop down."""
    def broken():
        raise RuntimeError("registry gone")

    monkeypatch.setattr(session_health, "wait_edges", broken)
    code = cli.main(["ops", "health", "--json"])
    assert code == 0
    rows = json.loads(capsys.readouterr().out.splitlines()[-1])
    assert isinstance(rows, list)
