"""Under the line you type: which posture you are in, and what is still out.

The posture was in the banner, which scrolls away after the first turn,
and the hint row under the prompt said only which keys exist. So the one
fact that decides what the next command will do — whether it asks, or
edits, or asks nothing at all — was off screen for the whole session.

    acceptEdits · esc interrupt · shift+tab approval mode · /help

It costs no row: the hint row was always there.

And what is running had a row of its own, which named the work and how
long it had been out but not the handle to reach it. `/bash <id>` was the
way in, and the id was the one thing on screen that was missing:

    ⚙ pytest tests/ 6m40s [01e5b151] · +2 more

Deliberately left for later: walking the list with the arrow keys and
attaching to one. This is the half that makes the other half possible —
you cannot reach a job whose id you were never shown.
"""

from __future__ import annotations

import pytest

from delfin.agent import background_view as BV


def _view(shells=(), agents=(), watches=()):
    return {"shells": list(shells), "agents": list(agents),
            "watches": list(watches), "wakeups": [], "errors": []}


# -- the posture ------------------------------------------------------------

def test_the_hint_row_names_the_posture():
    from delfin.agent import repl as R
    line = R._box_hint("acceptEdits")
    assert line.startswith("acceptEdits"), line
    assert "esc interrupt" in line, "the keys are still there"


def test_bypass_is_named_too():
    from delfin.agent import repl as R
    assert R._box_hint("bypassPermissions").startswith("bypassPermissions")


def test_no_posture_leaves_the_row_as_it_was():
    from delfin.agent import repl as R
    assert R._box_hint("") == R._BOX_HINT


# -- the handle -------------------------------------------------------------

def test_a_running_shell_shows_the_id_you_reach_it_by():
    line = BV.status_line(
        _view(shells=[{"id": "01e5b151", "label": "pytest tests/",
                       "since": 0}]), now=400)
    assert "01e5b151" in line, line
    assert "pytest tests/" in line


def test_the_id_is_short_enough_to_read():
    line = BV.status_line(
        _view(shells=[{"id": "01e5b151c0ffee", "label": "x", "since": 0}]),
        now=60)
    assert "01e5b151c0ffee" not in line, "the full handle is not a glance"
    assert "01e5b151" in line


def test_nothing_out_still_costs_no_row():
    assert BV.status_line(_view()) == ""


def test_it_still_never_raises():
    assert BV.status_line({"shells": "not a list"}) == ""


def test_the_count_of_the_rest_survives():
    line = BV.status_line(
        _view(shells=[{"id": f"id{i:06d}", "label": f"job {i}", "since": 0}
                      for i in range(9)]), now=60)
    assert "+6 more" in line, line
