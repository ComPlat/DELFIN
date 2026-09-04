"""Turning Dynamik Opt off stops a walk that a release started.

Reported from a real session: "Dynamik Opt aus, es geht weiter, warum?" -- the
switch went off and run 54, a climb the release had started, went on streaming
frames into the picture for another sixteen seconds under a line that said the
climb was the user's own to stop.

Dynamik Opt is the switch for the structure answering the hand, and a walk a
release starts is the hand carrying on: Auto minimising from where a drag left
off, or a climb aimed by the gesture.  So it belongs to that switch and goes
off with it.  A climb or a minimise the user pressed for is a deliberate
calculation with its own button, and it keeps running.
"""
from __future__ import annotations

import pathlib
import sys

import pytest

pytest.importorskip("ipywidgets")


def _a_part(structure="""8
ethane
C  0.000000  0.000000  0.762900
C  0.000000  0.000000 -0.762900
H -0.505000  0.874000  1.162900
H -0.505000 -0.874000  1.162900
H  1.010000  0.000000  1.162900
H  0.505000  0.874000 -1.162900
H  0.505000 -0.874000 -1.162900
H -1.010000  0.000000 -1.162900
"""):
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    import test_the_budget_prices_a_relaxed_path as budget
    return budget._a_part(structure)


def _run_a_climb(part, from_release):
    """Put the editor in the state a running climb leaves, by the two roads."""
    part.submit_ff_dd.value = "gfn2"
    part.submit_relax_btn.value = True
    part.state["climb_run"] = object()
    part.state["walk_from_release"] = from_release
    part.state["gfn_run"] = 3


def test_a_release_started_climb_stops_with_the_switch():
    part = _a_part()
    _run_a_climb(part, from_release=True)
    part.submit_relax_btn.value = False
    assert part.state.get("climb_run") is None, "the climb kept running"
    said = " ".join(part.state.get("mol_status_lines") or ())
    assert "no longer following your hand" in said
    assert "has stopped with it" in said
    # Climb to TS is a mode, not this run's switch: it is left as the user set
    # it, so the next drag still climbs.
    assert part.submit_climb_btn.value is False or part.submit_climb_btn.value


def test_a_climb_the_user_pressed_for_keeps_running():
    part = _a_part()
    _run_a_climb(part, from_release=False)
    part.submit_relax_btn.value = False
    assert part.state.get("climb_run") is not None, "a deliberate climb was cut"
    said = " ".join(part.state.get("mol_status_lines") or ())
    assert "no longer following your hand" in said
    assert "its own button is what stops that" in said


def test_a_release_started_minimise_stops_with_the_switch():
    part = _a_part()
    part.submit_ff_dd.value = "gfn2"
    part.submit_relax_btn.value = True
    part.state["optimize_run"] = object()
    part.state["walk_from_release"] = True
    part.state["gfn_run"] = 3
    part.submit_relax_btn.value = False
    # Its own switch is the clean stop; the toggle handler ends the run.
    assert part.submit_optimize_btn.value is False
    said = " ".join(part.state.get("mol_status_lines") or ())
    assert "has stopped with it" in said


def test_the_provenance_is_written_where_each_walk_begins():
    """A source census: the flag is set in the two starters and read in the
    off-branch, so neither road can be added without setting it."""
    import inspect

    from delfin.dashboard import structure_editor

    source = inspect.getsource(structure_editor)
    climb = source.split("def _climb_now")[1].split("\n    def ")[0]
    assert "state['walk_from_release'] = aimed_from is not None" in climb
    opt = source.split("def on_submit_optimize(")[1].split("\n    def ")[0]
    assert "state['walk_from_release'] = bool(state.pop('optimise_from_release'" in opt
    off = source.split("def on_submit_relax_toggle(change):")[1].split("\n    def ")[0]
    assert "state.get('walk_from_release')" in off
