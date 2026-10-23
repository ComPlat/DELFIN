"""Letting go of an atom must not rebuild the editor's world.

Two things met here and both were mine.

**The picture.**  The release draws whatever the wall kept, because the page's
model is where the *cursor* was and the wall's is where the chemistry allowed.
Where the box already holds the kept geometry -- which at a release it almost
always does -- there was nothing to write, so it redrew by calling the host's
``update_view()``.  That is the host's "a structure changed" path however it is
flagged: it throws the RDKit perception away and builds it again, and it ends
in a whole new viewer -- a WebGL context disposed and another made, of which a
browser grants a handful.  One per release is what every drag was paying.

**The budget.**  Its anchor is guarded against everything an energy depends on,
which is right, and a mismatch *dropped* it, which is the dangerous way to be
right: the switch stays lit, the drag goes on costing what a priced drag costs,
and nothing is refused any more.  Survivable when the user moved the question.
Not survivable when the editor moved it -- ``auto M`` settles on a spin that
every later run then uses, so one Optimise takes the wall out of force with
nobody asking.

Measured on an ethane, one hydrogen pulled fourteen times by 0.45 A under GFN2
with the budget lit and a solvent picked after Set: the wall in force leaves
C-H at 1.088 A with all seven bonds and reads "+16.0 of 22.3 kcal/mol"; the
wall out of force leaves 1.307 A, because at the release nothing is kept and
what stays is the cursor's wish.  Pull further and that is a proton off.
"""
from __future__ import annotations

import pathlib
import re

import pytest

_EDITOR = (pathlib.Path(__file__).resolve().parents[1]
           / 'delfin' / 'dashboard' / 'structure_editor.py').read_text()


def _body(name):
    return _EDITOR.split(f'def {name}(')[1].split('\n    def ')[0]


def test_the_release_swaps_the_model_it_does_not_rebuild_the_viewer():
    """The cheap half exists and was written for exactly this."""
    body = _body('_keep_this_geometry')
    swap = body.index('_swap_the_model_js(xyz_document(rows, why))')
    host = body.index('update_view()')
    assert swap < host, 'the host is asked first, so nothing is ever swapped'
    # And the host stays as the fallback, because without a live viewer to
    # swap into there is nothing else to draw with.
    assert 'if swap:' in body
    assert '_run_manip_js(swap)' in body


def test_an_anchor_that_outlived_its_question_is_measured_again():
    """Rather than let go of, which is the wall going out with the light on."""
    body = _body('_keep_the_anchor_current')
    assert '_set_thermal_anchor(' in body
    assert 'submit_thermal_btn.value' in body
    assert '_the_anchor_outlived_its_question()' in body
    # Two ways a lit budget ends up with no zero it can use, and both are
    # the editor's doing rather than the user's, so both are put right the
    # same way.  It is *stale* -- auto M settled on a spin, or a solvent was
    # chosen after Set -- or it is *gone*, because an Optimise wrote its
    # result back and structure_changed dropped the anchor with the structure
    # it belonged to.  The second was left out at first, and measured: the
    # budget on, one Optimise all pressed, and the next grab ran unwalled
    # with the button still lit.
    told = _body('_the_anchor_outlived_its_question')
    assert "if state.get('thermal_e0') is None:" in told
    assert 'return bool(_current_xyz())' in told
    assert '_the_question_an_anchor_answers()' in told


def test_it_is_asked_where_the_wall_is_needed():
    """At the grab, which happens once and is where a hand arrives.

    Hooked at the presses that move the question instead, it would be right
    for the press that moved it and miss every other way -- and the way it
    goes wrong is the molecule coming apart.
    """
    grab = _EDITOR.split("if verb == 'grabbed':")[1]
    grab = grab.split('\n        if verb')[0]
    assert '_keep_the_anchor_current(' in grab


def test_the_line_no_longer_blames_the_method_for_a_multiplicity():
    """Four things can move the question and the sentence named all four.

    The method is the one this can tell apart for certain, so it is named on
    its own and the other three are named as the three they are.
    """
    body = _body('_thermal_note')
    assert "state.get('thermal_method') != str(submit_ff_dd.value)" in body
    assert 'its zero was measured ' in body
    assert '_server_label(str(state.get(' in body
    # And the other three are named as the three they are.
    assert 'charge, multiplicity or solvent' in body
    assert 'method, charge, multiplicity or solvent' not in body


def test_a_follow_that_has_answered_owns_the_box_wall_or_no_wall():
    """The release must not put the cursor's wish over a computed geometry.

    The wall's branch decides *which* geometry to keep.  Underneath it sits
    the plainer question: whether the last word on a drag may be a geometry
    nothing computed.  The page's payload is where the cursor was, and the
    release is the one message that carries it with no answer on its way to
    overwrite it -- so one frame undid every relaxed answer the drag had made.

    Measured on an ethane with the budget lit and the engine moved from GFN2
    to GFN-FF after Set, which takes the anchor out of force and leaves the
    switch on: through the drag the box held answers -- the hand asked
    x=-0.744 and the box held -1.076 -- and at the release it held -0.206,
    the wish to the digit.
    """
    body = _EDITOR.split('def on_submit_manip_sync(')[1].split('\n    def ')[0]
    assert 'pulling = bool((dragging or released)' in body
    assert 'lit = _thermal_live() or bool(submit_topology_btn.value)' in body
    assert "if pulling and lit and int(state.get('gfn_follow_steps') or 0):" \
        in body
    # Before the raw payload, or it would never be reached.
    assert body.index('if pulling and lit') < body.index('payload = header')


def test_and_leaves_alone_the_two_cases_where_the_hand_is_the_answer():
    """No wall switched on, a placing hand, and a follow that has not answered.

    Each of the three was measured.  With no wall the release is *meant* to
    leave the structure where the hand put it and the next grab carries on
    from there -- refusing it put the second grab 0.29 A behind the first.
    Under a placing hand the answer is laid back onto the cursor, so the box
    has to follow the hand or the drag does nothing at all.  And a follow that
    has answered nothing is the one case where the page's geometry is all
    there is.
    """
    body = _EDITOR.split('def on_submit_manip_sync(')[1].split('\n    def ')[0]
    guard = body.split('if pulling and lit and int(')[1].split('\n\n')[0]
    assert "state.get('gfn_follow_steps') or 0)" in guard
    # And the count belongs to this grab: it is put back to nought where a
    # follow starts, beside the flag the other half reads.
    start = _EDITOR.split("state['gfn_follow'] = True")[1][:200]
    assert "state['gfn_follow_steps'] = 0" in start


@pytest.mark.parametrize('marker', [
    'C-H 1.088 A',
    'C-H 1.307 A',
    '+16.0 of 22.3',
])
def test_the_measurement_is_written_down_where_the_fix_is(marker):
    """A number in a commit message outlives nothing; one here outlives the
    next person who wonders why the anchor is measured twice."""
    assert marker in _body('_keep_the_anchor_current'), marker
