"""A step is a step when it moves atoms, not when its number is large.

The editor says so when a structure under the hand answers with two
arrangements rather than one -- a fact about the molecule at that load, which
its business is to name rather than to filter out.  It judged that on the
driven coordinate's own number, with one threshold for an Angstrom of a bond
and a degree of a turn, on the reasoning that neither is near a hundredth in
either unit.

Both halves of that are wrong, and drags with a standing cursor say by how
much.  A *converged* relaxation disagrees with itself by 0.08 to 0.15 degrees
of a torsion -- eight to fifteen times a hundredth -- so a firing on a
converged answer was a firing on that noise.  And the number says nothing
about what the user sees: 0.11 degrees of one torsion is 0.0090 A of movement
and is invisible, while 0.144 degrees of another is 0.0837 A and is plain on
the screen.  A factor of nine for the same-sized number.

So both tests are read as the atom travel they imply, which is the one
currency the perception already weighs a bond, an angle and a torsion in --
see :data:`gfn_optimize._MOVED_ANGSTROM`.  The threshold itself does not move,
because in that currency it lands where it should: measured, the acetate's
alternating C-O steps 0.052 A and is named on 21 of 26 answers, the 0.0837 A
torsion is named, and the 0.0090 A one is not.

Nothing here refuses or clamps anything.  The detector adds a clause to a
status line and touches no geometry, so no drag that went anywhere before goes
anywhere different now.
"""
from __future__ import annotations

import math
import pathlib
import tempfile

import pytest

from delfin.dashboard import gfn_optimize as gfn

#: Four atoms in an L, so the torsion 0-1-2-3 has a lever arm worth measuring:
#: atom 0 stands 1.5 A off the 1-2 axis, and a degree about it therefore
#: carries the atom 1.5 * pi/180 = 0.026 A.
_CHAIN = """4
a torsion with a lever arm of 1.5 A
C   0.00  0.00  0.00
C   0.00  0.00  1.50
C   1.50  0.00  1.50
C   1.50  1.20  1.50
"""


def _an_editor(text=_CHAIN):
    """One structure editor over a coordinate box of its own."""
    pytest.importorskip('ipywidgets')
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    room = pathlib.Path(tempfile.mkdtemp())
    for name in ('calc', 'archive', 'office'):
        (room / name).mkdir()
    ctx = DashboardContext(calc_dir=room / 'calc', archive_dir=room / 'archive',
                           office_dir=room / 'office')
    ctx.run_js = lambda _script: None
    return structure_editor.build(
        ctx, state={}, coords_widget=widgets.Textarea(value=text),
        viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None)


def _entry(kind, atoms):
    return {'kind': kind, 'atoms': list(atoms), 'value': 0.0, 'mode': 'drag'}


# ---------------------------------------------------------------------------
# The currency
# ---------------------------------------------------------------------------

def test_a_distance_is_already_the_travel_it_implies():
    """An Angstrom of a bond is an Angstrom of an atom, and no conversion."""
    said = gfn.travel_between(_CHAIN, _entry('distance', [0, 1]), 1.50, 1.55)
    assert said == pytest.approx(0.05, abs=1e-9)


def test_a_turn_is_its_lever_arm_times_the_angle():
    """Which is what makes a degree comparable with an Angstrom at all."""
    said = gfn.travel_between(_CHAIN, _entry('dihedral', [0, 1, 2, 3]),
                              10.0, 11.0)
    assert said == pytest.approx(1.5 * math.radians(1.0), rel=1e-6)
    # And it scales with the angle, not with anything else.
    twice = gfn.travel_between(_CHAIN, _entry('dihedral', [0, 1, 2, 3]),
                               10.0, 12.0)
    assert twice == pytest.approx(2.0 * said, rel=1e-9)


def test_a_turn_across_the_seam_is_the_short_way_round():
    """350 degrees away is ten degrees away.

    Subtracted plainly, a torsion alternating either side of 180 degrees reads
    as most of a circle -- which is a step every threshold passes, on a
    structure that has not moved.
    """
    said = gfn.travel_between(_CHAIN, _entry('dihedral', [0, 1, 2, 3]),
                              179.0, -179.0)
    assert said == pytest.approx(1.5 * math.radians(2.0), rel=1e-6)


def test_a_coordinate_that_is_not_about_this_structure_says_nothing():
    """Rather than answering about atoms that are not there."""
    assert gfn.travel_between(_CHAIN, _entry('distance', [0, 99]),
                              1.0, 1.1) is None
    assert gfn.travel_between('', _entry('distance', [0, 1]), 1.0, 1.1) is None


# ---------------------------------------------------------------------------
# And what it decides
# ---------------------------------------------------------------------------

def _alternating(part, entry, one, two, answers=6):
    """Feed the detector two values over and over; hand back what it said."""
    said = []
    for n in range(answers):
        said.append(bool(part._two_arrangements(
            one if n % 2 == 0 else two, entry, _CHAIN)))
    return said


def test_the_same_numbers_are_two_arrangements_or_none_by_what_they_move():
    """The whole of the change, in one comparison.

    A flip of 0.05 in the driven coordinate's number: as a bond that is 0.05 A
    of an atom and is plainly two arrangements; as a turn about this lever it
    is 0.0013 A, which is a converged relaxation disagreeing with itself and
    is invisible on the screen.  Judged on the number alone, the two were the
    same answer.
    """
    part = _an_editor()
    stretch = _alternating(part, _entry('distance', [0, 1]), 1.50, 1.55)
    assert any(stretch), 'a 0.05 A flip of a bond is two arrangements'

    part = _an_editor()
    turn = _alternating(part, _entry('dihedral', [0, 1, 2, 3]), 10.00, 10.05)
    assert not any(turn), (
        '0.05 degrees about this lever is 0.0013 A of movement, which is the '
        'relaxation disagreeing with itself')

    # And the same turn made big enough to see is named again, so this is a
    # conversion and not a way of never speaking about torsions.
    part = _an_editor()
    plainly = _alternating(part, _entry('dihedral', [0, 1, 2, 3]), 10.0, 13.0)
    assert any(plainly), '3 degrees about this lever is 0.079 A and is visible'


def test_nothing_is_said_before_three_answers_in_a_row():
    """One ragged answer is not an announcement."""
    part = _an_editor()
    said = _alternating(part, _entry('distance', [0, 1]), 1.50, 1.55,
                        answers=6)
    assert said[:4] == [False, False, False, False], said
    assert said[4] is True, said


def test_a_coordinate_the_geometry_does_not_carry_says_nothing():
    """Rather than reading a step off atoms that are not there."""
    part = _an_editor()
    said = _alternating(part, _entry('distance', [0, 99]), 1.50, 1.55)
    assert not any(said), said
