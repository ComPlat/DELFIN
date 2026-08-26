"""A drag that shakes with nobody moving the mouse.

Reported from a real session, twice in seven minutes: "ich kann die maus halten
an einer stelle es zappelt" and "immer gleiches atom gezogen zappelt".  The
journal that came with it says what was happening, and it is not what any of
the obvious answers would have been.

**Not the thermal wall.**  123 wall events over the session, every one of them
carrying ``wall: null`` -- nothing was ever refused and nothing was ever put
back.  The row said "Within the budget" throughout.

**Not a replayed trajectory.**  The follow sends the whole path from index 0
on every answer, but the page skips what it has already seen, and the frame
index never went backwards.

**The hand's own force.**  It rises with how far the structure has fallen
behind the cursor, which is what pulling on something is like and is the whole
point of :func:`gfn_optimize.as_pushes`.  With no delay that is a spring; with
one it is an oscillator.  From the recorded drag -- 52 answers, the same atom
throughout, the same coordinate held for every one of them:

    2314.4s   pulling 36.7 of a possible 44
    2315.8s   pulling 87.0 of a possible 87
    2317.3s   pulling 44.2 of a possible 44
    2318.7s   pulling 63.3 of a possible 79
    2321.4s   pulling 12.2 of a possible 44

44 is what the slider was set to.  Every other answer went far above it and
pulled at the ceiling, because the structure had fallen behind; the answer
after that overshot and the force collapsed.  One round of that is one answer,
and an answer on this sixty-eight-atom system in solvent was 1.7 seconds.
"""

from __future__ import annotations

import pathlib
import tempfile

import pytest


_ETHANE = """8
ethane
C   -0.7560000000   0.0000000000   0.0000000000
C    0.7560000000   0.0000000000   0.0000000000
H   -1.1404000000   1.0206000000   0.0000000000
H   -1.1404000000  -0.5103000000   0.8838000000
H   -1.1404000000  -0.5103000000  -0.8838000000
H    1.1404000000   0.5103000000   0.8838000000
H    1.1404000000   0.5103000000  -0.8838000000
H    1.1404000000  -1.0206000000   0.0000000000
"""

#: The force the hand asked for, answer by answer, off the recorded session --
#: the drag that begins at 2297.8 s and runs for 52 answers with one atom held
#: and one coordinate driven.
_AS_RECORDED = [
    44.0, 44.0, 44.0, 44.0, 44.0, 44.0, 50.0, 74.0, 44.0, 87.0, 44.0, 79.0,
    44.0, 44.0, 44.0, 44.0, 44.0, 44.0, 44.0, 44.0, 44.0, 58.0, 44.0, 56.0,
    44.0, 72.0, 44.0, 47.0, 44.0, 46.0, 44.0, 45.0, 44.0, 50.0, 44.0, 47.0,
]


def _an_editor(text=_ETHANE):
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
    state = {}
    part = structure_editor.build(
        ctx, state=state, coords_widget=widgets.Textarea(value=text),
        viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None)
    return part, state


def _push(force, atoms=(0, 1), kind='distance', value=1.9):
    return {'kind': kind, 'atoms': list(atoms), 'mode': 'push',
            'force': force, 'value': value}


def _swing(values):
    """The average step-to-step change, which is what a shake is."""
    return (sum(abs(b - a) for a, b in zip(values, values[1:]))
            / max(1, len(values) - 1))


# ---------------------------------------------------------------------------
# The loop, closed
# ---------------------------------------------------------------------------

def test_a_force_that_doubles_and_halves_is_damped_towards_the_middle():
    part, state = _an_editor()
    state['gfn_follow_run'] = 7

    out = []
    for asked in _AS_RECORDED:
        got = part._steady_hand([_push(asked)])
        out.append(float(got[0]['force']))

    assert _swing(_AS_RECORDED) > 7.0, 'the recording really does jump'
    assert _swing(out) < 0.55 * _swing(_AS_RECORDED), (
        _swing(_AS_RECORDED), _swing(out))
    assert max(out) < max(_AS_RECORDED), 'and it never reaches the spike'
    assert min(out) >= min(_AS_RECORDED) - 1e-9, 'nor undershoots the floor'


def test_the_first_answer_of_a_drag_is_not_damped_at_all():
    """Beginning a drag has to be immediate. It is only carrying one on that
    is damped, because only then is there a previous answer to disagree
    with."""
    part, state = _an_editor()
    state['gfn_follow_run'] = 1
    got = part._steady_hand([_push(80.0)])
    assert got[0]['force'] == pytest.approx(80.0)


def test_a_sustained_pull_still_gets_there():
    """The cost of damping, stated: the hand reaches a force it is held at,
    it just takes a few answers about it."""
    part, state = _an_editor()
    state['gfn_follow_run'] = 2
    part._steady_hand([_push(10.0)])
    reached = [float(part._steady_hand([_push(90.0)])[0]['force'])
               for _ in range(4)]
    assert reached[0] < 60.0, 'not in one answer'
    assert reached[2] > 0.85 * 90.0, ('and in three', reached)


def test_each_coordinate_is_damped_on_its_own():
    """A force smoothed across two different pairs of atoms is a force about
    neither of them, and the coordinate the hand drives does change while a
    fragment is walked past a molecule."""
    part, state = _an_editor()
    state['gfn_follow_run'] = 3
    part._steady_hand([_push(20.0, atoms=(0, 1))])
    got = part._steady_hand([_push(80.0, atoms=(2, 3))])
    assert got[0]['force'] == pytest.approx(80.0), (
        'a pair nobody has driven starts where it is asked to')


def test_a_new_run_starts_the_hand_again():
    """Every press that draws over this structure takes a new run number, and
    a force left over from the drag before it is about a geometry nobody
    has."""
    part, state = _an_editor()
    state['gfn_follow_run'] = 4
    part._steady_hand([_push(10.0)])
    state['gfn_follow_run'] = 5
    got = part._steady_hand([_push(90.0)])
    assert got[0]['force'] == pytest.approx(90.0)


def test_nothing_else_about_the_push_is_touched():
    part, state = _an_editor()
    state['gfn_follow_run'] = 6
    first = _push(30.0, atoms=(4, 5), value=2.31)
    got = part._steady_hand([first])[0]
    assert got['atoms'] == [4, 5]
    assert got['kind'] == 'distance'
    assert got['value'] == pytest.approx(2.31)
    assert got['mode'] == 'push'
    assert first['force'] == 30.0, 'and the entry handed in is not written to'


def test_a_push_with_no_force_is_passed_through_untouched():
    part, state = _an_editor()
    state['gfn_follow_run'] = 8
    held = {'kind': 'distance', 'atoms': [0, 1], 'mode': 'fix', 'value': 1.53}
    assert part._steady_hand([held])[0] is held


def test_the_damping_is_in_the_drag_and_nowhere_else():
    """It belongs to the follow, which is the only place a force is worked out
    again and again against a structure that is still moving. A scan holds
    values, and a push in a scan is ramped on purpose."""
    from editor_source import EDITOR_SOURCE

    follow = EDITOR_SOURCE.split('contacts = _gfn.as_pushes(', 1)[1]
    assert '_steady_hand(' in follow[:400]
    assert EDITOR_SOURCE.count('_steady_hand(contacts)') == 1
