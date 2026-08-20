"""Redo, and the census of every action it can put back.

Undo was proved by going backwards.  Redo replays the same history entries
*forwards*, which is the direction nothing had ever exercised, and that is what
makes it an audit as well as a feature: an entry that fails to carry something
shows up the moment it is walked in both directions.  So every
structure-changing action here is captured, done, undone, redone, and compared
byte for byte -- the coordinates and everything held with them, which is what
:func:`_structure_marks` covers.

What was measured on the twenty-odd actions of that census, walked forwards for
the first time: all of them came back identical.  The entry format Undo was
given carries enough, which is worth knowing rather than assuming.

The semantics are Excel's and Word's, because a surprise costs more here than
anywhere else:

    Undo    puts the state it is leaving on the way forward, then restores the
            previous one.
    Redo    takes the newest state off the way forward, puts the step back on
            the way back, and restores it.
    Anything else clears the way forward.  This is the part people rely on
    without knowing they do -- undo three things, do a fourth, and the three
    are gone rather than waiting to be redone on top of a structure they were
    never part of.  It lives in _remember, which every structure-changing
    action already passes through.

Redo inherits everything Undo had to be taught, because it has the same
exposure: it goes through the same _restore, so it ends what is running, moves
the run number on so an answer in flight is refused, drops the frame the
picture stood on, and puts back the whole of the marks rather than only the
coordinates.
"""

from __future__ import annotations

import copy
import threading

import pytest

_WATER = '3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n'

_BENZENE = """12
benzene
C  1.3970  0.0000  0.0000
C  0.6985  1.2098  0.0000
C -0.6985  1.2098  0.0000
C -1.3970  0.0000  0.0000
C -0.6985 -1.2098  0.0000
C  0.6985 -1.2098  0.0000
H  2.4810  0.0000  0.0000
H  1.2405  2.1486  0.0000
H -1.2405  2.1486  0.0000
H -2.4810  0.0000  0.0000
H -1.2405 -2.1486  0.0000
H  1.2405 -2.1486  0.0000
"""

_TBP = """6
FeCl2Br3
Fe  0.000  0.000  0.000
Cl  0.000  0.000  2.200
Cl  0.000  0.000 -2.200
Br  2.300  0.000  0.000
Br -1.150  1.992  0.000
Br -1.150 -1.992  0.000
"""

_MARKS = ('bond_edits', 'hand_bonds', 'constraints', 'hyb_overrides',
          'poly_applied', 'poly_metal', 'poly_assignment',
          'poly_arrangements', 'poly_arrangement_index')
_EMPTY = {'bond_edits': {}, 'hand_bonds': {}, 'constraints': [],
          'hyb_overrides': {}, 'poly_applied': None, 'poly_metal': None,
          'poly_assignment': None, 'poly_arrangements': [],
          'poly_arrangement_index': 0}


@pytest.fixture
def editor(tmp_path):
    """A Submit tab holding one structure, and a fresh one per structure."""
    pytest.importorskip('ipywidgets')
    from delfin.dashboard import tab_submit
    from delfin.dashboard.context import DashboardContext

    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / 'calc',
                           archive_dir=tmp_path / 'archive',
                           office_dir=tmp_path / 'office')
    ctx.run_js = lambda _script: None

    def load(structure):
        _widget, refs = tab_submit.create_tab(ctx)
        refs['coords_widget'].value = structure
        return refs

    return load


def _rows(refs):
    return [line for line in refs['coords_widget'].value.splitlines()[2:]
            if line.strip()]


def _marks(refs):
    state = refs['editor_state']
    return {key: copy.deepcopy(state.get(key) or _EMPTY[key])
            for key in _MARKS}


def _cmd(refs, verb, payload='', serial=[0]):        # noqa: B006
    serial[0] += 1
    refs['submit_cmd_sync'].value = f'{verb}:{serial[0]}:{payload}'


def _browser_moved(refs, rows, note='DELFIN drag-end'):
    refs['submit_manip_sync'].value = f'{len(rows)}\n{note}\n' + '\n'.join(rows)


def _one_step_first(refs):
    """A recorded step, so Undo has somewhere other than the beginning to go."""
    rows = _rows(refs)
    parts = rows[0].split()
    parts[1] = f'{float(parts[1]) + 0.11:.4f}'
    _cmd(refs, 'grabbed', '0')
    _browser_moved(refs, [' '.join(parts)] + rows[1:])


def _undo(refs):
    refs['submit_manip_undo_btn'].click()


def _redo(refs):
    refs['submit_manip_redo_btn'].click()


def _round_trip(refs, action, what):
    """Do it, take it back, put it back, and say whether it is identical."""
    _one_step_first(refs)
    before, was = refs['coords_widget'].value, _marks(refs)
    action()
    acted, marked = refs['coords_widget'].value, _marks(refs)
    assert (acted != before) or (marked != was), (
        f'{what} did nothing, so nothing was measured')
    _undo(refs)
    assert refs['coords_widget'].value == before, (
        f'Undo did not take back {what}')
    assert _marks(refs) == was, f'Undo left something standing after {what}'
    assert not refs['submit_manip_redo_btn'].disabled, (
        f'there was nothing to redo after undoing {what}')
    _redo(refs)
    assert refs['coords_widget'].value == acted, (
        f'Redo did not put back the coordinates of {what}')
    assert _marks(refs) == marked, (
        f'Redo put back the coordinates of {what} and not what was held '
        f'with them')


# ---------------------------------------------------------------------------
# the census: every action Undo takes back, walked forwards as well
# ---------------------------------------------------------------------------
def test_the_drawing_edits_go_forwards_again(editor):
    """Placing, growing, retyping, deleting, cycling a bond order.  Each of
    them derives more than coordinates -- a placed carbon is typed sp3 and
    that typing is held -- so a Redo that put the geometry back and left the
    typing behind would be the same defect Undo had, read the other way."""
    for structure, verb, payload, what in (
            (_WATER, 'addatom', 'C,3,0,0', 'placing an atom'),
            (_WATER, 'grow', '0,C,1,1,0,0', 'growing an atom'),
            (_WATER, 'setelement', '0,S', 'changing an element'),
            (_BENZENE, 'delatoms', '11', 'deleting an atom'),
            (_BENZENE, 'bondcycle', '0,1', 'cycling a bond order')):
        refs = editor(structure)
        _round_trip(refs, lambda v=verb, p=payload: _cmd(refs, v, p), what)


def test_the_bond_buttons_go_forwards_again(editor):
    """Drawing a bond and cutting one, with Adjust H both ways: the bond edit
    is not a coordinate, so this is the mark half of an entry being walked
    forwards."""
    for adjust in (True, False):
        refs = editor(_BENZENE)
        refs['submit_adjust_h_btn'].value = adjust

        def draw():
            refs['submit_pick_sync'].value = '0,3'
            refs['submit_bond_btn'].click()

        _round_trip(refs, draw, f'drawing a bond with Adjust H {adjust}')

    refs = editor(_BENZENE)

    def cut():
        refs['submit_pick_sync'].value = '0,1'
        refs['submit_unbond_btn'].click()

    _round_trip(refs, cut, 'cutting a bond')


def test_the_three_that_move_atoms_in_the_browser_go_forwards_again(editor):
    """A drag, setting an internal coordinate, exchanging two ligands.  These
    arrive as a geometry the browser made, so the entry is all there is."""
    refs = editor(_WATER)

    def drag():
        rows = _rows(refs)
        parts = rows[0].split()
        parts[2] = f'{float(parts[2]) + 0.31:.4f}'
        _cmd(refs, 'grabbed', '0')
        _browser_moved(refs, [' '.join(parts)] + rows[1:])

    _round_trip(refs, drag, 'a drag')

    refs = editor(_WATER)

    def turn():
        refs['submit_pick_sync'].value = '0,1'
        refs['submit_internal_btn'].value = True
        refs['submit_internal_value'].value = 1.40
        _browser_moved(refs, ['O 0.00 0.00 0.00', 'H 1.40 0.0 0.0',
                              'H -0.24 0.93 0.0'])

    _round_trip(refs, turn, 'setting an internal coordinate')

    refs = editor(_TBP)

    def swap():
        refs['editor_state']['poly_metal'] = 0
        refs['submit_pick_sync'].value = '1,3'
        refs['submit_swap_btn'].click()
        _browser_moved(refs, ['Fe  0.000  0.000  0.000',
                              'Cl  2.300  0.000  0.000',
                              'Cl  0.000  0.000 -2.200',
                              'Br  0.000  0.000  2.200',
                              'Br -1.150  1.992  0.000',
                              'Br -1.150 -1.992  0.000'])

    _round_trip(refs, swap, 'exchanging two ligands')


def test_the_ones_that_move_nothing_themselves_go_forwards_again(editor):
    """A polyhedron, turning it, typing an atom, holding a value, letting go
    of one, hardening a pull into a fix.  None of them is a coordinate and
    every one of them decides where the field takes the structure next, which
    is why they are entries at all."""
    refs = editor(_TBP)
    refs['submit_pick_sync'].value = '0'
    _round_trip(refs, lambda: setattr(refs['submit_poly_dd'], 'value', 'tbp'),
                'putting a polyhedron on')

    refs = editor(_TBP)
    refs['submit_pick_sync'].value = '0'
    refs['submit_poly_dd'].value = 'tbp'
    _round_trip(refs, lambda: refs['submit_poly_turn_btn'].click(),
                'turning the polyhedron')

    refs = editor(_BENZENE)

    def type_it():
        refs['submit_pick_sync'].value = '0'
        refs['submit_hyb_dd'].value = 'sp3'

    _round_trip(refs, type_it, 'typing an atom sp3')

    refs = editor(_WATER)

    def hold():
        refs['submit_pick_sync'].value = '0,1'
        refs['submit_internal_value'].value = 1.10
        refs['submit_hold_btn'].click()

    _round_trip(refs, hold, 'holding a value')

    refs = editor(_WATER)
    refs['submit_pick_sync'].value = '0,1'
    refs['submit_internal_value'].value = 1.10
    refs['submit_hold_btn'].click()

    def release():
        refs['submit_constraint_dd'].value = 'c0'
        refs['submit_constraint_del'].click()

    _round_trip(refs, release, 'letting go of a held value')

    refs = editor(_WATER)
    refs['submit_pick_sync'].value = '0,1'
    refs['submit_internal_value'].value = 1.10
    refs['submit_hold_btn'].click()

    def harden():
        refs['submit_constraint_dd'].value = 'c0'
        refs['submit_hold_mode'].value = 'fix'

    _round_trip(refs, harden, 'holding a value exactly')


def test_reset_goes_forwards_again(editor):
    """Reset is the largest step there is -- everything held on the structure
    goes with it -- so it is the one where an incomplete entry would show
    hardest in either direction."""
    refs = editor(_WATER)
    _cmd(refs, 'addatom', 'C,3,0,0')
    refs['submit_pick_sync'].value = '0,1'
    refs['submit_internal_value'].value = 1.10
    refs['submit_hold_btn'].click()
    built, held = refs['coords_widget'].value, _marks(refs)

    refs['submit_reset_btn'].click()
    reset_to, reset_marks = refs['coords_widget'].value, _marks(refs)
    assert len(_rows(refs)) == 3

    _undo(refs)
    assert refs['coords_widget'].value == built
    assert _marks(refs) == held

    _redo(refs)
    assert refs['coords_widget'].value == reset_to, (
        'Redo did not put the Reset back')
    assert _marks(refs) == reset_marks, (
        'Redo put the geometry back and left what Reset cleared standing')


# ---------------------------------------------------------------------------
# the shape of the pair
# ---------------------------------------------------------------------------
def test_a_run_of_undos_and_redos_returns_to_where_it_started(editor):
    """Three actions, three presses back, three presses forward, and the
    structure and the history are what they were.  The history depth matters
    as much as the coordinates: a Redo that appends something that merely
    looks like the entry Undo took off leaves the next Undo returning what the
    last one restored, which reads as a button that has stopped working."""
    refs = editor(_WATER)
    state = refs['editor_state']
    _one_step_first(refs)
    marks = []
    for element in ('C', 'N', 'O'):
        _cmd(refs, 'addatom', f'{element},3,0,0')
        marks.append(refs['coords_widget'].value)
    depth = len(state['history'])

    for _ in range(3):
        _undo(refs)
    for _ in range(3):
        _redo(refs)

    assert refs['coords_widget'].value == marks[-1], (
        'three back and three forward is not where it started')
    assert len(state['history']) == depth, (
        f'the history is {len(state["history"])} deep and was {depth}')
    assert not state['structure_undo'], 'there is still a way forward'

    # And the way back still works from there, one action at a time.
    for expected in (marks[1], marks[0]):
        _undo(refs)
        assert refs['coords_widget'].value == expected


def test_a_new_action_makes_what_was_undone_unreachable(editor):
    """The rule everybody relies on and nobody names.  Undo twice, do
    something else, and the two that were taken back are gone -- not waiting
    to be redone on top of a structure they were never part of."""
    refs = editor(_WATER)
    _one_step_first(refs)
    _cmd(refs, 'addatom', 'C,3,0,0')
    _cmd(refs, 'addatom', 'N,0,3,0')
    _undo(refs)
    _undo(refs)
    assert refs['editor_state']['structure_undo'], 'nothing was put aside'

    _cmd(refs, 'addatom', 'S,0,0,3')
    branched = refs['coords_widget'].value
    assert not refs['editor_state']['structure_undo'], (
        'the way forward survived a new action')
    assert refs['submit_manip_redo_btn'].disabled, (
        'the button still offers a way forward that is gone')

    _redo(refs)
    assert refs['coords_widget'].value == branched, (
        'a branch that should be unreachable was walked into')
    assert 'Nothing to redo' in refs['mol_status'].value


def test_the_walk_back_to_the_loaded_structure_can_be_walked_forwards(editor):
    """The one Undo that takes nothing off the stack: the last press lands on
    the structure as it arrived and leaves the first entry standing, because
    that entry is what Reset goes back to.  Redo has to put nothing back
    there either, or the history grows an entry per round trip."""
    refs = editor(_WATER)
    state = refs['editor_state']
    _cmd(refs, 'addatom', 'C,3,0,0')
    built = refs['coords_widget'].value
    depth = len(state['history'])

    _undo(refs)
    assert len(_rows(refs)) == 3
    _redo(refs)
    assert refs['coords_widget'].value == built
    assert len(state['history']) == depth, (
        'a round trip through the beginning grew the history')


def test_redo_is_offered_only_when_there_is_a_way_forward(editor):
    """A control that is lit and does nothing is worse than one that is not
    there.  It is greyed on a fresh structure, lit by an Undo, and greyed
    again once the way forward has been walked to its end."""
    refs = editor(_WATER)
    assert refs['submit_manip_redo_btn'].disabled
    _one_step_first(refs)
    assert refs['submit_manip_redo_btn'].disabled
    _undo(refs)
    assert not refs['submit_manip_redo_btn'].disabled
    _redo(refs)
    assert refs['submit_manip_redo_btn'].disabled


def test_the_keyboard_reaches_redo(editor):
    """Ctrl-Shift-Z and Ctrl-Y both arrive as one verb on the command bridge,
    because the browser keeps no way forward of its own -- its snapshot stack
    is cleared by every re-render -- so a redo is always this history's."""
    from delfin.dashboard import molecule_viewer

    script = molecule_viewer.SUBMIT_MANIP_BOOTSTRAP_JS
    assert "pushCommandToPython(undoneKeys[r], 'redo', 'structure')" in script
    assert "key === 'y' || key === 'Y'" in script

    refs = editor(_WATER)
    _one_step_first(refs)
    was = refs['coords_widget'].value
    _undo(refs)
    _cmd(refs, 'redo', 'structure')
    assert refs['coords_widget'].value == was, (
        'the command bridge does not reach Redo')


# ---------------------------------------------------------------------------
# what Redo had to inherit
# ---------------------------------------------------------------------------
def test_a_running_calculation_cannot_write_over_a_redo(editor):
    """The defect Undo was given _stop_what_is_running for, read forwards.

    xtb is minimising a geometry that Redo has just replaced, and it writes
    its answer into the box a second or two later -- so the press looks as
    though it did nothing.  Measured here with a stand-in engine that blocks
    until the test lets it go, so the press lands squarely in the middle of
    the run.  Redo goes through the same _restore, so it moves the run number
    on and the answer is refused.
    """
    from delfin.dashboard import gfn_optimize as gfn

    refs = editor(_WATER)
    _one_step_first(refs)
    dragged = refs['coords_widget'].value
    _cmd(refs, 'addatom', 'C,3,0,0')
    built = refs['coords_widget'].value

    started, release = threading.Event(), threading.Event()

    def stand_in(_xyz, _method, **_kw):
        started.set()
        release.wait(10)
        return {'ok': True, 'converged': True, 'frames': [], 'energy': -5.0,
                'xyz': '3\nrelaxed\nO 9.0 0.0 0.0\nH 9.96 0.0 0.0\n'
                       'H 8.76 0.93 0.0\n'}

    was = (gfn.optimize_with_gfn, gfn.find_binary, gfn.relax_steps)
    gfn.optimize_with_gfn = stand_in
    gfn.find_binary = lambda *_a, **_k: '/bin/true'
    gfn.relax_steps = lambda *_a, **_k: {'ok': True}
    try:
        refs['submit_ff_dd'].value = 'gfn2'
        refs['submit_optimize_btn'].value = True
        assert started.wait(10), 'the stand-in engine was never called'

        # Pressing Optimise is itself an action, so it clears the way
        # forward -- which is the rule and not an accident.  Two presses
        # back put one there again with the run still blocked in the
        # stand-in, which is the state this is about.
        _undo(refs)
        _undo(refs)
        assert refs['coords_widget'].value == dragged
        _redo(refs)
        assert refs['coords_widget'].value == built, 'Redo did not take'
        assert refs['editor_state']['optimize_run'] is None, (
            'the run still owns the coordinate box')

        release.set()
        thread = refs['editor_state'].get('optimize_thread')
        if thread is not None:
            thread.join(10)
        assert refs['coords_widget'].value == built, (
            'the finished run wrote its answer over what Redo put back')
        assert not refs['submit_optimize_btn'].value, (
            'the switch stayed lit over a run that had been stopped')
    finally:
        (gfn.optimize_with_gfn, gfn.find_binary, gfn.relax_steps) = was


def test_a_minimisation_waiting_to_start_stands_down_when_the_structure_does(
        editor):
    """Letting go with Auto on arms a minimisation a third of a second later,
    and a third of a second is long enough to press Undo in.  Nothing
    cancelled the wait: it woke up and ran on the geometry that had just been
    put back, so the press looked as though it had done nothing.  Every timer
    in the editor belongs to the structure it was started for, and the
    generation counter is how the rest of them already say so.
    """
    from editor_source import EDITOR_SOURCE as source

    for arm in ('def _arm_gfn_restart', 'def _arm_gfn_optimise'):
        body = source.split(arm)[1].split('\n    def ')[0]
        assert "generation = state.get('gfn_generation')" in body, (
            f'{arm} does not note which structure it is waiting for')
        assert "if state.get('gfn_generation') != generation:" in body, (
            f'{arm} runs whatever happened while it waited')
    stop = source.split('def _stop_what_is_running')[1].split('\n    def ')[0]
    assert '_gfn_new_generation()' in stop, (
        'taking a structure back does not move the generation on')
