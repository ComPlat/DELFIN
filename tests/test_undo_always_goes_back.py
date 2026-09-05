"""Undo takes back the thing that was just done, whatever it was.

Driven through the real Submit tab -- the widgets the browser writes to, the
buttons the user presses -- rather than by reading the source: every action
here captures the coordinate box, runs, presses Undo, and compares byte for
byte. Each action is preceded by one recorded step so that the "nothing left
to undo, here is the structure as it was loaded" fallback cannot make a
missing step look like a working one; that fallback is what hid most of this.

What was measured before the change, on the same twenty-four actions:

    restored            a drag, changing an element, holding a value,
                        retuning one, switching the relaxation on
    took back the wrong Bond, Unbond, Set, Swap, Turn, choosing a polyhedron,
    step                letting go of a held value, changing pull to fix,
                        typing an atom sp2, Reset
    restored the        placing an atom, growing one, deleting atoms,
    geometry only       cycling a bond order -- the typing the edit derived
                        stayed behind
    lost outright       Undo during a running optimisation: the geometry came
                        back and xtb wrote its answer over it seconds later

The rule the fixes follow: an action that changes the structure -- its
coordinates, or anything that decides where the field will take them next --
is one entry in the history, recorded before it happens.
"""

from __future__ import annotations

import copy
import threading

import pytest

_WATER = '3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n'

#: C0 and C3 sit across the ring, 2.795 A apart and not bonded, so a bond
#: drawn between them is a real edit whichever way Adjust H is set.
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

#: A trigonal-bipyramidal iron: two kinds of vertex, so its ligands can be
#: arranged more than one way and Turn has somewhere to step to.
_TBP = """6
FeCl2Br3
Fe  0.000  0.000  0.000
Cl  0.000  0.000  2.200
Cl  0.000  0.000 -2.200
Br  2.300  0.000  0.000
Br -1.150  1.992  0.000
Br -1.150 -1.992  0.000
"""

#: Everything about a structure the coordinates do not say, and that an action
#: can change: what is held, what was bonded by hand, how an atom is typed and
#: which polyhedron its ligands sit on.
_MARKS = ('bond_edits', 'hand_bonds', 'constraints', 'hyb_overrides',
          'poly_applied', 'poly_metal', 'poly_assignment',
          'poly_arrangements', 'poly_arrangement_index')
_EMPTY = {'bond_edits': {}, 'hand_bonds': {}, 'constraints': [],
          'hyb_overrides': {}, 'poly_applied': None, 'poly_metal': None,
          'poly_assignment': None, 'poly_arrangements': [],
          'poly_arrangement_index': 0}


@pytest.fixture
def editor(tmp_path):
    """A Submit tab holding one structure, and a fresh one per structure.

    One tab per case rather than one per test: writing a structure a tab is
    already showing does not fire the observer that loads it, so a second case
    in the same tab would go on working on the first one's molecule -- with
    the first one's history behind it.
    """
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
    """A gesture the browser could not finish on its own."""
    serial[0] += 1
    refs['submit_cmd_sync'].value = f'{verb}:{serial[0]}:{payload}'


def _browser_moved(refs, rows, note='DELFIN drag-end'):
    """The geometry the page sends back after it has moved the atoms."""
    refs['submit_manip_sync'].value = f'{len(rows)}\n{note}\n' + '\n'.join(rows)


def _one_step_first(refs):
    """A recorded step, so Undo has somewhere other than the beginning to go.

    Without this every check passes: with a single entry in the history Undo
    hands back the structure as it was loaded, which in a one-action test is
    the same text the action started from.
    """
    rows = _rows(refs)
    parts = rows[0].split()
    parts[1] = f'{float(parts[1]) + 0.11:.4f}'
    _cmd(refs, 'grabbed', '0')
    _browser_moved(refs, [' '.join(parts)] + rows[1:])


def _undo(refs):
    refs['submit_manip_undo_btn'].click()


def _goes_back(refs, action):
    """Run *action*, press Undo, and say whether everything came back."""
    _one_step_first(refs)
    before, marks = refs['coords_widget'].value, _marks(refs)
    depth = len(refs['editor_state']['history'])
    action()
    changed = (refs['coords_widget'].value != before) or (_marks(refs) != marks)
    recorded = len(refs['editor_state']['history']) > depth
    _undo(refs)
    return {
        'changed': changed,
        'recorded': recorded,
        'coords': refs['coords_widget'].value == before,
        'marks': _marks(refs) == marks,
    }


def _assert_undone(refs, action, what):
    outcome = _goes_back(refs, action)
    assert outcome['changed'], f'{what} did nothing, so nothing was measured'
    assert outcome['recorded'], f'{what} is not a step in the history'
    assert outcome['coords'], f'Undo did not put the coordinates back after {what}'
    assert outcome['marks'], f'Undo did not put everything else back after {what}'


# ---------------------------------------------------------------------------
# the drawing edits, which went through the structural path and were steps
# ---------------------------------------------------------------------------
def test_a_drawn_atom_comes_back_out_with_the_typing_it_brought(editor):
    """The geometry always came back; the hybridisation the edit derived did
    not. Placing a carbon types it sp3 and holds that, so an Undo that put the
    coordinates back and left the typing standing handed back a structure the
    field pulls into the shape of an atom that is no longer there."""
    refs = editor(_WATER)
    _assert_undone(refs, lambda: _cmd(refs, 'addatom', 'C,3,0,0'),
                   'placing an atom')


def test_every_other_drawing_edit_comes_back_too(editor):
    """One path -- growing, retyping, deleting, cycling a bond order -- and
    each of them one entry."""
    for structure, verb, payload, what in (
            (_WATER, 'grow', '0,C,1,1,0,0', 'growing an atom'),
            (_WATER, 'setelement', '0,S', 'changing an element'),
            (_BENZENE, 'delatoms', '11', 'deleting an atom'),
            (_BENZENE, 'bondcycle', '0,1', 'cycling a bond order')):
        refs = editor(structure)
        _assert_undone(refs, lambda v=verb, p=payload: _cmd(refs, v, p), what)


# ---------------------------------------------------------------------------
# the bond buttons, which were a step only when Adjust H was on
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('adjust', [True, False])
def test_a_bond_drawn_by_hand_is_a_step_either_way(editor, adjust):
    """With Adjust H on, the bond goes through the structural path and that
    path recorded the step. With it off nothing did: the bond was drawn, the
    field was told about it, and Undo walked straight past to whatever had
    been done before -- which is worse than a button that says it has nothing
    to take back, because it silently takes back something else."""
    refs = editor(_BENZENE)
    refs['submit_adjust_h_btn'].value = adjust

    def draw():
        refs['submit_pick_sync'].value = '0,3'
        refs['submit_bond_btn'].click()

    _assert_undone(refs, draw, f'drawing a bond with Adjust H {adjust}')


def test_cutting_a_bond_is_a_step(editor):
    """Unbond never took the structural path at all, whatever Adjust H said,
    so cutting a bond was never a step."""
    refs = editor(_BENZENE)

    def cut():
        refs['submit_pick_sync'].value = '0,1'
        refs['submit_unbond_btn'].click()

    _assert_undone(refs, cut, 'cutting a bond')


# ---------------------------------------------------------------------------
# the three that move atoms in the browser
# ---------------------------------------------------------------------------
def test_set_is_a_step(editor):
    """Set turns the selection by hand: the browser moves the atoms and hands
    the geometry back. Nothing recorded the structure it started from, so
    turning a dihedral and pressing Undo took back the step before the turn
    and left the turn standing."""
    refs = editor(_WATER)

    def turn():
        refs['submit_pick_sync'].value = '0,1'
        refs['submit_internal_btn'].value = True
        refs['submit_internal_value'].value = 1.40
        _browser_moved(refs, ['O 0.00 0.00 0.00', 'H 1.40 0.0 0.0',
                              'H -0.24 0.93 0.0'])

    _assert_undone(refs, turn, 'setting an internal coordinate')


def test_a_sweep_of_arrow_keys_is_one_step_and_not_two_hundred(editor):
    """With Set on, every repeat of an arrow key comes through the same
    handler -- half a degree of a dihedral at a time. One entry per press
    fills a two-hundred-entry history in two seconds of holding the key and
    makes Undo creep back through the sweep half a degree at a time, so a run
    of the same gesture is one step, ended by any other action. Which is what
    a run of typing is in a text editor."""
    refs = editor(_WATER)
    _one_step_first(refs)
    depth = len(refs['editor_state']['history'])

    refs['submit_pick_sync'].value = '0,1'
    refs['submit_internal_btn'].value = True
    for step in range(20):
        refs['submit_internal_value'].value = 1.00 + 0.01 * step
        _browser_moved(refs, ['O 0.00 0.00 0.00',
                              f'H {1.00 + 0.01 * step:.2f} 0.0 0.0',
                              'H -0.24 0.93 0.0'])

    assert len(refs['editor_state']['history']) - depth == 1, (
        'a sweep is one thing the user did, so it is one step back')


def test_swapping_two_ligands_is_a_step(editor):
    """The line the button prints says "Undo puts them back", and it did not:
    the exchange happens in the browser and arrives as a geometry, and nothing
    recorded the arrangement it replaced."""
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

    _assert_undone(refs, swap, 'exchanging two ligands')


# ---------------------------------------------------------------------------
# and the ones that move nothing themselves, but decide where the field goes
# ---------------------------------------------------------------------------
def test_the_polyhedron_and_turning_it_are_two_steps(editor):
    """Neither moves a coordinate itself. Putting a polyhedron on pulls the
    donors onto its vertices and Turn says which vertex each belongs on, and
    the field walks them there over the next second -- so there was nothing
    for Undo to find, and the ligands went on moving into the arrangement it
    was meant to take back."""
    refs = editor(_TBP)
    refs['submit_pick_sync'].value = '0'
    _assert_undone(refs, lambda: setattr(refs['submit_poly_dd'], 'value', 'tbp'),
                   'putting a polyhedron on')

    refs = editor(_TBP)
    refs['submit_pick_sync'].value = '0'
    refs['submit_poly_dd'].value = 'tbp'
    _assert_undone(refs, lambda: refs['submit_poly_turn_btn'].click(),
                   'turning the polyhedron')
    assert refs['editor_state']['poly_applied'] == 'tbp', (
        'one press is one step: the polyhedron itself was not taken back')


def test_typing_an_atom_is_a_step(editor):
    """A type is what the field builds its angles from, so typing a carbon sp2
    flattens the centre around it as surely as pulling on it would."""
    refs = editor(_BENZENE)

    def type_it():
        refs['submit_pick_sync'].value = '0'
        refs['submit_hyb_dd'].value = 'sp3'

    _assert_undone(refs, type_it, 'typing an atom sp3')


def test_letting_go_of_a_held_value_is_a_step(editor):
    """Setting one was a step and releasing it was not, so Undo could take
    back the holding and never the letting go."""
    refs = editor(_WATER)
    refs['submit_pick_sync'].value = '0,1'
    refs['submit_internal_value'].value = 1.10
    refs['submit_hold_btn'].click()

    def release():
        refs['submit_constraint_dd'].value = 'c0'
        refs['submit_constraint_del'].click()

    _assert_undone(refs, release, 'letting go of a held value')


def test_changing_a_hold_from_pull_to_fix_is_a_step(editor):
    """A pull settles at a compromise and a fix is met exactly, so this moves
    the structure -- and it armed the relaxation that moves it without leaving
    anywhere to come back to."""
    refs = editor(_WATER)
    refs['submit_pick_sync'].value = '0,1'
    refs['submit_internal_value'].value = 1.10
    refs['submit_hold_btn'].click()

    def harden():
        refs['submit_constraint_dd'].value = 'c0'
        refs['submit_hold_mode'].value = 'fix'

    _assert_undone(refs, harden, 'holding a value exactly')


# ---------------------------------------------------------------------------
# Reset
# ---------------------------------------------------------------------------
def test_reset_can_be_taken_back(editor):
    """It was the one action in the editor with no way back at all. An hour of
    work went at a press, and the history was cut to its first entry -- the
    structure Reset had just gone to -- so Undo landed on it and reported that
    there was nothing more to take back."""
    refs = editor(_WATER)
    _cmd(refs, 'addatom', 'C,3,0,0')
    built = refs['coords_widget'].value
    refs['submit_pick_sync'].value = '0,1'
    refs['submit_internal_value'].value = 1.10
    refs['submit_hold_btn'].click()
    held = list(refs['editor_state']['constraints'])

    refs['submit_reset_btn'].click()
    assert len(_rows(refs)) == 3, 'Reset did not go back to what was loaded'
    assert not refs['editor_state']['constraints']

    _undo(refs)
    assert refs['coords_widget'].value == built, 'the structure did not come back'
    assert refs['editor_state']['constraints'] == held, (
        'what was held at the moment of the Reset came back with it')


def test_five_actions_are_five_presses_back_in_order(editor):
    """The whole point, end to end: a drag, a drawn atom, a held value, a
    bond, a Reset -- and five presses of Undo walk back through them one at a
    time and in the order they were done. Before, the same five presses were
    a Reset that could not be taken back followed by four presses that undid
    the drag, the atom and then nothing at all."""
    refs = editor(_WATER)
    state = refs['editor_state']
    loaded = refs['coords_widget'].value

    _one_step_first(refs)
    dragged = refs['coords_widget'].value
    _cmd(refs, 'addatom', 'C,3,0,0')
    built = refs['coords_widget'].value
    refs['submit_pick_sync'].value = '0,1'
    refs['submit_internal_value'].value = 1.10
    refs['submit_hold_btn'].click()
    refs['submit_adjust_h_btn'].value = False
    refs['submit_pick_sync'].value = '0,4'
    refs['submit_bond_btn'].click()
    bonds = len(state['bond_edits'])
    refs['submit_reset_btn'].click()
    assert refs['coords_widget'].value == loaded

    _undo(refs)                          # the reset
    assert refs['coords_widget'].value == built
    assert len(state['bond_edits']) == bonds and len(state['constraints']) == 1

    _undo(refs)                          # the bond
    assert len(state['bond_edits']) == bonds - 1

    _undo(refs)                          # the held value
    assert not state['constraints']
    assert refs['coords_widget'].value == built

    _undo(refs)                          # the drawn atom
    assert refs['coords_widget'].value == dragged

    _undo(refs)                          # the drag
    assert refs['coords_widget'].value == loaded

    _undo(refs)                          # and it stays at the beginning
    assert refs['coords_widget'].value == loaded
    assert 'as it was loaded' in refs['mol_status'].value


def test_reset_still_reaches_the_beginning_afterwards(editor):
    """Reset being undoable must not cost it what it is for: from anywhere in
    a session it goes to the structure as it arrived, and the entry that is
    never dropped is what it goes to."""
    refs = editor(_WATER)
    loaded = refs['coords_widget'].value
    for element in ('C', 'N', 'O'):
        _cmd(refs, 'addatom', f'{element},3,0,0')
    refs['submit_reset_btn'].click()
    assert refs['coords_widget'].value.splitlines()[0] == '3'
    refs['submit_reset_btn'].click()
    assert [row.split()[0] for row in _rows(refs)] == [
        row.split()[0] for row in loaded.splitlines()[2:] if row.strip()]


def test_stepping_to_another_isomer_starts_its_own_history(editor):
    """A step between isomers is quiet -- it draws the structure itself rather
    than letting the write do it -- so the host's "this is a structure I have
    not seen" path never ran and the entries from the isomer before it stayed.
    Undo then put the previous isomer's geometry over this one, and the two
    need not even have the same number of atoms."""
    refs = editor(_WATER)
    assert refs['editor_state']['history'][0]['coords'] == _WATER

    refs['editor_state']['isomers'] = [
        ('O 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n', 3, 'one'),
        ('\n'.join(_BENZENE.splitlines()[2:]) + '\n', 12, 'two'),
    ]
    refs['editor_state']['isomer_index'] = 0
    refs['isomer_next_btn'].click()

    assert len(_rows(refs)) == 12, 'the second isomer is not on screen'
    named = [entry['what'] for entry in refs['editor_state']['history']]
    assert named == ['the structure as it was loaded'], (
        'the history of the isomer before it came along')
    here = refs['coords_widget'].value
    _undo(refs)
    assert refs['coords_widget'].value == here, (
        'Undo reached back into a different molecule')


# ---------------------------------------------------------------------------
# and the one that was not a missing entry at all
# ---------------------------------------------------------------------------
def test_a_running_optimisation_cannot_write_over_an_undo(editor):
    """xtb was minimising the geometry from before the press and wrote its
    answer into the box a few seconds later, over the structure Undo had just
    put back. From outside that is a button that does nothing: the coordinates
    come back and then quietly go again, and the only clue is that what is
    left is an optimisation's geometry rather than the one that was there.

    Measured here with a stand-in engine that blocks until the test lets it
    go, so the press lands squarely in the middle of the run. The run number
    is what settles it -- it is what tells a superseded calculation from the
    current one, and Undo moves it on.
    """
    from delfin.dashboard import gfn_optimize as gfn

    refs = editor(_WATER)
    _one_step_first(refs)
    started, release = threading.Event(), threading.Event()

    def stand_in(_xyz, _method, **_kw):
        started.set()
        release.wait(10)
        return {'ok': True, 'converged': True, 'frames': [],
                'energy': -5.0,
                'xyz': '3\nrelaxed\nO 9.0 0.0 0.0\nH 9.96 0.0 0.0\n'
                       'H 8.76 0.93 0.0\n'}

    was = (gfn.optimize_with_gfn, gfn.find_binary, gfn.relax_steps)
    gfn.optimize_with_gfn = stand_in
    gfn.find_binary = lambda *_a, **_k: '/bin/true'
    gfn.relax_steps = lambda *_a, **_k: {'ok': True}
    try:
        refs['submit_ff_dd'].value = 'gfn2'
        before = refs['coords_widget'].value
        refs['submit_optimize_btn'].value = True
        assert started.wait(10), 'the stand-in engine was never called'

        _undo(refs)
        assert refs['coords_widget'].value == before, 'Undo did not take'
        assert refs['editor_state']['optimize_run'] is None, (
            'the run still owns the coordinate box')

        release.set()
        thread = refs['editor_state'].get('optimize_thread')
        if thread is not None:
            thread.join(10)
        assert refs['coords_widget'].value == before, (
            'the finished run wrote its answer over the structure Undo restored')
        assert not refs['submit_optimize_btn'].value, (
            'the switch stayed lit over a run that had been stopped')
    finally:
        (gfn.optimize_with_gfn, gfn.find_binary, gfn.relax_steps) = was


def test_a_superseded_run_may_not_write_the_coordinate_box(editor):
    """The rule underneath the test above, on its own: a geometry computed
    under one run number is refused once the number has moved on. Every writer
    that answers from a thread carries the number it started under."""
    from editor_source import EDITOR_SOURCE as source

    write = source.split('def _write_coords')[1].split('\n    def ')[0]
    assert "if run is not None and int(state.get('gfn_run', 0)) != int(run):" in write
    # The three that answer from a thread and would otherwise land late.
    assert "drawn=True, run=run)" in source, 'the settle does not carry its run'
    assert "drawn=played[0], run=run_id)" in source, (
        'the optimisation does not carry its run')
    assert "run=state.get('scan_frame_run'))" in source, (
        'the scan does not carry its run')
