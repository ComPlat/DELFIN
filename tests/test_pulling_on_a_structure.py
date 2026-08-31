"""A load that deforms, and a load that only moves the molecule.

The whole method rests on one distinction: a molecule in free space has six
directions it can move in for nothing, and a load with a component along any of
them does not deform it -- it accelerates it.  E'(x) = E(x) - sum F.x is then
unbounded below, and a minimiser walks the structure off to infinity without
bending a single bond.

So the load is projected onto what is left, and the arithmetic for that is
worth driving rather than reading: it is six vectors, it is done at every level
because the rotations turn with the molecule, and getting it wrong looks like a
calculation that simply never converges.
"""

from __future__ import annotations

import numpy as np
import pytest

from delfin.dashboard import under_load as load


_WATER = np.array([[0.000, 0.000, 0.000],
                   [0.960, 0.000, 0.000],
                   [-0.240, 0.930, 0.000]])


def test_one_arrow_pulls_the_atom_and_leans_on_the_rest():
    """What a single arrow really asks for.

    Pull one atom of a water north with 10 kcal/mol/A and the load has a net
    push: left as drawn, the molecule travels north for ever and nothing bends.
    Projected, a fifth of it goes -- measured, 0.209 of the norm -- and what
    survives is the same pull with the other two atoms leaning back against it.
    That is the deformation the drawing meant, and it is the only part of it
    that has an answer.
    """
    one = np.zeros((3, 3))
    one[0] = (0.0, 10.0, 0.0)
    assert not np.allclose(one.sum(axis=0), 0.0)      # it pushes the molecule

    kept, share = load.deforming_part(one, _WATER)
    assert share == pytest.approx(0.209, abs=0.01), share

    # No net push left, and none of it turns the molecule either.
    assert np.allclose(kept.sum(axis=0), 0.0, atol=1e-9)
    torque = np.cross(_WATER - _WATER.mean(axis=0), kept).sum(axis=0)
    assert np.allclose(torque, 0.0, atol=1e-9)

    # The atom that was pulled is still pulled the way it was drawn, and the
    # other two now carry the reaction.
    assert kept[0][1] > 0.5 * one[0][1], kept
    assert kept[1][1] < 0 and kept[2][1] < 0, kept


def test_two_arrows_pulling_apart_are_all_deformation():
    """A pair of equal and opposite forces is the classic mechanochemical
    load, and there is nothing rigid in it to remove: it neither pushes the
    molecule anywhere nor spins it, it only stretches what is between."""
    pair = np.zeros((3, 3))
    pair[1] = (10.0, 0.0, 0.0)
    pair[2] = (-10.0, 0.0, 0.0)
    # Placed so the pair's line passes through the centre, or the couple has a
    # torque and part of it is a spin.
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    kept, share = load.deforming_part(pair, coords)
    assert share < 1e-6, share
    assert np.allclose(kept, pair, atol=1e-9)


def test_a_couple_that_spins_is_taken_out():
    """Two forces that cancel as a push but not as a turn.  Nothing moves off,
    and the molecule spins for ever instead -- which is the same failure with a
    different shape."""
    couple = np.zeros((3, 3))
    couple[1] = (0.0, 10.0, 0.0)
    couple[2] = (0.0, -10.0, 0.0)
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    assert np.allclose(couple.sum(axis=0), 0.0)          # no net push
    kept, share = load.deforming_part(couple, coords)
    assert share > 0.9, share                            # all of it was a spin
    # And what survives has neither a push nor a turn left in it.
    assert np.allclose(kept.sum(axis=0), 0.0, atol=1e-9)
    torque = np.cross(coords - coords.mean(axis=0), kept).sum(axis=0)
    assert np.allclose(torque, 0.0, atol=1e-9)


def test_the_six_directions_are_six_and_orthonormal():
    ways = load.rigid_directions(_WATER)
    assert ways.shape == (6, 9)
    assert np.allclose(ways @ ways.T, np.eye(6), atol=1e-9)

    # A diatomic has five: it cannot be turned about its own axis.
    two = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.1]])
    assert load.rigid_directions(two).shape == (5, 6)
    # And a single atom has three.
    assert load.rigid_directions(np.zeros((1, 3))).shape == (3, 3)


def test_a_force_in_the_units_the_editor_draws_in():
    """Everything the editor shows a force in is kcal/mol/A -- the push ramps
    in them and A_BOND_HOLDS is quoted in them -- so a drawn arrow is in them
    too, and this is the one place it becomes atomic units."""
    from delfin.dashboard import climb

    assert load.FORCE_IN_KCAL_PER_ANGSTROM == pytest.approx(
        climb.HARTREE_IN_KCAL / climb.BOHR)
    # Which is 1185.8, and a bond holds against 110 of them -- so the ramp's
    # far end at 240 is comfortably past what a bond survives, and its near
    # end at 8 is a nudge.
    assert load.FORCE_IN_KCAL_PER_ANGSTROM == pytest.approx(1185.8, abs=0.1)


def test_putting_it_back_where_it_started():
    """A minimiser takes finite steps, so a long relaxation drifts a little
    even under a balanced load.  Left alone, the frames of a walk swim across
    the screen and every one of them is a different orientation of the same
    answer."""
    turned = _WATER @ np.array([[0.0, -1.0, 0.0],
                                [1.0, 0.0, 0.0],
                                [0.0, 0.0, 1.0]]) + np.array([3.0, -2.0, 1.0])
    back = load._turned_onto(turned, _WATER)
    assert np.allclose(back, _WATER, atol=1e-9)

    # And the shape itself is never touched: a structure that really has been
    # deformed comes back deformed, only oriented.
    bent = _WATER.copy()
    bent[2] += (0.0, 0.30, 0.0)
    assert not np.allclose(load._turned_onto(bent, _WATER), _WATER, atol=1e-3)


_needs_xtb = pytest.mark.skipif(
    not __import__('delfin.dashboard.climb', fromlist=['x']).have_fast_gradients()
    and __import__('delfin.dashboard.climb',
                   fromlist=['x'])._gfn.find_xtb() is None,
    reason='no xtb to take gradients from')

_ETHANE = """8
ethane
C   0.000000  0.000000  0.765000
C   0.000000  0.000000 -0.765000
H   0.000000  1.019000  1.163000
H  -0.882000 -0.510000  1.163000
H   0.882000 -0.510000  1.163000
H   0.000000 -1.019000 -1.163000
H   0.882000  0.510000 -1.163000
H  -0.882000  0.510000 -1.163000
"""


@_needs_xtb
def test_an_ethane_pulled_apart_says_which_bond_gave_and_when():
    """The whole point, end to end, and it finds a number the rest of this
    codebase measured somewhere else.

    Two arrows on the two carbons, straight apart.  Nobody says which bond to
    watch -- that is the reason for not driving a coordinate -- and the ramp
    reports what gave and between which two loads.

    Measured here under GFN2, ten levels from 8 to 240 kcal/mol/A, 0.7 s and
    113 gradients on this machine:

        force  8.0   11.7   17.0   24.9   36.3   52.9   77.2   112.7
        C-C    1.534 1.540  1.549  1.564  1.586  1.625  1.699  3.299
        E      +0.00 +0.06  +0.19  +0.49  +1.19  +2.91  +7.78  +122.4

    It holds at 77 and is broken by 113, and
    :data:`delfin.dashboard.gfn_optimize.A_BOND_HOLDS` -- 110 kcal/mol/A, put
    there from a different measurement entirely -- is what a bond withstands.
    The two agree without ever having been shown each other.
    """
    from delfin.dashboard import gfn_optimize as gfn

    got = load.walk_under_load(
        _ETHANE,
        [{'atom': 0, 'vector': (0.0, 0.0, 1.0)},
         {'atom': 1, 'vector': (0.0, 0.0, -1.0)}],
        'gfn2', steps=10)

    assert got['ok'], got['status']
    assert len(got['points']) >= 8, got['points']

    # A pair pulling straight apart is pure deformation: it neither moves the
    # molecule nor spins it, so nothing had to be taken out of it.
    assert got['rigid'] < 1e-6, got['rigid']

    # The bond stretches before it goes, and the energy rises with it.
    forces = [one['force'] for one in got['points']]
    assert forces == sorted(forces)
    rising = [one['energy'] for one in got['points'][:6]]
    assert rising == sorted(rising), rising

    gave = got['gave']
    assert gave is not None, 'nothing was reported as having given way'
    assert 'C1-C2' in gave['said'], gave
    assert 'break' in gave['said'].lower(), gave
    # And it gave where a bond gives.
    assert gave['held'] < gfn.A_BOND_HOLDS < gave['broke'], (gave,
                                                             gfn.A_BOND_HOLDS)


@_needs_xtb
def test_the_ramp_stops_once_it_has_its_answer():
    """Past the break the load is only drawing two pieces apart: every further
    level moves them by the trust region and says nothing new."""
    got = load.walk_under_load(
        _ETHANE,
        [{'atom': 0, 'vector': (0.0, 0.0, 1.0)},
         {'atom': 1, 'vector': (0.0, 0.0, -1.0)}],
        'gfn2', steps=20)
    assert got['gave'] is not None
    # One level past the break to confirm it stays broken, and then it stops --
    # well short of the twenty that were offered.
    assert len(got['points']) < 20, len(got['points'])
    broke_at = next(i for i, one in enumerate(got['points'])
                    if one['force'] >= got['gave']['broke'] - 1e-6)
    assert len(got['points']) - broke_at <= 2, (broke_at, len(got['points']))


def _an_editor():
    pytest.importorskip('ipywidgets')
    import pathlib
    import tempfile

    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    water = '3\nwater\nO 0.0 0.0 0.0\nH 0.757 0.586 0.0\nH -0.757 0.586 0.0\n'
    room = pathlib.Path(tempfile.mkdtemp())
    for name in ('calc', 'archive', 'office'):
        (room / name).mkdir()
    sent: list = []
    ctx = DashboardContext(calc_dir=room / 'calc', archive_dir=room / 'archive',
                           office_dir=room / 'office')
    ctx.run_js = sent.append
    state: dict = {}
    part = structure_editor.build(
        ctx, state=state, coords_widget=widgets.Textarea(value=water),
        viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None, get_smiles_charge=lambda *a, **k: None)
    part._replace_mol_output_view(water)
    part.submit_ff_dd.value = 'gfn2'
    return part, state, sent


def test_an_arrow_is_hung_kept_listed_and_taken_off_again():
    """The load outlives the picture it was drawn on, so it is kept by atom
    number and not by anything the viewer owns: optimise, step to another
    frame, come back, and it is still the load that was set."""
    part, state, sent = _an_editor()

    part.submit_load_btn.value = True
    assert any("'load'" in one or '"load"' in one for one in sent), (
        'the page was never put into the mode')

    # The page reporting a drag, the way the gesture does.
    sent.clear()
    part.submit_cmd_sync.value = 'load:1:1,0.0,1.0,0.0'
    part.submit_cmd_sync.value = 'load:2:2,0.0,-1.0,0.0'
    assert state['loads'] == [{'atom': 1, 'vector': (0.0, 1.0, 0.0)},
                              {'atom': 2, 'vector': (0.0, -1.0, 0.0)}]
    assert 'setLoads(' in ''.join(sent), 'the page was never told'
    assert [label for label, _v in part.submit_load_dd.options] == [
        'H2  (1.00)', 'H3  (1.00)']

    # A second arrow on an atom that already has one replaces it rather than
    # adding to it: two forces on one atom are one force, and a user aiming
    # again means the new aim.
    part.submit_cmd_sync.value = 'load:4:1,1.0,0.0,0.0'
    assert len(state['loads']) == 2, state['loads']
    assert dict(state['loads'][-1])['vector'] == (1.0, 0.0, 0.0)

    # A press that went nowhere is a press, not an arrow.
    part.submit_cmd_sync.value = 'load:5:0,0.0,0.0,0.0'
    assert all(one['atom'] != 0 for one in state['loads']), state['loads']

    # And the right button takes one off.
    part.submit_cmd_sync.value = 'unload:6:1'
    assert [one['atom'] for one in state['loads']] == [2]


def test_the_scan_can_pull_instead_of_driving():
    """A load arms nothing -- there is no coordinate, which is what it is for
    -- so the press has to answer to the arrows as well as to the legs, or the
    one mode that needs no arming would be the one with no way to start."""
    part, state, _sent = _an_editor()

    assert [label for label, _v in part.submit_scan_how.options] == [
        'push with a force', 'walk the value', 'pull along the arrows']

    # Nothing armed and nothing drawn: no press.
    part.submit_scan_how.value = 'load'
    part._refresh_scan()
    assert part.submit_scan_run_btn.layout.display == 'none'

    part.submit_cmd_sync.value = 'load:1:1,0.0,1.0,0.0'
    part._refresh_scan()
    assert part.submit_scan_run_btn.layout.display == ''
    assert not part.submit_scan_run_btn.disabled

    # And pressing it without a method that can take a gradient says so
    # rather than starting something that cannot answer.
    part.submit_ff_dd.value = 'uff'
    part._pull_along_the_arrows()
    assert 'needs xtb' in part.mol_status.value, part.mol_status.value


def test_the_arrows_belong_to_the_scan_and_nothing_else_changed():
    """The mode exists inside the one press it is for.

    Offered under exactly the condition the Scan press is offered under, so a
    user who never scans never meets it -- and Select, Manipulate and Dynamik
    Opt behave as they always did.  A load without a method that can take a
    gradient is a load nothing can be done with.
    """
    part, state, _sent = _an_editor()

    # A browser force field has no gradient to pull against, and the arrows go
    # with the Scan press they belong to.
    part.submit_ff_dd.value = 'uff'
    assert part.submit_load_btn.layout.display == 'none'
    assert part.submit_scan_btn.layout.display == 'none'

    part.submit_ff_dd.value = 'gfn2'
    assert part.submit_load_btn.layout.display == ''
    assert part.submit_scan_btn.layout.display == ''

    # And a mode cannot be left standing when the row it lives on has gone.
    part.submit_load_btn.value = True
    part.submit_ff_dd.value = 'uff'
    assert part.submit_load_btn.value is False

    # The other three modes are exactly the four-button mutex they were: one
    # at a time, and pressing one puts the others down.
    part.submit_ff_dd.value = 'gfn2'
    part.submit_manip_btn.value = True
    assert part.submit_load_btn.value is False
    part.submit_load_btn.value = True
    assert part.submit_manip_btn.value is False
    assert part.submit_select_btn.value is False
    assert part.submit_draw_btn.value is False


def test_the_first_arrow_says_what_the_scan_is_for():
    """Read as "the box already says load", the mode could never be reached:
    the box that would say it is only on the row once something is armed, and
    the only thing that could arm it is the arrow."""
    part, state, _sent = _an_editor()

    assert part.submit_scan_how.value == 'push'
    part.submit_cmd_sync.value = 'load:1:1,0.0,1.0,0.0'
    assert part.submit_scan_how.value == 'load'
    assert part.submit_scan_run_btn.layout.display == ''
