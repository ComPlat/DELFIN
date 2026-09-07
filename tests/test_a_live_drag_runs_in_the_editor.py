"""A live drag, end to end through the built editor, not just the physics.

test_a_live_hand_steers_one_path proves :func:`climb.steer`; this proves the
wiring around it runs -- the follow worker takes the live branch, prices the
answer against the budget, writes the moved structure back into the box and
sends it out on the frame channel, with no exception anywhere in between.

It builds a real editor (the way the other follow tests do), configures a live
grab, feeds one drag-follow answer, and waits for the background worker.
"""
from __future__ import annotations

import json
import pathlib
import tempfile
import time

import numpy as np
import pytest

from delfin.dashboard import climb as C
from delfin.dashboard import gfn_optimize as gfn


_needs_xtb = pytest.mark.skipif(
    not C.have_fast_gradients() and gfn.find_xtb() is None,
    reason='no xtb to take gradients from')


_BUTANE = """14
butane, relaxed under GFN2
C  -1.54258393  -0.35536195   0.40104530
C  -0.70744036   0.40824399  -0.61910669
C   0.72604655  -0.11291589  -0.68783670
C   1.56472464   0.65918240  -1.69865454
H  -1.11190925  -0.25935801   1.39569002
H  -1.58190998  -1.41271644   0.14722368
H  -2.55983132   0.02872430   0.42912843
H  -0.69261017   1.46842827  -0.35515790
H  -1.16943788   0.31997548  -1.60535233
H   0.71098686  -1.17069122  -0.96132958
H   1.18541349  -0.03392440   0.30040534
H   1.61247655   1.71252240  -1.42993563
H   1.13116952   0.58052009  -2.69354809
H   2.57900527   0.26797097  -1.73457132
"""


def _an_editor(text):
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


def _coords(xyz):
    return np.array(gfn.coordinates_of(xyz), dtype=float).reshape(-1, 3)


def _live_grab(part):
    """Put the editor in a live grab under GFN2 with Dynamik Opt on."""
    part.submit_ff_dd.value = 'gfn2'
    part.submit_relax_btn.value = True
    # The method's own refresh decides the options; the live hand is among
    # them under xtb, and here it is chosen.
    part.submit_hand_dd.options = [('pull with a force', 'pull'),
                                   ('move the atom', 'move'),
                                   ('live dynamics', 'live')]
    part.submit_hand_dd.value = 'live'


def _wait(state, *, steps=1, timeout=40.0):
    end = time.time() + timeout
    while time.time() < end:
        if (int(state.get('gfn_follow_steps') or 0) >= steps
                and not state.get('gfn_follow_busy')):
            return True
        time.sleep(0.02)
    return False


@_needs_xtb
def test_a_live_answer_moves_the_atom_and_writes_it_back():
    part, state = _an_editor(_BUTANE)
    _live_grab(part)

    P = _coords(_BUTANE)
    wish = P.copy()
    wish[3] += (0.0, 0.0, 0.6)          # the cursor leads the terminal carbon
    rows = [f'{s} {r[0]:.6f} {r[1]:.6f} {r[2]:.6f}'
            for s, r in zip('CCCCHHHHHHHHHH', wish)]
    wish_xyz = f'14\nDELFIN drag-follow held=3\n' + '\n'.join(rows) + '\n'

    part._gfn_follow_step(wish_xyz, [3])
    assert _wait(state), state.get('gfn_last_status')

    # It ran the live engine, not an error path.
    assert int(state['gfn_follow_steps']) >= 1
    said = state.get('gfn_last_status') or ''
    assert 'stopped on an error' not in said, said
    assert 'steers the drag' in said, said

    # The box holds a moved, well-formed structure -- same 14 atoms, and the
    # terminal carbon has followed the hand somewhere new.
    box = part.coords_widget.value
    body = [line for line in box.splitlines()[2:] if line.strip()]
    assert len(body) == 14, box
    Q = _coords(box)
    assert np.linalg.norm(Q[3] - P[3]) > 0.02, 'the atom did not move'
    # And it walked, not jumped: the carried-forward step is bounded.
    assert np.linalg.norm(Q - P, axis=1).max() < 1.0, 'an atom jumped'

    # And it went out on the frame channel for the viewer to draw.
    payload = json.loads(part.submit_gfn_frame.value)
    assert payload.get('follow') == 1, payload
    assert payload.get('frames'), payload


@_needs_xtb
def test_a_live_drag_holds_to_the_thermal_budget():
    """The live hand is a force hand, so the budget prices it: with a low
    ceiling a long reach is held back rather than followed all the way, and the
    line says the budget it is resting against."""
    part, state = _an_editor(_BUTANE)
    _live_grab(part)
    # A budget from this structure, at a temperature low enough that stretching
    # a C-C most of the way apart is not affordable.
    part.submit_thermal_btn.value = True
    part.submit_temperature.value = 150.0
    part._set_thermal_anchor(relax=False)

    P = _coords(_BUTANE)
    said_lines = []
    # Drive the terminal carbon straight out along the C2-C3 bond, answer after
    # answer, each wish a little further than the last.
    for k in range(10):
        here = _coords(part.coords_widget.value)
        axis = here[3] - here[2]
        axis = axis / np.linalg.norm(axis)
        wish = here.copy()
        wish[3] = here[3] + axis * 0.25 * (k + 1)
        rows = [f'{s} {r[0]:.6f} {r[1]:.6f} {r[2]:.6f}'
                for s, r in zip('CCCCHHHHHHHHHH', wish)]
        wish_xyz = f'14\nDELFIN drag-follow held=3\n' + '\n'.join(rows) + '\n'
        state['gfn_follow_steps'] = 0
        part._gfn_follow_step(wish_xyz, [3])
        assert _wait(state), state.get('gfn_last_status')
        said_lines.append(state.get('gfn_last_status') or '')

    # The budget is on screen, priced against the anchor at the set
    # temperature, every answer.
    assert any('kcal/mol' in one and '150' in one for one in said_lines), \
        said_lines[-1]
    # And the bond did not run all the way out: the wall held the terminal
    # carbon back short of the last wish, which was 2.5 A of stretch.
    final = _coords(part.coords_widget.value)
    stretched = np.linalg.norm(final[3] - final[2])
    began = np.linalg.norm(P[3] - P[2])
    assert stretched < began + 2.0, (began, stretched)


# A relaxed cis-2-butene, carbons 0-1-2-3, C=C torsion near cis.
_CIS_BUTENE = """12
cis-2-butene, relaxed under GFN2
C   1.56133168  -0.30153550   0.34941920
C   0.58962547   0.81664449   0.15552991
C  -0.69758682   0.72956008  -0.14104061
C  -1.50476993  -0.50896586  -0.35699709
H   2.38089185  -0.20604371  -0.36277491
H   1.10146901  -1.27583353   0.21959173
H   1.98910506  -0.24943338   1.35057641
H   1.02420044   1.80112840   0.27987105
H  -1.26489964   1.64626528  -0.24753152
H  -1.93691854  -0.49814727  -1.35756982
H  -0.91382776  -1.41217448  -0.24475249
H  -2.32872082  -0.54156452   0.35577813
"""

_SYMS12 = 'CCCCHHHHHHHHHH'[:12]


def _dihedral(P, a, b, c, d):
    import math
    b0, b1, b2 = P[a] - P[b], P[c] - P[b], P[d] - P[c]
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.degrees(math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def _rot_about(P, i, ax, deg):
    import math
    a, b = P[ax[0]], P[ax[1]]
    u = (b - a) / np.linalg.norm(b - a)
    v = P[i] - a
    th = math.radians(deg)
    return (a + v * math.cos(th) + np.cross(u, v) * math.sin(th)
            + u * np.dot(u, v) * (1 - math.cos(th)))


@_needs_xtb
def test_the_drive_hand_turns_a_double_bond_from_cis_towards_trans():
    """The drive hand, end to end through the built editor: pick the four
    atoms of the C=C torsion, drag a terminal carbon, and the torsion is
    forced across its barrier -- what neither the pull (drives the wrong
    coordinate) nor the live hand (takes the softest way) could do.  A strong
    hand crosses 90 degrees, the twisted top; the coordinate is what is
    driven, not the atom, so nothing here is tuned to butene.
    """
    part, state = _an_editor(_CIS_BUTENE)
    part.submit_ff_dd.value = 'gfn2'
    part.submit_relax_btn.value = True
    part.submit_hand_dd.options = [('pull with a force', 'pull'),
                                   ('move the atom', 'move'),
                                   ('live dynamics', 'live'),
                                   ('drive coordinate', 'drive')]
    part.submit_hand_dd.value = 'drive'
    part.submit_pull_slider.value = 3.0          # a strong hand, to cross
    state['picked'] = [0, 1, 2, 3]               # the C=C torsion

    started = _dihedral(_coords(_CIS_BUTENE), 0, 1, 2, 3)
    for _ in range(45):
        P = _coords(part.coords_widget.value)
        wish = P.copy()
        wish[3] = _rot_about(P, 3, (1, 2), 20.0)  # drag C3 along the arc
        rows = [f'{s} {r[0]:.6f} {r[1]:.6f} {r[2]:.6f}'
                for s, r in zip(_SYMS12, wish)]
        wish_xyz = f'12\nDELFIN drag-follow held=3\n' + '\n'.join(rows) + '\n'
        state['gfn_follow_steps'] = 0
        part._gfn_follow_step(wish_xyz, [3])
        assert _wait(state), state.get('gfn_last_status')

    reached = abs(_dihedral(_coords(part.coords_widget.value), 0, 1, 2, 3))
    assert abs(started) < 30.0, started
    assert reached > 100.0, f'the drive did not cross the barrier: {reached:.0f}'
    said = state.get('gfn_last_status') or ''
    assert 'drives the dihedral' in said, said


@_needs_xtb
def test_the_drive_hand_works_under_a_group_drag():
    """The real browser case: the coordinate atoms are selected, and grabbing
    one drags the whole selection as a group -- which would leave the
    coordinate unchanged if every atom moved together.  The page names the one
    atom it actually grabbed (``gfn_grabbed``), and the drive reads its target
    from that atom alone, so the coordinate still moves.  This proves the fix
    that makes the drive hand usable with a selection, not only with a lone
    grabbed atom.
    """
    part, state = _an_editor(_CIS_BUTENE)
    part.submit_ff_dd.value = 'gfn2'
    part.submit_relax_btn.value = True
    part.submit_hand_dd.options = [('pull with a force', 'pull'),
                                   ('move the atom', 'move'),
                                   ('live dynamics', 'live'),
                                   ('drive coordinate', 'drive')]
    part.submit_hand_dd.value = 'drive'
    part.submit_pull_slider.value = 3.0
    state['picked'] = [0, 1, 2, 3]

    for _ in range(20):
        P = _coords(part.coords_widget.value)
        # A group drag: all four selected atoms move together (here rotated as
        # a block about the C=C), which alone changes no internal coordinate.
        moved = P.copy()
        for a in (0, 1, 2, 3):
            moved[a] = _rot_about(P, a, (1, 2), 15.0)
        rows = [f'{s} {r[0]:.6f} {r[1]:.6f} {r[2]:.6f}'
                for s, r in zip(_SYMS12, moved)]
        wish_xyz = f'12\nDELFIN drag-follow held=0,1,2,3\n' + '\n'.join(rows) + '\n'
        state['gfn_grabbed'] = 3          # the page says which one was taken
        state['gfn_follow_steps'] = 0
        part._gfn_follow_step(wish_xyz, [0, 1, 2, 3])
        assert _wait(state), state.get('gfn_last_status')

    reached = abs(_dihedral(_coords(part.coords_widget.value), 0, 1, 2, 3))
    assert reached > 40.0, f'the group drag did not drive the torsion: {reached:.0f}'


@_needs_xtb
def test_the_drive_hand_ramps_the_coordinate_on_the_wheel():
    """The drive hand's real gesture: pick the coordinate, then the mouse
    wheel ramps it.  Each notch advances an accumulating target, and the
    coordinate is driven towards it over its barrier -- which the drag could
    not do, because grabbing a picked atom moves the whole selection as a block
    and leaves the coordinate unchanged.  Here the wheel turns a cis double
    bond across 90 degrees towards trans.
    """
    part, state = _an_editor(_CIS_BUTENE)
    part.submit_ff_dd.value = 'gfn2'
    part.submit_relax_btn.value = True
    part.submit_hand_dd.options = [('pull with a force', 'pull'),
                                   ('move the atom', 'move'),
                                   ('live dynamics', 'live'),
                                   ('drive coordinate', 'drive')]
    part.submit_hand_dd.value = 'drive'
    part.submit_pull_slider.value = 3.0
    state['picked'] = [0, 1, 2, 3]

    started = abs(_dihedral(_coords(_CIS_BUTENE), 0, 1, 2, 3))
    for _ in range(30):
        state['gfn_follow_steps'] = 0
        part._drive_wheel(1)                 # one notch of the wheel
        assert _wait(state), state.get('gfn_last_status')

    reached = abs(_dihedral(_coords(part.coords_widget.value), 0, 1, 2, 3))
    assert started < 30.0, started
    assert reached > 90.0, f'the wheel did not drive across the barrier: {reached:.0f}'
    # The target accumulated notch by notch.
    assert state.get('drive_wheel_target') is not None


@_needs_xtb
def test_the_drive_wheel_needs_a_picked_coordinate():
    """Without two-to-four atoms picked there is no coordinate to ramp, so the
    wheel says so rather than doing something arbitrary."""
    part, state = _an_editor(_CIS_BUTENE)
    part.submit_ff_dd.value = 'gfn2'
    part.submit_relax_btn.value = True
    part.submit_hand_dd.options = [('pull with a force', 'pull'),
                                   ('move the atom', 'move'),
                                   ('live dynamics', 'live'),
                                   ('drive coordinate', 'drive')]
    part.submit_hand_dd.value = 'drive'
    state['picked'] = []                     # nothing picked
    part._drive_wheel(1)
    said = ' '.join(str(one) for one in (state.get('mol_status_lines') or ()))
    assert 'pick 2 atoms' in said.lower(), said
    assert state.get('drive_wheel_target') is None


@_needs_xtb
def test_the_wheel_drives_even_when_the_hand_is_counted_gone():
    """The real-browser bug: the frame player sends `gfnfree` whenever the
    mouse is not on an atom -- which is always, on the wheel -- and that clears
    `gfn_follow`, so `_hand_gone()` reads True.  A steer told to stop when the
    hand is gone then breaks on its first step and takes none: the status says
    "drives the distance" while the coordinate never moves.  The wheel's steer
    must not be stoppable that way.  Here `gfn_follow` is False throughout and
    the wheel still moves the coordinate.
    """
    import time as _t
    part, state = _an_editor(_CIS_BUTENE)
    part.submit_ff_dd.value = 'gfn2'
    part.submit_relax_btn.value = True
    part.submit_hand_dd.options = [('pull with a force', 'pull'),
                                   ('move the atom', 'move'),
                                   ('live dynamics', 'live'),
                                   ('drive coordinate', 'drive')]
    part.submit_hand_dd.value = 'drive'
    part.submit_pull_slider.value = 3.0
    state['picked'] = [0, 1, 2, 3]
    state['thermal_was'] = _CIS_BUTENE
    state['drive_wheel_target'] = _dihedral(_coords(_CIS_BUTENE), 0, 1, 2, 3) + 40.0
    state['gfn_follow'] = False          # as gfnfree leaves it -- hand "gone"

    before = _dihedral(_coords(_CIS_BUTENE), 0, 1, 2, 3)
    driven = {'kind': 'dihedral', 'atoms': [0, 1, 2, 3],
              'target': state['drive_wheel_target']}
    part._live_answer(_CIS_BUTENE, [0, 1, 2, 3], _t.perf_counter(),
                      'gfn2', 'GFN2-xTB', 0, 0, None, None, driven=driven)
    after = _dihedral(_coords(part.coords_widget.value), 0, 1, 2, 3)
    assert abs(after - before) > 3.0, (
        f'the wheel took no steps with the hand counted gone: {before:.0f} -> '
        f'{after:.0f}')


def _methyl_sum(P):
    import math
    def ang(i, j, k):
        u, v = P[i] - P[j], P[k] - P[j]
        c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        return math.degrees(math.acos(max(-1.0, min(1.0, c))))
    return ang(2, 0, 3) + ang(2, 0, 4) + ang(3, 0, 4)


_ETHANE8 = """8
ethane
C 0 0 0.765
C 0 0 -0.765
H 0 1.019 1.163
H -0.882 -0.510 1.163
H 0.882 -0.510 1.163
H 0 -1.019 -1.163
H 0.882 0.510 -1.163
H -0.882 0.510 -1.163
"""


@_needs_xtb
def test_the_rest_reacts_fully_when_the_wheel_stops():
    """Beim Anhalten fertig reagieren: a notch gives the rest only a few steps
    to keep up, so mid-scroll it lags above the relaxed path; when the wheel
    stops, the settle holds the driven coordinate where it was taken and lets
    the rest relax the whole way into the reacted geometry.  Proven by the
    energy: the settle lowers it (the rest found a better configuration) while
    the driven coordinate stays where the wheel left it.
    """
    base = gfn.optimize_with_gfn(_ETHANE8, 'gfn2', optimise=True,
                                 max_steps=200, timeout=120)['xyz']
    part, state = _an_editor(base)
    part.submit_ff_dd.value = 'gfn2'
    part.submit_relax_btn.value = True
    part.submit_hand_dd.options = [('pull with a force', 'pull'),
                                   ('move the atom', 'move'),
                                   ('live dynamics', 'live'),
                                   ('drive coordinate', 'drive')]
    part.submit_hand_dd.value = 'drive'
    state['picked'] = [0, 1]                 # the C-C bond

    for k in range(20):                      # drive it well apart
        state['gfn_follow_steps'] = 0
        part.submit_cmd_sync.value = f'drivewheel:{k}:1'
        assert _wait(state), state.get('gfn_last_status')
    timer = state.get('drive_settle_timer')  # the debounce timer this armed
    if timer is not None:
        timer.cancel()                       # fire it ourselves, deterministically

    def _energy(xyz):
        return gfn.optimize_with_gfn(xyz, 'gfn2', optimise=False,
                                     etemp=1000.0, timeout=60)['energy']

    P0 = _coords(part.coords_widget.value)
    cc0 = float(np.linalg.norm(P0[0] - P0[1]))
    e0 = _energy(part.coords_widget.value)

    state['gfn_follow_steps'] = 0
    part._drive_settle(state.get('drive_settle_serial', 0))
    assert _wait(state), state.get('gfn_last_status')

    P1 = _coords(part.coords_widget.value)
    cc1 = float(np.linalg.norm(P1[0] - P1[1]))
    e1 = _energy(part.coords_widget.value)

    # The rest relaxed -- the settle found a lower-energy configuration...
    assert e1 < e0 - 1e-4, (e0, e1)
    # ...while the driven coordinate was held where the wheel left it.
    assert abs(cc1 - cc0) < 0.2, (cc0, cc1)
    # And the rest actually moved (the hydrogens, not the held pair).
    assert np.linalg.norm(P1[2:] - P0[2:], axis=1).max() > 0.02, 'rest stood still'
