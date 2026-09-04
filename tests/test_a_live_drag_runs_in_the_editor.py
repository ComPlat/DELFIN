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
