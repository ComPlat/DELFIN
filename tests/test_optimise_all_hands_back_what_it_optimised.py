"""Optimise all writes the geometry it optimised into the box, every method.

The picture and the box are fed by two different paths.  A method that streams
its trajectory -- every GFN one does -- plays the optimiser's own walk into the
viewer, and the hand-over to the box was held back so that redrawing would not
rebuild the viewer under the playback.  Held back with it was the hand-over
itself, and in a tab that keeps one structure the box is the only place there
is: the coordinates, Copy XYZ and SUBMIT JOB all kept what the user pressed on
while the line above them said the frames had been optimised and priced a
geometry nobody could reach.

Measured live on a pulled propane under GFN2: "Optimised 2 of 2 frame(s) ...
E = -10.500828 Eh" over a box 0.4751 A from that geometry and 67.08 kcal/mol
above it.  The same press under UFF, which streams nothing, wrote the box.

Stubbed here rather than run, because what is being checked is which path the
answer takes and not what xtb computes.
"""
from __future__ import annotations

import pathlib
import re
import time

import pytest

from delfin.dashboard.context import DashboardContext

_SOURCE = (pathlib.Path(__file__).resolve().parents[1]
           / 'delfin' / 'dashboard' / 'structure_editor.py').read_text()

#: A body, the way a conversion makes one: atom lines and nothing else.
_PROPANE = '\n'.join([
    'C   -1.267300   -0.235800    0.000000',
    'C    0.000000    0.630000    0.000000',
    'C    1.950000   -0.400000    0.000000',
    'H   -2.160000    0.390000    0.000000',
    'H   -1.300000   -0.880000    0.880000',
])
_METHANOL = '\n'.join([
    'C    0.000000    0.000000    0.000000',
    'O    1.550000    0.000000    0.000000',
    'H   -0.400000    1.030000    0.000000',
])


def _moved(body):
    """The same atoms somewhere else, as an optimiser would hand them back."""
    rows = []
    for line in body.splitlines():
        one = line.split()
        rows.append('%s %.6f %.6f %.6f' % (one[0], float(one[1]) + 0.31,
                                           float(one[2]), float(one[3])))
    return '\n'.join(rows)


@pytest.fixture
def submit_tab(tmp_path):
    """The Submit tab, which is the host that keeps exactly one structure.

    The ORCA Builder takes the whole set through ``offer_structures`` and was
    never affected; this one supplies none, so the box is all it has.
    """
    pytest.importorskip('ipywidgets')
    from delfin.dashboard import tab_submit

    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / 'calc',
                           archive_dir=tmp_path / 'archive',
                           office_dir=tmp_path / 'office')
    ctx.run_js = lambda _script: None
    _widget, refs = tab_submit.create_tab(ctx)
    return refs


def _answers_with(streams):
    """An optimiser that hands back a whole XYZ document, as xtb does."""
    def fake(xyz, method, **kw):
        from delfin.dashboard import gfn_optimize as gfn

        rows = gfn.atom_lines(xyz)
        body = _moved('\n'.join(rows))
        frames = [[0.0] * (3 * len(rows))] * 3
        if streams and kw.get('on_frames'):
            kw['on_frames'](frames)
        return {
            'ok': True,
            # A document, comment line and all -- what xtb really returns.
            'xyz': '%d\n energy: -10.5 gnorm: 0.0002 xtb: 6.7.1\n%s\n'
                   % (len(rows), body),
            'energy': -10.5, 'converged': True, 'status': 'converged in 0.1 s',
            'frames': frames if streams else None,
        }
    return fake


def _press_optimise_all(refs, monkeypatch, method, streams):
    from delfin.dashboard import tab_submit

    monkeypatch.setattr(tab_submit._gfn, 'find_binary',
                        lambda _m=None: '/x/xtb')
    monkeypatch.setattr(tab_submit._gfn, 'optimize_with_gfn',
                        _answers_with(streams))
    state = refs['editor_state']
    refs['submit_ff_dd'].value = method
    # The box first: writing it is a structure the host has not seen, and that
    # is what puts the editor back to its defaults -- including the isomers.
    refs['coords_widget'].value = '5\npropane, pulled\n%s\n' % _PROPANE
    state['isomers'] = [(_PROPANE, 5, 'propane'), (_METHANOL, 3, 'methanol')]
    state['isomer_index'] = 0
    before = refs['coords_widget'].value
    refs['submit_optimize_all_btn'].value = True
    until = time.time() + 30
    while time.time() < until and refs['submit_optimize_all_btn'].value:
        time.sleep(0.02)
    time.sleep(0.2)
    return state, before


def _body_of(text):
    from delfin.dashboard import gfn_optimize as gfn
    return '\n'.join(gfn.atom_lines(text))


def test_a_streamed_optimisation_still_writes_the_box(submit_tab, monkeypatch):
    """The defect, at the method it happened under."""
    state, before = _press_optimise_all(submit_tab, monkeypatch, 'gfn2', True)
    box = submit_tab['coords_widget'].value
    line = re.sub('<[^>]+>', '', submit_tab['mol_status'].value)
    assert 'Optimised' in line, line
    assert box != before, 'the box still holds what the user pressed on'
    assert _body_of(box) == _body_of(state['isomers'][0][0]), box
    # And what Copy XYZ and the editor's own buttons read is the same one.
    held = (state.get('current_xyz_for_copy') or {}).get('content')
    assert held and _body_of(held) == _body_of(box), held


def test_and_a_method_that_streams_nothing_behaves_the_same_way(submit_tab,
                                                                monkeypatch):
    """One button, one meaning.  This half always worked; it is the control."""
    state, before = _press_optimise_all(submit_tab, monkeypatch, 'uff', False)
    box = submit_tab['coords_widget'].value
    assert box != before
    assert _body_of(box) == _body_of(state['isomers'][0][0]), box


def test_stepping_onto_an_optimised_structure_gives_a_valid_file(submit_tab,
                                                                 monkeypatch):
    """The list holds bare bodies; an optimisation wrote documents into it.

    Each use of the list puts its own header on, so the box ended up with a
    header saying 11 above thirteen more lines -- xtb's own count and its
    "energy: ... gnorm: ..." among them.  That is the text Copy XYZ and Submit
    hand on.
    """
    _press_optimise_all(submit_tab, monkeypatch, 'gfn2', True)
    for _step in range(3):
        submit_tab['isomer_next_btn'].click()
        rows = [one for one in submit_tab['coords_widget'].value.splitlines()
                if one.strip()]
        assert len(rows) >= 3, rows
        assert rows[0].split() and rows[0].split()[0].isdigit(), rows[0]
        # The header, its comment, and that many atoms -- nothing else.
        assert int(rows[0].split()[0]) == len(rows) - 2, rows[:3]
        assert (_body_of(submit_tab['coords_widget'].value)
                == '\n'.join(rows[2:]))


def test_only_the_drawing_is_held_back_while_a_trajectory_plays():
    """And only for the structure the playback is standing on.

    Redrawing rebuilds the viewer, which is what used to tear the playback
    down milliseconds after it started.  So the guard belongs on that one call
    and not on the hand-over around it -- which is the whole of this defect.
    """
    body = _SOURCE.split('def _show_isomer_at_index(')[1]
    body = body.split('\n    def ')[0]
    guard = 'if not drawn:\n            _replace_mol_output_view(xyz_data)'
    assert guard in body
    # Everything else runs either way: the box, what Copy reads, and the cache
    # Submit reads.
    for must in ("state['current_xyz_for_copy']",
                 "state['converted_xyz_cache']['xyz']",
                 'coords_widget.value = ('):
        assert must in body, must
    # And the press names the frame whose trajectory was played.
    press = _SOURCE.split('def on_submit_optimize(')[1].split('\n    def ')[0]
    assert '_offer_isomers(results, drawn=0 if played[0] else None)' in press
