"""A scan leaves the places along its walk that are worth coming back to.

A scan crosses three structures anyone would want again: where it started, the
highest point it went over, and where it came to.  It already works all three
out -- the free-energy mode takes a Hessian at each of them, and two of them
are handed to the path finder as its ends -- and then dropped them on the
floor.  What was left in the box was the last one, and there was no way back
to the other two short of running the scan again the other way.

They go into the history, in the order they were reached, so **Undo walks back
through them and Redo walks forward again**.  That shape rather than a control
of its own, for three reasons:

  * it is what was asked for -- "durch Undo ... wieder zu springen";
  * there is nothing new to press, and the standing rule on this editor is
    that fewer buttons is better;
  * it comes out right with Redo, which is being built at the same time: the
    walk back and the walk forward are the same list read in two directions.

It does not cost the invariant Undo was given.  Every structure-changing
action is still exactly one entry recorded before it happens; these are that
one action's own stations, put in afterwards from the geometries rather than
from the box, because by the time a scan can name them it is standing at the
far end of its own walk.

Measured, on a butane torsion walked from 180 to 60 degrees with GFN-FF in
twelve steps -- the ordinary case, one hill and a well on the far side of it:

    box after the scan          Scanned                       60 deg
    Undo                        the highest point the walk crossed
    Undo                        where the walk started
    Undo                        the structure from before the scan, byte
                                for byte
    Redo, Redo, Redo            the same three, forwards

The wording never assumes what kind of system is being computed: where it
started, the highest point it crossed, where it came to.
"""

from __future__ import annotations

import shutil
import time

import pytest

_needs_xtb = pytest.mark.skipif(not shutil.which('xtb'),
                                reason='xtb not installed')

#: Butane, anti.  A torsion walked from here crosses one eclipsed maximum and
#: falls into the gauche well, which is a scan with all three landmarks in it
#: and cheap enough under GFN-FF to run in a test.
_BUTANE = """14
butane
C -1.9026 -0.2415  0.0000
C -0.6178  0.5776  0.0000
C  0.6178 -0.5776  0.0000
C  1.9026  0.2415  0.0000
H -1.9741 -0.8776  0.8834
H -1.9741 -0.8776 -0.8834
H -2.7756  0.4128  0.0000
H -0.5892  1.2266  0.8794
H -0.5892  1.2266 -0.8794
H  0.5892 -1.2266  0.8794
H  0.5892 -1.2266 -0.8794
H  1.9741  0.8776  0.8834
H  1.9741  0.8776 -0.8834
H  2.7756 -0.4128  0.0000
"""


@pytest.fixture
def editor(tmp_path):
    pytest.importorskip('ipywidgets')
    from delfin.dashboard import tab_submit
    from delfin.dashboard.context import DashboardContext

    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / 'calc',
                           archive_dir=tmp_path / 'archive',
                           office_dir=tmp_path / 'office')
    ctx.run_js = lambda _script: None
    _widget, refs = tab_submit.create_tab(ctx)
    refs['coords_widget'].value = _BUTANE
    return refs


def _comment(refs):
    return refs['coords_widget'].value.splitlines()[1]


def _run_the_scan(refs, steps=10):
    refs['submit_ff_dd'].value = 'gfnff'
    # A torsion is walked to a value; a push is a force between two atoms and
    # is refused on anything but a distance.
    refs['submit_scan_how'].value = 'hold'
    refs['submit_pick_sync'].value = '0,1,2,3'
    refs['submit_internal_value'].value = 180.0
    # An end rather than a direction, which is the third answer in the box
    # that used to be direction-only with a checkbox beside it.
    refs['submit_scan_way'].value = 'to'
    refs['submit_scan_to'].value = 60.0
    refs['submit_scan_steps'].value = steps
    refs['submit_scan_whole'].value = True
    refs['submit_scan_btn'].click()
    before = refs['coords_widget'].value
    refs['submit_scan_run_btn'].click()
    began = time.time()
    while refs['editor_state'].get('scan_run') and time.time() - began < 240:
        time.sleep(0.05)
    time.sleep(0.3)
    return before


@_needs_xtb
def test_undo_walks_back_through_the_scans_landmarks(editor):
    """Where it came to, the highest point it crossed, where it started, and
    then the structure from before any of it -- one press each, in that order.

    The last of those has to come back byte for byte: a scan is still one
    action, and its stations must not cost the way back past it.
    """
    refs = editor
    before = _run_the_scan(refs)
    assert 'Scanned' in _comment(refs), 'the scan wrote nothing'
    came_to = refs['coords_widget'].value

    refs['submit_manip_undo_btn'].click()
    assert 'highest point' in _comment(refs), (
        f'the first press back reached {_comment(refs)!r}')
    summit = refs['coords_widget'].value

    refs['submit_manip_undo_btn'].click()
    assert 'where the walk started' in _comment(refs), (
        f'the second press back reached {_comment(refs)!r}')
    started = refs['coords_widget'].value

    refs['submit_manip_undo_btn'].click()
    assert refs['coords_widget'].value == before, (
        'the walk back does not reach the structure from before the scan')

    # Three distinct places, or there was nothing to walk through.
    assert len({came_to, summit, started}) == 3


@_needs_xtb
def test_redo_walks_forward_through_them_again(editor):
    """The same list read the other way, which is the whole reason the
    landmarks are entries rather than a control of their own."""
    refs = editor
    _run_the_scan(refs)
    reached = [refs['coords_widget'].value]
    for _ in range(3):
        refs['submit_manip_undo_btn'].click()
        reached.append(refs['coords_widget'].value)

    for expected in reversed(reached[:-1]):
        refs['submit_manip_redo_btn'].click()
        assert refs['coords_widget'].value == expected, (
            f'the walk forward reached {_comment(refs)!r} instead')


@_needs_xtb
def test_the_landmarks_are_the_scans_own_geometries(editor):
    """Not something reconstructed afterwards: where the walk started is the
    first point it relaxed, and it is the same structure the path finder is
    given as one of its two ends."""
    refs = editor
    _run_the_scan(refs)
    ends = refs['editor_state'].get('scan_ends')
    assert ends, 'the scan left no ends'
    began_at = [line.split()[1:] for line in ends[0].splitlines()[2:]
                if line.strip()]

    refs['submit_manip_undo_btn'].click()          # the highest point
    refs['submit_manip_undo_btn'].click()          # where it started
    here = [line.split()[1:] for line in
            refs['coords_widget'].value.splitlines()[2:] if line.strip()]
    assert len(here) == len(began_at)
    for mine, theirs in zip(here, began_at):
        for a, b in zip(mine, theirs):
            assert abs(float(a) - float(b)) < 1e-4, (
                'the landmark is not the geometry the scan walked from')


def test_a_landmark_that_repeats_the_one_before_it_is_not_a_step():
    """A walk that only ever went downhill has its highest point at its first
    step, so those two landmarks are one geometry.  Recorded twice, a press of
    Undo would appear to do nothing at all -- which is the one thing Undo may
    never look like."""
    from editor_source import EDITOR_SOURCE as source

    body = source.split('def _remember_landmark')[1].split('\n    def ')[0]
    assert "if history and history[-1].get('coords') == entry['coords']:" in body
    assert 'return False' in body


def test_the_wording_names_no_chemistry():
    """Every kind of system is computed here, so a scan's landmarks are named
    by where they are on the walk and never by what a particular chemistry
    would call them."""
    from editor_source import EDITOR_SOURCE as source

    said = source.split('def _done(final=walked)')[1].split('\n                schedule_ui_update')[0]
    for word in ('educt', 'reactant', 'product', 'transition state', 'TS'):
        assert word.lower() not in said.lower(), (
            f'the scan names its landmarks after {word}')
    assert 'where it started' in said
    assert 'the highest point it crossed' in said
