"""A walk of N even steps knows its summit to half a step and no better.

What it reports is the highest point it *sampled*, which is at or below the
real one and moves when somebody changes the number of steps.  Eight steps over
1.4 A is a summit known to 0.09 A, and a barrier that depends on the grid it
was measured on is not a barrier -- it is a property of the grid.

So the walk finishes by locating its own top: a parabola through the three
points that bracket it, one extra relaxation a round.  On a smooth top that
converges quadratically, which is why it is a parabola and not a bisection --
one point a round instead of two, for a better answer.

The push needs none of it.  It has :func:`_across`, which prices the crossing
it fell through in one step, and that is the same idea for the case where the
sampling failed completely rather than merely coarsely.
"""

from __future__ import annotations

import pathlib
import tempfile

import pytest

from editor_source import EDITOR_SOURCE as SOURCE


def _an_editor():
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
        ctx, state={}, coords_widget=widgets.Textarea(value=''),
        viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None)


def test_the_vertex_of_three_points():
    """The arithmetic, driven rather than read.

    Three points at 0, *share* and 1, and the top of the parabola through
    them.  Two answers are refusals and both matter more than the fit: a curve
    that opens upward has no top between the points at all, and one that is
    nearly straight puts its vertex far outside the bracket -- which is what
    the numbers look like on the far side of a discontinuity, where a fitted
    top would be a confident answer about a path that is not there.
    """
    vertex = _an_editor()._the_vertex

    # Symmetric: the top is exactly in the middle, and knowing that is the
    # whole reason the middle point is not already the answer.
    assert vertex(0.5, 0.0, 1.0, 0.0) == pytest.approx(0.5)
    # Leaning: the far end is lower, so the top sits nearer the near one.
    assert vertex(0.5, 0.0, 1.0, -1.0) == pytest.approx(0.4167, abs=1e-3)
    assert vertex(0.5, -1.0, 1.0, 0.0) == pytest.approx(0.5833, abs=1e-3)
    # A valley between two peaks has no top between them.
    assert vertex(0.5, 1.0, 0.0, 1.0) is None
    # A straight line has one nowhere near.
    assert vertex(0.5, 0.0, 0.5, 1.0) is None
    # And any three points of a parabola give its top back, wherever the
    # middle one sits.  A walk that stopped at its next minimum has uneven
    # spacing at the end, so the summit can be bracketed by two points of
    # different reach -- assumed to be a half, that is a top in the wrong
    # place with no sign that anything was assumed.
    def peak(x, top):
        return -(x - top) ** 2

    for top in (0.30, 0.50, 0.72):
        for share in (0.25, 0.50, 0.60):
            got = vertex(share, peak(0.0, top), peak(share, top),
                         peak(1.0, top))
            assert got == pytest.approx(top, abs=1e-9), (top, share, got)


def test_where_the_narrowing_sits_in_the_walk():
    """Order is most of what makes it right.

    Before the walk back, or the return leg retraces points the walk never
    held; the refinement's points lie *between* two driven ones and belong to
    the profile rather than to the path that was driven.

    Before the free energies and the landmarks, because both are taken at the
    summit and the summit is what has just moved -- three Hessians on the
    geometry that was second best is the expensive way of answering the wrong
    question.
    """
    walk = SOURCE.split('def on_submit_scan_run')[1].split('\n    def ')[0]
    assert '_narrow_on_the_top(path, summit)' in walk

    at = walk.index('_narrow_on_the_top(path, summit)')
    assert at < walk.index('walking it back'), 'the return leg would retrace it'
    assert at < walk.index('three Hessians'), 'the free energies moved with it'
    assert at < walk.index('_remember_landmark'), 'the landmarks moved with it'

    body = SOURCE.split('def _narrow_on_the_top')[1].split('\n            def ')[0]
    # It is the walk's, not the push's: a push has _across for the same job.
    assert 'if (not pushing' in walk[:at][-400:] or 'not pushing' in walk[at-400:at]
    # The driven values are untouched -- the walk back reads them.
    assert 'drove' not in body
    # And the frames are not added to the playback: they lie between two of
    # its frames, and a trajectory that steps backwards shows nothing.
    assert 'shown' not in body
    assert "state.get('scan_stop')" in body, 'Stop has to reach it'
