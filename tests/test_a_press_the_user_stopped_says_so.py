"""Stop is a result, and the one line the user reads called it an optimisation.

Everything else about a stopped press is honest.  The geometry kept is the
frame the picture stood on, and the comment written into the coordinate box by
the same function reads "stopped at the frame on screen".  The status row said
"Optimised with GFN2-xTB." over it.

That matters because the box is what Copy and Submit read.  Measured on a
stretched C20H42 under GFN2 by the swarm that found this, the press stopped at
frame 11 of 44 leaves a geometry 2.84 A and +400 kcal/mol from the minimum;
stopped at frame 2 it is +2404; and where the page never reported which frame
the picture stood on, the box is byte-identical to the input after twelve
frames had already been drawn -- all three under the same sentence.

It got through because the third kind of unfinished run had nobody to speak
for it.  The split into `failures` and `unfinished` keys the second on the
words "before converging", which is what a cycle limit says and not what a
Stop says, so a stopped run fell into neither list and `done == count` made it
a success.  A press knows it was stopped without reading its own prose --
though only while it is still running: by the time the row is written the run
number has moved on whatever happened, so the answer has to be remembered
rather than asked for again.
"""
from __future__ import annotations

import pathlib
import tempfile
import time

import pytest

from delfin.dashboard import gfn_optimize as gfn

_needs_xtb = pytest.mark.skipif(gfn.find_xtb() is None,
                                reason='no xtb to optimise with')

#: A hexane with every C-C stretched to 1.72 A, so a relaxation has a long way
#: to go and there are frames to stop in the middle of.
#:
#: Two structures, because the two questions want opposite things of them.
#:
#: Letting a press finish has to be quick, so that one is small.  Stopping a
#: press in the middle has to be *impossible to lose*: a small molecule under
#: load can be finished before the Stop arrives, and a test that races the
#: optimiser passes alone and fails in the suite -- which is how this one first
#: went in.  So the stopped case is given a chain long enough that GFN2 cannot
#: be done with it, and it never has to converge: the press is stopped at the
#: first frame and the whole test is over in under a second.
_SMALL = """8
ethane, pulled apart
C  0.000000  0.000000  0.900000
C  0.000000  0.000000 -0.900000
H -0.600000  1.030000  1.300000
H -0.600000 -1.030000  1.300000
H  1.200000  0.000000  1.300000
H  0.600000  1.030000 -1.300000
H  0.600000 -1.030000 -1.300000
H -1.200000  0.000000 -1.300000
"""

_LONG = "60\nan alkane chain, every C-C stretched\n" + "\n".join(
    f"C {1.72 * i:10.5f} {0.35 * (i % 2):10.5f}   0.00000" for i in range(20)
) + "\n" + "\n".join(
    f"H {1.72 * i:10.5f} {0.35 * (i % 2) + 1.09:10.5f}   0.00000"
    for i in range(20)
) + "\n" + "\n".join(
    f"H {1.72 * i:10.5f} {0.35 * (i % 2) - 1.09:10.5f}   0.00000"
    for i in range(20)
) + "\n"


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
    return structure_editor.build(
        ctx, state={}, coords_widget=widgets.Textarea(value=text),
        viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None)


def _pressed(text, stop_at=1, seconds=120):
    """Press Optimise, optionally letting go of it at a given frame."""
    part = _an_editor(text)
    part.submit_ff_dd.value = 'gfn2'
    before = part.coords_widget.value
    part.submit_optimize_btn.value = True
    began = time.time()
    if stop_at is not None:
        while (int(part.state.get('gfn_push_end') or 0) < stop_at
               and time.time() - began < seconds):
            time.sleep(0.005)
        part.submit_optimize_btn.value = False
    # The press puts its own switch back up when it is done, which is the one
    # signal that means "this press is over" whichever way it ended.
    while part.submit_optimize_btn.value and time.time() - began < seconds:
        time.sleep(0.02)
    time.sleep(0.3)
    said = str(part.state.get('gfn_last_status') or part.mol_status.value)
    return said, gfn.largest_shift(before, part.coords_widget.value)


@_needs_xtb
def test_a_finished_press_is_still_reported_as_one():
    """Otherwise the fix below would be a way of never saying "Optimised"."""
    said, moved = _pressed(_SMALL, stop_at=None)
    assert 'Optimised with' in said, said
    assert 'stopped' not in said, said
    assert moved > 0.1, (
        f'the press moved the structure {moved:.4f} A, so it never ran')


@_needs_xtb
def test_a_stopped_press_does_not_call_itself_an_optimisation():
    """The row and the coordinate box used to say opposite things.

    The box holds what the user stopped on, and its own comment says so.  The
    row is what they read.
    """
    said, moved = _pressed(_LONG)
    assert 'stopped' in said, said
    assert 'Optimised with' not in said, said
    # And it really was stopped early: the structure has not arrived anywhere.
    assert moved < 0.1, (
        f'the press had already moved {moved:.4f} A, so it was not stopped '
        f'in the middle of anything')
