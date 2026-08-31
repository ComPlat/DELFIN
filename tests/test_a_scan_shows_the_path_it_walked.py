"""A finished scan draws the profile it walked, and no longer than that.

A relaxed scan is the editor's instrument for walking a reaction: it drives
one or several internal coordinates and prices every point.  It reported
itself in a sentence and in nothing else -- a grep for "plot" and "matplotlib"
over ``structure_editor.py`` returned nothing -- and a sentence is the wrong
shape for what a person reads off a profile at a glance: where the barrier is,
whether it is one barrier or three, whether the walk ended in a minimum or was
still climbing, how flat the top is.  Every number needed was already in hand
when the verdict was written.

Six real walks were made with xtb while this was built, and each is a shape
the picture has to survive:

* A driven Diels-Alder approach, GFN2, C0-C10 from 3.26 to 1.55 A in twenty
  points and 144 s: a flat approach, +6.1 kcal/mol at 2.27, and a fall to
  -67.9 in the product well.
* Butane's central torsion walked the whole way round, twenty-four points in
  238 s: the anti minimum it starts in, two gauche wells and three eclipsed
  tops, where the sentence names one barrier.
* The same approach pushed with a force instead of driven, where the path is
  read off the structure rather than dictated and the crossing the push falls
  through is priced with points of its own -- so the walk has more points than
  it has steps and they are not evenly spaced.
* A walk stopped by hand after 8 of 20 points, which ends *at* its highest
  point: the top and the end are one point, and drawing both markers put two
  labels on one pixel.
* A walk past what the temperature can pay for, where an open circle marks the
  structure the box was handed and the dashed line says why it was not the one
  the walk ended on.
* A two-leg concerted walk, drawn against the first leg with the other named
  under the picture: there is no single x when two coordinates are walked
  together.

What the pictures showed, and what was changed because of it, is recorded in
the comments of :mod:`delfin.dashboard.scan_profile`.  Two of them are here as
tests because they are properties rather than impressions: a walk that closes
a bond drives its coordinate downwards, and drawn on an ascending axis it
reads backwards -- the product on the left, the structure the walk started
from on the right, against the sentence beside it; and a thermal ceiling 22
kcal/mol above a torsion whose whole profile is 4 kcal/mol flattens the
profile into a line, so the ceiling is drawn where the data leaves room for it
and nowhere else.

The rest is the lifetime.  A picture is a claim about what was measured, and
beside a structure the walk never visited it is the same class of lie as a
control offering something impossible.  It goes when a run starts drawing over
the structure -- a hand dragging, a settle, an optimisation, a climb, a saddle
search, a chain, the next scan -- and when a geometry that is not the scan's
own arrives in the box.  It does not go for the run numbers claimed to *stop*
things: Undo and Reset take one to reset the page's player, and Undo walks
back through the scan's own landmarks, which are the walk the picture is of.
"""

from __future__ import annotations

import base64
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import time

import pytest

from tests.editor_source import EDITOR_SOURCE


_needs_xtb = pytest.mark.skipif(not shutil.which('xtb'),
                                reason='xtb not installed')

#: A walk that closes a bond: it starts far out, crosses, and falls into a
#: well.  The numbers are the real Diels-Alder walk above, rounded.
_CLOSING = [(3.26, 0.0), (3.08, 0.3), (2.90, 1.0), (2.72, 2.1), (2.54, 3.9),
            (2.36, 5.6), (2.27, 6.0), (2.18, -17.8), (2.00, -33.7),
            (1.82, -49.5), (1.64, -61.4), (1.55, -63.8)]

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


def _figure(points=_CLOSING, **drawn):
    pytest.importorskip('matplotlib')
    from delfin.dashboard import scan_profile
    drawn.setdefault('x_label', 'C0-C10 (A)')
    drawn.setdefault('y_label', 'kcal/mol above the start')
    drawn.setdefault('title', 'C0-C10, walked')
    return scan_profile.figure(points, **drawn)


def _an_editor(text=_BUTANE):
    """One structure editor over a coordinate box of its own.

    The real part, driven the way the tab drives it: what is on screen and
    what a hook does to it are facts about the widgets, not about the source.
    """
    pytest.importorskip('ipywidgets')
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    room = pathlib.Path(tempfile.mkdtemp())
    for name in ('calc', 'archive', 'office'):
        (room / name).mkdir()
    ctx = DashboardContext(calc_dir=room / 'calc',
                           archive_dir=room / 'archive',
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


def _draw_one(part, points=_CLOSING):
    """Put a picture on the page the way a finished scan does.

    With a path and legs of its own rather than by running xtb: what is being
    measured here is when the picture goes away, and a two-minute walk would
    say nothing more about that than twelve numbers do.  Both halves, in the
    order the scan uses them -- drawn on the worker, shown in the turn.
    """
    part._show_scan_profile(part._scan_profile_html(
        points, [{'kind': 'distance', 'atoms': [0, 1], 'from': 3.26,
                  'to': 1.55, 'steps': 12}], False,
        began=points[0], kept=None))


def _shown(widget):
    return (widget.layout.display or '') != 'none'


def test_the_points_a_scan_computed_are_drawn_as_points():
    """A scan has a dozen points and a curve through them would claim more.

    The difference is exactly at the top, which is where a reader looks: a
    smooth interpolation through twelve points puts a maximum where nothing
    was computed.  So the marks are the points, and the line between them is
    thin and says only in what order they were walked.
    """
    fig = _figure()
    ax = fig.axes[0]
    drawn = [c for c in ax.collections if len(c.get_offsets())
             == len(_CLOSING)]
    assert drawn, 'the computed points are not on the picture as points'
    joins = [line for line in ax.lines
             if len(line.get_xdata()) == len(_CLOSING)]
    assert joins, 'the order they were walked in is not shown'
    assert joins[0].get_linewidth() <= 1.5


def test_the_walk_reads_in_the_direction_it_was_made():
    """A scan that closes a bond drives its coordinate downwards.

    Drawn on an ascending axis the walk reads backwards: measured on the
    Diels-Alder approach from 3.26 to 1.55 A, the product stood on the left
    and the structure the walk started from on the right, against every
    reaction profile ever drawn and against the sentence beside it.
    """
    closing = _figure(_CLOSING).axes[0]
    left, right = closing.get_xlim()
    assert left > right, ('a walk that drives its coordinate down was\n'
                          'drawn backwards')

    opening = _figure([(x, y) for x, y in reversed(_CLOSING)]).axes[0]
    left, right = opening.get_xlim()
    assert left < right, 'a walk outwards must not be turned round'


def test_the_ceiling_is_drawn_where_the_data_leaves_room_for_it():
    """The thermal ceiling is a line, and only where it is on this picture.

    Forced into the range, a 22 kcal/mol budget over a torsion whose whole
    profile is 4 kcal/mol flattens the profile into a line -- and the profile
    is what the picture is for.  Where the ceiling does not fit, the sentence
    beside the picture is what says the path is open; the axes stay the data's.
    """
    flat = [(0.0, 0.0), (60.0, 3.4), (120.0, 0.9), (180.0, 4.1)]
    near = _figure(flat, ceiling=3.0, ceiling_label='3.0 kcal/mol at 60 K')
    far = _figure(flat, ceiling=22.4, ceiling_label='22.4 kcal/mol at 298 K')

    def ceiling_lines(fig):
        ax = fig.axes[0]
        return [line for line in ax.lines
                if len(set(line.get_ydata())) == 1 and len(line.get_xdata())]

    assert ceiling_lines(near), 'a ceiling inside the picture is not drawn'
    assert not ceiling_lines(far), 'a ceiling nobody can see squashed the walk'
    assert far.axes[0].get_ylim()[1] < 10.0, far.axes[0].get_ylim()


def test_what_a_mark_means_is_never_carried_by_its_colour_alone():
    """The two landmarks are the pair a colour-blind reader cannot separate.

    Red for the highest point and green for where the walk was left is the
    house colouring -- the trajectory plot in the Calculations tab labels its
    peak and its end point in the same two -- so they differ in shape as well,
    they are named in the legend, and their numbers are written on the picture
    in ink rather than in the colour of the mark they belong to.
    """
    fig = _figure(started=(3.26, 0.0), top=(2.27, 6.0), ended=(1.55, -63.8))
    ax = fig.axes[0]
    labels = [text.get_text() for text in ax.get_legend().get_texts()]
    # The three the scan puts into the undo history, so the picture names the
    # same places Undo steps back through -- which is what the line under it
    # tells the reader.
    assert 'where it started' in labels
    assert 'the highest point' in labels
    assert 'where it ended' in labels
    shapes = {tuple(map(tuple, c.get_paths()[0].vertices.tolist()))
              for c in ax.collections if len(c.get_offsets()) == 1}
    assert len(shapes) == 3, 'the landmarks are told apart by colour alone'
    written = [text for text in ax.texts
               if text.get_text() in ('+6.0', '-63.8')]
    assert len(written) == 2, [t.get_text() for t in ax.texts]
    assert {text.get_color() for text in written} == {'#263238'}


def test_a_walk_left_at_its_highest_point_is_marked_once():
    """A walk that was stopped, or that ran out of steps still climbing.

    Measured on a Diels-Alder approach interrupted after 8 of 20 points: the
    triangle and the square were drawn on the same pixel with their two labels
    on top of each other.  There is one point there, and what it is, is that
    the walk was left there.
    """
    climbing = [(3.26, 0.0), (3.08, 0.3), (2.90, 1.0), (2.72, 1.5)]
    ax = _figure(climbing, top=(2.72, 1.5), ended=(2.72, 1.5),
                 ended_label='where you stopped it').axes[0]
    labels = [text.get_text() for text in ax.get_legend().get_texts()]
    assert 'where you stopped it' in labels
    assert 'the highest point' not in labels
    assert len([t for t in ax.texts if t.get_text() == '+1.5']) == 1


def test_the_picture_is_made_without_a_display_and_without_pyplot():
    """It is drawn under Voila, from the scan's own background thread.

    Two things follow, and neither is visible from a figure that was built in
    a process where something else has already imported matplotlib: there is
    no display, so the backend has to be a file backend; and ``pyplot`` keeps
    one global figure registry that threads share, so the object API is used
    instead and no figure is ever handed to it.  Asked of a fresh interpreter
    with no DISPLAY, which is the only way to measure either.
    """
    pytest.importorskip('matplotlib')
    root = str(pathlib.Path(__file__).resolve().parent.parent)
    said = subprocess.run(
        [sys.executable, '-c',
         'import sys;'
         'from delfin.dashboard import scan_profile;'
         'out = scan_profile.profile_html('
         '    [(1.0, 0.0), (1.2, 4.0), (1.4, 1.0)],'
         '    x_label="C0-C1 (A)", y_label="kcal/mol", title="walked");'
         'print(len(out), "matplotlib.pyplot" in sys.modules)'],
        capture_output=True, text=True, timeout=300,
        env={'PATH': '/usr/bin:/bin', 'PYTHONPATH': root, 'HOME': '/tmp'},
        cwd=root)
    assert said.returncode == 0, said.stderr[-2000:]
    size, pyplot = said.stdout.split()
    assert int(size) > 1000, said.stdout
    assert pyplot == 'False', 'the picture went through global pyplot state'


def test_the_picture_is_a_png_the_page_can_show_on_its_own():
    """Carried as a data URI, the way the trajectory plot in the Calculations
    tab is: a strict page and an overlay outside the tab both have to be able
    to draw it with nothing to fetch."""
    pytest.importorskip('matplotlib')
    from delfin.dashboard import scan_profile

    drawn = scan_profile.profile_html(
        _CLOSING, x_label='C0-C10 (A)', y_label='kcal/mol above the start',
        title='C0-C10, walked',
        note='Undo steps back through the marked points.')
    match = re.search(r"src='data:image/png;base64,([A-Za-z0-9+/=]+)'", drawn)
    assert match, drawn[:400]
    png = base64.b64decode(match.group(1))
    assert png[:8] == b'\x89PNG\r\n\x1a\n'
    # It scales with the column it is in: the editor is built in a tab, in the
    # ORCA Builder and in an overlay as wide as the screen.
    assert 'width:100%' in drawn.replace(' ', '')
    assert 'Undo steps back through the marked points.' in drawn


def test_the_profile_takes_the_structure_s_place_rather_than_a_row_of_its_own():
    """One or the other, never both, and nothing until there is a walk.

    It had a row under the viewer, and that was wrong in the one way that
    matters: after a scan the panel was two panels, the structure shrank to
    make room, and a user who wanted to go on manipulating the molecule had a
    graph in the way.  The viewer is a fixed height and the profile is a
    full-width picture; side by side they made the panel twice as tall and the
    structure half as useful.  Both are about the same walk and only one of
    them is wanted at a time, so there is a switch.

    An editor that has never scanned is laid out exactly as it was.
    """
    part, _state = _an_editor()
    assert not _shown(part.submit_scan_plot_btn)
    assert part.submit_scan_plot.value == ''
    assert _shown(part.mol_output), 'the structure is what is on screen'
    # It travels into fullscreen without being told to: it is inside the
    # viewer's own box, which is what fullscreen takes, so the picture that
    # replaces the structure goes wherever the structure goes.
    assert 'submit-scan-plot' in part.submit_scan_plot._dom_classes

    from delfin.dashboard import tab_submit
    source = pathlib.Path(tab_submit.__file__).read_text(encoding='utf-8')
    stack = source.split('mol_viewer_stack = widgets.Box(')[1].split(')')[0]
    assert 'submit_scan_plot' in stack, stack
    assert 'mol_output' in stack

    from delfin.dashboard import tab_orca_builder
    builder = pathlib.Path(
        tab_orca_builder.__file__).read_text(encoding='utf-8')
    other = builder.split('orca_mol_stack = widgets.Box(')[1].split(')')[0]
    assert 'submit_scan_plot' in other, other

    # And the switch is on the toolbar rather than with the picture, because
    # that row is always in view: put with the picture it went below the fold
    # with it, and there was no way back to the structure at all.
    assert part.submit_scan_plot_btn in part.submit_manip_toolbar.children


def test_a_run_over_the_structure_takes_the_picture_away():
    """Every walker that draws over this structure moves the run counter.

    The follow under a hand, the settle, an optimisation, a climb, a saddle
    search, a chain, and the next scan: one line at the counter is the whole
    of "something else is moving the structure now".  Taken away when the run
    starts rather than when it ends, because a profile beside a structure
    being pulled about has already stopped being true.

    The counter also moves for the opposite reason -- Undo, Reset and the
    Stops claim a number the page has never seen so that the player drops
    what it was playing -- and those draw nothing.  They are named apart at
    the counter, and a picture that went away on a press of Undo is what said
    the distinction had to be made.
    """
    for walker in ('follow', 'settle', 'optimise', 'climb', 'scan', 'saddle',
                   'chain'):
        part, state = _an_editor()
        _draw_one(part)
        assert _shown(part.submit_scan_plot_btn)
        assert 'base64' in part.submit_scan_plot.value
        part._note_the_run(int(state.get('gfn_run', 0)) + 1, walker)
        assert not _shown(part.submit_scan_plot_btn), walker
        assert part.submit_scan_plot.value == '', walker
        assert not state.get('scan_plot'), walker

    for stopping in ('press', 'abandoned'):
        part, state = _an_editor()
        _draw_one(part)
        part._note_the_run(int(state.get('gfn_run', 0)) + 1, stopping)
        assert _shown(part.submit_scan_plot_btn), stopping
        assert state.get('scan_plot'), stopping


def test_the_picture_stands_only_while_the_box_holds_the_walk():
    """A geometry that arrives without a run of its own is caught at the box.

    Undo, Reset, a structure someone loads, a drag writing where it ended.
    The scan's own landmarks keep it -- they are the walk, and Undo stepping
    back to the highest point it crossed is a way of reading the picture, not
    of leaving it -- and anything else takes it away.
    """
    part, state = _an_editor()
    box = part.coords_widget
    rows = [line for line in _BUTANE.splitlines()[2:] if line.strip()]
    from delfin.dashboard.structure_editor import xyz_document

    _draw_one(part)
    # The landmark Undo reaches for: the same molecule, the scan's own claim.
    box.value = xyz_document(
        rows, 'Scanned: the highest point the walk crossed')
    assert _shown(part.submit_scan_plot_btn), ('Undo took away the picture of '
                                           'the walk it was stepping through')

    box.value = xyz_document(rows, 'Edited in DELFIN viewer')
    assert not _shown(part.submit_scan_plot_btn)

    # And another molecule wearing the scan's comment is not this walk either.
    _draw_one(part)
    assert _shown(part.submit_scan_plot_btn)
    box.value = xyz_document(['O 0.0 0.0 0.0', 'H 0.76 0.59 0.0',
                              'H -0.76 0.59 0.0'], 'Scanned')
    assert not _shown(part.submit_scan_plot_btn)


def test_undo_walks_the_scans_own_landmarks_without_taking_the_picture():
    """The press the picture points at must not be the press that clears it.

    A scan puts three geometries into the undo history and the profile marks
    all three, so Undo is how a reader gets from the picture to the structure
    at a point on it.  Measured on a real twenty-four-point torsion before
    this was right: the walk finished, the picture stood, and one press of
    Undo took it away over "Scanned: the highest point the walk crossed" --
    the very walk the picture is of.

    The cause was the run counter.  Undo claims a number to stop whatever was
    walking and to give the page's player one it has never seen, and that
    claim looked exactly like an optimisation starting.  It is not: it draws
    nothing.  What arrives in the box decides this case instead.
    """
    part, state = _an_editor()
    from delfin.dashboard.structure_editor import xyz_document
    rows = [line for line in _BUTANE.splitlines()[2:] if line.strip()]

    # The three entries a finished scan leaves, in the order it leaves them.
    part._remember('a scan of C0-C1')
    part._remember_landmark(_BUTANE, 'the scan, on from where it started',
                            'Scanned: where the walk started')
    part._remember_landmark(_BUTANE.replace('-1.9026', '-1.9126'),
                            'the scan, on from the highest point it crossed',
                            'Scanned: the highest point the walk crossed')
    part.coords_widget.value = xyz_document(
        rows, 'Scanned to the next minimum')
    _draw_one(part)
    assert _shown(part.submit_scan_plot_btn)

    part.submit_manip_undo_btn.click()          # to the highest point crossed
    assert part.coords_widget.value.splitlines()[1].startswith('Scanned')
    assert _shown(part.submit_scan_plot_btn), 'Undo took away the picture of the '\
                                          'walk it was stepping through'
    part.submit_manip_undo_btn.click()          # to where the walk started
    assert _shown(part.submit_scan_plot_btn)
    part.submit_manip_undo_btn.click()          # out of the scan altogether
    assert not part.coords_widget.value.splitlines()[1].startswith('Scanned')
    assert not _shown(part.submit_scan_plot_btn)


def test_the_picture_is_drawn_once_off_the_turn_and_shown_after_the_answer():
    """Not a point at a time, not in the turn, and never before the sentence.

    Measured with the interpreter the dashboard runs on: 0.17 s to build the
    figure and encode it, and 40 to 66 kB of PNG.  The walks it was measured
    against cost 6 to 10 s a point under GFN2 -- the twenty-point Diels-Alder
    approach took 144 s, the twenty-four-point torsion 238 s -- so a picture
    per point would cost a few per cent of the time and push a megabyte of
    pictures nobody waits for down the channel the frame stream needs.  And it
    is drawn and shown only after the turn that carries the answer has been
    handed over, because the sentence is the answer and a picture must never
    be able to delay or replace it.

    Measured, and caught by six of the scan tests when it was the other way
    round: matplotlib's first import is 0.3 s on an idle machine and over a
    second on a loaded one, and with the figure built before the answer was
    scheduled that whole delay fell between the walk ending and the verdict
    reaching the row -- the tests found "step 14 of 20" still standing after
    the run had finished.  The figure is still built on the thread that did
    the walking rather than in the interface's turn, for the same reason
    turned round: a second spent inside a callback is a second in which the
    dashboard does not answer.
    """
    body = EDITOR_SOURCE.split('def on_submit_scan_run(')[1]
    body = body.split('\n    def _said_modes(')[0]
    assert body.count('_scan_profile_html(') == 1
    assert body.count('_show_scan_profile,') == 1
    # Nothing of the picture happens before the answer is on its way: not the
    # drawing, and not the showing.
    assert body.index('schedule_ui_update(_done)') < body.index(
        '_scan_profile_html(')
    assert body.index('schedule_ui_update(_done)') < body.index(
        '_show_scan_profile,')
    # And the figure is built outside the turn: it is the argument, worked
    # out here, and only the write is scheduled.
    assert ('schedule_ui_update(\n                    _show_scan_profile,'
            in body)
    # And it cannot take the answer away: the whole of the drawing is wrapped.
    drawing = EDITOR_SOURCE.split('def _scan_profile_html(')[1]
    drawing = drawing.split('\n    def ')[0]
    assert 'except Exception' in drawing


@_needs_xtb
def test_a_real_scan_leaves_its_profile_on_the_page():
    """The whole path, on the real editor: arm a leg, press, get a picture.

    Butane's central C-C under GFN-FF in eight points, which is seconds.  What
    is asserted is what the user sees: a switch that was not there before the
    press is there after it, the structure is still on screen, and pressing
    the switch puts a PNG of the walk in its place.
    """
    part, state = _an_editor()
    part.submit_ff_dd.value = 'gfnff'
    part.submit_scan_how.value = 'hold'
    part.submit_pick_sync.value = '1,2'
    part.submit_internal_value.value = 1.53
    part.submit_scan_way.value = 'to'
    part.submit_scan_to.value = 2.20
    part.submit_scan_steps.value = 8
    part.submit_scan_whole.value = True
    part.submit_scan_add_btn.click()
    assert not _shown(part.submit_scan_plot_btn)

    part.submit_scan_run_btn.click()
    began = time.time()
    while state.get('scan_run') and time.time() - began < 240:
        time.sleep(0.05)
    time.sleep(0.3)

    said = ' '.join(state.get('mol_status_lines') or ())
    assert 'The scan walked' in said, said
    assert _shown(part.submit_scan_plot_btn), said
    # The walk ends with the structure it reached on screen, which is what
    # the user goes on working with. The profile is there for the asking.
    assert _shown(part.mol_output)
    assert not _shown(part.submit_scan_plot)
    assert part.submit_scan_plot_btn.description == 'Show the profile'

    part.submit_scan_plot_btn.value = True
    assert _shown(part.submit_scan_plot)
    assert not _shown(part.mol_output), 'one or the other, never both'
    # And the line that lies on the picture goes with the structure. It costs a
    # molecule nothing -- there are empty corners -- but a profile has none:
    # its bottom edge is the axis, the numbers and the caption, and a line over
    # that hides the part drawn to be read.
    assert not _shown(part.mol_status)
    assert not _shown(part.mol_status_fs)
    assert part.submit_scan_plot_btn.description == 'Back to the structure'

    match = re.search(r"base64,([A-Za-z0-9+/=]+)'",
                      part.submit_scan_plot.value)
    assert match, part.submit_scan_plot.value[:300]
    assert base64.b64decode(match.group(1))[:8] == b'\x89PNG\r\n\x1a\n'


_MEASURE_JS = r"""
(names) => {
  const plot = document.querySelector('.submit-scan-plot');
  const stack = document.querySelector('.delfin-structure-viewer-stack');
  if (!plot || !stack) return {error: 'nothing to measure'};
  const img = plot.querySelector('img');
  const box = plot.getBoundingClientRect();
  const view = stack.getBoundingClientRect();
  const picture = img ? img.getBoundingClientRect() : null;
  const column = plot.parentElement.getBoundingClientRect();
  const viewer = document.querySelector('.submit-mol-output');
  return {
    shown: getComputedStyle(plot).display !== 'none',
    viewerShown: viewer
        ? getComputedStyle(viewer).display !== 'none' : null,
    overflows: picture ? picture.right > column.right + 0.5 : null,
    wider: picture ? picture.width > column.width + 0.5 : null,
    viewer: Math.round(view.height),
    plotHeight: Math.round(box.height),
    inOverlay: !!plot.closest('.delfin-structure-fs-overlay'),
    pageScrollsSideways:
        document.documentElement.scrollWidth > innerWidth + 1,
  };
}
"""


@pytest.mark.parametrize('width', (1920, 1280))
def test_the_profile_lays_out_in_the_structure_s_place_in_a_browser(width):
    """Measured in chromium, because no source says where a box ends up.

    Two things it has to be: in the structure's place rather than under it,
    and inside its column -- the image is 936 px wide as it is drawn and the
    column in the Submit tab is narrower than that at every window width, so
    it is the ``width: 100%`` on the image that keeps the page from scrolling
    sideways.  And in the fullscreen overlay it has to travel with the
    structure it is about.
    """
    pytest.importorskip('ipywidgets')
    playwright = pytest.importorskip(
        'playwright.sync_api', reason='needs a browser to lay the row out')
    pytest.importorskip('matplotlib')
    from delfin.dashboard import scan_profile
    from delfin.dashboard.molecule_viewer import (
        structure_viewer_fullscreen_bootstrap_js,
        structure_viewer_fullscreen_kind_js,
    )
    toolbar_test = pytest.importorskip('test_the_toolbar_stays_on_the_screen')

    stylesheet = toolbar_test._widget_stylesheet()
    tab_widget, exports = toolbar_test._build_tab()
    exports['submit_manip_toolbar'].layout.display = 'flex'
    # The state a finished scan leaves with the switch pressed: the picture
    # on screen and the structure out of the way.
    exports['submit_scan_plot'].value = scan_profile.profile_html(
        _CLOSING, x_label='C0-C10 (A)', y_label='kcal/mol above the start',
        title='C0-C10, walked', top=(2.27, 6.0), ended=(1.55, -63.8))
    exports['submit_scan_plot'].layout.display = ''
    exports['submit_scan_plot'].layout.display = ''
    exports['mol_output'].layout.display = 'none'

    from delfin.dashboard.molecule_viewer import (
        STRUCTURE_VIEWER_FULLSCREEN_CSS,
    )
    document = (
        '<!doctype html><html><head><meta charset="utf-8"><style>'
        + stylesheet + '\n' + STRUCTURE_VIEWER_FULLSCREEN_CSS
        + '\nhtml, body { margin: 0; padding: 0; }\n'
        + '</style></head><body>' + toolbar_test._render(tab_widget)
        + '<script>' + structure_viewer_fullscreen_bootstrap_js() + '\n'
        + structure_viewer_fullscreen_kind_js(
            'submit', 'submit-scope-', ['_submitMolViewerByScope'])
        + '</script></body></html>')

    with playwright.sync_playwright() as engine:
        try:
            browser = engine.chromium.launch()
        except Exception as exc:                     # no browser on this box
            pytest.skip(f'chromium unavailable: {exc}')
        try:
            page = browser.new_page(viewport={'width': width, 'height': 900})
            page.set_content(document)
            page.wait_for_timeout(150)
            embedded = page.evaluate(_MEASURE_JS, None)
            page.evaluate("() => document.querySelector("
                          "'.delfin-structure-fullscreen-btn').click()")
            page.wait_for_timeout(150)
            enlarged = page.evaluate(_MEASURE_JS, None)
            page.close()
        finally:
            browser.close()

    assert 'error' not in embedded, embedded
    assert embedded['shown'], embedded
    assert embedded['viewerShown'] is False, (
        'one or the other, never both')
    assert not embedded['overflows'] and not embedded['wider'], embedded
    assert not embedded['pageScrollsSideways'], embedded
    assert embedded['plotHeight'] > 40, embedded

    assert 'error' not in enlarged, enlarged
    assert enlarged['inOverlay'], 'the profile stayed on the page'


def test_the_line_comes_back_with_the_structure():
    """Nothing is lost by taking it away: what it was saying is the verdict of
    the walk the profile is of, and the profile carries its own caption."""
    part, _state = _an_editor()
    part._show_scan_profile("<img src='x'/>")
    assert _shown(part.mol_status), 'a finished walk still reports'

    part.submit_scan_plot_btn.value = True
    assert not _shown(part.mol_status)
    part.submit_scan_plot_btn.value = False
    assert _shown(part.mol_status), 'and it is back, unchanged'

    # And a picture that is dropped while it is showing must not leave the
    # panel with neither of them in it.
    part.submit_scan_plot_btn.value = True
    part._scan_plot_drop()
    assert _shown(part.mol_output) and _shown(part.mol_status)


def _a_walk_of(part, count=7):
    """A finished walk, left the way a real one leaves itself.

    On its last point.  That matters to every test below it: a walk ends by
    writing the structure it reached into the box, and the slider's whole
    rule is which of the walk's points the box is holding.
    """
    def geo(d):
        return f"2\nstep\nC 0.000 0.000 0.000\nO {d:.3f} 0.000 0.000\n"
    points = [(1.13 + 0.1 * i, i * 3.4, geo(1.13 + 0.1 * i))
              for i in range(count)]
    part.coords_widget.value = points[-1][2]
    part.state['scan_walk'] = {
        'points': points, 'method': 'gfn2', 'charge': 0, 'uhf': 0,
        'solvent': '', 'solvation_model': '',
        'structure': part._structure_fingerprint(part.coords_widget.value)}
    part.state['walk_points_stale'] = False
    part._refresh_the_walk_points()
    return points


def test_the_points_of_a_walk_can_be_stepped_through():
    """A scan keeps every geometry it computed and only the two ends could be
    reached. The steps between them are where the chemistry is -- the bond
    half broken, the point the profile marks as the top -- and they cost
    minutes to compute and nothing to look at."""
    part, _state = _an_editor()
    assert not _shown(part.submit_walk_at), 'nothing walked yet'

    points = _a_walk_of(part)
    assert _shown(part.submit_walk_at)
    assert part.submit_walk_at.max == len(points) - 1

    part.submit_walk_at.value = 3
    said = ' '.join(part.state.get('mol_status_lines') or ())
    assert 'Point 4 of 7' in said, said
    assert '1.43' in said and '+10.2' in said


def test_stepping_through_a_walk_goes_to_the_point():
    """What is on screen is what the presses act on. Anything else is a
    structure nobody can compute anything about.

    Measured from a real session on a 64-atom system: the slider was stepped
    to point 10 near the top of a barrier and Show the shape pressed, and the
    answer was "a minimum, not a transition state"; then point 11, the same
    press, and the same free energy to the second decimal for a different
    structure. Both Hessians had run on what was in the box, which was the
    minimum the scan ended on.
    """
    part, _state = _an_editor()
    points = _a_walk_of(part)

    for n in (1, 4, 2, 6):
        part.submit_walk_at.value = n
        assert (part._geometry_key(part.coords_widget.value)
                == part._geometry_key(points[n][2])), n
        assert part.coords_widget.value.splitlines()[1].lower().startswith(
            'scanned'), 'and it is still the walk the profile is of'


def test_a_run_of_stepping_is_one_press_of_undo():
    """Nothing is lost by going there, and going back is not twelve presses:
    a run of stepping is one step in the history, the way a sweep of an arrow
    key is one rather than two hundred."""
    part, state = _an_editor()
    _a_walk_of(part)
    before = part.coords_widget.value
    depth = len(state.get('history') or ())

    for n in (1, 4, 2, 6):
        part.submit_walk_at.value = n
    assert len(state.get('history') or ()) == depth + 1

    part.on_submit_manip_undo()
    assert part.coords_widget.value == before


def test_the_points_belong_to_the_molecule_they_were_walked_on():
    part, _state = _an_editor()
    _a_walk_of(part)
    assert _shown(part.submit_walk_at)
    part.coords_widget.value = "3\nother\nO 0 0 0\nH 1 0 0\nH 0 1 0\n"
    part._refresh_the_walk_points()
    assert not _shown(part.submit_walk_at)


def test_the_slider_opens_on_the_point_the_box_is_holding():
    """Which is the end, for a walk that has just finished: it leaves the
    structure it reached on screen, and a row reading "point 1" beside a
    viewer showing point twenty is the control contradicting the picture.

    Undo steps back through the walk's own landmarks -- where it started, the
    highest point it crossed -- and those are geometries of the walk too, so
    the slider follows the box to them.
    """
    part, _state = _an_editor()
    points = _a_walk_of(part)
    assert part.submit_walk_at.value == len(points) - 1

    part.coords_widget.value = points[2][2]
    part._refresh_the_walk_points()
    assert part.submit_walk_at.value == 2


def test_stepping_leaves_the_slider_where_the_user_put_it():
    """Looking at a point deliberately does not write the box, so a refresh
    that happens while somebody is reading the walk must not pull the slider
    back to the point the box is standing on."""
    part, _state = _an_editor()
    _a_walk_of(part)
    part.submit_walk_at.value = 1
    part._refresh_the_walk_points()
    assert part.submit_walk_at.value == 1


def test_reading_the_walk_does_not_end_the_picture_of_it():
    """Reported from a real session: the slider moved one notch and the
    switch standing next to it vanished.

    The two are one walk shown two ways -- as a shape and as structures.
    Looking at a point hands the page one geometry that walk itself computed;
    it writes no box and takes no undo step, so the structure the profile is
    a claim about has not moved, and the claim is still true.
    """
    part, state = _an_editor()
    _a_walk_of(part)
    _draw_one(part)
    assert _shown(part.submit_scan_plot_btn) and _shown(part.submit_walk_at)

    part.submit_walk_at.value = 3
    assert _shown(part.submit_scan_plot_btn), (
        'one nudge of the slider took the profile away')
    assert state.get('scan_plot')
    assert 'base64' in part.submit_scan_plot.value


def test_going_on_from_a_point_takes_the_walk_off_the_row():
    """A trajectory describes the structure it was walked on. Once something
    has drawn over that structure there is nothing under what is on screen to
    step through, and the slider goes the way the picture does.

    It has to go at the moment the run *starts*, not when it writes: a run
    writes the box only when it finishes, so in between the geometry in the
    box is still the walk's own while the structure on screen is already
    somewhere else. That is why a refresh must not bring it back.
    """
    part, state = _an_editor()
    _a_walk_of(part)
    assert _shown(part.submit_walk_at)

    part._note_the_run(int(state.get('gfn_run', 0)) + 1, 'settle')
    assert not _shown(part.submit_walk_at)
    part._refresh_the_walk_points()
    assert not _shown(part.submit_walk_at), (
        'the box still holds the walk, but the structure has moved on')


def test_a_geometry_the_walk_never_computed_takes_it_off_too():
    part, _state = _an_editor()
    points = _a_walk_of(part)
    part.coords_widget.value = "2\nedited\nC 0.000 0.000 0.000\nO 1.421 0.000 0.000\n"
    part._refresh_the_walk_points()
    assert not _shown(part.submit_walk_at), points[0]


def test_the_presses_that_draw_nothing_leave_the_walk_alone():
    """Undo, Reset and the Stops claim a run number so that the page drops
    what it was playing; looking at a point claims one so the page shows the
    one frame handed to it. None of the three draws over the structure."""
    for name in ('press', 'abandoned', 'look'):
        part, state = _an_editor()
        _a_walk_of(part)
        part._note_the_run(int(state.get('gfn_run', 0)) + 1, name)
        assert _shown(part.submit_walk_at), name


def test_a_driven_walk_keeps_the_points_it_paid_for():
    """Which is why the box is read as geometry rather than as a comment.

    The comment would have been cheaper, and it is what the picture is kept
    honest by, but a walk does not sign the box the same way every time: a
    driven scan ends on "Driven until the bonds were made and broken", and a
    rule reading for "Scanned" alone would take the points away from exactly
    the walks that cost the most to compute.
    """
    from delfin.dashboard.structure_editor import xyz_document

    part, _state = _an_editor()
    points = _a_walk_of(part)
    rows = [line for line in points[-1][2].splitlines()[2:] if line.strip()]
    part.coords_widget.value = xyz_document(
        rows, 'Driven until the bonds were made and broken')
    part._refresh_the_walk_points()
    assert _shown(part.submit_walk_at)
    assert part.submit_walk_at.value == len(points) - 1


def test_a_point_asked_for_brings_the_structure_back_in_front():
    """The profile and the structure share one corner of the panel and only
    one of them is in it at a time, so a geometry asked for while the picture
    is up would otherwise be sent to a viewer nobody can see."""
    part, _state = _an_editor()
    _a_walk_of(part)
    _draw_one(part)
    part.submit_scan_plot_btn.value = True
    assert not _shown(part.mol_output)

    part.submit_walk_at.value = 2
    assert part.submit_scan_plot_btn.value is False
    assert _shown(part.mol_output), 'the structure is in front again'
    assert _shown(part.submit_scan_plot_btn), 'and the picture is still there'


def test_the_slider_stands_beside_the_profile_switch():
    """The two are one question asked twice: the profile is the walk as a
    shape, the slider is the walk as structures."""
    part, _state = _an_editor()
    row = list(part.submit_manip_toolbar.children)
    assert part.submit_walk_at in row
    assert row.index(part.submit_walk_at) == row.index(
        part.submit_scan_plot_btn) + 1


# ---------------------------------------------------------------------------
# The end of the stream
# ---------------------------------------------------------------------------
#
# Reported three times in one afternoon, from two installs and two Python
# versions, with the same twelve numbers each time: "the scan result jumps
# when I grab it".  The journal says why.  A walk of sixteen atoms delivers
# its fifteenth and last frame at 7.88 s and the box holds the structure it
# ended on at 8.04 s -- and at 11.05 s the picture is still standing on frame
# 3, because a walk is paced to be watched and the calculation was faster
# than the watching.
#
# That lag is not the defect; it is the point of pacing.  The defect is that
# the page was never told the run had ended.  Both optimisers say so through
# the final write of _stream_frames; the scan writes its own payloads, one
# per point, and every one of the fifteen carried final: null -- against the
# optimisation two seconds earlier in the same journal, whose last write
# carried final: 1.  Without it the player cannot tell a run that finished
# from one that stopped, and those mean opposite things about a queue.


def _payload_of(part):
    import json
    return json.loads(part.submit_gfn_frame.value or '{}')


def test_a_walk_that_has_ended_says_so_on_the_channel():
    part, state = _an_editor()
    run = part._note_the_run(int(state.get('gfn_run', 0)) + 1, 'scan')
    state['gfn_run'] = run

    part._close_the_frames(run)
    said = _payload_of(part)
    assert said.get('final') == 1, said
    assert said.get('run') == run


def test_the_marker_carries_nothing_and_moves_nothing():
    """A write that also claimed to deliver the tail would tell the page it
    had shown frames it has not, and the rest of the walk would be skipped --
    which is the opposite of what this is for."""
    part, state = _an_editor()
    run = part._note_the_run(int(state.get('gfn_run', 0)) + 1, 'scan')
    state['gfn_run'] = run

    part._close_the_frames(run)
    said = _payload_of(part)
    assert said.get('frames') == []
    assert said.get('from') == 0
    assert said.get('follow') == 1, (
        'the scan streams as a follow, and the marker must not say otherwise')


def test_a_walk_that_was_superseded_closes_nothing():
    """The same rule every write on this channel answers to. A scan whose run
    was taken over while its last point was being priced would otherwise
    close a stream belonging to whatever replaced it."""
    part, state = _an_editor()
    run = part._note_the_run(int(state.get('gfn_run', 0)) + 1, 'scan')
    state['gfn_run'] = run
    part.submit_gfn_frame.value = ''

    part._note_the_run(run + 1, 'optimise')
    state['gfn_run'] = run + 1
    part._close_the_frames(run)
    assert part.submit_gfn_frame.value == ''


def test_every_way_out_of_a_walk_closes_its_stream():
    """Including the one that walked nothing: what this closes is the run,
    not the result."""
    from editor_source import EDITOR_SOURCE

    body = EDITOR_SOURCE.split('def _done(final=walked):', 1)[1]
    closes = body.index('_close_the_frames(')
    walked_nothing = body.index('if not path:')
    assert closes < walked_nothing, (
        'a walk that ended early left its run open on the page')


def test_the_view_panel_is_not_a_wall_across_the_picture():
    """The panel is not its controls.

    It sits over the picture at ``width: max-content`` -- which is the width of
    its widest row laid out *unbroken* -- with a cap at the width of the
    viewer.  The sliders in it are 250 px, so on a narrow viewer the box spans
    most of the picture, and every press that landed on its padding was a press
    the structure never saw: a rubber band could not be started at all.

    Measured in chromium on the real stylesheet, at a 500 px viewer, a press a
    short way down the middle of the picture:

        two buttons on the bottom row   panel 262 px, 52% of the width -> slider
        three buttons (Orders added)    panel 300 px, 60%             -> slider
        three buttons, panel passed through                           -> the picture

    So the wall was there before Orders widened it; the third button is what
    pushed it past the width people were working at.  The controls take their
    own presses back, which is what the panel is for.
    """
    playwright = pytest.importorskip(
        'playwright.sync_api', reason='needs a browser to lay the panel out')
    from delfin.dashboard.molecule_viewer import STRUCTURE_VIEWER_FULLSCREEN_CSS

    css = STRUCTURE_VIEWER_FULLSCREEN_CSS
    assert '.delfin-structure-view-over' in css

    def page(width):
        return (
            "<!doctype html><html><head><style>body{margin:0}"
            ".row{display:flex;gap:4px;flex-flow:row wrap;align-items:center}"
            ".slider{width:250px;height:24px;background:#ccd}"
            + css +
            "\n#stack{position:relative!important;height:560px;"
            f"width:{width}px!important;max-width:none!important}}"
            "</style></head><body>"
            "<div id='stack' class='delfin-structure-viewer-stack'>"
            "<div class='delfin-structure-view-over'>"
            "<div class='slider'></div><div class='slider'></div>"
            "<div class='slider'></div><div class='slider'></div>"
            "<div class='row'>"
            "<button style='height:30px;width:104px'>Dyn. bonds</button>"
            "<button style='height:30px;width:86px'>Orders</button>"
            "<button style='height:30px;width:90px'>Centre</button>"
            "</div></div></div></body></html>")

    with playwright.sync_playwright() as engine:
        try:
            browser = engine.chromium.launch()
        except Exception as exc:                     # no browser on this box
            pytest.skip(f'chromium unavailable: {exc}')
        try:
            view = browser.new_page(viewport={'width': 1600, 'height': 900})
            for width in (500, 640, 900):
                view.set_content(page(width))
                got = view.evaluate("""() => {
                  const st = document.getElementById('stack');
                  const p = document.querySelector('.delfin-structure-view-over');
                  const r = p.getBoundingClientRect(), s = st.getBoundingClientRect();
                  const hit = (fx, fy) => {
                    const el = document.elementFromPoint(s.left + s.width * fx,
                                                         s.top + s.height * fy);
                    return el ? (el.id || el.className || el.tagName) : 'nothing';
                  };
                  return {wide: r.width > s.width * 0.5,
                          picture: hit(0.5, 0.12), control: hit(0.8, 0.20)};
                }""")
                # Wherever the panel reaches, the picture is still reachable.
                assert got['picture'] == 'stack', (width, got)
                # And a press on a control in it is still that control's.
                assert got['control'] == 'BUTTON', (width, got)
            # The narrow case is the one that was broken, so it has to be the
            # one that is actually narrow: a check that never covers half the
            # picture is not testing the wall.
            view.set_content(page(500))
            covers = view.evaluate("""() => {
              const st = document.getElementById('stack');
              const p = document.querySelector('.delfin-structure-view-over');
              return p.getBoundingClientRect().width
                     / st.getBoundingClientRect().width; }""")
            assert covers > 0.5, covers
        finally:
            browser.close()
