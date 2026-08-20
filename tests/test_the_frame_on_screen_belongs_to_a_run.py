"""The picture is what gets edited, so the number that names it must be true.

The page counts the frames it has drawn and hands that count over when an atom
is picked up; the kernel reads it as a count and keeps ``walked[shown - 1]`` as
the structure.  Two things were wrong with that, and both were measured on a
real editor driven through the page's own messages, with a real browser playing
the frames and xtb computing them.

**The count was one short.**  The branch that draws the first frame of a run
has nothing to interpolate from, and it drew without counting.  Measured in
Chromium over ten frames written one at a time, the way a scan writes them: the
picture stood on the tenth and the page reported nine.  So every grab during a
run kept the frame *before* the one on screen -- and a grab while the picture
stood on the first frame of a run reported zero, which the kernel reads as "no
frame at all" and answers with wherever the calculation had got to instead.
Measured end to end on a strained butane under GFN-FF, taking hold while the
picture stood on the first frame drawn: the geometry drawn and the geometry the
kernel went on to hold were 3.234 A apart, largest atom.  With the count fixed,
0.000 A.

**The count named no walk.**  It travelled alone, and a count on its own is a
plausible index into any path there is.  After a scan the page is still
standing on the scan's trajectory: a minimisation that starts a moment later
claims a run number the page has not seen yet, and a hand arriving in that
window handed it "frame 29" -- which it applied to its own path, a geometry
nobody had seen.  The run travels with the number now and is checked against
the run whose frames are being cut, so a number left over from a walk that has
ended can only be refused, never misread.
"""

from __future__ import annotations

import pytest

from editor_source import EDITOR_SOURCE as SOURCE


def _watcher():
    return SOURCE.split('def _install_gfn_frame_watcher')[1].split('\n    def ')[0]


def test_every_frame_the_picture_draws_is_counted():
    """Including the first of a run, which is the one that drew without
    counting.  There are three places a frame leaves the queue -- the whole
    queue at once, however many steps are due, and the single first frame --
    and every one of them has to move the count on, or the number means a
    different thing depending on where in the run it is read."""
    watcher = _watcher()
    # Three places a frame leaves the queue, and each of them counts it.
    assert watcher.count('play.shown=(play.shown||0)+1') == 2, (
        'the two single-frame paths do not both count')
    assert 'play.shown=(play.shown||0)+play.queue.length' in watcher, (
        'taking the whole queue at once does not count what it took')
    first = watcher.split(
        'play.last=play.queue.shift(); show(null,play.last,1)')[1]
    assert 'play.shown=(play.shown||0)+1' in first[:200], (
        'the first frame of a run is drawn without being counted')


def test_the_grab_says_which_walk_the_frame_belongs_to():
    """A count is only a frame if the run it counts along comes with it."""
    watcher = _watcher()
    assert 'String(play.shown||0)+","+(play.run||0)' in watcher, (
        'the grab hands over a count with no walk attached to it')
    assert '" of run "+(play.run||0)' in watcher, (
        'the stop report names no run either')


def test_the_kernel_reads_the_frame_only_out_of_the_walk_it_names():
    """The whole of the rule, in the one function that answers it, and the
    two places that ask it -- a run that was stopped and a run a hand cut
    off.  Both used to index whatever number happened to be lying about."""
    reader = SOURCE.split('def _the_shown_frame_of')[1] \
                   .split('\n    def ')[0]
    assert "state.get('gfn_shown_run') != int(run or 0)" in reader
    # The arithmetic that turns the count into a frame lives in one place,
    # which is not this one: three copies of an off-by-one is three chances
    # to hand back a geometry nobody was looking at.
    cutter = SOURCE.split('def _frame_as_xyz(')[1].split('\n    def ')[0]
    assert 'walked[shown - 1]' in cutter
    handler = SOURCE.split(
        'def on_submit_optimize(change=None, every_frame=False)')[1] \
        .split('\n    def ')[0]
    assert handler.count('_the_shown_frame_of(run_id)') == 2, (
        'one of the two readers takes the count without asking whose walk '
        'it counts along')
    assert 'walked[shown - 1]' not in handler, (
        'a reader still indexes the path itself')


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
    refs['coords_widget'].value = (
        '3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n')
    return refs


def _says(refs, verb, payload, serial=[0]):          # noqa: B006
    serial[0] += 1
    refs['submit_cmd_sync'].value = f'{verb}:{serial[0]}:{payload}'


def test_a_grab_records_the_frame_and_the_walk_it_belongs_to(editor):
    """Driven the way the browser drives it: the message, not the internals."""
    state = editor['editor_state']
    _says(editor, 'gfngrab', '29,7')
    assert state['gfn_shown_frame'] == 29
    assert state['gfn_shown_run'] == 7


def test_a_stop_report_records_both_as_well(editor):
    """The other way the page names a frame -- the switch going up rather than
    a hand arriving -- and it used to be parsed by taking the last word of the
    sentence, which is now the run."""
    state = editor['editor_state']
    _says(editor, 'gfnplay', 'stopped at frame 12 of run 4')
    assert state['gfn_shown_frame'] == 12
    assert state['gfn_shown_run'] == 4


def test_a_page_that_names_no_walk_is_refused_rather_than_guessed_at(editor):
    """An older page sends the count alone.  Both are dropped rather than one
    of them kept: an index into the wrong path is worse than no index at all,
    because the wrong one is acted on and the missing one is not."""
    state = editor['editor_state']
    _says(editor, 'gfngrab', '12,3')
    _says(editor, 'gfngrab', '12')
    assert state.get('gfn_shown_frame') is None
    assert state.get('gfn_shown_run') is None


def test_a_walk_that_ends_takes_its_frame_number_with_it(editor):
    """Starting a run drops whatever the last one left behind, and so does
    taking a structure back.  Belt as well as braces now that the run travels
    with the number, and both are worth keeping: the number is written by the
    page and the guard is the only thing standing between a stale one and a
    path it does not describe."""
    state = editor['editor_state']
    _says(editor, 'gfngrab', '29,7')
    editor['submit_manip_undo_btn'].click()
    assert state.get('gfn_shown_frame') is None
    assert state.get('gfn_shown_run') is None
