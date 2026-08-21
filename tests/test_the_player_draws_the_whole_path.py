"""The trajectory player, driven in a real browser.

Everything else about the player is checked by reading the script it ships.
That catches a rule that was removed and misses a rule that does not do what
it says, and both of those were in here: the kernel sent a fixed tail of eight
frames however many had been computed, and the page threw away everything past
a twenty-frame queue.  Neither is visible in the source without knowing what
the numbers cost, and both showed up the moment the thing was run --
thirty-five of sixty frames drawn, the oldest of every burst missing.

So this loads the real script into Chromium against a stubbed viewer that
records every frame it is asked to draw, and writes the frame field the way
ipywidgets writes it: the value set directly, no event.  What is measured is
the path from the kernel's write to the picture.

Skipped where there is no browser, which is most machines and the CI runner.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from delfin.dashboard.context import DashboardContext


PAGE = """<!doctype html><html><body>
<div class="__SCOPE__">
  <div class="submit-gfn-frame"><textarea></textarea></div>
  <div class="submit-cmd-sync"><input></div>
  <button class="submit-optimize-switch mod-active">Optimize</button>
</div>
<script>
window.__drawn = [];
window._submitManipStateByScope = {};
window.__delfinSubmitManip = {
  setPositions: function(scope, out, held) {
    window.__drawn.push(out.slice());
    return true;
  },
  pushXyz: function(){ return true; }
};
</script>
</body></html>"""


def _player_script():
    """The whole self-executing player block, as the page receives it."""
    from delfin.dashboard import tab_submit

    folder = Path(tempfile.mkdtemp())
    for name in ("calc", "archive", "office"):
        (folder / name).mkdir()
    ctx = DashboardContext(calc_dir=folder / "calc",
                           archive_dir=folder / "archive",
                           office_dir=folder / "office")
    ctx.run_js = lambda _script: None
    tab_submit.create_tab(ctx)
    startup = "\n".join(ctx.init_js_parts)
    at = startup.index("window.__delfinGfnPlay")
    return startup[startup.rindex("(function(){", 0, at):]


@pytest.fixture(scope="module")
def browser():
    pytest.importorskip("ipywidgets")
    playwright = pytest.importorskip("playwright.sync_api")
    with playwright.sync_playwright() as p:
        try:
            engine = p.chromium.launch()
        except Exception as exc:                      # no browser installed
            pytest.skip(f"chromium not available: {exc}")
        yield engine
        engine.close()


def _play(browser, writes, settle_ms, pace=None, switch_off_after=None):
    """Write those payloads into the frame field and report what was drawn.

    *switch_off_after* takes the Optimise switch up after that many writes,
    which is what the kernel does the moment the last round lands -- while the
    picture is still walking the path.
    """
    script = _player_script()
    scope = script.split('var scope=')[1].split(';')[0].strip().strip('"')
    page = browser.new_page()
    try:
        page.set_content(PAGE.replace("__SCOPE__", scope))
        page.evaluate(script)
        if pace is not None:
            page.evaluate("([s, ms]) => {window.__delfinGfnPlay[s].pace = ms;}",
                          [scope, pace])
        for written, (payload, wait) in enumerate(writes, 1):
            page.evaluate(
                """([sel, text]) => {
                    document.querySelector(sel).value = text;
                }""",
                [".submit-gfn-frame textarea", json.dumps(payload)])
            page.wait_for_timeout(wait)
            if switch_off_after is not None and written == switch_off_after:
                page.evaluate("""() => {
                    document.querySelector('.submit-optimize-switch')
                        .classList.remove('mod-active');
                }""")
        page.wait_for_timeout(settle_ms)
        return {
            "frames": page.evaluate("window.__drawn"),
            "shown": page.evaluate("s => window.__delfinGfnPlay[s].shown", scope),
            "queue": page.evaluate("s => window.__delfinGfnPlay[s].queue.length",
                                   scope),
        }
    finally:
        page.close()


def _bursts(total, size):
    """The path, written the way the kernel writes it: each window starts
    where the previous window started, so every frame goes out twice."""
    path = [[float(i)] * 9 for i in range(total)]
    writes, start, end = [], 0, 0
    for cut in range(size, total + 1, size):
        writes.append(({"run": 1, "from": start, "frames": path[start:cut]}, 120))
        start, end = end, cut
    writes.append(({"run": 1, "from": start, "frames": path[start:],
                    "final": 1}, 60))
    return writes


def _reached(drawn):
    """Which points of the path were drawn.  Every component of frame i is i,
    so a drawn frame names itself; the interpolated positions in between fall
    between two integers and are not counted."""
    return {round(f[0]) for f in drawn if abs(f[0] - round(f[0])) < 1e-6}


def test_every_frame_of_the_path_reaches_the_picture(browser):
    """Thirty-five of sixty, before this.

    The kernel sent a fixed tail of eight frames per write, and the page kept
    a queue of twenty and dropped the front of it.  Two independent rules,
    each defensible on its own, and between them the viewer showed a sample of
    the optimisation rather than the optimisation.
    """
    # Forty milliseconds is twenty-five a second, the top of the control.
    # Asking for faster than the screen draws could only mean showing fewer
    # frames, and the control does not reach there -- so the promise that
    # every frame is drawn is tested at the fastest that can be asked for.
    out = _play(browser, _bursts(60, 20), settle_ms=6000, pace=40)
    reached = _reached(out["frames"])

    missing = sorted(set(range(60)) - reached)
    assert not missing, f"never drawn: {missing}"
    assert out["queue"] == 0, "the path was not played to the end"
    assert out["frames"][-1][0] == 59.0, "the picture stops short of the end"


def test_the_pace_control_sets_the_rate(browser):
    """Frames a second, and it is the user's number that decides.

    The player used to speed up whenever the queue grew, which is exactly what
    a slow setting asks it not to do.
    """
    path = [[float(i)] * 9 for i in range(120)]
    one = [({"run": 1, "from": 0, "frames": path, "final": 1}, 50)]

    fast = _play(browser, one, settle_ms=1500, pace=20)
    slow = _play(browser, one, settle_ms=1500, pace=200)

    # 1000/20 = 50 a second against 1000/200 = 5, both inside a wide band:
    # the animation frame is 60 Hz, so the fast end is bounded by that and
    # neither end is exact on a shared machine.
    assert 30 <= fast["shown"] <= 90, fast["shown"]
    assert 3 <= slow["shown"] <= 12, slow["shown"]
    assert fast["shown"] > slow["shown"] * 3, (
        f"the pace did not pace: {fast['shown']} against {slow['shown']}")


def test_a_hand_on_the_structure_cuts_the_run_where_the_picture_stands(browser):
    """Taking hold mid-playback keeps the frame on screen and drops the rest.

    The frames in front of the picture were computed for a structure the hand
    is now changing, and nobody has seen them.  Dropping them was already
    right; saying *where* was not -- the message went out as "gfngrab:3:" with
    no frame on it, so the kernel knew a hand had arrived and not where the
    picture stood.  The geometry the user had taken hold of then lived only in
    the browser until some later drag pushed it back, and taking hold and
    letting go without moving left the box and the picture apart.
    """
    script = _player_script()
    scope = script.split('var scope=')[1].split(';')[0].strip().strip('"')
    path = [[float(i)] * 9 for i in range(120)]

    page = browser.new_page()
    try:
        page.set_content(PAGE.replace("__SCOPE__", scope))
        page.evaluate(script)
        page.evaluate("([s, ms]) => {window.__delfinGfnPlay[s].pace = ms;}",
                      [scope, 60])
        page.evaluate("([sel, t]) => {document.querySelector(sel).value = t;}",
                      [".submit-gfn-frame textarea",
                       json.dumps({"run": 1, "from": 0, "frames": path,
                                   "final": 1})])
        # Waited for rather than slept through.  A fixed window assumes the
        # browser gets the processor within it, and on a machine running the
        # rest of the suite it does not -- the picture had drawn nothing and
        # the test read that as the player being broken.
        page.wait_for_function(
            "s => window.__delfinGfnPlay[s].shown > 3", arg=scope, timeout=30000)

        queued = page.evaluate("s => window.__delfinGfnPlay[s].queue.length",
                               scope)
        shown = page.evaluate("s => window.__delfinGfnPlay[s].shown", scope)
        assert queued > 20, "the calculation has to be ahead for this to mean anything"
        assert 0 < shown < 120, shown
        standing = page.evaluate("window.__drawn[window.__drawn.length - 1][0]")

        # A hand arrives on the structure.  movedEnough, because that is what
        # tells a drag from a tap now: a press that has not passed the slop has
        # moved nothing, and the player leaves it alone rather than cutting the
        # run under it.
        page.evaluate("""s => {window._submitManipStateByScope[s] =
            {drag: {kind: 'translate', targets: [1], movedEnough: true}};}""",
                      scope)
        page.wait_for_function(
            "s => window.__delfinGfnPlay[s].queue.length === 0",
            arg=scope, timeout=30000)

        assert page.evaluate(
            "s => window.__delfinGfnPlay[s].queue.length", scope) == 0, (
            "what was computed past the picture has to go")
        # And the picture does not jump forward to where the calculation got to.
        after = page.evaluate("window.__drawn[window.__drawn.length - 1][0]")
        assert abs(after - standing) < 2.0, (
            f"the picture moved on after the grab: {standing} -> {after}")

        said = page.evaluate(
            "document.querySelector('.submit-cmd-sync input').value")
        verb, _serial, payload = said.split(":", 2)
        assert verb == "gfngrab"
        # The count and the walk it counts along.  The count alone is a
        # plausible index into any path there is, and after a scan it was
        # being read into somebody else's.
        assert payload == f"{shown},1", (
            f"the kernel was told {payload!r}, the picture stood on {shown} "
            f"of run 1")
    finally:
        page.close()


def test_the_count_the_page_reports_is_the_number_of_frames_it_has_drawn(
        browser):
    """A count, not an index, and the kernel reads it as ``walked[shown-1]``.

    The branch that draws the first frame of a run has nothing to interpolate
    from and drew without counting, so the number was one short for the whole
    of the run after it.  Measured here over ten frames written one at a time,
    the way a scan writes them: the picture stood on the tenth of the path and
    the page reported nine, so a grab kept the ninth -- and a grab while the
    picture stood on the first frame reported zero, which the kernel reads as
    "no frame at all" and answers with wherever the calculation had got to.

    Frames name themselves, so the picture can be read off the last draw:
    every component of frame i is i, and the count has to be i + 1.
    """
    writes = [({"run": 1, "from": i, "follow": 1,
                "frames": [[float(i)] * 9]}, 150) for i in range(10)]
    out = _play(browser, writes, settle_ms=2000, pace=55)

    assert out["queue"] == 0, "the path was not played to the end"
    standing = out["frames"][-1][0]
    assert standing == 9.0, f"the picture stands on frame {standing}"
    assert out["shown"] == 10, (
        f"ten frames were drawn and the page reports {out['shown']}")


def test_a_finished_run_plays_its_path_out_after_the_switch_goes_up(browser):
    """The run being finished does not mean the picture is.

    The kernel turns the Optimise switch off the moment the last round lands,
    and the player abandoned whatever was still queued when it saw that go up.
    While the queue was capped at twenty and a backlog was played at speed,
    almost nothing was left by then.  Keeping the whole path and walking it at
    a pace anyone can watch made that the ordinary case: measured here at
    twelve frames a second, eight of sixty drawn and fifty-two thrown away the
    instant the switch moved -- the picture showing the first jerk of an
    optimisation and never arriving at the minimum it had just computed.

    A finished run and a stopped one are different things.  The last write of
    a finished run is marked; a stop arrives as a halt, which clears the queue
    on purpose.
    """
    path = [[float(i)] * 9 for i in range(60)]
    out = _play(browser,
                [({"run": 1, "from": 0, "frames": path, "final": 1}, 600)],
                settle_ms=7000, pace=83, switch_off_after=1)

    assert out["queue"] == 0, "it has to finish walking"
    assert _reached(out["frames"]) == set(range(60)), "the whole path"
    assert out["frames"][-1][0] == 59.0, (
        "the picture has to arrive at the geometry the run reached")


def test_the_mark_survives_the_write_that_carries_it(browser):
    """A run short enough to finish in one message carries the mark on the
    same write that starts the run, and the run-change reset ran after it --
    so the mark was cleared by the very write that brought it, and the path
    was abandoned exactly as before.
    """
    script = _player_script()
    body = script.split("if(run!==play.run){")[1]
    reset_at = body.index("play.complete=0")
    mark_at = body.index("data.final) play.complete=1")
    assert reset_at < mark_at, (
        "the mark has to be set after the run reset, not before it")


def test_every_setting_of_the_pace_makes_a_difference(browser):
    """Most of the slider did nothing.

    The loop runs on the animation frame and took one step per frame, resetting
    its clock to now each time -- so the remainder was thrown away and the rate
    snapped to whole divisors of the screen.  Measured: 62.7 a second, then
    31.3, and nothing in between.  Every setting from 32 to 59 drew 31.3, which
    is most of the top half of the control.

    The clock moves on by exactly the steps taken now, so the remainder carries
    and the average is the rate that was asked for.
    """
    path = [[float(i)] * 9 for i in range(400)]
    one = [({"run": 1, "from": 0, "frames": path, "final": 1}, 50)]

    rates = {}
    for wanted in (5, 10, 17, 25):
        out = _play(browser, one, settle_ms=1500,
                    pace=max(1, round(1000 / wanted)))
        rates[wanted] = out["shown"] / 1.5

    # Each step up is a real step up, not a repeat of the one below it.
    for lower, higher in ((5, 10), (10, 17), (17, 25)):
        assert rates[higher] > rates[lower] * 1.15, (
            f"{higher}/s drew {rates[higher]:.1f} against {lower}/s at "
            f"{rates[lower]:.1f} -- the setting did nothing")


def test_the_pace_reaches_down_to_one_frame_every_ten_seconds(browser):
    """Slow is the useful end: the calculation runs on ahead while the picture
    walks, and taking hold keeps the frame that is on screen.

    Half a frame a second was the floor and it was not slow enough to watch a
    coordination sphere rearrange a step at a time.
    """
    path = [[float(i)] * 9 for i in range(60)]
    one = [({"run": 1, "from": 0, "frames": path, "final": 1}, 50)]

    # The count includes the frame the run opens on, which is drawn as soon as
    # it arrives and has nothing to be paced against -- so three paced steps
    # read as four.
    out = _play(browser, one, settle_ms=6000, pace=2000)   # half a frame a second
    assert 3 <= out["shown"] <= 5, (
        f"half a frame a second over six seconds is three, not {out['shown']}")
    assert out["queue"] > 50, "and the rest is still waiting to be walked"

    # And the floor: one every ten seconds, so six seconds moves nothing on
    # past the frame the run opened on.
    crawl = _play(browser, one, settle_ms=6000, pace=10000)
    assert crawl["shown"] == 1, crawl["shown"]
    assert crawl["queue"] >= 58, "the path is all still ahead"


def test_the_number_that_travels_for_a_held_atom_is_its_index(browser):
    """A serial and an index are different numbers, and both are plausible.

    pushXyzToPython names the atoms the hand is on with ffIndicesOf, so what
    reaches the kernel is a position in getAtoms -- not the atom's serial.  The
    thermal wall was keyed by serial and looked up nothing wherever the two
    differ, which is the ordinary case: a benzene could still be pulled open
    while the budget underneath said it could not.

    Measured rather than read, because this is a contract between two files and
    reading either one on its own agrees with itself.
    """
    from delfin.dashboard import tab_submit

    scope = "sc"
    page = browser.new_page()
    try:
        page.set_content(
            f'<!doctype html><html><body><div class="{scope}">'
            '<div class="submit-manip-sync"><textarea></textarea></div></div>'
            '<script>'
            'window._submitManipStateByScope={};'
            # Serials that are deliberately not the indices, as 3Dmol hands
            # them out once a model has been rebuilt.
            'window.__ATOMS=[{serial:101,elem:"C",x:0,y:0,z:0},'
            '{serial:102,elem:"C",x:1.65,y:0,z:0}];'
            'var MODEL={selectedAtoms:function(){return window.__ATOMS;}};'
            f'window._submitMolViewerByScope={{"{scope}":'
            '{getModel:function(){return MODEL;},render:function(){}}};'
            '</script></body></html>')
        page.evaluate(tab_submit.submit_manip_bootstrap_js())

        # Split here rather than in the page: a newline escape inside a
        # Python string is a newline before the browser ever sees it, and a
        # newline inside a JavaScript string literal is a syntax error.
        pushed = page.evaluate("""(s) => {
            window._submitManipStateByScope[s] =
                {drag: {kind: 'translate', targets: [102]}};
            window.__delfinSubmitManip.pushXyz(s, 'drag-follow');
            return document.querySelector('.submit-manip-sync textarea').value;
        }""", scope)
        note = pushed.splitlines()[1]

        assert "held=1" in note, note
        assert "held=102" not in note, (
            "the kernel is given the index, so the wall must be keyed by it")
    finally:
        page.close()


def test_the_wall_reaches_the_page_and_can_be_taken_down(browser):
    """It was sent with run_js and never arrived.

    run_js clears its output before displaying, so of two calls in quick
    succession only the second survives -- and the kernel writes the wall
    about ten times a second while a drag is running.  The budget said the
    structure was fifty-five kcal/mol past what it can spend and the ring came
    apart anyway, because the holding-back never left the kernel.

    It travels as a widget value the page reads on its own clock now, the same
    way the trajectory does, in a field of its own: the frame field is one
    slot too, and a wall written between two frames would be gone before the
    page looked.
    """
    import json as _json

    from delfin.dashboard import tab_submit

    script = _player_script()
    scope = script.split('var scope=')[1].split(';')[0].strip().strip('"')
    page = browser.new_page()
    try:
        page.set_content(
            f'<!doctype html><html><body><div class="{scope}">'
            '<div class="submit-gfn-frame"><textarea></textarea></div>'
            '<div class="submit-gfn-wall"><input></div>'
            '<div class="submit-cmd-sync"><input></div>'
            '<button class="submit-optimize-switch mod-active">O</button></div>'
            '<script>window._submitManipStateByScope={};'
            'window.__ATOMS=[{serial:101,elem:"C",x:0,y:0,z:0},'
            '{serial:102,elem:"C",x:1.65,y:0,z:0}];'
            'var MODEL={selectedAtoms:function(){return window.__ATOMS;}};'
            f'window._submitMolViewerByScope={{"{scope}":'
            '{getModel:function(){return MODEL;},render:function(){}}};'
            '</script></body></html>')
        page.evaluate(tab_submit.submit_manip_bootstrap_js())
        page.evaluate(script)

        def wall():
            return page.evaluate(
                "s => (window._submitManipStateByScope[s] || {}).thermalWall"
                " || null", scope)

        def write(payload):
            page.evaluate("([sel, t]) => {document.querySelector(sel).value = t;}",
                          [".submit-gfn-wall input", _json.dumps(payload)])

        assert wall() is None

        # Keyed by index, which is the number the kernel is given.
        write({"n": 1, "wall": {"1": [1.65, 0.0, 0.0]}})
        page.wait_for_function(
            "s => !!(window._submitManipStateByScope[s] || {}).thermalWall",
            arg=scope, timeout=15000)
        assert wall() == {"1": [1.65, 0, 0]}

        # And down again on the far side of a barrier, where the energy falls
        # and the structure is back inside its budget.
        write({"n": 2, "wall": None})
        page.wait_for_function(
            "s => !(window._submitManipStateByScope[s] || {}).thermalWall",
            arg=scope, timeout=15000)
        assert wall() is None
    finally:
        page.close()


def _drag_page(browser, scope_holder):
    """A page with the manipulation script live and its mouse handlers bound.

    onViewerReady is what attaches the overlay and binds them; without it the
    events go nowhere and every assertion about a drag passes because nothing
    moved at all.  That is how the first version of this test read.
    """
    from delfin.dashboard import tab_submit

    scope = "sc"
    scope_holder.append(scope)
    page = browser.new_page()
    page.set_content(
        f'<!doctype html><html><body><div class="{scope}">'
        '<div class="submit-mol-output" '
        'style="position:relative;width:600px;height:600px">'
        '<canvas width="600" height="600"></canvas></div>'
        '<div class="submit-manip-sync"><textarea></textarea></div></div>'
        '<script>window._submitManipStateByScope={};'
        'window.__ATOMS=[{serial:101,elem:"C",x:0,y:0,z:0},'
        '{serial:102,elem:"C",x:1.379,y:0,z:0}];'
        'var MODEL={selectedAtoms:function(){return window.__ATOMS;}};'
        f'window._submitMolViewerByScope={{"{scope}":'
        '{getModel:function(){return MODEL;},render:function(){}}};'
        '</script></body></html>')
    page.evaluate(tab_submit.submit_manip_bootstrap_js())
    page.evaluate("(s) => {window.__delfinSubmitManip.onViewerReady("
                  "s, document.querySelector('.submit-mol-output'));}", scope)
    return page


def test_the_hand_cannot_outrun_the_calculation(browser):
    """A mouse crosses an angstrom in a hundred milliseconds.

    The energy behind it arrives about ten times a second, so a drag that is
    simply applied jumps from one geometry to the next and the structures in
    between are never evaluated.  Measuring only where the hand ends up cannot
    say whether the path was possible: the ring is open before anything can
    object, and what a budget then refuses is a distortion that has already
    happened.

    Measured here on the real drag path -- the same mousemove handler the page
    binds -- with twenty events of two pixels, which is what an ordinary drag
    is made of:

        no leash          +1.20 A     the ring comes apart
        leash of 0.10 A   +0.06 A     it waits for the calculation

    The mark advances with every answer that is inside the budget, so in use
    the atom keeps moving -- at the rate the path can be checked, which is
    what walking it means.
    """
    holder = []
    page = _drag_page(browser, holder)
    scope = holder[0]
    try:
        def where():
            return round(page.evaluate("window.__ATOMS[1].x"), 3)

        def start_at(value):
            page.evaluate("v => {window.__ATOMS[1].x = v;}", value)

        def leash(mark, reach):
            page.evaluate(
                "([s, m, r]) => {window.__delfinSubmitManip.setThermalWall("
                "s, m === null ? null : {'1': [m, 0, 0]}, r);}",
                [scope, mark, reach])

        def drag(pixels_each, times):
            page.evaluate("""([s, d, n]) => {
                var st = window._submitManipStateByScope[s];
                st.dragSensitivity = 1;
                var X = 100;
                st.drag = {kind: 'translate', targets: [102],
                           movedEnough: true, snapshotted: true,
                           lastX: X, lastY: 100, startX: X, startY: 100};
                for (var i = 0; i < n; i++) {
                    X += d;
                    window.dispatchEvent(new MouseEvent('mousemove',
                        {clientX: X, clientY: 100, bubbles: true}));
                }
            }""", [scope, pixels_each, times])

        # Nothing holding it: the hand takes the atom wherever it likes.
        start_at(1.379)
        leash(None, 0)
        drag(2, 20)
        loose = where()
        assert loose > 2.3, (
            f"the drag has to be able to tear it open, or this proves nothing: "
            f"{loose}")

        # On the leash, with the mark standing still because no answer has
        # come back yet: it advances to the edge and waits.
        start_at(1.379)
        leash(1.379, 0.10)
        drag(2, 20)
        held = where()
        assert held - 1.379 <= 0.101, f"the leash let go: {held}"
        assert held > 1.379, "and it must still follow the hand as far as it may"

        # Back is always possible: a structure that has run out of budget has
        # to be undoable, and that is the one thing needed from there.
        before = where()
        drag(-2, 20)
        assert where() < before, f"the way back was blocked: {before} -> {where()}"
    finally:
        page.close()
