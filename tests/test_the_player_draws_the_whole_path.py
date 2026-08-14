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


def _play(browser, writes, settle_ms, pace=None):
    """Write those payloads into the frame field and report what was drawn."""
    script = _player_script()
    scope = script.split('var scope=')[1].split(';')[0].strip().strip('"')
    page = browser.new_page()
    try:
        page.set_content(PAGE.replace("__SCOPE__", scope))
        page.evaluate(script)
        if pace is not None:
            page.evaluate("([s, ms]) => {window.__delfinGfnPlay[s].pace = ms;}",
                          [scope, pace])
        for payload, wait in writes:
            page.evaluate(
                """([sel, text]) => {
                    document.querySelector(sel).value = text;
                }""",
                [".submit-gfn-frame textarea", json.dumps(payload)])
            page.wait_for_timeout(wait)
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
    out = _play(browser, _bursts(60, 20), settle_ms=6000, pace=8)
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
        page.wait_for_timeout(900)

        queued = page.evaluate("s => window.__delfinGfnPlay[s].queue.length",
                               scope)
        shown = page.evaluate("s => window.__delfinGfnPlay[s].shown", scope)
        assert queued > 20, "the calculation has to be ahead for this to mean anything"
        assert 0 < shown < 120, shown
        standing = page.evaluate("window.__drawn[window.__drawn.length - 1][0]")

        # A hand arrives on the structure.
        page.evaluate("""s => {window._submitManipStateByScope[s] =
            {drag: {kind: 'translate', targets: [1]}};}""", scope)
        page.wait_for_timeout(500)

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
        assert payload == str(shown), (
            f"the kernel was told {payload!r}, the picture stood on {shown}")
    finally:
        page.close()
