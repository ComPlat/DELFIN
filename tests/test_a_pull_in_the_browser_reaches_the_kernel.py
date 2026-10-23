"""A real press, a real canvas, a real mouse, and what the kernel is told.

Everything else about the hand is measured on the kernel's side, where the
arithmetic is: what the budget allows, what xtb answers, what the wall keeps.
None of that says the gesture works.  Between the mouse and the kernel there
is a raycaster, an overlay that has to take the press without stealing the
scene's rotation, a force field that decides how much of the wish the atoms
get, and one hidden textarea that is the whole channel back.  This drives that
half in a browser, with real mouse events on the real canvas.

What it pins down:

  the press finds an atom      through the shipped raycaster, not through a
                               projection written for the test
  the pull moves nothing here  the page records the wish and lets the field
                               answer -- see ffPullWant.  A page that moved
                               the atom itself would be showing travel the
                               chemistry never agreed to
  the kernel is told           through the .submit-manip-sync input, in the
                               format the kernel parses, carrying a geometry
                               that actually moved
  a still mouse asks nothing   no move, no report.  Every answer the kernel
                               gives is about a wish somebody made; a page
                               that re-asked on a timer would have the
                               structure recomputed under a hand that is not
                               moving, which is the shape of the shaking this
                               editor has spent a long time removing

Skipped where there is no browser, which is most machines.
"""

from __future__ import annotations

from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
_BUNDLE = _REPO / 'delfin' / 'dashboard' / 'static' / '3Dmol-min.js'
_SCOPE = 'submit-scope-browsertest'

_ETHANE = """8
ethane
C  0.000000  0.000000  0.000000
C  1.530000  0.000000  0.000000
H -0.360000  1.020000  0.000000
H -0.360000 -0.510000  0.884000
H -0.360000 -0.510000 -0.884000
H  1.890000  0.510000  0.884000
H  1.890000  0.510000 -0.884000
H  1.890000 -1.020000  0.000000
"""

#: Set up the editor over a viewer of our own, the way the tab does it.
_SETUP = """
(args) => {
  const [scope, xyz, terms] = args;
  const host = document.getElementById('host');
  const v = window.__delfinCreateViewer(host, {backgroundColor: 'white'});
  window._submitMolViewerByScope = window._submitMolViewerByScope || {};
  window._submitMolViewerByScope[scope] = v;
  v.addModel(xyz, 'xyz');
  v.setStyle({}, {stick: {radius: 0.15}, sphere: {scale: 0.28}});
  v.zoomTo(); v.render();
  const M = window.__delfinSubmitManip;
  M.onViewerReady(scope, host);
  M.setMode(scope, 'manipulate');
  M.setLiveHand(scope, true);
  M.setPullStrength(scope, 0.5, null);
  const ff = M.setForceField(scope, terms);
  M.setOptimizerStrength(scope, 3);
  M.setSettleOnRelease(scope, false);
  window.__reports = [];
  document.getElementById('sync').addEventListener('input', function () {
    window.__reports.push(this.value);
  });
  return ff && ff.ok !== false;
}
"""

_LOOK = """
(scope) => {
  const st = window._submitManipStateByScope[scope];
  const v = window._submitMolViewerByScope[scope];
  return {
    kind: st && st.drag ? st.drag.kind : null,
    share: st ? st.pullShare : null,
    xs: v.getModel().selectedAtoms({}).map((a) => +a.x.toFixed(3)),
    reports: window.__reports.length,
  };
}
"""


def _blocks():
    from delfin.dashboard.molecule_forcefield_js import molecule_ff_bootstrap_js
    from delfin.dashboard.molecule_viewer import (
        RIGHT_MOUSE_TRANSLATE_PATCH_JS, submit_manip_bootstrap_js,
        vendored_3dmol_js)

    return [vendored_3dmol_js(), RIGHT_MOUSE_TRANSLATE_PATCH_JS,
            molecule_ff_bootstrap_js(), submit_manip_bootstrap_js()]


def _grab(page, box):
    """The first point on the canvas where a press starts a drag.

    Scanned rather than projected.  Which pixel is over an atom is the
    raycaster's answer and nobody else's, and a projection written here to
    agree with it would be testing itself: measured, a hand-written
    modelToScreen put the carbon 62 pixels away from where the press had to
    land.
    """
    for gy in range(6, int(box['h']), 14):
        for gx in range(6, int(box['w']), 14):
            at = (box['l'] + gx, box['t'] + gy)
            page.mouse.move(*at)
            page.mouse.down()
            kind = page.evaluate(
                "(s) => {const st = window._submitManipStateByScope[s];"
                "return st && st.drag ? st.drag.kind : null;}", _SCOPE)
            page.mouse.up()
            if kind:
                return at
    return None


def test_a_pull_on_the_canvas_reaches_the_kernel_and_a_still_mouse_does_not():
    playwright = pytest.importorskip(
        'playwright.sync_api', reason='needs a browser to press an atom')
    if not _BUNDLE.is_file():
        pytest.skip('vendored 3Dmol bundle is not present')
    from delfin.dashboard.molecule_forcefield import export_forcefield_terms

    terms = export_forcefield_terms(_ETHANE, method='uff')
    assert terms.get('ok'), terms

    html = (
        f"<div class='{_SCOPE}'>"
        "<div id='host' style='width:640px;height:480px;position:relative'>"
        "</div>"
        "<div class='submit-manip-sync'><textarea id='sync'></textarea></div>"
        "<div class='submit-manip-status'></div></div>"
        # type="module" is the host's own scope: `this` undefined, and a
        # top-level var does not become a global.  See
        # test_viewer_scripts_in_host_scope.
        + ''.join(f'<script type="module">{one}</script>' for one in _blocks()))

    errors: list[str] = []
    with playwright.sync_playwright() as p:
        try:
            browser = p.chromium.launch(args=['--use-gl=swiftshader'])
        except Exception as exc:
            pytest.skip(f'chromium unavailable: {exc}')
        try:
            page = browser.new_page(viewport={'width': 900, 'height': 700})
            page.on('pageerror', lambda e: errors.append(str(e)))
            page.set_content(html)
            page.wait_for_timeout(3000)
            assert page.evaluate(_SETUP, [_SCOPE, _ETHANE, terms]) is True
            page.wait_for_timeout(400)

            box = page.evaluate(
                "() => {const b = document.getElementById('host')"
                ".getBoundingClientRect();"
                "return {l: b.left, t: b.top, w: b.width, h: b.height};}")
            at = _grab(page, box)
            assert at is not None, 'no press anywhere on the canvas took an atom'

            page.mouse.move(*at)
            page.mouse.down()
            page.wait_for_timeout(120)
            assert page.evaluate(_LOOK, _SCOPE)['kind'] == 'translate'

            was = page.evaluate(_LOOK, _SCOPE)['xs']
            for step in range(1, 11):
                page.mouse.move(at[0] - 12 * step, at[1])
                page.wait_for_timeout(60)
            pulled = page.evaluate(_LOOK, _SCOPE)

            # The hand is a pull, and the atoms have moved -- but they have
            # moved because the field answered, not because the page put them
            # there.  What the page contributes is the wish.
            assert pulled['share'] == pytest.approx(0.5)
            assert pulled['xs'] != was
            assert pulled['reports'] >= 1, 'the kernel was never told'

            # And a mouse that is not moving asks for nothing.  Every answer
            # the kernel gives is about a wish somebody made.
            page.wait_for_timeout(1500)
            assert page.evaluate(_LOOK, _SCOPE)['reports'] == pulled['reports']

            page.mouse.up()
            page.wait_for_timeout(400)
            said = page.evaluate('() => window.__reports')
            assert said, 'the release told the kernel nothing'
            assert 'drag-end' in said[-1].splitlines()[1], said[-1]
            # And what it hands over is a geometry that moved, in the format
            # the kernel parses: the release is the last word on a drag, and a
            # release that reported the structure it started from would undo
            # the whole gesture.
            rows = [line for line in said[-1].splitlines()[2:] if line.strip()]
            assert len(rows) == 8, said[-1]
            assert [float(one.split()[1]) for one in rows] != [
                float(one.split()[1]) for one in _ETHANE.splitlines()[2:]
                if one.strip()]
        finally:
            browser.close()
    assert not errors, errors
