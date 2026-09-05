"""One reading, for the atoms that are picked now, and gone when they are not.

The reported symptom, on a saddle being examined in the editor: "das feld
davon sieht man manchmal mehrere und sie bleiben auch nachdem die markierungen
weg sind und auch manchmal uebereinander ... ich will die immer nur zu
ausgewaehlten atomen sehen."

The box itself was always right about what it said.  What was wrong is who
owned it.  It is built once per viewer and kept on that viewer's state, and a
re-render -- a new structure, a scan point, a frame the player lands on --
runs :js:func:`onViewerReady`, which reset the state and set ``measureBox`` to
null.  Null is not removal.  The element stayed in the corner of the picture
holding the last distance it had been given, owned by nothing, so nothing ever
hid it again when the picks went away; and the next pick built a second box on
top of it.  Enough re-renders, enough boxes.

Driven here rather than read off the source, because what is being claimed is
about the DOM after two renders and not about a line of JavaScript.
"""
from __future__ import annotations

from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
_BUNDLE = _REPO / "delfin" / "dashboard" / "static" / "3Dmol-min.js"

_ETHANE = """8
ethane
C  0.000000  0.000000  0.762900
C  0.000000  0.000000 -0.762900
H -0.505000  0.874000  1.162900
H -0.505000 -0.874000  1.162900
H  1.010000  0.000000  1.162900
H  0.505000  0.874000 -1.162900
H  0.505000 -0.874000 -1.162900
H -1.010000  0.000000 -1.162900
"""


def _blocks():
    """Every block the page gets, in the order the page gets them.

    Fewer is not a smaller test but a different one: the viewer factory is
    installed by the mouse patch, so a page without it has no picture to
    measure in.
    """
    from delfin.dashboard.molecule_forcefield_js import molecule_ff_bootstrap_js
    from delfin.dashboard.molecule_viewer import (
        RIGHT_MOUSE_TRANSLATE_PATCH_JS,
        measurement_bootstrap_js,
        structure_viewer_fullscreen_bootstrap_js,
        submit_manip_bootstrap_js,
        vendored_3dmol_js,
    )

    return [vendored_3dmol_js(), RIGHT_MOUSE_TRANSLATE_PATCH_JS,
            molecule_ff_bootstrap_js(), submit_manip_bootstrap_js(),
            measurement_bootstrap_js(),
            structure_viewer_fullscreen_bootstrap_js()]


def test_a_re_rendered_picture_does_not_collect_measure_boxes():
    """Two renders, two picks, and still exactly one box in the picture.

    Before the fix this ended with two, the older one showing the distance it
    had been given before the re-render and answering to nobody.
    """
    playwright = pytest.importorskip(
        "playwright.sync_api", reason="needs a browser to drive the viewer")
    if not _BUNDLE.is_file():
        pytest.skip("vendored 3Dmol bundle is not present")

    page_source = (
        "<div id='host' style='width:320px;height:240px;position:relative'>"
        "</div>")
    page_source += "".join(
        f'<script type="module">{source}</script>' for source in _blocks())

    errors: list[str] = []
    with playwright.sync_playwright() as p:
        try:
            browser = p.chromium.launch(args=["--use-gl=swiftshader"])
        except Exception as exc:
            pytest.skip(f"chromium unavailable: {exc}")
        page = browser.new_page(viewport={"width": 800, "height": 600})
        page.on("pageerror", lambda e: errors.append(str(e)))
        page.set_content(page_source)
        page.wait_for_timeout(2500)
        got = page.evaluate(
            """async (xyz) => {
                // The box is filled on the next drawn frame, not on the call
                // that picks the atoms: picking schedules a redraw and the
                // redraw is what writes the reading.  So each step here waits
                // for two frames, which is what the picture itself waits for.
                const drawn = () => new Promise((go) =>
                    requestAnimationFrame(() => requestAnimationFrame(go)));
                const host = document.getElementById('host');
                const view = window.__delfinCreateViewer(
                    host, {backgroundColor: 'white'});
                view.addModel(xyz, 'xyz');
                view.setStyle({}, {stick: {}});
                view.zoomTo();
                view.render();
                window._submitMolViewerByScope = {trial: view};
                const manip = window.__delfinSubmitManip;
                const count = () => host.querySelectorAll(
                    '.submit-manip-measure-box').length;
                const shown = () => Array.from(host.querySelectorAll(
                    '.submit-manip-measure-box')
                ).filter((b) => b.style.display !== 'none').length;
                manip.onViewerReady('trial', host);
                manip.setPicks('trial', [0, 1]);
                await drawn();
                const first = {boxes: count(), showing: shown()};
                // What a new structure, a scan point or a landed frame does.
                manip.onViewerReady('trial', host);
                manip.setPicks('trial', [0, 1]);
                await drawn();
                const second = {boxes: count(), showing: shown()};
                // And with nothing picked, nothing is on show.
                manip.setPicks('trial', []);
                await drawn();
                const empty = {boxes: count(), showing: shown()};
                return {first, second, empty};
            }""",
            _ETHANE,
        )
        browser.close()

    assert not errors, f"the viewer threw: {errors[:3]}"
    assert got["first"] == {"boxes": 1, "showing": 1}, got
    # The one that was reported: a second render used to leave the first box
    # standing and build another beside it.
    assert got["second"] == {"boxes": 1, "showing": 1}, got
    # And a picture with nothing picked shows no reading at all.
    assert got["empty"]["showing"] == 0, got
