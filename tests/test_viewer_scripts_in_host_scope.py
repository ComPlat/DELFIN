"""Every shipped viewer script must survive the host's evaluation scope.

``ctx.run_js`` hands scripts to JupyterLab/Voila through
``display(Javascript(...))``, which evaluates them with ``this`` undefined and
without giving top-level ``var`` declarations global scope. A ``<script>`` tag
does neither of those things, so a harness built on one cannot see this class
of failure -- and did not: the vendored 3Dmol bundle begins with ``root = this``
and ends by assigning ``root["3Dmol"]``, so it threw

    Cannot set properties of undefined (setting '3Dmol')

which aborted the rest of the start-up script and left every viewer empty.

This test loads each block the way the host does, in a module scope, and
checks that the globals the dashboard depends on actually exist afterwards.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
_BUNDLE = _REPO / "delfin" / "dashboard" / "static" / "3Dmol-min.js"

# Globals other code calls by name. A block that throws part way through
# silently takes the ones defined below its failure point with it.
_REQUIRED_GLOBALS = (
    "$3Dmol",
    "$3Dmolpromise",
    "__delfinCreateViewer",
    "__delfinCoalesce",
    "__delfinOnViewerCreated",
    "__delfinDisposeViewer",
    "__delfinRenderViewerPng",
    "__delfinFF",
    "__delfinSubmitManip",
    "__delfinMeasure",
)

_XYZ = """5
methane
C  0.000000  0.000000  0.000000
H  0.629118  0.629118  0.629118
H -0.629118 -0.629118  0.629118
H -0.629118  0.629118 -0.629118
H  0.629118 -0.629118 -0.629118
"""


def _script_blocks():
    from delfin.dashboard.molecule_forcefield_js import molecule_ff_bootstrap_js
    from delfin.dashboard.molecule_viewer import (
        RIGHT_MOUSE_TRANSLATE_PATCH_JS,
        measurement_bootstrap_js,
        structure_viewer_fullscreen_bootstrap_js,
        submit_manip_bootstrap_js,
        vendored_3dmol_js,
    )

    return [
        ("vendored 3Dmol", vendored_3dmol_js()),
        ("mouse patch", RIGHT_MOUSE_TRANSLATE_PATCH_JS),
        ("force field", molecule_ff_bootstrap_js()),
        ("editor", submit_manip_bootstrap_js()),
        ("measurement", measurement_bootstrap_js()),
        ("fullscreen", structure_viewer_fullscreen_bootstrap_js()),
    ]


def test_every_shipped_script_survives_a_module_scope():
    playwright = pytest.importorskip(
        "playwright.sync_api", reason="needs a browser to evaluate the scripts"
    )
    if not _BUNDLE.is_file():
        pytest.skip("vendored 3Dmol bundle is not present")

    blocks = _script_blocks()
    # type="module" reproduces the host's scope: `this` is undefined and a
    # top-level `var` does not become a global.
    page_source = "<div id='host' style='width:320px;height:240px;position:relative'></div>"
    page_source += "".join(
        f'<script type="module">{source}\n'
        f'window.__delfinBlocksRun=(window.__delfinBlocksRun||0)+1;</script>'
        for _name, source in blocks
    )

    errors: list[str] = []
    with playwright.sync_playwright() as p:
        try:
            browser = p.chromium.launch(args=["--use-gl=swiftshader"])
        except Exception as exc:  # no browser installed in this environment
            pytest.skip(f"chromium unavailable: {exc}")
        page = browser.new_page(viewport={"width": 800, "height": 600})
        page.on("pageerror", lambda e: errors.append(str(e)))
        page.set_content(page_source)
        page.wait_for_timeout(2500)

        finished = page.evaluate("() => window.__delfinBlocksRun || 0")
        present = page.evaluate(
            "(names) => names.map((n) => [n, typeof window[n]])",
            list(_REQUIRED_GLOBALS),
        )
        built = page.evaluate(
            """(xyz) => {
                try {
                    window.__delfinViewerPixelRatio = 2;
                    const v = window.__delfinCreateViewer(
                        document.getElementById('host'), {backgroundColor: 'white'});
                    v.addModel(xyz, 'xyz');
                    v.setStyle({}, {stick: {}});
                    v.zoomTo();
                    v.render();
                    return {ok: true,
                            atoms: v.getModel().selectedAtoms({}).length};
                } catch (e) {
                    return {ok: false, error: String(e)};
                }
            }""",
            _XYZ,
        )
        browser.close()

    assert not errors, f"scripts threw in a module scope: {errors[:3]}"
    assert finished == len(blocks), (
        f"only {finished} of {len(blocks)} script blocks ran to completion"
    )
    missing = [name for name, kind in present if kind == "undefined"]
    assert not missing, f"globals missing after load: {missing}"
    assert built.get("ok"), f"no viewer could be built: {built.get('error')}"
    assert built["atoms"] == 5, json.dumps(built)
