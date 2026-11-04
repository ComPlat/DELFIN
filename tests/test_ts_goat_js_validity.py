"""The generated TS->GOAT viewer JavaScript must be syntactically
valid.

Regression for the P0 white-viewer bug of 2026-11-02: a stick-restyle
removal left dangling ``viewer.setStyle({serial: serial},`` calls and
orphan ``bt`` references inside three JS generators, so the whole
``(function(){...})()`` failed to parse -- the viewer stayed blank and
no panel message appeared, because the syntax error only reached the
browser console.  Every ts_goat JS generator is now checked with
``node --check``; skipped when node is unavailable.
"""

import shutil
import subprocess

import pytest

node = shutil.which("node")

pytestmark = pytest.mark.skipif(
    node is None, reason="node is not installed on this host")


def _check(js: str, name: str) -> None:
    import tempfile
    from pathlib import Path

    with tempfile.NamedTemporaryFile(
            "w", suffix=".js", delete=False) as fh:
        fh.write(js)
        path = fh.name
    try:
        result = subprocess.run(
            [node, "--check", path], capture_output=True, text=True)
    finally:
        Path(path).unlink(missing_ok=True)
    assert result.returncode == 0, (
        f"{name} produced syntactically invalid JavaScript: "
        f"{result.stderr.strip()[:400]}")


def test_bootstrap_js_is_valid():
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_bootstrap_js
    _check(
        ts_goat_viewer_bootstrap_js(
            "vid-test", arrows=[0.1, 0.2, 0.3] * 3),
        "ts_goat_viewer_bootstrap_js")


def test_step_js_is_valid():
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_step_js
    _check(ts_goat_viewer_step_js("vid-test", 1.0), "ts_goat_viewer_step_js")


def test_frames_js_is_valid():
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_frames_js
    _check(
        ts_goat_viewer_frames_js("vid-test", running=True),
        "ts_goat_viewer_frames_js")


def test_center_js_is_valid():
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_center_js
    _check(ts_goat_viewer_center_js("vid-test"), "ts_goat_viewer_center_js")


def test_bootstrap_without_arrows_is_valid():
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_bootstrap_js
    _check(
        ts_goat_viewer_bootstrap_js("vid-test", arrows=None),
        "ts_goat_viewer_bootstrap_js(no arrows)")



def test_highlight_js_is_valid():
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_highlights_js
    _check(
        ts_goat_viewer_highlights_js("vid-test", [0, 1, 8]),
        "ts_goat_viewer_highlights_js")
