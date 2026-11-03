"""A rebuilt chat keeps its 3D viewers instead of building new ones.

Input: the chat's 3D card script and a transcript Python replaces as a
whole on every update. Measured in the real page (2026-10-22): a session
with three 3D cards created 15 canvases in 40 s -- five transcript
updates, every viewer built again each time, each with a new WebGL
context. Chrome keeps about 16 contexts and drops the oldest; over a
55-minute field session with many cards (report 18:42, KIT cluster) the
page froze and was reloaded.

After: a viewer that was built is kept (by its stage id and a hash over
the whole coordinate text), and a rebuilt card that shows
the same molecule gets that viewer's stage moved into it -- canvas,
WebGL context and frame included -- with its controls wired to it again.
A card whose molecule changed is built anew. A kept viewer that no card
claims within ten seconds gives its WebGL context back, and at most 12
detached ones are held, below Chrome's limit. The real-page measurement
is in the commit that makes the change; this pins the mechanism in the
shipped script.
"""
from __future__ import annotations

from pathlib import Path


def _script() -> str:
    from delfin.dashboard import tab_agent
    src = Path(tab_agent.__file__).read_text(encoding="utf-8")
    # The 3D card script, from its install guard to its entry point's end.
    start = src.index("window.__delfinMol3DInstalled = true;")
    end = src.index("window.__delfinMol3D = function(el)", start)
    return src[start:end + 6000]


def test_a_built_viewer_is_kept_by_stage_and_coordinates():
    js = _script()
    assert "__delfinMol3DKept" in js
    assert "stage.id" in js and "data-mol3d" in js
    # The whole text decides, not a prefix: two conformers share the first
    # lines of their coordinates.
    assert "hashOf(d)" in js and "d.slice(" not in js


def test_an_unclaimed_viewer_gives_its_context_back():
    js = _script()
    assert "KEEP_MAX = 12" in js
    assert "WEBGL_lose_context" in js and "loseContext()" in js
    assert "release(k)" in js


def test_a_kept_viewer_is_moved_into_the_rebuilt_card():
    js = _script()
    assert "replaceChild(kept.stage, stage)" in js
    assert "wireTrajectory(card, kept.stage, kept.viewer" in js
    assert "wireFullscreen(card, kept.viewer" in js
