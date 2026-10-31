"""A chat rebuilt while an answer streams keeps its height.

Input: the transcript HTML, which Python replaces as a whole on every
streaming update, and its 3D cards. Before (measured in the real page,
2026-10-22): every 3D card came back with a fixed 380 px stage and was
sized by the browser script to 80% of the chat frame (619 px at
1920x1080) only after 3Dmol had built it. With three cards above the
end, the transcript shrank by 3 x 239 = 717 px and grew back on every
update, about every 1.5 s while an answer streamed: a reader at the end
saw the chat jump up and back.
After: the script keeps the stage height it computed on the chat frame
(``--delfin-mol3d-h``), which is not replaced with the transcript, and a
card's stage starts at that height (380 px before any card was built).
"""
from __future__ import annotations

from pathlib import Path


def test_a_new_stage_starts_at_the_height_the_last_one_was_given(tmp_path):
    from delfin.dashboard import tab_agent
    xyz = tmp_path / "w.xyz"
    xyz.write_text("3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n")
    html = tab_agent._mol3d_card_html(Path(xyz), "w.xyz")
    stage = html.split('class="delfin-mol3d-stage"', 1)[1].split(">", 1)[0]
    assert "height:var(--delfin-mol3d-h, 380px)" in stage, stage


def test_the_script_keeps_the_height_on_the_frame():
    from delfin.dashboard import tab_agent
    src = Path(tab_agent.__file__).read_text(encoding="utf-8")
    build = src[src.index("function build(card, stage)"):]
    build = build[:build.index("$3Dmol.createViewer")]
    assert "setProperty('--delfin-mol3d-h'" in build
