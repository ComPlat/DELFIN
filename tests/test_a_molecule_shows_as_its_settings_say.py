"""A molecule card shows as the viewer settings say, folded turn or not.

Input: the agent chat (delfin/dashboard/tab_agent.py). Three defects:
a finished turn folded its molecule card away with its tool steps; the
inline card drew a fixed ball-and-stick instead of the global viewer
style; the fullscreen view called setStyle before adding the model, so
3Dmol styled no atoms and drew its default thin lines. The browser
script is checked as text: no node is needed to run these tests.
"""
from __future__ import annotations

import pathlib
import re

from delfin.dashboard import tab_agent as T

_SRC = (pathlib.Path(__file__).resolve().parents[1]
        / "delfin" / "dashboard" / "tab_agent.py").read_text()

_CARD = ('<div class="delfin-artifact delfin-chat-mol3d" '
         'data-mol3d="x">mol</div>')


def _render(m):
    return f"<p>{m['content']}</p>"


def test_a_folded_turn_keeps_its_molecule_card_visible():
    msgs = [
        {"role": "user", "content": "show it"},
        {"role": "tool", "content": '<span class="tool-name">Read</span>'},
        {"role": "tool", "content": _CARD},
        {"role": "tool", "content": '<span class="tool-name">Bash</span>'},
        {"role": "assistant", "content": "the answer",
         "turn_fold": {"duration_s": 3, "outcome": "ok"}},
    ]
    parts = T._render_transcript(msgs, _render)
    fold = [p for p in parts if "delfin-chat-turn-steps" in p]
    assert len(fold) == 1
    assert "tool-name\">Read" in fold[0] and "tool-name\">Bash" in fold[0]
    assert "delfin-chat-mol3d" not in fold[0]
    top = [p for p in parts if "delfin-chat-mol3d" in p]
    assert len(top) == 1
    assert parts.index(top[0]) == parts.index(fold[0]) + 1
    assert parts[-1] == "<p>the answer</p>"


def test_a_turn_of_only_display_results_draws_no_empty_fold():
    msgs = [
        {"role": "user", "content": "show it"},
        {"role": "tool", "content": _CARD},
        {"role": "assistant", "content": "the answer",
         "turn_fold": {"duration_s": 3, "outcome": "ok"}},
    ]
    parts = T._render_transcript(msgs, _render)
    assert not any("delfin-chat-turn-steps" in p for p in parts)
    assert any("delfin-chat-mol3d" in p for p in parts)


def _script():
    a = _SRC.index("function build(card, stage)")
    b = _SRC.index("function closeFs()", a)
    return _SRC[a:b]


def _fullscreen():
    s = _script()
    return s[s.index("function wireFullscreen("):]


def test_the_inline_card_draws_the_configured_style():
    build = _script()[:_script().index("function wireTrajectory(")]
    assert "{stick: {radius: 0.12}" not in build
    assert build.count("setStyle({}, look.style)") == 3
    assert "createViewer(stage, look.config)" in build


def test_the_fullscreen_viewer_uses_the_same_look():
    fs = _fullscreen()
    assert "createViewer(big, look.config)" in fs
    assert "{backgroundColor: 'white'}" not in fs


def test_the_fullscreen_style_is_set_after_the_model():
    fs = _fullscreen()
    # Every branch adds its model, then styles it; a setStyle before the
    # first addModel styles an empty viewer (3Dmol's default lines).
    first_add = re.search(r"v2\.addModel", fs).start()
    first_style = fs.index("v2.setStyle(")
    assert first_add < first_style
    assert fs.count("v2.setStyle({}, look.style)") == 3


def test_the_card_carries_the_global_viewer_config(monkeypatch, tmp_path):
    from delfin.dashboard import molecule_viewer as MV
    prof = dict(MV.get_viewer_profile())
    prof["viewer_config_js"] = '{"backgroundColor":"#123456"}'
    monkeypatch.setattr(MV, "get_viewer_profile", lambda: prof)
    xyz = tmp_path / "w.xyz"
    xyz.write_text("3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n")
    card = T._mol3d_card_html(xyz, "w.xyz")
    assert card is not None
    assert 'data-mol3d-config="{&quot;backgroundColor&quot;:&quot;#123456&quot;}"' in card
