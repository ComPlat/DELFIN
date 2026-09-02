"""Every hotkey went to the dashboard, and Paste did nothing at all.

Ketcher binds ``keydown``, and ``copy``, ``cut`` and ``paste``, on the document
*inside* its frame.  Those handlers run only while that document holds the
focus.  Nothing in DELFIN ever gave it to them -- the frame carried no
tabindex, no ``.focus()`` was ever called on it -- so pressing ``t`` over the
canvas typed into the dashboard, and ``navigator.clipboard.read()``, which the
bundle calls in three places, threw *"Document is not focused"*.

One cause with two names: the keyboard did nothing and Paste did nothing.  The
same fix answers both, and it is also what makes a structure pasted out of
ChemDraw arrive -- the bundle already reads ``chemical/x-cdx`` and
``chemical/x-cdxml``; it simply never saw the paste.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.dashboard import ketcher
from delfin.dashboard.context import DashboardContext
from editor_source import SUBMIT_SOURCE


# ---------------------------------------------------------------------------
# the frame
# ---------------------------------------------------------------------------
def test_the_frame_can_hold_the_focus_and_carries_the_clipboard():
    frame = ketcher.frame_html("/voila/files/.delfin/ketcher/index.html?v=3.17.0")

    assert "tabindex='0'" in frame, "a frame with no tabindex cannot be focused"
    assert "allow='clipboard-read; clipboard-write'" in frame
    assert "/voila/files/.delfin/ketcher/index.html?v=3.17.0" in frame


def test_the_frame_the_structure_editor_shows_is_that_one():
    """One implementation, so a fix to it is a fix everywhere it is shown."""
    maker = SUBMIT_SOURCE.split("def _draw_frame_html")[1].split("\n    def ")[0]

    assert "_ketcher.frame_html" in maker
    assert "<iframe" not in maker, "no second copy of the frame to fall behind"


def test_the_frame_is_as_tall_as_it_was():
    """The panel is not the place to make the editor bigger: the tab is."""
    assert "height:560px" in ketcher.frame_html("x", height="560px")
    assert "height:72vh" in ketcher.frame_html("x", height="72vh")


# ---------------------------------------------------------------------------
# the focus
# ---------------------------------------------------------------------------
def test_the_focus_is_handed_over_when_the_pointer_is_on_the_editor():
    js = ketcher.focus_js(".submit-ketcher-frame.scope-1")

    assert "contentWindow.focus()" in js
    assert "addEventListener('pointerenter',reach)" in js
    assert "addEventListener('pointerdown',reach,true)" in js
    assert ".submit-ketcher-frame.scope-1" in js


def test_the_focus_is_never_taken_by_a_timer_or_on_load():
    """A frame that takes the focus by itself a couple of seconds after it
    appears is what scrolls the page away under the user -- measured at +3 s,
    pane 0 to 654, which is what the scroll hold exists to undo.  And taking
    it while somebody is typing in a field of the dashboard would be a worse
    bug than the one being fixed.

    The timer that is here looks for the *element*: a widget's value is set
    from the kernel and the DOM catches up afterwards.  It calls ``bind``,
    never ``reach``.
    """
    js = ketcher.focus_js(".frame")

    assert "setTimeout(bind,150)" in js
    assert "setTimeout(reach" not in js
    assert "addEventListener('load'" not in js
    for line in js.splitlines():
        if "setTimeout" in line or "setInterval" in line:
            assert "reach" not in line, line


def test_binding_twice_binds_once():
    """It travels with every question asked of the editor, so it arrives many
    times over a session."""
    js = ketcher.focus_js(".frame")

    assert "__delfinKetcherFocus" in js
    assert "if(host.__delfinKetcherFocus) return;" in js


def test_it_is_sent_with_the_scroll_hold_rather_than_after_it():
    """``run_js`` clears its output before displaying the next script, so two
    calls in a row can mean the first is thrown away before the browser has
    run it.  That is why the answer channel is a widget value, and it is why
    these two travel joined."""
    opener = SUBMIT_SOURCE.split("def on_submit_draw_open")[1].split("\n    def ")[0]

    assert "_run_manip_js(_KETCHER_SCROLL_HOLD_JS + _KETCHER_FOCUS_JS)" in opener
    assert opener.count("_run_manip_js(") == 1, "one script, not two"


# ---------------------------------------------------------------------------
# in the tab
# ---------------------------------------------------------------------------
@pytest.fixture
def editor(tmp_path, monkeypatch):
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_submit

    for name in ("calc", "archive", "office", "home"):
        (tmp_path / name).mkdir()
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(tmp_path))
    monkeypatch.setattr(pathlib.Path, "home", lambda: tmp_path / "home")
    folder = tmp_path / ".delfin" / "ketcher"
    folder.mkdir(parents=True)
    (folder / "index.html").write_text("<html></html>")
    (folder / ".delfin-ketcher-version").write_text("3.17.0")
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    sent = []
    ctx.run_js = sent.append
    _widget, refs = tab_submit.create_tab(ctx)
    refs["sent_js"] = sent
    return refs


def test_opening_the_editor_hands_it_the_keyboard(editor):
    editor["sent_js"].clear()

    editor["submit_draw_open_btn"].value = True

    frame = editor["submit_draw_frame"].value
    assert "tabindex='0'" in frame
    assert "clipboard-read" in frame
    script = "\n".join(editor["sent_js"])
    assert "contentWindow.focus()" in script
    assert "pointerenter" in script


def test_each_editor_focuses_its_own_frame(editor, tmp_path, monkeypatch):
    """Both tabs keep their drawing frame outside the scope container, so the
    class on the element is the only way to tell one editor's frame from the
    other's -- which is how TO SMILES once answered into the wrong tab."""
    editor["submit_draw_open_btn"].value = True
    scope = editor["submit_scope_id"]

    script = "\n".join(editor["sent_js"])
    assert f".submit-ketcher-frame.{scope}" in script
    assert ".submit-ketcher-frame'" not in script, "never 'the' frame"
