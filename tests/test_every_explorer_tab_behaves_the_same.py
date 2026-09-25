"""Calculations, Archive, Office and Home are one browser, so they behave alike.

Viewing a document was the Office tab's privilege: only there did a PDF scroll
through instead of turning one page at a time, only there did a Word document
render at all, only there did the frame take the height of the pane, and only
there did picking a row leave the viewport where it was. The other three tabs
are the same browser on another folder, so the difference was an accident of
which flag the code read.

What still differs is writing: the Backups folder belongs beside a document,
not inside a calculation.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard import tab_home
from delfin.dashboard.context import DashboardContext


def _ctx(tmp_path: Path) -> DashboardContext:
    for folder in ("calc", "archive", "office"):
        (tmp_path / folder).mkdir(exist_ok=True)
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    ctx.run_js = lambda script: None
    return ctx


def _listing(refs) -> list[str]:
    refs["calc_list_directory"]()
    return [str(o) for o in refs["calc_file_list"].options]


# -- one behaviour in every clone -----------------------------------------

def test_a_document_is_viewed_the_same_way_everywhere():
    """The flags the viewer reads must not depend on which folder it sits on."""
    source = inspect.getsource(browser.create_tab)
    assert "_DOC_VIEW_FEEL = True" in source

    for site in ("continuous=_DOC_VIEW_FEEL",                 # a PDF scrolls through
                 "height_px=None if _DOC_VIEW_FEEL",           # and takes the pane's height
                 "office=_DOC_VIEW_FEEL",                      # a sheet is a grid
                 "suffix == '.docx' and _DOC_VIEW_FEEL",       # a Word file renders
                 "[calc_fullscreen_btn] if _DOC_VIEW_FEEL"):   # and can fill the window
        assert site in source, site


def test_writing_still_belongs_to_the_office_tab():
    """A Backups folder beside a document is right; inside a calculation it is not."""
    source = inspect.getsource(browser.create_tab)
    assert "OFFICE_BACKUP_DIR if _OFFICE_DOC_FEEL else None" in source
    assert "versioned=_OFFICE_DOC_FEEL" in source


def test_picking_a_row_does_not_drop_the_user_at_the_far_end():
    """Selecting a row anchored at the far end, so the viewport holds still.

    Outside the office feel the selection ran to the last column and the view
    was revealed there -- pick row 12 and land at the end of the sheet. The
    grid decides that from data-office, which the browser now always sets.
    """
    from delfin.dashboard import spreadsheet_view as sheet

    source = inspect.getsource(sheet)
    assert "var OFFICE = wrap.dataset.office === '1';" in source
    assert "anchor = {r: ri, c: colCount()};" in source      # the far end
    assert "moveTo(ri, 1, true, true);" in source            # and the view holds

    # and the flag reaches the markup the grid reads it from
    assert 'data-office="{"1" if office else "0"}"' in source


# -- the dotfile toggle ----------------------------------------------------

def test_the_calculations_browser_shows_everything_and_offers_no_toggle(tmp_path):
    """A calculations folder is not mostly dot folders, so there is nothing to hide."""
    ctx = _ctx(tmp_path)
    (ctx.calc_dir / "a_run").mkdir()
    (ctx.calc_dir / ".delfin_last_run.json").write_text("{}", encoding="utf-8")

    _widget, refs = browser.create_tab(ctx)

    assert any(".delfin_last_run.json" in o for o in _listing(refs))
    row = refs["calc_nav_selection_row"] if "calc_nav_selection_row" in refs else None
    if row is not None:
        assert refs["calc_dotfiles_btn"] not in row.children


def test_in_the_home_tab_the_toggle_hides_and_shows_again(tmp_path, monkeypatch):
    fake_home = tmp_path / "home"
    (fake_home / ".cache").mkdir(parents=True)
    (fake_home / "messwerte").mkdir()
    monkeypatch.setattr(Path, "home", staticmethod(lambda: fake_home))

    _widget, refs = tab_home.create_tab(_ctx(tmp_path))
    button = refs["calc_dotfiles_btn"]

    assert button.value is False
    assert not any(".cache" in o for o in _listing(refs))

    button.value = True
    assert any(".cache" in o for o in _listing(refs))
    assert any("messwerte" in o for o in _listing(refs))           # the rest is still there

    button.value = False
    assert not any(".cache" in o for o in _listing(refs))


def test_the_home_tab_starts_without_the_dot_folders(tmp_path, monkeypatch):
    fake_home = tmp_path / "home"
    (fake_home / ".cache").mkdir(parents=True)
    (fake_home / "messwerte").mkdir()
    monkeypatch.setattr(Path, "home", staticmethod(lambda: fake_home))

    _widget, refs = tab_home.create_tab(_ctx(tmp_path))

    assert refs["calc_dotfiles_btn"].value is False
    shown = _listing(refs)
    assert any("messwerte" in o for o in shown)
    assert not any(".cache" in o for o in shown)


# -- icons -----------------------------------------------------------------

@pytest.mark.parametrize("name, icon", [
    ("a_paper.pdf", "📕"),
    ("a_table.xlsx", "📊"),
    ("a_table.csv", "📊"),
    ("a_letter.docx", "📃"),
    ("a_picture.jpg", "🖼"),
    ("an_archive.zip", "🗜"),
    ("a_talk.pptx", "📽"),
    ("a_structure.xyz", "🔬"),
])
def test_a_file_type_is_recognisable_by_its_icon(tmp_path, name, icon):
    ctx = _ctx(tmp_path)
    (ctx.calc_dir / name).write_bytes(b"x")

    _widget, refs = browser.create_tab(ctx)
    row = [o for o in _listing(refs) if name in o]

    assert row, f"{name} was not listed"
    assert row[0].startswith(icon), row[0]


# -- the splitter ----------------------------------------------------------

def test_the_splitter_is_dragged_with_a_pointer_not_only_a_mouse():
    """A finger and a pen produce pointer events; mouse events they do not."""
    source = inspect.getsource(browser.create_tab)
    assert "addEventListener('pointerdown'" in source
    assert "setPointerCapture" in source           # the drag survives leaving the strip
    assert "touch-action:none" in source           # and a finger drags instead of scrolling
    assert "addEventListener('mousedown', function(e) {{\n                e.preventDefault();\n                document.addEventListener('mousemove'" not in source
