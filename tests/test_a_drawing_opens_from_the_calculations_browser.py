"""A .ket is JSON and a .mol is a table of numbers.

Shown as text -- which is what the browser's default branch does with anything
it has no viewer for -- neither is something a chemist can read.  The one
program that can read them is already in the dashboard, so a drawing
double-clicked in the Calculations tab goes there and is editable.

The Ketcher tab is registered additively, and a registered tab that failed to
build is marked unavailable, so this has to work when there is no Ketcher tab
at all: the text preview is what there is then, and it is said rather than
shown silently.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.dashboard import ketcher
from delfin.dashboard.context import DashboardContext


@pytest.fixture
def dashboard(tmp_path, monkeypatch):
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_calculations_browser, tab_ketcher

    for name in ("calc", "archive", "office", "home"):
        (tmp_path / name).mkdir()
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(tmp_path))
    monkeypatch.setattr(pathlib.Path, "home", lambda: tmp_path / "home")
    folder = tmp_path / ".delfin" / "ketcher"
    folder.mkdir(parents=True)
    (folder / "index.html").write_text("<html></html>")
    (folder / ".delfin-ketcher-version").write_text("3.17.0")

    kept = tmp_path / "calc" / "Ketcher"
    kept.mkdir()
    (kept / "aspirin.ket").write_text('{"root":{"nodes":[]}}')
    (kept / "notes.txt").write_text("not a drawing")

    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    sent = []
    ctx.run_js = sent.append
    ctx.tabs_widget = type("Tabs", (), {"selected_index": None})()
    ctx.tab_indices = {"Ketcher": 7}

    _kt, ketcher_refs = tab_ketcher.create_tab(ctx)
    _cb, browser_refs = tab_calculations_browser.create_tab(ctx)
    browser_refs["sheet_state"]["current_path"] = "Ketcher"
    browser_refs["calc_list_directory"]()
    return {"ctx": ctx, "browser": browser_refs, "ketcher": ketcher_refs,
            "sent": sent, "tmp_path": tmp_path}


def test_a_drawing_is_marked_as_one_in_the_listing(dashboard):
    """A pencil, so it does not read as one more text file in a folder of
    output."""
    listed = list(dashboard["browser"]["calc_file_list"].options)

    assert "✏️ aspirin.ket" in listed
    assert any(name.endswith("notes.txt") and "✏" not in name
               for name in listed), "and only the drawings are marked"


def test_opening_it_puts_it_in_the_editor_and_goes_there(dashboard):
    dashboard["sent"].clear()

    dashboard["browser"]["calc_open_input"].value = "✏️ aspirin.ket"

    script = "\n".join(dashboard["sent"])
    assert "setMolecule" in script
    assert '{\\"root\\":{\\"nodes\\":[]}}' in script, "the drawing itself"
    assert dashboard["ctx"].tabs_widget.selected_index == 7, "and it is shown"
    assert "aspirin.ket" in dashboard["ketcher"]["ketcher_panel"].status.value


def test_it_is_not_shown_as_text(dashboard):
    """Which is what the default branch would have done with it."""
    dashboard["browser"]["calc_open_input"].value = "✏️ aspirin.ket"

    shown = dashboard["browser"]["calc_content_area"].value or ""
    assert '"root"' not in shown, "the JSON of their own structure"
    assert "aspirin.ket" in dashboard["browser"]["calc_file_info"].value


def test_it_can_be_reached_after_the_file_is_already_selected(dashboard):
    """Double-clicking opens it; this is for getting back to it without
    opening it again."""
    dashboard["browser"]["calc_update_options_dropdown"]("✏️ aspirin.ket")

    assert dashboard["browser"]["calc_options_dropdown"].options == (
        "(Options)", "Open in Ketcher")

    dashboard["sent"].clear()
    dashboard["browser"]["calc_file_list"].value = ("✏️ aspirin.ket",)
    dashboard["browser"]["calc_options_dropdown"].value = "Open in Ketcher"

    assert "setMolecule" in "\n".join(dashboard["sent"])


def test_a_dashboard_with_no_ketcher_tab_says_so_rather_than_doing_nothing(
        dashboard):
    """It is registered additively and can be hidden or fail to build.  A
    double-click that silently does nothing is worse than a text preview."""
    dashboard["ctx"].ketcher_refs = {}
    dashboard["sent"].clear()

    dashboard["browser"]["calc_open_input"].value = ""
    dashboard["browser"]["calc_open_input"].value = "✏️ aspirin.ket"

    shown = dashboard["browser"]["calc_content_area"].value or ""
    assert "setMolecule" not in "\n".join(dashboard["sent"])
    assert "no Ketcher tab" in shown, "a double-click that does nothing is worse"
    assert "nodes" in shown, "and the file itself is shown, as it was before"


def test_the_archive_and_office_tabs_find_it_too(tmp_path, monkeypatch):
    """They are the Calculations browser over a different root, made with
    ``dataclasses.replace`` -- a second context sharing this dict -- and they
    are built *before* the registered tabs are.  Rebinding the attribute would
    leave them holding the empty dict they were made with, and a drawing
    opened from Archive would say there was no Ketcher tab while it was
    sitting two tabs along."""
    pytest.importorskip("ipywidgets")
    import dataclasses

    from delfin.dashboard import tab_ketcher

    (tmp_path / "calc").mkdir()
    (tmp_path / "archive").mkdir()
    (tmp_path / "home").mkdir()
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(tmp_path))
    monkeypatch.setattr(pathlib.Path, "home", lambda: tmp_path / "home")
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive")
    ctx.run_js = lambda _script: None

    # The clone is made first, exactly as create_dashboard makes it.
    archive_ctx = dataclasses.replace(ctx, calc_dir=ctx.archive_dir)
    tab_ketcher.create_tab(ctx)

    assert callable(archive_ctx.ketcher_refs.get("open_drawing"))


def test_the_tab_publishes_itself_because_the_registry_keeps_only_the_widget(
        dashboard):
    """``registered_tab_specs`` takes ``result[0]`` and throws the refs away,
    so a tab that has to be reachable has to put itself somewhere reachable."""
    from delfin.dashboard import tab_registry

    assert "delfin.dashboard.tab_ketcher" in tab_registry._BUILTIN_DYNAMIC_TABS
    assert callable(dashboard["ctx"].ketcher_refs.get("open_drawing"))
    assert "ketcher" in [tab.id for tab in tab_registry.registered_tabs()]


def test_what_the_editor_opens_and_what_the_browser_offers_are_one_list():
    """Two lists would drift, and the one that drifted would be the one that
    made a file look unopenable."""
    for suffix in ketcher.DRAWING_SUFFIXES:
        assert ketcher.is_drawing(f"drawn{suffix}") is True
