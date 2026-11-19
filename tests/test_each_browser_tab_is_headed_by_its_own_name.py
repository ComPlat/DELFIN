"""Archive, Home and Office are headed by their own names.

The three are the calculations browser on another root. Only Office had a
heading of its own; Archive and Home were headed "📂 Calculations
Browser", so a reader in the Archive tab was told they were somewhere
else. Each clone now passes its tab's title, and the heading is that.
"""
from __future__ import annotations

from types import SimpleNamespace

import pytest

from delfin.dashboard import tab_archive_statistics, tab_home, tab_office
from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.context import DashboardContext


@pytest.mark.parametrize("title, heading", [
    ("Archive", "<h3>🗄️ Archive</h3>"),
    ("Home", "<h3>🏠 Home</h3>"),
    ("Office", "<h3>📁 Office</h3>"),
])
def test_a_titled_clone_is_headed_by_its_title(title, heading):
    ctx = SimpleNamespace(browser_title=title)
    assert browser._browser_heading(ctx, is_office_tab=False) == heading


def test_without_a_title_the_root_decides_as_before():
    assert (browser._browser_heading(SimpleNamespace(), is_office_tab=False)
            == "<h3>📂 Calculations Browser</h3>")
    assert (browser._browser_heading(SimpleNamespace(), is_office_tab=True)
            == "<h3>📁 Office</h3>")


def test_a_title_is_escaped():
    ctx = SimpleNamespace(browser_title="<b>x</b>")
    assert "<b>" not in browser._browser_heading(ctx, is_office_tab=False)


def test_each_clone_passes_its_own_title():
    """The wrappers say which tab they are; the main context does not."""
    import inspect
    for module, title in ((tab_archive_statistics, "Archive"),
                          (tab_home, "Home"), (tab_office, "Office")):
        assert f"browser_title = '{title}'" in inspect.getsource(module.create_tab)


@pytest.fixture
def ctx(tmp_path, monkeypatch):
    monkeypatch.setattr(browser, "load_remote_archive_enabled", lambda: False)
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    return DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )


def _heading(widget):
    for child in widget.children:
        value = getattr(child, "value", "")
        if isinstance(value, str) and value.startswith("<h3>"):
            return value
    raise AssertionError("no heading")


def test_the_archive_tab_is_headed_archive(ctx):
    archive, _ = tab_archive_statistics.create_tab(ctx)
    calc, _ = browser.create_tab(ctx)
    assert _heading(archive) == "<h3>🗄️ Archive</h3>"
    assert _heading(calc) == "<h3>📂 Calculations Browser</h3>"
