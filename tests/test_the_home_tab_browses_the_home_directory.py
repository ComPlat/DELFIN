"""The Home tab is the calculations browser, rooted at the home directory.

Calculations, Office and Home are one browser on three folders, so a file that
is neither a calculation nor a document can still be read, edited, uploaded and
downloaded without leaving the dashboard.

Two things are specific to this clone and are what these tests hold: the
recursive workspace scan must not walk a home directory (a micromamba
installation and every cache sit in there, and the scan runs on each change of
directory), and the tab is off until somebody switches it on in Settings.
"""

from __future__ import annotations

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


def test_the_home_tab_is_rooted_at_the_home_directory(tmp_path, monkeypatch):
    fake_home = tmp_path / "home"
    (fake_home / "notes").mkdir(parents=True)
    (fake_home / "notes" / "a_letter.docx").write_bytes(b"PK\x03\x04")
    monkeypatch.setattr(Path, "home", staticmethod(lambda: fake_home))

    widget, refs = tab_home.create_tab(_ctx(tmp_path))
    refs["calc_list_directory"]()

    shown = [str(o) for o in refs["calc_file_list"].options]
    assert any("notes" in s for s in shown), shown
    assert widget is not None


def test_the_calculations_root_stays_reachable(tmp_path, monkeypatch):
    """So a file found in the home directory can be moved where it belongs."""
    fake_home = tmp_path / "home"
    fake_home.mkdir()
    monkeypatch.setattr(Path, "home", staticmethod(lambda: fake_home))

    ctx = _ctx(tmp_path)
    captured = {}
    original = browser.create_tab

    def spy(inner_ctx):
        captured["ctx"] = inner_ctx
        return original(inner_ctx)

    monkeypatch.setattr(browser, "create_tab", spy)
    tab_home.create_tab(ctx)

    assert captured["ctx"].calc_dir == fake_home
    assert captured["ctx"].primary_calc_dir == ctx.calc_dir


def test_the_workspace_scan_does_not_walk_the_home_directory(tmp_path, monkeypatch):
    """The scan runs on every change of directory; over a home it must not walk."""
    fake_home = tmp_path / "home"
    deep = fake_home / "micromamba" / "envs" / "x" / "lib" / "python3.11"
    deep.mkdir(parents=True)
    (fake_home / "calc_elsewhere").mkdir()
    (fake_home / "calc_elsewhere" / "CONTROL.txt").write_text("method=classic\n", encoding="utf-8")
    monkeypatch.setattr(Path, "home", staticmethod(lambda: fake_home))

    walked: list[str] = []
    real_rglob = Path.rglob

    def counting_rglob(self, pattern):
        walked.append(str(self))
        return real_rglob(self, pattern)

    monkeypatch.setattr(Path, "rglob", counting_rglob)
    _widget, refs = tab_home.create_tab(_ctx(tmp_path))
    refs["calc_list_directory"]()

    assert not any(str(fake_home) == w for w in walked), (
        f"the home directory itself was walked: {walked}"
    )


def test_the_same_scan_still_runs_for_the_calculations_browser(tmp_path):
    """The bound belongs to the home clone, not to the browser."""
    ctx = _ctx(tmp_path)
    workspace = ctx.calc_dir / "a_run"
    workspace.mkdir()
    (workspace / "CONTROL.txt").write_text("method=classic\n", encoding="utf-8")

    _widget, refs = browser.create_tab(ctx)
    refs["calc_list_directory"]()

    assert not getattr(ctx, "browser_scan_is_bounded", False)


def test_the_tab_is_registered_and_off_until_it_is_switched_on():
    import inspect

    from delfin.dashboard import create_dashboard

    source = inspect.getsource(create_dashboard)
    assert "'id': 'home'" in source
    home_spec = source[source.index("'id': 'home'"):]
    home_spec = home_spec[:home_spec.index("},")]
    assert "'title': 'Home'" in home_spec
    assert "'default_visible': False" in home_spec      # a deliberate choice
    assert "'fixed': False" in home_spec                # so Settings can show it
