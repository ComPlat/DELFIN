"""The download control was dead for a folder, and silent about it.

A field report says a download on a folder "im Root" does not respond, or
only sporadically, and that no message says whether it is forbidden.

There is no permission rule involved. A disabled ipywidgets Button fires
no click, so ``calc_on_download`` -- and every message it owns -- never
ran. The button's state was recomputed in only two places: at the end of
a directory listing, and on the path that opens a FILE. Selecting a folder
returns early on both open paths, so at the explorer root a folder left
the button grey unless a file had been clicked first in the same listing.
That ordering is the "sporadically".

Three changes: every selection change re-evaluates the control, a control
that is genuinely disabled carries the reason, and a folder is measured
before its archive is built instead of after.

Universal: one browser instance over ``tmp_path``, no host paths.
"""

from __future__ import annotations

import pytest

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.context import DashboardContext


@pytest.fixture
def refs(tmp_path):
    """A browser at a root holding one folder and one file."""
    root = tmp_path / "calc"
    (root / "a_folder").mkdir(parents=True)
    (root / "a_folder" / "inner.txt").write_text("x\n", encoding="utf-8")
    (root / "note.txt").write_text("hello\n", encoding="utf-8")
    ctx = DashboardContext(calc_dir=root, archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    _widget, refs = browser.create_tab(ctx)
    refs["calc_list_directory"]()
    return refs


def _option_for(refs, name):
    for opt in refs["calc_file_list"].options:
        value = opt[1] if isinstance(opt, tuple) else opt
        if name in str(value) or name in str(opt):
            return value
    raise AssertionError(f"{name} not in {refs['calc_file_list'].options}")


def _select(refs, name):
    """Select one entry the way a click does, observers included."""
    value = _option_for(refs, name)
    refs["calc_file_list"].value = (value,)
    refs["calc_on_selection_change"]({"new": (value,)})


# ---------------------------------------------------------------------------
# The control follows the selection, whatever kind it is
# ---------------------------------------------------------------------------

def test_selecting_a_folder_at_the_root_enables_the_download(refs):
    _select(refs, "a_folder")
    btn = refs["calc_download_btn"]
    assert btn.disabled is False
    assert btn.description == "Download ZIP"


def test_selecting_a_file_at_the_root_enables_the_download(refs):
    _select(refs, "note.txt")
    btn = refs["calc_download_btn"]
    assert btn.disabled is False
    assert btn.description == "Download"


def test_the_folder_works_without_a_file_clicked_first(refs):
    """The reported ordering. A folder as the FIRST action used to leave the
    button grey; only a prior file click made it work."""
    _select(refs, "a_folder")
    assert refs["calc_download_btn"].disabled is False


def test_the_folder_still_works_after_a_file_was_clicked(refs):
    """The ordering that already worked must keep working."""
    _select(refs, "note.txt")
    _select(refs, "a_folder")
    btn = refs["calc_download_btn"]
    assert btn.disabled is False
    assert btn.description == "Download ZIP"


def test_going_back_to_a_file_switches_the_label_back(refs):
    _select(refs, "a_folder")
    _select(refs, "note.txt")
    assert refs["calc_download_btn"].description == "Download"


# ---------------------------------------------------------------------------
# A disabled control is never silent
# ---------------------------------------------------------------------------

def test_a_disabled_download_says_why(refs):
    """At the root with nothing selected there is no target. The button
    stays grey -- but a reader can now find out why without guessing that
    it is a refusal."""
    refs["calc_file_list"].value = ()
    refs["calc_update_download_btn"]()
    btn = refs["calc_download_btn"]
    assert btn.disabled is True
    assert btn.tooltip
    assert "select" in btn.tooltip.lower()


def test_an_enabled_download_names_its_target(refs):
    _select(refs, "a_folder")
    assert "a_folder" in refs["calc_download_btn"].tooltip


def test_a_click_with_no_target_writes_a_message(refs):
    """Belt and braces: if the control is ever reachable with no target,
    the handler says so rather than returning in silence."""
    refs["calc_file_list"].value = ()
    refs["calc_download_status"].value = ""
    refs["calc_on_download"](None)
    assert refs["calc_download_status"].value


# ---------------------------------------------------------------------------
# A folder is measured before it is built
# ---------------------------------------------------------------------------

def test_a_download_reports_that_it_started(refs):
    _select(refs, "a_folder")
    refs["calc_on_download"](None)
    status = refs["calc_download_status"].value
    assert "Download started" in status
    assert ".zip" in status


def test_a_folder_over_the_limit_is_refused_before_it_is_zipped(tmp_path,
                                                                monkeypatch):
    """The limit used to be applied to the finished payload, so the kernel
    built the whole archive first and then threw it away -- the browser went
    unresponsive with no message.

    The oversized file is sparse: ``st_size`` reports 40 MB while almost no
    disk is used, so the test is fast and needs nothing from the host.
    """
    import shutil as _shutil

    root = tmp_path / "calc"
    big = root / "big_folder"
    big.mkdir(parents=True)
    with open(big / "sparse.bin", "wb") as fh:
        fh.truncate(40 * 1024 * 1024)

    built = []
    real_make_archive = _shutil.make_archive
    monkeypatch.setattr(_shutil, "make_archive",
                        lambda *a, **k: (built.append(a),
                                         real_make_archive(*a, **k))[1])

    ctx = DashboardContext(calc_dir=root, archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    _widget, refs = browser.create_tab(ctx)
    refs["calc_list_directory"]()
    value = _option_for(refs, "big_folder")
    refs["calc_file_list"].value = (value,)
    refs["calc_on_selection_change"]({"new": (value,)})
    assert refs["calc_download_btn"].disabled is False

    refs["calc_on_download"](None)
    assert built == [], "the archive was built despite the size being known"
    status = refs["calc_download_status"].value
    assert "too large" in status.lower()
    assert "before compression" in status


def test_the_size_estimate_counts_the_files_under_a_folder(tmp_path):
    root = tmp_path / "t"
    (root / "sub").mkdir(parents=True)
    (root / "a.bin").write_bytes(b"x" * 10)
    (root / "sub" / "b.bin").write_bytes(b"y" * 5)
    ctx = DashboardContext(calc_dir=tmp_path, archive_dir=tmp_path,
                           office_dir=tmp_path)
    ctx.run_js = lambda script: None
    _widget, refs = browser.create_tab(ctx)
    assert refs["_calc_tree_size"](root) == 15


def test_the_estimate_stops_once_the_limit_is_passed(tmp_path):
    """A home directory is not worth walking to the end to learn that it
    is too big."""
    root = tmp_path / "t"
    root.mkdir()
    for i in range(20):
        (root / f"f{i}.bin").write_bytes(b"z" * 100)
    ctx = DashboardContext(calc_dir=tmp_path, archive_dir=tmp_path,
                           office_dir=tmp_path)
    ctx.run_js = lambda script: None
    _widget, refs = browser.create_tab(ctx)
    assert refs["_calc_tree_size"](root, cap=250) <= 2000
    assert refs["_calc_tree_size"](root, cap=250) > 250


def test_the_estimate_does_not_raise_on_a_broken_link(tmp_path):
    root = tmp_path / "t"
    root.mkdir()
    (root / "good.bin").write_bytes(b"x" * 3)
    (root / "dangling").symlink_to(tmp_path / "nothing-here")
    ctx = DashboardContext(calc_dir=tmp_path, archive_dir=tmp_path,
                           office_dir=tmp_path)
    ctx.run_js = lambda script: None
    _widget, refs = browser.create_tab(ctx)
    assert refs["_calc_tree_size"](root) == 3
