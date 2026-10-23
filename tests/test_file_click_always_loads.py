"""Clicking a file loads it -- also when it is the file that is already selected.

The browser's file list reported *selection changes*, and a click on the entry
that is already selected changes nothing.  So the click produced no message, the
file did not load, and the way out was to deselect it and select it again.  The
list now reports the click itself, with a counter, so no click is dropped as
"same value".
"""

from __future__ import annotations

import pytest

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.context import DashboardContext


@pytest.fixture
def tab(tmp_path):
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    (tmp_path / "calc" / "note.txt").write_text("line\n" * 40, encoding="utf-8")

    sent: list[str] = []
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    ctx.run_js = sent.append
    _widget, refs = browser.create_tab(ctx)
    refs["calc_list_directory"]()
    return refs, sent, ctx


def _label_of(file_list, name):
    for option in file_list.options:
        text = option[1] if isinstance(option, tuple) else option
        if name in str(text):
            return text
    raise AssertionError(f"{name} is not in the list: {file_list.options}")


def test_a_click_on_the_selected_file_still_loads_it(tab):
    refs, sent, _ctx = tab
    file_list = refs["calc_file_list"]
    label = _label_of(file_list, "note.txt")

    file_list.value = (label,)               # first click: selection changes
    assert sent, "selecting the file did not open it"

    # Second click on the same entry: the selection is unchanged, so only the
    # click bridge can carry it.
    sent.clear()
    refs["calc_open_input"].value = f"{label}\x1f2"

    assert sent, "clicking the already selected file loaded nothing"


def test_the_bridge_re_arms_itself(tab):
    refs, _sent, ctx = tab
    label = _label_of(refs["calc_file_list"], "note.txt")

    refs["calc_open_input"].value = f"{label}\x1f1"

    # Left standing, the next click with the same counter would be swallowed.
    assert refs["calc_open_input"].value == ""


def test_a_folder_click_does_not_open_anything(tab, tmp_path):
    refs, sent, _ctx = tab
    sent.clear()

    refs["calc_open_input"].value = "(no such entry)\x1f1"

    assert not sent, "a non-file entry must not paint the content area"


def test_the_click_is_reported_by_the_browser_with_a_counter(tab):
    """The page must send a value that differs for every click."""
    _refs, sent, ctx = tab
    page = "\n".join(list(ctx.init_js_parts) + sent)

    assert "calc-cmd-open" in page, "the click is never reported to the kernel"
    assert "_openSeq" in page, "without a counter a repeated click is dropped"
