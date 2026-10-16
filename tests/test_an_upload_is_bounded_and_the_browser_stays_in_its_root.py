"""Uploads into the Calculations browser are bounded; navigation stays in the root.

Driven through the real tab (``tab_calculations_browser.create_tab``) on a
tmp calc directory: the hidden upload widget receives a value the way the
browser front end sets it, and a double click is the command text the page
writes into the hidden ``calc-cmd-dblclick`` input.

* A ZIP whose unpacked size exceeds ``DELFIN_ZIP_MAX_UNPACKED_BYTES``, or
  whose entry count exceeds ``DELFIN_ZIP_MAX_ENTRIES``, leaves nothing
  unpacked; a file over ``DELFIN_UPLOAD_MAX_BYTES`` is not written. The
  sample archive is 64 KiB of zeros, compressed to a few hundred bytes.
* A double click on a link that points out of the root does not enter it.
"""

from __future__ import annotations

import io
import pathlib
import zipfile

import pytest

from delfin.dashboard.context import DashboardContext


def _walk(widget):
    yield widget
    for child in getattr(widget, "children", ()) or ():
        yield from _walk(child)


@pytest.fixture
def browser(tmp_path, monkeypatch):
    widgets = pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_calculations_browser

    for name in ("calc", "archive", "office", "home"):
        (tmp_path / name).mkdir()
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(tmp_path))
    monkeypatch.setattr(pathlib.Path, "home", lambda: tmp_path / "home")
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    ctx.run_js = lambda *a, **k: None
    tab, refs = tab_calculations_browser.create_tab(ctx)
    found = list(_walk(tab))
    upload = next(w for w in found if isinstance(w, widgets.FileUpload))
    dblclick = next(w for w in found if isinstance(w, widgets.Text)
                    and "calc-cmd-dblclick" in getattr(w, "_dom_classes", ()))
    return {"tmp": tmp_path, "calc": tmp_path / "calc", "refs": refs,
            "upload": upload, "dblclick": dblclick}


def _zip(entries) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, data in entries:
            archive.writestr(name, data)
    return buffer.getvalue()


def _send(upload, name: str, content: bytes) -> None:
    import datetime

    upload.value = ({"name": name, "type": "", "size": len(content),
                     "content": memoryview(content),
                     "last_modified": datetime.datetime.now(datetime.timezone.utc)},)


def _files_below(path: pathlib.Path):
    return [p for p in path.rglob("*") if p.is_file()]


def test_a_zip_that_unpacks_past_the_limit_leaves_nothing(browser, monkeypatch):
    monkeypatch.setenv("DELFIN_ZIP_MAX_UNPACKED_BYTES", "4096")
    _send(browser["upload"], "big.zip", _zip([("zeros.bin", b"\0" * 65536)]))
    assert not any(p.name == "zeros.bin" for p in _files_below(browser["calc"]))


def test_a_zip_with_too_many_entries_leaves_nothing(browser, monkeypatch):
    monkeypatch.setenv("DELFIN_ZIP_MAX_ENTRIES", "5")
    _send(browser["upload"], "many.zip",
          _zip([(f"f{i}.txt", b"x") for i in range(10)]))
    assert _files_below(browser["calc"]) == []


def test_a_file_over_the_upload_limit_is_not_written(browser, monkeypatch):
    monkeypatch.setenv("DELFIN_UPLOAD_MAX_BYTES", "1000")
    _send(browser["upload"], "big.txt", b"y" * 2000)
    assert _files_below(browser["calc"]) == []


def test_an_upload_within_the_limits_still_lands(browser, monkeypatch):
    monkeypatch.setenv("DELFIN_UPLOAD_MAX_BYTES", "100000")
    monkeypatch.setenv("DELFIN_ZIP_MAX_UNPACKED_BYTES", "100000")
    _send(browser["upload"], "small.zip",
          _zip([("job/a.xyz", b"1\n\nH 0 0 0\n")]))
    _send(browser["upload"], "note.txt", b"hello")
    names = sorted(p.name for p in _files_below(browser["calc"]))
    assert names == ["a.xyz", "note.txt"]


def test_a_double_click_does_not_follow_a_link_out_of_the_root(browser):
    outside = browser["tmp"] / "outside"
    outside.mkdir()
    (outside / "secret.txt").write_text("x", encoding="utf-8")
    (browser["calc"] / "inside").mkdir()
    (browser["calc"] / "away").symlink_to(outside, target_is_directory=True)
    refs = browser["refs"]
    refs["calc_list_directory"]()
    labels = list(refs["calc_file_list"].options)
    away = next(label for label in labels if label.endswith("away"))
    inside = next(label for label in labels if label.endswith("inside"))

    browser["dblclick"].value = away
    assert refs["sheet_state"]["current_path"] == ""

    browser["dblclick"].value = inside
    assert refs["sheet_state"]["current_path"] == "inside"


# -- the limits themselves -------------------------------------------------------

def test_the_limits_read_their_variables_and_fall_back_to_defaults(monkeypatch):
    from delfin.dashboard import upload_limits as UL

    monkeypatch.setenv(UL.UPLOAD_ENV, "123")
    monkeypatch.setenv(UL.UNPACKED_ENV, "not a number")
    monkeypatch.setenv(UL.ENTRIES_ENV, "-4")
    assert UL.max_upload_bytes() == 123
    assert UL.max_unpacked_bytes() == UL.DEFAULT_UNPACKED_BYTES
    assert UL.max_entries() == UL.DEFAULT_ENTRIES
    with pytest.raises(UL.UploadRefused, match=UL.UPLOAD_ENV):
        UL.check_upload_size(124, "x")
    UL.check_upload_size(123, "x")


def test_a_chunk_announcement_cannot_reserve_unbounded_memory():
    from delfin.dashboard import upload_limits as UL

    with pytest.raises(UL.UploadRefused):
        UL.check_chunk_plan(UL.MAX_CHUNKS + 1, 0)
    UL.check_chunk_plan(10, 100)


def test_a_new_file_is_never_written_through_an_existing_link(tmp_path):
    from delfin.dashboard import upload_limits as UL

    target = tmp_path / "elsewhere.txt"
    target.write_text("kept", encoding="utf-8")
    link = tmp_path / "planted"
    link.symlink_to(target)
    with pytest.raises(FileExistsError):
        UL.write_new_bytes(link, b"new")
    assert target.read_text(encoding="utf-8") == "kept"


def test_a_budget_overrun_leaves_no_partial_file(tmp_path):
    from delfin.dashboard import upload_limits as UL

    destination = tmp_path / "out.bin"
    with pytest.raises(UL.UploadRefused):
        UL.write_new_file(destination, io.BytesIO(b"\0" * 5000),
                          UL.Budget(1000))
    assert not destination.exists()


def test_a_file_linked_from_outside_the_root_is_not_shown(browser):
    outside = browser["tmp"] / "outside"
    outside.mkdir()
    (outside / "secret.txt").write_text("outside-marker-7d1", encoding="utf-8")
    (browser["calc"] / "leak.txt").symlink_to(outside / "secret.txt")
    (browser["calc"] / "plain.txt").write_text("inside-marker-3c9",
                                               encoding="utf-8")
    refs = browser["refs"]
    refs["calc_list_directory"]()
    labels = list(refs["calc_file_list"].options)

    state = refs["sheet_state"]

    def shown():
        return (refs["calc_content_area"].value + str(state.get("file_content") or "")
                + str(state.get("selected_file_path") or ""))

    refs["calc_open_input"].value = next(l for l in labels if l.endswith("leak.txt"))
    assert "outside-marker-7d1" not in shown()
    assert "secret.txt" not in shown() and "leak.txt" not in str(
        state.get("selected_file_path") or "")

    refs["calc_open_input"].value = next(l for l in labels if l.endswith("plain.txt"))
    assert "plain.txt" in str(state.get("selected_file_path") or "")
