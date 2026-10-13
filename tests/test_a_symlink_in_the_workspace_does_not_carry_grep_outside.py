"""A symlink in the workspace does not carry grep_file outside it.

The read gate asks before a file outside the allowed folders is read.
``grep_file`` gates the folder it searches and then opened every file the
walk found there, symlinks included: a link in the workspace pointing at a
file elsewhere was read line by line with no dialog, while reading the same
link by name was asked about. A link is a second path whose target nobody
was asked about, so the search does not follow it out of the folder.

Driven through the document executor's public ``execute``, with a confirm
callback that records every question and answers no.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor

SECRET = "MARKER-OUTSIDE-7f3a"


@pytest.fixture
def box(tmp_path, monkeypatch):
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "home")
    (tmp_path / "home").mkdir()
    ws = tmp_path / "ws"
    (ws / "sub").mkdir(parents=True)
    out = tmp_path / "outside"
    out.mkdir()
    (out / "notes.txt").write_text(f"{SECRET} line outside\n", encoding="utf-8")
    (ws / "a.txt").write_text(f"{SECRET} line inside\n", encoding="utf-8")
    asked: list = []
    perms = KitToolPermissions(workspace=ws, mode="default",
                               task_session_id="grep-symlink")
    perms.confirm_callback = lambda *a, **k: (asked.append(a) or False)
    return ws, out, perms, asked


def _grep(perms, path="."):
    return str(_doc_executor.execute(
        "grep_file", {"pattern": "MARKER", "path": path}, perms))


def test_a_file_link_to_outside_is_not_read_by_a_folder_search(box):
    ws, out, perms, asked = box
    (ws / "link.txt").symlink_to(out / "notes.txt")
    res = _grep(perms)
    assert "line outside" not in res
    assert "a.txt:1:" in res                      # the search itself still works
    assert "1 symlink(s) pointing outside" in res  # and says what it left out


def test_the_same_link_in_a_subfolder_is_not_read_either(box):
    ws, out, perms, _ = box
    (ws / "sub" / "deep.txt").symlink_to(out / "notes.txt")
    assert "line outside" not in _grep(perms)
    assert "line outside" not in _grep(perms, "sub")


def test_a_relative_link_that_climbs_out_is_not_read(box):
    ws, out, perms, _ = box
    (ws / "sub" / "up.txt").symlink_to(Path("..") / ".." / "outside" / "notes.txt")
    assert "line outside" not in _grep(perms)


def test_a_chain_of_links_ending_outside_is_not_read(box):
    ws, out, perms, _ = box
    (ws / "second.txt").symlink_to(out / "notes.txt")
    (ws / "first.txt").symlink_to(ws / "second.txt")
    assert "line outside" not in _grep(perms)


def test_a_link_that_stays_inside_the_folder_is_still_searched(box):
    ws, _, perms, _ = box
    (ws / "alias.txt").symlink_to(ws / "a.txt")
    res = _grep(perms)
    assert "alias.txt:1:" in res
    assert "pointing outside" not in res


def test_an_innocent_name_linked_to_a_secret_file_inside_reads_nothing(box):
    ws, _, perms, _ = box
    (ws / ".env").write_text(f"{SECRET}=in-the-env-file\n", encoding="utf-8")
    (ws / "readme_copy.txt").symlink_to(ws / ".env")
    assert "in-the-env-file" not in _grep(perms)


def test_naming_the_link_directly_asks_and_a_no_reads_nothing(box):
    ws, out, perms, asked = box
    (ws / "link.txt").symlink_to(out / "notes.txt")
    res = _grep(perms, "link.txt")
    assert "line outside" not in res
    assert len(asked) == 1, asked


def test_a_dangling_link_does_not_break_the_search(box):
    ws, out, perms, _ = box
    (ws / "gone.txt").symlink_to(out / "never_there.txt")
    assert "a.txt:1:" in _grep(perms)


def test_a_linked_folder_is_not_walked(box):
    """os.walk does not descend into a symlinked directory; this pins it,
    because the file rule above would not catch files under one."""
    ws, out, perms, _ = box
    (ws / "outdir").symlink_to(out, target_is_directory=True)
    assert "line outside" not in _grep(perms)
