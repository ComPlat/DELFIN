"""A patch that removes a file must not take the only copy with it.

``apply_patch`` read every pre-image into memory BEFORE applying, then
journalled by re-reading each file afterwards — with ``except OSError:
continue`` around it. For a deletion that read is exactly the one that
fails, so the branch that skipped the journal was the branch holding the
last copy of the file. The user was told ``status: ok``.

A rename hit it from the other side: the diff scanner returned only the
old path (the first group of ``diff --git a/old b/new``), so the removed
file went down the same skip and the new path was never journalled at
all.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from delfin.agent import change_journal as cj
from delfin.agent import patch_apply as pa
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


@pytest.fixture
def home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "home")
    (tmp_path / "home").mkdir()
    return tmp_path / "home"


@pytest.fixture
def gitws(tmp_path):
    if shutil.which("git") is None:
        pytest.skip("git not available")
    d = tmp_path / "ws"
    d.mkdir()
    subprocess.run(["git", "init", "-q", str(d)], check=True)
    subprocess.run(["git", "-C", str(d), "config", "user.email", "a@b.c"],
                   check=True)
    subprocess.run(["git", "-C", str(d), "config", "user.name", "t"], check=True)
    return d


def _perms(ws: Path, sid: str) -> KitToolPermissions:
    p = KitToolPermissions(workspace=ws, mode="default")
    p.task_session_id = sid
    return p


_DELETE_DIFF = (
    "diff --git a/doomed.txt b/doomed.txt\n"
    "deleted file mode 100644\n"
    "--- a/doomed.txt\n"
    "+++ /dev/null\n"
    "@@ -1,2 +0,0 @@\n"
    "-one\n"
    "-two\n"
)


def test_a_deleted_file_is_journalled_and_comes_back(home, gitws):
    target = gitws / "doomed.txt"
    target.write_text("one\ntwo\n", encoding="utf-8")

    ex = _DocToolExecutor()
    out = ex._execute_apply_patch({"diff": _DELETE_DIFF}, _perms(gitws, "s-del"))
    assert '"status": "ok"' in out, out
    assert not target.exists()

    records = cj.list_changes("s-del")
    assert records, "the only copy of the file was dropped on the floor"
    assert records[-1]["deleted"] is True
    assert records[-1]["path"] == str(target)

    res = cj.revert("s-del", scope="last", workspace=gitws)
    assert res["reverted"] == [str(target)], res
    assert target.read_text(encoding="utf-8") == "one\ntwo\n"


def test_a_recreated_file_is_not_overwritten_by_the_undo(home, gitws):
    target = gitws / "doomed.txt"
    target.write_text("one\ntwo\n", encoding="utf-8")
    ex = _DocToolExecutor()
    ex._execute_apply_patch({"diff": _DELETE_DIFF}, _perms(gitws, "s-del2"))

    target.write_text("the user made a new one\n", encoding="utf-8")
    res = cj.revert("s-del2", scope="last", workspace=gitws)
    assert res["reverted"] == []
    assert "exists again" in res["conflicts"][0]["reason"]
    assert target.read_text(encoding="utf-8") == "the user made a new one\n"


def test_undoing_the_undo_removes_the_file_again(home, gitws):
    target = gitws / "doomed.txt"
    target.write_text("one\ntwo\n", encoding="utf-8")
    ex = _DocToolExecutor()
    ex._execute_apply_patch({"diff": _DELETE_DIFF}, _perms(gitws, "s-del3"))

    res = cj.revert("s-del3", scope="last", workspace=gitws)
    assert target.exists()
    undo_seq = res["undo_seqs"][0]

    back = cj.revert("s-del3", scope="turn", turn_seqs=[undo_seq],
                     workspace=gitws)
    assert back["reverted"] == [str(target)], back
    assert not target.exists()


def test_a_rename_names_both_sides(gitws):
    """The scanner is what decides which paths get a pre-image."""
    diff = (
        "diff --git a/old.txt b/new.txt\n"
        "similarity index 100%\n"
        "rename from old.txt\n"
        "rename to new.txt\n"
    )
    assert pa._files_in_diff(diff) == ["old.txt", "new.txt"]


def test_a_renamed_file_can_be_put_back(home, gitws):
    old = gitws / "old.txt"
    old.write_text("keep me\n", encoding="utf-8")
    subprocess.run(["git", "-C", str(gitws), "add", "old.txt"], check=True)
    subprocess.run(["git", "-C", str(gitws), "commit", "-qm", "x"], check=True)

    diff = (
        "diff --git a/old.txt b/new.txt\n"
        "similarity index 100%\n"
        "rename from old.txt\n"
        "rename to new.txt\n"
    )
    ex = _DocToolExecutor()
    out = ex._execute_apply_patch({"diff": diff}, _perms(gitws, "s-ren"))
    assert '"status": "ok"' in out, out
    assert not old.exists() and (gitws / "new.txt").exists()

    paths = {r["path"] for r in cj.list_changes("s-ren")}
    assert str(old) in paths, "the removed side was never journalled"
    assert str(gitws / "new.txt") in paths, "the new side was never journalled"

    res = cj.revert("s-ren", scope="session", workspace=gitws)
    assert res["conflicts"] == [], res
    assert old.read_text(encoding="utf-8") == "keep me\n"
    assert not (gitws / "new.txt").exists()


def test_a_created_file_still_reverts_by_deletion(home, gitws):
    diff = (
        "diff --git a/fresh.txt b/fresh.txt\n"
        "new file mode 100644\n"
        "--- /dev/null\n"
        "+++ b/fresh.txt\n"
        "@@ -0,0 +1 @@\n"
        "+hello\n"
    )
    ex = _DocToolExecutor()
    out = ex._execute_apply_patch({"diff": diff}, _perms(gitws, "s-new"))
    assert '"status": "ok"' in out, out
    assert (gitws / "fresh.txt").exists()

    res = cj.revert("s-new", scope="last", workspace=gitws)
    assert res["reverted"] == [str(gitws / "fresh.txt")], res
    assert not (gitws / "fresh.txt").exists()
