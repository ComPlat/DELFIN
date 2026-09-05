"""Undo has to return the file's BYTES, or refuse — never something else.

Two ways it did neither, both reproduced against the real modules:

* The journal stored the write path's NORMALISED text (LF, no BOM,
  decoded) and hashed that, while ``revert`` verifies against the file's
  RAW bytes. A CRLF csv and a cp1252 .tex both came back "content
  changed since the agent's edit (post-hash mismatch) — not overwriting"
  when nobody had touched them: undo was structurally unavailable for
  every Windows-authored file in the workspace, and the message blamed a
  third party for it.

* ``apply_patch`` read both the pre-image and the post-image with
  ``decode("utf-8", errors="replace")``. Because BOTH sides were damaged
  the same way, the post-hash guard matched and the restore fired:
  ``b"caf\\xe9 alpha\\n"`` came back as ``b"caf\\xef\\xbf\\xbd alpha\\n"``
  and was reported as reverted.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from delfin.agent import change_journal as cj
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


@pytest.fixture
def home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "home")
    (tmp_path / "home").mkdir()
    return tmp_path / "home"


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "ws"
    d.mkdir()
    return d


def _perms(ws: Path, sid: str) -> KitToolPermissions:
    p = KitToolPermissions(workspace=ws, mode="default")
    p.task_session_id = sid
    return p


def _write_through_the_tool(ws: Path, sid: str, name: str, content: str) -> str:
    """Drive the real write executor, read baseline included."""
    ex = _DocToolExecutor()
    perms = _perms(ws, sid)
    target = ws / name
    if target.exists():
        perms.read_tracker[str(target)] = target.stat().st_mtime
    return ex._execute_write_file({"path": name, "content": content}, perms)


# ---------------------------------------------------------------------------
# The file's own convention survives the round trip
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name,before,after_text", [
    ("windows.csv", b"a,b\r\nc,d\r\n", "a,b\nc,D\n"),
    ("bom.txt", "﻿first\nsecond\n".encode("utf-8"), "first\nSECOND\n"),
    ("latin.tex", b"caf\xe9 au lait\n", "café au lait, please\n"),
    ("plain.txt", b"one\ntwo\n", "one\nTWO\n"),
])
def test_a_write_can_be_undone_whatever_the_file_is(
        home, ws, name, before, after_text):
    (ws / name).write_bytes(before)
    result = _write_through_the_tool(ws, "s-conv", name, after_text)
    assert "error" not in result[:20], result
    assert (ws / name).read_bytes() != before        # the write happened

    out = cj.revert("s-conv", scope="last", workspace=ws)
    assert out["conflicts"] == [], out["conflicts"]
    assert out["reverted"] == [str(ws / name)]
    assert (ws / name).read_bytes() == before, (
        "the undo did not put the original bytes back")


def test_the_conflict_guard_still_refuses_a_file_the_user_changed(home, ws):
    (ws / "win.csv").write_bytes(b"a,b\r\nc,d\r\n")
    _write_through_the_tool(ws, "s-guard", "win.csv", "a,b\nc,D\n")
    (ws / "win.csv").write_bytes(b"the user typed this\r\n")

    out = cj.revert("s-guard", scope="last", workspace=ws)
    assert out["reverted"] == []
    assert "content changed since the agent's edit" in out["conflicts"][0]["reason"]
    assert (ws / "win.csv").read_bytes() == b"the user typed this\r\n"


# ---------------------------------------------------------------------------
# apply_patch: bytes it cannot decode go down the byte-exact path
# ---------------------------------------------------------------------------

def _git_ws(ws: Path) -> bool:
    if shutil.which("git") is None:
        return False
    subprocess.run(["git", "init", "-q", str(ws)], check=True)
    return True


def test_a_patch_to_a_non_utf8_file_undoes_byte_for_byte(home, ws):
    if not _git_ws(ws):
        pytest.skip("git not available")
    target = ws / "cafe.txt"
    original = b"caf\xe9 alpha\nbeta\n"
    target.write_bytes(original)

    ex = _DocToolExecutor()
    diff = "--- a/cafe.txt\n+++ b/cafe.txt\n@@ -2 +2 @@\n-beta\n+gamma\n"
    result = ex._execute_apply_patch({"diff": diff}, _perms(ws, "s-patch"))
    assert '"status": "ok"' in result, result
    assert target.read_bytes() == b"caf\xe9 alpha\ngamma\n"

    out = cj.revert("s-patch", scope="last", workspace=ws)
    assert out["reverted"] == [str(target)], out
    assert target.read_bytes() == original, (
        "the undo wrote the replacement character over the user's bytes")


def test_an_ordinary_patch_still_applies_and_undoes(home, ws):
    if not _git_ws(ws):
        pytest.skip("git not available")
    target = ws / "note.txt"
    target.write_text("one\ntwo\n", encoding="utf-8")

    ex = _DocToolExecutor()
    diff = "--- a/note.txt\n+++ b/note.txt\n@@ -2 +2 @@\n-two\n+TWO\n"
    result = ex._execute_apply_patch({"diff": diff}, _perms(ws, "s-plain"))
    assert '"status": "ok"' in result, result
    assert target.read_text(encoding="utf-8") == "one\nTWO\n"

    out = cj.revert("s-plain", scope="last", workspace=ws)
    assert out["reverted"] == [str(target)]
    assert target.read_text(encoding="utf-8") == "one\ntwo\n"


def test_a_pre_image_that_cannot_be_stored_exactly_is_refused(home, ws):
    """A caller that read the file with errors="replace" hands over text
    with U+FFFD in it. Writing that back would change bytes the agent
    never touched, and the post-hash guard cannot catch it because the
    damage is in the pre-image."""
    p = ws / "note.tex"
    p.write_bytes(b"caf\xe9 x\n")
    lossy_pre = b"caf\xe9 alpha\n".decode("utf-8", errors="replace")
    cj.record_change("s-lossy", tool="notebook_edit", path=p,
                     old_text=lossy_pre, new_text="café x\n")

    rec = cj.list_changes("s-lossy")[-1]
    assert rec["lossy"] is True
    out = cj.revert("s-lossy", scope="last", workspace=ws)
    assert out["reverted"] == []
    assert "could not be stored" in out["skipped"][0]["reason"]
    assert p.read_bytes() == b"caf\xe9 x\n"


def test_an_ordinary_pre_image_is_not_flagged_lossy(home, ws):
    p = ws / "plain.txt"
    p.write_text("new\n", encoding="utf-8")
    cj.record_change("s-fine", tool="write_file", path=p,
                     old_text="old\n", new_text="new\n")
    assert "lossy" not in cj.list_changes("s-fine")[-1]


def test_the_python_applier_refuses_a_file_it_cannot_decode(ws):
    """No git in the workspace: the fallback used to raise
    UnicodeDecodeError out of the tool call, and rewriting the file from
    a lossy decode would have been worse."""
    from delfin.agent import patch_apply as pa

    (ws / "cafe.txt").write_bytes(b"caf\xe9 alpha\nbeta\n")
    out = pa.apply_patch(
        ws, "--- a/cafe.txt\n+++ b/cafe.txt\n@@ -2 +2 @@\n-beta\n+gamma\n",
        prefer="py")
    assert out["status"] == "check_failed"
    assert "UTF-8" in out["error"]
    assert (ws / "cafe.txt").read_bytes() == b"caf\xe9 alpha\nbeta\n"
