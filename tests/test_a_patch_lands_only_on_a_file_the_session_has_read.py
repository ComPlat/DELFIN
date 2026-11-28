"""A patch is held to the read-before-write rule the other write tools hold.

``write_file`` refuses an existing file it has not read ("refusing to
overwrite ... without a prior read_file"), ``edit_file`` and
``multi_edit`` refuse without a read baseline, ``notebook_edit`` refuses
without ``notebook_read``, and a shell redirect is refused by the same
rule (``_bash_target_needs_a_read``, added after an agent destroyed a
concurrent user edit with ``sed -i``).

``apply_patch`` ran none of it. Nothing in ``_execute_apply_patch``,
``patch_apply.py`` or ``patch_targets.py`` reads ``perms.read_tracker``,
so the one tool a model uses to change code by context could rewrite a
file the session never opened, and land on top of a change another
process made after the session read it. The same content, reached
through the tool the agent happened to pick.

The rule is the same one, so the exemption is too: a file the patch
CREATES has nothing to clobber and needs no baseline, and ``check_only``
writes nothing -- it stays a dry run, usable in plan mode.

Driven through the real dispatch (``_doc_executor.execute``, the call
delfin/agent/engine.py:978 makes for every tool call of a turn), so the
gate, the mode pin and the executor all run. Measured red on main
(7fafac914) by name: the three refusal cases below applied their patch
instead of refusing it.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor

ORIGINAL = "one\ntwo\nthree\n"
PATCHED = "one\nTWO\nthree\n"

# Changes line 2; the context lines are 1 and 3.
PATCH_MODIFY = (
    "--- a/f.txt\n+++ b/f.txt\n@@ -1,3 +1,3 @@\n one\n-two\n+TWO\n three\n"
)
PATCH_NEW = "--- /dev/null\n+++ b/new.txt\n@@ -0,0 +1,1 @@\n+hello\n"
# Two files in one diff: a.txt is read first in the mixed case, b.txt is not.
PATCH_TWO_FILES = (
    "--- a/a.txt\n+++ b/a.txt\n@@ -1 +1 @@\n-a\n+A\n"
    "--- a/b.txt\n+++ b/b.txt\n@@ -1 +1 @@\n-b\n+B\n"
)


@pytest.fixture
def ws(tmp_path) -> Path:
    d = tmp_path / "ws"
    d.mkdir()
    return d


def _perms(ws: Path, sid: str = "s-patch", mode: str = "acceptEdits"):
    p = KitToolPermissions(workspace=str(ws), mode=mode)
    p.task_session_id = sid
    return p


def _patch(perms, diff: str, **extra):
    """One apply_patch call through the product's tool dispatch."""
    args = {"diff": diff, **extra}
    return json.loads(_doc_executor.execute("apply_patch", args, perms))


def _read(perms, rel: str):
    return _doc_executor.execute("read_file", {"path": rel}, perms)


def _stale_later(ws: Path, rel: str, text: str) -> None:
    """Another process writes the file, deterministically later.

    The staleness rule allows 1e-3 s; a real second writer can land
    inside one clock tick, so the mtime is pinned rather than trusted to
    wall-clock granularity.
    """
    f = ws / rel
    f.write_text(text, encoding="utf-8")
    st = f.stat()
    import os
    os.utime(f, ns=(int((st.st_atime + 10) * 1e9),
                    int((st.st_mtime + 10) * 1e9)))


# ---------------------------------------------------------------------------
# controls: what must keep working, unchanged, before and after
# ---------------------------------------------------------------------------

def test_a_patch_to_a_file_the_session_read_lands(ws):
    (ws / "f.txt").write_text(ORIGINAL)
    perms = _perms(ws)
    _read(perms, "f.txt")
    out = _patch(perms, PATCH_MODIFY)
    assert out.get("status") == "ok", out
    assert (ws / "f.txt").read_text() == PATCHED


def test_a_patch_that_creates_a_file_needs_no_read(ws):
    """The exemption write_file follows: nothing to clobber, nothing to
    have read. A fix that refuses this broke the tool's main use.

    The newline the applier leaves off the created file is its own
    behaviour and not this rule's business, so the content is compared
    stripped; a refusal, or a file that never appears, still fails here.
    """
    perms = _perms(ws)
    out = _patch(perms, PATCH_NEW)
    assert out.get("status") == "ok", out
    assert (ws / "new.txt").read_text().strip() == "hello"


def test_a_dry_run_needs_no_read_and_writes_nothing(ws):
    """check_only changes nothing, so it needs no baseline -- and must stay
    allowed where writing is not (plan mode)."""
    (ws / "f.txt").write_text(ORIGINAL)
    perms = _perms(ws)
    out = _patch(perms, PATCH_MODIFY, check_only=True)
    assert "error" not in out, out
    assert (ws / "f.txt").read_text() == ORIGINAL


# ---------------------------------------------------------------------------
# the rule
# ---------------------------------------------------------------------------

def test_a_blind_patch_is_refused_and_writes_nothing(ws):
    """RED ON MAIN: the patch applied and the original bytes were gone.

    Fails when apply_patch lets a diff through without consulting the
    session's read baseline -- which is the case where the agent has
    never seen the file it is overwriting.
    """
    (ws / "f.txt").write_text(ORIGINAL)
    perms = _perms(ws)
    out = _patch(perms, PATCH_MODIFY)
    assert "error" in out, out
    assert "without a prior read_file" in out["error"]
    assert (ws / "f.txt").read_text() == ORIGINAL


def test_a_patch_over_a_change_made_after_the_read_is_refused(ws):
    """RED ON MAIN: the patch's own context lines still matched (the other
    process edited a line outside the hunk), so it applied and destroyed
    that change.

    The refusal must be the staleness one, not "unread" and not a hunk
    failure -- the agent has to learn the file moved under it.
    """
    moved = ORIGINAL + "four\nfive\nSIX-BY-OTHERS\n"
    (ws / "f.txt").write_text(ORIGINAL + "four\nfive\nsix\n")
    perms = _perms(ws)
    _read(perms, "f.txt")
    _stale_later(ws, "f.txt", moved)

    out = _patch(perms, PATCH_MODIFY)
    assert "error" in out, out
    assert "modified since last read_file" in out["error"]
    assert (ws / "f.txt").read_text() == moved


def test_a_multi_file_patch_refuses_when_one_target_is_unread(ws):
    """RED ON MAIN: both files were patched.

    One read file does not license a second the session never opened, and
    the answer must name the unread one; nothing in the diff is written.
    """
    (ws / "a.txt").write_text("a\n")
    (ws / "b.txt").write_text("b\n")
    perms = _perms(ws)
    _read(perms, "a.txt")

    out = _patch(perms, PATCH_TWO_FILES)
    assert "error" in out, out
    assert "b.txt" in out["error"]
    assert (ws / "a.txt").read_text() == "a\n"
    assert (ws / "b.txt").read_text() == "b\n"


def test_the_staging_path_refuses_a_blind_patch(tmp_path, monkeypatch):
    """diff_approval does not write -- it stages a change for a person to
    approve. Staging a change to a file nobody read asks them to approve
    something they cannot check, so the rule holds before staging too.

    The pending store is Path.home()/.delfin/pending, so HOME is moved --
    to a SIBLING of the workspace, because a workspace that is an ancestor
    of $HOME is refused by _is_forbidden_workspace_root (9548) and the
    first version of this test tripped exactly that.
    """
    fake_home = tmp_path / "fakehome"
    (fake_home / ".delfin" / "pending").mkdir(parents=True)
    monkeypatch.setattr(Path, "home", lambda: fake_home)
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "f.txt").write_text(ORIGINAL)
    perms = _perms(ws, mode="default")
    perms.mode = "diff_approval"

    out = _patch(perms, PATCH_MODIFY)
    assert "error" in out, out
    assert "without a prior read_file" in out["error"]
    assert (ws / "f.txt").read_text() == ORIGINAL
