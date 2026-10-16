"""A write into shared lessons or a memory store is asked in every mode.

Input: write_file calls through the real executor (``_doc_executor.
execute``) in default, acceptEdits and bypassPermissions,
with a confirm callback that records the preview and answers no, and
with no callback at all. Targets: ``delfin/agent/pack/knowledge/`` (the
global lessons), ``.delfin/knowledge/`` of a repository (its project
lessons), and, with the workspace at a temporary HOME, the stores under
``~/.delfin``: the per-project and user-wide memory, the private
lessons, the swarm findings and the knowledge review queue.

Output: whether the callback was asked, the preview it was shown, the
tool's answer, and the file on disk.

Semantics pinned here: each target is asked in every mode (bypass
included) with a preview that names the store, is refused where nobody
can be asked, and is left unwritten after a no. A plain file in the
same workspace is written unasked in acceptEdits (the control that the
gate is not simply asking everything). The cage's read-only binds
(``_caged_readonly_paths``) include the existing stores.

All paths are under ``tmp_path``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _doc_executor

MODES = ("default", "acceptEdits", "bypassPermissions")

REPO_TARGETS = (
    "delfin/agent/pack/knowledge/lesson_x.md",
    ".delfin/knowledge/lesson_y.md",
    "sub/project/.delfin/knowledge/lesson_z.md",
)
HOME_TARGETS = (
    ".delfin/projects/-ws/memory/project_note.md",
    ".delfin/memory/user_note.md",
    ".delfin/knowledge/private/lessons.json",
    ".delfin/knowledge/queue/entry.json",
    ".delfin/swarm/run1/entries.jsonl",
)


@pytest.fixture
def home(tmp_path, monkeypatch):
    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setattr(Path, "home", lambda: h)
    monkeypatch.delenv("DELFIN_KNOWLEDGE_QUEUE", raising=False)
    import delfin.user_settings as us
    monkeypatch.setattr(us, "load_settings", lambda *a, **k: {})
    return h


def _write(ws: Path, rel: str, mode: str, callback, *, root=None,
           extra=()):
    """write_file of a NEW file ``(root or ws)/rel``; returns the answer
    and the file's text afterwards ("" when it was not written)."""
    perms = KitToolPermissions(workspace=ws, mode=mode,
                               confirm_callback=callback,
                               extra_workspace_dirs=tuple(extra))
    target = (root or ws) / rel
    target.parent.mkdir(parents=True, exist_ok=True)
    arg = rel if root is None else str(target)
    out = _doc_executor.execute(
        "write_file", {"path": arg, "content": "changed\n"},
        permissions=perms)
    return out, (target.read_text(encoding="utf-8")
                 if target.exists() else "")


def _ws(tmp_path) -> Path:
    ws = tmp_path / "ws"
    ws.mkdir(exist_ok=True)
    return ws


def _asked():
    seen: list[str] = []

    def cb(_name, _args, preview):
        seen.append(preview)
        return False
    return cb, seen


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("rel", REPO_TARGETS)
def test_a_lesson_write_is_asked_in_every_mode(tmp_path, home, mode, rel):
    ws = _ws(tmp_path)
    cb, seen = _asked()
    out, after = _write(ws, rel, mode, cb)
    assert seen, f"{mode}: no question for {rel}: {out}"
    assert "KNOWLEDGE" in seen[0] or "MEMORY" in seen[0], seen[0][:200]
    assert after == ""


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("rel", HOME_TARGETS)
def test_a_memory_store_write_is_asked_in_every_mode(tmp_path, home, mode,
                                                     rel):
    """The home directory is never a workspace; the stores are reachable
    through a granted directory, here ``--add-dir ~/.delfin``."""
    cb, seen = _asked()
    out, after = _write(_ws(tmp_path), rel, mode, cb, root=home,
                        extra=(home / ".delfin",))
    assert seen, f"{mode}: no question for ~/{rel}: {out}"
    assert "KNOWLEDGE" in seen[0] or "MEMORY" in seen[0], seen[0][:200]
    assert after == ""


def test_with_nobody_to_ask_the_write_is_refused(tmp_path, home):
    ws = _ws(tmp_path)
    out, after = _write(ws, REPO_TARGETS[0], "bypassPermissions", None)
    assert "refused" in out.lower(), out
    assert after == ""
    out, after = _write(ws, HOME_TARGETS[0], "bypassPermissions", None,
                        root=home, extra=(home / ".delfin",))
    assert "refused" in out.lower(), out
    assert after == ""


def test_a_plain_file_is_written_unasked_in_accept_edits(tmp_path, home):
    ws = _ws(tmp_path)
    cb, seen = _asked()
    _out, after = _write(ws, "notes/plain.md", "acceptEdits", cb)
    assert not seen
    assert after == "changed\n"


def test_the_cage_binds_the_existing_stores_read_only(tmp_path, home):
    ws = tmp_path / "ws"
    for rel in REPO_TARGETS[:2]:
        (ws / rel).parent.mkdir(parents=True, exist_ok=True)
    for rel in (".delfin/memory", ".delfin/knowledge",
                ".delfin/projects/-ws/memory"):
        (home / rel).mkdir(parents=True, exist_ok=True)
    perms = KitToolPermissions(workspace=ws, mode="bypassPermissions")
    ro = set(A._caged_readonly_paths(perms))
    import os
    for p in (ws / "delfin/agent/pack/knowledge", ws / ".delfin/knowledge",
              home / ".delfin/memory", home / ".delfin/knowledge",
              home / ".delfin/projects/-ws/memory"):
        assert os.path.realpath(p) in ro, p
