"""Undo journal: pre-image capture, listing, and conflict-safe revert."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import change_journal as cj


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


def _edit(ws: Path, name: str, old: str | None, new: str) -> Path:
    """Simulate the agent's write path: file ends up with ``new`` on disk."""
    p = ws / name
    p.write_text(new, encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# record + list
# ---------------------------------------------------------------------------

def test_record_and_list_round_trip(home, ws):
    p = _edit(ws, "a.txt", "v1\n", "v2\n")
    r1 = cj.record_change("s1", tool="write_file", path=p,
                          old_text="v1\n", new_text="v2\n")
    r2 = cj.record_change("s1", tool="edit_file", path=p,
                          old_text="v2\n", new_text="v3\n")
    assert r1 is not None and r2 is not None
    assert r1["seq"] == 1 and r2["seq"] == 2
    assert r1["post_hash"] == r2["pre_hash"]  # chained edits share the hash

    records = cj.list_changes("s1")
    assert [r["seq"] for r in records] == [1, 2]  # newest last
    assert records[0]["tool"] == "write_file"
    assert records[0]["created"] is False
    assert records[0]["truncated"] is False
    assert records[0]["path"] == str(p)

    # Pre-image is stored verbatim and content-addressed.
    sdir = home / ".delfin" / "undo" / "s1"
    pre = sdir / records[0]["pre_file"]
    assert pre.read_text(encoding="utf-8") == "v1\n"
    assert records[0]["pre_file"].startswith("000001-")

    assert cj.list_changes("s1", last_n=1) == [records[1]]
    assert cj.list_changes("nope") == []


def test_list_skips_torn_lines(home, ws):
    p = _edit(ws, "a.txt", "v1", "v2")
    cj.record_change("s1", tool="write_file", path=p,
                     old_text="v1", new_text="v2")
    journal = home / ".delfin" / "undo" / "s1" / "journal.jsonl"
    with journal.open("a", encoding="utf-8") as fh:
        fh.write('{"seq": 99, "truncat')  # torn mid-append
    assert [r["seq"] for r in cj.list_changes("s1")] == [1]


# ---------------------------------------------------------------------------
# revert: last / created / conflicts
# ---------------------------------------------------------------------------

def test_revert_last_restores_exact_bytes(home, ws):
    original = "line1\nline2\r\n\ttab\n"   # mixed newlines on purpose
    p = _edit(ws, "a.txt", original, "clobbered\n")
    cj.record_change("s1", tool="write_file", path=p,
                     old_text=original, new_text="clobbered\n")

    res = cj.revert("s1", scope="last")
    assert res["reverted"] == [str(p)]
    assert res["conflicts"] == [] and res["skipped"] == []
    # Raw-byte comparison: read_text would collapse \r\n and mask a
    # newline-translation bug in the restore path.
    assert p.read_bytes() == original.encode("utf-8")


def test_created_file_revert_deletes(home, ws):
    p = _edit(ws, "new.txt", None, "fresh\n")
    cj.record_change("s1", tool="write_file", path=p,
                     old_text=None, new_text="fresh\n")

    res = cj.revert("s1", scope="last")
    assert res["reverted"] == [str(p)]
    assert not p.exists()

    # Second revert of the same entry: file already absent → skipped.
    res2 = cj.revert("s1", scope="last")
    assert res2["reverted"] == []
    assert res2["skipped"] and "absent" in res2["skipped"][0]["reason"]


def test_post_hash_mismatch_is_conflict_file_untouched(home, ws):
    p = _edit(ws, "a.txt", "v1", "v2")
    cj.record_change("s1", tool="write_file", path=p,
                     old_text="v1", new_text="v2")
    p.write_text("user edited this afterwards", encoding="utf-8")

    res = cj.revert("s1", scope="last")
    assert res["reverted"] == []
    assert len(res["conflicts"]) == 1
    assert res["conflicts"][0]["path"] == str(p)
    assert "hash mismatch" in res["conflicts"][0]["reason"]
    assert p.read_text(encoding="utf-8") == "user edited this afterwards"


def test_created_file_conflict_when_user_changed_it(home, ws):
    p = _edit(ws, "new.txt", None, "fresh\n")
    cj.record_change("s1", tool="write_file", path=p,
                     old_text=None, new_text="fresh\n")
    p.write_text("user work in the created file", encoding="utf-8")

    res = cj.revert("s1", scope="last")
    assert res["reverted"] == []
    assert len(res["conflicts"]) == 1
    assert p.exists()  # never deleted when the content moved on


# ---------------------------------------------------------------------------
# revert: session chain order + turn scope + workspace guard
# ---------------------------------------------------------------------------

def test_session_revert_unwinds_chain_newest_first(home, ws):
    p = _edit(ws, "a.txt", "v1", "v3")  # disk ends at v3
    cj.record_change("s1", tool="edit_file", path=p,
                     old_text="v1", new_text="v2")
    cj.record_change("s1", tool="edit_file", path=p,
                     old_text="v2", new_text="v3")

    res = cj.revert("s1", scope="session")
    # Both entries revert (newest first: v3→v2, then v2→v1).
    assert res["reverted"] == [str(p), str(p)]
    assert res["conflicts"] == []
    assert p.read_text(encoding="utf-8") == "v1"


def test_turn_scope_reverts_only_given_seqs(home, ws):
    pa = _edit(ws, "a.txt", "a1", "a2")
    pb = _edit(ws, "b.txt", "b1", "b2")
    cj.record_change("s1", tool="edit_file", path=pa,
                     old_text="a1", new_text="a2")
    r2 = cj.record_change("s1", tool="edit_file", path=pb,
                          old_text="b1", new_text="b2")

    res = cj.revert("s1", scope="turn", turn_seqs=[r2["seq"]])
    assert res["reverted"] == [str(pb)]
    assert pb.read_text(encoding="utf-8") == "b1"
    assert pa.read_text(encoding="utf-8") == "a2"  # untouched

    # Empty / unknown seqs → nothing happens.
    assert cj.revert("s1", scope="turn", turn_seqs=[999]) == {
        "reverted": [], "conflicts": [], "skipped": []}


def test_workspace_guard_skips_outside_paths(home, ws, tmp_path):
    outside = tmp_path / "elsewhere"
    outside.mkdir()
    p = _edit(outside, "x.txt", "v1", "v2")
    cj.record_change("s1", tool="write_file", path=p,
                     old_text="v1", new_text="v2")

    res = cj.revert("s1", scope="last", workspace=ws)
    assert res["reverted"] == []
    assert res["skipped"] == [
        {"path": str(p), "reason": "outside workspace"}]
    assert p.read_text(encoding="utf-8") == "v2"

    # Without the workspace restriction the same entry reverts.
    res2 = cj.revert("s1", scope="last")
    assert res2["reverted"] == [str(p)]


def test_unknown_scope_is_reported_not_raised(home):
    res = cj.revert("s1", scope="bogus")
    assert res["reverted"] == [] and res["conflicts"] == []
    assert "unknown scope" in res["skipped"][0]["reason"]


# ---------------------------------------------------------------------------
# size guard, traversal safety, cap pruning, never-raises
# ---------------------------------------------------------------------------

def test_oversize_pre_image_stored_truncated_and_revert_refuses(home, ws):
    big = "x" * (cj.MAX_PRE_IMAGE_BYTES + 100)
    p = _edit(ws, "big.txt", big, "small now")
    rec = cj.record_change("s1", tool="write_file", path=p,
                           old_text=big, new_text="small now")
    assert rec is not None
    entry = cj.list_changes("s1")[-1]
    assert entry["truncated"] is True
    pre = home / ".delfin" / "undo" / "s1" / entry["pre_file"]
    assert pre.stat().st_size <= cj.MAX_PRE_IMAGE_BYTES

    res = cj.revert("s1", scope="last")
    assert res["reverted"] == []
    assert res["skipped"] and "truncated" in res["skipped"][0]["reason"]
    assert p.read_text(encoding="utf-8") == "small now"  # untouched


def test_session_id_traversal_is_neutralized(home, ws):
    p = _edit(ws, "a.txt", "v1", "v2")
    rec = cj.record_change("../../evil", tool="write_file", path=p,
                           old_text="v1", new_text="v2")
    assert rec is not None
    undo_root = home / ".delfin" / "undo"
    # Everything stays under the undo root; nothing escaped upward.
    assert not (home / ".delfin" / "evil").exists()
    assert not (home.parent / "evil").exists()
    subdirs = [d.name for d in undo_root.iterdir()]
    assert len(subdirs) == 1
    import re
    assert re.fullmatch(r"[A-Za-z0-9_-]+", subdirs[0])
    assert "/" not in subdirs[0] and ".." not in subdirs[0]
    assert cj.list_changes("../../evil")[0]["seq"] == 1


def test_cap_pruning_drops_oldest_and_unlinks_pre_images(home, ws, monkeypatch):
    monkeypatch.setattr(cj, "MAX_ENTRIES_PER_SESSION", 3)
    p = ws / "a.txt"
    files = []
    for i in range(5):
        p.write_text(f"v{i + 1}", encoding="utf-8")
        cj.record_change("s1", tool="edit_file", path=p,
                         old_text=f"v{i}", new_text=f"v{i + 1}")
        files.append(cj.list_changes("s1")[-1]["pre_file"])

    records = cj.list_changes("s1")
    assert [r["seq"] for r in records] == [3, 4, 5]
    sdir = home / ".delfin" / "undo" / "s1"
    assert not (sdir / files[0]).exists()
    assert not (sdir / files[1]).exists()
    for kept in files[2:]:
        assert (sdir / kept).exists()
    # Seq numbering continues past the prune.
    r = cj.record_change("s1", tool="edit_file", path=p,
                         old_text="v5", new_text="v6")
    assert r["seq"] == 6


def test_record_never_raises_on_unwritable_store(home, ws, monkeypatch):
    import os

    p = _edit(ws, "a.txt", "v1", "v2")
    original_write = cj._atomic_write

    def boom(path, text):
        raise OSError("disk full")

    monkeypatch.setattr(cj, "_atomic_write", boom)
    assert cj.record_change("s1", tool="write_file", path=p,
                            old_text="v1", new_text="v2") is None
    monkeypatch.setattr(cj, "_atomic_write", original_write)

    # A read-only undo root must not raise either (chmod is advisory for
    # root, so only assert the None contract when it can actually block).
    if os.geteuid() != 0:
        root = home / ".delfin"
        root.mkdir(exist_ok=True)
        (root / "undo").mkdir(exist_ok=True)
        (root / "undo").chmod(0o500)
        try:
            assert cj.record_change("s2", tool="write_file", path=p,
                                    old_text="v1", new_text="v2") is None
        finally:
            (root / "undo").chmod(0o700)
        assert cj.list_changes("s2") == []
