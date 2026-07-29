"""Pending-edit store (diff-approval mode): staging, conflict-safe
approval, rejection, caps, sanitization, and dashboard rendering."""

from __future__ import annotations

import json
import os
import stat
from pathlib import Path

import pytest

from delfin.agent import change_journal as cj
from delfin.agent import pending_changes as pc


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


def _stage_edit(sid, path, old, new, tool="edit_file", note=""):
    return pc.stage(sid, tool=tool, path=path,
                    old_text=old, new_text=new, note=note)


# ---------------------------------------------------------------------------
# stage
# ---------------------------------------------------------------------------

def test_stage_returns_record_with_diff(home, ws):
    p = ws / "a.txt"
    p.write_text("v1\n", encoding="utf-8")
    rec = _stage_edit("s1", p, "v1\n", "v2\n", note="bump")
    assert "error" not in rec
    assert rec["id"] == "1" and rec["seq"] == 1
    assert rec["status"] == "pending"
    assert rec["created"] is False
    assert rec["created_at"]
    assert rec["note"] == "bump"
    assert "-v1" in rec["diff"] and "+v2" in rec["diff"]
    # File on disk is untouched — staging never applies.
    assert p.read_text(encoding="utf-8") == "v1\n"


def test_stage_create_file_marks_created(home, ws):
    p = ws / "new.txt"
    rec = _stage_edit("s1", p, None, "hello\n", tool="write_file")
    assert "error" not in rec
    assert rec["created"] is True
    assert rec["old_file"] == ""
    assert "+hello" in rec["diff"]
    assert not p.exists()


def test_stage_refuses_identical_content(home, ws):
    rec = _stage_edit("s1", ws / "a.txt", "same\n", "same\n")
    assert "error" in rec


def test_stage_refuses_second_pending_for_same_path(home, ws):
    p = ws / "a.txt"
    assert "error" not in _stage_edit("s1", p, "v1\n", "v2\n")
    rec = _stage_edit("s1", p, "v1\n", "v3\n")
    assert "error" in rec
    assert "already pending" in rec["error"]
    # A different path still stages fine.
    assert "error" not in _stage_edit("s1", ws / "b.txt", None, "x\n")


def test_stage_refuses_oversize_image(home, ws, monkeypatch):
    monkeypatch.setattr(pc, "MAX_IMAGE_BYTES", 16)
    rec = _stage_edit("s1", ws / "a.txt", "v1\n", "x" * 64)
    assert "error" in rec


def test_stage_refuses_when_queue_full(home, ws, monkeypatch):
    monkeypatch.setattr(pc, "MAX_PENDING_PER_SESSION", 2)
    assert "error" not in _stage_edit("s1", ws / "a.txt", None, "a\n")
    assert "error" not in _stage_edit("s1", ws / "b.txt", None, "b\n")
    rec = _stage_edit("s1", ws / "c.txt", None, "c\n")
    assert "error" in rec and "full" in rec["error"]


def test_stage_prunes_only_terminal_records(home, ws, monkeypatch):
    monkeypatch.setattr(pc, "MAX_TERMINAL_RECORDS", 1)
    _stage_edit("s1", ws / "a.txt", None, "a\n")
    _stage_edit("s1", ws / "b.txt", None, "b\n")
    assert pc.reject("s1", "1")["status"] == "rejected"
    assert pc.reject("s1", "2")["status"] == "rejected"
    # Next stage triggers the prune of the oldest terminal record.
    _stage_edit("s1", ws / "c.txt", None, "c\n")
    ids = {r["seq"] for r in pc._read_journal("s1")}
    assert 1 not in ids            # oldest rejected pruned
    assert {2, 3} <= ids           # newer terminal + pending survive


# ---------------------------------------------------------------------------
# list_pending / get
# ---------------------------------------------------------------------------

def test_list_pending_filters_and_orders(home, ws):
    _stage_edit("s1", ws / "a.txt", None, "a\n")
    _stage_edit("s1", ws / "b.txt", None, "b\n")
    assert [r["id"] for r in pc.list_pending("s1")] == ["1", "2"]
    pc.reject("s1", "1")
    assert [r["id"] for r in pc.list_pending("s1")] == ["2"]


def test_get_accepts_int_and_padded_ids(home, ws):
    _stage_edit("s1", ws / "a.txt", None, "a\n")
    assert pc.get("s1", "1")["id"] == "1"
    assert pc.get("s1", 1)["id"] == "1"
    assert pc.get("s1", "000001")["id"] == "1"
    assert pc.get("s1", "999") is None
    assert pc.get("other", "1") is None


# ---------------------------------------------------------------------------
# approve
# ---------------------------------------------------------------------------

def test_approve_applies_edit_and_records_undo(home, ws):
    p = ws / "a.txt"
    p.write_text("v1\n", encoding="utf-8")
    _stage_edit("s1", p, "v1\n", "v2\n")
    res = pc.approve("s1", "1", workspace=ws)
    assert res["status"] == "applied"
    assert p.read_text(encoding="utf-8") == "v2\n"
    # Undo-journal hand-off happened inside approve().
    undo = cj.list_changes("s1")
    assert len(undo) == 1
    assert undo[0]["path"] == str(p)
    assert res["undo_seq"] == undo[0]["seq"]
    # Record is terminal, images gone, queue empty.
    assert pc.get("s1", "1")["status"] == "applied"
    assert pc.list_pending("s1") == []
    sdir = home / ".delfin" / "pending" / "s1"
    assert not list(sdir.glob("*-old.txt")) and not list(sdir.glob("*-new.txt"))


def test_approve_creates_new_file_with_parents(home, ws):
    p = ws / "sub" / "dir" / "new.txt"
    _stage_edit("s1", p, None, "made\n", tool="write_file")
    res = pc.approve("s1", "1", workspace=ws)
    assert res["status"] == "applied"
    assert p.read_text(encoding="utf-8") == "made\n"
    # Created entries revert by deletion via the undo journal.
    assert cj.list_changes("s1")[0]["created"] is True


def test_approve_conflict_when_content_changed(home, ws):
    p = ws / "a.txt"
    p.write_text("v1\n", encoding="utf-8")
    _stage_edit("s1", p, "v1\n", "v2\n")
    p.write_text("user-edit\n", encoding="utf-8")
    res = pc.approve("s1", "1", workspace=ws)
    assert res["status"] == "conflict"
    assert "changed since staging" in res["reason"]
    assert p.read_text(encoding="utf-8") == "user-edit\n"   # untouched
    # Stays pending, conflict reason persisted for the dashboard.
    rec = pc.get("s1", "1")
    assert rec["status"] == "pending"
    assert "changed since staging" in rec["last_conflict"]
    assert cj.list_changes("s1") == []                       # no undo entry


def test_approve_conflict_when_created_file_appeared(home, ws):
    p = ws / "new.txt"
    _stage_edit("s1", p, None, "agent\n")
    p.write_text("user\n", encoding="utf-8")
    res = pc.approve("s1", "1", workspace=ws)
    assert res["status"] == "conflict"
    assert p.read_text(encoding="utf-8") == "user\n"


def test_approve_conflict_when_file_missing(home, ws):
    p = ws / "a.txt"
    p.write_text("v1\n", encoding="utf-8")
    _stage_edit("s1", p, "v1\n", "v2\n")
    p.unlink()
    res = pc.approve("s1", "1", workspace=ws)
    assert res["status"] == "conflict"
    assert "missing" in res["reason"]


def test_approve_refuses_target_outside_workspace(home, ws, tmp_path):
    outside = tmp_path / "elsewhere.txt"
    outside.write_text("v1\n", encoding="utf-8")
    _stage_edit("s1", outside, "v1\n", "v2\n")
    res = pc.approve("s1", "1", workspace=ws)
    assert res["status"] == "conflict"
    assert "outside the workspace" in res["reason"]
    assert outside.read_text(encoding="utf-8") == "v1\n"


def test_approve_preserves_permission_bits(home, ws):
    p = ws / "a.sh"
    p.write_text("v1\n", encoding="utf-8")
    os.chmod(p, 0o750)
    _stage_edit("s1", p, "v1\n", "v2\n")
    assert pc.approve("s1", "1", workspace=ws)["status"] == "applied"
    assert stat.S_IMODE(p.stat().st_mode) == 0o750
    assert p.read_text(encoding="utf-8") == "v2\n"


def test_approve_not_found_and_double_approve(home, ws):
    assert pc.approve("s1", "7", workspace=ws)["status"] == "not_found"
    p = ws / "a.txt"
    p.write_text("v1\n", encoding="utf-8")
    _stage_edit("s1", p, "v1\n", "v2\n")
    assert pc.approve("s1", "1", workspace=ws)["status"] == "applied"
    assert pc.approve("s1", "1", workspace=ws)["status"] == "already_applied"


# ---------------------------------------------------------------------------
# reject / bulk operations
# ---------------------------------------------------------------------------

def test_reject_marks_record_and_removes_images(home, ws):
    p = ws / "a.txt"
    p.write_text("v1\n", encoding="utf-8")
    _stage_edit("s1", p, "v1\n", "v2\n")
    res = pc.reject("s1", "1")
    assert res["status"] == "rejected"
    assert p.read_text(encoding="utf-8") == "v1\n"
    assert pc.get("s1", "1")["status"] == "rejected"
    assert pc.reject("s1", "1")["status"] == "already_rejected"
    sdir = home / ".delfin" / "pending" / "s1"
    assert not list(sdir.glob("*-old.txt")) and not list(sdir.glob("*-new.txt"))


def test_approve_all_reports_applied_and_conflicts(home, ws):
    a = ws / "a.txt"
    a.write_text("v1\n", encoding="utf-8")
    _stage_edit("s1", a, "v1\n", "v2\n")
    b = ws / "b.txt"
    b.write_text("v1\n", encoding="utf-8")
    _stage_edit("s1", b, "v1\n", "v2\n")
    b.write_text("user\n", encoding="utf-8")     # provoke a conflict on b
    res = pc.approve_all("s1", workspace=ws)
    assert [r["path"] for r in res["applied"]] == [str(a)]
    assert [r["path"] for r in res["conflicts"]] == [str(b)]
    assert res["errors"] == []
    assert a.read_text(encoding="utf-8") == "v2\n"
    assert b.read_text(encoding="utf-8") == "user\n"
    assert [r["id"] for r in pc.list_pending("s1")] == ["2"]  # conflict stays


def test_reject_all(home, ws):
    _stage_edit("s1", ws / "a.txt", None, "a\n")
    _stage_edit("s1", ws / "b.txt", None, "b\n")
    res = pc.reject_all("s1")
    assert res["rejected"] == ["1", "2"]
    assert pc.list_pending("s1") == []


# ---------------------------------------------------------------------------
# store hygiene: sanitized sids, permissions, never-raises
# ---------------------------------------------------------------------------

def test_session_id_is_sanitized_no_traversal(home, ws):
    rec = _stage_edit("../../etc", ws / "a.txt", None, "a\n")
    assert "error" not in rec
    root = home / ".delfin" / "pending"
    (sdir,) = list(root.iterdir())
    assert sdir.name == "______etc"
    assert ".." not in sdir.name and "/" not in sdir.name
    # The same raw sid resolves to the same store on read.
    assert [r["id"] for r in pc.list_pending("../../etc")] == ["1"]


def test_store_files_are_private(home, ws):
    _stage_edit("s1", ws / "a.txt", "v1\n", "v2\n")
    sdir = home / ".delfin" / "pending" / "s1"
    assert stat.S_IMODE(sdir.stat().st_mode) == 0o700
    for f in sdir.iterdir():
        assert stat.S_IMODE(f.stat().st_mode) == 0o600


def test_torn_journal_lines_are_skipped(home, ws):
    _stage_edit("s1", ws / "a.txt", None, "a\n")
    journal = home / ".delfin" / "pending" / "s1" / "pending.jsonl"
    with journal.open("a", encoding="utf-8") as fh:
        fh.write("{torn json\n\nnot json at all\n")
        fh.write(json.dumps({"no_seq": True}) + "\n")
    assert [r["id"] for r in pc.list_pending("s1")] == ["1"]
    # A later stage continues the seq chain past the garbage.
    rec = _stage_edit("s1", ws / "b.txt", None, "b\n")
    assert rec["seq"] == 2


def test_public_api_never_raises_without_store(home, ws):
    # No session dir exists for any of these.
    assert pc.list_pending("ghost") == []
    assert pc.get("ghost", "1") is None
    assert pc.approve("ghost", "1", workspace=ws)["status"] == "not_found"
    assert pc.reject("ghost", "1")["status"] == "not_found"
    assert pc.approve_all("ghost", workspace=ws) == {
        "applied": [], "conflicts": [], "errors": []}
    assert pc.reject_all("ghost") == {"rejected": [], "errors": []}
    assert isinstance(pc.render_pending("ghost"), str)


def test_stage_reports_error_instead_of_raising(home, ws, monkeypatch):
    # Force a filesystem failure inside stage(): journal append fails.
    monkeypatch.setattr(pc, "_journal_path",
                        lambda sid: home / "missing" / "pending.jsonl")
    rec = _stage_edit("s1", ws / "a.txt", None, "a\n")
    assert "error" in rec


# ---------------------------------------------------------------------------
# render_pending
# ---------------------------------------------------------------------------

def test_render_pending_empty(home):
    out = pc.render_pending("s1")
    assert "No pending changes" in out


def test_render_pending_compact_markdown(home, ws):
    p = ws / "a.txt"
    p.write_text("v1\n", encoding="utf-8")
    _stage_edit("s1", p, "v1\n", "v2\n", note="tweak")
    _stage_edit("s1", ws / "new.txt", None, "x\n", tool="write_file")
    out = pc.render_pending("s1")
    assert "2 awaiting approval" in out
    assert "**#1**" in out and "**#2**" in out
    assert "edit_file (edit)" in out
    assert "write_file (new file)" in out
    assert "```diff" in out
    assert "note: tweak" in out
    assert "/approve <id>" in out and "/reject <id>" in out


def test_render_pending_truncates_long_queue_and_diffs(home, ws):
    for i in range(4):
        _stage_edit("s1", ws / f"f{i}.txt", None,
                    "\n".join(f"line{j}" for j in range(40)) + "\n")
    out = pc.render_pending("s1", max_items=2, max_diff_lines=5)
    assert "... and 2 more" in out
    assert "more diff lines)" in out


def test_render_pending_shows_last_conflict(home, ws):
    p = ws / "a.txt"
    p.write_text("v1\n", encoding="utf-8")
    _stage_edit("s1", p, "v1\n", "v2\n")
    p.write_text("user\n", encoding="utf-8")
    assert pc.approve("s1", "1", workspace=ws)["status"] == "conflict"
    out = pc.render_pending("s1")
    assert "last conflict:" in out
