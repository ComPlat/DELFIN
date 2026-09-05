"""An undo overwrites files too — so it has to be recorded like one.

Reproduced before the fix: ``revert`` wrote over the user's file and
deleted created files without journalling, auditing or tracing any of
it. Running undo twice then blamed the user — "content changed since the
agent's edit" — for the change the agent's own undo had made, and the
version it had replaced existed nowhere.

The other half of the same gap: entries pruned at the per-session cap
were forgotten outright, so undoing one returned three empty lists,
byte-identical to "there was nothing to undo". The user asked for their
file back and was told, in the same words as a no-op, that nothing
happened.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import audit_log as al
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


def _record(sid: str, path: Path, old: str, new: str) -> dict:
    """The write path's order: write first, journal the pre-image after."""
    path.write_text(new, encoding="utf-8")
    return cj.record_change(sid, tool="write_file", path=path,
                            old_text=old, new_text=new)


# ---------------------------------------------------------------------------
# The undo is journalled, and can itself be undone
# ---------------------------------------------------------------------------

def test_the_undo_is_written_back_to_the_journal(home, ws):
    p = ws / "a.txt"
    _record("s1", p, "user\n", "agent\n")

    res = cj.revert("s1", scope="last", workspace=ws)
    assert res["reverted"] == [str(p)]
    assert p.read_text(encoding="utf-8") == "user\n"

    records = cj.list_changes("s1")
    assert [r["tool"] for r in records] == ["write_file", "undo"]
    assert records[-1]["undo_of"] == records[0]["seq"]
    assert records[0]["undone"]["by_seq"] == records[-1]["seq"]


def test_the_agents_version_is_not_lost_by_the_undo(home, ws):
    p = ws / "a.txt"
    _record("s2", p, "user\n", "agent\n")
    res = cj.revert("s2", scope="last", workspace=ws)

    back = cj.revert("s2", scope="turn", turn_seqs=res["undo_seqs"],
                     workspace=ws)
    assert back["reverted"] == [str(p)], back
    assert p.read_text(encoding="utf-8") == "agent\n", (
        "the version the undo replaced existed nowhere")


def test_a_second_undo_does_not_blame_the_user(home, ws):
    p = ws / "a.txt"
    _record("s3", p, "user\n", "agent\n")
    cj.revert("s3", scope="last", workspace=ws)

    again = cj.revert("s3", scope="last", workspace=ws)
    assert again["conflicts"] == [], (
        "the user was accused of the agent's own undo")
    assert again["reverted"] == []
    assert "already been undone" in again["skipped"][0]["reason"]
    assert p.read_text(encoding="utf-8") == "user\n"


def test_repeated_undo_walks_back_through_the_session(home, ws):
    a, b = ws / "a.txt", ws / "b.txt"
    _record("s4", a, "a-user\n", "a-agent\n")
    _record("s4", b, "b-user\n", "b-agent\n")

    cj.revert("s4", scope="last", workspace=ws)
    assert b.read_text(encoding="utf-8") == "b-user\n"
    assert a.read_text(encoding="utf-8") == "a-agent\n"

    cj.revert("s4", scope="last", workspace=ws)
    assert a.read_text(encoding="utf-8") == "a-user\n"


def test_a_session_undo_does_not_re_apply_its_own_undo(home, ws):
    p = ws / "a.txt"
    _record("s5", p, "user\n", "agent\n")
    cj.revert("s5", scope="session", workspace=ws)
    assert p.read_text(encoding="utf-8") == "user\n"

    cj.revert("s5", scope="session", workspace=ws)
    assert p.read_text(encoding="utf-8") == "user\n", (
        "undoing the session a second time put the agent's version back")


def test_the_undo_tool_records_the_undo_in_the_audit_log(home, ws):
    p = ws / "a.txt"
    ex = _DocToolExecutor()
    perms = KitToolPermissions(workspace=ws, mode="default")
    perms.task_session_id = "s6"
    perms.read_tracker[str(p)] = 0.0
    p.write_text("user\n", encoding="utf-8")
    perms.read_tracker[str(p)] = p.stat().st_mtime
    ex._execute_write_file({"path": "a.txt", "content": "agent\n"}, perms)

    out = json.loads(ex._execute_undo_changes({"scope": "last"}, perms))
    assert out["reverted"] == [str(p)], out

    log = home / ".delfin" / "audit.log"
    records = al.read_last_n(20, log_path=log)
    undos = [r for r in records if r.get("tool") == "undo_changes"]
    assert undos, "the one write nobody asked for twice has no audit record"
    assert undos[-1]["path"] == str(p)


def test_the_turn_cursor_is_kept_per_session(home, ws):
    """One executor serves every parallel sub-agent; a single flat list
    let a turn-scoped undo revert a sibling's writes."""
    ex = _DocToolExecutor()
    out = []
    for sid, name in (("sess-a", "a.txt"), ("sess-b", "b.txt")):
        perms = KitToolPermissions(workspace=ws, mode="default")
        perms.task_session_id = sid
        target = ws / name
        target.write_text("user\n", encoding="utf-8")
        perms.read_tracker[str(target)] = target.stat().st_mtime
        ex._execute_write_file({"path": name, "content": "agent\n"}, perms)
        out.append((sid, target))

    a_seqs = ex._turn_seqs_for("sess-a")
    b_seqs = ex._turn_seqs_for("sess-b")
    assert a_seqs and b_seqs
    assert len(a_seqs) == 1 and len(b_seqs) == 1

    perms_b = KitToolPermissions(workspace=ws, mode="default")
    perms_b.task_session_id = "sess-b"
    json.loads(ex._execute_undo_changes({"scope": "turn"}, perms_b))
    assert (ws / "b.txt").read_text(encoding="utf-8") == "user\n"
    assert (ws / "a.txt").read_text(encoding="utf-8") == "agent\n", (
        "the turn undo reached into a sibling session's files")


# ---------------------------------------------------------------------------
# The cap: forgetting is not the same as nothing to undo
# ---------------------------------------------------------------------------

def test_a_pruned_entry_says_the_pre_image_was_dropped(home, ws, monkeypatch):
    monkeypatch.setattr(cj, "MAX_ENTRIES_PER_SESSION", 2)
    p = ws / "a.txt"
    first = _record("s7", p, "v0\n", "v1\n")
    for i in range(1, 4):
        _record("s7", p, f"v{i}\n", f"v{i + 1}\n")

    res = cj.revert("s7", scope="turn", turn_seqs=[first["seq"]], workspace=ws)
    assert res["reverted"] == []
    assert "dropped at the" in res["skipped"][0]["reason"]

    nothing = cj.revert("s7", scope="turn", turn_seqs=[9999], workspace=ws)
    assert res["skipped"] != nothing["skipped"], (
        "a dropped pre-image reads exactly like 'there was nothing to undo'")


def test_the_record_that_pushed_an_entry_out_says_so(home, ws, monkeypatch):
    monkeypatch.setattr(cj, "MAX_ENTRIES_PER_SESSION", 2)
    p = ws / "a.txt"
    rec = None
    for i in range(3):
        rec = _record("s8", p, f"v{i}\n", f"v{i + 1}\n")
    assert rec["dropped"] == [1], rec


def test_an_ordinary_record_reports_no_drop(home, ws):
    p = ws / "a.txt"
    rec = _record("s9", p, "v0\n", "v1\n")
    assert rec["dropped"] == []
