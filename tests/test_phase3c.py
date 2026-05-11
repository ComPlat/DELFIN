"""Tests for Phase 3c: task store, audit log, web tools (safety side)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor
from delfin.agent import audit_log as al
from delfin.agent import web_tools as wt


@pytest.fixture
def workspace(tmp_path) -> Path:
    ws = tmp_path / "ws"
    ws.mkdir()
    return ws


# ---------------------------------------------------------------------------
# Task store
# ---------------------------------------------------------------------------


def test_task_create_then_list(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="default")
    out = json.loads(_doc_executor.execute("task_create", {
        "subject": "Wire up BoTorch", "description": "ax + gp surrogate",
    }, perms))
    assert out["status"] == "created"
    tid = out["task"]["id"]
    out2 = json.loads(_doc_executor.execute(
        "task_list", {}, perms))
    assert out2["count"] == 1
    assert out2["tasks"][0]["id"] == tid


def test_task_update_status(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="default")
    tid = json.loads(_doc_executor.execute("task_create", {
        "subject": "X",
    }, perms))["task"]["id"]
    out = json.loads(_doc_executor.execute("task_update", {
        "task_id": tid, "status": "in_progress",
    }, perms))
    assert out["task"]["status"] == "in_progress"
    out2 = json.loads(_doc_executor.execute("task_update", {
        "task_id": tid, "status": "completed",
    }, perms))
    assert out2["task"]["status"] == "completed"


def test_task_persists_across_reload(workspace):
    """The store reads JSON from disk so two callers see the same state."""
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _doc_executor.execute("task_create", {"subject": "Persist me"}, perms)
    # Drop the cached store to force a re-read.
    from delfin.agent import agent_tasks
    agent_tasks._STORES.clear()
    out = json.loads(_doc_executor.execute("task_list", {}, perms))
    assert out["count"] == 1
    assert out["tasks"][0]["subject"] == "Persist me"


def test_tasks_are_scoped_to_session(workspace):
    perms_a = KitToolPermissions(
        workspace=workspace, mode="default", task_session_id="session-a"
    )
    perms_b = KitToolPermissions(
        workspace=workspace, mode="default", task_session_id="session-b"
    )
    _doc_executor.execute("task_create", {"subject": "Only A"}, perms_a)
    _doc_executor.execute("task_create", {"subject": "Only B"}, perms_b)

    out_a = json.loads(_doc_executor.execute("task_list", {}, perms_a))
    out_b = json.loads(_doc_executor.execute("task_list", {}, perms_b))

    assert out_a["count"] == 1
    assert out_a["tasks"][0]["subject"] == "Only A"
    assert out_b["count"] == 1
    assert out_b["tasks"][0]["subject"] == "Only B"


def test_task_update_unknown_id(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="default")
    out = json.loads(_doc_executor.execute("task_update", {
        "task_id": 999, "status": "completed",
    }, perms))
    assert "error" in out
    assert "not found" in out["error"]


def test_task_create_rejects_empty_subject(workspace):
    perms = KitToolPermissions(workspace=workspace, mode="default")
    out = json.loads(_doc_executor.execute("task_create", {
        "subject": "   ",
    }, perms))
    assert "error" in out


# ---------------------------------------------------------------------------
# Audit log
# ---------------------------------------------------------------------------


def test_audit_writes_jsonl(workspace, tmp_path, monkeypatch):
    log_path = tmp_path / "audit.log"
    monkeypatch.setattr(al, "_default_log_path", lambda: log_path)

    perms = KitToolPermissions(workspace=workspace, mode="default")
    f = workspace / "a.py"
    f.write_text("x = 1\n")
    _doc_executor.execute("read_file", {"path": str(f)}, perms)
    _doc_executor.execute("edit_file", {
        "path": str(f), "old_string": "x = 1", "new_string": "x = 2",
    }, perms)

    assert log_path.exists()
    lines = [json.loads(l) for l in log_path.read_text().splitlines() if l]
    assert any(r["tool"] == "edit_file" and r["decision"] == "ok" for r in lines)


def test_audit_logs_denials(workspace, tmp_path, monkeypatch):
    log_path = tmp_path / "audit.log"
    monkeypatch.setattr(al, "_default_log_path", lambda: log_path)
    perms = KitToolPermissions(workspace=workspace, mode="default")
    _doc_executor.execute("bash", {
        "command": "rm -rf /tmp/fake",
        "description": "evil",
    }, perms)
    lines = [json.loads(l) for l in log_path.read_text().splitlines() if l]
    assert any(r["tool"] == "bash" and r["decision"] == "denied" for r in lines)


def test_audit_rotates_weekly(tmp_path, monkeypatch):
    """Old logs from a previous ISO week get renamed to audit-YYYY-Www.log."""
    log_path = tmp_path / "audit.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    # Manually write an old record with last-year timestamp.
    old = {"ts": "2020-01-01T00:00:00Z", "tool": "edit_file",
           "decision": "ok", "session_id": "", "mode": "default"}
    log_path.write_text(json.dumps(old) + "\n")

    monkeypatch.setattr(al, "_default_log_path", lambda: log_path)
    al.append({"tool": "x", "decision": "ok", "mode": "default",
               "session_id": "", "pid": 1})
    # Old log should be renamed.
    rotated = list(tmp_path.glob("audit-2020*.log"))
    assert rotated, f"no rotated log; dir contents = {list(tmp_path.iterdir())}"


# ---------------------------------------------------------------------------
# Web tools — safety only (no live network needed)
# ---------------------------------------------------------------------------


def test_web_fetch_blocks_localhost():
    out = wt.web_fetch("http://localhost:8080/")
    assert "error" in out
    assert "deny-list" in out["error"]


def test_web_fetch_blocks_metadata():
    out = wt.web_fetch("http://169.254.169.254/latest/meta-data/")
    assert "error" in out


def test_web_fetch_blocks_internal_tld():
    out = wt.web_fetch("http://kit.internal/secrets")
    assert "error" in out


def test_web_fetch_blocks_rfc1918():
    out = wt.web_fetch("http://192.168.1.1/")
    assert "error" in out


def test_web_fetch_rejects_non_http():
    out = wt.web_fetch("file:///etc/passwd")
    assert "error" in out
    assert "http(s)" in out["error"]


def test_web_search_rejects_empty_query():
    out = wt.web_search("")
    assert "error" in out
