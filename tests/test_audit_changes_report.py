"""Audit-log changes report: session stamping + aggregation + formatting.

"What did you change?" must be answered from the audit record, not from
model memory. These tests pin the read side introduced for that:

* append(..., session_id=...) stamps records that carry no session
  (and never overwrites one they already carry);
* build_changes_report() groups writes per path, lists commands with
  their cwd, includes denied actions, filters by session, and skips
  torn/malformed lines without raising;
* format_changes_report() renders paths, counts and commands for chat.
"""

import json

from delfin.agent import audit_log


def _log(tmp_path):
    return tmp_path / "audit.log"


# ---------------------------------------------------------------------------
# Session stamping
# ---------------------------------------------------------------------------

def test_append_stamps_session_id(tmp_path):
    p = _log(tmp_path)
    rec = audit_log.make_record(tool="write_file", decision="ok",
                                path="a.py")
    rec.pop("session_id", None)          # simulate a legacy caller
    audit_log.append(rec, log_path=p, session_id="sess-42")
    out = audit_log.read_last_n(10, log_path=p)
    assert len(out) == 1
    assert out[0]["session_id"] == "sess-42"


def test_append_keeps_existing_session_id(tmp_path):
    p = _log(tmp_path)
    rec = audit_log.make_record(tool="bash", decision="ok",
                                command="ls", session_id="original")
    audit_log.append(rec, log_path=p, session_id="other")
    out = audit_log.read_last_n(10, log_path=p)
    assert out[0]["session_id"] == "original"


def test_make_record_roundtrip_session_id(tmp_path):
    p = _log(tmp_path)
    rec = audit_log.make_record(tool="edit_file", decision="ok",
                                path="b.py", session_id="sess-rt")
    audit_log.append(rec, log_path=p)
    out = audit_log.read_last_n(10, log_path=p)
    assert out[0]["session_id"] == "sess-rt"
    assert out[0]["tool"] == "edit_file"
    assert out[0]["path"] == "b.py"


# ---------------------------------------------------------------------------
# Report aggregation
# ---------------------------------------------------------------------------

def _seed_session(p, sid="sess-1"):
    """Two writes + one edit on the same file, one bash, one denied write."""
    for tool in ("write_file", "write_file", "edit_file"):
        audit_log.append(
            audit_log.make_record(tool=tool, decision="ok",
                                  path="delfin/core.py", session_id=sid),
            log_path=p)
    audit_log.append(
        audit_log.make_record(tool="bash", decision="ok",
                              command="pytest tests/ -q", session_id=sid,
                              extra={"cwd": "tests"}),
        log_path=p)
    audit_log.append(
        audit_log.make_record(tool="write_file", decision="denied",
                              path="/etc/passwd", session_id=sid),
        log_path=p)


def test_report_groups_writes_per_path(tmp_path):
    p = _log(tmp_path)
    _seed_session(p)
    rep = audit_log.build_changes_report("sess-1", log_path=p)
    assert len(rep["files_written"]) == 1
    entry = rep["files_written"][0]
    assert entry["path"] == "delfin/core.py"
    assert entry["count"] == 3
    assert "write_file" in entry["tool"] and "edit_file" in entry["tool"]


def test_report_lists_commands_with_cwd(tmp_path):
    p = _log(tmp_path)
    _seed_session(p)
    rep = audit_log.build_changes_report("sess-1", log_path=p)
    assert rep["commands"] == [
        {"command": "pytest tests/ -q", "cwd": "tests"}]


def test_report_includes_denied(tmp_path):
    p = _log(tmp_path)
    _seed_session(p)
    rep = audit_log.build_changes_report("sess-1", log_path=p)
    assert len(rep["denied"]) == 1
    assert rep["denied"][0]["tool"] == "write_file"
    assert rep["denied"][0]["target"] == "/etc/passwd"
    # A denied write must not count as a file written.
    assert all(f["path"] != "/etc/passwd" for f in rep["files_written"])


def test_report_includes_permissions_persisted(tmp_path):
    p = _log(tmp_path)
    audit_log.append(
        audit_log.make_record(
            tool="remember_permission", decision="ok", session_id="s",
            persistence={"pattern": "pytest *", "scope": "workspace"}),
        log_path=p)
    rep = audit_log.build_changes_report("s", log_path=p)
    assert rep["permissions_persisted"][0]["tool"] == "remember_permission"
    assert rep["permissions_persisted"][0]["persistence"]["pattern"] == "pytest *"


def test_report_window_counts_records(tmp_path):
    p = _log(tmp_path)
    _seed_session(p)
    rep = audit_log.build_changes_report("sess-1", log_path=p)
    assert rep["window"]["records"] == 5
    assert rep["window"]["from_ts"] != ""
    assert rep["window"]["to_ts"] != ""


# ---------------------------------------------------------------------------
# Session filter / windowing
# ---------------------------------------------------------------------------

def test_report_filters_by_session(tmp_path):
    p = _log(tmp_path)
    _seed_session(p, sid="sess-a")
    audit_log.append(
        audit_log.make_record(tool="write_file", decision="ok",
                              path="other.py", session_id="sess-b"),
        log_path=p)
    rep = audit_log.build_changes_report("sess-b", log_path=p)
    assert [f["path"] for f in rep["files_written"]] == ["other.py"]
    assert rep["commands"] == []
    assert rep["denied"] == []
    assert rep["window"]["records"] == 1


def test_report_no_session_uses_all_records(tmp_path):
    p = _log(tmp_path)
    _seed_session(p, sid="sess-a")
    audit_log.append(
        audit_log.make_record(tool="write_file", decision="ok",
                              path="other.py", session_id="sess-b"),
        log_path=p)
    rep = audit_log.build_changes_report(log_path=p)
    paths = {f["path"] for f in rep["files_written"]}
    assert paths == {"delfin/core.py", "other.py"}


def test_report_last_n_caps_records(tmp_path):
    p = _log(tmp_path)
    for i in range(10):
        audit_log.append(
            audit_log.make_record(tool="write_file", decision="ok",
                                  path=f"f{i}.py", session_id="s"),
            log_path=p)
    rep = audit_log.build_changes_report("s", last_n=3, log_path=p)
    assert rep["window"]["records"] == 3
    assert [f["path"] for f in rep["files_written"]] == [
        "f7.py", "f8.py", "f9.py"]


def test_report_since_ts_filter(tmp_path):
    p = _log(tmp_path)
    _seed_session(p)
    # since_ts far in the future excludes everything; far past keeps all.
    rep = audit_log.build_changes_report(
        since_ts="2999-01-01T00:00:00Z", log_path=p)
    assert rep["window"]["records"] == 0
    rep = audit_log.build_changes_report(
        since_ts="2000-01-01T00:00:00Z", log_path=p)
    assert rep["window"]["records"] == 5


# ---------------------------------------------------------------------------
# Robustness
# ---------------------------------------------------------------------------

def test_malformed_lines_skipped(tmp_path):
    p = _log(tmp_path)
    _seed_session(p)
    with p.open("a", encoding="utf-8") as fh:
        fh.write('{"ts": "2026-07-29T10:00:00Z", "tool": "wri')  # torn line
        fh.write("\nnot json at all\n")
    audit_log.append(
        audit_log.make_record(tool="write_file", decision="ok",
                              path="after.py", session_id="sess-1"),
        log_path=p)
    rep = audit_log.build_changes_report("sess-1", log_path=p)
    paths = {f["path"] for f in rep["files_written"]}
    assert "after.py" in paths and "delfin/core.py" in paths
    assert rep["window"]["records"] == 6      # 5 seeded + 1 after the tear


def test_report_never_raises_on_missing_log(tmp_path):
    rep = audit_log.build_changes_report(
        "any", log_path=tmp_path / "does-not-exist.log")
    assert rep["files_written"] == []
    assert rep["window"]["records"] == 0


def test_report_never_raises_on_binary_garbage(tmp_path):
    p = _log(tmp_path)
    p.write_bytes(b"\x00\xff\xfe garbage \x9c\n\x01\n")
    rep = audit_log.build_changes_report(log_path=p)
    assert rep["window"]["records"] == 0


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

def test_format_contains_paths_counts_and_commands(tmp_path):
    p = _log(tmp_path)
    _seed_session(p)
    text = audit_log.format_changes_report(
        audit_log.build_changes_report("sess-1", log_path=p))
    assert "delfin/core.py" in text
    assert "×3" in text
    assert "pytest tests/ -q" in text
    assert "cwd: tests" in text
    assert "/etc/passwd" in text
    assert "Files written:" in text
    assert "Commands run:" in text


def test_format_empty_report(tmp_path):
    text = audit_log.format_changes_report(
        audit_log.build_changes_report(
            "nobody", log_path=tmp_path / "missing.log"))
    assert "No recorded changes" in text


def test_format_never_raises_on_junk():
    assert isinstance(audit_log.format_changes_report({}), str)
    assert isinstance(
        audit_log.format_changes_report({"files_written": "junk"}), str)
    assert isinstance(audit_log.format_changes_report(None), str)


def test_report_is_json_serializable(tmp_path):
    p = _log(tmp_path)
    _seed_session(p)
    rep = audit_log.build_changes_report("sess-1", log_path=p)
    json.dumps(rep)
