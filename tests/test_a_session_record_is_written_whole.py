"""The open-sessions list and a presence record are written whole, and a
session open on another login node is not reopened here.

Review 2026-09-16 of the session list: both files were written with
truncate-then-write, so a reader on another node (the home directory is
shared) could see an empty file -- the whole list of open sessions gone;
and the list was not host-scoped, so node B reopened node A's conversation
and both saved over it.
"""
import json
import os

from delfin.agent import state_paths as SP
from delfin.dashboard import agent_sessions as AS


def test_an_atomic_write_leaves_the_file_whole_and_owner_only(tmp_path):
    target = tmp_path / "state.json"
    target.write_text("old")
    SP.write_text_atomic(target, "new")
    assert target.read_text() == "new"
    assert not [p for p in tmp_path.iterdir() if p.name.startswith(".state.json.")]
    assert oct(target.stat().st_mode & 0o777) == "0o600"


def test_a_failed_write_keeps_the_old_file(tmp_path, monkeypatch):
    target = tmp_path / "state.json"
    target.write_text("old")

    def boom(*a, **k):
        raise OSError("disk")
    monkeypatch.setattr(os, "replace", boom)
    try:
        SP.write_text_atomic(target, "new")
    except OSError:
        pass
    assert target.read_text() == "old"
    assert not [p for p in tmp_path.iterdir() if p.name.startswith(".state.json.")]


def test_open_sessions_carry_their_host_and_only_this_hosts_come_back(tmp_path, monkeypatch):
    monkeypatch.setattr(AS, "_OPEN_SESSIONS_PATH", tmp_path / "open.json")
    monkeypatch.setattr(AS, "_host", lambda: "node-a")
    AS.save_open_sessions([{"session_id": "s1", "workspace": "/w"}])
    data = json.loads((tmp_path / "open.json").read_text())
    assert data["sessions"][0]["host"] == "node-a"
    data["sessions"].append({"session_id": "s2", "workspace": "/w", "host": "node-b"})
    data["sessions"].append({"session_id": "s3", "workspace": "/w"})
    (tmp_path / "open.json").write_text(json.dumps(data))
    assert [r["session_id"] for r in AS.load_open_sessions()] == ["s1", "s3"]


def test_the_presence_record_is_written_whole():
    import inspect
    from delfin.agent import session_presence as P
    assert "write_text_atomic(" in inspect.getsource(P.announce)
