"""Tests for session resume helpers (resume_latest, cleanup_old_sessions, find_sessions)."""

from __future__ import annotations

import os
import tempfile
import time
from pathlib import Path

import pytest

from delfin.agent import session_store as SS


@pytest.fixture
def tmp_sessions_dir(monkeypatch):
    with tempfile.TemporaryDirectory() as d:
        path = Path(d)
        monkeypatch.setattr(SS, "_SESSIONS_DIR", path)
        yield path


def _save(session_id: str, *, mode="quick", title="t", msg="hi") -> None:
    SS.save_session(
        session_id=session_id,
        mode=mode,
        title=title,
        chat_messages=[{"role": "user", "content": msg}],
    )


def test_resume_latest_returns_newest(tmp_sessions_dir):
    _save("a", title="Alpha")
    time.sleep(0.05)
    _save("b", title="Beta")
    time.sleep(0.05)
    _save("c", title="Gamma")
    latest = SS.resume_latest()
    assert latest is not None
    assert latest["session_id"] == "c"


def test_resume_latest_none_when_empty(tmp_sessions_dir):
    assert SS.resume_latest() is None


def test_resume_latest_respects_max_age(tmp_sessions_dir):
    _save("old", title="Old")
    # Backdate: rewrite updated_at older than 1 hour
    p = tmp_sessions_dir / "old.json"
    import json
    data = json.loads(p.read_text())
    data["updated_at"] = time.time() - 7200
    p.write_text(json.dumps(data))
    assert SS.resume_latest(max_age_s=3600) is None
    assert SS.resume_latest(max_age_s=None) is not None


def test_cleanup_old_sessions(tmp_sessions_dir):
    _save("recent")
    _save("ancient")
    # Backdate 'ancient' file mtime by 60 days
    p = tmp_sessions_dir / "ancient.json"
    old_t = time.time() - 60 * 86_400
    os.utime(p, (old_t, old_t))
    n = SS.cleanup_old_sessions(max_age_days=30)
    assert n == 1
    assert (tmp_sessions_dir / "recent.json").exists()
    assert not (tmp_sessions_dir / "ancient.json").exists()


def test_find_sessions_by_title(tmp_sessions_dir):
    _save("a", title="Refactor logging system")
    _save("b", title="Fix audit log bug")
    _save("c", title="Unrelated")
    matches = SS.find_sessions("audit")
    ids = [s["session_id"] for s in matches]
    assert "b" in ids
    assert "c" not in ids


def test_find_sessions_by_first_message(tmp_sessions_dir):
    _save("x", title="generic", msg="please refactor the auth module")
    _save("y", title="generic", msg="run the build")
    matches = SS.find_sessions("auth")
    ids = [s["session_id"] for s in matches]
    assert "x" in ids
    assert "y" not in ids


def test_find_sessions_empty_query_returns_all(tmp_sessions_dir):
    _save("a")
    _save("b")
    matches = SS.find_sessions("")
    assert len(matches) == 2
