"""The resume table says whether a session still lives.

A saved record carries model, provider, host and a heartbeat since the
2026-09-19 fix in ``session_store`` — but ``cli_resume`` kept its own row
reader and dropped all of it: the terminal's ``--sessions`` table showed
neither when a session was last seen nor whether its process was still
running. With several parallel sessions on several nodes that is exactly
the question the reader asks.
"""

from __future__ import annotations

import json
import os
import time

from delfin.agent import cli_resume as CR


def _write(d, sid, **fields):
    data = {"session_id": sid, "created_at": time.time() - 3600,
            "updated_at": time.time() - 60, "model": "kit",
            "chat_messages": [{"role": "user", "content": "job"}]}
    data.update(fields)
    (d / f"{sid}.json").write_text(json.dumps(data))


def test_rows_carry_the_heartbeat_and_liveness(tmp_path):
    _write(tmp_path, "live-1", pid=os.getpid(),
           proc_start=CR.proc_identity.process_start(os.getpid()),
           host=CR.proc_identity.socket.gethostname(),
           heartbeat_at=time.time() - 5)
    rows = CR.list_sessions(sessions_dir=tmp_path)
    assert rows[0]["heartbeat_at"] > 0
    assert rows[0]["alive"] is True


def test_a_dead_pid_is_marked_dead(tmp_path):
    _write(tmp_path, "dead-1", pid=999_999_999, proc_start="12345",
           host="", heartbeat_at=time.time() - 5)
    rows = CR.list_sessions(sessions_dir=tmp_path)
    assert rows[0]["alive"] is False


def test_an_old_record_without_the_fields_says_unknown(tmp_path):
    _write(tmp_path, "old-1")   # no pid, no heartbeat
    rows = CR.list_sessions(sessions_dir=tmp_path)
    assert rows[0]["alive"] is None
    assert rows[0]["heartbeat_at"] == 0


def test_the_table_shows_last_seen_and_life(tmp_path):
    _write(tmp_path, "t-live", pid=os.getpid(),
           proc_start=CR.proc_identity.process_start(os.getpid()),
           host="", heartbeat_at=time.time() - 3)
    _write(tmp_path, "t-old")   # predates the fields
    text = CR.render_sessions(CR.list_sessions(sessions_dir=tmp_path))
    assert "last seen" in text          # the column is named
    assert "LIVE" in text               # the live one is marked
    assert "unknown" in text            # the old one degrades, not crashes


def test_a_stranger_on_the_pid_is_not_alive(tmp_path):
    _write(tmp_path, "t-stranger", pid=os.getpid(), proc_start="1",
           host="", heartbeat_at=time.time())
    rows = CR.list_sessions(sessions_dir=tmp_path)
    assert rows[0]["alive"] is False
