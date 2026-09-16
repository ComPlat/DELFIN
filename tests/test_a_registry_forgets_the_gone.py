"""Presence records of sessions long gone are removed, and one take of the
session inbox is bounded.

Review 2026-09-16: a crashed kernel never withdraws its presence record, so
the directory grew by one file per crash and every refresh read them all;
and take() handed every waiting message to one prompt with no cap.
"""
import json
import os
import time

from delfin.agent import session_messages as SM
from delfin.agent import session_presence as P


def _record(tmp_path, monkeypatch, key, *, age_s, host="elsewhere", pid=1):
    monkeypatch.setattr(P, "_DIR", tmp_path)
    rec = {"key": key, "host": host, "pid": pid, "updated_at": time.time() - age_s,
           "workspace": "/w", "title": key}
    (tmp_path / f"{key}.json").write_text(json.dumps(rec))
    return tmp_path / f"{key}.json"


def test_a_record_long_stale_is_removed_and_a_fresh_one_kept(tmp_path, monkeypatch):
    old = _record(tmp_path, monkeypatch, "old", age_s=P._REAP_AFTER_S + 60)
    fresh = _record(tmp_path, monkeypatch, "fresh", age_s=5)
    stale = _record(tmp_path, monkeypatch, "stale", age_s=P._STALE_S + 60)
    keys = [r["key"] for r in P.open_sessions()]
    assert keys == ["fresh"]
    assert not old.exists(), "long gone: reaped"
    assert stale.exists(), "merely stale: skipped, its session may come back"
    assert fresh.exists()


def test_a_dead_process_on_this_host_is_reaped_at_once(tmp_path, monkeypatch):
    import socket
    dead = _record(tmp_path, monkeypatch, "dead", age_s=5, host=socket.gethostname(), pid=2 ** 22 - 7)
    P.open_sessions()
    assert not dead.exists()


def test_one_take_is_bounded_and_says_what_it_dropped(tmp_path, monkeypatch):
    monkeypatch.setattr(SM, "_inbox", lambda key: tmp_path / f"{key}.jsonl")
    for i in range(SM._MAX_TAKE + 5):
        SM.send("b", f"m{i}", from_key="a")
    got = SM.take("b")
    assert len(got) == SM._MAX_TAKE + 1
    assert "5 earlier message(s) were not delivered" in got[0]["text"]
    assert got[-1]["text"] == f"m{SM._MAX_TAKE + 4}"
    assert SM.take("b") == []
