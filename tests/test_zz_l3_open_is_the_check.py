"""L3 (operator build): what the index opens is decided on the open file.

A check on the path followed by a read of the path can be raced: the
sessions directory is writable by the agent. The read itself must refuse.
"""

import os
import time

from delfin.agent import session_index as SI


def test_a_symlink_is_not_followed_even_to_a_regular_file(tmp_path):
    target = tmp_path / "other_account.json"
    target.write_text('{"title": "SECRET-TITLE"}')
    link = tmp_path / "s.json"
    link.symlink_to(target)
    assert SI._safe_read(link) is None
    assert SI._is_regular_file(link) is False


def test_a_fifo_returns_at_once(tmp_path):
    fifo = tmp_path / "s.jsonl"
    os.mkfifo(fifo)
    started = time.monotonic()
    assert SI._safe_read(fifo) is None
    assert time.monotonic() - started < 2.0


def test_a_directory_and_a_missing_file_yield_nothing(tmp_path):
    assert SI._safe_read(tmp_path) is None
    assert SI._safe_read(tmp_path / "absent.json") is None


def test_the_cap_is_on_what_is_read_not_only_on_what_stat_said(tmp_path):
    f = tmp_path / "s.json"
    f.write_bytes(b"x" * 100)
    assert SI._safe_read(f, max_bytes=99) is None
    data, st = SI._safe_read(f, max_bytes=100)
    assert len(data) == 100 and st.st_size == 100


def test_a_session_id_is_a_stem_never_a_path():
    for bad in ("../x", "a/b", "", ".hidden", "a" * 200, "x\ny", None):
        assert SI._valid_session_id(bad) is False, bad
    for good in ("sess-1", "2026-10-11_abc.def", "A1"):
        assert SI._valid_session_id(good) is True, good
    assert SI.index_session("../../etc/passwd") is False


def test_checkpoints_and_links_are_not_sessions(tmp_path, monkeypatch):
    from delfin.agent import session_store
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path)
    (tmp_path / "real.json").write_text("{}")
    (tmp_path / "real.turn.json").write_text("{}")
    (tmp_path / "link.json").symlink_to(tmp_path / "real.json")
    os.mkfifo(tmp_path / "pipe.json")
    (tmp_path / "sub.json").mkdir()
    assert SI._session_ids_on_disk() == ["real"]
