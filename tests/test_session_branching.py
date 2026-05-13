"""Tests for session-branching primitives: fork_session +
session_lineage + session_children."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import session_store as ss


@pytest.fixture
def fake_home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    monkeypatch.setattr(
        ss, "_SESSIONS_DIR", tmp_path / ".delfin" / "agent_sessions"
    )
    return tmp_path


def _make(sid: str, title: str = "", **extras):
    return ss.save_session(
        sid,
        mode=extras.pop("mode", "solo"),
        chat_messages=extras.pop("chat_messages",
                                  [{"role": "user", "content": "hi"}]),
        title=title or sid,
        **extras,
    )


def test_fork_records_parent(fake_home):
    _make("root")
    new_sid = ss.fork_session("root")
    assert new_sid is not None
    forked = ss.load_session(new_sid)
    assert forked is not None
    assert forked["forked_from"] == "root"
    # Source is untouched
    source = ss.load_session("root")
    assert source is not None
    assert "forked_from" not in source or not source.get("forked_from")


def test_session_lineage_walks_back_to_root(fake_home):
    """root → a → b → c; lineage of c returns 4 records ending at root."""
    _make("root")
    a = ss.fork_session("root", new_id="a")
    b = ss.fork_session(a, new_id="b")
    c = ss.fork_session(b, new_id="c")
    chain = ss.session_lineage(c)
    sids = [r["session_id"] for r in chain]
    assert sids == ["c", "b", "a", "root"]


def test_session_lineage_handles_missing_parent_gracefully(fake_home):
    """If a referenced parent has been deleted, lineage stops at the
    last reachable session — no crash."""
    _make("orphan_parent_root")
    new_sid = ss.fork_session("orphan_parent_root", new_id="orphan_child")
    ss.delete_session("orphan_parent_root")
    chain = ss.session_lineage(new_sid)
    assert chain[0]["session_id"] == new_sid
    # Walks only as far as available
    assert len(chain) == 1


def test_session_lineage_single_session(fake_home):
    _make("alone")
    chain = ss.session_lineage("alone")
    assert len(chain) == 1
    assert chain[0]["session_id"] == "alone"


def test_session_lineage_breaks_cycles_safely(fake_home, monkeypatch):
    """Defence against pathological data: a forked_from that points at
    itself (or a previously-seen ID) must not loop forever."""
    _make("loop")
    # Hand-craft a self-loop
    sessions_dir = ss._ensure_dir()
    p = sessions_dir / "loop.json"
    data = json.loads(p.read_text(encoding="utf-8"))
    data["forked_from"] = "loop"
    p.write_text(json.dumps(data), encoding="utf-8")
    chain = ss.session_lineage("loop")
    assert len(chain) == 1


def test_session_children_lists_forks(fake_home):
    _make("trunk")
    ss.fork_session("trunk", new_id="branch_a")
    ss.fork_session("trunk", new_id="branch_b")
    ss.fork_session("trunk", new_id="branch_c")
    kids = ss.session_children("trunk")
    sids = {k["session_id"] for k in kids}
    assert sids == {"branch_a", "branch_b", "branch_c"}


def test_session_children_returns_empty_for_leaf(fake_home):
    _make("isolated")
    assert ss.session_children("isolated") == []


def test_fork_session_title_suffix_idempotent(fake_home):
    """Forking a fork should not append the (fork) suffix twice."""
    _make("a", title="My Session")
    b = ss.fork_session("a", new_id="b")
    c = ss.fork_session(b, new_id="c")
    bd = ss.load_session(b)
    cd = ss.load_session(c)
    # Both forks share the same suffix-count
    assert bd["title"].count("(fork)") == 1
    assert cd["title"].count("(fork)") == 1
