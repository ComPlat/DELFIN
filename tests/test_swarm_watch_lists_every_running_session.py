"""Swarm Watch without a roster lists every running session, on every node.

A wave or swarm started by a script runs in a SLURM job on a compute node;
the dashboard knows no roster for it and the tab said "No swarm to watch"
while eleven sessions ran. Without a roster the tab now lists every open
session from the presence records in the shared HOME, with its host. A
session whose heartbeat stopped (its job ended) is left out. The read skips
the per-worktree git walk (seconds per call) and runs at most every 30 s.
"""
from __future__ import annotations

from delfin.dashboard import tab_swarm as T

NOW = 1_000_000.0


def _rec(name, host, **kw):
    base = {"name": name, "session_id": name + "-id", "host": host,
            "branch": f"agent/{name}", "open_approvals": [], "idle_s": 30,
            "last_turn": {"model": "kit.glm"}, "pid_alive": True,
            "paused": False, "inbox": 0}
    base.update(kw)
    return base


def test_sessions_of_other_nodes_are_listed_with_their_host():
    recs = [_rec("nacht-s1", "uc3n001.localdomain"),
            _rec("nacht-s5", "uc3n034.localdomain")]
    rows = T.running_session_rows(recs, now=NOW, settings={},
                                  beats={"nacht-s1": NOW - 20, "nacht-s5": NOW - 40})
    assert [(r["name"], r["host"]) for r in rows] == [
        ("nacht-s1", "uc3n001"), ("nacht-s5", "uc3n034")]
    assert rows[0]["model"] == "kit.glm"


def test_a_session_whose_heartbeat_stopped_is_left_out():
    recs = [_rec("alive", "n1"), _rec("ended", "n2")]
    rows = T.running_session_rows(
        recs, now=NOW, settings={},
        beats={"alive": NOW - 30, "ended": NOW - T.HEARTBEAT_STALE_S - 1})
    assert [r["name"] for r in rows] == ["alive"]


def test_a_pause_flag_without_a_session_is_no_row():
    recs = [{"name": "nacht-s30", "session_id": "", "host": None}]
    assert T.running_session_rows(recs, now=NOW, settings={}, beats={}) == []


def test_rows_are_escaped_and_flag_open_approvals():
    rows = [{"name": "<b>x</b>", "host": "n1", "branch": "b", "model": "m",
             "state": "working", "idle": "1m", "approvals": 2, "inbox": 0}]
    html = "".join(T.running_rows_html(rows))
    assert "&lt;b&gt;x&lt;/b&gt;" in html and "2 open" in html
    assert "<button" not in html and "<input" not in html


def test_the_tab_reads_without_the_git_walk(monkeypatch):
    seen = {}
    from delfin.agent import session_health

    def collect(settings=None, *, with_git=True):
        seen["with_git"] = with_git
        return []
    monkeypatch.setattr(session_health, "collect", collect)
    T.running_session_rows(now=NOW, settings={}, beats={})
    assert seen["with_git"] is False
