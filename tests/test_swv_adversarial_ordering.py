"""Adversarial tests for SWV (reviewer nacht-s17, branch agent/s17-swvr1).

These pin the ordering and state contracts of delfin/agent/swarm/observe.py
(member_feed and _state_of). Each test is a control that was RED before the
builder's fixes and GREEN at the merged head 000c4f022:

* F1 -- feed step ordering: member_feed must return steps newest-LAST
  (the audit log already delivers steps newest-last via steps_since, and the
  blanket entries.reverse() used to flip them to newest-FIRST). Commits (git
  log is newest-first) are inverted separately so the whole feed is newest
  last.
* F2 -- the brief's 6th member state "waiting for a message": a member whose
  session-health shape is `waiting_on_person` (no approval, not paused, no
  resume planned) must project as "waiting for a message", not "working".
* commit order -- commit entries stay newest-LAST relative to each other.

These adversarial tests are independent of the builder's own suite; they call
the module's public API (member_feed, _state_of) with the same hostile/edge
shapes that previously produced the wrong result.
"""
from __future__ import annotations

from delfin.agent.swarm import observe as _observe


def _oldest_first_step(offset=0, *, session="", **kw):
    # Mirrors audit_log.steps_since: rows oldest-first (newest last), so a
    # member_feed call must NOT reverse them back to newest-first. Accepts
    # member_feed's call shape: (since, session=..., log_path=...).
    rows = [
        {"ts": "2026-06-01T10:00:00", "tool": "read_file", "command": "a.py",
         "decision": ""},
        {"ts": "2026-06-01T10:00:05", "tool": "write_file", "command": "b.py",
         "decision": ""},
    ]
    return rows, 999


def _commits_newest_first(member_name: str):
    # Mirrors _default_commits via git log: NEWEST-first subjects. member_feed
    # must invert ONLY these (not the steps) to satisfy newest-last.
    return {member_name: ["most-recent-commit", "older-commit"]}


def _health_waiting_on_person():
    # The shape session_health.collect() builds for a member blocked on a
    # person/message (verdict waiting_on_person): no approval, not paused, no
    # resume planned, turn not flagged stopped.
    return {
        "open_approvals": [],
        "paused": False,
        "no_token_s": 30,
        "resume_planned": False,
        "idle_s": 200,
        "last_turn": {"stopped": False},
    }


def test_feed_tool_steps_newest_last():
    """F1: the feed's tool-call steps must be newest-LAST, never reversed."""
    feed = _observe.member_feed("s1", steps_fn=_oldest_first_step,
                                commits_fn=lambda s: {s: []})
    kinds = [e["text"] for e in feed if e["kind"] == "tool"]
    # Oldest (10:00:00 a.py) first, newest (10:00:05 b.py) last.
    assert kinds == ["a.py", "b.py"], kinds


def test_feed_commits_newest_last():
    """Commit entries stay newest-LAST relative to each other."""
    feed = _observe.member_feed("s1", steps_fn=_oldest_first_step,
                                commits_fn=_commits_newest_first)
    commits = [e["text"] for e in feed if e["kind"] == "commit"]
    # git log gave newest-first; the feed must present newest last.
    assert commits == ["older-commit", "most-recent-commit"], commits


def test_state_is_waiting_for_a_message_when_blocked_on_a_person():
    """F2 (QS route A): a member blocked on a person is 'waiting for a message'."""
    state = _observe._state_of(_health_waiting_on_person())
    assert state == "waiting for a message", state
