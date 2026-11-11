"""SWV phase 1: the read-only swarm overview observe contract (package SWV).

Two pure functions in delfin/agent/swarm/observe.py, called by the dashboard
(tab_swarm.py), the TUI (tui/src/swarm) and an SW3 lead — one data API:

  swarm_overview(source[, records=, steps_fn=]) -> one row per member
  member_feed(session, since[, records=, steps_fn=]) -> the member's steps

Data writers are EXISTING stores, never a second trace: audit_log.steps_since
(tool calls + refusals) and the session_health-shaped record (session_health
collect()'s shape: open_approvals / paused / idle_s / no_token_s /
last_turn.stopped / api_errors_in_a_row). observe.py may inject those writers
(records=, steps_fn=) so tests feed real-record fixtures; the production
wrappers call session_health.collect() and audit_log.steps_since.

Safety (task target 4, QS ruling): pure reads only — observing never prompts,
pauses, answers an approval or writes anything; a feed is shown only to the
session's own user (the ownership check lives in the service layer; this
module renders only the data it is given).

Red control: observe.py does not exist on the base commit, so every case here
fails red on the previous commit.
"""
from __future__ import annotations

import json
import time

import pytest

from delfin.agent import audit_log as _al
from delfin.agent.swarm import observe as _observe  # missing -> red


# ---------------------------------------------------------------------------
# real-record helpers
# ---------------------------------------------------------------------------

def _health_record(session: str, **over) -> dict:
    """One record shaped exactly as session_health.collect() produces it."""
    rec = {
        "name": session,
        "session_id": session,
        "title": "",
        "paused": False,
        "pid_alive": True,
        "open_approvals": [],
        "approval_opened_at": None,
        "idle_s": 5.0,
        "no_token_s": 0.0,
        "api_errors_in_a_row": 0,
        "open_work": False,
        "resume_planned": True,
        "started_at": time.time() - 3600,
        "last_turn": {"stopped": False},
    }
    rec.update(over)
    return rec


def _step(session: str, tool: str, command: str, *,
          decision: str = "", reason: str = "") -> dict:
    rec = {"session_id": session, "tool": tool, "command": command}
    if decision:
        rec["decision"] = decision
    if reason:
        rec["reason"] = reason
    return rec


@pytest.fixture
def audit_log(tmp_path):
    """Write real records through audit_log.append into a tmp log."""
    path = tmp_path / "audit.jsonl"
    def write(*records):
        for r in records:
            _al.append(r, log_path=path)
        return path
    return write


# ---------------------------------------------------------------------------
# swarm_overview
# ---------------------------------------------------------------------------

def _spec(packages):
    return {"goal": "x", "packages": packages}


def test_overview_one_row_per_member_with_all_fields(audit_log):
    source = _spec([{"name": "swv", "roles": {"builder": True, "reviewer": True}}])
    records = [_health_record("swv-builder"),
               _health_record("swv-reviewer")]
    rows = _observe.swarm_overview(source, records=records, steps_fn=lambda *a, **k: ([], 0))
    assert isinstance(rows, list)
    by_name = {r["role"]: r for r in rows}
    assert set(by_name) == {"builder", "reviewer"}
    for row in rows:
        for field in ("package", "role", "state", "model", "commits_on_branch",
                      "last_tool_call", "last_commit_subject",
                      "last_activity_age_s", "pending_approvals_count"):
            assert field in row, f"missing {field}"
    assert by_name["builder"]["package"] == "swv"
    assert by_name["builder"]["state"] == row_state(records[0])


def row_state(record: dict) -> str:
    """The SWV state that maps from a session_health record (mirrors the
    module's own mapping; the overview uses it)."""
    if record.get("open_approvals"):
        return "NEEDS APPROVAL"
    if record.get("paused"):
        return "paused"
    # no_token_s + last_turn.stopped mirror session_health.classify's axes;
    # the overview keeps the same priority: approval > paused > ended > model
    turn = record.get("last_turn") or {}
    if turn.get("stopped") is True and not record.get("resume_planned"):
        return "ended"
    if (record.get("no_token_s") or 0) >= 10 * 60:
        return "waiting for model"
    return "working"


def test_overview_state_map_all_six_cases(audit_log):
    source = _spec([{"name": "p", "roles": {
        "a": True, "b": True, "c": True, "d": True, "e": True, "f": True}}])
    records = [
        _health_record("p-a"),                                        # working
        _health_record("p-b", no_token_s=11 * 60),                    # waiting for model
        _health_record("p-c", idle_s=200, resume_planned=False,       # ended (stopped turn, nothing resumes)
                       last_turn={"stopped": True}),
        _health_record("p-d", open_approvals=["r1"]),                 # NEEDS APPROVAL
        _health_record("p-e", paused=True),                           # paused
        _health_record("p-f", idle_s=200, resume_planned=False,       # waiting for a message
                       last_turn={"stopped": False}),
    ]
    rows = _observe.swarm_overview(source, records=records,
                                   steps_fn=lambda *a, **k: ([], 0))
    states = {r["role"]: r["state"] for r in rows}
    assert states["a"] == "working"
    assert states["b"] == "waiting for model"
    assert states["c"] == "ended"
    assert states["d"] == "NEEDS APPROVAL"
    assert states["e"] == "paused"
    assert states["f"] == "waiting for a message"


def test_overview_last_tool_call_and_activity_age_from_audit_log(audit_log):
    source = _spec([{"name": "swv", "roles": {"builder": True}}])
    rec = _health_record("swv-builder")
    steps = {"swv-builder": [
        {"ts": "2026-06-01T10:00:05", "tool": "run_tests",
         "command": "tests/test_swv_observe.py", "decision": ""},
        {"ts": "2026-06-01T10:00:10", "tool": "write_file",
         "command": "delfin/x.py", "decision": "denied",
         "reason": "outside scope"},
    ]}
    def steps_fn(offset=0, *, session="", denied_only=False, **kw):
        all_rows = sorted(steps.get(session, []), key=lambda s: s["ts"])
        return all_rows[offset:], len(all_rows)

    rows = _observe.swarm_overview(source, records=[rec], steps_fn=steps_fn)
    row = rows[0]
    assert row["last_tool_call"]["tool"] == "write_file"
    assert "delfin/x.py" in row["last_tool_call"]["arg"]
    assert row["pending_approvals_count"] == 0
    assert isinstance(row["last_activity_age_s"], (int, float))


def test_overview_pending_approvals_count(audit_log):
    source = _spec([{"name": "p", "roles": {"builder": True}}])
    records = [_health_record("p-builder", open_approvals=["r1", "r2"])]
    rows = _observe.swarm_overview(source, records=records,
                                   steps_fn=lambda *a, **k: ([], 0))
    assert rows[0]["pending_approvals_count"] == 2


def test_overview_roster_source_plain_file(audit_log, tmp_path):
    """An operator wave has no spec: a plain roster file (name -> package/role)
    passed explicitly maps to the same rows."""
    roster = tmp_path / "roster.json"
    roster.write_text(json.dumps({
        "ops-1": {"package": "ops", "role": "builder"},
        "ops-2": {"package": "ops", "role": "reviewer"},
    }), encoding="utf-8")
    records = [_health_record("ops-1"), _health_record("ops-2")]
    rows = _observe.swarm_overview(str(roster), records=records,
                                   steps_fn=lambda *a, **k: ([], 0))
    assert len(rows) == 2
    by_name = {r["role"]: r for r in rows}
    assert by_name["builder"]["package"] == "ops"
    assert by_name["reviewer"]["package"] == "ops"


def test_overview_roster_source_plain_dict(audit_log):
    """The same roster handed over as an object (the RPC's ``source``)."""
    roster = {"ops-1": {"package": "ops", "role": "builder"}}
    rows = _observe.swarm_overview(roster, records=[_health_record("ops-1")],
                                   steps_fn=lambda *a, **k: ([], 0))
    assert [(r["name"], r["package"], r["role"]) for r in rows] == [
        ("ops-1", "ops", "builder")]


def test_overview_spec_file_source(audit_log, tmp_path):
    """A swarm.yaml path gives its members, not an empty table."""
    spec = tmp_path / "swarm.yaml"
    spec.write_text(json.dumps({"packages": [
        {"name": "alpha", "roles": {"builder": True, "reviewer": True}}]}),
        encoding="utf-8")
    rows = _observe.swarm_overview(str(spec), records=[],
                                   steps_fn=lambda *a, **k: ([], 0))
    assert sorted(r["name"] for r in rows) == ["alpha-builder", "alpha-reviewer"]


def test_overview_missing_member_never_lies(audit_log):
    """A member with no session record still has a row, its state 'no
    session': absence proves neither work nor a pause -- never invent a
    liveness."""
    source = _spec([{"name": "swv", "roles": {"builder": True}}])
    rows = _observe.swarm_overview(source, records=[],  # no record at all
                                   steps_fn=lambda *a, **k: ([], 0))
    assert len(rows) == 1
    assert rows[0]["state"] == "no session"


# ---------------------------------------------------------------------------
# member_feed
# ---------------------------------------------------------------------------

def test_feed_steps_newest_last_with_result_status(audit_log):
    path = audit_log(
        _step("swv-builder", "bash", "git status", decision=""),
        _step("swv-builder", "write_file", "delfin/x.py", decision="denied",
              reason="outside scope"),
        _step("swv-builder", "run_tests", "tests/t.py", decision=""),
    )
    feed = _observe.member_feed("swv-builder", since=0,
                                records=[_health_record("swv-builder")],
                                steps_fn=_al.steps_since, log_path=path)
    # newest last -- pin the REAL order, not a vacuous kind check. The
    # steps are written into the log oldest-first (bash, write_file,
    # run_tests) and steps_since returns them the same way, so the newest
    # step (run_tests) must be the LAST feed entry. The old check
    # `kinds[-1] in ("tool", "refusal")` was always true and never caught
    # the blanket entries.reverse() that flipped the steps to newest-first
    # (F1), mixing them with the commits.
    step_texts = [e["text"] for e in feed if e["kind"] in ("tool", "refusal")]
    assert step_texts[0] == "git status"    # oldest step first
    assert step_texts[-1] == "tests/t.py"   # newest step last
    assert all(e["kind"] in ("tool", "refusal", "commit") for e in feed)
    # the denial is present and marked as a refusal
    denied = [e for e in feed if e["kind"] == "refusal"]
    assert denied and "outside scope" in denied[0]["text"]


def test_overview_without_records_reads_the_real_sessions(monkeypatch):
    """The dashboard and the TUI call swarm_overview(source) with no records:
    it must read the live sessions, not treat every member as missing (the
    whole running wave showed 'no session')."""
    monkeypatch.setattr(_observe, "_default_records",
                        lambda: [_health_record("ops-1")])
    rows = _observe.swarm_overview({"ops-1": {"package": "ops", "role": "builder"}},
                                   steps_fn=lambda *a, **k: ([], 0),
                                   commits_fn=lambda name: {})
    assert rows[0]["state"] == "working"
