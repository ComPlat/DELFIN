"""Package SW1 phase 1: the swarm status JSON contract (SW2 pins on it).

Control: tests/test_sw1_status_schema.py, 4 cases, all red on the
previous commit (ModuleNotFoundError: delfin.agent.swarm.status).
"""

from __future__ import annotations

from delfin.agent.swarm import status


def test_schema_version_is_pinned_and_report_validates():
    assert status.STATUS_SCHEMA_VERSION == 1
    assert status.validate_report(status.STATUS_EXAMPLE) == []


def test_broken_reports_are_caught():
    bad = dict(status.STATUS_EXAMPLE)
    bad.pop("goal")
    assert any("goal" in p for p in status.validate_report(bad))
    bad2 = {"schema_version": 2, "goal": "x", "packages": []}
    assert status.validate_report(bad2)


def test_unknown_phase_is_rejected():
    bad = dict(status.STATUS_EXAMPLE)
    bad["packages"] = [dict(bad["packages"][0])]
    bad["packages"][0]["phase"] = "banana"
    assert any("phase" in p for p in status.validate_report(bad))


def test_events_reach_subscribers_and_are_cleaned_up():
    seen = []
    off = status.subscribe(seen.append)
    status.emit_event({"kind": "phase_change", "package": "core",
                       "phase": "in_review"})
    off()
    status.emit_event({"kind": "phase_change", "package": "core",
                       "phase": "qs_ok"})
    assert len(seen) == 1
    assert seen[0]["kind"] == "phase_change"


def test_a_subscriber_raising_does_not_break_others():
    got = []

    def boom(_):
        raise RuntimeError("bad subscriber")

    off1 = status.subscribe(boom)
    off2 = status.subscribe(got.append)
    status.emit_event({"kind": "member_stalled"})
    off1()
    off2()
    assert got == [{"kind": "member_stalled"}]
