"""Tests for delfin.dashboard.job_events (D5 — proactive job notifications)."""
from __future__ import annotations

from delfin.dashboard import job_events


def _job(job_id, status, name="", job_dir=""):
    """Build a plain-dict JobInfo replacement for tests."""
    return {
        "job_id": job_id, "status": status, "name": name, "job_dir": job_dir,
    }


# ---------------------------------------------------------------------------
# diff_job_states
# ---------------------------------------------------------------------------

def test_diff_no_change_returns_empty():
    snap = {"1": _job("1", "RUNNING")}
    assert job_events.diff_job_states(snap, snap) == []


def test_diff_status_change_emits_event():
    before = {"1": _job("1", "RUNNING", name="Cas1")}
    after = {"1": _job("1", "COMPLETED", name="Cas1")}
    events = job_events.diff_job_states(before, after)
    assert len(events) == 1
    e = events[0]
    assert e.job_id == "1"
    assert e.from_status == "RUNNING"
    assert e.to_status == "COMPLETED"
    assert e.is_terminal is True
    assert e.transition == "RUNNING→COMPLETED"


def test_diff_new_job_emits_event_with_empty_from_status():
    before = {}
    after = {"42": _job("42", "PENDING", name="new")}
    events = job_events.diff_job_states(before, after)
    assert len(events) == 1
    e = events[0]
    assert e.from_status == ""
    assert e.to_status == "PENDING"
    assert e.transition == "PENDING"


def test_diff_disappeared_job_does_not_emit_event():
    """Jobs missing from the new snapshot are silently dropped."""
    before = {"1": _job("1", "RUNNING")}
    after = {}
    events = job_events.diff_job_states(before, after)
    assert events == []


def test_diff_handles_none_inputs():
    """None inputs are treated as empty maps."""
    assert job_events.diff_job_states(None, None) == []
    assert job_events.diff_job_states({}, None) == []
    # None as ``before`` → all jobs in ``after`` are seen for the first time.
    # That emits events with from_status="" — they're filtered later by
    # is_notable_event.
    out = job_events.diff_job_states(None, {"1": _job("1", "RUNNING")})
    assert len(out) == 1
    assert out[0].from_status == ""


def test_diff_skips_jobs_with_no_status():
    after = {"1": _job("1", "")}
    assert job_events.diff_job_states({}, after) == []


def test_diff_works_with_dataclass_inputs():
    """Real JobInfo dataclasses are supported."""
    from delfin.dashboard.backend_base import JobInfo
    before = {"1": JobInfo(job_id="1", status="RUNNING", name="x")}
    after = {"1": JobInfo(job_id="1", status="FAILED", name="x")}
    events = job_events.diff_job_states(before, after)
    assert len(events) == 1
    assert events[0].to_status == "FAILED"


def test_diff_multiple_simultaneous_transitions():
    before = {
        "1": _job("1", "RUNNING"),
        "2": _job("2", "RUNNING"),
        "3": _job("3", "PENDING"),
    }
    after = {
        "1": _job("1", "COMPLETED"),
        "2": _job("2", "FAILED"),
        "3": _job("3", "RUNNING"),
    }
    events = job_events.diff_job_states(before, after)
    assert len(events) == 3
    transitions = {e.job_id: e.to_status for e in events}
    assert transitions == {
        "1": "COMPLETED", "2": "FAILED", "3": "RUNNING",
    }


# ---------------------------------------------------------------------------
# is_notable_event
# ---------------------------------------------------------------------------

def test_notable_running_to_completed():
    e = job_events.JobEvent(
        job_id="1", name="x", job_dir="", from_status="RUNNING",
        to_status="COMPLETED", is_terminal=True,
    )
    assert job_events.is_notable_event(e) is True


def test_notable_pending_to_failed():
    e = job_events.JobEvent(
        job_id="1", name="x", job_dir="", from_status="PENDING",
        to_status="FAILED", is_terminal=True,
    )
    assert job_events.is_notable_event(e) is True


def test_notable_running_to_timeout():
    e = job_events.JobEvent(
        job_id="1", name="x", job_dir="", from_status="RUNNING",
        to_status="TIMEOUT", is_terminal=True,
    )
    assert job_events.is_notable_event(e) is True


def test_not_notable_pending_to_running():
    """Job started — not chat-worthy."""
    e = job_events.JobEvent(
        job_id="1", name="x", job_dir="", from_status="PENDING",
        to_status="RUNNING", is_terminal=False,
    )
    assert job_events.is_notable_event(e) is False


def test_not_notable_first_time_seen_terminal():
    """We never saw the job running, so seeing it as COMPLETED is not news."""
    e = job_events.JobEvent(
        job_id="1", name="x", job_dir="", from_status="",
        to_status="COMPLETED", is_terminal=True,
    )
    assert job_events.is_notable_event(e) is False


def test_not_notable_first_time_seen_pending():
    e = job_events.JobEvent(
        job_id="1", name="x", job_dir="", from_status="",
        to_status="PENDING", is_terminal=False,
    )
    assert job_events.is_notable_event(e) is False


# ---------------------------------------------------------------------------
# format_job_event_message
# ---------------------------------------------------------------------------

def test_format_completed_includes_followup():
    e = job_events.JobEvent(
        job_id="123", name="Cas_red", job_dir="/calc/Cas_red",
        from_status="RUNNING", to_status="COMPLETED", is_terminal=True,
    )
    msg = job_events.format_job_event_message(e)
    assert "Cas_red" in msg
    assert "123" in msg
    assert "/calc/Cas_red" in msg
    assert "✅" in msg
    # Suggests an analysis or recalc-check
    assert "Energien" in msg or "Recalc" in msg


def test_format_failed_includes_recalc_suggestion():
    e = job_events.JobEvent(
        job_id="9", name="x", job_dir="",
        from_status="RUNNING", to_status="FAILED", is_terminal=True,
    )
    msg = job_events.format_job_event_message(e)
    assert "❌" in msg
    assert "fehlgeschlagen" in msg
    assert "Recalc" in msg


def test_format_timeout_suggests_time_limit_bump():
    e = job_events.JobEvent(
        job_id="9", name="long_job", job_dir="",
        from_status="RUNNING", to_status="TIMEOUT", is_terminal=True,
    )
    msg = job_events.format_job_event_message(e)
    assert "⏱️" in msg
    assert "Zeitlimit" in msg or "time_limit" in msg


def test_format_cancelled_no_action_suggestion():
    """Cancelled = user did it themselves, no proactive suggestion."""
    e = job_events.JobEvent(
        job_id="9", name="x", job_dir="",
        from_status="RUNNING", to_status="CANCELLED", is_terminal=True,
    )
    msg = job_events.format_job_event_message(e)
    assert "🚫" in msg
    assert "abgebrochen" in msg


def test_format_unknown_terminal_uses_generic_message():
    e = job_events.JobEvent(
        job_id="9", name="x", job_dir="",
        from_status="RUNNING", to_status="WEIRD", is_terminal=False,
    )
    msg = job_events.format_job_event_message(e)
    # Falls through to generic "ℹ️ ... transition"
    assert "RUNNING→WEIRD" in msg


def test_format_uses_job_id_when_name_missing():
    e = job_events.JobEvent(
        job_id="42", name="", job_dir="",
        from_status="RUNNING", to_status="COMPLETED", is_terminal=True,
    )
    msg = job_events.format_job_event_message(e)
    assert "42" in msg


# ---------------------------------------------------------------------------
# snapshot_jobs
# ---------------------------------------------------------------------------

def test_snapshot_empty_input():
    assert job_events.snapshot_jobs([]) == {}
    assert job_events.snapshot_jobs(None) == {}


def test_snapshot_keys_jobs_by_id():
    from delfin.dashboard.backend_base import JobInfo
    jobs = [
        JobInfo(job_id="1", name="a", status="RUNNING"),
        JobInfo(job_id="2", name="b", status="PENDING"),
    ]
    snap = job_events.snapshot_jobs(jobs)
    assert set(snap.keys()) == {"1", "2"}
    assert snap["1"].name == "a"


def test_snapshot_skips_jobs_without_id():
    """A JobInfo missing job_id is silently dropped."""
    from delfin.dashboard.backend_base import JobInfo
    jobs = [JobInfo(job_id=None, name="bad")]
    assert job_events.snapshot_jobs(jobs) == {}


# ---------------------------------------------------------------------------
# End-to-end: realistic scenario
# ---------------------------------------------------------------------------

def test_e2e_typical_workflow():
    """Submit → run → finish → diff yields one notable event."""
    # Snapshot 1: just submitted
    snap1 = {
        "100": _job("100", "PENDING", name="prod_run", job_dir="/c/prod"),
    }
    # Snapshot 2: running
    snap2 = {
        "100": _job("100", "RUNNING", name="prod_run", job_dir="/c/prod"),
    }
    # Snapshot 3: failed
    snap3 = {
        "100": _job("100", "FAILED", name="prod_run", job_dir="/c/prod"),
    }
    e1 = job_events.diff_job_states(snap1, snap2)
    e2 = job_events.diff_job_states(snap2, snap3)
    notable_first = [e for e in e1 if job_events.is_notable_event(e)]
    notable_second = [e for e in e2 if job_events.is_notable_event(e)]
    # PENDING→RUNNING is not notable, RUNNING→FAILED is notable
    assert notable_first == []
    assert len(notable_second) == 1
    msg = job_events.format_job_event_message(notable_second[0])
    assert "fehlgeschlagen" in msg
    assert "prod_run" in msg
