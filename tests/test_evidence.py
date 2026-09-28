"""Controls for delfin/agent/evidence.py (package 8: verify_evidence).

Each test names the behaviour the wave contract demands: an evidence
entry is CHECKED, never believed. Red on the previous commit (module
did not exist).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from delfin.agent import evidence as ev  # noqa: E402


# ---------------------------------------------------------------------------
# test-kind evidence
# ---------------------------------------------------------------------------

def test_test_evidence_green_run_ok(tmp_path):
    f = tmp_path / "tests" / "test_x.py"
    f.parent.mkdir()
    f.write_text("def test_y():\n    pass\n")
    runs = [{"command": "tests/test_x.py", "exit_code": 0,
             "status": "ok", "passed": 1, "failed": 0}]
    ok, detail = ev.verify_evidence(
        {"kind": "test", "ref": "tests/test_x.py::test_y"},
        workspace=tmp_path, runs=runs)
    assert ok is True, detail


def test_test_evidence_missing_node_fails(tmp_path):
    ok, detail = ev.verify_evidence(
        {"kind": "test", "ref": "tests/nope.py::test_y"}, workspace=tmp_path)
    assert ok is False and "not found" in detail


def test_test_evidence_red_run_fails(tmp_path):
    f = tmp_path / "tests" / "test_x.py"
    f.parent.mkdir()
    f.write_text("def test_y():\n    pass\n")
    runs = [{"command": "tests/test_x.py", "exit_code": 1,
             "status": "failed", "passed": 0, "failed": 1}]
    ok, detail = ev.verify_evidence(
        {"kind": "test", "ref": "tests/test_x.py::test_y"},
        workspace=tmp_path, runs=runs)
    assert ok is False and "not green" in detail


def test_test_evidence_no_run_record_fails(tmp_path):
    f = tmp_path / "tests" / "test_x.py"
    f.parent.mkdir()
    f.write_text("def test_y():\n    pass\n")
    ok, detail = ev.verify_evidence(
        {"kind": "test", "ref": "tests/test_x.py::test_y"},
        workspace=tmp_path, runs=[])
    assert ok is False and "no green run" in detail


# ---------------------------------------------------------------------------
# job-kind evidence
# ---------------------------------------------------------------------------

class _Job:
    """Stand-in shaped like the job objects DELFIN's list_jobs returns."""
    def __init__(self, job_id, state="RUNNING", status="ok"):
        self.job_id = job_id
        self.state = state
        self.status = status


def test_job_evidence_running_ok():
    ok, detail = ev.verify_evidence(
        {"kind": "job", "ref": "12345"},
        list_jobs=lambda: [_Job("12345", state="RUNNING")])
    assert ok is True, detail


def test_job_evidence_failed_state_fails():
    """A job id merely being listed is not evidence: a FAILED job must
    be rejected, not believed. Red before the freshness fix."""
    ok, detail = ev.verify_evidence(
        {"kind": "job", "ref": "12345"},
        list_jobs=lambda: [_Job("12345", state="FAILED")])
    assert ok is False and "state" in detail


def test_job_evidence_completed_state_fails():
    """A COMPLETED job is history, not a live run: must be rejected."""
    ok, detail = ev.verify_evidence(
        {"kind": "job", "ref": "12345"},
        list_jobs=lambda: [_Job("12345", state="COMPLETED")])
    assert ok is False and "state" in detail


def test_job_evidence_cancelled_state_fails():
    ok, detail = ev.verify_evidence(
        {"kind": "job", "ref": "12345"},
        list_jobs=lambda: [_Job("12345", state="CANCELLED")])
    assert ok is False and "state" in detail


def test_job_evidence_status_failed_fails():
    """If only a `status` field is populated, a failed one rejects too."""
    ok, detail = ev.verify_evidence(
        {"kind": "job", "ref": "12345"},
        list_jobs=lambda: [_Job("12345", state="", status="failed")])
    assert ok is False and "state" in detail


def test_job_evidence_unknown_state_rejects_fail_closed():
    """An unknown/empty state is unconfirmed, not accepted."""
    ok, detail = ev.verify_evidence(
        {"kind": "job", "ref": "12345"},
        list_jobs=lambda: [_Job("12345", state="WEIRD")])
    assert ok is False


# ---------------------------------------------------------------------------
# non-git workspaces: a green run cannot be judged, only marked
# ---------------------------------------------------------------------------

def test_test_evidence_non_git_marked_unconfirmed(tmp_path, monkeypatch):
    """Outside git nothing is stamped, so nothing can be judged -- the
    green run must come back marked unconfirmed, not believed."""
    f = tmp_path / "tests" / "test_x.py"
    f.parent.mkdir()
    f.write_text("def test_y():\n    pass\n")
    runs = [{"command": "tests/test_x.py", "exit_code": 0,
             "status": "ok", "passed": 1, "failed": 0}]
    # make the workspace look non-git: fingerprint returns no commit
    import delfin.agent.evidence_freshness as ef
    monkeypatch.setattr(ef, "fingerprint", lambda ws: {"commit": ""})
    ok, detail = ev.verify_evidence(
        {"kind": "test", "ref": "tests/test_x.py::test_y"},
        workspace=tmp_path, runs=runs)
    assert ok is True and "unconfirmed" in detail


def test_unknown_kind_fails():
    ok, detail = ev.verify_evidence({"kind": "gutfeeling", "ref": "x"})
    assert ok is False and "unknown evidence kind" in detail


def test_evidence_object_with_attrs_supported(tmp_path):
    """The contract's Evidence dataclass must work, not only dicts."""
    class E:  # stand-in shaped like skill_proposals.Evidence
        kind, ref, detail, verified_at = "test", "tests/a.py::t", "", ""
    f = tmp_path / "tests" / "a.py"
    f.parent.mkdir()
    f.write_text("def t():\n    pass\n")
    runs = [{"command": "tests/a.py", "exit_code": 0,
             "status": "ok", "passed": 1, "failed": 0}]
    ok, _ = ev.verify_evidence(E(), workspace=tmp_path, runs=runs)
    assert ok is True


# ---------------------------------------------------------------------------
# integration: staleness reaches the verify_evidence verdict
# ---------------------------------------------------------------------------

def test_stale_green_run_rejected_at_verify_level(tmp_path, monkeypatch):
    """A recorded green run whose tree changed since the run must be
    REJECTED by verify_evidence itself, not only by is_stale. Uses a
    fake git-shaped workspace: fingerprint says a commit exists, the
    recorded run's stamp names a different (older) commit."""
    f = tmp_path / "tests" / "test_x.py"
    f.parent.mkdir()
    f.write_text("def test_y():\n    pass\n")
    runs = [{"command": "tests/test_x.py", "exit_code": 0,
             "status": "ok", "passed": 1, "failed": 0,
             "fingerprint": {"commit": "oldsha", "dirty": {}}}]
    import delfin.agent.evidence_freshness as ef
    monkeypatch.setattr(ef, "fingerprint", lambda ws: {
        "commit": "newsha", "dirty": {}})
    ok, detail = ev.verify_evidence(
        {"kind": "test", "ref": "tests/test_x.py::test_y"},
        workspace=tmp_path, runs=runs)
    assert ok is False and "stale" in detail
