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
