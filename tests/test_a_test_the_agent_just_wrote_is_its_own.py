"""The test-tamper gate does not fire on a test file the session itself created.

Driven 2026-09-16: an agent wrote a NEW tests/test_delfin_doctor.py, ran it
(red: its fixture kept the real PATH), fixed the fixture, and every fix came
back with "TEST-TAMPER GATE ... your final answer MUST state why the test
was wrong". Existing tests, and a test file that was overwritten rather
than created, stay guarded.
"""
from delfin.agent import api_client as A


def _red(evidence, red):
    A._observe_test_evidence(
        evidence, red, "run_tests", {"target": "tests/test_delfin_doctor.py"},
        '{"status": "failed", "summary": {"passed": 1, "failed": 1}, '
        '"failures": [{"node_id": "tests/test_delfin_doctor.py::test_x"}]}')
    assert red


def test_a_created_test_file_is_the_agents_own(monkeypatch):
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    evidence, red, created = [], set(), set()
    A._observe_test_evidence(evidence, red, "write_file",
                             {"path": "tests/test_delfin_doctor.py", "content": "x"},
                             "File created: tests/test_delfin_doctor.py\n\n+x", created=created)
    _red(evidence, red)
    note = A._observe_test_evidence(evidence, red, "edit_file",
                                    {"path": "tests/test_delfin_doctor.py"},
                                    "File edited", created=created)
    assert note == ""


def test_an_existing_test_is_still_guarded(monkeypatch):
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    evidence, red, created = [], set(), set()
    A._observe_test_evidence(evidence, red, "write_file",
                             {"path": "tests/test_delfin_doctor.py", "content": "x"},
                             "File overwritten: tests/test_delfin_doctor.py", created=created)
    _red(evidence, red)
    note = A._observe_test_evidence(evidence, red, "edit_file",
                                    {"path": "tests/test_delfin_doctor.py"},
                                    "File edited", created=created)
    assert "TEST-TAMPER GATE" in note


def test_a_diff_from_dev_null_counts_as_created(monkeypatch):
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    evidence, red, created = [], set(), set()
    diff = "--- /dev/null\n+++ b/tests/test_new.py\n@@ -0,0 +1 @@\n+def test(): pass\n"
    A._observe_test_evidence(evidence, red, "apply_patch", {"diff": diff}, "Patch applied", created=created)
    assert "tests/test_new.py" in created


def test_without_the_created_set_the_old_behaviour_holds(monkeypatch):
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    evidence, red = [], set()
    _red(evidence, red)
    note = A._observe_test_evidence(evidence, red, "edit_file",
                                    {"path": "tests/test_delfin_doctor.py"}, "File edited")
    assert "TEST-TAMPER GATE" in note
