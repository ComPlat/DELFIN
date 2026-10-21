"""A multi-part task is done only when every part has evidence.

Welle 11, package B, phase 2. Wave 10 (s14): "fertig" was reported with
1 of 5 parts done. ``check_completion_claim`` judged the subject as one
blob — whichever branch fired first decided, and the other parts were
invisible. When the description enumerates parts (numbered ``1.`` /
``1)``, lettered ``a)`` / ``a.``), each part is now judged as a task of
its own, and one part without evidence holds the whole task back — the
note names the missing part.

Scope decision, stated so it can be challenged: bare dash/bullet lists
are NOT parts. In this repo's task descriptions a dash list is the
shape of context notes ("use pytest -q", "do not push") far more often
than of work parts, and test-vocabulary inside such a note would nag
every task that mentions how to test. Numbered and lettered lists are
the work-part shape (the wave-10 failure used numbered parts).
"""

from __future__ import annotations

from delfin.agent.task_evidence import check_completion_claim


def _check(subject, description="", **kw):
    return check_completion_claim(subject, description, **kw)


DESCRIPTION = (
    "Fix the ordering bug.\n"
    "1. Correct the sort in delfin/report.py\n"
    "2. Add the regression test tests/test_report_order.py\n"
    "3. Write the summary as report.md"
)


def test_one_written_part_of_three_does_not_complete_the_task():
    changes = [{"path": "delfin/report.py", "ts": 100.0, "created": False}]
    check = _check("Fix the ordering bug", DESCRIPTION, changes=changes,
                   window_start=0.0)
    assert check["verdict"] == "unmet", check


def test_the_note_names_the_missing_part():
    changes = [{"path": "delfin/report.py", "ts": 100.0, "created": False}]
    check = _check("Fix the ordering bug", DESCRIPTION, changes=changes,
                   window_start=0.0)
    assert "tests/test_report_order.py" in check["note"], check


def test_the_note_names_the_part_number():
    changes = [{"path": "delfin/report.py", "ts": 100.0, "created": False}]
    check = _check("Fix the ordering bug", DESCRIPTION, changes=changes,
                   window_start=0.0)
    assert "2" in check["note"], check


def test_all_parts_written_verifies():
    changes = [
        {"path": "delfin/report.py", "ts": 100.0, "created": False},
        {"path": "tests/test_report_order.py", "ts": 101.0,
         "created": True},
        {"path": "report.md", "ts": 102.0, "created": True},
    ]
    check = _check("Fix the ordering bug", DESCRIPTION, changes=changes,
                   window_start=0.0)
    assert check["verdict"] == "verified", check


def test_every_part_checkable_and_evidenced_when_subject_carries_a_path():
    """A subject path beside enumerated parts: every piece needs its own
    evidence, the subject is not swallowed by the parts."""
    changes = [
        {"path": "delfin/report.py", "ts": 101.0, "created": False},
        {"path": "tests/test_report_order.py", "ts": 102.0,
         "created": True},
    ]
    subject = "Fix the ordering bug in delfin/report.py"
    check = _check(subject, DESCRIPTION.replace("\n3. Write the summary"
                                                " as report.md", ""),
                   changes=changes, window_start=0.0)
    assert check["verdict"] == "verified", check


def test_prose_numbering_without_claims_stays_unchecked():
    """A numbered list of context sentences (no paths, no format words,
    no test/calc vocabulary) changes nothing: the aggregate is whatever
    the single-subject check says."""
    desc = ("Notes for the write-up:\n"
            "1. The method was validated against the 2026-09-22 run\n"
            "2. The protocol follows the wave-10 decision\n"
            "3. Numbers carry their units")
    check = _check("Summarise the methodology", desc, changes=[],
                   window_start=0.0)
    assert check["verdict"] == "unchecked", check


def test_a_single_numbered_item_is_not_multi_part():
    changes = [{"path": "delfin/report.py", "ts": 100.0, "created": False}]
    desc = "1. Correct the sort in delfin/report.py"
    check = _check("Fix the ordering bug", desc, changes=changes,
                   window_start=0.0)
    assert check["verdict"] == "verified", check


def test_lettered_parts_are_detected():
    desc = ("a) Correct the sort in delfin/report.py\n"
            "b) Add the regression test tests/test_report_order.py")
    changes = [{"path": "delfin/report.py", "ts": 100.0, "created": False}]
    check = _check("Fix the ordering bug", desc, changes=changes,
                   window_start=0.0)
    assert check["verdict"] == "unmet", check
    assert "tests/test_report_order.py" in check["note"], check


def test_lettered_parts_with_evidence_verify():
    desc = ("a) Correct the sort in delfin/report.py\n"
            "b) Add the regression test tests/test_report_order.py")
    changes = [
        {"path": "delfin/report.py", "ts": 100.0, "created": False},
        {"path": "tests/test_report_order.py", "ts": 101.0,
         "created": True},
    ]
    check = _check("Fix the ordering bug", desc, changes=changes,
                   window_start=0.0)
    assert check["verdict"] == "verified", check


def test_a_dash_bullet_list_is_not_judged_as_parts():
    """The documented scope: dash bullets are context notes here. A
    bullet list beside one written file must not turn the task unmet."""
    desc = ("Use the gate for tests:\n"
            "- pytest via the gate, never directly\n"
            "- do not push, the integrator collects")
    changes = [{"path": "delfin/report.py", "ts": 100.0, "created": False}]
    check = _check("Fix the ordering bug", desc, changes=changes,
                   window_start=0.0)
    assert check["verdict"] == "verified", check


def test_a_test_part_needs_a_test_run():
    """A part that says "run the tests" is evidence-demanding in its own
    right: with an empty test ledger (no runs recorded) the task is not
    done. With no ledger passed at all the part is unchecked -- the
    honest unknown -- and so is the task."""
    desc = ("1. Correct the sort in delfin/report.py\n"
            "2. Run the test suite")
    changes = [{"path": "delfin/report.py", "ts": 100.0, "created": False}]
    check = _check("Fix the ordering bug", desc, changes=changes,
                   tests=[], window_start=0.0)
    assert check["verdict"] == "unmet", check
    assert "2" in check["note"], check


def test_a_test_part_with_a_green_run_verifies():
    desc = ("1. Correct the sort in delfin/report.py\n"
            "2. Run the test suite")
    changes = [{"path": "delfin/report.py", "ts": 100.0, "created": False}]
    tests = [{"tool": "run_tests", "command": "tests/", "exit_code": 0,
              "status": "ok", "passed": 4, "failed": 0, "ts": 100.0}]
    check = _check("Fix the ordering bug", desc, changes=changes,
                   tests=tests, window_start=0.0)
    assert check["verdict"] == "verified", check


def test_no_description_no_parts_no_change():
    """The pinned old behaviour: without enumerated parts the check is
    exactly what it was."""
    check = _check("Fix the ordering bug in delfin/report.py",
                   changes=[], window_start=0.0)
    assert check["verdict"] == "unmet", check
    assert check["kind"] == "path_unwritten", check
