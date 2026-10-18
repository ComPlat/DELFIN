"""CI1 adversarial tests (nacht-s6, reviewer): rollup/log edge cases.

These attack summarize/failing_tests/next_step contracts the builder's
17 red cases do not cover: conclusions beyond SUCCESS/FAILURE, missing
keys, case sensitivity, boundary attempts, hostile shapes.
They run against delfin.agent.ci_watch (built in phase 2).
"""
import pytest

import delfin.agent.ci_watch as ci_watch


def _rollup(entries, head="abc1234567890"):
    return {"headRefOid": head, "statusCheckRollup": entries}


def _check(name, status="COMPLETED", conclusion="SUCCESS"):
    return {"name": name, "status": status, "conclusion": conclusion}


# ---------------------------------------------------------------- summarize

@pytest.mark.parametrize("conclusion", ["CANCELLED", "TIMED_OUT",
                                        "STARTUP_FAILURE", "ACTION_REQUIRED"])
def test_summarize_non_success_conclusions_are_not_green(conclusion):
    entries = [_check("lint"), _check("pytest", conclusion=conclusion)]
    data = ci_watch.summarize(_rollup(entries))
    assert data["state"] == "red"
    assert "pytest" in data["failing"]


def test_summarize_skipped_is_not_a_failure():
    entries = [_check("lint", conclusion="SKIPPED")]
    data = ci_watch.summarize(_rollup(entries))
    assert data["state"] == "green"


def test_summarize_queued_is_running():
    entries = [_check("pytest", status="QUEUED")]
    assert ci_watch.summarize(_rollup(entries))["state"] == "running"


def test_summarize_missing_rollup_key_is_none_not_crash():
    data = ci_watch.summarize({"headRefOid": "abc1234567890"})
    assert data["state"] == "none"
    assert data["failing"] == []


def test_summarize_none_rollup_json_is_none_state():
    assert ci_watch.summarize(None)["state"] == "none"


def test_summarize_missing_head_oid_reports_it():
    data = ci_watch.summarize({"statusCheckRollup": [_check("lint")]})
    assert data["state"] == "green"
    assert data["head_sha"] in (None, "")


def test_summarize_duplicate_failing_names_deduped():
    entries = [_check("pytest", conclusion="FAILURE"),
               _check("pytest", conclusion="FAILURE")]
    data = ci_watch.summarize(_rollup(entries))
    assert data["failing"] == ["pytest"]


def test_summarize_stale_failing_list_is_empty():
    entries = [_check("pytest", conclusion="FAILURE")]
    data = ci_watch.summarize(_rollup(entries, head="old"),
                              current_head="new")
    assert data["failing"] == []


# ------------------------------------------------------------- failing_tests

def test_failing_tests_indented_and_crlf_lines_match():
    log = "  FAILED tests/test_x.py::test_a - boom\r\n" \
          "\tERROR tests/test_x.py::test_b - nope\r\n"
    out = ci_watch.failing_tests(log)
    assert any("test_a" in line for line in out)
    assert any("test_b" in line for line in out)


def test_failing_tests_unicode_lookalike_not_matched():
    # A Cyrillic A (U+0410) masquerading as the latin A must NOT match FAILED.
    log = "FАILED tests/test_x.py::test_evil - looks like FAILED"
    assert ci_watch.failing_tests(log) == []


def test_failing_tests_empty_log():
    assert ci_watch.failing_tests("") == []


def test_failing_tests_none_log():
    assert ci_watch.failing_tests(None) == []


# ----------------------------------------------------------------- next_step

def test_next_step_boundary_attempts_two_is_repair():
    assert ci_watch.next_step("red", 2, {"max_repairs": 3})[0] == "repair"


def test_next_step_missing_limits_key_uses_default():
    step, _ = ci_watch.next_step("red", 0, {})
    assert step in ("repair", "wait")


def test_next_step_negative_attempts_is_repair_not_crash():
    assert ci_watch.next_step("red", -1, {"max_repairs": 3})[0] == "repair"


def test_next_step_unknown_state_asks_person():
    assert ci_watch.next_step("bogus", 0, {"max_repairs": 3})[0] == "ask_person"


def test_next_step_returns_reason_string_always():
    for state in ("green", "red", "running", "none"):
        step, reason = ci_watch.next_step(state, 0, {"max_repairs": 3})
        assert isinstance(reason, str) and reason
