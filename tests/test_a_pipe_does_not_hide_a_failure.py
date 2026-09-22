"""A filtered pipe whose exit code hides a failure is named.

A pipe exits with its last command's code. Reports 20260915-090037 and
-091158 ran ``pytest ... 2>&1 | tail -3`` and got exit code 0 with the
failures on screen; report 20260915-132613 pushed through ``| tail -3``,
where a rejected push exits 0 too.
"""

from __future__ import annotations

import json

from delfin.agent.api_client import _pipe_exit_note


def _result(output, exit_code=0):
    return json.dumps({"exit_code": exit_code, "stdout": output, "stderr": ""})


def test_failing_tests_behind_tail_are_named():
    note = _pipe_exit_note(
        {"command": "python -m pytest tests/test_x.py -q 2>&1 | tail -3"},
        _result("FAILED tests/test_x.py::test_a - assert 1 == 2\n"
                "1 failed, 4 passed in 0.31s\n"))
    assert "exit code 0" in note and "pipefail" in note


def test_a_rejected_push_behind_tail_is_named():
    note = _pipe_exit_note(
        {"command": "git push origin main 2>&1 | tail -3"},
        _result(" ! [rejected]        main -> main (fetch first)\n"
                "error: failed to push some refs to 'github.com:o/r.git'\n"))
    assert note


def test_a_green_run_behind_tail_stays_quiet():
    assert _pipe_exit_note(
        {"command": "python -m pytest tests/test_x.py -q 2>&1 | tail -3"},
        _result("5 passed in 0.20s\n")) == ""


def test_pipefail_keeps_the_real_exit_code():
    assert _pipe_exit_note(
        {"command": "set -o pipefail; python -m pytest -q 2>&1 | tail -3"},
        _result("1 failed, 4 passed in 0.31s\n")) == ""


def test_a_failure_the_exit_code_already_reports_gets_no_note():
    assert _pipe_exit_note(
        {"command": "python -m pytest -q 2>&1 | tail -3"},
        _result("1 failed, 4 passed in 0.31s\n", exit_code=1)) == ""
    assert _pipe_exit_note(
        {"command": "python -m pytest -q"},
        _result("1 failed, 4 passed in 0.31s\n")) == ""
