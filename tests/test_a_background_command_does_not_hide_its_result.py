"""A background command that hides its own result is told so.

Report 20260915-125310 started the test suite as
``pytest -q -m "not slow" 2>&1 | tail -30`` -- bash_output showed nothing for
thirty minutes, because tail prints only at the end -- and then as
``pytest ... > log; echo "EXIT=$?" >> log``, whose exit code was echo's: the
suite failed 25 tests and the job finished with exit code 0.
"""

from __future__ import annotations

from delfin.agent.api_client import _background_command_note


def test_a_tail_pipe_is_named():
    note = _background_command_note(
        {"command": 'python -m pytest -q -m "not slow" 2>&1 | tail -30'})
    assert "| tail" in note and "bash_output" in note


def test_a_trailing_echo_is_named():
    note = _background_command_note({"command": (
        'python -m pytest -q > pytest_fast.log 2>&1; '
        'echo "EXIT=$?" >> pytest_fast.log')})
    assert "echo" in note and "exit_code" in note


def test_a_plain_command_gets_no_note():
    assert _background_command_note(
        {"command": 'python -m pytest -q -m "not slow"'}) == ""
    assert _background_command_note(
        {"command": "echo start > log && python run.py >> log 2>&1"}) == ""
