"""Three ways a delegate's report passed a check it should have failed.

The cross-check exists so a parent does not take a delegate's word for
its own work. Each of these let the word through.

**A discarded write counted as a change.** The bash write pattern matches
any redirect, so ``python -c '...' > /dev/null`` -- a command whose whole
point is to throw the output away -- registered as file-writing evidence.
A claim of "updated the config" was then reported as backed.

**A failing suite backed a passing claim.** A test-result claim was
confirmed as soon as ANY test-running call was in the trace. The recorded
output was consulted only to look for a claimed NUMBER. So a delegate
that ran pytest, got "3 failed, 10 passed", and reported "all tests pass"
came back supported -- the verifier actively endorsing the false claim it
was built to catch.

**An attempted read grounded a citation.** Evidence paths are taken from
each call's INPUT arguments, whatever the call returned. A read that was
refused, or that found nothing, still put the path in the observed set --
which grounds the delegate's file citations and, because delegate
evidence is merged upward, tells the PARENT it observed a file nobody
read.

The module's principle is kept: missing bookkeeping is never treated as
evidence of fabrication. A path is dropped only when the recorded output
is there and says the read failed. When nothing was recorded, nothing is
concluded.
"""

from __future__ import annotations

import pytest

from delfin.agent import subagents as sa


def _ev(calls):
    return sa.collect_report_evidence({"final_text": "x", "tool_calls": calls})


# ---------------------------------------------------------------------------
# A discarded write is not a change
# ---------------------------------------------------------------------------

def test_a_redirect_to_dev_null_is_not_a_write():
    ev = _ev([{"name": "bash", "input": {"command": "pytest -q > /dev/null"},
               "output": "ok"}])
    assert ev["writes"] == []


@pytest.mark.parametrize("sink", ["/dev/null", "/dev/stdout", "/dev/stderr"])
def test_no_ephemeral_sink_counts(sink):
    ev = _ev([{"name": "bash", "input": {"command": f"echo hi > {sink}"},
               "output": ""}])
    assert ev["writes"] == [], sink


def test_a_real_redirect_still_counts():
    ev = _ev([{"name": "bash", "input": {"command": "echo hi > out.txt"},
               "output": ""}])
    assert ev["writes"]


def test_a_write_tool_still_counts():
    ev = _ev([{"name": "write_file", "input": {"path": "a.py"},
               "output": "written"}])
    assert ev["writes"]


def test_a_mutation_claim_is_not_backed_by_a_discarded_write():
    v = sa.verify_subagent_report({
        "final_text": "I updated the config.",
        "tool_calls": [{"name": "bash",
                        "input": {"command": "cat cfg > /dev/null"},
                        "output": ""}],
    })
    assert any(u["kind"] == "completion" for u in v["unsupported"])


# ---------------------------------------------------------------------------
# A failing suite does not back a passing claim
# ---------------------------------------------------------------------------

def _run(output: str) -> dict:
    return sa.verify_subagent_report({
        "final_text": "All tests pass.",
        "tool_calls": [{"name": "bash", "input": {"command": "pytest -q"},
                        "output": output}],
    })


def test_a_pass_claim_over_a_failing_run_is_unsupported():
    v = _run("3 failed, 10 passed in 2.1s")
    assert any(u["kind"] == "test_result" for u in v["unsupported"])


def test_the_reason_names_the_failure():
    v = _run("3 failed, 10 passed in 2.1s")
    why = " ".join(u.get("why", "") for u in v["unsupported"])
    assert "fail" in why.lower()


@pytest.mark.parametrize("out", [
    "3 failed, 10 passed in 2.1s",
    "FAILED tests/test_a.py::test_b - AssertionError",
    "1 error in 0.3s",
    "=== 2 failed, 1 error ===",
])
def test_every_shape_of_failure_is_caught(out):
    assert any(u["kind"] == "test_result" for u in _run(out)["unsupported"]), out


def test_a_genuinely_passing_run_is_still_supported():
    v = _run("120 passed, 3 skipped in 8.4s")
    assert any(s["kind"] == "test_result" for s in v["supported"])
    assert not any(u["kind"] == "test_result" for u in v["unsupported"])


def test_a_run_with_no_recorded_output_is_left_alone():
    """Missing bookkeeping is not evidence of fabrication."""
    v = sa.verify_subagent_report({
        "final_text": "All tests pass.",
        "tool_calls": [{"name": "bash", "input": {"command": "pytest -q"}}],
    })
    assert not any(u["kind"] == "test_result" for u in v["unsupported"])


def test_the_word_failed_inside_a_passing_summary_is_not_a_failure():
    """'0 failed' is a pass, and so is 'expected failures'."""
    v = _run("50 passed, 0 failed in 3.0s")
    assert not any(u["kind"] == "test_result" for u in v["unsupported"])


# ---------------------------------------------------------------------------
# An attempted read does not ground a citation
# ---------------------------------------------------------------------------

def test_a_refused_read_does_not_count_as_observed():
    ev = _ev([{"name": "read_file", "input": {"path": "secret.py"},
               "output": '{"error": "denied by the folder lock"}'}])
    assert "secret.py" not in ev["files"]


def test_an_empty_read_does_not_count_as_observed():
    ev = _ev([{"name": "read_file", "input": {"path": "gone.py"},
               "output": "(file not found)"}])
    assert "gone.py" not in ev["files"]


def test_a_successful_read_still_counts():
    ev = _ev([{"name": "read_file", "input": {"path": "a.py"},
               "output": "1  import os"}])
    assert "a.py" in ev["files"]


def test_a_call_with_no_recorded_output_still_counts():
    """The payload trace carries inputs without outputs; concluding failure
    from an absent record would discard every path in that whole path."""
    ev = _ev([{"name": "read_file", "input": {"path": "a.py"}}])
    assert "a.py" in ev["files"]


def test_a_citation_is_not_grounded_by_the_attempt():
    v = sa.verify_subagent_report({
        "final_text": "The limit is set in config/limits.py.",
        "tool_calls": [{"name": "read_file",
                        "input": {"path": "config/limits.py"},
                        "output": '{"error": "denied"}'}],
    })
    assert any(u["kind"] == "file_reference" for u in v["unsupported"])


def test_the_worktree_diff_is_still_evidence():
    """A diff summary reports what actually changed inside the child."""
    ev = sa.collect_report_evidence({
        "final_text": "x",
        "tool_calls": [],
        "worktree": {"diff_summary": " delfin/agent/engine.py | 4 ++--"},
    })
    assert any("engine.py" in f for f in ev["files"])
