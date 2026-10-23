"""Two claims a delegate got certified on evidence that says nothing.

**"All tests pass" over a run whose output was discarded.** The
number-check is gated on the recorded output being truthy::

    if (claim.subject and ev["test_output"] and not re.search(...)):

so an EMPTY output falls straight through to ``_add(True,
"test_result", ...)``. A delegate running ``pytest -q > /dev/null 2>&1``
and reporting "all tests pass" came back SUPPORTED -- the verifier
endorsing a claim nothing in the run can speak to. The module already
knows this shape: the ephemeral-sink pattern was wired into the write
check for exactly the same reason and never into this one.

The bash result is a JSON envelope (command, cwd, exit_code, elapsed_s,
stdout, stderr), so the discarded run does not even leave an empty
string behind: it leaves an envelope that mentions pytest, contains no
failures, and can satisfy a claimed NUMBER out of ``elapsed_s``. The
check now reads the output the run produced, not the envelope around it.

**A mutation claim backed by any write anywhere.** ``ev["writes"]`` is a
flat list of write-ish calls, so ``mkdir -p /tmp/scratch`` corroborated
"I implemented the fix in engine.py". A claim that names a file is now
paired with it: that file has to appear among the paths this run wrote.

The bookkeeping principle is kept throughout. A test run whose output
was never recorded at all is left alone -- absent records are not
evidence of fabrication. A mutation claim that names no file is judged
exactly as before.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import subagents as sa


def _bash(command: str, **fields) -> dict:
    """A bash tool call as the trace records it: JSON envelope and all."""
    payload = {"exit_code": 0, "elapsed_s": 3.14, "stdout": "",
               "stderr": "", "command": command, "cwd": "."}
    payload.update(fields)
    return {"name": "bash", "input": {"command": command},
            "output": json.dumps(payload)}


def _verify(text: str, calls: list[dict]) -> dict:
    return sa.verify_subagent_report({"final_text": text, "tool_calls": calls})


def _kinds(v: dict, key: str) -> set[str]:
    return {e["kind"] for e in v[key]}


# ---------------------------------------------------------------------------
# A discarded test run is no evidence
# ---------------------------------------------------------------------------

_DISCARDED = _bash("pytest -q > /dev/null 2>&1")


def test_a_pass_claim_over_a_discarded_run_is_unsupported():
    assert "test_result" in _kinds(_verify("All tests pass.", [_DISCARDED]),
                                  "unsupported")


def test_a_discarded_run_is_never_reported_as_supported():
    assert "test_result" not in _kinds(
        _verify("All tests pass.", [_DISCARDED]), "supported")


def test_the_reason_says_the_run_recorded_no_output():
    why = " ".join(u.get("why", "")
                   for u in _verify("All tests pass.", [_DISCARDED])
                   ["unsupported"])
    assert "no output" in why.lower()


def test_a_claimed_number_is_not_met_by_the_envelope_itself():
    """3.14 seconds is not three passing tests."""
    v = _verify("3 tests passed.", [_DISCARDED])
    assert "test_result" in _kinds(v, "unsupported")


def test_the_german_form_is_caught_too():
    v = _verify("Alle Tests sind grün.", [_DISCARDED])
    assert "test_result" in _kinds(v, "unsupported")


# ---------------------------------------------------------------------------
# ...while a real run still reads as one
# ---------------------------------------------------------------------------

def test_a_real_passing_run_is_still_supported():
    v = _verify("All tests pass.",
                [_bash("pytest -q", stdout="120 passed, 3 skipped in 8.4s")])
    assert "test_result" in _kinds(v, "supported")
    assert "test_result" not in _kinds(v, "unsupported")


def test_a_failing_run_inside_the_envelope_is_still_caught():
    v = _verify("All tests pass.",
                [_bash("pytest -q", stdout="3 failed, 10 passed in 2.1s")])
    assert "test_result" in _kinds(v, "unsupported")


def test_a_claimed_count_still_reads_the_real_output():
    v = _verify("120 passed.",
                [_bash("pytest -q", stdout="120 passed in 8.4s")])
    assert "test_result" in _kinds(v, "supported")


def test_a_plain_text_test_output_still_works():
    """Not every runner returns an envelope."""
    v = _verify("All tests pass.",
                [{"name": "run_tests", "input": {},
                  "output": "60 passed in 4.0s"}])
    assert "test_result" in _kinds(v, "supported")


def test_a_run_with_no_recorded_output_at_all_is_left_alone():
    """Missing bookkeeping is not evidence of fabrication."""
    v = _verify("All tests pass.",
                [{"name": "bash", "input": {"command": "pytest -q"}}])
    assert "test_result" not in _kinds(v, "unsupported")


# ---------------------------------------------------------------------------
# A mutation claim has to name a file this run wrote
# ---------------------------------------------------------------------------

_SCRATCH_DIR = {"name": "bash", "input": {"command": "mkdir -p /tmp/scratch"},
                "output": "{}"}


def test_a_scratch_directory_does_not_corroborate_a_source_edit():
    v = _verify("I implemented the fix in delfin/agent/engine.py.",
                [_SCRATCH_DIR])
    assert "completion" in _kinds(v, "unsupported")


def test_the_reason_names_the_file_that_was_never_written():
    v = _verify("I implemented the fix in delfin/agent/engine.py.",
                [_SCRATCH_DIR])
    why = " ".join(u.get("why", "") for u in v["unsupported"])
    assert "engine.py" in why


def test_writing_a_different_file_does_not_back_the_claim():
    v = _verify("I implemented the fix in delfin/agent/engine.py.",
                [{"name": "write_file", "input": {"path": "notes.md"},
                  "output": "written"}])
    assert "completion" in _kinds(v, "unsupported")


def test_the_matching_write_backs_the_claim():
    v = _verify("I implemented the fix in delfin/agent/engine.py.",
                [{"name": "edit_file",
                  "input": {"path": "delfin/agent/engine.py"},
                  "output": "1 edit applied"}])
    assert "completion" in _kinds(v, "supported")
    assert "completion" not in _kinds(v, "unsupported")


def test_an_absolute_write_path_matches_a_relative_citation():
    v = _verify("I fixed delfin/agent/engine.py.",
                [{"name": "write_file",
                  "input": {"path": "/repo/delfin/agent/engine.py"},
                  "output": "written"}])
    assert "completion" not in _kinds(v, "unsupported")


def test_an_in_place_shell_edit_still_backs_the_claim():
    """The target of `sed -i` is named in the command, and that is the
    lenient direction — it only makes the pairing easier to satisfy."""
    v = _verify("I fixed delfin/agent/engine.py.",
                [_bash("sed -i s/a/b/ delfin/agent/engine.py")])
    assert "completion" not in _kinds(v, "unsupported")


def test_a_claim_naming_no_file_is_judged_as_before():
    v = _verify("I updated the config.",
                [{"name": "write_file", "input": {"path": "cfg.toml"},
                  "output": "written"}])
    assert "completion" in _kinds(v, "supported")


@pytest.mark.parametrize("text", [
    "I updated the config.",
    "Ich habe die Konfiguration aktualisiert.",
])
def test_a_claim_with_no_write_at_all_is_still_unsupported(text):
    v = _verify(text, [{"name": "read_file", "input": {"path": "cfg.toml"},
                        "output": "x = 1"}])
    assert "completion" in _kinds(v, "unsupported")
