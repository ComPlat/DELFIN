"""Regression tests for submit_delfin.sh wrapper behavior.

These tests guard against set -e bugs in the post-ORCA cleanup phase that
previously caused the wrapper to die silently, leaving incomplete ORCA
output files in the submit directory and SLURM reporting FAILED(1) even
when ORCA itself finished normally. Root causes fixed in commit 437106f:

1. `wait "$_DELFIN_PID"` triggered set -e on ORCA failure, killing the
   wrapper before the final sync could run.
2. `_postprocess_nmr_if_needed` used `grep -qi 'EPRNMR' || return`, and
   `return` without an argument propagated grep's non-zero exit code,
   which re-triggered set -e at the caller.
"""

import subprocess
from pathlib import Path


_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "delfin"
    / "submit_templates"
    / "submit_delfin.sh"
)


def _run_bash(script: str) -> subprocess.CompletedProcess:
    """Run a bash snippet and return the CompletedProcess (no raise on error)."""
    return subprocess.run(
        ["bash", "-c", script],
        capture_output=True,
        text=True,
        timeout=10,
    )


def test_wait_failure_does_not_kill_wrapper_post_orca():
    """Reproduces the ORCA-FAILED case: wait returns non-zero, wrapper must
    still reach the completion marker (and propagate the real exit code)."""
    script = """
    set -euo pipefail

    # Simulate local_runner that fails (ORCA crashed)
    (exit 1) &
    _DELFIN_PID=$!

    # The fix: set +e before wait so errexit doesn't kill the wrapper
    set +e
    wait "$_DELFIN_PID"
    EXIT_CODE=$?

    echo "Exit Code: $EXIT_CODE"
    echo "Copied N result files (job-end, profile=full)"
    echo "DELFIN WRAPPER COMPLETED (exit $EXIT_CODE)"
    exit $EXIT_CODE
    """
    result = _run_bash(script)
    assert result.returncode == 1, "wrapper must propagate ORCA's exit code"
    assert "Exit Code: 1" in result.stdout
    assert "Copied N result files" in result.stdout, \
        "final sync must run even when ORCA failed"
    assert "DELFIN WRAPPER COMPLETED (exit 1)" in result.stdout, \
        "completion marker must appear even when ORCA failed"


def test_old_broken_wait_behavior_is_what_we_fixed():
    """Confirms the pre-fix code path actually reproduced the bug, so this
    test battery is guarding a real regression."""
    script = """
    set -euo pipefail

    (exit 1) &
    _DELFIN_PID=$!

    # OLD broken: no set +e — set -e kills the script at wait
    wait "$_DELFIN_PID"
    EXIT_CODE=$?

    echo "Copied N result files"
    echo "DELFIN WRAPPER COMPLETED"
    exit $EXIT_CODE
    """
    result = _run_bash(script)
    assert result.returncode == 1
    # The broken version never reaches any of these lines
    assert "Copied N result files" not in result.stdout
    assert "DELFIN WRAPPER COMPLETED" not in result.stdout


def test_postprocess_nmr_returns_zero_for_non_nmr_orca(tmp_path):
    """The real _postprocess_nmr_if_needed function must return 0 when the
    ORCA input is not an NMR job, otherwise set -e at the caller will
    kill the wrapper before the final sync."""
    inp_file = tmp_path / "job.inp"
    inp_file.write_text("! B3LYP def2-SVP Opt\n* xyz 0 1\nH 0 0 0\n*\n")

    # Extract and run the actual function from submit_delfin.sh
    script = f"""
    set -euo pipefail

    # Source just the function definition by extracting it from the script.
    eval "$(awk '/^_postprocess_nmr_if_needed\\(\\) \\{{/,/^\\}}$/' {_SCRIPT})"

    export DELFIN_MODE=orca
    export DELFIN_INP_FILE={inp_file}
    export PYTHON_BIN=/usr/bin/python3

    cd {tmp_path}

    # If _postprocess_nmr_if_needed returns non-zero, set -e would kill
    # the shell here. The test passes only if we reach the echo below.
    _postprocess_nmr_if_needed
    echo "post-process returned without aborting the wrapper"
    """
    result = _run_bash(script)
    assert result.returncode == 0, \
        f"_postprocess_nmr_if_needed triggered set -e. stderr={result.stderr}"
    assert "post-process returned without aborting the wrapper" in result.stdout


def test_postprocess_nmr_returns_zero_for_missing_inp(tmp_path):
    """Even with a missing .inp, the function must return cleanly."""
    script = f"""
    set -euo pipefail

    eval "$(awk '/^_postprocess_nmr_if_needed\\(\\) \\{{/,/^\\}}$/' {_SCRIPT})"

    export DELFIN_MODE=orca
    export DELFIN_INP_FILE=/nonexistent/path.inp
    export PYTHON_BIN=/usr/bin/python3

    cd {tmp_path}
    _postprocess_nmr_if_needed
    echo "OK"
    """
    result = _run_bash(script)
    assert result.returncode == 0
    assert "OK" in result.stdout


def test_submit_script_has_set_plus_e_before_wait():
    """Guard against accidental removal of the set +e fix."""
    content = _SCRIPT.read_text()
    # Find the wait line and confirm set +e appears within a reasonable
    # window before it (same section, not 200 lines away)
    lines = content.splitlines()
    wait_idx = next(
        i for i, line in enumerate(lines)
        if 'wait "$_DELFIN_PID"' in line
    )
    window = "\n".join(lines[max(0, wait_idx - 10):wait_idx])
    assert "set +e" in window, \
        "set +e must be set shortly before `wait $_DELFIN_PID` to prevent " \
        "set -e from killing the wrapper on ORCA failure"


def test_submit_script_postprocess_nmr_has_explicit_return_0():
    """Guard against reintroducing the `grep || return` pattern."""
    content = _SCRIPT.read_text()
    # The function body should use `if ... then return 0` form, not the
    # fragile `[ ... ] && return` or `grep ... || return` idioms.
    func_start = content.index("_postprocess_nmr_if_needed()")
    func_end = content.index("\n}\n", func_start)
    func_body = content[func_start:func_end]
    assert "|| return\n" not in func_body, \
        "`|| return` without explicit 0 propagates the failing command's " \
        "exit code and kills the wrapper via set -e"
    assert "return 0" in func_body, \
        "_postprocess_nmr_if_needed must explicitly return 0 on early exits"
