"""A4 red controls -- job-log failure naming (D2-D4) beyond the builder's own test.

Reviewer s13, on the builder's two follow-up commits (617ec42bf, 4b1f9034e).
The builder's tests/test_a4_e2e_jobfailures.py pins the new behaviour with

    assert msgs.lower() != "does not end"

which joins the FULL finding sentences and compares them to the bare phrase
"does not end" -- a comparison that is always True.  That assertion would
stay green even if the _JOB_FAILURES scan were removed entirely, as long as
some other error finding made ok False.  So the new behaviour (naming WHICH
of D2/D3/D4 killed the run) is not actually pinned.

These controls tighten that:

- A01-A03: each live-failure tail must flip ok=False AND the finding message
  must contain the specific failure CATEGORY (mpi / orca_main / progext /
  external_exit), not merely any error finding.  This is the concrete promise
  of 617ec42bf and the builder's own loose assertion does not detect it.
- H01: the deny-list scan reads every *.err in the run dir (checker line
  ~139, ``run_dir.glob("*.err")``) and matches bare substrings, not ORCA-
  anchored context.  A HEALTHY run whose directory happens to hold an
  unrelated .err file carrying the generic strings must STILL pass (an
  OCCUPIER or xTB stderr dropped next to a finished GOAT job is plausible).
  A deny-list (blacklist) that fires on a good run is worse than one that
  misses, so this is the sshd of the tester's local problem: name it.

Runs under ``check_goat_e2e`` only, on synthetic files built in the test.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.mlp_tools.goat_e2e_check import _JOB_FAILURES, check_goat_e2e

# --- a minimal, well-formed ensemble (real-reader format) -------------------
GOOD_ENSEMBLE = """3
-76.4000000000 converged=true
C  0.0  0.0  0.0
C  1.5  0.0  0.0
O  2.2  1.2  0.0
-76.3990000000 converged=true
C  0.0  0.0  0.0
C  1.5  0.0  0.0
O  2.2  1.2  0.0
-76.3980000000 converged=true
C  0.0  0.0  0.0
C  1.5  0.0  0.0
O  2.2  1.2  0.0
"""

GOOD_GLOBAL_MINIMUM = """3
C  0.0  0.0  0.0
C  1.5  0.0  0.0
O  2.2  1.2  0.0
"""

# a clean .out with the ORCA termination marker and no failure line
CLEAN_OUT = """...geometry optimization...
****ORCA TERMINATED NORMALLY****
"""

# the live-run failure tails, cited in tests/test_a4_e2e_jobfailures.py
D2 = "ORCA finished by error termination in orca_util\n[file orca_tools/qcmsg.cpp, line 394]: aborting the run\n"
D3 = 'ERROR: File "orca-mlips-orb" not found.\n[file orca_main/main_input_check.cpp, line 10783]: Error (ORCA_MAIN): ... aborting the run\n'
D4 = "Calling orca-mlips-orb ... --backend orb --device cuda\nError (ORCA/EXT): external program failed with exit code 512\n"

# the concrete categories the deny-list needles actually fire for each tail
# (the checker names what a needle's tuple second element is).  D2's needle
# set is "ORCA finished by error termination" (orca_util) + "aborting the
# run" (abort); the "not enough slots" MPI text lives in the job .err in the
# real run and is not part of the D2 .out tail.
TAIL_CATEGORY = {
    "A01_D2_mpi_slots": "abort",
    "A02_D3_progext_not_found": "progext",
    "A03_D4_external_exit_code": "external_exit",
}


def _write_healthy_run(tmp_path: Path, extra_err: list[str] = ()) -> Path:
    """A healthy GOAT run directory: good ensemble, clean .out, optional extras.

    ``extra_err``: additional ``*.err`` files written into the SAME directory,
    as an unrelated job's stderr dropped next to a finished GOAT job would be.
    """
    run = tmp_path / "run"
    run.mkdir()
    (run / "extopt.out").write_text(CLEAN_OUT, encoding="utf8")
    (run / "extopt.finalensemble.xyz").write_text(GOOD_ENSEMBLE, encoding="utf8")
    (run / "extopt.globalminimum.xyz").write_text(GOOD_GLOBAL_MINIMUM, encoding="utf8")
    for i, err in enumerate(extra_err):
        (run / f"unrelated_{i}.err").write_text(err, encoding="utf8")
    return run


@pytest.mark.parametrize("name,tail", [
    ("A01_D2_mpi_slots", D2),
    ("A02_D3_progext_not_found", D3),
    ("A03_D4_external_exit_code", D4),
])
def test_failure_tail_names_its_category(tmp_path, name, tail):
    """The finding must NAME the concrete D2/D3/D4 category, not just error."""
    run = _write_healthy_run(tmp_path)
    (run / "extopt.out").write_text(CLEAN_OUT + tail, encoding="utf8")
    rep = check_goat_e2e(run, base_name="extopt")
    assert rep["ok"] is False
    msgs = " ".join(f.get("message", "") for f in rep["findings"]).lower()
    assert TAIL_CATEGORY[name] in msgs, (
        f"finding does not name the {TAIL_CATEGORY[name]} category: {msgs}")

    # every tail must also be found by the deny-list scan (it is driven by it)
    found = [tag for tag, _cat in _JOB_FAILURES if tag.lower() in msgs]
    assert found, f"no known needle matched the {name} tail"


def test_h01_healthy_run_survives_a_benign_matching_err(tmp_path):
    """A good run must still pass when an unrelated .err carries a needle word.

    The scan ``run_dir.glob("*.err")`` matches bare substrings from a deny-
    list.  A healthy finished GOAT run must not be marked failed because a
    neighbour file in the same directory contains the words the deny-list
    hunts (e.g. the string ``aborting the run`` or ``not enough slots`` in a
    different tool's diagnostic).  This is the classic deny-list false
    positive, and it is the direction the operator's A4 checklist flags:
    "which class does the check mean to allow".
    """
    run = _write_healthy_run(
        tmp_path,
        extra_err=[
            # a benign sentence that contains the literal needle substrings
            "safety monitor: aborting the run was never triggered on node 3; "
            "report: not enough slots available anywhere on this partition; "
            f"Error (ORCA_MAIN): (no such message today)\n",
        ],
    )
    rep = check_goat_e2e(run, base_name="extopt")
    assert rep["ok"] is True, (f"healthy run wrongly failed: "
                               f"{[f.get('message') for f in rep['findings']]}")


def test_h02_only_an_error_level_finding_flips_ok(tmp_path):
    """Warnings (GPU evidence absent, stale marker) never flip ok=False."""
    run = _write_healthy_run(tmp_path)
    rep = check_goat_e2e(run, base_name="extopt")
    hard = [f for f in rep["findings"] if f.get("level") == "error"]
    assert rep["ok"] is True
    assert hard == []
