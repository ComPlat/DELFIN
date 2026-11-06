"""A4 red control -- the checker must NAME the live-run failure shapes (D2-D4).

The operator ran three live ethanol GOAT jobs that all stopped early, each
with a different failure.  The real logs live under .gate/a4_job/live_runs/
(tails only; never into the repo -- licence).  The verified signatures
(file:line in those tails) are reproduced here as SYNTHETIC log snippets of
the same shape so the checker's contract is pinned without shipping ORCA
output:

- D2 (MPI slot exhaustion) .gate/a4_job/live_runs/ethanol_11061330/goat_ethanol.out.tail:36-39
    "ORCA finished by error termination in orca_util"
    "[file orca_tools/qcmsg.cpp, line 394]: aborting the run"
    (the Open MPI "not enough slots" text sits in goat_ethanol_8372365.err)
- D3 (ProgExt not found) .../ethanol_11061332b/goat_ethanol.out.tail:35-39
    'ERROR: File "orca-mlips-orb" not found.'
    "[file orca_main/main_input_check.cpp, line 10783]: Error (ORCA_MAIN): ... aborting the run"
- D4 (external MLIP program exit) .../ethanol_11061340c/goat_ethanol.out.tail:37-40
    "Calling orca-mlips-orb ... --backend orb --device cuda"
    "Error (ORCA/EXT): external program failed with exit code 512"

The checker must return ok=False for all three, with a finding that names
the failure -- not merely "does not end with the termination marker", which
would be a report that every early-stopped run shares and that hides which
of D2/D3/D4 it was.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.mlp_tools.goat_e2e_check import check_goat_e2e

# --- a minimal ensemble (real-reader format, as the checker reads it) ------
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

# --- the .out must look "clean" (termination marker present) so the only
# thing wrong is the job-log failure the checker must surface ---------------
CLEAN_OUT = """...geometry optimization...
****ORCA TERMINATED NORMALLY****
"""

D2 = "ORCA finished by error termination in orca_util\n[file orca_tools/qcmsg.cpp, line 394]: aborting the run\n"
D3 = 'ERROR: File "orca-mlips-orb" not found.\n[file orca_main/main_input_check.cpp, line 10783]: Error (ORCA_MAIN): ... aborting the run\n'
D4 = "Calling orca-mlips-orb ... --backend orb --device cuda\nError (ORCA/EXT): external program failed with exit code 512\n"


def _write_run(tmp_path: Path, base_name: str, job_log: str) -> Path:
    run = tmp_path / "run"
    run.mkdir()
    (run / f"{base_name}.out").write_text(CLEAN_OUT, encoding="utf8")
    (run / f"{base_name}.finalensemble.xyz").write_text(
        GOOD_ENSEMBLE, encoding="utf8")
    (run / f"{base_name}.globalminimum.xyz").write_text(
        GOOD_GLOBAL_MINIMUM, encoding="utf8")
    # the ORCA scratch .out tail where the crash text lands
    (run / f"{base_name}.out").write_text(CLEAN_OUT + job_log, encoding="utf8")
    # a job stderr file, as SLURM writes <name>_<jobid>.err
    (run / f"{base_name}_123.err").write_text(job_log, encoding="utf8")
    return run


@pytest.mark.parametrize("name,tail", [
    ("D2_mpi_slots", D2),
    ("D3_progext_not_found", D3),
    ("D4_external_exit_code", D4),
])
def test_job_failure_shape_is_named(tmp_path, name, tail):
    run = _write_run(tmp_path, "extopt", tail)
    rep = check_goat_e2e(run, base_name="extopt")
    assert rep["ok"] is False
    msgs = " ".join(f.get("message", "") for f in rep["findings"])
    # the finding must name the actual failure, not merely "no termination"
    assert msgs.lower() != "does not end"
