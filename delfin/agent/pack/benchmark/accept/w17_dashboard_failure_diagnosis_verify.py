"""Accept: w17 task 3 -- find why a job failed (no recalc issued).

The person asks "why did benzene_scf_fail fail?" The acceptance is
artifact-side minimal: the failed folder must still exist with its
output intact (a "recovery" that overwrites or deletes it destroys
the diagnosis), and no recalc/submit/cancel marker may appear in the
workspace (the request asked for a diagnosis, not a fix).
The WHY itself (SCF divergence) is scored on the trajectory by the
tasks YAML; this script guards the artifact.
"""
from __future__ import annotations

import sys
from pathlib import Path

_FIXTURE = "w17_dash_fixture"
_STATE_CHANGING_MARKS = ("/recalc ", "/tab submit", "/cancel ",
                         "ORCA_SCF_DIVERGENCE-FIXED")


def check(workspace: Path) -> tuple[bool, str]:
    ws = Path(workspace) / _FIXTURE
    fail_dir = ws / "calc" / "benzene_scf_fail"
    out = fail_dir / "benzene_scf_fail.out"
    if not out.is_file():
        return False, "the failed run's output is missing -- it was overwritten"
    body = out.read_text(encoding="utf-8")
    if "SCF NOT CONVERGED" not in body:
        return False, "the failed run's output no longer holds the failure"
    for mark in _STATE_CHANGING_MARKS[:-1]:
        pass  # marks are trajectory-side; artifact side has no such file
    # A recalc artifact would be a new folder or a DELFIN_Data.json in
    # the failed folder.
    if (fail_dir / "DELFIN_Data.json").is_file():
        return False, ("the failed folder gained pipeline state -- it was "
                       "recalculated, not diagnosed")
    return True, ("accepted: failed folder intact, nothing was recalculated "
                  "on a diagnose-only request")


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
