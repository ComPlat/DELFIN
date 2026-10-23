#!/usr/bin/env python3
"""Three calculations in three states, for a question that is three
independent questions and one answer.

"Look at each of these runs on its own and tell me where they stand" is
the shape of a fan-out: every folder is a self-contained audit, nothing
in one depends on another, and the answer is the synthesis. It is the
case `orchestrate` exists for (parallel stage, barrier, synthesis stage)
and, before this task, no benchmark ever asked for it.

The three states are chosen so that a lazy read is wrong:

  runs/succeeded   converged, exit code 0, an energy in the output
  runs/failed      SCF did not converge, exit code 1025 -- an exit-code
                   file exists, so the index says completed: true, and
                   a reader who stops there calls it a result
  runs/running     a run log and no exit code: still going, or crashed

Exit non-zero on any failure: a task whose precondition was not built is
reported as unmeasured, never as a model failure.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

_CONTROL = """functional = {functional}
basis_set = def2-SVP
solvent = DMF
charge = 0
"""

_OUT_OK = """
                                 *****************
                                 * O   R   C   A *
                                 *****************

Your calculation uses the libint2 library
-------------------------------------------------------------------------------
                       SCF CONVERGED AFTER  14 CYCLES
-------------------------------------------------------------------------------
FINAL SINGLE POINT ENERGY      -613.4171289012

                    ***  THE OPTIMIZATION HAS CONVERGED  ***

                             ****ORCA TERMINATED NORMALLY****
TOTAL RUN TIME: 0 days 0 hours 41 minutes 12 seconds 300 msec
"""

_OUT_FAILED = """
                                 *****************
                                 * O   R   C   A *
                                 *****************

-------------------------------------------------------------------------------
                       SCF NOT CONVERGED AFTER 125 CYCLES
-------------------------------------------------------------------------------
Error: SCF not converged in 125 iterations.
       The wavefunction is not converged. Increase MaxIter or use a different
       SCF algorithm (e.g. ! SlowConv / ! KDIIS / ! SOSCF).

ORCA finished by error termination in SCF
"""

_RUN_LOG_RUNNING = """[2026-09-10 05:12:03] delfin: job started (ORCA 6.0.1, 16 cores)
[2026-09-10 05:12:04] delfin: step 1/3 optimisation
[2026-09-10 06:40:51] delfin: step 1/3 optimisation still running (cycle 38)
"""


def _folder(root: Path, name: str, functional: str) -> Path:
    d = root / name
    d.mkdir(parents=True)
    (d / "CONTROL.txt").write_text(_CONTROL.format(functional=functional))
    (d / "run.inp").write_text(f"! {functional} def2-SVP Opt\n* xyz 0 1\nO 0 0 0\nH 0 0 1\nH 0 1 0\n*\n")
    return d


def build(root: Path) -> None:
    d = _folder(root, "succeeded", "PBE0")
    (d / "run.out").write_text(_OUT_OK)
    (d / ".exit_code_0").write_text("")
    (d / "DELFIN_Data.json").write_text(json.dumps({
        "functional": "PBE0", "basis_set": "def2-SVP", "solvent": "DMF",
        "electronic_energy": -613.4171289012, "completed": True,
    }, indent=1))

    d = _folder(root, "failed", "B3LYP")
    (d / "run.out").write_text(_OUT_FAILED)
    (d / ".exit_code_1025").write_text("")
    (d / "delfin_run.log").write_text(
        "[2026-09-09 22:01:10] delfin: job started\n"
        "[2026-09-09 22:47:33] delfin: ORCA error termination in SCF (exit 1025)\n")

    d = _folder(root, "running", "PBE0")
    (d / "delfin_run.log").write_text(_RUN_LOG_RUNNING)


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    root = ws / "runs"
    if root.exists():
        print(f"refusing to overwrite {root}", file=sys.stderr)
        return 1
    build(root)
    # The precondition, checked: the indexer must see all three and must
    # tell them apart, or the task measures the fixture and not the model.
    try:
        from delfin.doc_server.calc_indexer import _scan_calc_dir
        recs = {r["rel_path"].split("/")[-1]: r for r in _scan_calc_dir(root, "calc", quiet=True)}
    except Exception as exc:                                   # pragma: no cover
        print(f"the indexer could not read the fixture: {exc}", file=sys.stderr)
        return 1
    want = {"succeeded": "succeeded", "failed": "failed", "running": "running or crashed"}
    for name, prefix in want.items():
        got = str(recs.get(name, {}).get("outcome", ""))
        if not got.startswith(prefix):
            print(f"{name}: outcome {got!r}, expected to start with {prefix!r}", file=sys.stderr)
            return 1
    print(f"prepared {root} — three runs: succeeded, failed (exit 1025), running")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
