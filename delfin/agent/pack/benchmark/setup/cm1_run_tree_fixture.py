#!/usr/bin/env python3
"""Fixture for the CM1 (wave 26) Code Mode benchmark family.

Builds a run tree of 20 job folders under <workspace>/cm1/run_tree/,
each folder holding a host-neutral recorded ORCA .out excerpt (final
energy line, status markers, method header) plus a small CONTROL.txt
excerpt. Deterministic values (no host names, no user names). Refuses
to overwrite an existing cm1/ tree. No subprocess, no scheduler, no
network.
"""
from __future__ import annotations

import sys
from pathlib import Path

# (folder, method, dispersion, basis, energy_eh, status)
_JOBS = [
    ("job_01", "r2SCAN-3c", "D4", "def2-mTZVP", -402.512341, "converged"),
    ("job_02", "r2SCAN-3c", "D4", "def2-mTZVP", -402.518220, "converged"),
    ("job_03", "r2SCAN-3c", "D4", "def2-mTZVP", -402.514890, "converged"),
    ("job_04", "r2SCAN-3c", "D4", "def2-mTZVP", -402.509912,
     "FAILED: SCF not converged"),
    ("job_05", "r2SCAN-3c", "D4", "def2-mTZVP", -402.521477, "converged"),
    ("job_06", "B3LYP", "D3BJ", "def2-SVP", -401.778210, "converged"),
    ("job_07", "B3LYP", "D3BJ", "def2-SVP", -401.781904, "converged"),
    ("job_08", "r2SCAN-3c", "D4", "def2-mTZVP", -402.516000, "converged"),
    ("job_09", "B3LYP", "D3BJ", "def2-SVP", -401.779801, "converged"),
    ("job_10", "r2SCAN-3c", "D4", "def2-mTZVP", -402.519921, "converged"),
    ("job_11", "r2SCAN-3c", "D4", "def2-mTZVP", -402.513350, "converged"),
    ("job_12", "r2SCAN-3c", "D4", "def2-mTZVP", -402.514771, "converged"),
    ("job_13", "B3LYP", "D3BJ", "def2-SVP", -401.775612, "FAILED: walltime"),
    ("job_14", "r2SCAN-3c", "D4", "def2-mTZVP", -402.520100, "converged"),
    ("job_15", "r2SCAN-3c", "D4", "def2-mTZVP", -402.517850, "converged"),
    ("job_16", "B3LYP", "D3BJ", "def2-SVP", -401.780450, "converged"),
    ("job_17", "r2SCAN-3c", "D4", "def2-mTZVP", -402.515290, "converged"),
    ("job_18", "r2SCAN-3c", "D4", "def2-mTZVP", -402.516730, "converged"),
    ("job_19", "r2SCAN-3c", "D4", "def2-mTZVP", -402.514410, "converged"),
    ("job_20", "r2SCAN-3c", "D4", "def2-mTZVP", -402.513980, "converged"),
]

_OUT_TEMPLATE = """\
* Info: ORCA calculation from a recorded run tree (fixture, host-neutral)
Program Version: ORCA 5.0.4 (recorded excerpt)
Solvent:          CPCM, water
Functional:       {method}
Dispersion:       {disp} ( Grimme's D4 )
Basis:            {basis}

...
-------------------------
TOTAL ENERGY SUMMARY
-------------------------
Final Single Point Energy ...  {energy:.6f} Eh
-------------------------
* Status: {status}
"""


def build(ws: Path) -> None:
    root = ws / "cm1" / "run_tree"
    if root.exists():
        print(f"refusing to overwrite {root}", file=sys.stderr)
        raise SystemExit(1)
    root.mkdir(parents=True)
    for folder, method, disp, basis, energy, status in _JOBS:
        job = root / folder
        job.mkdir()
        (job / "job.out").write_text(
            _OUT_TEMPLATE.format(method=method, disp=disp, basis=basis,
                                 energy=energy, status=status),
            encoding="utf-8")
        (job / "CONTROL.txt").write_text(
            "# CONTROL.txt (recorded excerpt)\n"
            f"method = {method}/{basis} D4\n"
            "solvent = water\n",
            encoding="utf-8")
    (root / "README.txt").write_text(
        "CM1 benchmark fixture: 20 recorded job folders, host-neutral. "
        "Two jobs failed (SCF, job_04; walltime, job_13); two methods "
        "are present (r2SCAN-3c/def2-mTZVP D4: 15 jobs, "
        "B3LYP/def2-SVP D3BJ: 5 jobs). Never compare total energies "
        "across methods.\n", encoding="utf-8")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    build(ws)
    print(f"prepared {ws / 'cm1' / 'run_tree'} — 20 job folders, "
          "2 failed runs, 2 methods")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
