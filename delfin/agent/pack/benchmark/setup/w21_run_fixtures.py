#!/usr/bin/env python3
"""Fixtures for the wave-21 real-run workflow tasks (package WB0).

Like w21_build_fixtures.py but for the three tasks that RUN the built
workflow for real (tiny, <= 6 atoms, through DELFIN's submit path).
Each folder carries start.xyz + brief.txt + a CONTROL.txt the agent's
run uses; the CONTROL pins charge/multiplicity so CB0A's
charge_spin_consistent can check the run against it.

Exit 1 without writing anything when any folder already exists.
Exit 2 on a missing argument.
"""
from __future__ import annotations

import sys
from pathlib import Path

# (folder, n_atoms, header, coords, brief, control)
_TASKS = [
    ("w21_run_methanol", 6, "CH4O charge=0 mult=1",
     "C 0.0 0.0 0.0\nO 1.43 0.0 0.0"
     "\nH 0.36 0.0 0.98\nH 0.36 0.0 -0.98\nH 1.79 0.0 -0.94\nH 1.79 0.0 0.94",
     "Optimise methanol in the gas phase and confirm the geometry is a "
     "true minimum; report the thermochemical correction at 298 K.",
     "workflow = classic_opt_freq\ncharge = 0\nmult = 1\n"
     "method = B3LYP\nbasis = def2-SVP\nsolvent = \n"),
    ("w21_run_water_ir", 3, "H2O charge=0 mult=1",
     "O 0.0 0.0 0.0\nH 0.76 0.59 0.0\nH -0.76 0.59 0.0",
     "Compute the harmonic frequencies and IR intensities of water and "
     "report the frequency of the bending mode.",
     "workflow = classic_opt_freq\ncharge = 0\nmult = 1\n"
     "method = B3LYP\nbasis = def2-SVP\nsolvent = \n"),
    ("w21_run_hcl_spectrum", 2, "HCl charge=0 mult=1",
     "H 0.0 0.0 0.0\nCl 1.27 0.0 0.0",
     "Compute the infrared spectrum of hydrogen chloride and report the "
     "stretch frequency.",
     "workflow = classic_opt_freq\ncharge = 0\nmult = 1\n"
     "method = B3LYP\nbasis = def2-SVP\nsolvent = \n"),
]


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: w21_run_fixtures.py <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    existing = [name for name, _, _, _, _, _ in _TASKS
                if (ws / name).exists()]
    if existing:
        print(f"refusing to overwrite {existing}", file=sys.stderr)
        return 1
    for name, natoms, header, coords, brief, control in _TASKS:
        d = ws / name
        d.mkdir(parents=True)
        (d / "start.xyz").write_text(
            f"{natoms}\n{header}\n{coords}\n", encoding="utf-8")
        (d / "brief.txt").write_text(brief + "\n", encoding="utf-8")
        (d / "CONTROL.txt").write_text(control, encoding="utf-8")
    print(f"prepared {len(_TASKS)} w21 run fixtures under {ws}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
