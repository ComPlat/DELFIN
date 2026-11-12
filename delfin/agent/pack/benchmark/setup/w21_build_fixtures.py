#!/usr/bin/env python3
"""Fixtures for the wave-21 workflow-BUILDING tasks (package WB0).

Input: the workspace root (argv[1]). Output: w21_build_<k>/ with a start
structure (start.xyz) and brief.txt, one folder per task, following the
convention of w20_rest.py. The tasks ask the agent to BUILD A WORKFLOW
(a DELFIN pipeline spec) for a chemist's request — no calculation is
graded by the acceptance script; workflow_grader.grade() judges the spec
the agent writes into spec.json / answer.txt in the task folder.

Exit 1 without writing anything when any folder already exists.
Exit 2 on a missing argument.
"""
from __future__ import annotations

import sys
from pathlib import Path

# (folder, n_atoms, header, coords, brief) — one entry per task.
_TASKS = [
    ("w21_build_conformers", 9, "C2H6O charge=0 mult=1",
     "C 0.0 0.0 0.0\nC 1.43 0.0 0.0\nO -0.7 1.2 0.0"
     "\nH 0.36 0.0 0.98\nH 0.36 0.0 -0.98\nH 1.79 0.0 -0.94"
     "\nH 1.79 0.0 0.94\nH -0.7 1.9 0.62\nH -0.7 1.9 -0.62",
     "butanol's parent ethanol is flexible. Before anything expensive I "
     "need its conformational landscape mapped cheaply, and then the "
     "thermochemistry of the lowest conformer at one consistent level of "
     "theory."),
    ("w21_build_redox", 5, "CH2O charge=1 mult=2",
     "C 0.0 0.0 0.0\nO 0.0 0.0 1.22\nH 0.94 0.0 -0.61"
     "\nH -0.94 0.0 -0.61",
     "I need a one-electron potential for the formaldehyde cation couple "
     "in acetonitrile. Both members of the couple must be computed the "
     "same way — I will compare them."),
    ("w21_build_ts", 5, "CH3+ H-abstraction coords",
     "C 0.0 0.0 0.0\nH 1.08 0.0 0.0\nH -0.54 0.93 0.0"
     "\nH -0.54 -0.93 0.0\nH 0.0 0.0 1.30",
     "Find the transition state for hydrogen abstraction from methane by "
     "a bare proton, and I need it verified — a saddle point claim alone "
     "is worthless to me. Methane is CH4, the proton comes in along z."),
    ("w21_build_spin", 6, "FeL5 dummy ligand field",
     "Fe 0.0 0.0 0.0\nN 0.0 0.0 2.0\nN 0.0 2.0 0.0\nN 2.0 0.0 0.0"
     "\nN -2.0 0.0 0.0\nN 0.0 -2.0 0.0",
     "A small iron(II) pentamine complex: is the high-spin or the "
     "low-spin state lower? Rank the spin states and tell me the "
     "splitting. Fe is Fe(II), the five N ligands are neutral donors."),
    ("w21_build_cycle", 4, "NH3 in water",
     "N 0.0 0.0 0.12\nH 0.0 0.94 -0.05\nH 0.81 -0.47 -0.05"
     "\nH -0.81 -0.47 -0.05",
     "Give me the free energy of transfer of ammonia from the gas phase "
     "into water, as a thermodynamic cycle: the solvation contribution at "
     "one level of theory plus the gas-phase reference, both members "
     "computed consistently."),
    ("w21_build_spectrum", 5, "CH2O vib/rot IR",
     "C 0.0 0.0 0.0\nO 0.0 0.0 1.22\nH 0.94 0.0 -0.61"
     "\nH -0.94 0.0 -0.61",
     "I want the infrared spectrum of formaldehyde in water: the harmonic "
     "frequencies, their IR intensities, and which bands dominate the "
     "spectrum."),
]


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: w21_build_fixtures.py <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    existing = [name for name, _, _, _, _ in _TASKS if (ws / name).exists()]
    if existing:
        print(f"refusing to overwrite {existing}", file=sys.stderr)
        return 1
    for name, natoms, header, coords, brief in _TASKS:
        d = ws / name
        d.mkdir(parents=True)
        (d / "start.xyz").write_text(
            f"{natoms}\n{header}\n{coords}\n", encoding="utf-8")
        (d / "brief.txt").write_text(brief + "\n", encoding="utf-8")
    print(f"prepared {len(_TASKS)} w21 build fixtures under {ws}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
