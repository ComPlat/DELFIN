"""Setup for the OR2 formaldehyde tasks (w20-or2-01/02).

Writes one start geometry of formaldehyde (H2CO, 4 atoms, charge 0)
plus a pre-optimised copy for the TD-DFT task, into the task's own
workspace folder. Nothing outside it is touched.

MAX_ATOMS documents the small-systems limit this package is graded on.
"""

from __future__ import annotations

import pathlib

MAX_ATOMS = 4


def run(workspace: pathlib.Path) -> dict:
    """Write formaldehyde geometries under <workspace>/or2/formaldehyde/."""
    d = workspace / "or2" / "formaldehyde"
    d.mkdir(parents=True, exist_ok=True)

    # C at origin, O on +x, H2 symmetric on +y/z-distorted C2v plane
    start = d / "start.xyz"
    start.write_text(
        "4\n"
        "formaldehyde start geometry (C2v, charge 0)\n"
        "C  0.000000  0.000000  0.000000\n"
        "O  1.210000  0.000000  0.000000\n"
        "H -0.650000  0.930000  0.000000\n"
        "H -0.650000 -0.930000  0.000000\n",
        encoding="utf-8")

    opt = d / "formaldehyde_opt.xyz"
    opt.write_text(
        "4\n"
        "formaldehyde pre-optimised geometry (B3LYP/def2-SVP quality placeholder)\n"
        "C  0.000000  0.000000  0.000000\n"
        "O  1.208500  0.000000  0.000000\n"
        "H -0.654500  0.935800  0.000000\n"
        "H -0.654500 -0.935800  0.000000\n",
        encoding="utf-8")

    return {"folder": str(d), "files": ["start.xyz",
                                        "formaldehyde_opt.xyz"]}
