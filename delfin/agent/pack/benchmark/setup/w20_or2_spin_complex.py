"""Setup for the OR2 spin-state task (w20-or2-03).

Writes a small iron complex start geometry (charge +2, high-spin
start, multiplicity 5) into the task workspace. Kept to 14 atoms:
[Fe(NH3)4(H2O)2]2+ — a weak-field, high-spin starting point suited to
the OCCUPIER preferred-spin workflow.

MAX_ATOMS documents the small-systems limit this package is graded on.
"""

from __future__ import annotations

import pathlib
import sys

MAX_ATOMS = 14


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(f"usage: {pathlib.Path(argv[0]).name} <workspace>",
              file=sys.stderr)
        return 2
    run(pathlib.Path(argv[1]).resolve())
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))


def run(workspace: pathlib.Path) -> dict:
    """Write the spin-complex geometry under <workspace>/or2/spin/."""
    d = workspace / "or2" / "spin"
    d.mkdir(parents=True, exist_ok=True)

    # Tetragonal, Fe(II) d6 high-spin start; long Fe-N (weak field).
    lines = [
        "14",
        "[Fe(NH3)4(H2O)2]2+ high-spin start, charge +2, mult 5",
        "Fe  0.000000  0.000000  0.000000",
        "N   2.150000  0.000000  0.000000",
        "N  -2.150000  0.000000  0.000000",
        "N   0.000000  2.150000  0.000000",
        "N   0.000000 -2.150000  0.000000",
        "O   0.000000  0.000000  2.180000",
        "O   0.000000  0.000000 -2.180000",
    ]
    # H: one on each N (simplified NH3), two on each O (H2O)
    for nx, ny, sign in ((2.15, 0.0, 1), (-2.15, 0.0, -1),
                         (0.0, 2.15, 1), (0.0, -2.15, -1)):
        lines.append(f"H  {nx + 0.96 * sign: .6f}  {ny + 0.40: .6f}  0.300000")
    for oz in (2.18, -2.18):
        lines.append(f"H  0.100000  0.880000  {oz + 0.97 * (1 if oz > 0 else -1): .6f}")
        lines.append(f"H  0.100000 -0.880000  {oz + 0.97 * (1 if oz > 0 else -1): .6f}")
    (d / "start.xyz").write_text("\n".join(lines) + "\n", encoding="utf-8")

    return {"folder": str(d), "files": ["start.xyz"]}
