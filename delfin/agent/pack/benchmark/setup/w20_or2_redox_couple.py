"""Setup for the OR2 redox task (w20-or2-04).

Writes optimised structures of a small organic redox couple into the
task workspace: p-phenylenediamine-like fragment reduced to 18 atoms by
methyl-free ring (p-phenylenediamine, C6H8N2 = 16 atoms; the radical
cation is the 17th electron — same atoms, charge +1).

MAX_ATOMS documents the small-systems limit this package is graded on.
"""

from __future__ import annotations

import pathlib
import sys

MAX_ATOMS = 16


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
    """Write neutral + oxidised structures under <workspace>/or2/redox/."""
    d = workspace / "or2" / "redox"
    d.mkdir(parents=True, exist_ok=True)

    # Planar p-phenylenediamine, ring on xy, N at para positions.
    frame = [
        ("C", 1.210, 0.697, 0.0), ("C", 1.210, -0.697, 0.0),
        ("C", 0.000, 1.395, 0.0), ("C", 0.000, -1.395, 0.0),
        ("C", -1.210, 0.697, 0.0), ("C", -1.210, -0.697, 0.0),
        ("N", 0.000, 2.565, 0.0), ("N", 0.000, -2.565, 0.0),
    ]
    hydrogens = [
        ("H", 2.140, 1.235, 0.0), ("H", 2.140, -1.235, 0.0),
        ("H", -2.140, 1.235, 0.0), ("H", -2.140, -1.235, 0.0),
        ("H", 0.870, 3.130, 0.0), ("H", -0.870, 3.130, 0.0),
        ("H", 0.870, -3.130, 0.0), ("H", -0.870, -3.130, 0.0),
    ]
    atoms = frame + hydrogens  # 16 atoms
    body = "".join(f"{s} {x: .6f} {y: .6f} {z: .6f}\n"
                   for s, x, y, z in atoms)
    (d / "neutral.xyz").write_text(
        f"16\np-phenylenediamine neutral, charge 0\n{body}", encoding="utf-8")
    (d / "oxidised.xyz").write_text(
        f"16\np-phenylenediamine radical cation, charge +1\n{body}",
        encoding="utf-8")

    return {"folder": str(d), "files": ["neutral.xyz", "oxidised.xyz"]}
