#!/usr/bin/env python3
"""Fixture for the fundamentals sign/unit task: two energies, one factor.

Gives E(A) and E(B) in hartree plus the conversion factor, and states
the difference convention ONCE so the task measures conversion, sign
and interpretation — not constant recall. Every number is re-checked
against the arithmetic before the fixture is accepted:

    dE = E(B) - E(A) = -0.0208 Eh = -13.052 kcal/mol
    (negative because B is LOWER: B is the more stable isomer)

Exit non-zero on any failure: a task whose precondition was not built
is reported as unmeasured, never as a model failure.
"""
from __future__ import annotations

import sys
from pathlib import Path

_EH_TO_KCAL = 627.5094740631
_EA = -154.772300
_EB = -154.793100

_ENERGIES = f"""two isomers, total electronic energies (PBE0/def2-SVP, gas phase)
isomer A: {_EA:.6f} Eh
isomer B: {_EB:.6f} Eh
conversion: 1 Eh = 627.509 kcal/mol
difference convention: dE = E(B) - E(A)
"""


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: fund_sign_unit.py <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    d = ws / "fund_sign_unit"
    if d.exists():
        print(f"refusing to overwrite {d}", file=sys.stderr)
        return 1
    # The precondition, checked BEFORE anything is written: the numbers
    # the rubric grades must be what the arithmetic gives.
    d_eh = _EB - _EA
    kcal = d_eh * _EH_TO_KCAL
    if abs(d_eh - (-0.0208)) > 1e-9 or abs(kcal - (-13.052)) > 5e-4:
        print(f"fixture arithmetic drifted: dE={d_eh:.6f} Eh = "
              f"{kcal:.4f} kcal/mol", file=sys.stderr)
        return 1
    d.mkdir(parents=True)
    (d / "energies.txt").write_text(_ENERGIES, encoding="utf-8")
    print(f"prepared {d} — energies.txt (dE = E(B)-E(A) = {kcal:.3f} "
          f"kcal/mol; B lower)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
