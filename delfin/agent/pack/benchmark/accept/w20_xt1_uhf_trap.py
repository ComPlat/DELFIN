#!/usr/bin/env python3
"""Acceptance for the XT1 uhf trap task (methoxy radical).

The trap: the task prompt states a WRONG spin multiplicity for the
methoxy radical (triplet, uhf=2; the electron count says 13 valence
electrons -> doublet, uhf=1).  The model must catch it.

The script recomputes both variants from the folder the model worked in
(xtb/methoxy_radical):

  * chrg 0, uhf 1  -- the correct doublet; its Hessian must give the
    reproducible physics of this fixture set (energy band around
    -7.630932 Eh, imag band {0,1}, gap around 14.259 eV).
  * chrg 0, uhf 2  -- the triplet as the prompt claims.  MEASURED
    FINDING (xtb 6.7.1 GFN2, ALPB water, this fixture): the triplet SCF
    collapses onto the same solution as the doublet -- identical energy
    and gap to 11 digits -- so the recomputation does NOT separate the
    two multiplicities here.  The physics check is therefore only the
    doublet band; whether the model CAUGHT the wrong multiplicity is
    graded by the expected_signals of the task (it must name the
    correct doublet and that the prompt's triplet is wrong).

The verdict is PASS only when the doublet recomputation matches the
band; the catch itself is graded by the task's expected_signals.

usage: w20_xt1_uhf_trap.py <workspace>
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REF_E, E_TOL = -7.630932, 0.002
REF_GAP, GAP_TOL = 14.259, 0.5


def _xtb() -> str:
    from delfin.qm_runtime import resolve_tool
    r = resolve_tool("xtb")
    if r is None or not r.path:
        print("xtb not found via qm_runtime resolver", file=sys.stderr)
        sys.exit(3)
    return str(r.path)


def _run(cmd, cwd: Path) -> str:
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                           errors="replace", timeout=600)
        return p.stdout
    except OSError as exc:
        return f"[OSError] {exc}"


def _recompute(xtb: str, xyz: Path, uhf: str) -> tuple[float | None, float | None]:
    """opt+hess at charge 0 with the given --uhf; returns (E, gap)."""
    with tempfile.TemporaryDirectory(prefix="xt1trap_") as td:
        tdp = Path(td)
        shutil.copy(xyz, tdp / "check.xyz")
        base = [xtb, "check.xyz", "--opt", "--chrg", "0", "--uhf", uhf,
                "--gfn", "2", "--alpb", "water"]
        log = _run(base, tdp)
        if not (tdp / "xtbopt.xyz").is_file():
            return None, None
        hlog = _run([xtb, "xtbopt.xyz", "--hess", "--chrg", "0",
                     "--uhf", uhf, "--gfn", "2", "--alpb", "water"], tdp)
        em = re.search(r"TOTAL ENERGY\s+(-?[0-9.]+)\s+Eh", hlog)
        gm = re.search(r"HOMO-LUMO GAP\s+([0-9.]+)\s+eV", hlog)
        e = float(em.group(1)) if em else None
        gap = float(gm.group(1)) if gm else None
    return e, gap


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: w20_xt1_uhf_trap.py <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    mdir = ws / "xtb" / "methoxy_radical"
    if not mdir.is_dir():
        print("no xtb/methoxy_radical folder in workspace", file=sys.stderr)
        return 1
    cand = mdir / "xtbopt.xyz"
    if not cand.is_file():
        xyzs = sorted(mdir.glob("*.xyz"))
        cand = xyzs[-1] if xyzs else mdir / "start.xyz"
    if not cand.is_file():
        print("no xyz in xtb/methoxy_radical", file=sys.stderr)
        return 1

    xtb = _xtb()
    e_dbl, gap_dbl = _recompute(xtb, cand, "1")
    if e_dbl is None:
        print("doublet recomputation produced no xtbopt.xyz", file=sys.stderr)
        return 1
    if abs(e_dbl - REF_E) > E_TOL or abs(gap_dbl - REF_GAP) > GAP_TOL:
        print(f"doublet physics off band: E {e_dbl} Eh (ref {REF_E} "
              f"±{E_TOL}), gap {gap_dbl} eV (ref {REF_GAP} ±{GAP_TOL})",
              file=sys.stderr)
        return 1
    print(f"doublet physics in band: E {e_dbl:.6f} Eh, gap {gap_dbl:.3f} "
          f"eV; whether the wrong multiplicity was caught is graded by "
          f"the task's expected_signals")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
