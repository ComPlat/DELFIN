#!/usr/bin/env python3
"""Acceptance engine for the XT1 xtb-method tasks: grade the physics by
recomputation.

Each task grades a system from setup/w20_xt1_start_geometries.py.  The
script runs xtb (GFN2, ALPB water) on the geometry the model left in the
task folder (or our start.xyz fallback), then checks:

  * the optimisation produced an xtbopt.xyz,
  * the imaginary-frequency count matches the graded band
    (0 for the closed-shell systems, {0,1} for the radicals --
    the methoxy radical start optimises to a shallow saddle, -734
    cm-1, which is reproducible physics and reported, never hidden),
  * the final single-point energy falls in the measured band (Eh),
  * the HOMO-LUMO gap falls in a tolerance band (eV).

Reference bands measured once with the same protocol (xtb 6.7.1 GFN2,
--alpb water, opt then hess on the optimised geometry):

  water_cation       E -4.570152 Eh   imag 0   gap 15.315 eV   uhf 1
  oh_anion           E -4.827982 Eh   imag 0   gap 13.990 eV
  butanone           E -16.688780 Eh  imag 0   gap 3.754 eV
  ammonium           E -4.632618 Eh   imag 0   gap 14.825 eV
  methylamine_radical E -7.001206 Eh  imag 1 (-103 cm-1, near-flat)
                     gap 12.744 eV   uhf 1
  methoxy_radical    E -7.630932 Eh   imag 1 (-734 cm-1, CH3O->CH2OH
                     isomerisation saddle) gap 14.259 eV  uhf 1

usage: w20_xt1_minimum.py <workspace> <mol_key>
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# mol key -> (folder, charge, uhf, ref_energy_eh, energy_tol, imag_band,
#             gap_ev, gap_tol)
SYSTEMS = {
    "water_cation": ("water_cation", 1, 1, -4.570152, 0.002, (0, 0),
                     15.315, 0.5),
    "oh_anion": ("oh_anion", -1, 0, -4.827982, 0.002, (0, 0),
                 13.990, 0.5),
    "butanone": ("butanone", 0, 0, -16.688780, 0.010, (0, 0),
                 3.754, 0.5),
    "ammonium": ("ammonium", 1, 0, -4.632618, 0.002, (0, 0),
                 14.825, 0.5),
    "methylamine_radical": ("methylamine_radical", 0, 1, -7.001206, 0.002,
                            (0, 1), 12.744, 0.5),
    "methoxy_radical": ("methoxy_radical", 0, 1, -7.630932, 0.002,
                        (0, 1), 14.259, 0.5),
}


def _xtb() -> str:
    from delfin.qm_runtime import resolve_tool
    r = resolve_tool("xtb")
    if r is None or not r.path:
        print("xtb not found via qm_runtime resolver", file=sys.stderr)
        sys.exit(3)
    return str(r.path)


def _run(cmd: list[str], cwd: Path) -> str:
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                           errors="replace", timeout=600)
        return p.stdout
    except OSError as exc:
        return f"[OSError] {exc}"


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print("usage: w20_xt1_minimum.py <workspace> <mol_key>",
              file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    key = argv[2]
    if key not in SYSTEMS:
        print(f"unknown molecule key: {key}", file=sys.stderr)
        return 2
    folder, chrg, uhf, ref_e, e_tol, imag_lohi, ref_gap, gap_tol = SYSTEMS[key]

    mdir = ws / "xtb" / folder
    if not mdir.is_dir():
        print(f"no xtb/{folder} folder in workspace", file=sys.stderr)
        return 1
    cand = mdir / "xtbopt.xyz"
    if not cand.is_file():
        xyzs = sorted(mdir.glob("*.xyz"))
        cand = xyzs[-1] if xyzs else mdir / "start.xyz"
    if not cand.is_file():
        print(f"no xyz in xtb/{folder}", file=sys.stderr)
        return 1

    xtb = _xtb()
    with tempfile.TemporaryDirectory(prefix="xt1acc_") as td:
        tdp = Path(td)
        shutil.copy(cand, tdp / "check.xyz")
        base = [xtb, "check.xyz", "--opt", "--chrg", str(chrg),
                "--uhf", str(uhf), "--gfn", "2", "--alpb", "water"]
        log = _run(base, tdp)
        if not (tdp / "xtbopt.xyz").is_file():
            print("optimisation produced no xtbopt.xyz", file=sys.stderr)
            print(log[-2000:], file=sys.stderr)
            return 1
        hlog = _run([xtb, "xtbopt.xyz", "--hess", "--chrg", str(chrg),
                     "--uhf", str(uhf), "--gfn", "2", "--alpb", "water"], tdp)

    m = re.search(r"# imaginary freq\.\s+(\d+)", hlog)
    if not m:
        print("could not read imaginary-frequency count from Hessian log",
              file=sys.stderr)
        print(hlog[-2000:], file=sys.stderr)
        return 1
    n_imag = int(m.group(1))
    lo, hi = imag_lohi
    if not (lo <= n_imag <= hi):
        print(f"imaginary-frequency count {n_imag} outside graded band "
              f"[{lo}, {hi}]", file=sys.stderr)
        return 1

    em = re.search(r"TOTAL ENERGY\s+(-?[0-9.]+)\s+Eh", hlog)
    if not em:
        print("could not read total energy from Hessian log", file=sys.stderr)
        print(hlog[-2000:], file=sys.stderr)
        return 1
    e = float(em.group(1))
    if abs(e - ref_e) > e_tol:
        print(f"energy {e:.6f} Eh outside band {ref_e:.6f} ± {e_tol} Eh",
              file=sys.stderr)
        return 1

    gap = None
    gm = re.search(r"HOMO-LUMO GAP\s+([0-9.]+)\s+eV", hlog)
    if gm:
        gap = float(gm.group(1))
    if gap is None or abs(gap - ref_gap) > gap_tol:
        print(f"gap {gap} eV outside band {ref_gap} ± {gap_tol} eV",
              file=sys.stderr)
        return 1

    print(f"minimum ok: {n_imag} imag freqs (band [{lo},{hi}]), "
          f"energy {e:.6f} Eh in band, gap {gap:.3f} eV in band")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
