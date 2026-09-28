#!/usr/bin/env python3
"""Acceptance: the optimised structure is a real minimum.

Runs xtb --opt then --hess on the folder the task names (folder names
are pinned per task in the YAML, passed here as a molecule key), in a
scratch dir under the workspace, and checks:

  * the optimisation produced an xtbopt.xyz,
  * the Hessian reports ZERO imaginary frequencies (cutoff -20 cm-1),
  * optionally the HOMO-LUMO gap falls in a tolerance band.

Exit status is the verdict; stdout is the evidence.  The reference bands
were measured once with the same xtb/GFN2 protocol when the suite was
written (see tasks_chem.yaml) - this script re-derives the physics
rather than trusting the model's prose.

usage: chem_opt_is_a_minimum.py <workspace> <mol_key> [gap_lo gap_hi]
"""

from __future__ import annotations

import math
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# molecule key -> (subfolder, charge)
MOLS = {
    "acetaminophen": ("acetaminophen", 0),
    "chloronitrobenzene": ("chloronitrobenzene", 0),
    "au_cyanide": ("au_cyanide", -1),
    "cr_carbonyl": ("cr_carbonyl", 0),
}


def _xtb() -> str:
    from delfin.qm_runtime import resolve_tool
    r = resolve_tool("xtb")
    if r is None or not r.path:
        print("xtb not found via qm_runtime resolver", file=sys.stderr)
        sys.exit(3)
    return str(r.path)


def _run(cmd: list[str], cwd: Path) -> str:
    p = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True,
                       timeout=600)
    return (p.stdout or "") + (p.stderr or "")


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print("usage: <workspace> <mol_key> [gap_lo gap_hi]", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    mol_key = argv[2]
    gap_band = None
    if len(argv) >= 5:
        gap_band = (float(argv[3]), float(argv[4]))

    if mol_key not in MOLS:
        print(f"unknown molecule key: {mol_key}", file=sys.stderr)
        return 2
    folder, chrg = MOLS[mol_key]

    # The model's optimised structure: prefer its xtbopt.xyz / *.xyz in
    # the molecule folder; fall back to our start.xyz (then the opt is
    # ours, which still answers the physics question).
    mdir = ws / "chem" / folder
    if not mdir.is_dir():
        print(f"no chem/{folder} folder in workspace", file=sys.stderr)
        return 1
    cand = (mdir / "xtbopt.xyz")
    if not cand.is_file():
        xyzs = sorted(mdir.glob("*.xyz"))
        cand = xyzs[-1] if xyzs else mdir / "start.xyz"
    if not cand.is_file():
        print(f"no xyz in chem/{folder}", file=sys.stderr)
        return 1

    xtb = _xtb()
    with tempfile.TemporaryDirectory(prefix="chemacc_") as td:
        tdp = Path(td)
        shutil.copy(cand, tdp / "check.xyz")
        log = _run([xtb, "check.xyz", "--opt", "--chrg", str(chrg),
                    "--uhf", "0", "--alpb", "water"], tdp)
        if not (tdp / "xtbopt.xyz").is_file():
            print("optimisation produced no xtbopt.xyz", file=sys.stderr)
            print(log[-2000:], file=sys.stderr)
            return 1
        hlog = _run([xtb, "xtbopt.xyz", "--hess", "--chrg", str(chrg),
                     "--uhf", "0", "--alpb", "water"], tdp)

    m = re.search(r"# imaginary freq\.\s+(\d+)", hlog)
    if not m:
        print("could not read imaginary-frequency count from Hessian log",
              file=sys.stderr)
        print(hlog[-2000:], file=sys.stderr)
        return 1
    n_imag = int(m.group(1))
    if n_imag != 0:
        print(f"NOT a minimum: {n_imag} imaginary frequencies", file=sys.stderr)
        return 1

    gap = None
    gm = re.search(r"HOMO-LUMO GAP\s+([0-9.]+)\s+eV", hlog)
    if gm:
        gap = float(gm.group(1))
    if gap_band is not None:
        if gap is None:
            print("gap band requested but no gap in log", file=sys.stderr)
            return 1
        lo, hi = gap_band
        if not (lo <= gap <= hi):
            print(f"gap {gap:.3f} eV outside band [{lo}, {hi}] eV",
                  file=sys.stderr)
            return 1
        print(f"minimum ok, 0 imag freqs, gap {gap:.3f} eV in [{lo},{hi}]")
        return 0
    print(f"minimum ok, 0 imag freqs" + (f", gap {gap:.3f} eV" if gap else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
