"""GFN-FF energy evaluation for conformer ranking (license-clean, metal-capable).

UFF (via OpenBabel) gives unusable energies on transition-metal complexes — measured
spreads of 170-300 kcal/mol and frequent numeric blow-ups (1e16) where the physical
conformer window is a few kcal/mol.  GFN-FF (Grimme's generic force field, xtb) is
parametrised for the whole periodic table including metals and returns physical
relative conformer energies (measured 8 kcal/mol spread on the same pool UFF put at
301).  This module shells out to the ``xtb`` CLI (``--gfnff --sp``) and parses the
total energy; results are cached by coordinate hash + charge so a repeated frame is
free.  Everything here is OPT-IN (the caller gates on DELFIN_CONF_GFNFF_RANK); with
the flag off this module is never imported.

Charge note: for ranking the conformers OF ONE isomer the total charge is constant,
so even an approximate charge ranks them self-consistently; the caller may pass the
true net charge when it has it.
"""
from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple

_HARTREE_TO_KCAL = 627.5094740631

# resolve the xtb binary once (PATH, then common conda/micromamba locations)
def _find_xtb() -> Optional[str]:
    cand = os.environ.get("DELFIN_XTB_PATH")
    if cand and os.path.exists(cand):
        return cand
    w = shutil.which("xtb")
    if w:
        return w
    for p in (
        os.path.expanduser("~/micromamba/bin/xtb"),
        "/opt/xtb/bin/xtb",
        os.path.expanduser("~/miniconda3/bin/xtb"),
        os.path.expanduser("~/anaconda3/bin/xtb"),
    ):
        if os.path.exists(p):
            return p
    return None


_XTB = _find_xtb()
_CACHE: Dict[Tuple[str, int], Optional[float]] = {}


def available() -> bool:
    """True if an xtb binary was found (GFN-FF ranking is possible)."""
    return _XTB is not None


def _natoms(xyz_block: str) -> int:
    return sum(1 for ln in xyz_block.splitlines() if len(ln.split()) == 4)


def gfnff_energy(xyz_block: str, charge: int = 0,
                 uhf: int = 0, timeout: float = 120.0) -> Optional[float]:
    """Return the GFN-FF total energy of ``xyz_block`` in kcal/mol (or None on any
    failure).  ``xyz_block`` is a header-less ``Sym x y z`` block (the canonical
    DELFIN format).  Cached by (coordinate-hash, charge)."""
    if _XTB is None:
        return None
    key = (hashlib.sha256(xyz_block.encode()).hexdigest(), int(charge))
    if key in _CACHE:
        return _CACHE[key]
    val: Optional[float] = None
    try:
        na = _natoms(xyz_block)
        if na < 2:
            _CACHE[key] = None
            return None
        with tempfile.TemporaryDirectory() as td:
            fp = os.path.join(td, "conf.xyz")
            with open(fp, "w") as fh:
                fh.write(f"{na}\n\n{xyz_block}\n")
            cmd = [_XTB, fp, "--gfnff", "--sp",
                   "--chrg", str(int(charge)), "--uhf", str(int(uhf))]
            res = subprocess.run(cmd, capture_output=True, text=True,
                                 timeout=timeout, cwd=td,
                                 env={**os.environ, "OMP_NUM_THREADS": "1"})
            for line in res.stdout.splitlines():
                # "          | TOTAL ENERGY              -42.123456789 Eh   |"
                if "TOTAL ENERGY" in line:
                    parts = line.split()
                    for i, tok in enumerate(parts):
                        if tok == "ENERGY" and i + 1 < len(parts):
                            try:
                                val = float(parts[i + 1]) * _HARTREE_TO_KCAL
                            except ValueError:
                                val = None
                            break
                    break
    except Exception:
        val = None
    _CACHE[key] = val
    return val


def rerank(scored: List[Tuple], top_m: int, charge: int = 0,
           uhf: int = 0) -> Optional[List[Tuple]]:
    """Re-rank a UFF-scored candidate list with GFN-FF.

    ``scored`` is the Welle-5o ``[(uff_energy, cand_xyz, tag, coords), ...]`` list,
    pre-sorted-ish by UFF.  Take the ``top_m`` lowest-UFF candidates (the cheap
    shortlist), recompute their energy with GFN-FF, and return the SAME tuple shape
    with the GFN-FF energy substituted, sorted ascending.  Returns None if GFN-FF is
    unavailable or every shortlist evaluation failed (caller then keeps UFF order)."""
    if _XTB is None or not scored:
        return None
    shortlist = sorted(scored, key=lambda t: t[0])[:max(1, top_m)]
    out: List[Tuple] = []
    for tup in shortlist:
        cand_xyz = tup[1]
        e = gfnff_energy(cand_xyz, charge=charge, uhf=uhf)
        if e is None:
            continue
        out.append((e,) + tuple(tup[1:]))
    if not out:
        return None
    out.sort(key=lambda t: t[0])
    return out
