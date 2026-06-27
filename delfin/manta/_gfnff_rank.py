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

# Selectable ranking Hamiltonian (DELFIN_CONF_RANK_METHOD).  All run through the SAME
# xtb CLI / energy parse, so switching is just the flag.  GFN-FF is the fast default
# (force field, ~physical conformer energies on metals); GFN2-xTB is the accurate
# semiempirical tight-binding (better ranking, ~10-50x slower); GFN1/GFN0 in between.
_METHOD_FLAGS = {
    "gfnff": ["--gfnff"],
    "gfn2": ["--gfn", "2"],
    "gfn1": ["--gfn", "1"],
    "gfn0": ["--gfn", "0"],
}


def _resolve_method() -> str:
    # Default GFN2-xTB: GFN-FF is too crude to RANK transition-metal-complex
    # conformers — measured on ABAKOE (W(V)) its energy minimum is GFN2's
    # near-maximum (482 kcal/mol), so a GFN-FF ranking keeps garbage.  GFN2 gives a
    # physical low-energy ladder.  Set DELFIN_CONF_RANK_METHOD=gfnff for a fast
    # (cruder) draft pass.
    m = os.environ.get("DELFIN_CONF_RANK_METHOD", "gfn2").strip().lower()
    return m if m in _METHOD_FLAGS else "gfn2"

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
_CACHE: Dict[Tuple[str, int, str], Optional[float]] = {}


def available() -> bool:
    """True if an xtb binary was found (GFN-FF ranking is possible)."""
    return _XTB is not None


def _natoms(xyz_block: str) -> int:
    return sum(1 for ln in xyz_block.splitlines() if len(ln.split()) == 4)


def gfnff_energy(xyz_block: str, charge: int = 0, uhf: int = 0,
                 timeout: float = 120.0, method: Optional[str] = None) -> Optional[float]:
    """Return the total energy of ``xyz_block`` in kcal/mol under the selected xtb
    Hamiltonian (or None on any failure).  ``xyz_block`` is a header-less
    ``Sym x y z`` block (the canonical DELFIN format).  ``method`` overrides the
    DELFIN_CONF_RANK_METHOD env (gfnff | gfn2 | gfn1 | gfn0).  Cached by
    (coordinate-hash, charge, method)."""
    if _XTB is None:
        return None
    meth = (method or _resolve_method())
    if meth not in _METHOD_FLAGS:
        meth = "gfnff"
    key = (hashlib.sha256(xyz_block.encode()).hexdigest(), int(charge), meth)
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
            cmd = [_XTB, fp] + _METHOD_FLAGS[meth] + ["--sp",
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


def gfnff_optimize(xyz_block: str, charge: int = 0, uhf: int = 0,
                   method: Optional[str] = None,
                   timeout: float = 600.0):
    """Geometry OPTIMIZATION via ``xtb --opt`` (default GFN2).

    Unlike :func:`gfnff_energy` (single-point ranking), this RELAXES the
    geometry — used to polish the top-ranked structures into best-possible
    coordinates.  Returns ``(optimized_xyz_block, energy_kcal_or_None)`` or
    ``None`` on any failure (caller then keeps the unrelaxed structure).
    Single-threaded (OMP=1) so each optimization is reproducible.
    """
    if _XTB is None:
        return None
    meth = (method or _resolve_method())
    if meth not in _METHOD_FLAGS:
        meth = "gfn2"
    try:
        na = _natoms(xyz_block)
        if na < 2:
            return None
        with tempfile.TemporaryDirectory() as td:
            fp = os.path.join(td, "conf.xyz")
            with open(fp, "w") as fh:
                fh.write(f"{na}\n\n{xyz_block}\n")
            cmd = [_XTB, fp] + _METHOD_FLAGS[meth] + ["--opt",
                   "--chrg", str(int(charge)), "--uhf", str(int(uhf))]
            res = subprocess.run(cmd, capture_output=True, text=True,
                                 timeout=timeout, cwd=td,
                                 env={**os.environ, "OMP_NUM_THREADS": "1"})
            opt_fp = os.path.join(td, "xtbopt.xyz")
            if not os.path.exists(opt_fp):
                return None
            with open(opt_fp) as fh:
                lines = fh.read().splitlines()
            if len(lines) < 3:
                return None
            try:
                n = int(lines[0].strip())
            except ValueError:
                return None
            opt_xyz = "\n".join(lines[2:2 + n])
            # energy: xtbopt.xyz comment line is e.g. " energy: -42.1 gnorm: ..."
            energy = None
            toks = lines[1].split()
            for i, tok in enumerate(toks):
                if tok.lower().startswith("energy") and i + 1 < len(toks):
                    try:
                        energy = float(toks[i + 1]) * _HARTREE_TO_KCAL
                    except ValueError:
                        energy = None
                    break
            if energy is None:
                for line in res.stdout.splitlines():
                    if "TOTAL ENERGY" in line:
                        parts = line.split()
                        for i, t in enumerate(parts):
                            if t == "ENERGY" and i + 1 < len(parts):
                                try:
                                    energy = float(parts[i + 1]) * _HARTREE_TO_KCAL
                                except ValueError:
                                    energy = None
                                break
                        break
            return (opt_xyz, energy)
    except Exception:
        return None


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
