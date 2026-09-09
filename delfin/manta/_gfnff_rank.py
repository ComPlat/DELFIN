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
import logging
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

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
    # g-xTB is not one of the family above even though the flag looks like it.
    # It ships as a statically linked xtb of its own, and the ordinary xtb
    # beside it **accepts** ``--gxtb`` and silently runs GFN2 -- measured
    # identical to the last digit.  So it needs its own binary and a check that
    # the binary is really it; see _find_gxtb and _binary_is_really_gxtb.
    "gxtb": ["--gxtb"],
}

#: Methods that need the separate g-xTB build rather than the ordinary xtb.
_NEEDS_GXTB = frozenset({"gxtb"})


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


def _find_gxtb() -> Optional[str]:
    """The xtb build that actually has g-xTB in it.

    Looked for under two names because the two exist in the wild: DELFIN's own
    installer lays it down as ``xtb-gxtb`` beside the ordinary one, while the
    upstream release installs as plain ``gxtb`` (that is what is on this
    machine, in /usr/local/bin).  A resolver that knew only the first name
    reported "not found" on a system where it was installed.
    """
    override = os.environ.get("DELFIN_GXTB_BINARY")
    if override and os.path.exists(override):
        return override
    try:
        from delfin.qm_runtime import get_qm_tools_bin_dir

        candidate = os.path.join(get_qm_tools_bin_dir(), "xtb-gxtb")
        if os.path.exists(candidate) and os.access(candidate, os.X_OK):
            return candidate
    except Exception:                                     # noqa: BLE001
        pass
    for name in ("xtb-gxtb", "gxtb"):
        found = shutil.which(name)
        if found:
            return found
    return None


_XTB = _find_xtb()
_GXTB = _find_gxtb()
_CACHE: Dict[Tuple[str, int, str], Optional[float]] = {}
_GXTB_VERIFIED: Dict[str, bool] = {}

#: Two hydrogens, far enough apart to be a molecule and small enough that both
#: Hamiltonians answer in well under a second.  Used only to tell the two
#: programs apart, never for chemistry.
_PROBE_XYZ = "H 0.000000 0.000000 0.000000\nH 0.000000 0.000000 0.740000\n"


def _binary_is_really_gxtb(binary: str, timeout: float = 60.0) -> bool:
    """Whether ``binary`` runs g-xTB when asked to, or quietly runs GFN2.

    The failure this guards against is silent by construction: an ordinary xtb
    given ``--gxtb`` prints no warning, exits zero, and returns the GFN2 energy.
    A ranking built on that is a GFN2 ranking wearing a g-xTB label, and nothing
    in the output says so.

    So the two are asked the same question and their answers compared.  A real
    g-xTB and a GFN2 disagree on any molecule; the same number twice means one
    program answered twice.
    """
    cached = _GXTB_VERIFIED.get(binary)
    if cached is not None:
        return cached
    answers = []
    for flags in (_METHOD_FLAGS["gxtb"], _METHOD_FLAGS["gfn2"]):
        answers.append(_probe_energy(binary, flags, timeout))
    good = (answers[0] is not None and answers[1] is not None
            and abs(answers[0] - answers[1]) > 1e-8)
    if not good:
        logger.warning(
            "%s does not appear to be a g-xTB build: --gxtb and --gfn 2 return "
            "%s and %s. An ordinary xtb accepts --gxtb and runs GFN2 silently, "
            "so g-xTB ranking is refused rather than mislabelled.",
            binary, answers[0], answers[1])
    _GXTB_VERIFIED[binary] = good
    return good


def _probe_energy(binary: str, flags, timeout: float) -> Optional[float]:
    """One total energy from ``binary`` on the probe molecule, or None."""
    try:
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, "probe.xyz")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("2\n\n" + _PROBE_XYZ)
            done = subprocess.run(
                [binary, path, *flags, "--norestart"],
                cwd=folder, capture_output=True, text=True, timeout=timeout,
                env={**os.environ, "OMP_NUM_THREADS": "1"})
            for line in reversed((done.stdout or "").splitlines()):
                if "TOTAL ENERGY" in line:
                    for token in line.split():
                        try:
                            return float(token)
                        except ValueError:
                            continue
    except Exception:                                     # noqa: BLE001
        return None
    return None


def binary_for(method: Optional[str] = None) -> Optional[str]:
    """Whichever program the chosen ranking method needs, or None.

    Returns None for g-xTB when no verified g-xTB build is present -- the
    caller is then expected to fall back and say so, rather than rank with
    something else under the g-xTB name.
    """
    name = (method or _resolve_method()).strip().lower()
    if name in _NEEDS_GXTB:
        if _GXTB is None:
            return None
        return _GXTB if _binary_is_really_gxtb(_GXTB) else None
    return _XTB


def available(method: Optional[str] = None) -> bool:
    """True if the binary the chosen method needs was found (and verified)."""
    if method is None:
        return _XTB is not None
    return binary_for(method) is not None


def _natoms(xyz_block: str) -> int:
    return sum(1 for ln in xyz_block.splitlines() if len(ln.split()) == 4)


def gfnff_energy(xyz_block: str, charge: int = 0, uhf: int = 0,
                 timeout: float = 120.0, method: Optional[str] = None,
                 solvent: str = "") -> Optional[float]:
    """Return the total energy of ``xyz_block`` in kcal/mol under the selected xtb
    Hamiltonian (or None on any failure).  ``solvent`` adds ALPB implicit
    solvation, so the screen ranks in the same medium the optimisations and
    everything downstream run in; an unknown solvent name makes xtb fail and
    the frame simply gets no energy rather than a gas-phase one.  ``xyz_block`` is a header-less
    ``Sym x y z`` block (the canonical DELFIN format).  ``method`` overrides the
    DELFIN_CONF_RANK_METHOD env (gfnff | gfn2 | gfn1 | gfn0).  Cached by
    (coordinate-hash, charge, method)."""
    meth = (method or _resolve_method())
    if meth not in _METHOD_FLAGS:
        # Do not quietly substitute a Hamiltonian.  GFN-FF used to be the
        # fallback here, and on a transition-metal complex it is not a slightly
        # worse GFN2: measured on ABAKOE (W(V)) its energy minimum sits at
        # GFN2's near-maximum, 482 kcal/mol away.  A ranking produced that way
        # is not a rougher ranking, it is a different one, and nothing in the
        # result says which method produced it.
        logger.error(
            "Unknown ranking method %r; expected one of %s. Refusing to "
            "substitute another Hamiltonian -- this frame gets no energy.",
            meth, ", ".join(sorted(_METHOD_FLAGS)))
        return None
    binary = binary_for(meth)
    if binary is None:
        return None
    solvent_name = str(solvent or "").strip().lower()
    key = (hashlib.sha256(xyz_block.encode()).hexdigest(), int(charge), meth,
           solvent_name)
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
            cmd = [binary, fp] + _METHOD_FLAGS[meth] + ["--sp",
                   "--chrg", str(int(charge)), "--uhf", str(int(uhf))]
            if solvent_name and solvent_name not in ("gas", "none", "vacuum"):
                cmd += ["--alpb", solvent_name]
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
    meth = (method or _resolve_method())
    if meth not in _METHOD_FLAGS:
        meth = "gfn2"
    binary = binary_for(meth)
    if binary is None:
        return None
    try:
        na = _natoms(xyz_block)
        if na < 2:
            return None
        with tempfile.TemporaryDirectory() as td:
            fp = os.path.join(td, "conf.xyz")
            with open(fp, "w") as fh:
                fh.write(f"{na}\n\n{xyz_block}\n")
            cmd = [binary, fp] + _METHOD_FLAGS[meth] + ["--opt",
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


_Z_GFF = {s: i for i, s in enumerate(
    ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu "
     "Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba "
     "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi"
     ).split(), 1)}


def _n_electrons(xyz_block: str, charge: int) -> int:
    tot = 0
    for line in xyz_block.splitlines():
        p = line.split()
        if len(p) >= 4:
            tot += _Z_GFF.get(p[0], 0)
    return tot - int(charge)


def gfnff_optimize_autospin(xyz_block: str, charge: int = 0,
                            method: Optional[str] = None, timeout: float = 600.0):
    """Geometry-opt SCANNING multiplicity (the determined-spin/D3 path): try UHF in
    {parity, parity+2, parity+4} (parity = n_electrons % 2) and return the LOWEST-energy
    (ground-state multiplicity) result.  Essential for OPEN-SHELL transition metals where
    a fixed uhf=0 gives a wrong GFN2 energy -> wrong ranking / wrong global minimum.
    Falls back to the parity attempt if all energies fail.  Returns ``(xyz, energy)`` or None."""
    try:
        par = _n_electrons(xyz_block, charge) % 2
    except Exception:
        par = 0
    best = None
    for u in (par, par + 2, par + 4):
        r = gfnff_optimize(xyz_block, charge=charge, uhf=u, method=method, timeout=timeout)
        if r and r[1] is not None and (best is None or r[1] < best[1]):
            best = r
    if best is None:
        return gfnff_optimize(xyz_block, charge=charge, uhf=par, method=method, timeout=timeout)
    return best


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
