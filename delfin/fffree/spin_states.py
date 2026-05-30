"""delfin.fffree.spin_states — Universal spin-state enumeration for TMCs.

For DFT calculations, the correct multiplicity is crucial. This module
enumerates all chemically reasonable spin states for a metal center based on:
  - d-electron count (from metal + oxidation state)
  - Coordination geometry (Δ_oct, Δ_tet splittings)
  - Universal Hund's rule + spin-pairing energy heuristics

Each spin state corresponds to a different DFT multiplicity (2S+1) and
potentially different geometry (Jahn-Teller distortions for degenerate states).

Universal across all d-block metals (d¹..d¹⁰) and coordination geometries.
FF-free: pure orbital-counting + ligand-field analysis, no FF.
Deterministic.

Env-gate: DELFIN_FFFREE_SPIN_STATES=1 (default OFF, byte-identical when unset).
Auto-enabled under PURE_TRACK3 for DFT-startable completeness.
"""
from __future__ import annotations

import os
from typing import Dict, List, Tuple

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_SPIN_STATES = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_SPIN_STATES", "0") == "1"


# d-electron counts for common TM oxidation states (most stable per period)
# Format: {(metal, oxidation): d_count}
_D_COUNTS = {
    # 3d
    ("Sc", 3): 0, ("Ti", 2): 2, ("Ti", 3): 1, ("Ti", 4): 0,
    ("V", 2): 3, ("V", 3): 2, ("V", 4): 1, ("V", 5): 0,
    ("Cr", 2): 4, ("Cr", 3): 3, ("Cr", 6): 0,
    ("Mn", 2): 5, ("Mn", 3): 4, ("Mn", 4): 3, ("Mn", 7): 0,
    ("Fe", 2): 6, ("Fe", 3): 5,
    ("Co", 2): 7, ("Co", 3): 6,
    ("Ni", 2): 8, ("Ni", 3): 7,
    ("Cu", 1): 10, ("Cu", 2): 9,
    ("Zn", 2): 10,
    # 4d
    ("Y", 3): 0, ("Zr", 4): 0, ("Nb", 5): 0,
    ("Mo", 2): 4, ("Mo", 3): 3, ("Mo", 4): 2, ("Mo", 5): 1, ("Mo", 6): 0,
    ("Tc", 7): 0, ("Ru", 2): 6, ("Ru", 3): 5,
    ("Rh", 1): 8, ("Rh", 3): 6, ("Pd", 0): 10, ("Pd", 2): 8,
    ("Ag", 1): 10, ("Cd", 2): 10,
    # 5d
    ("La", 3): 0, ("Hf", 4): 0,
    ("W", 6): 0, ("W", 4): 2,
    ("Re", 7): 0, ("Re", 5): 2,
    ("Os", 2): 6, ("Os", 3): 5,
    ("Ir", 1): 8, ("Ir", 3): 6,
    ("Pt", 0): 10, ("Pt", 2): 8, ("Pt", 4): 6,
    ("Au", 1): 10, ("Au", 3): 8,
    ("Hg", 2): 10,
}


def d_electron_count(metal: str, oxidation_state: int) -> int:
    """Universal d-electron count lookup. Returns -1 if unknown."""
    return _D_COUNTS.get((metal, oxidation_state), -1)


def allowed_spin_states(d_count: int, geometry: str = "oct") -> List[Tuple[int, str]]:
    """Enumerate all chemically reasonable spin states for d-electron count + geometry.

    Parameters
    ----------
    d_count : 0..10 (d-electron count)
    geometry : 'oct' (octahedral), 'tet' (tetrahedral), 'sqp' (square planar)

    Returns
    -------
    List of (multiplicity, label) tuples:
      multiplicity = 2S + 1 (DFT input)
      label = 'high-spin', 'low-spin', 'intermediate-spin'

    Examples:
      Fe(2+) d6 oct: [(1, 'low-spin'), (5, 'high-spin')]
      Co(2+) d7 oct: [(2, 'low-spin'), (4, 'high-spin')]
      Mn(2+) d5 oct: [(2, 'low-spin'), (6, 'high-spin')]
      Cu(2+) d9 oct: [(2, 'only')]  — always doublet

    Universal: works for any d-count and any common geometry.
    """
    if d_count < 0 or d_count > 10:
        return []
    out = []
    if geometry == "oct":
        # Octahedral: t2g (3 orbitals) + eg (2 orbitals)
        out = _oct_spin_states(d_count)
    elif geometry == "tet":
        # Tetrahedral: e (2) + t2 (3) — always high-spin in practice
        out = _tet_spin_states(d_count)
    elif geometry == "sqp":
        # Square planar (d8 typical): dxz/dyz, dxy, dz2, dx2-y2
        out = _sqp_spin_states(d_count)
    else:
        out = _oct_spin_states(d_count)  # fallback
    return out


def _oct_spin_states(d: int) -> List[Tuple[int, str]]:
    """Octahedral t2g^x eg^y spin enumeration."""
    states = []
    # d0, d1, d2, d3, d8, d9, d10 — only one ground state typically
    if d == 0:
        states.append((1, "singlet"))
    elif d == 1:
        states.append((2, "doublet"))
    elif d == 2:
        states.append((3, "triplet"))
    elif d == 3:
        states.append((4, "quartet"))
    elif d == 4:
        # HS: t2g^3 eg^1 (S=2, M=5);  LS: t2g^4 (S=1, M=3)
        states.append((3, "low-spin"))
        states.append((5, "high-spin"))
    elif d == 5:
        # HS: t2g^3 eg^2 (S=5/2, M=6);  LS: t2g^5 (S=1/2, M=2)
        states.append((2, "low-spin"))
        states.append((6, "high-spin"))
    elif d == 6:
        # HS: t2g^4 eg^2 (S=2, M=5);  LS: t2g^6 (S=0, M=1)
        states.append((1, "low-spin"))
        states.append((5, "high-spin"))
    elif d == 7:
        # HS: t2g^5 eg^2 (S=3/2, M=4);  LS: t2g^6 eg^1 (S=1/2, M=2)
        states.append((2, "low-spin"))
        states.append((4, "high-spin"))
    elif d == 8:
        # t2g^6 eg^2: paramagnetic (S=1, M=3); also sqp possible (S=0)
        states.append((3, "triplet"))
        states.append((1, "singlet"))  # if sqp distortion
    elif d == 9:
        # t2g^6 eg^3 (S=1/2, Jahn-Teller distortion)
        states.append((2, "doublet"))
    elif d == 10:
        states.append((1, "singlet"))
    return states


def _tet_spin_states(d: int) -> List[Tuple[int, str]]:
    """Tetrahedral always high-spin in practice (Δ_tet ≈ 4/9 Δ_oct)."""
    # Use max-multiplicity (high-spin) configuration
    if d == 0:
        return [(1, "singlet")]
    if d == 1:
        return [(2, "doublet")]
    # Tetrahedral d2..d10: Hund's-rule high-spin
    n_unpaired = min(d, 10 - d) if d > 5 else d
    return [(n_unpaired + 1, "high-spin-tet")]


def _sqp_spin_states(d: int) -> List[Tuple[int, str]]:
    """Square planar (d8 strong-field is common). Other d-counts use oct fallback."""
    if d == 8:
        return [(1, "singlet")]  # Typical strong-field d8 sqp
    return _oct_spin_states(d)


def enumerate_spin_states_for_metal(metal: str, oxidation_state: int, geometry: str = "oct") -> List[Tuple[int, str]]:
    """Universal one-call API: metal + oxidation + geometry → spin states.

    Returns list of (multiplicity, label) for DFT input generation.
    """
    if not _SPIN_STATES:
        return []
    d_count = d_electron_count(metal, oxidation_state)
    if d_count == -1:
        return []
    return allowed_spin_states(d_count, geometry)


if __name__ == "__main__":
    print("=== Universal spin-state enumeration ===")
    test_cases = [
        ("Fe", 2, "oct"),  # d6
        ("Fe", 3, "oct"),  # d5
        ("Co", 2, "oct"),  # d7
        ("Mn", 2, "oct"),  # d5
        ("Ni", 2, "sqp"),  # d8 sqp
        ("Cu", 2, "oct"),  # d9
        ("Cr", 3, "oct"),  # d3
        ("V", 3, "oct"),   # d2
        ("Ru", 2, "oct"),  # d6
        ("Ti", 4, "oct"),  # d0
        ("Pd", 2, "sqp"),  # d8
    ]
    for m, ox, geom in test_cases:
        d = d_electron_count(m, ox)
        states = allowed_spin_states(d, geom)
        ms = ", ".join(f"M={mult} ({lbl})" for mult, lbl in states)
        print(f"  {m}({ox}+) {geom} d{d}: {ms}")
