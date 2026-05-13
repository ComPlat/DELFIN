"""Per-element van der Waals radii for DELFIN clash detection.

Used by Baustein 5 (post-optimizer) and related geometry-validation passes
to derive distance thresholds for steric-clash and H-H crowding checks.

References
----------
Bondi, A. (1964). "van der Waals Volumes and Radii."
    J. Phys. Chem. 68 (3): 441-451.
    Primary source for main-group elements (H, C, N, O, F, P, S, Cl, ...).

Alvarez, S. (2013). "A cartography of the van der Waals territories."
    Dalton Trans. 42: 8617-8636.
    Extension for transition metals, lanthanides, and actinides.

Notes
-----
- Metals are normally NOT checked against each other for vdW clashes
  (they act as anchors). For genuine M-M bonded systems, use
  ``_METAL_METAL_BOND_LENGTHS`` from ``delfin.smiles_converter`` instead.
- The default ``H_H_FACTOR`` (0.75) is intentionally more permissive than
  the heavy-atom factor (0.85), because Bondi H-H sum (2.40 A) is known
  to over-call mild proximity in crowded ligand spheres.
"""

from __future__ import annotations

from typing import Dict

# ---------------------------------------------------------------------------
# Primary radii table (Angstrom).
# Main-group: Bondi (1964). Transition metals / Ln / An: Alvarez (2013).
# ---------------------------------------------------------------------------
VDW_RADII: Dict[str, float] = {
    # Main-group (Bondi 1964)
    "H":  1.20, "He": 1.40, "Li": 1.82, "Be": 1.53, "B":  1.92, "C":  1.70,
    "N":  1.55, "O":  1.52, "F":  1.47, "Ne": 1.54, "Na": 2.27, "Mg": 1.73,
    "Al": 1.84, "Si": 2.10, "P":  1.80, "S":  1.80, "Cl": 1.75, "Ar": 1.88,
    "K":  2.75, "Ca": 2.31, "Ga": 1.87, "Ge": 2.11, "As": 1.85, "Se": 1.90,
    "Br": 1.85, "Kr": 2.02, "Rb": 3.03, "Sr": 2.49, "In": 1.93, "Sn": 2.17,
    "Sb": 2.06, "Te": 2.06, "I":  1.98, "Xe": 2.16, "Cs": 3.43, "Ba": 2.68,
    "Tl": 1.96, "Pb": 2.02, "Bi": 2.07, "Po": 1.97, "At": 2.02, "Rn": 2.20,

    # 3d transition metals (Alvarez 2013)
    "Sc": 2.11, "Ti": 1.87, "V":  1.79, "Cr": 1.89, "Mn": 1.97, "Fe": 1.94,
    "Co": 1.92, "Ni": 1.84, "Cu": 1.86, "Zn": 2.10,
    # 4d transition metals
    "Y":  2.19, "Zr": 1.86, "Nb": 2.07, "Mo": 2.09, "Tc": 2.09, "Ru": 2.07,
    "Rh": 1.95, "Pd": 2.02, "Ag": 2.03, "Cd": 2.30,
    # Lanthanides
    "La": 2.43, "Ce": 2.42, "Pr": 2.40, "Nd": 2.39, "Sm": 2.38, "Eu": 2.40,
    "Gd": 2.38, "Tb": 2.37, "Dy": 2.35, "Ho": 2.33, "Er": 2.32, "Tm": 2.30,
    "Yb": 2.40, "Lu": 2.21,
    # 5d transition metals
    "Hf": 2.23, "Ta": 2.22, "W":  2.18, "Re": 2.16, "Os": 2.16, "Ir": 2.13,
    "Pt": 2.13, "Au": 2.14, "Hg": 2.23,
    # Actinides (limited Alvarez coverage)
    "Th": 2.43, "U":  1.86,
}

# Default scaling factors used by clash-threshold helpers.
DEFAULT_FACTOR: float = 0.85          # heavy-atom clash gate
DEFAULT_H_H_FACTOR: float = 0.75      # looser gate for H/H pairs
_COVALENT_FALLBACK_SCALE: float = 1.2  # vdW ~ 1.2 * covalent if unknown
_LAST_RESORT_RADIUS: float = 1.80      # used if even covalent is missing


def _covalent_fallback(symbol: str) -> float:
    """Return ``1.2 * covalent_radius`` for ``symbol``; ``_LAST_RESORT_RADIUS`` if unknown.

    Imported lazily to avoid a module-import cycle with ``smiles_converter``.
    """
    try:
        from delfin.smiles_converter import _COVALENT_RADII  # type: ignore
    except Exception:
        return _LAST_RESORT_RADIUS
    r_cov = _COVALENT_RADII.get(symbol)
    if r_cov is None:
        return _LAST_RESORT_RADIUS
    return float(r_cov) * _COVALENT_FALLBACK_SCALE


def get_vdw_radius(symbol: str) -> float:
    """Return the Bondi/Alvarez vdW radius (Angstrom) for ``symbol``.

    If the element is not tabulated, fall back to ``1.2 * covalent_radius``
    (using ``delfin.smiles_converter._COVALENT_RADII``); if even that is
    missing, return ``1.80`` as a last-resort default.
    """
    if not symbol:
        return _LAST_RESORT_RADIUS
    r = VDW_RADII.get(symbol)
    if r is not None:
        return float(r)
    return _covalent_fallback(symbol)


def get_vdw_sum(sym1: str, sym2: str) -> float:
    """Return the sum of two vdW radii (Angstrom)."""
    return get_vdw_radius(sym1) + get_vdw_radius(sym2)


def get_clash_threshold(
    sym1: str,
    sym2: str,
    factor: float = DEFAULT_FACTOR,
    h_h_factor: float = DEFAULT_H_H_FACTOR,
) -> float:
    """Distance below which ``(sym1, sym2)`` is considered clashing.

    H-H pairs use the more permissive ``h_h_factor`` (default 0.75) instead
    of ``factor`` (default 0.85), because the strict Bondi 2.04 A H-H
    threshold over-calls mild proximity in crowded ligands.
    """
    if sym1 == "H" and sym2 == "H":
        return get_vdw_sum(sym1, sym2) * h_h_factor
    return get_vdw_sum(sym1, sym2) * factor


def is_clash(
    sym1: str,
    sym2: str,
    distance: float,
    factor: float = DEFAULT_FACTOR,
    h_h_factor: float = DEFAULT_H_H_FACTOR,
) -> bool:
    """Return ``True`` iff ``distance`` is below the clash threshold for the pair."""
    return distance < get_clash_threshold(sym1, sym2, factor=factor, h_h_factor=h_h_factor)


if __name__ == "__main__":  # pragma: no cover - quick sanity check
    print("VDW_RADII size:", len(VDW_RADII))
    samples = ["H", "C", "N", "O", "F", "Fe", "Cu", "Pt", "U", "Xx"]
    print("Sample radii (A):")
    for s in samples:
        print(f"  {s:>3s}: {get_vdw_radius(s):.3f}")

    print("\nPair thresholds (factor=0.85, h_h=0.75):")
    pairs = [("H", "H"), ("C", "C"), ("C", "H"), ("N", "O"), ("Cu", "N"), ("Fe", "O")]
    for a, b in pairs:
        thr = get_clash_threshold(a, b)
        print(f"  {a}-{b}: sum={get_vdw_sum(a, b):.3f}  thr={thr:.3f}")

    print("\nis_clash checks:")
    cases = [
        ("H", "H", 2.10),  # 2.10 > 1.80 -> no clash (h_h=0.75 -> 1.80)
        ("H", "H", 1.70),  # below 1.80 -> clash
        ("C", "C", 3.50),  # above 2.89 -> no clash (pi-stack OK)
        ("C", "C", 2.50),  # below 2.89 -> clash
        ("C", "H", 2.30),  # threshold ~ 2.465 -> clash
        ("C", "H", 2.60),  # above ~ 2.465 -> no clash
    ]
    for a, b, d in cases:
        print(f"  {a}-{b} @ {d:.2f} A -> clash={is_clash(a, b, d)}")
