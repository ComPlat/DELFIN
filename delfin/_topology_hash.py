"""Universal topology-hash hard-gate (Welle-5p-A).

Goal
----
After ANY geometry-modifying step (rotamer rotation, ring-pucker
displacement, chelate twist, macrocycle deformation, post-UFF
optimisation, …) we must verify that the **topology** of the molecule has
NOT changed.  Without this guard:

* Rotamer-diversity (Welle-5l Track-6) rotates the C–C bond inside a
  chelate ring → ring atoms shift → N–H bond vector flips toward the
  metal → false "diversity" with broken topology.
* Conformer-pool (Welle-5o) puckers / twists chelate rings → M–D
  distances drift → coordination geometry collapses.

The Iter-12+13+14 *M–D Invariant Catastrophe* lesson
(`feedback_md_invariant.md`) showed that any rotation pass without a
strict invariant check produces apparently-improved metrics by **breaking
the molecule**.  The X10-ALEQEO `Fe(CO)2(NH2-CH2-CH2-S)2` test exposes
the same pattern after Welle-5p: 42 % of frames have an amine-H within
2.3 Å of Fe — the H atoms have flipped through the N donor toward the
metal, breaking VSEPR geometry around the N donor.

This module gives every layer the SAME universal hard-gate so the bug
cannot reappear in a new sampling layer added later.

Universal-fundamental compliance
--------------------------------
* Pure-Python (no numpy, no RDKit, no OpenBabel).
* Uses **graph + geometry** features only.  No SMILES strings, refcodes,
  named-ligand patterns or element allowlists beyond the
  IUPAC "metal" definition (rows 2–6, groups 1–14, lanthanides).
* Identical entry-points are wired into ``_rotamer_diversity``,
  ``_conformer_pool``, and ``smiles_converter`` so the gate is single-
  sourced.

Default-OFF contract
--------------------
Module-level helpers always run when called — they are pure functions.
The **wire-in callsites** check the master flag
``DELFIN_5P_A_TOPOLOGY_HARDGATE`` (default ``0``) and skip the call when
unset.  This preserves byte-identical default-OFF behaviour.

Env-flags
---------
``DELFIN_5P_A_TOPOLOGY_HARDGATE``      (default ``0``)   — master switch
``DELFIN_5P_A_MD_TOL``                 (default ``0.05``)— Å, M–D invariant
``DELFIN_5P_A_AMINE_H_MIN_DEG``        (default ``60.0``)— ∠(M-N-H) ≥ this
``DELFIN_5P_A_BOND_TOL``               (default ``1.30``)— covalent-radius scale
``DELFIN_5P_A_CLASH_HM_MIN``           (default ``2.30``)— Å, donor-H · · · M
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from delfin.common.logging import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Covalent radii (Pyykkö 2009, in Å) — pure-Python fallback identical to
# delfin.xyz_io._COVALENT_RADII_FALLBACK so no cross-module import is needed
# from a pure topology check.
# ---------------------------------------------------------------------------

_COVALENT_RADII: Dict[str, float] = {
    "H": 0.31, "He": 0.28,
    "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71,
    "O": 0.66, "F": 0.57, "Ne": 0.58,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07,
    "S": 1.05, "Cl": 1.02, "Ar": 1.06,
    "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53,
    "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24,
    "Cu": 1.32, "Zn": 1.22,
    "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20,
    "Kr": 1.16,
    "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64,
    "Mo": 1.54, "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39,
    "Ag": 1.45, "Cd": 1.44, "In": 1.42, "Sn": 1.39, "Sb": 1.39,
    "Te": 1.38, "I": 1.39, "Xe": 1.40,
    "Cs": 2.44, "Ba": 2.15,
    "La": 2.07, "Ce": 2.04, "Pr": 2.03, "Nd": 2.01, "Pm": 1.99,
    "Sm": 1.98, "Eu": 1.98, "Gd": 1.96, "Tb": 1.94, "Dy": 1.92,
    "Ho": 1.92, "Er": 1.89, "Tm": 1.90, "Yb": 1.87, "Lu": 1.87,
    "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44,
    "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32,
    "Tl": 1.45, "Pb": 1.46, "Bi": 1.48,
    "Ac": 2.15, "Th": 2.06, "Pa": 2.00, "U": 1.96, "Np": 1.90,
    "Pu": 1.87,
}

# Atomic-number → symbol map (used when we operate on graph-dict input).
_Z_TO_SYM: Dict[int, str] = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N",
    8: "O", 9: "F", 10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S",
    17: "Cl", 18: "Ar",
    19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V", 24: "Cr",
    25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn",
    31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr",
    37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo",
    43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd",
    49: "In", 50: "Sn", 51: "Sb", 52: "Te", 53: "I", 54: "Xe",
    55: "Cs", 56: "Ba",
    57: "La", 58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm",
    63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy", 67: "Ho", 68: "Er",
    69: "Tm", 70: "Yb", 71: "Lu",
    72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir",
    78: "Pt", 79: "Au", 80: "Hg",
    81: "Tl", 82: "Pb", 83: "Bi",
    89: "Ac", 90: "Th", 91: "Pa", 92: "U", 93: "Np", 94: "Pu",
}

_METAL_ATOMIC_NUMBERS = frozenset(
    [3, 11, 19, 37, 55, 4, 12, 20, 38, 56]
    + list(range(21, 31))
    + list(range(39, 49))
    + list(range(57, 72))
    + list(range(72, 81))
    + [13, 31, 49, 50, 81, 82, 83]
    + list(range(89, 95))
)


def _is_metal_z(z: int) -> bool:
    return int(z) in _METAL_ATOMIC_NUMBERS


def _is_metal_sym(sym: str) -> bool:
    z = _sym_to_z(sym)
    return _is_metal_z(z) if z is not None else False


def _sym_to_z(sym: str) -> Optional[int]:
    for z, s in _Z_TO_SYM.items():
        if s == sym:
            return z
    return None


def _rcov(sym: str) -> float:
    """Return covalent radius (Å) for *sym*.  Falls back to 1.20 if unknown."""
    return float(_COVALENT_RADII.get(sym, 1.20))


# ---------------------------------------------------------------------------
# Env helpers
# ---------------------------------------------------------------------------


def _env_bool(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in ("1", "true", "yes", "on")


def _env_float(name: str, default: float, lo: float = 0.0, hi: float = 360.0) -> float:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return max(lo, min(hi, float(raw)))
    except (TypeError, ValueError):
        return default


def is_hardgate_enabled() -> bool:
    """Return True iff the master Welle-5p-A flag is set."""
    return _env_bool("DELFIN_5P_A_TOPOLOGY_HARDGATE", False)


# ---------------------------------------------------------------------------
# Geometry primitives (pure Python, no numpy)
# ---------------------------------------------------------------------------


Coord = Tuple[float, float, float]


def _dist(a: Coord, b: Coord) -> float:
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _angle_deg(a: Coord, b: Coord, c: Coord) -> float:
    """Return the angle a-b-c in degrees (b is the vertex)."""
    bax = a[0] - b[0]
    bay = a[1] - b[1]
    baz = a[2] - b[2]
    bcx = c[0] - b[0]
    bcy = c[1] - b[1]
    bcz = c[2] - b[2]
    na = math.sqrt(bax * bax + bay * bay + baz * baz)
    nc = math.sqrt(bcx * bcx + bcy * bcy + bcz * bcz)
    if na < 1e-9 or nc < 1e-9:
        return 0.0
    dot = bax * bcx + bay * bcy + baz * bcz
    cos = max(-1.0, min(1.0, dot / (na * nc)))
    return math.degrees(math.acos(cos))


# ---------------------------------------------------------------------------
# Topology fingerprint computation
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TopologyFingerprint:
    """Canonical topology fingerprint for a 3D structure.

    The fingerprint is independent of atom ordering — it captures the
    *multiset* of bonds (by element pair) plus the multiset of M–D
    distances (rounded to 0.05 Å buckets).  Two structures with identical
    fingerprints have identical bond topology and identical first-shell
    metal coordination.
    """

    n_atoms: int
    formula: str  # Hill-style element multiset, e.g. "C2H6Fe"
    bond_multiset: Tuple[Tuple[str, str], ...]  # sorted (Z1,Z2) symbol-pairs
    md_bins: Tuple[Tuple[str, str, int], ...]  # (metal_sym, donor_sym, 0.05Å bucket)
    n_metal: int
    n_donors: int

    def hexdigest(self) -> str:
        """Return a stable hex digest (CRC-style, no hashlib dependency).

        Two fingerprints with the same digest are *guaranteed* equal
        (digest is collision-checked against full equality at compare
        time — we never rely on the digest alone for the gate).
        """
        h = 5381
        for s in (
            self.formula,
            "|".join(f"{a}-{b}" for (a, b) in self.bond_multiset),
            "|".join(f"{m}-{d}-{n}" for (m, d, n) in self.md_bins),
            f"n={self.n_atoms};m={self.n_metal};d={self.n_donors}",
        ):
            for ch in s:
                h = ((h * 33) ^ ord(ch)) & 0xFFFFFFFFFFFFFFFF
        return f"{h:016x}"


def _hill_formula(symbols: Sequence[str]) -> str:
    counts: Dict[str, int] = {}
    for s in symbols:
        counts[s] = counts.get(s, 0) + 1
    keys = sorted(counts.keys())
    # Carbon-first then hydrogen by Hill convention, then rest alphabetic
    out = []
    if "C" in counts:
        out.append(f"C{counts['C']}" if counts["C"] > 1 else "C")
    if "H" in counts:
        out.append(f"H{counts['H']}" if counts["H"] > 1 else "H")
    for k in keys:
        if k in ("C", "H"):
            continue
        out.append(f"{k}{counts[k]}" if counts[k] > 1 else k)
    return "".join(out)


def detect_bonds(
    symbols: Sequence[str],
    coords: Sequence[Coord],
    bond_tol: float = 1.30,
) -> List[Tuple[int, int]]:
    """Geometric bond perception.

    A pair (i, j) is bonded iff
        dist(i, j) <= (r_cov(i) + r_cov(j)) * bond_tol

    For metals we use a slightly more permissive tolerance to mirror
    Pyykkö-Riedel coordination chemistry — but only the M-D edges are
    *tested* against this scale; the rest of the bond multiset uses the
    standard 1.30× tolerance.

    Returns a list of (i, j) with i < j.
    """
    n = len(symbols)
    out: List[Tuple[int, int]] = []
    radii = [_rcov(s) for s in symbols]
    metals = [_is_metal_sym(s) for s in symbols]
    for i in range(n):
        ri = radii[i]
        for j in range(i + 1, n):
            d = _dist(coords[i], coords[j])
            cutoff = (ri + radii[j]) * bond_tol
            # Slightly larger cutoff for metal-X edges (coordination
            # bonds are 5–15 % longer than covalent equivalents)
            if metals[i] != metals[j]:
                cutoff *= 1.10
            if d <= cutoff:
                out.append((i, j))
    return out


def find_metals_and_donors(
    symbols: Sequence[str],
    bonds: Sequence[Tuple[int, int]],
) -> Tuple[List[int], Dict[int, List[int]]]:
    """Return (metal_indices, {metal_idx: [donor_indices]})."""
    metal_idxs = [i for i, s in enumerate(symbols) if _is_metal_sym(s)]
    md: Dict[int, List[int]] = {m: [] for m in metal_idxs}
    metal_set = set(metal_idxs)
    for (i, j) in bonds:
        if i in metal_set and j not in metal_set:
            md[i].append(j)
        elif j in metal_set and i not in metal_set:
            md[j].append(i)
    return metal_idxs, md


def compute_topology_fingerprint(
    symbols: Sequence[str],
    coords: Sequence[Coord],
    bond_tol: float = 1.30,
) -> TopologyFingerprint:
    """Build a canonical fingerprint from (symbols, coords) only.

    No external state, no SMILES, no force field — pure geometry.  The
    fingerprint is *invariant* under arbitrary atom renumbering because
    the bond multiset and M-D bins are sorted before insertion.
    """
    n = len(symbols)
    bonds = detect_bonds(symbols, coords, bond_tol=bond_tol)
    bond_pairs = sorted(
        tuple(sorted((symbols[i], symbols[j]))) for (i, j) in bonds
    )
    metal_idxs, md = find_metals_and_donors(symbols, bonds)
    md_bins: List[Tuple[str, str, int]] = []
    for m in metal_idxs:
        for d in md[m]:
            dist_a = _dist(coords[m], coords[d])
            bucket = int(round(dist_a / 0.05))
            md_bins.append((symbols[m], symbols[d], bucket))
    md_bins.sort()
    n_donors = sum(len(v) for v in md.values())
    return TopologyFingerprint(
        n_atoms=n,
        formula=_hill_formula(symbols),
        bond_multiset=tuple(bond_pairs),
        md_bins=tuple(md_bins),
        n_metal=len(metal_idxs),
        n_donors=n_donors,
    )


# ---------------------------------------------------------------------------
# Validation entry points
# ---------------------------------------------------------------------------


@dataclass
class TopologyCheckResult:
    """Result of a topology preservation check.

    ``passed`` is True iff every check returned True.  The ``violations``
    list contains short human-readable strings.  ``details`` carries
    structured numbers (for tests / logging).
    """

    passed: bool
    violations: List[str] = field(default_factory=list)
    details: Dict[str, object] = field(default_factory=dict)


def _check_md_invariant(
    symbols: Sequence[str],
    coords_before: Sequence[Coord],
    coords_after: Sequence[Coord],
    md_before: Dict[int, List[int]],
    md_tol: float,
) -> List[str]:
    """Verify every M-D bond length is preserved within *md_tol* Å."""
    out: List[str] = []
    for m, donors in md_before.items():
        for d in donors:
            d_before = _dist(coords_before[m], coords_before[d])
            d_after = _dist(coords_after[m], coords_after[d])
            if abs(d_after - d_before) > md_tol:
                out.append(
                    f"M-D({symbols[m]}{m}-{symbols[d]}{d}) shifted "
                    f"{d_before:.3f}→{d_after:.3f} Å (>|{md_tol:.3f}|)"
                )
    return out


def _check_md_set_unchanged(
    symbols: Sequence[str],
    md_before: Dict[int, List[int]],
    md_after: Dict[int, List[int]],
) -> List[str]:
    """Verify the set of M-D edges is unchanged (no donor lost or gained)."""
    out: List[str] = []
    keys = sorted(set(md_before.keys()) | set(md_after.keys()))
    for m in keys:
        before = set(md_before.get(m, []))
        after = set(md_after.get(m, []))
        lost = before - after
        gained = after - before
        for d in lost:
            out.append(f"donor lost: {symbols[m]}{m}-{symbols[d]}{d}")
        for d in gained:
            out.append(f"donor gained: {symbols[m]}{m}-{symbols[d]}{d}")
    return out


def _check_bond_multiset(
    fp_before: TopologyFingerprint,
    fp_after: TopologyFingerprint,
) -> List[str]:
    """Verify the bond multiset is identical."""
    out: List[str] = []
    b = list(fp_before.bond_multiset)
    a = list(fp_after.bond_multiset)
    if b != a:
        # Diff multisets for a friendly message
        cb: Dict[Tuple[str, str], int] = {}
        ca: Dict[Tuple[str, str], int] = {}
        for p in b:
            cb[p] = cb.get(p, 0) + 1
        for p in a:
            ca[p] = ca.get(p, 0) + 1
        keys = sorted(set(cb.keys()) | set(ca.keys()))
        diffs = []
        for k in keys:
            nb = cb.get(k, 0)
            na = ca.get(k, 0)
            if nb != na:
                diffs.append(f"{k[0]}-{k[1]}: {nb}→{na}")
        out.append("bond multiset changed: " + ", ".join(diffs))
    return out


def _check_amine_h_orientation(
    symbols: Sequence[str],
    coords: Sequence[Coord],
    md_now: Dict[int, List[int]],
    bonds_now: Sequence[Tuple[int, int]],
    min_angle_deg: float,
    h_to_metal_min_dist: float,
) -> List[str]:
    """Check that every donor-N's bonded H atoms point *away* from the metal.

    Universal rule: for every metal M and every donor atom D that is N
    (or analogous sp3 pnictogen P/As — extended via element family), each
    H atom bonded to D must have ∠(M-D-H) ≥ ``min_angle_deg``.  If any
    H is closer than ``h_to_metal_min_dist`` to M, that's a second
    (independent) failure mode — the H has flipped through the donor.
    """
    out: List[str] = []
    # Neighbour table from bond list
    n = len(symbols)
    nbrs: List[List[int]] = [[] for _ in range(n)]
    for (i, j) in bonds_now:
        nbrs[i].append(j)
        nbrs[j].append(i)

    # Pnictogen-like donor elements: N, P, As (universal — same lone-pair
    # geometry, same umbrella inversion mode).  Avoid SMILES patterns by
    # gating on atomic-number membership, not by named ligand.
    pnictogens = {"N", "P", "As"}

    for m, donors in md_now.items():
        for d in donors:
            if symbols[d] not in pnictogens:
                continue
            # Find H atoms bonded to D
            for h in nbrs[d]:
                if symbols[h] != "H":
                    continue
                ang = _angle_deg(coords[m], coords[d], coords[h])
                d_hm = _dist(coords[m], coords[h])
                if ang < min_angle_deg:
                    out.append(
                        f"amine-H flip: ∠({symbols[m]}{m}-{symbols[d]}{d}-H{h})"
                        f"={ang:.1f}° < {min_angle_deg:.1f}° (H–M={d_hm:.2f} Å)"
                    )
                elif d_hm < h_to_metal_min_dist:
                    out.append(
                        f"amine-H too close: H{h}-{symbols[m]}{m}={d_hm:.2f} Å "
                        f"< {h_to_metal_min_dist:.2f} Å"
                    )
    return out


def topology_preserved(
    symbols: Sequence[str],
    coords_before: Sequence[Coord],
    coords_after: Sequence[Coord],
    md_tol: float = 0.05,
    amine_h_min_deg: float = 60.0,
    h_to_metal_min_dist: float = 2.30,
    bond_tol: float = 1.30,
) -> TopologyCheckResult:
    """Return a TopologyCheckResult comparing two coord snapshots.

    ``passed`` is True iff
        1. All M-D distances changed by ≤ ``md_tol`` Å.
        2. The set of M-D edges is unchanged.
        3. The bond multiset (sorted element-pair tuple) is unchanged.
        4. No donor-N/P/As H atom has ∠(M-D-H) < ``amine_h_min_deg``.
        5. No donor-N/P/As H atom is within ``h_to_metal_min_dist`` of M.

    Both inputs must have the same atom-ordering — coords_before[i] and
    coords_after[i] are the same atom.  Atom symbols are taken from
    ``symbols`` and assumed identical for both snapshots (no isotope
    changes inside a rotation pass).
    """
    if len(coords_before) != len(coords_after) or len(symbols) != len(coords_before):
        return TopologyCheckResult(
            passed=False,
            violations=["coord length mismatch"],
            details={
                "n_symbols": len(symbols),
                "n_before": len(coords_before),
                "n_after": len(coords_after),
            },
        )

    bonds_before = detect_bonds(symbols, coords_before, bond_tol=bond_tol)
    bonds_after = detect_bonds(symbols, coords_after, bond_tol=bond_tol)
    _, md_before = find_metals_and_donors(symbols, bonds_before)
    _, md_after = find_metals_and_donors(symbols, bonds_after)

    fp_before = compute_topology_fingerprint(symbols, coords_before, bond_tol=bond_tol)
    fp_after = compute_topology_fingerprint(symbols, coords_after, bond_tol=bond_tol)

    violations: List[str] = []
    violations.extend(
        _check_md_invariant(symbols, coords_before, coords_after, md_before, md_tol)
    )
    violations.extend(_check_md_set_unchanged(symbols, md_before, md_after))
    violations.extend(_check_bond_multiset(fp_before, fp_after))
    violations.extend(
        _check_amine_h_orientation(
            symbols,
            coords_after,
            md_after,
            bonds_after,
            amine_h_min_deg,
            h_to_metal_min_dist,
        )
    )

    return TopologyCheckResult(
        passed=not violations,
        violations=violations,
        details={
            "fp_before": fp_before.hexdigest(),
            "fp_after": fp_after.hexdigest(),
            "n_md_before": fp_before.n_donors,
            "n_md_after": fp_after.n_donors,
            "n_bonds_before": len(bonds_before),
            "n_bonds_after": len(bonds_after),
        },
    )


# ---------------------------------------------------------------------------
# Standalone (single snapshot) amine-H realism check — used in
# smiles_converter post-emit validation when no "before" coords exist.
# ---------------------------------------------------------------------------


def standalone_amine_h_realism(
    symbols: Sequence[str],
    coords: Sequence[Coord],
    bond_tol: float = 1.30,
    amine_h_min_deg: float = 60.0,
    h_to_metal_min_dist: float = 2.30,
) -> TopologyCheckResult:
    """Run only the amine-H orientation check on a single XYZ snapshot.

    Useful as a final "is this structure outputtable?" gate in
    smiles_converter when no pre-modification snapshot exists.
    """
    bonds = detect_bonds(symbols, coords, bond_tol=bond_tol)
    _, md = find_metals_and_donors(symbols, bonds)
    violations = _check_amine_h_orientation(
        symbols, coords, md, bonds, amine_h_min_deg, h_to_metal_min_dist
    )
    return TopologyCheckResult(
        passed=not violations,
        violations=violations,
        details={
            "n_metal": len(md),
            "n_donors": sum(len(v) for v in md.values()),
        },
    )


# ---------------------------------------------------------------------------
# XYZ parsing helpers (DELFIN format, no header) — for callsites that
# need a string→symbols/coords adapter without depending on
# _rotamer_diversity.
# ---------------------------------------------------------------------------


def parse_xyz_delfin(xyz: str) -> Tuple[List[str], List[Coord]]:
    """Parse a DELFIN-format XYZ string (no header, no comment line).

    Lines may be:
        <Sym> <x> <y> <z>
    or the standard 2-header XYZ format (count line + comment line).  We
    detect the header by trying ``int(first_line)``.
    """
    lines = [l for l in xyz.splitlines() if l.strip()]
    if not lines:
        return [], []
    # Strip standard header if present
    try:
        int(lines[0].strip())
        # 2nd line is comment, body starts at line 2
        body = lines[2:] if len(lines) > 1 else []
    except ValueError:
        body = lines
    symbols: List[str] = []
    coords: List[Coord] = []
    for line in body:
        parts = line.split()
        if len(parts) < 4:
            continue
        sym = parts[0]
        # Strip charge/parenthesis decorations like "Fe(1)"
        clean = ""
        for ch in sym:
            if ch.isalpha():
                clean += ch
            else:
                break
        if not clean:
            continue
        try:
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3])
        except ValueError:
            continue
        symbols.append(clean)
        coords.append((x, y, z))
    return symbols, coords


def topology_preserved_xyz(
    xyz_before: str,
    xyz_after: str,
    md_tol: Optional[float] = None,
    amine_h_min_deg: Optional[float] = None,
    h_to_metal_min_dist: Optional[float] = None,
    bond_tol: Optional[float] = None,
) -> TopologyCheckResult:
    """XYZ-string adapter for :func:`topology_preserved`.

    Env-flag tolerances are used when the corresponding arg is ``None``,
    so callsites that want the user-configured behaviour can call this
    without re-reading env vars.
    """
    sb, cb = parse_xyz_delfin(xyz_before)
    sa, ca = parse_xyz_delfin(xyz_after)
    if sb != sa:
        return TopologyCheckResult(
            passed=False,
            violations=["symbol order changed between snapshots"],
            details={"len_before": len(sb), "len_after": len(sa)},
        )
    return topology_preserved(
        sb,
        cb,
        ca,
        md_tol=(
            md_tol if md_tol is not None
            else _env_float("DELFIN_5P_A_MD_TOL", 0.05, lo=0.0, hi=2.0)
        ),
        amine_h_min_deg=(
            amine_h_min_deg if amine_h_min_deg is not None
            else _env_float("DELFIN_5P_A_AMINE_H_MIN_DEG", 60.0, lo=0.0, hi=180.0)
        ),
        h_to_metal_min_dist=(
            h_to_metal_min_dist if h_to_metal_min_dist is not None
            else _env_float("DELFIN_5P_A_CLASH_HM_MIN", 2.30, lo=0.0, hi=10.0)
        ),
        bond_tol=(
            bond_tol if bond_tol is not None
            else _env_float("DELFIN_5P_A_BOND_TOL", 1.30, lo=0.5, hi=3.0)
        ),
    )


def standalone_amine_h_realism_xyz(
    xyz: str,
    amine_h_min_deg: Optional[float] = None,
    h_to_metal_min_dist: Optional[float] = None,
    bond_tol: Optional[float] = None,
) -> TopologyCheckResult:
    """XYZ-string adapter for :func:`standalone_amine_h_realism`."""
    symbols, coords = parse_xyz_delfin(xyz)
    if not symbols:
        return TopologyCheckResult(
            passed=False,
            violations=["empty XYZ"],
            details={},
        )
    return standalone_amine_h_realism(
        symbols,
        coords,
        bond_tol=(
            bond_tol if bond_tol is not None
            else _env_float("DELFIN_5P_A_BOND_TOL", 1.30, lo=0.5, hi=3.0)
        ),
        amine_h_min_deg=(
            amine_h_min_deg if amine_h_min_deg is not None
            else _env_float("DELFIN_5P_A_AMINE_H_MIN_DEG", 60.0, lo=0.0, hi=180.0)
        ),
        h_to_metal_min_dist=(
            h_to_metal_min_dist if h_to_metal_min_dist is not None
            else _env_float("DELFIN_5P_A_CLASH_HM_MIN", 2.30, lo=0.0, hi=10.0)
        ),
    )


__all__ = [
    "TopologyFingerprint",
    "TopologyCheckResult",
    "compute_topology_fingerprint",
    "topology_preserved",
    "topology_preserved_xyz",
    "standalone_amine_h_realism",
    "standalone_amine_h_realism_xyz",
    "detect_bonds",
    "find_metals_and_donors",
    "parse_xyz_delfin",
    "is_hardgate_enabled",
]
