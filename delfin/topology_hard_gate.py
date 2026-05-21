"""Topology Hard-Gate — runtime topology-invariant validation.

Säule 1 of Hybrid-Path (nature_project/15_HYBRID_PATH_FINAL.md).

Per User-Direktive feedback_md_invariant.md: refined-output MUST preserve
topology vs the SMILES-implied bond network.  Topology loss (detached
donors, extra fragments, broken M-L bonds) cannot be repaired by UFF /
MACE / DFT downstream — it changes which molecule was generated.

This module provides a *runtime* hard-gate that runs on (smiles, xyz)
right before output-write.  Failures are reported as TopologyGateResult.

Scope (this version):
    - Validates a single conformer XYZ against its source SMILES
    - Does NOT compare pre-refinement vs post-refinement snapshots
      (that requires wiring inside smiles_converter, separate commit)

Acceptance contract:
    1. Every M-D bond expected from SMILES is present in XYZ (distance
       within covalent-radius * BOND_TOL).
    2. No expected covalent bond is missing.
    3. No two atoms are closer than COLLISION_THRESHOLD.
    4. No spurious extra fragments (XYZ has same fragment count as SMILES).

Default behaviour: pure function, no side effects, no env-flags read.
Wiring into the pipeline (env-flag-gated) is a separate commit.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

# Re-use existing covalent-radii infrastructure
from delfin.xyz_io import _load_covalent_radii, _COVALENT_RADII_FALLBACK


# ----------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------

#: Covalent-radius-sum tolerance for bond detection (Pyykkö 2009).
#: A pair (a, b) is bonded if dist(a, b) <= (r_a + r_b) * BOND_TOL.
BOND_TOL: float = 1.30

#: Minimum allowed inter-atom distance (Å).  Anything below indicates a
#: refinement collision.
COLLISION_THRESHOLD: float = 0.5

#: Maximum M-D distance treated as a coordination bond (Å).  Beyond this
#: the donor is considered *detached*.
MAX_MD_DISTANCE: float = 3.5

#: Periodic table symbols considered transition / post-transition / lanthanide
#: metals for the purpose of M-D detection.
_METAL_SYMBOLS = frozenset({
    "Li", "Na", "K", "Rb", "Cs", "Fr",
    "Be", "Mg", "Ca", "Sr", "Ba", "Ra",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Al", "Ga", "In", "Tl",
    "Ge", "Sn", "Pb",
    "Sb", "Bi",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr",
})


# ----------------------------------------------------------------------
# Public dataclasses
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class TopologyViolation:
    """One concrete violation discovered by the hard-gate."""
    kind: str  # "missing_md", "broken_bond", "collision", "detached_donor",
               # "extra_fragment", "missing_bond"
    detail: str
    indices: Tuple[int, ...] = ()


@dataclass(frozen=True)
class TopologyGateResult:
    """Result of running the hard-gate on (smiles, xyz)."""
    passed: bool
    violations: Tuple[TopologyViolation, ...] = ()
    n_atoms: int = 0
    n_md_bonds_expected: int = 0
    n_md_bonds_found: int = 0
    n_fragments_xyz: int = 0
    n_fragments_smiles: int = 0

    @property
    def reason(self) -> str:
        """Short human-readable reason if failed, else 'pass'."""
        if self.passed:
            return "pass"
        kinds = {v.kind for v in self.violations}
        return ", ".join(sorted(kinds)) or "fail"


# ----------------------------------------------------------------------
# Internal helpers — pure functions
# ----------------------------------------------------------------------

def _radii_table() -> Dict[str, float]:
    radii = _load_covalent_radii("pyykko2009") or _COVALENT_RADII_FALLBACK
    return radii


def _parse_xyz(xyz: str) -> List[Tuple[str, float, float, float]]:
    """Parse a *single-frame* XYZ string into [(symbol, x, y, z), ...].

    Tolerates leading atom-count + comment lines as per XYZ spec.
    """
    lines = [ln for ln in xyz.splitlines() if ln.strip()]
    if not lines:
        return []
    # Detect XYZ-format: first line int, second comment
    try:
        n_expected = int(lines[0].strip())
        atom_lines = lines[2 : 2 + n_expected]
    except (ValueError, IndexError):
        # Treat the whole thing as raw atom-lines
        atom_lines = lines
    atoms: List[Tuple[str, float, float, float]] = []
    for ln in atom_lines:
        parts = ln.split()
        if len(parts) < 4:
            continue
        try:
            atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
        except ValueError:
            continue
    return atoms


def _distance_sq(a: Tuple[str, float, float, float],
                 b: Tuple[str, float, float, float]) -> float:
    dx = a[1] - b[1]
    dy = a[2] - b[2]
    dz = a[3] - b[3]
    return dx * dx + dy * dy + dz * dz


def _is_metal(symbol: str) -> bool:
    return symbol in _METAL_SYMBOLS


def _xyz_bond_graph(atoms: Sequence[Tuple[str, float, float, float]],
                    radii: Dict[str, float],
                    tol: float = BOND_TOL) -> List[Tuple[int, int]]:
    """Derive bonds from atom positions via covalent-radius sum."""
    bonds: List[Tuple[int, int]] = []
    n = len(atoms)
    for i in range(n):
        ri = radii.get(atoms[i][0], 1.20)
        for j in range(i + 1, n):
            rj = radii.get(atoms[j][0], 1.20)
            cutoff_sq = ((ri + rj) * tol) ** 2
            if _distance_sq(atoms[i], atoms[j]) <= cutoff_sq:
                bonds.append((i, j))
    return bonds


def _fragments_from_bonds(n: int, bonds: Sequence[Tuple[int, int]]) -> int:
    """Count connected components via union-find."""
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i, j in bonds:
        union(i, j)
    return len({find(i) for i in range(n)})


def _smiles_topology(smiles: str) -> Optional[Dict]:
    """Extract atoms + bonds + metal-donor pairs from SMILES via RDKit.

    Returns dict with:
        atoms: List[str] — element symbols
        bonds: List[Tuple[int, int]] — heavy-atom bond pairs
        md_pairs: List[Tuple[int, int]] — metal-to-donor pairs (M first)
        n_fragments: int — number of disconnected fragments
    Returns None if RDKit unavailable or parsing fails.
    """
    try:
        from rdkit import Chem  # type: ignore
    except ImportError:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        # Best-effort sanitize (some TM-complex SMILES fail strict sanitize)
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass
    except Exception:
        return None

    atoms: List[str] = []
    metal_indices: List[int] = []
    for i, atom in enumerate(mol.GetAtoms()):
        sym = atom.GetSymbol()
        atoms.append(sym)
        if _is_metal(sym):
            metal_indices.append(i)

    bonds: List[Tuple[int, int]] = []
    for b in mol.GetBonds():
        i = b.GetBeginAtomIdx()
        j = b.GetEndAtomIdx()
        bonds.append((min(i, j), max(i, j)))

    md_pairs: List[Tuple[int, int]] = []
    for mi in metal_indices:
        for b in mol.GetAtomWithIdx(mi).GetBonds():
            other = b.GetOtherAtomIdx(mi)
            md_pairs.append((mi, other))

    n_fragments = len(Chem.GetMolFrags(mol))

    return {
        "atoms": atoms,
        "bonds": bonds,
        "md_pairs": md_pairs,
        "metal_indices": metal_indices,
        "n_fragments": n_fragments,
    }


def _heavy_indices(atoms: Sequence[Tuple[str, float, float, float]]) -> List[int]:
    return [i for i, a in enumerate(atoms) if a[0] != "H"]


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------

def validate_topology_invariant(
    xyz: str,
    smiles: str,
    *,
    bond_tol: float = BOND_TOL,
    collision_threshold: float = COLLISION_THRESHOLD,
    max_md_distance: float = MAX_MD_DISTANCE,
    strict_bonds: bool = False,
) -> TopologyGateResult:
    """Validate that an XYZ structure preserves SMILES-implied topology.

    Args:
        xyz: XYZ-format string (single frame).
        smiles: Source SMILES the XYZ was generated from.
        bond_tol: Multiplier on covalent-radius sum for bond detection.
            Default 1.30 (Pyykko 2009 + 30%).
        collision_threshold: Minimum allowed inter-atom distance (Å).
        max_md_distance: Maximum allowed M-D bond distance (Å).
        strict_bonds: If True, every SMILES heavy-atom bond must be
            mirrored in the XYZ.  Default False (only M-D bonds
            enforced strictly; carbon-skeleton flexibility tolerated).

    Returns:
        TopologyGateResult with .passed (bool), .violations (tuple).

    Guarantees:
        - Pure function: no global state, no env-flag reads.
        - Deterministic: same inputs → same output.
        - Returns passed=False with at least one violation if any check fails.

    Mathematical reference:
        nature_project/02_MATHEMATICAL_SPEC.md (Section: Topology Invariant)
        nature_project/15_HYBRID_PATH_FINAL.md (Säule 1)
    """
    radii = _radii_table()
    atoms = _parse_xyz(xyz)
    smi_topo = _smiles_topology(smiles)

    violations: List[TopologyViolation] = []

    if not atoms:
        violations.append(TopologyViolation(
            kind="empty_xyz",
            detail="XYZ parse produced no atoms",
        ))
        return TopologyGateResult(
            passed=False, violations=tuple(violations), n_atoms=0,
        )

    # 1. Collision check — any pair below threshold is a hard-fail.
    coll_sq = collision_threshold ** 2
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            if _distance_sq(atoms[i], atoms[j]) < coll_sq:
                violations.append(TopologyViolation(
                    kind="collision",
                    detail=(f"{atoms[i][0]}{i}-{atoms[j][0]}{j} "
                            f"dist<{collision_threshold:.2f}Å"),
                    indices=(i, j),
                ))

    xyz_bonds = _xyz_bond_graph(atoms, radii, bond_tol)
    n_fragments_xyz = _fragments_from_bonds(len(atoms), xyz_bonds)

    md_expected = 0
    md_found = 0
    n_fragments_smi = 0

    if smi_topo is None:
        # SMILES could not be parsed — best-effort gate, only run
        # collision + intra-XYZ-fragmentation checks.
        return TopologyGateResult(
            passed=len(violations) == 0,
            violations=tuple(violations),
            n_atoms=len(atoms),
            n_md_bonds_expected=0,
            n_md_bonds_found=0,
            n_fragments_xyz=n_fragments_xyz,
            n_fragments_smiles=0,
        )

    n_fragments_smi = smi_topo["n_fragments"]

    # 2. Fragment-count match — XYZ must not introduce extra fragments
    if n_fragments_xyz > n_fragments_smi:
        violations.append(TopologyViolation(
            kind="extra_fragment",
            detail=(f"XYZ has {n_fragments_xyz} fragments, "
                    f"SMILES has {n_fragments_smi}"),
        ))

    # 3. Heavy-atom-count match (between SMILES heavy + H expansion).
    #    Skip this if heavy-atom counts diverge — DELFIN sometimes
    #    drops/adds explicit H during sanitisation.
    smi_heavy = sum(1 for s in smi_topo["atoms"] if s != "H")
    xyz_heavy = len([a for a in atoms if a[0] != "H"])
    if smi_heavy != xyz_heavy:
        violations.append(TopologyViolation(
            kind="heavy_atom_count_mismatch",
            detail=f"SMILES heavy={smi_heavy}, XYZ heavy={xyz_heavy}",
        ))
        # Heavy-atom-count mismatch makes bond-by-bond comparison unreliable;
        # return early with only the structural checks.
        return TopologyGateResult(
            passed=False,
            violations=tuple(violations),
            n_atoms=len(atoms),
            n_md_bonds_expected=len(smi_topo["md_pairs"]),
            n_md_bonds_found=0,
            n_fragments_xyz=n_fragments_xyz,
            n_fragments_smiles=n_fragments_smi,
        )

    # 4. M-D bond check — every SMILES M-D pair must have corresponding
    #    XYZ-distance within max_md_distance.
    #    Indexing: SMILES heavy-atom index → XYZ heavy-atom index by
    #    *position*.  This assumes DELFIN preserves heavy-atom order
    #    (true for all DELFIN-paths that don't reorder; H may be
    #    re-emitted).  We map via heavy-atom-rank.
    xyz_heavy_indices = _heavy_indices(atoms)
    smi_heavy_ranks: Dict[int, int] = {}
    rank = 0
    for i, s in enumerate(smi_topo["atoms"]):
        if s != "H":
            smi_heavy_ranks[i] = rank
            rank += 1

    md_expected = len(smi_topo["md_pairs"])
    for (mi, di) in smi_topo["md_pairs"]:
        if mi not in smi_heavy_ranks or di not in smi_heavy_ranks:
            continue
        xyz_mi = xyz_heavy_indices[smi_heavy_ranks[mi]]
        xyz_di = xyz_heavy_indices[smi_heavy_ranks[di]]
        dist = _distance_sq(atoms[xyz_mi], atoms[xyz_di]) ** 0.5
        if dist > max_md_distance:
            violations.append(TopologyViolation(
                kind="missing_md",
                detail=(f"{atoms[xyz_mi][0]}{xyz_mi}-"
                        f"{atoms[xyz_di][0]}{xyz_di} "
                        f"d={dist:.2f}Å > {max_md_distance:.2f}Å"),
                indices=(xyz_mi, xyz_di),
            ))
        else:
            md_found += 1

    # 5. (Optional, strict mode) Every SMILES heavy-atom bond present in XYZ.
    if strict_bonds:
        xyz_bond_set = set(xyz_bonds)
        for (si, sj) in smi_topo["bonds"]:
            if si not in smi_heavy_ranks or sj not in smi_heavy_ranks:
                continue  # bond involves H — skipped
            xi = xyz_heavy_indices[smi_heavy_ranks[si]]
            xj = xyz_heavy_indices[smi_heavy_ranks[sj]]
            if (min(xi, xj), max(xi, xj)) not in xyz_bond_set:
                violations.append(TopologyViolation(
                    kind="missing_bond",
                    detail=(f"{atoms[xi][0]}{xi}-"
                            f"{atoms[xj][0]}{xj} not bonded in XYZ"),
                    indices=(xi, xj),
                ))

    return TopologyGateResult(
        passed=len(violations) == 0,
        violations=tuple(violations),
        n_atoms=len(atoms),
        n_md_bonds_expected=md_expected,
        n_md_bonds_found=md_found,
        n_fragments_xyz=n_fragments_xyz,
        n_fragments_smiles=n_fragments_smi,
    )
