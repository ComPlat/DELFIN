"""delfin.grip_io — XYZ I/O + atom-mapping utilities for the ``delfin-grip`` CLI.

This module is the *user-facing data layer* for the GRIP standalone tool:

* :func:`read_xyz` / :func:`write_xyz` — multi-frame XYZ parsing/emission
* :func:`perceive_bonds_from_xyz` — build a RDKit RWMol with bonds detected
  from 3D coordinates, using ``rdkit.Chem.rdDetermineBonds`` as the primary
  perception engine (with an OpenBabel fallback)
* :func:`match_smiles_to_xyz` — robust mapping from SMILES-atom-idx ->
  XYZ-atom-idx (handles equivalent atoms via substructure-match + a
  Hungarian fallback)
* :func:`detect_metal_and_donors` — heuristic metal/donor detection so the
  CLI works without explicit ``--metal`` / ``--donors`` flags

The module is *pure I/O + perception* — it never calls into GRIP itself.
It is deliberately a thin layer so callers can use it independently of the
CLI (e.g. notebooks, RDKit-plugins).

Design contract
---------------
* Deterministic — no RNG, sorted iteration order, stable element keys.
* Honours ``PYTHONHASHSEED=0`` for the dict-keyed atom maps.
* Robust to user mistakes (bad XYZ format, missing comment line, mixed
  element capitalisation, trailing whitespace).
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

_LOG = logging.getLogger(__name__)

__all__ = [
    "Frame",
    "read_xyz",
    "write_xyz",
    "perceive_bonds_from_xyz",
    "match_smiles_to_xyz",
    "detect_metal_and_donors",
    "TRANSITION_METAL_Z",
    "LANTHANIDE_Z",
    "ACTINIDE_Z",
]

# Z numbers for transition metals + lanthanides + actinides.  Used by the
# metal-detector heuristic.  Set rather than range so membership-check is O(1).
TRANSITION_METAL_Z = frozenset(
    list(range(21, 31))   # Sc..Zn
    + list(range(39, 49))   # Y..Cd
    + list(range(57, 81))   # La..Hg  (covers Hf..Hg incl. La..Lu via lanthanides)
    + list(range(89, 113))  # Ac..Cn  (covers actinides + super-heavy)
)
LANTHANIDE_Z = frozenset(range(57, 72))
ACTINIDE_Z = frozenset(range(89, 104))


# ---------------------------------------------------------------------------
# Frame dataclass — one snapshot of N atoms
# ---------------------------------------------------------------------------
@dataclass
class Frame:
    """One XYZ frame: elements + coordinates + comment line.

    ``coordinates`` is always a ``(N, 3)`` float64 ndarray.  ``elements`` is
    a list of length N with element symbols (canonical capitalisation
    "Pd", "Cl", ...).  ``comment`` is the second line of the XYZ block as
    a stripped string.
    """

    elements: List[str] = field(default_factory=list)
    coordinates: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), dtype=np.float64))
    comment: str = ""

    @property
    def n_atoms(self) -> int:
        return len(self.elements)

    def __post_init__(self):
        self.coordinates = np.asarray(self.coordinates, dtype=np.float64)
        if self.coordinates.ndim == 1 and self.coordinates.size == 0:
            self.coordinates = np.zeros((0, 3), dtype=np.float64)
        if self.coordinates.ndim == 1:
            if self.coordinates.size % 3 != 0:
                raise ValueError(
                    f"Frame.coordinates flat length {self.coordinates.size} not div by 3"
                )
            self.coordinates = self.coordinates.reshape(-1, 3)
        if self.coordinates.shape[0] != len(self.elements):
            raise ValueError(
                f"Frame: elements ({len(self.elements)}) and coordinates "
                f"({self.coordinates.shape[0]}) length mismatch"
            )


# ---------------------------------------------------------------------------
# Element canonicalisation
# ---------------------------------------------------------------------------
def _canonical_element(token: str) -> str:
    """Normalise an element token to canonical capitalisation.

    XYZ files in the wild contain ``H``, ``h``, ``CL``, ``cl``, ``Cl``,
    ``11`` (atomic number), etc.  We accept all of these.
    """
    t = token.strip()
    if not t:
        raise ValueError("empty element token")
    if t.isdigit():
        # Atomic number — convert via RDKit.
        try:
            from rdkit.Chem import GetPeriodicTable
            pt = GetPeriodicTable()
            return pt.GetElementSymbol(int(t))
        except Exception as exc:  # pragma: no cover
            raise ValueError(f"cannot resolve atomic number {t}: {exc!r}")
    # Standard alpha symbol.  Capitalise correctly: 'cl' -> 'Cl', 'NA' -> 'Na'
    return t[0].upper() + t[1:].lower() if len(t) > 1 else t.upper()


# ---------------------------------------------------------------------------
# read_xyz / write_xyz
# ---------------------------------------------------------------------------
def read_xyz(path: Union[str, Path]) -> List[Frame]:
    """Read a (possibly multi-frame) XYZ file into a list of :class:`Frame`.

    Recognises the standard XYZ format:

    .. code-block::

        <N>
        <comment>
        <element>  <x>  <y>  <z>
        ...

    Multi-frame trajectories are detected by repeating ``<N>`` blocks.
    Empty files raise ``ValueError`` with a clear message.  Lines beyond
    ``2 + N`` per frame are interpreted as the next frame's header.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to the XYZ file.

    Returns
    -------
    list of :class:`Frame`
        One Frame per snapshot.  At least one is always returned (or
        ``ValueError`` raised).
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"XYZ file not found: {p}")

    text = p.read_text(encoding="utf-8", errors="replace")
    if not text.strip():
        raise ValueError(f"XYZ file is empty: {p}")
    lines = text.splitlines()

    frames: List[Frame] = []
    i = 0
    while i < len(lines):
        # Skip leading blank lines between frames
        while i < len(lines) and not lines[i].strip():
            i += 1
        if i >= len(lines):
            break

        # Header (atom count)
        try:
            n_atoms = int(lines[i].strip().split()[0])
        except (ValueError, IndexError):
            raise ValueError(
                f"XYZ parse error in {p}: line {i+1!r}: expected atom count, "
                f"got {lines[i]!r}"
            )
        if n_atoms <= 0:
            raise ValueError(
                f"XYZ parse error in {p}: non-positive atom count {n_atoms} "
                f"at line {i+1}"
            )

        # Comment
        comment = lines[i + 1].strip() if i + 1 < len(lines) else ""

        # Atom block
        block_start = i + 2
        block_end = block_start + n_atoms
        if block_end > len(lines):
            raise ValueError(
                f"XYZ parse error in {p}: frame at line {i+1} declares "
                f"{n_atoms} atoms but only {len(lines) - block_start} remain"
            )
        elements: List[str] = []
        coords = np.zeros((n_atoms, 3), dtype=np.float64)
        for j in range(n_atoms):
            row = lines[block_start + j].split()
            if len(row) < 4:
                raise ValueError(
                    f"XYZ parse error in {p}: line {block_start + j + 1}: "
                    f"expected 4 tokens, got {len(row)}"
                )
            elements.append(_canonical_element(row[0]))
            try:
                coords[j] = (float(row[1]), float(row[2]), float(row[3]))
            except ValueError as exc:
                raise ValueError(
                    f"XYZ parse error in {p}: line {block_start + j + 1}: "
                    f"{exc}"
                )

        frames.append(Frame(elements=elements, coordinates=coords, comment=comment))
        i = block_end

    if not frames:
        raise ValueError(f"XYZ file contains no frames: {p}")
    return frames


def write_xyz(
    path: Union[str, Path],
    frames: Union[Frame, Sequence[Frame]],
    comment_prefix: str = "",
) -> None:
    """Write one or more :class:`Frame`s to an XYZ file.

    Parameters
    ----------
    path : str or Path
        Destination path.
    frames : Frame or sequence of Frame
        Frames to emit.  A single Frame is wrapped automatically.
    comment_prefix : str
        Prepended to each frame's comment line (with a trailing space if
        non-empty).  Useful for tagging: ``comment_prefix="GRIP refined"``
        produces ``"GRIP refined original-comment"``.
    """
    if isinstance(frames, Frame):
        frames = [frames]
    frames = list(frames)
    if not frames:
        raise ValueError("write_xyz: no frames to write")

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    lines: List[str] = []
    for frame in frames:
        n = frame.n_atoms
        comment = frame.comment
        if comment_prefix:
            comment = (comment_prefix + " " + comment).strip()
        lines.append(str(n))
        lines.append(comment)
        for elem, (x, y, z) in zip(frame.elements, frame.coordinates):
            lines.append(f"{elem:<3s} {x:15.8f} {y:15.8f} {z:15.8f}")
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Bond perception
# ---------------------------------------------------------------------------
def perceive_bonds_from_xyz(
    elements: Sequence[str],
    coords: np.ndarray,
    charge: int = 0,
    *,
    use_huckel: bool = False,
    embed_chirality: bool = True,
):
    """Build a RDKit ``RWMol`` with bonds + bond orders + hybridisation
    perceived from the 3D coordinates.

    Uses :func:`rdkit.Chem.rdDetermineBonds.DetermineBonds` as the primary
    engine (covers organics + most main-group + many transition metals via
    distance-based connectivity).  Falls back to a distance-only graph
    built from covalent-radius sums when ``DetermineBonds`` raises (typical
    failure mode: f-block + open-shell TM with charge mismatch).

    Parameters
    ----------
    elements : sequence of str
        Element symbols, length N.
    coords : ndarray (N, 3)
        Cartesian coordinates in Ångström.
    charge : int, default 0
        Net molecular charge.  Required by ``DetermineBonds`` for proper
        bond-order assignment.
    use_huckel : bool, default False
        Pass through to ``DetermineBonds`` to enable Hückel aromaticity
        perception (slower; only useful for π-systems).
    embed_chirality : bool, default True
        Run ``Chem.AssignStereochemistryFrom3D`` after bond perception.

    Returns
    -------
    rdkit.Chem.RWMol
        Mol with the same atom order as the input (atom index ``i`` in
        the mol corresponds to ``elements[i]`` / ``coords[i]``).
    """
    from rdkit import Chem
    from rdkit.Chem import RWMol

    coords = np.asarray(coords, dtype=np.float64)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(f"coords must be (N, 3), got shape {coords.shape}")
    if coords.shape[0] != len(elements):
        raise ValueError(
            f"elements ({len(elements)}) and coords ({coords.shape[0]}) length mismatch"
        )

    rwmol = RWMol()
    for sym in elements:
        rwmol.AddAtom(Chem.Atom(sym))
    conf = Chem.Conformer(len(elements))
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, (float(x), float(y), float(z)))
    rwmol.AddConformer(conf, assignId=True)

    mol = rwmol.GetMol()
    try:
        from rdkit.Chem import rdDetermineBonds
        rdDetermineBonds.DetermineBonds(
            mol, charge=int(charge), useHueckel=bool(use_huckel)
        )
    except Exception as exc:
        _LOG.warning(
            "rdDetermineBonds.DetermineBonds failed (%r); falling back to "
            "distance-only connectivity (no bond orders).", exc,
        )
        mol = _distance_only_connectivity(elements, coords)

    if embed_chirality:
        try:
            Chem.AssignStereochemistryFrom3D(mol)
        except Exception:
            pass

    return RWMol(mol)


def _distance_only_connectivity(
    elements: Sequence[str], coords: np.ndarray
):
    """Fallback bond perception using covalent-radius sums.

    Two atoms are connected when their distance is below
    ``(r_cov_i + r_cov_j) * 1.3``.  No bond orders or aromaticity are
    assigned (everything is a SINGLE bond).  Sufficient as a graph for
    GRIP fragment detection, which only cares about connectivity +
    element identity.
    """
    from rdkit import Chem
    from rdkit.Chem import GetPeriodicTable, RWMol

    pt = GetPeriodicTable()
    coords = np.asarray(coords, dtype=np.float64)
    n = len(elements)

    rwmol = RWMol()
    for sym in elements:
        rwmol.AddAtom(Chem.Atom(sym))

    # Per-atom covalent radius (fallback 0.75 for unknown).
    radii: List[float] = []
    for sym in elements:
        try:
            radii.append(float(pt.GetRcovalent(sym)))
        except Exception:
            radii.append(0.75)

    # O(N^2) edges; fine for typical sub-1000-atom user inputs.
    for i in range(n):
        for j in range(i + 1, n):
            dij = float(np.linalg.norm(coords[i] - coords[j]))
            cutoff = 1.3 * (radii[i] + radii[j])
            if dij < cutoff:
                try:
                    rwmol.AddBond(i, j, Chem.BondType.SINGLE)
                except Exception:
                    pass

    conf = Chem.Conformer(n)
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, (float(x), float(y), float(z)))
    rwmol.AddConformer(conf, assignId=True)
    return rwmol.GetMol()


# ---------------------------------------------------------------------------
# SMILES <-> XYZ atom mapping
# ---------------------------------------------------------------------------
def match_smiles_to_xyz(
    smiles: str,
    xyz_mol,
    xyz_coords: Optional[np.ndarray] = None,
    *,
    strict_element: bool = True,
) -> Optional[Dict[int, int]]:
    """Match SMILES atom indices to XYZ atom indices.

    Three strategies (in order):

    1. **Substructure match** — ``Chem.GetSubstructMatch`` on the SMILES
       mol against the XYZ mol.  Handles renaming and most reorderings.
    2. **Hungarian assignment** — when the substructure match returns
       empty or atom-counts disagree (typical: SMILES has implicit H, XYZ
       has explicit H, OR SMILES is a slight isomer of XYZ).  Builds a
       cost matrix from element-type mismatch + 3D-environment-similarity
       and solves with :func:`scipy.optimize.linear_sum_assignment`.
    3. **Identity** — last-resort fallback when both above fail and atom
       counts match: return ``{i: i}`` with an informational warning.

    Parameters
    ----------
    smiles : str
        Input SMILES.  Hydrogens are added explicitly via
        ``Chem.AddHs(MolFromSmiles(smiles))`` for the comparison.
    xyz_mol : rdkit RWMol
        Bond-perceived mol from :func:`perceive_bonds_from_xyz`.
    xyz_coords : ndarray, optional
        Coordinates used by strategy 2.  Defaults to the conformer in
        ``xyz_mol``.
    strict_element : bool, default True
        If True, mappings whose mapped elements disagree are rejected
        (returns ``None``).  Set to False to allow heuristic matches with
        a warning.

    Returns
    -------
    dict[int, int] or None
        ``{smiles_atom_idx: xyz_atom_idx}``.  ``None`` when no acceptable
        mapping was found.
    """
    from rdkit import Chem

    smiles_mol = Chem.MolFromSmiles(smiles)
    if smiles_mol is None:
        _LOG.warning("match_smiles_to_xyz: cannot parse SMILES %r", smiles)
        return None
    smiles_mol = Chem.AddHs(smiles_mol)

    n_smi = smiles_mol.GetNumAtoms()
    n_xyz = xyz_mol.GetNumAtoms()

    # --- Strategy 1: substructure match ------------------------------------
    # Try both directions: SMILES-in-XYZ (typical when SMILES omits some
    # weakly-bound ligands), and XYZ-in-SMILES (typical when XYZ has
    # incomplete connectivity).
    try:
        match = xyz_mol.GetSubstructMatch(smiles_mol)
        if match and len(match) == n_smi:
            mapping = {i: int(match[i]) for i in range(n_smi)}
            if _validate_mapping_elements(smiles_mol, xyz_mol, mapping, strict_element):
                return mapping
    except Exception as exc:
        _LOG.debug("substructure match (SMILES->XYZ) failed: %r", exc)

    try:
        match = smiles_mol.GetSubstructMatch(xyz_mol)
        if match and len(match) == n_xyz:
            # invert mapping
            mapping = {int(match[i]): i for i in range(n_xyz)}
            if _validate_mapping_elements(smiles_mol, xyz_mol, mapping, strict_element):
                return mapping
    except Exception as exc:
        _LOG.debug("substructure match (XYZ->SMILES) failed: %r", exc)

    # --- Strategy 2: Hungarian assignment ----------------------------------
    if n_smi == n_xyz:
        try:
            mapping = _hungarian_mapping(
                smiles_mol, xyz_mol, xyz_coords
            )
            if mapping is not None and _validate_mapping_elements(
                smiles_mol, xyz_mol, mapping, strict_element
            ):
                return mapping
        except Exception as exc:
            _LOG.debug("Hungarian fallback failed: %r", exc)

    # --- Strategy 3: identity (only if counts equal AND elements match) ---
    if n_smi == n_xyz:
        mapping = {i: i for i in range(n_smi)}
        if _validate_mapping_elements(smiles_mol, xyz_mol, mapping, strict_element=True):
            _LOG.info(
                "match_smiles_to_xyz: substructure + Hungarian failed; "
                "falling back to identity mapping (SMILES and XYZ atom order "
                "happen to align)."
            )
            return mapping

    _LOG.warning(
        "match_smiles_to_xyz: no valid mapping (SMILES %d atoms, XYZ %d "
        "atoms)", n_smi, n_xyz,
    )
    return None


def _validate_mapping_elements(
    smiles_mol, xyz_mol, mapping: Dict[int, int], strict_element: bool
) -> bool:
    """Sanity-check that mapped atoms have matching element symbols."""
    mismatches = 0
    for smi_i, xyz_i in mapping.items():
        try:
            s_sym = smiles_mol.GetAtomWithIdx(int(smi_i)).GetSymbol()
            x_sym = xyz_mol.GetAtomWithIdx(int(xyz_i)).GetSymbol()
        except Exception:
            return False
        if s_sym != x_sym:
            mismatches += 1
    if mismatches == 0:
        return True
    if strict_element:
        return False
    _LOG.warning(
        "match_smiles_to_xyz: %d element mismatches in mapping (proceeding "
        "anyway, strict_element=False)", mismatches,
    )
    return True


def _hungarian_mapping(
    smiles_mol, xyz_mol, xyz_coords: Optional[np.ndarray]
) -> Optional[Dict[int, int]]:
    """Hungarian-assignment fallback.

    Cost = (element-mismatch penalty) + (3D-neighborhood-distance penalty).
    Same-element atoms get a 0-penalty; cross-element atoms a 100-penalty
    (effectively forbidden).  The 3D distance term uses a heuristic
    embedding from RDKit's ``EmbedMolecule`` for the SMILES side.
    """
    try:
        from scipy.optimize import linear_sum_assignment
        from rdkit.Chem import AllChem
    except Exception:
        return None

    n_smi = smiles_mol.GetNumAtoms()
    n_xyz = xyz_mol.GetNumAtoms()
    if n_smi != n_xyz:
        return None

    syms_smi = [smiles_mol.GetAtomWithIdx(i).GetSymbol() for i in range(n_smi)]
    syms_xyz = [xyz_mol.GetAtomWithIdx(i).GetSymbol() for i in range(n_xyz)]

    if xyz_coords is None:
        try:
            conf = xyz_mol.GetConformer()
            xyz_coords = np.array(
                [list(conf.GetAtomPosition(i)) for i in range(n_xyz)],
                dtype=np.float64,
            )
        except Exception:
            xyz_coords = None

    # Embed the SMILES side for a structural cost (deterministic seed).
    try:
        smi_copy = type(smiles_mol)(smiles_mol)  # shallow clone
        AllChem.EmbedMolecule(smi_copy, randomSeed=42)
        conf_s = smi_copy.GetConformer()
        smi_coords = np.array(
            [list(conf_s.GetAtomPosition(i)) for i in range(n_smi)],
            dtype=np.float64,
        )
    except Exception:
        smi_coords = None

    # Build cost matrix.
    cost = np.full((n_smi, n_xyz), 100.0, dtype=np.float64)
    for i in range(n_smi):
        for j in range(n_xyz):
            if syms_smi[i] != syms_xyz[j]:
                continue
            base = 0.0
            if smi_coords is not None and xyz_coords is not None:
                # use environment fingerprint: sorted distances to other atoms
                d_s = np.sort(np.linalg.norm(smi_coords - smi_coords[i], axis=1))
                d_x = np.sort(np.linalg.norm(xyz_coords - xyz_coords[j], axis=1))
                base = float(np.linalg.norm(d_s - d_x))
            cost[i, j] = base

    try:
        row_ind, col_ind = linear_sum_assignment(cost)
    except Exception:
        return None

    mapping = {int(r): int(c) for r, c in zip(row_ind, col_ind)}
    # Reject if any row got assigned to a forbidden column.
    for r, c in mapping.items():
        if cost[r, c] >= 99.0:
            return None
    return mapping


# ---------------------------------------------------------------------------
# Metal + donor heuristic detection
# ---------------------------------------------------------------------------
def detect_metal_and_donors(
    mol,
    coords: Optional[np.ndarray] = None,
    *,
    distance_multiplier: float = 1.3,
) -> Tuple[Optional[int], List[int]]:
    """Heuristically identify the metal atom and its donor atoms.

    Strategy:

    1. ``metal_idx`` = the first atom whose Z is in
       :data:`TRANSITION_METAL_Z` (or lanthanide/actinide).  If multiple
       metals exist, returns the FIRST one (most common single-metal case);
       multi-metal complexes need the explicit ``--metal`` flag.
    2. ``donors`` = atoms within ``distance_multiplier * (r_cov_metal +
       r_cov_atom)`` of the metal.  Hydrogens are skipped (donors are
       typically heavy atoms).
    3. When the perceived mol has explicit metal-bonded edges, those are
       preferred over the distance heuristic (covers ``rdDetermineBonds``
       successful cases).

    Parameters
    ----------
    mol : RDKit Mol-like
        Bond-perceived mol.
    coords : ndarray (N, 3), optional
        Atomic coordinates.  Falls back to the mol's conformer.
    distance_multiplier : float
        Coordination cutoff multiplier on the covalent-radius sum.

    Returns
    -------
    (metal_idx, donors) : (int or None, list of int)
        Sorted list of donor indices.  Returns ``(None, [])`` when no
        transition-metal atom is present (organic-only input).
    """
    from rdkit.Chem import GetPeriodicTable

    pt = GetPeriodicTable()
    n = mol.GetNumAtoms()

    if coords is None:
        try:
            conf = mol.GetConformer()
            coords = np.array(
                [list(conf.GetAtomPosition(i)) for i in range(n)],
                dtype=np.float64,
            )
        except Exception:
            coords = None

    # --- 1. Locate the metal atom (first transition-metal Z encountered) ---
    metal_idx: Optional[int] = None
    for i in range(n):
        try:
            z = int(mol.GetAtomWithIdx(i).GetAtomicNum())
        except Exception:
            continue
        if z in TRANSITION_METAL_Z:
            metal_idx = i
            break
    if metal_idx is None:
        return None, []

    # --- 2. Donors from explicit metal-bonded edges ------------------------
    donors: List[int] = []
    try:
        m_atom = mol.GetAtomWithIdx(metal_idx)
        for nbr in m_atom.GetNeighbors():
            ni = int(nbr.GetIdx())
            if nbr.GetSymbol() == "H":
                continue
            donors.append(ni)
    except Exception:
        donors = []

    # --- 3. Fallback: distance-based -------------------------------------
    if not donors and coords is not None:
        try:
            m_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
            r_m = float(pt.GetRcovalent(m_sym))
        except Exception:
            r_m = 1.40
        for i in range(n):
            if i == metal_idx:
                continue
            try:
                sym_i = mol.GetAtomWithIdx(i).GetSymbol()
            except Exception:
                continue
            if sym_i == "H":
                continue
            try:
                r_i = float(pt.GetRcovalent(sym_i))
            except Exception:
                r_i = 0.75
            d = float(np.linalg.norm(coords[i] - coords[metal_idx]))
            if d < distance_multiplier * (r_m + r_i):
                donors.append(i)

    return int(metal_idx), sorted(set(donors))
