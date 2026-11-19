"""Heavy-atom RMSD between poses of the same molecule, without alignment.

Input: two coordinate sets of the same molecule (a pose and its reference)
plus the atom mapping that pairs their heavy atoms.  Output: the RMSD over
the heavy atoms, in Angstrom, in the coordinates as they stand.

Semantics: the docking frame is the frame.  Re-aligning the two coordinate
sets before the RMSD answers "how differently do they sit in space" and
gives a lower number for every pose a rigid docking program places, so the
pose quality a redocking check asks about would be hidden — no
superposition happens here.

The mapping is where symmetry lives: a benzene ring, a methyl group, two
equal substituents on one phenyl — swapping their equal atoms changes the
coordinates without changing the molecule, so the reported RMSD is the
smallest one over the atom pairings the molecule's own symmetry allows.
The pairings come from RDKit's canonical atom ranking
(``Chem.CanonicalRankAtoms`` with ``breakTies=False``), which gives equal
ranks exactly to the atoms the molecule itself cannot tell apart; within
each equal-rank class the pairing that minimises the RMSD is chosen
(classes are independent, so this is exact).  Hydrogens are excluded:
their positions are placed by the preparation, not by the pose, and
including them would put a number into a question about the heavy-atom
geometry.

A caller that already has its own atom mapping (a docking program that
writes the same atom order in every pose, for instance) can pass it as
``pairs`` and skip the symmetry search.

Rejected alternatives: taking the RMSD over one fixed atom order (a pose
written with its ring atoms in the other order reads as far away while
sitting on top of the reference); aligning first (hides pose differences,
see above); including hydrogens (measures the hydrogen placement, not the
pose).
"""

from __future__ import annotations

import math
from itertools import permutations
from typing import Dict, List, Optional, Sequence, Tuple

from delfin.tools.bio._errors import StructureFormatError
from delfin.tools.bio._model import SmallMolecule

_RDKIT = None


def _rdkit():
    global _RDKIT
    if _RDKIT is None:
        try:
            import rdkit as module
            import rdkit.Chem  # noqa: F401 — pulls the Chem submodule in
        except ImportError as exc:
            raise StructureFormatError(
                "pose RMSD needs RDKit, which is not importable "
                f"({exc}); install the 'rdkit' package (BSD-3 licence)"
            ) from None
        _RDKIT = module
    return _RDKIT


def _is_same_molecule(a: SmallMolecule, b: SmallMolecule) -> bool:
    """The same molecule, judged on the heavy-atom graph.

    Both sides go through RDKit as a bond graph with element symbols; the
    comparison is on canonical SMILES, so the atom order the two files
    wrote is not a factor.
    """
    return _canonical_smiles(a) == _canonical_smiles(b)


def _canonical_smiles(molecule: SmallMolecule) -> str:
    rdkit = _rdkit()
    mol = _rdmol(molecule, heavy_only=True)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return rdkit.Chem.MolToSmiles(mol)


def _rdmol(molecule: SmallMolecule, *, heavy_only: bool):
    """An RDKit graph from a :class:`SmallMolecule`; hydrogens dropped when
    *heavy_only*.  Bond order 4 (SDF aromatic) becomes an aromatic bond."""
    rdkit = _rdkit()
    if molecule.molblock:
        from delfin.tools.bio._prepare import _from_small_molecule
        mol = _from_small_molecule(molecule)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        return rdkit.Chem.RemoveAllHs(mol) if heavy_only else mol
    heavy = [(i, a) for i, a in enumerate(molecule.atoms)
             if not (heavy_only and a.is_hydrogen)]
    renumber = {old: new for new, (old, _) in enumerate(heavy)}
    rw = rdkit.Chem.RWMol()
    for _, atom in heavy:
        rw.AddAtom(rdkit.Chem.Atom(atom.element))
    for i, j, order in molecule.bonds:
        if int(i) in renumber and int(j) in renumber:
            bond_type = {1: rdkit.Chem.BondType.SINGLE,
                         2: rdkit.Chem.BondType.DOUBLE,
                         3: rdkit.Chem.BondType.TRIPLE}.get(
                             int(order), rdkit.Chem.BondType.AROMATIC)
            rw.AddBond(renumber[int(i)], renumber[int(j)], bond_type)
    mol = rw.GetMol()
    rdkit.Chem.SanitizeMol(mol)
    return mol


def heavy_atom_pairs(pose: SmallMolecule,
                     reference: SmallMolecule) -> Tuple[Tuple[int, int], ...]:
    """Pair the pose's heavy atoms onto the reference's, symmetry-aware.

    Atoms are grouped by their canonical rank in the heavy-atom graph;
    within each equal-rank class the pairing that minimises the RMSD is
    chosen, so the result is the best mapping the molecule's symmetry
    allows.  Raises :class:`StructureFormatError` when the two are not the
    same molecule, or when the heavy-atom counts differ.
    """
    if not _is_same_molecule(pose, reference):
        raise StructureFormatError(
            "the two structures are not the same molecule: canonical "
            f"heavy-atom SMILES {_canonical_smiles(pose)!r} vs "
            f"{_canonical_smiles(reference)!r}"
        )
    pose_heavy = [i for i, a in enumerate(pose.atoms) if not a.is_hydrogen]
    ref_heavy = [i for i, a in enumerate(reference.atoms)
                 if not a.is_hydrogen]
    if len(pose_heavy) != len(ref_heavy):
        raise StructureFormatError(
            f"the two structures hold {len(pose_heavy)} and "
            f"{len(ref_heavy)} heavy atoms; the same molecule must have "
            "the same number"
        )
    rdkit = _rdkit()
    pose_mol = _rdmol(pose, heavy_only=True)
    ref_mol = _rdmol(reference, heavy_only=True)
    pose_ranks = list(rdkit.Chem.CanonicalRankAtoms(pose_mol,
                                                    breakTies=False))
    ref_ranks = list(rdkit.Chem.CanonicalRankAtoms(ref_mol,
                                                   breakTies=False))
    pose_class: Dict[int, List[int]] = {}
    for index, rank in zip(pose_heavy, pose_ranks):
        pose_class.setdefault(rank, []).append(index)
    ref_class: Dict[int, List[int]] = {}
    for index, rank in zip(ref_heavy, ref_ranks):
        ref_class.setdefault(rank, []).append(index)
    if set(pose_class) != set(ref_class):
        raise StructureFormatError(
            "the two structures are not the same molecule: their "
            "heavy-atom symmetry classes differ"
        )

    def distance(i: int, j: int) -> float:
        a, b = pose.atoms[i], reference.atoms[j]
        return ((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)

    pairs: List[Tuple[int, int]] = []
    for rank in sorted(pose_class):
        left = pose_class[rank]
        right = ref_class[rank]
        if len(left) != len(right):
            raise StructureFormatError(
                f"symmetry class {rank} holds {len(left)} heavy atoms in "
                f"the pose and {len(right)} in the reference; the two are "
                "not the same molecule"
            )
        best = None
        best_cost = None
        for order in permutations(right):
            cost = sum(distance(i, j) for i, j in zip(left, order))
            if best_cost is None or cost < best_cost:
                best_cost = cost
                best = order
        pairs.extend(zip(left, best))
    return tuple(pairs)


def pose_rmsd(pose: SmallMolecule, reference: SmallMolecule,
              pairs: Optional[Sequence[Tuple[int, int]]] = None) -> float:
    """The heavy-atom RMSD between two poses, in Angstrom, unaligned.

    ``pairs`` comes from :func:`heavy_atom_pairs` unless the caller already
    has its own mapping (one pose atom index onto one reference atom index,
    heavy atoms only).  Raises :class:`StructureFormatError` when the two
    are different molecules.
    """
    if pairs is None:
        pairs = heavy_atom_pairs(pose, reference)
    if not pairs:
        raise StructureFormatError(
            "the two structures have no heavy atoms to compare"
        )
    total = 0.0
    for i, j in pairs:
        a, b = pose.atoms[i], reference.atoms[j]
        total += (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2
    return math.sqrt(total / len(pairs))


__all__ = ["heavy_atom_pairs", "pose_rmsd"]
