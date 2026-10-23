"""Hapto fragment extraction, embedding and rigid alignment for hybrid (hapto plus sigma) complexes in the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.converter_flags import (
    _PIPELINE_SEEDS,
)
from delfin.manta.embed_timeout import (
    _embed_with_timeout,
)
from delfin.manta.ml_tables import (
    AllChem,
    Chem,
    Point3D,
    RDKIT_AVAILABLE,
    _METAL_SET,
    _get_ml_bond_length,
)

logger = get_logger("delfin.smiles_converter")


@dataclass
class _HybridHaptoFragment:
    """Metal-free fragment plus atom mappings for hybrid hapto assembly."""

    atom_indices: List[int]
    donor_atom_indices: List[int]
    metal_neighbor_indices: List[int]
    bridging_donor_indices: List[int]
    hapto_group_ids: List[int]
    anchor_atom_indices: List[int]
    original_to_fragment: Dict[int, int]
    fragment_to_original: Dict[int, int]
    fragment_mol: object
    use_scaffold_only: bool = False


@dataclass
class _HybridHaptoDecomposition:
    """Fragmented view of a metal complex with hapto metadata preserved."""

    metal_indices: List[int]
    donor_atom_indices: List[int]
    hapto_groups: List[Tuple[int, List[int]]]
    fragments: List[_HybridHaptoFragment]


@dataclass
class _PrimaryOrganometalModule:
    """Primary hapto-metal module with eligible non-hapto donor fragments."""

    metal_idx: int
    hapto_group_ids: List[int]
    donor_atom_indices: List[int]
    correlated_fragment_indices: List[int]
    terminal_fragment_indices: List[int]


def _copy_fragment_atom(atom, atom_map_num: int):
    """Create a detached copy of an atom with its valence/H metadata."""
    new_atom = Chem.Atom(atom.GetAtomicNum())
    new_atom.SetFormalCharge(atom.GetFormalCharge())
    new_atom.SetNumExplicitHs(atom.GetNumExplicitHs())
    new_atom.SetNoImplicit(atom.GetNoImplicit())
    new_atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons())
    new_atom.SetIsAromatic(atom.GetIsAromatic())
    new_atom.SetChiralTag(atom.GetChiralTag())
    new_atom.SetAtomMapNum(atom_map_num)
    try:
        new_atom.SetHybridization(atom.GetHybridization())
    except Exception:
        pass
    try:
        new_atom.SetIsotope(atom.GetIsotope())
    except Exception:
        pass
    return new_atom


def _extract_fragment_mol(mol, atom_indices: List[int]):
    """Copy a connected metal-free fragment out of ``mol``."""
    if not RDKIT_AVAILABLE or mol is None or not atom_indices:
        return None, {}, {}

    atom_set = set(atom_indices)
    rw = Chem.RWMol()
    original_to_fragment: Dict[int, int] = {}
    fragment_to_original: Dict[int, int] = {}

    for old_idx in atom_indices:
        old_atom = mol.GetAtomWithIdx(old_idx)
        new_idx = rw.AddAtom(_copy_fragment_atom(old_atom, old_idx + 1))
        original_to_fragment[old_idx] = new_idx
        fragment_to_original[new_idx] = old_idx

    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if begin_idx not in atom_set or end_idx not in atom_set:
            continue
        new_begin = original_to_fragment[begin_idx]
        new_end = original_to_fragment[end_idx]
        rw.AddBond(new_begin, new_end, bond.GetBondType())
        new_bond = rw.GetBondBetweenAtoms(new_begin, new_end)
        if new_bond is None:
            continue
        new_bond.SetIsAromatic(bond.GetIsAromatic())
        new_bond.SetIsConjugated(bond.GetIsConjugated())
        try:
            new_bond.SetBondDir(bond.GetBondDir())
            new_bond.SetStereo(bond.GetStereo())
            stereo_atoms = list(bond.GetStereoAtoms())
            if len(stereo_atoms) == 2:
                new_bond.SetStereoAtoms(
                    original_to_fragment[stereo_atoms[0]],
                    original_to_fragment[stereo_atoms[1]],
                )
        except Exception:
            pass

    frag = rw.GetMol()
    try:
        frag.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    try:
        Chem.FastFindRings(frag)
    except Exception:
        pass
    return frag, original_to_fragment, fragment_to_original


def _choose_fragment_anchor_atoms(
    mol,
    frag_atoms: List[int],
    donor_indices: set,
    hapto_groups: List[Tuple[int, List[int]]],
    hapto_atom_to_group: Dict[int, int],
) -> Tuple[List[int], List[int], List[int]]:
    """Choose chemically meaningful anchors for rigid fragment alignment."""
    frag_set = set(frag_atoms)
    frag_donors = sorted(idx for idx in frag_atoms if idx in donor_indices)
    frag_group_ids = sorted({hapto_atom_to_group[idx] for idx in frag_atoms if idx in hapto_atom_to_group})

    anchors: List[int] = []
    for gid in frag_group_ids:
        for atom_idx in hapto_groups[gid][1]:
            if atom_idx in frag_set and atom_idx not in anchors:
                anchors.append(atom_idx)

    for atom_idx in frag_donors:
        if atom_idx not in anchors:
            anchors.append(atom_idx)

    seed_atoms = list(anchors) if anchors else list(frag_donors) if frag_donors else list(frag_atoms)
    for atom_idx in seed_atoms:
        for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in frag_set and nbr_idx not in anchors and nbr.GetAtomicNum() > 1:
                anchors.append(nbr_idx)
                if len(anchors) >= 6:
                    break
        if len(anchors) >= 6:
            break

    if len(anchors) < 3:
        for atom_idx in frag_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() <= 1 or atom_idx in anchors:
                continue
            anchors.append(atom_idx)
            if len(anchors) >= 6:
                break

    return frag_donors, frag_group_ids, anchors


def _decompose_hapto_complex(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
) -> Optional[_HybridHaptoDecomposition]:
    """Split a metal complex into ligand fragments while keeping hapto metadata."""
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return None

    metal_indices = sorted(
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() in _METAL_SET
    )
    if not metal_indices:
        return None

    metal_set = set(metal_indices)
    donor_indices: set = set()
    bonds_to_remove: List[Tuple[int, int]] = []
    hapto_atom_to_group: Dict[int, int] = {}
    for group_idx, (_metal_idx, group_atoms) in enumerate(hapto_groups):
        for atom_idx in group_atoms:
            hapto_atom_to_group[atom_idx] = group_idx

    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if begin_idx in metal_set:
            donor_indices.add(end_idx)
            bonds_to_remove.append((begin_idx, end_idx))
        elif end_idx in metal_set:
            donor_indices.add(begin_idx)
            bonds_to_remove.append((begin_idx, end_idx))

    if not bonds_to_remove:
        return None

    rw = Chem.RWMol(mol)
    for begin_idx, end_idx in bonds_to_remove:
        if rw.GetBondBetweenAtoms(begin_idx, end_idx) is not None:
            rw.RemoveBond(begin_idx, end_idx)
    for metal_idx in sorted(metal_indices, reverse=True):
        rw.RemoveAtom(metal_idx)

    metal_free = rw.GetMol()
    try:
        metal_free.UpdatePropertyCache(strict=False)
    except Exception:
        pass

    old_to_new: Dict[int, int] = {}
    new_to_old: Dict[int, int] = {}
    removed_count = 0
    for old_idx in range(mol.GetNumAtoms()):
        if old_idx in metal_set:
            removed_count += 1
            continue
        new_idx = old_idx - removed_count
        old_to_new[old_idx] = new_idx
        new_to_old[new_idx] = old_idx

    try:
        frags = list(Chem.GetMolFrags(metal_free, asMols=False, sanitizeFrags=False))
    except Exception:
        frags = []

    fragments: List[_HybridHaptoFragment] = []
    for frag_new in frags:
        frag_atoms = sorted(new_to_old[idx] for idx in frag_new if idx in new_to_old)
        if not frag_atoms:
            continue
        fragment_mol, original_to_fragment, fragment_to_original = _extract_fragment_mol(mol, frag_atoms)
        if fragment_mol is None:
            continue
        frag_donors, frag_group_ids, anchors = _choose_fragment_anchor_atoms(
            mol,
            frag_atoms,
            donor_indices,
            hapto_groups,
            hapto_atom_to_group,
        )
        metal_neighbors = sorted({
            nbr.GetIdx()
            for atom_idx in frag_donors
            for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors()
            if nbr.GetSymbol() in _METAL_SET
        })
        bridging_donors = sorted(
            atom_idx
            for atom_idx in frag_donors
            if sum(
                1
                for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors()
                if nbr.GetSymbol() in _METAL_SET
            ) >= 2
        )
        use_scaffold_only = len(metal_neighbors) >= 2 or bool(bridging_donors)
        fragments.append(
            _HybridHaptoFragment(
                atom_indices=frag_atoms,
                donor_atom_indices=frag_donors,
                metal_neighbor_indices=metal_neighbors,
                bridging_donor_indices=bridging_donors,
                hapto_group_ids=frag_group_ids,
                anchor_atom_indices=anchors,
                original_to_fragment=original_to_fragment,
                fragment_to_original=fragment_to_original,
                fragment_mol=fragment_mol,
                use_scaffold_only=use_scaffold_only,
            )
        )

    if not fragments:
        return None

    fragments.sort(
        key=lambda frag: (
            frag.use_scaffold_only,
            -len(frag.hapto_group_ids),
            -len(frag.anchor_atom_indices),
            -len(frag.donor_atom_indices),
            -len(frag.atom_indices),
        )
    )
    return _HybridHaptoDecomposition(
        metal_indices=metal_indices,
        donor_atom_indices=sorted(donor_indices),
        hapto_groups=hapto_groups,
        fragments=fragments,
    )


def _embed_hybrid_fragment(fragment_mol):
    """Embed a detached ligand fragment while preserving its local chemistry."""
    if not RDKIT_AVAILABLE or fragment_mol is None:
        return None

    work = Chem.Mol(fragment_mol)
    work.RemoveAllConformers()

    for atom in work.GetAtoms():
        if atom.GetDegree() < 3 and atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        if (
            atom.GetSymbol() == 'O'
            and atom.GetFormalCharge() == -1
            and any(b.GetBondType() == Chem.BondType.DOUBLE for b in atom.GetBonds())
        ):
            atom.SetFormalCharge(0)

    try:
        work.UpdatePropertyCache(strict=False)
    except Exception:
        pass

    seeds = _PIPELINE_SEEDS[:3]
    embedded = False
    for seed in seeds:
        try:
            params = AllChem.ETKDGv3()
            params.randomSeed = seed
            params.useRandomCoords = True
            params.enforceChirality = False
            result = _embed_with_timeout(work, params)
        except Exception:
            result = -1
        if result == 0:
            embedded = True
            break

    if not embedded:
        try:
            fallback_params = AllChem.EmbedParameters()
            fallback_params.useRandomCoords = True
            fallback_params.randomSeed = 42
            result = _embed_with_timeout(work, fallback_params)
        except Exception:
            result = -1
        embedded = result == 0

    if not embedded:
        return None

    try:
        AllChem.UFFOptimizeMolecule(work, maxIters=200)
    except Exception:
        try:
            AllChem.MMFFOptimizeMolecule(work, maxIters=200)
        except Exception:
            pass

    return work


def _detect_primary_organometal_module(
    decomposition: Optional[_HybridHaptoDecomposition],
    hapto_groups: List[Tuple[int, List[int]]],
) -> Optional[_PrimaryOrganometalModule]:
    """Detect a primary hapto center with buildable non-hapto donor blocks."""
    if decomposition is None or not hapto_groups:
        return None

    hapto_metals = sorted({metal_idx for metal_idx, _grp in hapto_groups})
    if len(hapto_metals) != 1:
        return None

    metal_idx = hapto_metals[0]
    hapto_group_ids = [
        group_idx for group_idx, (group_metal_idx, _group_atoms) in enumerate(hapto_groups)
        if group_metal_idx == metal_idx
    ]
    if not hapto_group_ids:
        return None

    correlated_fragment_indices: List[int] = []
    terminal_fragment_indices: List[int] = []
    donor_atom_indices: List[int] = []
    for frag_idx, fragment in enumerate(decomposition.fragments):
        if metal_idx not in fragment.metal_neighbor_indices:
            continue
        if any(group_id in hapto_group_ids for group_id in fragment.hapto_group_ids):
            continue
        if fragment.use_scaffold_only:
            continue
        if any(neighbor_idx != metal_idx for neighbor_idx in fragment.metal_neighbor_indices):
            continue
        frag_donors = [
            donor_idx for donor_idx in fragment.donor_atom_indices
            if donor_idx not in donor_atom_indices
        ]
        if not frag_donors:
            continue
        donor_atom_indices.extend(frag_donors)
        if len(fragment.donor_atom_indices) >= 2:
            correlated_fragment_indices.append(frag_idx)
        elif len(fragment.donor_atom_indices) == 1:
            terminal_fragment_indices.append(frag_idx)

    if not correlated_fragment_indices and not terminal_fragment_indices:
        return None

    return _PrimaryOrganometalModule(
        metal_idx=metal_idx,
        hapto_group_ids=hapto_group_ids,
        donor_atom_indices=sorted(donor_atom_indices),
        correlated_fragment_indices=correlated_fragment_indices,
        terminal_fragment_indices=terminal_fragment_indices,
    )


def _primary_organometal_module_quality_ok(
    mol,
    decomposition: Optional[_HybridHaptoDecomposition],
    module: Optional[_PrimaryOrganometalModule],
) -> bool:
    """Cheap quality gate for the primary-metal organometal module path."""
    if not RDKIT_AVAILABLE or mol is None or decomposition is None or module is None:
        return False
    try:
        import numpy as np
    except ImportError:
        return False

    try:
        conf = mol.GetConformer(0)
    except Exception:
        return False

    metal_idx = int(module.metal_idx)
    metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
    metal_pos = np.array(conf.GetAtomPosition(metal_idx), dtype=float)

    donor_set = set(module.donor_atom_indices)
    for donor_idx in sorted(donor_set):
        donor_pos = np.array(conf.GetAtomPosition(donor_idx), dtype=float)
        target_len = float(_get_ml_bond_length(metal_sym, mol.GetAtomWithIdx(donor_idx).GetSymbol()))
        dist = float(np.linalg.norm(donor_pos - metal_pos))
        max_err = 0.70 if mol.GetAtomWithIdx(donor_idx).GetSymbol() == 'C' else 0.60
        if abs(dist - target_len) > max_err:
            return False

    inspected_atoms: set = set()
    for frag_idx in module.correlated_fragment_indices + module.terminal_fragment_indices:
        if not (0 <= int(frag_idx) < len(decomposition.fragments)):
            continue
        fragment = decomposition.fragments[int(frag_idx)]
        for atom_idx in fragment.atom_indices:
            if atom_idx in donor_set or atom_idx in inspected_atoms:
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
                continue
            inspected_atoms.add(atom_idx)
            atom_pos = np.array(conf.GetAtomPosition(atom_idx), dtype=float)
            dist = float(np.linalg.norm(atom_pos - metal_pos))
            min_allowed = max(
                1.65,
                0.90 * float(_get_ml_bond_length(metal_sym, atom.GetSymbol())),
            )
            if atom.GetIsAromatic():
                min_allowed = max(min_allowed, 1.85)
            elif atom.GetSymbol() in {'C', 'N', 'O'}:
                min_allowed = max(min_allowed, 1.75)
            if dist < min_allowed - 0.08:
                return False

    return True


def _rotation_matrix_from_vectors(vec_a, vec_b):
    """Return a rotation matrix that maps ``vec_a`` onto ``vec_b``."""
    import numpy as np

    a = np.asarray(vec_a, dtype=float)
    b = np.asarray(vec_b, dtype=float)
    norm_a = float(np.linalg.norm(a))
    norm_b = float(np.linalg.norm(b))
    if norm_a < 1e-12 or norm_b < 1e-12:
        return np.eye(3)

    a /= norm_a
    b /= norm_b
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    s = float(np.linalg.norm(v))

    if s < 1e-12:
        if c > 0.0:
            return np.eye(3)
        axis = np.array([1.0, 0.0, 0.0])
        if abs(a[0]) > 0.8:
            axis = np.array([0.0, 1.0, 0.0])
        v = np.cross(a, axis)
        v /= max(float(np.linalg.norm(v)), 1e-12)
        vx = np.array([
            [0.0, -v[2], v[1]],
            [v[2], 0.0, -v[0]],
            [-v[1], v[0], 0.0],
        ])
        return np.eye(3) + 2.0 * vx @ vx

    vx = np.array([
        [0.0, -v[2], v[1]],
        [v[2], 0.0, -v[0]],
        [-v[1], v[0], 0.0],
    ])
    return np.eye(3) + vx + vx @ vx * ((1.0 - c) / (s * s))


def _align_hybrid_fragment_onto_scaffold(
    scaffold_mol,
    fragment: _HybridHaptoFragment,
    embedded_fragment,
) -> bool:
    """Rigidly align an embedded fragment onto scaffold anchor coordinates."""
    if not RDKIT_AVAILABLE or scaffold_mol is None or embedded_fragment is None:
        return False
    try:
        import numpy as np
    except ImportError:
        return False

    try:
        scaffold_conf = scaffold_mol.GetConformer()
        fragment_conf = embedded_fragment.GetConformer()
    except Exception:
        return False

    anchor_pairs: List[Tuple[int, int]] = []
    for original_idx in fragment.anchor_atom_indices:
        fragment_idx = fragment.original_to_fragment.get(original_idx)
        if fragment_idx is None:
            continue
        scaffold_pos = scaffold_conf.GetAtomPosition(original_idx)
        target_vec = np.array([scaffold_pos.x, scaffold_pos.y, scaffold_pos.z], dtype=float)
        if not np.all(np.isfinite(target_vec)):
            continue
        anchor_pairs.append((fragment_idx, original_idx))

    if not anchor_pairs:
        return False

    source = []
    target = []
    for fragment_idx, original_idx in anchor_pairs:
        frag_pos = fragment_conf.GetAtomPosition(fragment_idx)
        source.append([frag_pos.x, frag_pos.y, frag_pos.z])
        scaf_pos = scaffold_conf.GetAtomPosition(original_idx)
        target.append([scaf_pos.x, scaf_pos.y, scaf_pos.z])

    src = np.asarray(source, dtype=float)
    tgt = np.asarray(target, dtype=float)
    all_coords = np.asarray(
        [
            [
                fragment_conf.GetAtomPosition(i).x,
                fragment_conf.GetAtomPosition(i).y,
                fragment_conf.GetAtomPosition(i).z,
            ]
            for i in range(embedded_fragment.GetNumAtoms())
        ],
        dtype=float,
    )

    if len(anchor_pairs) == 1:
        delta = tgt[0] - src[0]
        aligned = all_coords + delta
    elif len(anchor_pairs) == 2:
        rot = _rotation_matrix_from_vectors(src[1] - src[0], tgt[1] - tgt[0])
        aligned = (all_coords - src[0]) @ rot.T + tgt[0]
    else:
        src_centroid = src.mean(axis=0)
        tgt_centroid = tgt.mean(axis=0)
        covariance = (src - src_centroid).T @ (tgt - tgt_centroid)
        try:
            u, _s, vt = np.linalg.svd(covariance)
        except Exception:
            return False
        rot = vt.T @ u.T
        if float(np.linalg.det(rot)) < 0.0:
            vt[-1, :] *= -1.0
            rot = vt.T @ u.T
        aligned = (all_coords - src_centroid) @ rot.T + tgt_centroid

    # Keep scaffold-defining anchors exact and let the subsequent hapto-only
    # local relaxation repair nearby internal strain.
    for pair_idx, (fragment_idx, _original_idx) in enumerate(anchor_pairs):
        aligned[fragment_idx] = tgt[pair_idx]

    for fragment_idx, original_idx in fragment.fragment_to_original.items():
        scaffold_conf.SetAtomPosition(
            original_idx,
            Point3D(
                float(aligned[fragment_idx, 0]),
                float(aligned[fragment_idx, 1]),
                float(aligned[fragment_idx, 2]),
            ),
        )
    return True


def _hapto_primary_donor_indices(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
) -> set:
    """Return non-metal donors bound directly to hapto metal centers."""
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return set()

    hapto_metals = {metal_idx for metal_idx, _grp in hapto_groups}
    donors: set = set()
    for metal_idx in hapto_metals:
        for nbr in mol.GetAtomWithIdx(metal_idx).GetNeighbors():
            if nbr.GetSymbol() not in _METAL_SET:
                donors.add(nbr.GetIdx())
    return donors
