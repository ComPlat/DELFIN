"""ASE-based complex builder v2.

Workflow:
1. Parse a metal-complex SMILES and split into metals + ligand fragments.
2. Track the *actual coordinating atoms* from the original SMILES.
3. Build 3D ligand geometries.
4. Place ligands initially far from coordination sites and relax with
   ASE (elastic donor anchoring + soft steric repulsion + intraligand springs).
5. Enumerate coordination templates/permutations and keep the best candidate.

This is a pragmatic nudge-builder and not a final electronic-structure optimizer.
"""

from __future__ import annotations

import argparse
import itertools
import logging
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.constraints import FixAtoms
from ase.data import atomic_numbers, covalent_radii
from ase.io import write as ase_write
from ase.optimize import BFGS

from delfin.common.logging import configure_logging, get_logger
from delfin.smiles_converter import RDKIT_AVAILABLE, _METAL_SET, contains_metal

if RDKIT_AVAILABLE:
    from rdkit import Chem
    from rdkit.Chem import AllChem

logger = get_logger(__name__)


@dataclass
class MetalInfo:
    """Metal atom metadata parsed from the original complex SMILES."""

    old_index: int
    symbol: str
    charge: int


@dataclass
class LigandInfo:
    """Ligand fragment with 3D atoms and donor-to-metal links."""

    name: str
    smiles: str
    atoms: Atoms
    donor_links: List[Tuple[int, int]]  # (donor_atom_local_idx, metal_idx)
    bond_pairs: List[Tuple[int, int]]


@dataclass
class AnchorSpec:
    """One donor anchor that must be placed at a metal coordination site."""

    ligand_idx: int
    donor_local_idx: int
    metal_idx: int
    donor_symbol: str


@dataclass
class CandidateResult:
    """One optimized candidate structure."""

    tag: str
    atoms: Atoms
    energy: float
    clash_score: float


class ElasticNudgeCalculator(Calculator):
    """Simple elastic field for BUILD COMPLEX2 placement optimization."""

    implemented_properties = ["energy", "forces"]

    def __init__(
        self,
        bond_terms: Sequence[Tuple[int, int, float, float]],
        anchor_terms: Sequence[Tuple[int, np.ndarray, float]],
        repulsion_terms: Sequence[Tuple[int, int, float, float]],
    ):
        super().__init__()
        self.bond_terms = list(bond_terms)
        self.anchor_terms = list(anchor_terms)
        self.repulsion_terms = list(repulsion_terms)

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        positions = atoms.get_positions()
        forces = np.zeros_like(positions)
        energy = 0.0
        eps = 1.0e-12

        # Intraligand bond springs: preserve ligand geometry during nudge.
        for i, j, d0, k_bond in self.bond_terms:
            vec = positions[j] - positions[i]
            dist = float(np.linalg.norm(vec))
            if dist < eps:
                continue
            delta = dist - d0
            energy += 0.5 * k_bond * delta * delta
            fvec = (k_bond * delta / dist) * vec
            forces[i] += fvec
            forces[j] -= fvec

        # Donor anchors: force coordinating atoms toward coordination-site targets.
        for i, target, k_anchor in self.anchor_terms:
            vec = positions[i] - target
            dist = float(np.linalg.norm(vec))
            if dist < eps:
                continue
            energy += 0.5 * k_anchor * dist * dist
            forces[i] -= (k_anchor / dist) * vec * dist

        # Soft steric wall between atoms from different fragments.
        for i, j, cutoff, k_rep in self.repulsion_terms:
            vec = positions[j] - positions[i]
            dist = float(np.linalg.norm(vec))
            if dist < eps:
                dist = eps
            if dist >= cutoff:
                continue
            delta = cutoff - dist
            energy += 0.5 * k_rep * delta * delta
            direction = vec / dist
            fvec = k_rep * delta * direction
            forces[i] -= fvec
            forces[j] += fvec

        self.results["energy"] = float(energy)
        self.results["forces"] = forces


def _normalize(vec: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vec))
    if norm < 1.0e-12:
        return np.array([1.0, 0.0, 0.0], dtype=float)
    return vec / norm


def _rotation_matrix_from_vectors(v_from: np.ndarray, v_to: np.ndarray) -> np.ndarray:
    """Return rotation matrix that maps v_from to v_to."""
    a = _normalize(v_from)
    b = _normalize(v_to)
    cross = np.cross(a, b)
    dot = float(np.dot(a, b))
    cross_norm = float(np.linalg.norm(cross))

    if cross_norm < 1.0e-12:
        if dot > 0.0:
            return np.eye(3)
        # 180° rotation around any axis orthogonal to a.
        ortho = np.array([1.0, 0.0, 0.0], dtype=float)
        if abs(float(np.dot(a, ortho))) > 0.9:
            ortho = np.array([0.0, 1.0, 0.0], dtype=float)
        axis = _normalize(np.cross(a, ortho))
        return _axis_angle_rotation(axis, math.pi)

    axis = cross / cross_norm
    angle = math.atan2(cross_norm, dot)
    return _axis_angle_rotation(axis, angle)


def _axis_angle_rotation(axis: np.ndarray, angle: float) -> np.ndarray:
    x, y, z = axis
    c = math.cos(angle)
    s = math.sin(angle)
    one_c = 1.0 - c
    return np.array(
        [
            [c + x * x * one_c, x * y * one_c - z * s, x * z * one_c + y * s],
            [y * x * one_c + z * s, c + y * y * one_c, y * z * one_c - x * s],
            [z * x * one_c - y * s, z * y * one_c + x * s, c + z * z * one_c],
        ],
        dtype=float,
    )


def _coordination_templates(n_sites: int) -> List[Tuple[str, np.ndarray]]:
    """Return named unit-vector templates for coordination-site directions."""
    if n_sites <= 0:
        return []
    if n_sites == 1:
        return [("single", np.array([[1.0, 0.0, 0.0]], dtype=float))]
    if n_sites == 2:
        return [
            ("linear", np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]], dtype=float)),
        ]
    if n_sites == 3:
        v = np.array(
            [
                [1.0, 0.0, 0.0],
                [-0.5, math.sqrt(3.0) / 2.0, 0.0],
                [-0.5, -math.sqrt(3.0) / 2.0, 0.0],
            ],
            dtype=float,
        )
        return [("trigonal_planar", v)]
    if n_sites == 4:
        tetra = np.array(
            [[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]],
            dtype=float,
        )
        tetra = np.array([_normalize(row) for row in tetra], dtype=float)
        square = np.array(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]],
            dtype=float,
        )
        return [("tetrahedral", tetra), ("square_planar", square)]
    if n_sites == 5:
        tbp = np.array(
            [
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
                [1.0, 0.0, 0.0],
                [-0.5, math.sqrt(3.0) / 2.0, 0.0],
                [-0.5, -math.sqrt(3.0) / 2.0, 0.0],
            ],
            dtype=float,
        )
        sqp = np.array(
            [
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, -1.0, 0.0],
            ],
            dtype=float,
        )
        return [("trigonal_bipyramidal", tbp), ("square_pyramidal", sqp)]
    if n_sites == 6:
        octa = np.array(
            [
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
            ],
            dtype=float,
        )
        return [("octahedral", octa)]

    # Generic fallback for higher coordination numbers: quasi-uniform sphere points.
    points = []
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    for i in range(n_sites):
        z = 1.0 - (2.0 * i + 1.0) / n_sites
        radius = math.sqrt(max(0.0, 1.0 - z * z))
        theta = 2.0 * math.pi * i / phi
        points.append([radius * math.cos(theta), radius * math.sin(theta), z])
    return [(f"spherical_{n_sites}", np.array(points, dtype=float))]


def _ideal_metal_ligand_distance(metal_symbol: str, donor_symbol: str) -> float:
    """Simple idealized metal-donor distance from covalent radii."""
    z_metal = atomic_numbers.get(metal_symbol, 0)
    z_donor = atomic_numbers.get(donor_symbol, 0)
    r_metal = covalent_radii[z_metal] if z_metal > 0 else 1.40
    r_donor = covalent_radii[z_donor] if z_donor > 0 else 0.75
    return float(max(1.7, min(2.6, r_metal + r_donor + 0.20)))


def _get_covalent_radius(symbol: str) -> float:
    z = atomic_numbers.get(symbol, 0)
    if z <= 0:
        return 0.80
    return float(covalent_radii[z])


def _extract_complex_components(
    smiles: str,
) -> Tuple[Optional[List[MetalInfo]], Optional[List[Tuple["Chem.Mol", List[Tuple[int, int]], str]]], Optional[str]]:
    """Parse complex SMILES and preserve coordinating donor atoms."""
    if not RDKIT_AVAILABLE:
        return None, None, "RDKit not available"

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None, None, "Could not parse SMILES"

    metal_old_indices = [
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() in _METAL_SET
    ]
    if not metal_old_indices:
        return None, None, "No metal atom found in SMILES"
    metal_old_set = set(metal_old_indices)
    metal_old_to_local = {old_idx: i for i, old_idx in enumerate(metal_old_indices)}
    metals = []
    for old_idx in metal_old_indices:
        atom = mol.GetAtomWithIdx(old_idx)
        metals.append(
            MetalInfo(old_index=old_idx, symbol=atom.GetSymbol(), charge=atom.GetFormalCharge())
        )

    # Preserve original atom properties for robust ligand reconstruction.
    orig_props = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        orig_props[idx] = {
            "formal_charge": atom.GetFormalCharge(),
            "explicit_h": atom.GetNumExplicitHs(),
            "no_implicit": atom.GetNoImplicit(),
            "rad_e": atom.GetNumRadicalElectrons(),
        }

    edit_mol = Chem.RWMol(mol)
    donor_old_to_metal_olds: Dict[int, set] = {}
    bonds_to_remove = set()
    for bond in edit_mol.GetBonds():
        b = bond.GetBeginAtomIdx()
        e = bond.GetEndAtomIdx()
        b_is_metal = b in metal_old_set
        e_is_metal = e in metal_old_set
        if b_is_metal and e_is_metal:
            bonds_to_remove.add((min(b, e), max(b, e)))
        elif b_is_metal:
            donor_old_to_metal_olds.setdefault(e, set()).add(b)
            bonds_to_remove.add((min(b, e), max(b, e)))
        elif e_is_metal:
            donor_old_to_metal_olds.setdefault(b, set()).add(e)
            bonds_to_remove.add((min(b, e), max(b, e)))

    for b, e in sorted(bonds_to_remove):
        if edit_mol.GetBondBetweenAtoms(b, e) is not None:
            edit_mol.RemoveBond(b, e)

    for metal_old_idx in sorted(metal_old_indices, reverse=True):
        edit_mol.RemoveAtom(metal_old_idx)

    mol_no_metal = edit_mol.GetMol()

    # Index maps between original molecule and metal-removed molecule.
    idx_map_old_to_new: Dict[int, int] = {}
    removed = 0
    for old_idx in range(mol.GetNumAtoms()):
        if old_idx in metal_old_set:
            removed += 1
            continue
        idx_map_old_to_new[old_idx] = old_idx - removed

    try:
        for old_idx, props in orig_props.items():
            if old_idx in metal_old_set:
                continue
            new_idx = idx_map_old_to_new[old_idx]
            a = mol_no_metal.GetAtomWithIdx(new_idx)
            a.SetNumRadicalElectrons(int(props["rad_e"]))
            a.SetFormalCharge(int(props["formal_charge"]))
            if old_idx in donor_old_to_metal_olds:
                # Coordinating atoms should not get auto-hydrogenated after metal removal.
                a.SetNumExplicitHs(0)
                a.SetNoImplicit(True)
            else:
                a.SetNumExplicitHs(int(props["explicit_h"]))
                a.SetNoImplicit(bool(props["no_implicit"]))
            a.SetAtomMapNum(old_idx + 1)  # Keep original index mapping through fragment split.
    except Exception:
        pass

    try:
        mol_no_metal.UpdatePropertyCache(strict=False)
    except Exception:
        pass

    frag_maps: List[Tuple[int, ...]] = []
    frags = Chem.GetMolFrags(
        mol_no_metal,
        asMols=True,
        sanitizeFrags=False,
        fragsMolAtomMapping=frag_maps,
    )
    if not frags:
        return None, None, "No ligand fragments found after removing metals"

    ligands = []
    for frag_idx, frag_mol in enumerate(frags, start=1):
        donor_links: List[Tuple[int, int]] = []
        for atom in frag_mol.GetAtoms():
            atom_local_idx = atom.GetIdx()
            old_idx = atom.GetAtomMapNum() - 1
            if old_idx < 0:
                continue
            if old_idx not in donor_old_to_metal_olds:
                continue
            for metal_old_idx in sorted(donor_old_to_metal_olds[old_idx]):
                donor_links.append((atom_local_idx, metal_old_to_local[metal_old_idx]))

        frag_for_smiles = Chem.Mol(frag_mol)
        for atom in frag_for_smiles.GetAtoms():
            atom.SetAtomMapNum(0)
        frag_smiles = Chem.MolToSmiles(frag_for_smiles, canonical=True)
        ligands.append((frag_mol, donor_links, frag_smiles))

    return metals, ligands, None


def _embed_fragment_to_atoms(
    frag_mol: "Chem.Mol",
    donor_indices: Sequence[int],
    random_seed: int,
) -> Tuple[Optional[Atoms], Optional[List[Tuple[int, int]]], Optional[str]]:
    """Generate a 3D ligand geometry while preserving donor atom indexing."""
    mol = Chem.Mol(frag_mol)

    # Preserve donor atoms without implicit H addition.
    atoms_no_h = set(int(i) for i in donor_indices)
    for atom in mol.GetAtoms():
        if atom.GetNoImplicit():
            atoms_no_h.add(atom.GetIdx())

    try:
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES,
        )
    except Exception:
        pass
    try:
        mol.UpdatePropertyCache(strict=False)
    except Exception:
        pass

    mol_h = None
    try:
        if atoms_no_h:
            add_h_on = [i for i in range(mol.GetNumAtoms()) if i not in atoms_no_h]
            if add_h_on:
                mol_h = Chem.AddHs(mol, onlyOnAtoms=add_h_on)
            else:
                mol_h = Chem.Mol(mol)
        else:
            mol_h = Chem.AddHs(mol)
    except Exception:
        mol_h = Chem.Mol(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = int(random_seed)
    status = AllChem.EmbedMolecule(mol_h, params)
    if status != 0:
        status = AllChem.EmbedMolecule(mol_h, useRandomCoords=True, randomSeed=int(random_seed))
    if status != 0:
        return None, None, "Could not generate 3D coordinates"

    try:
        AllChem.UFFOptimizeMolecule(mol_h, maxIters=300)
    except Exception:
        try:
            AllChem.MMFFOptimizeMolecule(mol_h, maxIters=300)
        except Exception:
            pass

    conf = mol_h.GetConformer()
    symbols = []
    positions = []
    for i in range(mol_h.GetNumAtoms()):
        atom = mol_h.GetAtomWithIdx(i)
        pos = conf.GetAtomPosition(i)
        symbols.append(atom.GetSymbol())
        positions.append([pos.x, pos.y, pos.z])
    atoms = Atoms(symbols=symbols, positions=np.array(positions, dtype=float))
    bond_pairs = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol_h.GetBonds()]
    return atoms, bond_pairs, None


def _initial_metal_positions(metals: Sequence[MetalInfo]) -> np.ndarray:
    """Place metals with simple spacing as fixed coordination centers."""
    n = len(metals)
    if n == 1:
        return np.array([[0.0, 0.0, 0.0]], dtype=float)
    spacing = 3.0
    start = -0.5 * spacing * (n - 1)
    return np.array([[start + i * spacing, 0.0, 0.0] for i in range(n)], dtype=float)


def _anchor_label(anchor: AnchorSpec, ligands: Sequence[LigandInfo]) -> str:
    ligand = ligands[anchor.ligand_idx]
    return f"{ligand.smiles}:{anchor.donor_symbol}"


def _unique_anchor_permutations(
    anchors: Sequence[AnchorSpec],
    ligands: Sequence[LigandInfo],
    limit: int,
) -> List[Tuple[int, ...]]:
    """Generate assignment permutations; deduplicate by equivalent donor labels."""
    n = len(anchors)
    if n <= 1:
        return [tuple(range(n))]
    indices = list(range(n))
    seen_labels = set()
    out = []
    for perm in itertools.permutations(indices):
        label_signature = tuple(_anchor_label(anchors[p], ligands) for p in perm)
        if label_signature in seen_labels:
            continue
        seen_labels.add(label_signature)
        out.append(perm)
        if limit > 0 and len(out) >= limit:
            break
    return out


def _generate_target_sets(
    metals: Sequence[MetalInfo],
    ligands: Sequence[LigandInfo],
    anchors: Sequence[AnchorSpec],
    metal_positions: np.ndarray,
    max_structures: int,
    per_template_perm_limit: int,
) -> List[Tuple[str, Dict[int, np.ndarray]]]:
    """Create target-point assignments for donor anchors across coordination templates."""
    anchors_by_metal: Dict[int, List[int]] = {}
    for idx, anchor in enumerate(anchors):
        anchors_by_metal.setdefault(anchor.metal_idx, []).append(idx)

    options_per_metal = []
    for metal_idx in range(len(metals)):
        anchor_indices = anchors_by_metal.get(metal_idx, [])
        if not anchor_indices:
            continue

        metal_options = []
        anchor_specs_local = [anchors[i] for i in anchor_indices]
        templates = _coordination_templates(len(anchor_indices))
        for template_name, vectors in templates:
            perms = _unique_anchor_permutations(
                anchor_specs_local,
                ligands,
                limit=per_template_perm_limit,
            )
            for perm in perms:
                local_targets: Dict[int, np.ndarray] = {}
                for site_i, perm_pos in enumerate(perm):
                    global_anchor_idx = anchor_indices[perm_pos]
                    anchor = anchors[global_anchor_idx]
                    donor_symbol = anchor.donor_symbol
                    dist = _ideal_metal_ligand_distance(metals[metal_idx].symbol, donor_symbol)
                    target = metal_positions[metal_idx] + vectors[site_i] * dist
                    local_targets[global_anchor_idx] = target
                metal_options.append((template_name, local_targets))

        if not metal_options:
            continue
        options_per_metal.append((metal_idx, metal_options))

    if not options_per_metal:
        return []

    candidate_sets = []
    option_lists = [opts for _, opts in options_per_metal]
    metal_order = [m for m, _ in options_per_metal]

    for combo in itertools.product(*option_lists):
        target_map: Dict[int, np.ndarray] = {}
        tags = []
        for metal_idx, (template_name, local_targets) in zip(metal_order, combo):
            tags.append(f"M{metal_idx + 1}:{template_name}")
            target_map.update(local_targets)
        candidate_sets.append((" | ".join(tags), target_map))
        if max_structures > 0 and len(candidate_sets) >= max_structures:
            break

    return candidate_sets


def _build_candidate_geometry(
    metals: Sequence[MetalInfo],
    ligands: Sequence[LigandInfo],
    anchors: Sequence[AnchorSpec],
    target_map: Dict[int, np.ndarray],
    metal_positions: np.ndarray,
) -> Tuple[Atoms, List[int], List[int], List[int], List[Tuple[int, int]], List[Tuple[int, int]]]:
    """Assemble a combined Atoms object from one target assignment."""
    symbols: List[str] = []
    positions: List[np.ndarray] = []
    fragment_ids: List[int] = []
    metal_global_indices: List[int] = []
    ligand_starts: List[int] = []
    anchored_pairs: List[Tuple[int, int]] = []

    for i, metal in enumerate(metals):
        metal_global_indices.append(len(symbols))
        symbols.append(metal.symbol)
        positions.append(np.array(metal_positions[i], dtype=float))
        fragment_ids.append(i)

    # Prepare quick access: per ligand list of anchor indices.
    anchors_by_ligand: Dict[int, List[int]] = {}
    for anchor_idx, anchor in enumerate(anchors):
        anchors_by_ligand.setdefault(anchor.ligand_idx, []).append(anchor_idx)

    for lig_idx, ligand in enumerate(ligands):
        lig_atoms = ligand.atoms.copy()
        lig_pos = np.array(lig_atoms.get_positions(), dtype=float)
        lig_anchor_indices = anchors_by_ligand.get(lig_idx, [])

        if lig_anchor_indices:
            primary_anchor_idx = lig_anchor_indices[0]
            primary_anchor = anchors[primary_anchor_idx]
            donor_idx = primary_anchor.donor_local_idx
            target = target_map[primary_anchor_idx]
            metal_pos = metal_positions[primary_anchor.metal_idx]
            site_dir = _normalize(target - metal_pos)

            donor_pos = lig_pos[donor_idx]
            com = lig_atoms.get_center_of_mass()
            donor_to_com = _normalize(np.array(com, dtype=float) - donor_pos)
            rot = _rotation_matrix_from_vectors(donor_to_com, site_dir)
            centered = lig_pos - donor_pos
            rotated = centered @ rot.T

            # Start "from distance": donor begins outside the target shell.
            far_target = target + site_dir * (3.5 + 0.35 * lig_idx)
            translated = rotated + far_target
            lig_pos = translated

        lig_start = len(symbols)
        ligand_starts.append(lig_start)
        for atom_i, sym in enumerate(lig_atoms.get_chemical_symbols()):
            symbols.append(sym)
            positions.append(lig_pos[atom_i])
            fragment_ids.append(len(metals) + lig_idx)

        for donor_local_idx, metal_idx in ligand.donor_links:
            donor_global = lig_start + donor_local_idx
            metal_global = metal_global_indices[metal_idx]
            anchored_pairs.append((donor_global, metal_global))

    atoms = Atoms(symbols=symbols, positions=np.array(positions, dtype=float))
    return (
        atoms,
        fragment_ids,
        metal_global_indices,
        ligand_starts,
        anchored_pairs,
        [],  # reserved for future extensions
    )


def _prepare_force_field_terms(
    atoms: Atoms,
    ligands: Sequence[LigandInfo],
    anchors: Sequence[AnchorSpec],
    target_map: Dict[int, np.ndarray],
    ligand_starts: Sequence[int],
    fragment_ids: Sequence[int],
    anchored_pairs: Sequence[Tuple[int, int]],
) -> Tuple[List[Tuple[int, int, float, float]], List[Tuple[int, np.ndarray, float]], List[Tuple[int, int, float, float]]]:
    """Build force-field terms for the elastic nudge optimization."""
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()

    bond_terms: List[Tuple[int, int, float, float]] = []
    for lig_idx, ligand in enumerate(ligands):
        start = ligand_starts[lig_idx]
        for i_local, j_local in ligand.bond_pairs:
            i_global = start + i_local
            j_global = start + j_local
            d0 = float(np.linalg.norm(positions[j_global] - positions[i_global]))
            bond_terms.append((i_global, j_global, d0, 90.0))

    anchor_terms: List[Tuple[int, np.ndarray, float]] = []
    for anchor_idx, anchor in enumerate(anchors):
        donor_global = ligand_starts[anchor.ligand_idx] + anchor.donor_local_idx
        target = np.array(target_map[anchor_idx], dtype=float)
        anchor_terms.append((donor_global, target, 16.0))

    anchored_pair_set = {tuple(sorted(pair)) for pair in anchored_pairs}
    repulsion_terms: List[Tuple[int, int, float, float]] = []
    n_atoms = len(atoms)
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            if fragment_ids[i] == fragment_ids[j]:
                continue
            if (i, j) in anchored_pair_set:
                continue
            r_i = _get_covalent_radius(symbols[i])
            r_j = _get_covalent_radius(symbols[j])
            cutoff = max(1.35, 0.95 * (r_i + r_j))
            repulsion_terms.append((i, j, cutoff, 28.0))

    return bond_terms, anchor_terms, repulsion_terms


def _compute_clash_score(atoms: Atoms, fragment_ids: Sequence[int]) -> float:
    """Lower is better; sums only significant inter-fragment overlaps."""
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    score = 0.0
    n_atoms = len(atoms)
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            if fragment_ids[i] == fragment_ids[j]:
                continue
            dist = float(np.linalg.norm(positions[j] - positions[i]))
            cutoff = 0.90 * (_get_covalent_radius(symbols[i]) + _get_covalent_radius(symbols[j]))
            if dist < cutoff:
                score += cutoff - max(dist, 1.0e-8)
    return float(score)


def _write_xyz_and_start(atoms: Atoms, parent_dir: Path, xyz_name: str) -> None:
    out_xyz = parent_dir / xyz_name
    ase_write(out_xyz, atoms, format="xyz")
    lines = out_xyz.read_text(encoding="utf-8").splitlines()
    coords = lines[2:] if len(lines) >= 2 else []
    (parent_dir / "start.txt").write_text(
        "\n".join(coords).strip() + ("\n" if coords else ""),
        encoding="utf-8",
    )


def run_build_up2(
    smiles: str,
    builder_dir: Path,
    max_structures: int = 0,
    per_template_perm_limit: int = 0,
    max_steps: int = 320,
    fmax: float = 0.08,
    write_candidates: bool = False,
) -> bool:
    """Run BUILD COMPLEX2 workflow."""
    if not RDKIT_AVAILABLE:
        logger.error("RDKit is required for BUILD COMPLEX2")
        return False

    if not contains_metal(smiles):
        logger.error("SMILES does not contain a metal complex")
        return False

    metals, raw_ligands, error = _extract_complex_components(smiles)
    if error:
        logger.error(f"Failed to parse complex: {error}")
        return False
    assert metals is not None
    assert raw_ligands is not None

    # Keep only coordinating ligands (donor links explicitly present in SMILES).
    coordinating_raw_ligands = [entry for entry in raw_ligands if entry[1]]
    if not coordinating_raw_ligands:
        logger.error("No coordinating ligand atoms found in SMILES")
        return False
    if len(coordinating_raw_ligands) != len(raw_ligands):
        logger.info(
            "Skipping %d non-coordinating fragment(s); only coordinating ligands are built.",
            len(raw_ligands) - len(coordinating_raw_ligands),
        )

    builder_dir.mkdir(parents=True, exist_ok=True)

    ligands: List[LigandInfo] = []
    for lig_idx, (frag_mol, donor_links, frag_smiles) in enumerate(coordinating_raw_ligands, start=1):
        donor_indices = sorted({idx for idx, _ in donor_links})
        atoms, bond_pairs, emb_error = _embed_fragment_to_atoms(
            frag_mol=frag_mol,
            donor_indices=donor_indices,
            random_seed=1103 + lig_idx * 17,
        )
        if emb_error or atoms is None or bond_pairs is None:
            logger.error(f"Ligand {lig_idx} embedding failed: {emb_error or 'unknown error'}")
            return False

        # Keep only donor links that are valid in the embedded atom indexing.
        valid_donor_links = [(a_idx, m_idx) for (a_idx, m_idx) in donor_links if a_idx < len(atoms)]
        if not valid_donor_links:
            logger.error(
                "Ligand %d has no valid donor atom indices after embedding; aborting.",
                lig_idx,
            )
            return False

        ligand = LigandInfo(
            name=f"ligand{lig_idx}",
            smiles=frag_smiles,
            atoms=atoms,
            donor_links=valid_donor_links,
            bond_pairs=bond_pairs,
        )
        ligands.append(ligand)
        ase_write(builder_dir / f"{ligand.name}.xyz", atoms, format="xyz")

    anchors: List[AnchorSpec] = []
    for lig_idx, ligand in enumerate(ligands):
        for donor_local_idx, metal_idx in ligand.donor_links:
            donor_symbol = ligand.atoms[donor_local_idx].symbol
            anchors.append(
                AnchorSpec(
                    ligand_idx=lig_idx,
                    donor_local_idx=donor_local_idx,
                    metal_idx=metal_idx,
                    donor_symbol=donor_symbol,
                )
            )
    if not anchors:
        logger.error("No donor anchors found for BUILD COMPLEX2")
        return False

    logger.info("Metals: %s", ", ".join(f"{m.symbol}({m.charge:+d})" for m in metals))
    logger.info("Coordinating ligands: %d", len(ligands))
    logger.info("Donor anchors from SMILES: %d", len(anchors))

    metal_positions = _initial_metal_positions(metals)
    target_sets = _generate_target_sets(
        metals=metals,
        ligands=ligands,
        anchors=anchors,
        metal_positions=metal_positions,
        max_structures=max_structures,
        per_template_perm_limit=per_template_perm_limit,
    )
    if not target_sets:
        logger.error("Could not generate any coordination target sets")
        return False
    logger.info("Enumerating %d candidate coordination arrangements", len(target_sets))

    summary_lines = ["#candidate\ttag\tenergy\tclash"]
    best: Optional[CandidateResult] = None

    for cand_idx, (tag, target_map) in enumerate(target_sets, start=1):
        (
            atoms,
            fragment_ids,
            metal_global_indices,
            ligand_starts,
            anchored_pairs,
            _,
        ) = _build_candidate_geometry(
            metals=metals,
            ligands=ligands,
            anchors=anchors,
            target_map=target_map,
            metal_positions=metal_positions,
        )

        bond_terms, anchor_terms, repulsion_terms = _prepare_force_field_terms(
            atoms=atoms,
            ligands=ligands,
            anchors=anchors,
            target_map=target_map,
            ligand_starts=ligand_starts,
            fragment_ids=fragment_ids,
            anchored_pairs=anchored_pairs,
        )

        atoms.calc = ElasticNudgeCalculator(
            bond_terms=bond_terms,
            anchor_terms=anchor_terms,
            repulsion_terms=repulsion_terms,
        )
        atoms.set_constraint(FixAtoms(indices=metal_global_indices))

        try:
            dyn = BFGS(atoms, logfile=None)
            dyn.run(fmax=fmax, steps=max_steps)
            energy = float(atoms.get_potential_energy())
        except Exception as exc:
            logger.warning("Candidate %d failed during ASE optimization: %s", cand_idx, exc)
            continue

        clash = _compute_clash_score(atoms, fragment_ids)
        summary_lines.append(f"{cand_idx}\t{tag}\t{energy:.8f}\t{clash:.8f}")

        if write_candidates:
            cand_file = builder_dir / f"candidate_{cand_idx:03d}.xyz"
            ase_write(cand_file, atoms, format="xyz")

        cand = CandidateResult(tag=tag, atoms=atoms.copy(), energy=energy, clash_score=clash)
        if best is None:
            best = cand
        else:
            if (cand.clash_score, cand.energy) < (best.clash_score, best.energy):
                best = cand

    (builder_dir / "candidate_summary.tsv").write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )

    if best is None:
        logger.error("No valid candidate could be optimized")
        return False

    best_file = builder_dir / "best_complex2.xyz"
    ase_write(best_file, best.atoms, format="xyz")

    parent_dir = builder_dir.parent
    _write_xyz_and_start(best.atoms, parent_dir=parent_dir, xyz_name="build_complex2.xyz")
    logger.info("Best candidate: %s", best.tag)
    logger.info("Best energy: %.8f, clash: %.8f", best.energy, best.clash_score)
    logger.info("Wrote %s", best_file)
    logger.info("Wrote %s", parent_dir / "build_complex2.xyz")
    logger.info("Wrote %s", parent_dir / "start.txt")
    return True


def _read_smiles_from_input(input_file: Path) -> Tuple[Optional[str], Optional[str]]:
    if not input_file.exists():
        return None, f"Input file not found: {input_file}"
    try:
        content = input_file.read_text(encoding="utf-8")
    except Exception as exc:
        return None, f"Could not read input file: {exc}"

    for line in content.splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            return line, None
    return None, "No SMILES string found in input file"


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="delfin-build2",
        description="ASE-based complex builder (BUILD COMPLEX2)",
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        default="input.txt",
        help="Input file containing a complex SMILES (default: input.txt)",
    )
    parser.add_argument(
        "-d",
        "--directory",
        default="builder2",
        help="Output directory for intermediate files (default: builder2)",
    )
    parser.add_argument(
        "--max-structures",
        type=int,
        default=0,
        help="Maximum number of coordination candidates (0 = all, default: 0)",
    )
    parser.add_argument(
        "--perm-limit",
        type=int,
        default=0,
        help="Max unique donor-permutations per template (0 = all, default: 0)",
    )
    parser.add_argument(
        "--max-steps",
        type=int,
        default=320,
        help="Maximum ASE optimization steps per candidate (default: 320)",
    )
    parser.add_argument(
        "--fmax",
        type=float,
        default=0.08,
        help="ASE BFGS force threshold (default: 0.08 eV/A)",
    )
    parser.add_argument(
        "--write-candidate-files",
        action="store_true",
        help="Write per-candidate XYZ files (default: off)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )

    args = parser.parse_args(argv)
    configure_logging(level=logging.DEBUG if args.verbose else logging.INFO)

    smiles, error = _read_smiles_from_input(Path(args.input_file))
    if error:
        logger.error(error)
        return 1
    assert smiles is not None

    success = run_build_up2(
        smiles=smiles,
        builder_dir=Path(args.directory),
        max_structures=int(args.max_structures),
        per_template_perm_limit=int(args.perm_limit),
        max_steps=max(20, int(args.max_steps)),
        fmax=max(0.01, float(args.fmax)),
        write_candidates=bool(args.write_candidate_files),
    )
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
