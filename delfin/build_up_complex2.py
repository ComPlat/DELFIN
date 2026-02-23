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
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
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


_EMPIRICAL_ML_DISTANCES: Dict[Tuple[str, str], float] = {
    # Cadmium
    ("Cd", "N"): 2.25, ("Cd", "O"): 2.25, ("Cd", "S"): 2.55,
    ("Cd", "Se"): 2.65, ("Cd", "P"): 2.60, ("Cd", "Cl"): 2.55,
    ("Cd", "Br"): 2.65,
    # Zinc
    ("Zn", "N"): 2.05, ("Zn", "O"): 1.98, ("Zn", "S"): 2.30,
    ("Zn", "Se"): 2.40, ("Zn", "P"): 2.40, ("Zn", "Cl"): 2.30,
    # Copper(II)
    ("Cu", "N"): 2.00, ("Cu", "O"): 2.00, ("Cu", "S"): 2.30,
    ("Cu", "Cl"): 2.25, ("Cu", "Br"): 2.35,
    # Nickel
    ("Ni", "N"): 1.90, ("Ni", "O"): 1.85, ("Ni", "S"): 2.20,
    ("Ni", "P"): 2.20, ("Ni", "Cl"): 2.20,
    # Cobalt
    ("Co", "N"): 1.95, ("Co", "O"): 1.95, ("Co", "S"): 2.25,
    ("Co", "P"): 2.25, ("Co", "Cl"): 2.25,
    # Iron
    ("Fe", "N"): 1.95, ("Fe", "O"): 1.95, ("Fe", "S"): 2.25,
    ("Fe", "Cl"): 2.25,
    # Manganese
    ("Mn", "N"): 2.20, ("Mn", "O"): 2.15, ("Mn", "S"): 2.45,
    ("Mn", "Cl"): 2.35,
    # Platinum
    ("Pt", "N"): 2.05, ("Pt", "O"): 2.05, ("Pt", "S"): 2.30,
    ("Pt", "P"): 2.25, ("Pt", "Cl"): 2.30,
    # Palladium
    ("Pd", "N"): 2.05, ("Pd", "O"): 2.05, ("Pd", "S"): 2.30,
    ("Pd", "P"): 2.25, ("Pd", "Cl"): 2.30,
    # Ruthenium
    ("Ru", "N"): 2.10, ("Ru", "O"): 2.05, ("Ru", "S"): 2.35,
    ("Ru", "P"): 2.30, ("Ru", "Cl"): 2.35,
    # Rhodium / Iridium
    ("Rh", "N"): 2.10, ("Rh", "P"): 2.25, ("Rh", "Cl"): 2.35,
    ("Ir", "N"): 2.10, ("Ir", "P"): 2.25, ("Ir", "Cl"): 2.35,
}


def _ideal_metal_ligand_distance(metal_symbol: str, donor_symbol: str) -> float:
    """Return the idealized metal–donor bond distance in Å.

    Uses an empirical lookup table first; falls back to covalent-radii sum.
    """
    key = (metal_symbol, donor_symbol)
    if key in _EMPIRICAL_ML_DISTANCES:
        return _EMPIRICAL_ML_DISTANCES[key]
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


# Bondi Van der Waals radii (Å) – used for inter-fragment clash detection.
_VDW_RADII: Dict[str, float] = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
    "P": 1.80, "F": 1.47, "Cl": 1.75, "Br": 1.85, "I": 1.98,
    "Se": 1.90, "Te": 2.06,
    # Transition metals (Alvarez 2013 / Bondi extended)
    "Cd": 1.58, "Zn": 1.39, "Cu": 1.40, "Ni": 1.63, "Co": 1.70,
    "Fe": 1.56, "Mn": 1.97, "Cr": 1.66, "V": 1.79,
    "Pt": 1.75, "Pd": 1.63, "Ru": 1.46, "Rh": 1.45, "Ir": 1.41,
    "Ag": 1.72, "Au": 1.66, "Hg": 1.55,
}


def _get_vdw_radius(symbol: str) -> float:
    """Return Van der Waals radius in Å (Bondi values; default C-like 1.70)."""
    return _VDW_RADII.get(symbol, 1.70)


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
    metal_global_indices: Sequence[int] = (),
) -> Tuple[List[Tuple[int, int, float, float]], List[Tuple[int, np.ndarray, float]], List[Tuple[int, int, float, float]]]:
    """Build force-field terms for the elastic nudge optimization."""
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()

    # CRITICAL: use bond lengths from the ORIGINAL RDKit geometry, NOT from the
    # current (PSO-distorted) positions.  If PSO compressed a C-H to 0.88 Å,
    # using that as d0 would lock in the distortion through BFGS.
    bond_terms: List[Tuple[int, int, float, float]] = []
    for lig_idx, ligand in enumerate(ligands):
        start = ligand_starts[lig_idx]
        orig_pos = ligand.atoms.get_positions()
        for i_local, j_local in ligand.bond_pairs:
            i_global = start + i_local
            j_global = start + j_local
            d0 = float(np.linalg.norm(orig_pos[j_local] - orig_pos[i_local]))
            bond_terms.append((i_global, j_global, d0, 90.0))

    anchor_terms: List[Tuple[int, np.ndarray, float]] = []
    for anchor_idx, anchor in enumerate(anchors):
        donor_global = ligand_starts[anchor.ligand_idx] + anchor.donor_local_idx
        target = np.array(target_map[anchor_idx], dtype=float)
        anchor_terms.append((donor_global, target, 16.0))

    # Collect donor global indices so they are excluded from metal repulsion below.
    anchored_donor_set: set = {d for d, _m in anchored_pairs}
    anchored_pair_set = {tuple(sorted(pair)) for pair in anchored_pairs}
    metal_set = set(metal_global_indices)
    repulsion_terms: List[Tuple[int, int, float, float]] = []
    n_atoms = len(atoms)
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            if fragment_ids[i] == fragment_ids[j]:
                continue
            if (i, j) in anchored_pair_set:
                continue
            vdw_i = _get_vdw_radius(symbols[i])
            vdw_j = _get_vdw_radius(symbols[j])

            i_is_metal = i in metal_set
            j_is_metal = j in metal_set
            if i_is_metal or j_is_metal:
                # Metal ↔ non-donor repulsion: prevents ring C/H from pointing
                # toward the metal (e.g. pyridyl ring oriented backwards).
                # Donor atoms CAN approach legitimately → skip them here.
                non_metal_idx = j if i_is_metal else i
                if non_metal_idx in anchored_donor_set:
                    continue  # donor atom — no repulsion from metal
                # 85 % of VdW sum: Cd-C → ~2.79 Å, Cd-H → ~2.36 Å.
                cutoff = max(2.0, 0.85 * (vdw_i + vdw_j))
                repulsion_terms.append((i, j, cutoff, 60.0))
            else:
                # Standard inter-ligand repulsion.
                # 65 % of VdW sum: C-C → 2.21 Å, H-C → 1.89 Å, H-H → 1.56 Å.
                cutoff = max(1.5, 0.65 * (vdw_i + vdw_j))
                repulsion_terms.append((i, j, cutoff, 35.0))

    return bond_terms, anchor_terms, repulsion_terms


def _check_topology_integrity(
    atoms: Atoms,
    fragment_ids: Sequence[int],
    ligands: Sequence[LigandInfo],
    ligand_starts: Sequence[int],
    metal_global_indices: Sequence[int],
    anchored_pairs: Sequence[Tuple[int, int]],
) -> Tuple[bool, str]:
    """Verify that no inter-fragment "bonds" have formed and intra-ligand geometry is intact.

    Returns ``(ok, reason)`` where *ok* is True when the structure is valid.
    A candidate is rejected when:

    * Any two atoms from **different non-metal fragments** are within the
      covalent bonding distance (0.85 × sum of covalent radii) – this would
      indicate that ligands have "reacted" with each other.
    * Any intra-ligand bond has stretched to more than 1.6 × its original
      length – indicating topology corruption during BFGS.

    Metal ↔ donor bonds are explicitly allowed (they are in *anchored_pairs*).
    """
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    metal_set = set(metal_global_indices)
    anchored_set = {tuple(sorted(p)) for p in anchored_pairs}

    # 1. Inter-ligand covalent contact check.
    n = len(atoms)
    for i in range(n):
        if i in metal_set:
            continue
        for j in range(i + 1, n):
            if j in metal_set:
                continue
            if fragment_ids[i] == fragment_ids[j]:
                continue
            if tuple(sorted((i, j))) in anchored_set:
                continue
            dist = float(np.linalg.norm(positions[j] - positions[i]))
            bond_cutoff = 0.85 * (_get_covalent_radius(symbols[i]) + _get_covalent_radius(symbols[j]))
            if dist < bond_cutoff:
                return (
                    False,
                    f"inter-ligand bond formed: {symbols[i]}[{i}]-{symbols[j]}[{j}] "
                    f"at {dist:.3f} Å (cutoff {bond_cutoff:.3f} Å)",
                )

    # 2. Intra-ligand bond distortion check (stretch AND compress).
    for lig_idx, ligand in enumerate(ligands):
        start = ligand_starts[lig_idx]
        orig_pos = ligand.atoms.get_positions()
        curr_pos = positions[start : start + len(ligand.atoms)]
        for i_local, j_local in ligand.bond_pairs:
            orig_d = float(np.linalg.norm(orig_pos[j_local] - orig_pos[i_local]))
            curr_d = float(np.linalg.norm(curr_pos[j_local] - curr_pos[i_local]))
            if orig_d < 0.1:
                continue
            i_g = start + i_local
            j_g = start + j_local
            if curr_d > 1.6 * orig_d:
                return (
                    False,
                    f"intra-ligand bond stretched: {symbols[i_g]}[{i_g}]-{symbols[j_g]}[{j_g}] "
                    f"{orig_d:.3f}→{curr_d:.3f} Å (>1.6×)",
                )
            if curr_d < 0.65 * orig_d:
                return (
                    False,
                    f"intra-ligand bond compressed: {symbols[i_g]}[{i_g}]-{symbols[j_g]}[{j_g}] "
                    f"{orig_d:.3f}→{curr_d:.3f} Å (<0.65×)",
                )

    return True, "OK"


def _compute_clash_score(atoms: Atoms, fragment_ids: Sequence[int]) -> float:
    """Lower is better; sums inter-fragment VdW overlaps (75 % of VdW sum cutoff).

    Using Van-der-Waals radii ensures realistic detection: an H···C contact
    at 1.4 Å (impossible non-bonded) is correctly flagged, whereas the old
    covalent-radius cutoff (~1.0 Å) would have missed it entirely.
    """
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    score = 0.0
    n_atoms = len(atoms)
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            if fragment_ids[i] == fragment_ids[j]:
                continue
            dist = float(np.linalg.norm(positions[j] - positions[i]))
            cutoff = 0.75 * (_get_vdw_radius(symbols[i]) + _get_vdw_radius(symbols[j]))
            if dist < cutoff:
                score += cutoff - max(dist, 1.0e-8)
    return float(score)


# ── xTB single-point ranking ─────────────────────────────────────────────────


def _xtb_singlepoint(
    atoms: Atoms,
    charge: int,
    multiplicity: int,
    xtb_binary: str,
    nthreads: int = 1,
    gfn_level: int = 2,
) -> Optional[float]:
    """Run one xTB single-point and return the total energy in Hartree.

    Returns *None* if xTB fails or the binary is not found.
    """
    uhf = multiplicity - 1  # number of unpaired electrons
    with tempfile.TemporaryDirectory(prefix="delfin_xtb_") as tmpdir:
        xyz_path = os.path.join(tmpdir, "mol.xyz")
        ase_write(xyz_path, atoms, format="xyz")
        cmd = [
            xtb_binary, "mol.xyz",
            "--sp", "--gfn", str(gfn_level),
            "--chrg", str(charge),
            "--uhf", str(uhf),
            "--norestart",
        ]
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(nthreads)
        env["MKL_NUM_THREADS"] = str(nthreads)
        env["XTBPATH"] = os.path.dirname(xtb_binary)
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=180,
                cwd=tmpdir,
                env=env,
            )
            for line in result.stdout.splitlines():
                if "TOTAL ENERGY" in line:
                    parts = line.split()
                    for idx, tok in enumerate(parts):
                        if tok == "ENERGY" and idx + 1 < len(parts):
                            try:
                                return float(parts[idx + 1])
                            except ValueError:
                                pass
        except Exception:
            pass
    return None


def _xtb_optimize(
    atoms: "Atoms",
    charge: int,
    multiplicity: int,
    xtb_binary: str,
    fixed_indices: Optional[List[int]] = None,
    opt_level: str = "loose",
    gfn_level: int = 2,
    nthreads: int = 1,
) -> Optional[Tuple["Atoms", float]]:
    """Run xTB geometry optimization; return (optimized_atoms, energy_Ha) or None.

    The metal atoms listed in *fixed_indices* are constrained to their input
    positions via xTB's ``$fix`` block, preserving the coordination geometry
    while allowing the ligands to flex/deform freely.  This is the key
    advantage over BFGS: xTB knows about chemical bonding, so intra-ligand
    topology is automatically preserved and no inter-ligand bonds can form.
    """
    from ase.io import read as ase_read

    uhf = multiplicity - 1
    with tempfile.TemporaryDirectory(prefix="delfin_xtbopt_") as tmpdir:
        xyz_path = os.path.join(tmpdir, "mol.xyz")
        ase_write(xyz_path, atoms, format="xyz")

        cmd = [
            xtb_binary, "mol.xyz",
            "--opt", opt_level,
            "--gfn", str(gfn_level),
            "--chrg", str(charge),
            "--uhf", str(uhf),
            "--norestart",
        ]

        # Fix metal positions via xTB input file (atoms are 1-indexed in xTB).
        if fixed_indices:
            input_path = os.path.join(tmpdir, "xtb.inp")
            atom_list = ",".join(str(i + 1) for i in sorted(fixed_indices))
            with open(input_path, "w", encoding="utf-8") as fh:
                fh.write(f"$fix\n  atoms: {atom_list}\n$end\n")
            cmd += ["--input", "xtb.inp"]

        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(nthreads)
        env["MKL_NUM_THREADS"] = str(nthreads)

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600,
                cwd=tmpdir,
                env=env,
            )
            opt_xyz = os.path.join(tmpdir, "xtbopt.xyz")
            if not os.path.exists(opt_xyz):
                logger.debug(
                    "xTB opt: xtbopt.xyz not found; stderr tail: %s",
                    result.stderr[-400:] if result.stderr else "",
                )
                return None
            opt_atoms = ase_read(opt_xyz)
            # Parse total energy from xTB stdout.
            energy: Optional[float] = None
            for line in result.stdout.splitlines():
                if "TOTAL ENERGY" in line:
                    parts = line.split()
                    for tok_idx, tok in enumerate(parts):
                        if tok == "ENERGY" and tok_idx + 1 < len(parts):
                            try:
                                energy = float(parts[tok_idx + 1])
                            except ValueError:
                                pass
            if energy is None:
                logger.debug("xTB opt: could not parse energy from stdout.")
                return None
            return opt_atoms, energy
        except Exception as exc:
            logger.debug("xTB opt exception: %s", exc)
            return None


def _rank_candidates_by_xtb(
    candidates: List["CandidateResult"],
    total_charge: int,
    multiplicity: int,
    xtb_binary: str,
    n_parallel: int,
    gfn_level: int = 2,
) -> List["CandidateResult"]:
    """Re-rank *candidates* using parallel xTB single-points.

    Candidates for which xTB fails keep their original (elastic) energy.
    """
    logger.info(
        "Running xTB GFN%d single-points on %d candidates (%d parallel)...",
        gfn_level, len(candidates), n_parallel,
    )

    def _score(item: Tuple[int, "CandidateResult"]) -> Tuple[int, Optional[float]]:
        idx, cand = item
        e = _xtb_singlepoint(
            cand.atoms, charge=total_charge, multiplicity=multiplicity,
            xtb_binary=xtb_binary, gfn_level=gfn_level,
        )
        return idx, e

    xtb_energies: Dict[int, Optional[float]] = {}
    with ThreadPoolExecutor(max_workers=max(1, n_parallel)) as pool:
        futures = {pool.submit(_score, (i, c)): i for i, c in enumerate(candidates)}
        for fut in as_completed(futures):
            idx, e = fut.result()
            xtb_energies[idx] = e
            if e is not None:
                logger.debug("xTB candidate %d: %.8f Ha", idx, e)
            else:
                logger.debug("xTB candidate %d: failed, keeping elastic energy", idx)

    # Build re-ranked candidate list.
    result = []
    for i, cand in enumerate(candidates):
        e = xtb_energies.get(i)
        if e is not None:
            result.append(CandidateResult(
                tag=cand.tag, atoms=cand.atoms,
                energy=e, clash_score=cand.clash_score,
            ))
        else:
            result.append(cand)

    n_ok = sum(1 for e in xtb_energies.values() if e is not None)
    logger.info("xTB: %d/%d successful", n_ok, len(candidates))
    return result


# ── Swarm / Particle-Swarm-Optimisation (PSO) ligand placer ─────────────────


def _rvec_to_matrix(rvec: np.ndarray) -> np.ndarray:
    """Rodrigues rotation vector → 3×3 rotation matrix."""
    angle = float(np.linalg.norm(rvec))
    if angle < 1.0e-12:
        return np.eye(3, dtype=float)
    return _axis_angle_rotation(rvec / angle, angle)


def _matrix_to_rvec(R: np.ndarray) -> np.ndarray:
    """3×3 rotation matrix → Rodrigues rotation vector (inverse Rodrigues)."""
    cos_t = max(-1.0, min(1.0, (float(R[0, 0] + R[1, 1] + R[2, 2]) - 1.0) / 2.0))
    angle = math.acos(cos_t)
    if angle < 1.0e-12:
        return np.zeros(3, dtype=float)
    sin_t = math.sin(angle)
    if abs(sin_t) < 1.0e-8:
        # 180° rotation: extract axis from diagonal of R = 2*n*nᵀ - I.
        xx = max(0.0, (R[0, 0] + 1.0) / 2.0)
        yy = max(0.0, (R[1, 1] + 1.0) / 2.0)
        zz = max(0.0, (R[2, 2] + 1.0) / 2.0)
        x, y, z = math.sqrt(xx), math.sqrt(yy), math.sqrt(zz)
        if R[0, 1] < 0.0:
            y = -y
        if R[0, 2] < 0.0:
            z = -z
        axis = np.array([x, y, z], dtype=float)
        n = float(np.linalg.norm(axis))
        axis = axis / n if n > 1.0e-12 else np.array([1.0, 0.0, 0.0])
        return axis * angle
    axis = np.array(
        [R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]], dtype=float
    )
    axis /= 2.0 * sin_t
    return axis * angle


def _apply_rigid_pose(
    base_positions: np.ndarray,
    base_com: np.ndarray,
    tvec: np.ndarray,
    rvec: np.ndarray,
) -> np.ndarray:
    """Rotate *base_positions* around *base_com* by *rvec*, then translate by *tvec*."""
    R = _rvec_to_matrix(rvec)
    return (base_positions - base_com) @ R.T + base_com + tvec


@dataclass
class SwarmParticle:
    """One PSO particle = rigid-body poses for all ligands simultaneously.

    ``poses[k]`` = ``[tx, ty, tz, rx, ry, rz]`` for ligand *k*, where
    columns 0-2 are the COM translation and columns 3-5 are the Rodrigues
    rotation vector (axis × angle).
    """

    poses: np.ndarray      # shape (n_lig, 6)
    velocity: np.ndarray   # shape (n_lig, 6)
    pbest_poses: np.ndarray
    pbest_score: float = float("inf")


def _pso_fitness(
    poses: np.ndarray,
    base_pos_list: List[np.ndarray],
    base_coms: List[np.ndarray],
    vdw_radii_list: List[np.ndarray],
    anchors: Sequence[AnchorSpec],
    target_map: Dict[int, np.ndarray],
    w_anchor: float = 4.0,
    w_clash: float = 6.0,
) -> float:
    """Rigid-body PSO fitness = anchor distance² + inter-ligand VdW clash penalty.

    Clash cutoffs use 75 % of the sum of Van-der-Waals radii so that real
    steric overlaps (e.g. H···C < 2.2 Å) are strongly penalised.
    """
    n_lig = len(base_pos_list)
    world_pos = [
        _apply_rigid_pose(base_pos_list[k], base_coms[k], poses[k, :3], poses[k, 3:])
        for k in range(n_lig)
    ]

    score = 0.0
    for anchor_idx, anchor in enumerate(anchors):
        d = world_pos[anchor.ligand_idx][anchor.donor_local_idx] - target_map[anchor_idx]
        score += w_anchor * float(np.dot(d, d))

    for i in range(n_lig):
        for j in range(i + 1, n_lig):
            diff_ij = world_pos[i][:, None, :] - world_pos[j][None, :, :]  # (Ni,Nj,3)
            dists = np.linalg.norm(diff_ij, axis=2)                         # (Ni,Nj)
            # 75 % of VdW contact: C-C → 2.55 Å, H-C → 2.18 Å, H-H → 1.80 Å
            cutoffs = 0.75 * (vdw_radii_list[i][:, None] + vdw_radii_list[j][None, :])
            score += w_clash * float(np.maximum(0.0, cutoffs - dists).sum())

    return score


def _run_pso(
    ligands: Sequence[LigandInfo],
    anchors: Sequence[AnchorSpec],
    target_map: Dict[int, np.ndarray],
    metal_positions: np.ndarray,
    n_particles: int = 20,
    n_iter: int = 150,
    c1: float = 1.50,
    c2: float = 1.50,
    record_trajectory: bool = False,
    trajectory_stride: int = 10,
) -> Tuple[np.ndarray, float, List[np.ndarray]]:
    """Particle Swarm Optimisation over rigid-body ligand poses.

    Each particle represents one complete assignment of 6-DOF rigid-body
    poses to all ligands simultaneously.  The fitness combines the squared
    distance of each donor atom to its coordination target plus a soft
    inter-ligand VdW clash penalty.  Ligand internal geometry is preserved
    exactly (rigid-body moves only), so intra-ligand topology cannot be
    corrupted and ligands cannot react with each other.

    Returns ``(best_poses, best_score, trajectory)`` where *best_poses* has
    shape ``(n_lig, 6)`` and *trajectory* is a list of gbest_poses snapshots
    (empty when *record_trajectory* is False).
    """
    n_lig = len(ligands)
    base_pos_list = [lig.atoms.get_positions().copy() for lig in ligands]
    base_coms = [np.mean(bp, axis=0) for bp in base_pos_list]
    # Use VdW radii (not covalent) so that non-bonded clashes are detected at
    # realistic distances (H-C clash kicks in at ~2.2 Å, not 1.0 Å).
    vdw_radii_list = [
        np.array([_get_vdw_radius(s) for s in lig.atoms.get_chemical_symbols()])
        for lig in ligands
    ]

    # Group anchor indices by ligand for initial-pose computation.
    anchors_by_lig: Dict[int, List[Tuple[int, AnchorSpec]]] = {}
    for a_idx, anchor in enumerate(anchors):
        anchors_by_lig.setdefault(anchor.ligand_idx, []).append((a_idx, anchor))

    # Ideal initial pose per ligand.
    # ≥2 donors → orthogonal Procrustes (aligns ALL donors to their targets).
    # 1 donor  → orient donor-to-COM along site direction.
    ideal_tvecs: List[np.ndarray] = []
    ideal_rvecs: List[np.ndarray] = []
    for k in range(n_lig):
        lig_a = anchors_by_lig.get(k, [])
        if not lig_a:
            ideal_tvecs.append(np.zeros(3))
            ideal_rvecs.append(np.zeros(3))
            continue

        com = base_coms[k]

        if len(lig_a) >= 2:
            # Orthogonal Procrustes: minimise sum |R(D_i-COM)+COM+t - T_i|²
            donors_body = np.array(
                [base_pos_list[k][a.donor_local_idx] for _, a in lig_a], dtype=float
            )  # (n, 3) positions in body frame
            targets_world = np.array(
                [target_map[a_idx] for a_idx, _ in lig_a], dtype=float
            )  # (n, 3) coordination targets
            # Centre both sets around their means.
            D_c = donors_body.mean(axis=0)
            T_c = targets_world.mean(axis=0)
            B = donors_body - D_c          # centred body donors
            T = targets_world - T_c        # centred targets
            # Optimal rotation via SVD of cross-covariance (body → world).
            # We need R s.t. B @ R.T ≈ T  →  cross-cov H = B.T @ T
            H = B.T @ T
            U, _S, Vt = np.linalg.svd(H)
            R = Vt.T @ U.T
            # Correct reflection (det should be +1).
            if float(np.linalg.det(R)) < 0.0:
                Vt[-1, :] *= -1.0
                R = Vt.T @ U.T
            rvec = _matrix_to_rvec(R)
            # tvec: places centroid of rotated donors at centroid of targets.
            # world(D_i) = (D_i - COM) @ R.T + COM + tvec
            # mean_world = (D_c - COM) @ R.T + COM + tvec = T_c
            rotated_D_c = R @ (D_c - com) + com
            tvec = T_c - rotated_D_c
        else:
            # Single donor: orient so donor→COM points along site direction.
            a_idx, anchor = lig_a[0]
            target = target_map[a_idx]
            metal_pos = metal_positions[anchor.metal_idx]
            donor_base = base_pos_list[k][anchor.donor_local_idx]
            site_dir = _normalize(target - metal_pos)
            donor_to_com = _normalize(com - donor_base)
            R = _rotation_matrix_from_vectors(donor_to_com, site_dir)
            rvec = _matrix_to_rvec(R)
            rotated_donor = R @ (donor_base - com) + com
            tvec = target - rotated_donor

        ideal_tvecs.append(tvec)
        ideal_rvecs.append(rvec)

    # Unit-vector pointing FROM metal TOWARD each ligand's coordination centroid.
    # Used to place particles "beyond" the target so they approach from afar.
    approach_dirs: List[np.ndarray] = []
    for k in range(n_lig):
        lig_a = anchors_by_lig.get(k, [])
        if not lig_a:
            approach_dirs.append(np.array([1.0, 0.0, 0.0]))
            continue
        targets_world = np.array([target_map[a_idx] for a_idx, _ in lig_a], dtype=float)
        T_c = targets_world.mean(axis=0)
        metal_pos_k = metal_positions[lig_a[0][1].metal_idx]
        approach_dirs.append(_normalize(T_c - metal_pos_k))

    # v_max is scaled to allow particles to traverse the full approach distance
    # (starting distance + ideal-placement distance) in a reasonable number of steps.
    START_FAR_DIST = 8.0  # Å beyond the coordination target (away from metal)
    max_disp = max(
        (float(np.linalg.norm(ideal_tvecs[k])) for k in range(n_lig)),
        default=2.0,
    )
    v_max = max(1.5, 0.40 * (max_disp + START_FAR_DIST))

    rng = np.random.default_rng()

    def _make_particle(far_dist: float, noise_lateral: float, noise_r: float) -> SwarmParticle:
        """Create one particle with ligands placed *far_dist* Å beyond their targets.

        All particles start far from the metal so the trajectory shows a genuine
        collective approach from afar toward the coordination sites.
        """
        poses = np.zeros((n_lig, 6))
        for k in range(n_lig):
            # Place ligand beyond coordination site along approach direction;
            # add lateral scatter so particles don't all start at identical spots.
            lateral = rng.standard_normal(3) * noise_lateral
            poses[k, :3] = ideal_tvecs[k] + approach_dirs[k] * far_dist + lateral
            # Keep Procrustes/single-donor orientation but add small rotation noise.
            poses[k, 3:] = ideal_rvecs[k] + rng.standard_normal(3) * noise_r
        vel = rng.standard_normal(poses.shape) * (v_max * 0.1)
        sc = _pso_fitness(poses, base_pos_list, base_coms, vdw_radii_list, anchors, target_map)
        return SwarmParticle(poses=poses, velocity=vel, pbest_poses=poses.copy(), pbest_score=sc)

    # All particles start FAR from the metal; minimum 5 Å beyond the coordination
    # target so that gbest at t=0 is genuinely far and the trajectory shows the
    # collective approach from afar.  Never placed at the ideal (target) position.
    far_schedule   = [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 9.0, 8.0, 7.0, 6.0,
                       5.0, 10.0, 7.0, 6.0, 8.0, 9.0, 5.0, 7.0, 6.0, 8.0]
    lat_schedule   = [0.5,  0.8,  1.0, 1.2, 1.0, 0.8, 1.2, 1.0, 0.8, 0.6,
                       0.5,  1.0,  1.2, 0.8, 1.0, 1.2, 0.6, 0.8, 1.0, 0.5]
    rot_schedule   = [0.30, 0.40, 0.35, 0.50, 0.40, 0.30, 0.45, 0.35, 0.25, 0.40,
                       0.30, 0.50, 0.35, 0.25, 0.40, 0.45, 0.20, 0.35, 0.30, 0.40]
    particles: List[SwarmParticle] = []
    for i in range(n_particles):
        idx = i % len(far_schedule)
        particles.append(_make_particle(far_schedule[idx], lat_schedule[idx], rot_schedule[idx]))

    gbest_score = min(p.pbest_score for p in particles)
    gbest_poses = min(particles, key=lambda p: p.pbest_score).pbest_poses.copy()

    trajectory: List[np.ndarray] = []
    if record_trajectory:
        trajectory.append(gbest_poses.copy())

    w = 0.80
    w_decay = (0.35 / 0.80) ** (1.0 / max(1, n_iter))

    for _it in range(n_iter):
        # Evaluate and update personal bests.
        for p in particles:
            sc = _pso_fitness(p.poses, base_pos_list, base_coms, vdw_radii_list, anchors, target_map)
            if sc < p.pbest_score:
                p.pbest_score = sc
                p.pbest_poses = p.poses.copy()
            if sc < gbest_score:
                gbest_score = sc
                gbest_poses = p.poses.copy()

        if gbest_score < 1.0e-3:
            break  # well converged

        # Record swarm best every trajectory_stride iterations.
        if record_trajectory and (_it + 1) % trajectory_stride == 0:
            trajectory.append(gbest_poses.copy())

        # Velocity and position update.
        for p in particles:
            r1 = rng.random((n_lig, 6))
            r2 = rng.random((n_lig, 6))
            p.velocity = (
                w * p.velocity
                + c1 * r1 * (p.pbest_poses - p.poses)
                + c2 * r2 * (gbest_poses - p.poses)
            )
            p.velocity = np.clip(p.velocity, -v_max, v_max)
            p.poses += p.velocity

        w *= w_decay

    # Always record final frame.
    if record_trajectory:
        trajectory.append(gbest_poses.copy())

    return gbest_poses, gbest_score, trajectory


def _poses_to_candidate_geometry(
    ligands: Sequence[LigandInfo],
    metals: Sequence[MetalInfo],
    metal_positions: np.ndarray,
    poses: np.ndarray,
) -> Tuple[Atoms, List[int], List[int], List[int], List[Tuple[int, int]]]:
    """Assemble an :class:`~ase.Atoms` from PSO-optimised rigid-body poses."""
    base_pos_list = [lig.atoms.get_positions().copy() for lig in ligands]
    base_coms = [np.mean(bp, axis=0) for bp in base_pos_list]

    symbols: List[str] = []
    positions_out: List[np.ndarray] = []
    fragment_ids: List[int] = []
    metal_global_indices: List[int] = []
    ligand_starts: List[int] = []
    anchored_pairs: List[Tuple[int, int]] = []

    for i, metal in enumerate(metals):
        metal_global_indices.append(len(symbols))
        symbols.append(metal.symbol)
        positions_out.append(np.array(metal_positions[i], dtype=float))
        fragment_ids.append(i)

    for lig_idx, ligand in enumerate(ligands):
        world_pos = _apply_rigid_pose(
            base_pos_list[lig_idx], base_coms[lig_idx],
            poses[lig_idx, :3], poses[lig_idx, 3:],
        )
        lig_start = len(symbols)
        ligand_starts.append(lig_start)
        for sym, pos in zip(ligand.atoms.get_chemical_symbols(), world_pos):
            symbols.append(sym)
            positions_out.append(pos)
            fragment_ids.append(len(metals) + lig_idx)
        for donor_local_idx, metal_idx in ligand.donor_links:
            anchored_pairs.append((lig_start + donor_local_idx, metal_global_indices[metal_idx]))

    atoms = Atoms(symbols=symbols, positions=np.array(positions_out, dtype=float))
    return atoms, fragment_ids, metal_global_indices, ligand_starts, anchored_pairs


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
    use_swarm: bool = True,
    swarm_particles: int = 20,
    swarm_iterations: int = 150,
    swarm_refine_steps: int = 80,
    use_xtb: bool = False,
    xtb_binary: str = "xtb",
    xtb_charge: int = 0,
    xtb_multiplicity: int = 1,
    xtb_parallel: int = 4,
    xtb_gfn: int = 2,
    use_xtb_opt: bool = False,
    xtb_opt_level: str = "loose",
    xtb_handoff_dist: float = 3.5,
    write_swarm_trajectory: bool = False,
    trajectory_stride: int = 10,
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

    # Auto-detect total formal charge from SMILES for xTB calculations.
    # The user-supplied xtb_charge=0 (default) is overridden when SMILES
    # contains explicit charges (e.g. [Cd+2], deprotonated ligand atoms).
    _smiles_mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if _smiles_mol is not None:
        _auto_charge = sum(a.GetFormalCharge() for a in _smiles_mol.GetAtoms())
        if xtb_charge == 0 and _auto_charge != 0:
            logger.info(
                "Auto-detected total formal charge %+d from SMILES; "
                "using for xTB (override with --xtb-charge).",
                _auto_charge,
            )
            xtb_charge = _auto_charge
        elif xtb_charge != 0 and xtb_charge != _auto_charge:
            logger.warning(
                "Manual xtb_charge=%+d differs from SMILES formal charge %+d.",
                xtb_charge, _auto_charge,
            )

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
    best_by_tag: Dict[str, CandidateResult] = {}  # best per coordination topology
    all_candidates: List[CandidateResult] = []    # for optional xTB re-ranking

    for cand_idx, (tag, target_map) in enumerate(target_sets, start=1):
        swarm_used = False

        # ── PSO rigid-body placement ──────────────────────────────────────────
        if use_swarm:
            try:
                best_poses, pso_score, traj = _run_pso(
                    ligands=ligands,
                    anchors=anchors,
                    target_map=target_map,
                    metal_positions=metal_positions,
                    n_particles=swarm_particles,
                    n_iter=swarm_iterations,
                    record_trajectory=write_swarm_trajectory,
                    trajectory_stride=trajectory_stride,
                )
                logger.debug(
                    "Candidate %d PSO converged with score %.6f", cand_idx, pso_score
                )
                (
                    atoms,
                    fragment_ids,
                    metal_global_indices,
                    ligand_starts,
                    anchored_pairs,
                ) = _poses_to_candidate_geometry(
                    ligands=ligands,
                    metals=metals,
                    metal_positions=metal_positions,
                    poses=best_poses,
                )
                # Write PSO trajectory as multi-frame XYZ (one frame per stride).
                if write_swarm_trajectory and traj:
                    traj_path = builder_dir / f"swarm_traj_{cand_idx:03d}.xyz"
                    with open(traj_path, "w", encoding="utf-8") as traj_fh:
                        for frame_poses in traj:
                            frame_atoms, _, _, _, _ = _poses_to_candidate_geometry(
                                ligands=ligands,
                                metals=metals,
                                metal_positions=metal_positions,
                                poses=frame_poses,
                            )
                            ase_write(traj_fh, frame_atoms, format="xyz")
                swarm_used = True
            except Exception as exc:
                logger.warning(
                    "Candidate %d PSO failed (%s); falling back to BFGS.", cand_idx, exc
                )

        # ── Fallback: original geometry assembly ──────────────────────────────
        if not swarm_used:
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

        # ── BFGS refinement (always first — topology-preserving elastic FF) ─────
        # Bond spring terms keep intra-ligand topology intact; anchor terms pull
        # donors to coordination targets; repulsion prevents inter-ligand clashes.
        # This MUST run before xTB so that xTB starts from a clean geometry.
        refine_steps = swarm_refine_steps if swarm_used else max_steps
        bond_terms, anchor_terms, repulsion_terms = _prepare_force_field_terms(
            atoms=atoms,
            ligands=ligands,
            anchors=anchors,
            target_map=target_map,
            ligand_starts=ligand_starts,
            fragment_ids=fragment_ids,
            anchored_pairs=anchored_pairs,
            metal_global_indices=metal_global_indices,
        )

        atoms.calc = ElasticNudgeCalculator(
            bond_terms=bond_terms,
            anchor_terms=anchor_terms,
            repulsion_terms=repulsion_terms,
        )
        atoms.set_constraint(FixAtoms(indices=metal_global_indices))

        try:
            dyn = BFGS(atoms, logfile=None)
            dyn.run(fmax=fmax, steps=refine_steps)
            energy = float(atoms.get_potential_energy())
        except Exception as exc:
            logger.warning("Candidate %d BFGS refinement failed: %s", cand_idx, exc)
            continue

        # ── Topology integrity check after BFGS ──────────────────────────────
        topo_ok, topo_reason = _check_topology_integrity(
            atoms=atoms,
            fragment_ids=fragment_ids,
            ligands=ligands,
            ligand_starts=ligand_starts,
            metal_global_indices=metal_global_indices,
            anchored_pairs=anchored_pairs,
        )
        if not topo_ok:
            logger.warning(
                "Candidate %d rejected (topology violation): %s", cand_idx, topo_reason
            )
            summary_lines.append(f"{cand_idx}\t{tag}\tREJECTED\t{topo_reason}")
            continue

        # ── Optional xTB geometry polish (runs AFTER BFGS gives clean topology) ─
        # BFGS guarantees a valid starting geometry; xTB can then flex multidentate
        # ligands and let them deform to reach all coordination sites properly.
        if use_xtb_opt:
            _bfgs_pos = atoms.get_positions()
            _max_dd = 0.0
            for _ai, _anch in enumerate(anchors):
                _dg = ligand_starts[_anch.ligand_idx] + _anch.donor_local_idx
                _tgt = np.array(target_map[_ai], dtype=float)
                _max_dd = max(_max_dd, float(np.linalg.norm(_bfgs_pos[_dg] - _tgt)))
            if _max_dd < xtb_handoff_dist:
                _xtb_res = _xtb_optimize(
                    atoms=atoms,
                    charge=xtb_charge,
                    multiplicity=xtb_multiplicity,
                    xtb_binary=xtb_binary,
                    fixed_indices=list(metal_global_indices),
                    opt_level=xtb_opt_level,
                    gfn_level=xtb_gfn,
                )
                if _xtb_res is not None:
                    _xtb_atoms, _xtb_e = _xtb_res
                    # Verify xTB didn't break the topology.
                    _atoms_test = atoms.copy()
                    _atoms_test.set_positions(_xtb_atoms.get_positions())
                    _topo2_ok, _topo2_reason = _check_topology_integrity(
                        atoms=_atoms_test,
                        fragment_ids=fragment_ids,
                        ligands=ligands,
                        ligand_starts=ligand_starts,
                        metal_global_indices=metal_global_indices,
                        anchored_pairs=anchored_pairs,
                    )
                    if _topo2_ok:
                        atoms.set_positions(_xtb_atoms.get_positions())
                        energy = _xtb_e
                        logger.info(
                            "Candidate %d xTB polish OK, E = %.8f Ha", cand_idx, energy
                        )
                    else:
                        logger.warning(
                            "Candidate %d xTB polish broke topology (%s); keeping BFGS geometry.",
                            cand_idx, _topo2_reason,
                        )
                else:
                    logger.debug("Candidate %d xTB polish failed; keeping BFGS geometry.", cand_idx)
            else:
                logger.debug(
                    "Candidate %d: BFGS donors still %.2f Å from targets → skip xTB polish.",
                    cand_idx, _max_dd,
                )

        clash = _compute_clash_score(atoms, fragment_ids)
        summary_lines.append(f"{cand_idx}\t{tag}\t{energy:.8f}\t{clash:.8f}")

        if write_candidates:
            cand_file = builder_dir / f"candidate_{cand_idx:03d}.xyz"
            ase_write(cand_file, atoms, format="xyz")

        cand = CandidateResult(tag=tag, atoms=atoms.copy(), energy=energy, clash_score=clash)
        all_candidates.append(cand)
        # Preliminary ranking by elastic energy (overridden by xTB below if enabled).
        if best is None or (cand.clash_score, cand.energy) < (best.clash_score, best.energy):
            best = cand
        prev = best_by_tag.get(tag)
        if prev is None or (cand.clash_score, cand.energy) < (prev.clash_score, prev.energy):
            best_by_tag[tag] = cand

    (builder_dir / "candidate_summary.tsv").write_text(
        "\n".join(summary_lines) + "\n",
        encoding="utf-8",
    )

    if best is None:
        logger.error("No valid candidate could be optimized")
        return False

    # ── Optional xTB re-ranking ───────────────────────────────────────────────
    if use_xtb and all_candidates:
        resolved_xtb = shutil.which(xtb_binary) or xtb_binary
        if not os.path.isfile(resolved_xtb):
            logger.warning("xTB binary not found (%s); skipping re-ranking.", xtb_binary)
        else:
            ranked = _rank_candidates_by_xtb(
                candidates=all_candidates,
                total_charge=xtb_charge,
                multiplicity=xtb_multiplicity,
                xtb_binary=resolved_xtb,
                n_parallel=xtb_parallel,
                gfn_level=xtb_gfn,
            )
            # Re-build best / best_by_tag from xTB energies.
            best = None
            best_by_tag = {}
            xtb_lines = ["#candidate\ttag\txtb_energy_Ha\tclash"]
            for cand_idx_xtb, cand in enumerate(ranked, start=1):
                xtb_lines.append(
                    f"{cand_idx_xtb}\t{cand.tag}\t{cand.energy:.10f}\t{cand.clash_score:.8f}"
                )
                if best is None or (cand.clash_score, cand.energy) < (best.clash_score, best.energy):
                    best = cand
                prev = best_by_tag.get(cand.tag)
                if prev is None or (cand.clash_score, cand.energy) < (prev.clash_score, prev.energy):
                    best_by_tag[cand.tag] = cand
            (builder_dir / "xtb_ranking.tsv").write_text(
                "\n".join(xtb_lines) + "\n", encoding="utf-8"
            )
            logger.info("xTB re-ranking written to xtb_ranking.tsv")

    # ── Write outputs ─────────────────────────────────────────────────────────
    def _safe_tag(t: str) -> str:
        """Convert tag label to a safe filename stem."""
        return t.replace(":", "_").replace(" ", "").replace("|", "__")

    parent_dir = builder_dir.parent

    # One XYZ + start.txt per coordination topology in BOTH builder_dir AND parent_dir.
    # This ensures every distinct coordination geometry is available for further processing
    # (e.g. ORCA optimisation) without having to dig into builder2/.
    logger.info(
        "Saving best structure for each of the %d distinct coordination topology/topologies:",
        len(best_by_tag),
    )
    for tag_label, tag_cand in sorted(best_by_tag.items()):
        stem = _safe_tag(tag_label)
        # Inside builder2/ (archive)
        tag_file = builder_dir / f"best_{stem}.xyz"
        ase_write(tag_file, tag_cand.atoms, format="xyz")
        # In parent working dir (directly usable)
        xyz_name = f"build_complex2_{stem}.xyz"
        _write_xyz_and_start(
            tag_cand.atoms,
            parent_dir=parent_dir,
            xyz_name=xyz_name,
        )
        # Rename the generic start.txt to a topology-specific name.
        generic_start = parent_dir / "start.txt"
        topo_start = parent_dir / f"start_{stem}.txt"
        if generic_start.exists():
            generic_start.rename(topo_start)
        logger.info(
            "  [%s]  energy=%.6f  clash=%.6f  →  %s  /  %s",
            tag_label,
            tag_cand.energy,
            tag_cand.clash_score,
            xyz_name,
            topo_start.name,
        )

    # Overall best → canonical build_complex2.xyz + start.txt (backward compat).
    best_file = builder_dir / "best_complex2.xyz"
    ase_write(best_file, best.atoms, format="xyz")
    _write_xyz_and_start(best.atoms, parent_dir=parent_dir, xyz_name="build_complex2.xyz")

    logger.info(
        "Overall best: [%s]  energy=%.8f  clash=%.8f",
        best.tag, best.energy, best.clash_score,
    )
    logger.info("Wrote %s", best_file)
    logger.info("Wrote %s  (overall best)", parent_dir / "build_complex2.xyz")
    logger.info("Wrote %s  (overall best)", parent_dir / "start.txt")
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
        "--swarm",
        dest="use_swarm",
        action="store_true",
        default=True,
        help="Use PSO swarm optimizer for collision-free ligand placement (default: on)",
    )
    parser.add_argument(
        "--no-swarm",
        dest="use_swarm",
        action="store_false",
        help="Disable swarm optimizer; use BFGS only",
    )
    parser.add_argument(
        "--swarm-particles",
        type=int,
        default=20,
        help="Number of PSO particles (default: 20)",
    )
    parser.add_argument(
        "--swarm-iterations",
        type=int,
        default=150,
        help="Number of PSO iterations per candidate (default: 150)",
    )
    parser.add_argument(
        "--swarm-refine-steps",
        type=int,
        default=80,
        help="BFGS refinement steps after PSO (default: 80; 0 to skip)",
    )
    parser.add_argument(
        "--write-swarm-trajectory",
        action="store_true",
        default=False,
        help="Write multi-frame XYZ trajectory of the PSO swarm per candidate",
    )
    parser.add_argument(
        "--trajectory-stride",
        type=int,
        default=10,
        help="Record swarm trajectory every N iterations (default: 10)",
    )
    parser.add_argument(
        "--xtb-opt",
        dest="use_xtb_opt",
        action="store_true",
        default=False,
        help=(
            "Use xTB geometry optimization as the refinement step after PSO "
            "(instead of BFGS). Allows ligands to flex/deform internally; "
            "preserves topology via GFN2 bonding model (default: off)"
        ),
    )
    parser.add_argument(
        "--xtb-opt-level",
        default="loose",
        choices=["crude", "loose", "normal", "tight"],
        help="xTB optimization convergence level (default: loose)",
    )
    parser.add_argument(
        "--xtb-handoff-dist",
        type=float,
        default=3.5,
        help=(
            "Hand off to xTB opt only when all donors are within this distance "
            "[Å] of their coordination targets after PSO (default: 3.5)"
        ),
    )
    parser.add_argument(
        "--xtb",
        dest="use_xtb",
        action="store_true",
        default=False,
        help="Re-rank candidates with xTB single-points after PSO/BFGS (default: off)",
    )
    parser.add_argument(
        "--xtb-binary",
        default="xtb",
        help="Path or name of the xTB executable (default: xtb)",
    )
    parser.add_argument(
        "--xtb-charge",
        type=int,
        default=0,
        help="Total charge for xTB calculation (default: 0)",
    )
    parser.add_argument(
        "--xtb-multiplicity",
        type=int,
        default=1,
        help="Spin multiplicity for xTB (default: 1 = singlet)",
    )
    parser.add_argument(
        "--xtb-parallel",
        type=int,
        default=4,
        help="Number of parallel xTB jobs (default: 4)",
    )
    parser.add_argument(
        "--xtb-gfn",
        type=int,
        default=2,
        choices=[0, 1, 2],
        help="GFN level for xTB (default: 2)",
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
        use_swarm=bool(args.use_swarm),
        swarm_particles=max(2, int(args.swarm_particles)),
        swarm_iterations=max(10, int(args.swarm_iterations)),
        swarm_refine_steps=max(0, int(args.swarm_refine_steps)),
        use_xtb=bool(args.use_xtb),
        xtb_binary=str(args.xtb_binary),
        xtb_charge=int(args.xtb_charge),
        xtb_multiplicity=max(1, int(args.xtb_multiplicity)),
        xtb_parallel=max(1, int(args.xtb_parallel)),
        xtb_gfn=int(args.xtb_gfn),
        use_xtb_opt=bool(args.use_xtb_opt),
        xtb_opt_level=str(args.xtb_opt_level),
        xtb_handoff_dist=max(1.0, float(args.xtb_handoff_dist)),
        write_swarm_trajectory=bool(args.write_swarm_trajectory),
        trajectory_stride=max(1, int(args.trajectory_stride)),
    )
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
