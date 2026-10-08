"""Embedding strategies of the MANTA constructor: retries, the no-valence-check path, hydrogen geometry repair, clash-aware hydrogen placement, the manual metal embed and the unsanitised fallback.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
import threading
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.conformer_io import (
    _fix_zero_coord_hydrogens,
    _mol_to_xyz,
    _openbabel_generate_conformer_xyz,
    _xyz_to_rdkit_conformer,
)
from delfin.manta.converter_flags import (
    _EMBED_TIMEOUT,
    _delfin_env_int,
    _deterministic_embed_seed,
)
from delfin.manta.embed_timeout import (
    _embed_with_timeout,
    _make_random_embed_params,
)
from delfin.manta.geometry_quality import (
    _geometry_quality_score,
)
from delfin.manta.hapto_detect import (
    mol_from_smiles_rdkit,
)
from delfin.manta.metal_smiles import (
    _denormalize_metal_smiles,
    _normalize_metal_smiles,
)
from delfin.manta.ml_tables import (
    AllChem,
    Chem,
    OPENBABEL_AVAILABLE,
    RDKIT_AVAILABLE,
    STK_AVAILABLE,
    _METAL_SET,
    _PREFERRED_CN4_GEOMETRY,
    _get_ml_bond_length,
    _ml_bond_kind,
    _prefer_no_sanitize,
    stk,
)
from delfin.manta.openbabel_optimize import (
    _optimize_xyz_openbabel_safe,
)
from delfin.manta.stage_hooks import (
    _apply_f19_to_fallback_xyz,
)

logger = get_logger("delfin.smiles_converter")


def _try_multiple_strategies(
    smiles: str,
    output_path: Optional[str] = None,
    deterministic: bool = True,
) -> Tuple[Optional[str], Optional[str]]:
    """Try multiple parsing strategies for metal complexes.

    IMPORTANT: Does NOT convert bonds to dative - both neutral and charged
    SMILES notations should produce the SAME result (no extra H atoms).

    This function tries different approaches in order of preference:
    1. stk (if available) - often handles complex coordination well
    2. RDKit with original SMILES
    3. RDKit with normalized (charged) SMILES
    4. RDKit with denormalized (neutral) SMILES
    5. RDKit unsanitized fallback (no H addition)
    6. Open Babel 3D generation (only when ``deterministic`` is ``False``)

    Args:
        smiles: SMILES string to convert.
        output_path: Optional path to write the resulting XYZ file.
        deterministic: When ``True`` (default) the Open Babel strategy is
            skipped because Open Babel's ``make3D`` builder is seeded from
            the wall clock with no Python-side seed hook and would make the
            output non-reproducible.  A deterministic error is preferred over
            a non-reproducible geometry.

    Returns:
        Tuple of (xyz_content, error_message)
    """
    errors = []
    normalized = _normalize_metal_smiles(smiles)
    denormalized = _denormalize_metal_smiles(smiles)

    # Build list of SMILES variants to try (original + alternatives)
    smiles_variants = [smiles]
    if normalized and normalized not in smiles_variants:
        smiles_variants.append(normalized)
    if denormalized and denormalized not in smiles_variants:
        smiles_variants.append(denormalized)

    # Strategy 1: Try stk first (good for metal complexes)
    # stk handles H atoms correctly based on the SMILES notation.
    # stk.BuildingBlock internally runs ETKDG embedding which can hang
    # on complex metal ring systems → guard with timeout.
    if STK_AVAILABLE:
        for smi in smiles_variants:
            try:
                _stk_result = [None]

                def _stk_build(s=smi):
                    try:
                        bb = stk.BuildingBlock(s)
                        _stk_result[0] = bb.to_rdkit_mol()
                    except Exception:
                        _stk_result[0] = None

                _st = threading.Thread(target=_stk_build, daemon=True)
                _st.start()
                _st.join(timeout=_EMBED_TIMEOUT)
                if _st.is_alive():
                    errors.append(f"stk({smi[:30]}...): timed out")
                    continue
                mol = _stk_result[0]
                if mol is not None:
                    # Embed if needed
                    if mol.GetNumConformers() == 0:
                        params = AllChem.ETKDGv3()
                        params.randomSeed = 42
                        params.useRandomCoords = True
                        result = _embed_with_timeout(mol, params)
                        if result != 0:
                            # Fallback: permissive embed (no ETKDG knowledge)
                            result = _embed_with_timeout(
                                mol, _make_random_embed_params(42))
                        if result != 0:
                            continue

                    xyz_content = _mol_to_xyz(mol)
                    if output_path:
                        Path(output_path).write_text(xyz_content, encoding='utf-8')
                        logger.info(f"Converted SMILES to XYZ using stk: {output_path}")
                    return xyz_content, None
            except Exception as e:
                errors.append(f"stk({smi[:30]}...): {e}")

    # Strategy 2: RDKit with partial sanitize (no dative conversion)
    for smi in smiles_variants:
        try:
            mol, note = mol_from_smiles_rdkit(smi, allow_metal=True)
            if mol is not None:
                # Try to add H without modifying bond types
                try:
                    mol = Chem.AddHs(mol, addCoords=True)
                except Exception:
                    pass  # Continue without H if it fails

                params = AllChem.ETKDGv3()
                params.randomSeed = 42
                params.useRandomCoords = True
                result = _embed_with_timeout(mol, params)
                if result != 0:
                    # Fallback: permissive embed (no ETKDG knowledge)
                    result = _embed_with_timeout(
                        mol, _make_random_embed_params(42))
                if result == 0:
                    xyz_content = _mol_to_xyz(mol)
                    if output_path:
                        Path(output_path).write_text(xyz_content, encoding='utf-8')
                        logger.info(f"Converted SMILES to XYZ using RDKit: {output_path}")
                    return xyz_content, None
        except Exception as e:
            errors.append(f"RDKit({smi[:30]}...): {e}")

    # Strategy 3: Unsanitized fallback (preserves SMILES as-is)
    for smi in smiles_variants:
        xyz_content, err = _smiles_to_xyz_unsanitized_fallback(smi)
        if xyz_content:
            if output_path:
                Path(output_path).write_text(xyz_content, encoding='utf-8')
                logger.info(f"Converted SMILES to XYZ using unsanitized fallback: {output_path}")
            return xyz_content, None
        if err:
            errors.append(f"unsanitized({smi[:30]}...): {err}")

    # Strategy 4: Ultra-permissive (bypass all valence checks)
    for smi in smiles_variants:
        xyz_content, err = _smiles_to_xyz_no_valence_check(smi)
        if xyz_content:
            if output_path:
                Path(output_path).write_text(xyz_content, encoding='utf-8')
                logger.info(f"Converted SMILES to XYZ using no-valence-check fallback: {output_path}")
            return xyz_content, None
        if err:
            errors.append(f"no-valence({smi[:30]}...): {err}")

    # Strategy 5: Open Babel / Avogadro-like 3D generation.
    # Open Babel's make3D builder is wall-clock seeded with no Python seed
    # hook, so this strategy is non-reproducible.  Skip it in deterministic
    # mode — a deterministic error beats a geometry that changes every run
    # and poisons every downstream metric.
    if OPENBABEL_AVAILABLE and not deterministic:
        xyz_blocks, ob_err = _openbabel_generate_conformer_xyz(
            smiles, num_confs=1, deterministic=False
        )
        if xyz_blocks:
            xyz_content = xyz_blocks[0]
            if output_path:
                Path(output_path).write_text(xyz_content, encoding='utf-8')
                logger.info("Converted SMILES to XYZ using Open Babel: %s", output_path)
            return xyz_content, None
        if ob_err:
            errors.append(f"openbabel({smiles[:30]}...): {ob_err}")
    elif OPENBABEL_AVAILABLE:
        errors.append("openbabel: skipped (deterministic mode)")

    return None, f"All strategies failed: {'; '.join(errors)}"


def _smiles_to_xyz_no_valence_check(smiles: str) -> Tuple[Optional[str], Optional[str]]:
    """Ultra-permissive fallback: bypass all valence checks and geometry rules.

    This is for problematic metal complexes where RDKit's valence rules
    don't apply (unusual coordination numbers, complex ring systems, etc.).
    No H atoms are added - only what's in the SMILES.
    Uses minimal geometry knowledge to handle exotic structures.
    """
    if not RDKIT_AVAILABLE:
        return None, "RDKit is not installed"

    try:
        # Parse without sanitization
        try:
            params = Chem.SmilesParserParams()
            params.sanitize = False
            params.removeHs = False
            params.strictParsing = False
            mol = Chem.MolFromSmiles(smiles, params)
        except Exception:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)

        if mol is None:
            return None, "Failed to parse SMILES"

        # Disable all implicit H handling to avoid valence checks
        for atom in mol.GetAtoms():
            atom.SetNoImplicit(True)

        # Use minimal geometry requirements - this handles exotic metal complexes
        # that don't follow standard geometry rules
        embed_params = AllChem.EmbedParameters()
        embed_params.randomSeed = 42
        embed_params.useRandomCoords = True
        embed_params.ignoreSmoothingFailures = True
        embed_params.useExpTorsionAnglePrefs = False  # Don't enforce torsion angles
        embed_params.useBasicKnowledge = False  # Don't use standard geometry knowledge
        embed_params.enforceChirality = False

        result = _embed_with_timeout(mol, embed_params)

        if result != 0:
            # Try with different random seeds
            for seed in [123, 456, 789, 1000]:
                embed_params.randomSeed = seed
                result = _embed_with_timeout(mol, embed_params)
                if result == 0:
                    break

        if result != 0:
            return None, "Failed to generate 3D coordinates"

        # Try to add H atoms after embedding (coordinates will be placed)
        try:
            # Reset noImplicit for non-metal atoms to allow H calculation
            for atom in mol.GetAtoms():
                if atom.GetSymbol() not in _METAL_SET:
                    atom.SetNoImplicit(False)
            mol.UpdatePropertyCache(strict=False)
            mol = Chem.AddHs(mol, addCoords=True)
        except Exception:
            # If H addition fails, continue with what we have
            pass

        xyz_content = _mol_to_xyz(mol)
        return xyz_content, None

    except Exception as e:
        return None, f"No-valence-check error: {e}"


def _fix_h_geometry_via_smiles(xyz: str, smiles: str) -> str:
    """Apply ``_fix_h_geometry_universal`` using a SMILES-derived mol.

    Parses the SMILES, adds H, attempts to align the atom order with
    the XYZ (greedy element-aware), then runs the universal H-fix.
    Returns the original xyz unchanged if alignment fails.
    """
    if not RDKIT_AVAILABLE or not xyz or not smiles:
        return xyz
    try:
        # Parse the SMILES leniently (matches the converter's tolerance)
        mol = None
        for sanitize_flag in (True, False):
            try:
                if sanitize_flag:
                    mol = Chem.MolFromSmiles(smiles)
                else:
                    mol = Chem.MolFromSmiles(smiles, sanitize=False)
                if mol is not None:
                    break
            except Exception:
                continue
        if mol is None:
            return xyz
        try:
            mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass
        try:
            mol = Chem.AddHs(mol)
        except Exception:
            return xyz
        # Verify atom count matches (XYZ should already include H)
        n_atoms_xyz = sum(
            1 for line in xyz.strip().splitlines()
            if line.strip() and len(line.split()) >= 4
        )
        if n_atoms_xyz != mol.GetNumAtoms():
            return xyz
        return _fix_h_geometry_universal(xyz, mol)
    except Exception:
        return xyz


def _clash_aware_h_place_xyz(xyz: str, threshold: float = 1.0) -> str:
    """Iter-6 T1 — Per-H clash-aware re-rotation around its parent heavy atom.

    Ported from 81f8a1f's ``_build_multimetal_hapto_sequential`` Step 6b.
    For each H whose nearest-non-parent neighbour is closer than
    ``threshold`` (Angstrom), tries 5 trial rotations (60deg, 120deg,
    180deg, 240deg, 300deg) around an axis perpendicular to the H-parent
    bond, picking the rotation that maximises that nearest-neighbour
    distance.  Three passes (the H placement of one atom can free space
    for the next).

    Pure XYZ-text in/out -- no RDKit dependency.  The "parent" of an H is
    inferred as its closest non-H atom in the input geometry.

    Gated by ``DELFIN_H_CLASH_ROTATE`` (default 0).  Bit-exact passthrough
    when disabled or when the input has no H atoms.
    """
    if not xyz or not _delfin_env_int("DELFIN_H_CLASH_ROTATE", 0):
        return xyz
    try:
        import numpy as _np
    except Exception:
        return xyz
    try:
        lines_in = xyz.strip().splitlines()
        header_lines: List[str] = []
        coord_lines: List[str] = []
        for line in lines_in:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    float(parts[1]); float(parts[2]); float(parts[3])
                    coord_lines.append(line)
                    continue
                except ValueError:
                    pass
            header_lines.append(line)
        if not coord_lines:
            return xyz
        positions: List[List[float]] = []
        syms: List[str] = []
        for line in coord_lines:
            parts = line.split()
            positions.append([
                float(parts[1]), float(parts[2]), float(parts[3])
            ])
            syms.append(parts[0])
        n = len(positions)
        pos = _np.asarray(positions, dtype=float)
        h_indices = [i for i, s in enumerate(syms) if s == "H"]
        if not h_indices:
            return xyz
        # Parent = closest non-H atom for each H
        parent_of: Dict[int, int] = {}
        for i in h_indices:
            best_j = -1
            best_d = 1e9
            for j in range(n):
                if j == i or syms[j] == "H":
                    continue
                d = float(_np.linalg.norm(pos[i] - pos[j]))
                if d < best_d:
                    best_d = d
                    best_j = j
            if best_j >= 0:
                parent_of[i] = best_j
        angles_deg = [60.0, 120.0, 180.0, 240.0, 300.0]
        for _pass in range(3):
            moved_any = False
            for i in h_indices:
                p = parent_of.get(i, -1)
                if p < 0:
                    continue
                # Current nearest non-(self,parent) distance
                cur_d = 1e9
                for j in range(n):
                    if j == i or j == p:
                        continue
                    dd = float(_np.linalg.norm(pos[i] - pos[j]))
                    if dd < cur_d:
                        cur_d = dd
                if cur_d >= threshold:
                    continue
                bond_vec = pos[i] - pos[p]
                bond_len = float(_np.linalg.norm(bond_vec))
                if bond_len < 1e-6:
                    continue
                bond_dir = bond_vec / bond_len
                ref = _np.array([1.0, 0.0, 0.0])
                if abs(float(_np.dot(bond_dir, ref))) > 0.95:
                    ref = _np.array([0.0, 1.0, 0.0])
                axis = _np.cross(bond_dir, ref)
                axn = float(_np.linalg.norm(axis))
                if axn < 1e-6:
                    continue
                axis = axis / axn
                best_pos = pos[i].copy()
                best_d = cur_d
                for ang_deg in angles_deg:
                    a = float(ang_deg) * 3.141592653589793 / 180.0
                    ca = float(_np.cos(a))
                    sa = float(_np.sin(a))
                    rv = (bond_dir * ca
                          + _np.cross(axis, bond_dir) * sa
                          + axis * float(_np.dot(axis, bond_dir)) * (1.0 - ca))
                    rvn = float(_np.linalg.norm(rv))
                    if rvn < 1e-9:
                        continue
                    cand = pos[p] + (rv / rvn) * bond_len
                    cand_d = 1e9
                    for j in range(n):
                        if j == i or j == p:
                            continue
                        dd = float(_np.linalg.norm(cand - pos[j]))
                        if dd < cand_d:
                            cand_d = dd
                    if cand_d > best_d:
                        best_d = cand_d
                        best_pos = cand
                        if best_d >= threshold:
                            break
                if best_d > cur_d + 1e-6:
                    pos[i] = best_pos
                    moved_any = True
            if not moved_any:
                break
        out_coord_lines: List[str] = []
        for i in range(n):
            out_coord_lines.append(
                "%s %.6f %.6f %.6f" % (
                    syms[i], float(pos[i, 0]),
                    float(pos[i, 1]), float(pos[i, 2]),
                )
            )
        return "\n".join(header_lines + out_coord_lines) + "\n"
    except Exception:
        return xyz


def _fix_h_geometry_universal(xyz: str, mol) -> str:
    """Re-place H atoms via VSEPR rules (universal — neighbour count only).

    For every non-H atom with at least one H neighbour, replaces all H
    positions with ideal sp3-Td (109.5°) / sp2-trigonal (120°) / sp-
    linear (180°) positions derived from the heavy-atom positions.

    Why: OB UFF (used as final refinement on metal complexes) lacks
    proper transition-metal parameters and frequently places H atoms
    on orthogonal x/y/z axes (90°/180° H-C-H) — physically nonsensical
    but topology-intact, so it slips through the topology gate.  This
    function fixes it as a final post-process.

    Heavy-atom positions are NEVER touched.  Bond lengths use 1.09 Å
    (C-H), 1.01 Å (N-H), 0.96 Å (O-H), 1.34 Å (S-H), else 1.10 Å.
    """
    if not RDKIT_AVAILABLE or not xyz or mol is None:
        return xyz

    h_bond_len = {
        "C": 1.09, "N": 1.01, "O": 0.96, "S": 1.34, "P": 1.42,
        "B": 1.19, "Si": 1.48, "F": 0.92, "Cl": 1.27, "Br": 1.41,
    }

    try:
        # Parse xyz
        lines = [l for l in xyz.strip().splitlines() if l.strip()]
        positions: List[List[float]] = []
        syms: List[str] = []
        for line in lines:
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
                syms.append(parts[0])
            except ValueError:
                continue
        n = len(positions)
        if n != mol.GetNumAtoms():
            return xyz  # mismatch, give up

        # For each heavy atom with H neighbours, recompute H positions
        for atom in mol.GetAtoms():
            ci = atom.GetIdx()
            if syms[ci] == "H":
                continue
            h_nbrs = [
                n.GetIdx() for n in atom.GetNeighbors()
                if syms[n.GetIdx()] == "H"
            ]
            non_h_nbrs = [
                n.GetIdx() for n in atom.GetNeighbors()
                if syms[n.GetIdx()] != "H"
            ]
            if not h_nbrs:
                continue

            total = len(h_nbrs) + len(non_h_nbrs)
            if total < 2 or total > 4:
                continue

            # Skip if atom has unrecognized total nbr count
            if total >= 4:
                expected_angle = 109.47
            elif total == 3:
                expected_angle = 120.0
            else:
                expected_angle = 180.0

            cx, cy, cz = positions[ci]
            # Build axis from heavy nbrs back to centre
            ax = ay = az = 0.0
            for oi in non_h_nbrs:
                ox, oy, oz = positions[oi]
                ax += cx - ox
                ay += cy - oy
                az += cz - oz
            axis_mag = math.sqrt(ax * ax + ay * ay + az * az)
            if axis_mag < 1e-8 or not non_h_nbrs:
                # No reference axis — keep original H positions
                continue
            ax /= axis_mag; ay /= axis_mag; az /= axis_mag

            # Build orthonormal basis perpendicular to axis
            if abs(az) < 0.9:
                ux, uy, uz = -ay, ax, 0.0
            else:
                ux, uy, uz = 1.0, 0.0, 0.0
            dot = ux * ax + uy * ay + uz * az
            ux -= dot * ax; uy -= dot * ay; uz -= dot * az
            u_mag = math.sqrt(ux * ux + uy * uy + uz * uz)
            if u_mag < 1e-8:
                ux, uy, uz = 1.0, 0.0, 0.0; u_mag = 1.0
            ux /= u_mag; uy /= u_mag; uz /= u_mag
            vx = ay * uz - az * uy
            vy = az * ux - ax * uz
            vz = ax * uy - ay * ux

            cos_h = math.cos(math.pi - math.radians(expected_angle))
            sin_h = math.sin(math.pi - math.radians(expected_angle))
            n_h = len(h_nbrs)
            phase_step = 2 * math.pi / max(n_h, 1)
            bl = h_bond_len.get(syms[ci], 1.10)

            for k, hi in enumerate(h_nbrs):
                phase = k * phase_step
                cp, sp = math.cos(phase), math.sin(phase)
                perp_x = cp * ux + sp * vx
                perp_y = cp * uy + sp * vy
                perp_z = cp * uz + sp * vz
                dx = (cos_h * ax + sin_h * perp_x) * bl
                dy = (cos_h * ay + sin_h * perp_y) * bl
                dz = (cos_h * az + sin_h * perp_z) * bl
                positions[hi] = [cx + dx, cy + dy, cz + dz]

        # Reassemble xyz
        out_lines = []
        for i in range(n):
            x, y, z = positions[i]
            out_lines.append(f"{syms[i]:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
        return "\n".join(out_lines) + "\n"
    except Exception:
        return xyz


def _manual_metal_embed(smiles: str) -> Tuple[Optional[str], Optional[str]]:
    """Last-resort fallback: manually construct coordinates for metal complexes.

    Places the metal at the origin and coordinating atoms at idealized
    positions (octahedral for 6-coord, tetrahedral for 4-coord, etc.).
    Remaining ligand atoms are placed along bond directions using a simple
    BFS traversal.  The geometry is rough but usable as a GOAT/xTB input.

    For 4-coordinate metals, both tetrahedral and square-planar geometries are
    tried, OB UFF is applied to each, and the better-scoring geometry is kept.
    """
    if not RDKIT_AVAILABLE:
        return None, "RDKit is not installed"

    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None, "Failed to parse SMILES"

        # Suppress valence checks
        for atom in mol.GetAtoms():
            atom.SetNoImplicit(True)

        # Try to add H for non-metal atoms
        try:
            for atom in mol.GetAtoms():
                if atom.GetSymbol() not in _METAL_SET:
                    atom.SetNoImplicit(False)
            mol.UpdatePropertyCache(strict=False)
            mol = Chem.AddHs(mol)
        except Exception:
            pass

        # Find metal centers
        metal_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _METAL_SET]
        if not metal_indices:
            return None, "No metal found"

        # Idealized coordination vectors for common coordination numbers
        _COORD_VECTORS_BASE = {
            2: [(1, 0, 0), (-1, 0, 0)],
            3: [(1, 0, 0), (-0.5, 0.866, 0), (-0.5, -0.866, 0)],
            4: [(1, 1, 1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1)],  # tetrahedral
            5: [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -0.5, 0.866), (0, -0.5, -0.866)],
            6: [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)],
            7: [(0, 0, 1), (0, 0, -1), (1, 0, 0), (0.309, 0.951, 0),
                (-0.809, 0.588, 0), (-0.809, -0.588, 0), (0.309, -0.951, 0)],
            8: [(1.414, 0, 0.8), (0, 1.414, 0.8), (-1.414, 0, 0.8), (0, -1.414, 0.8),
                (1, 1, -0.8), (-1, 1, -0.8), (-1, -1, -0.8), (1, -1, -0.8)],
        }

        n_atoms = mol.GetNumAtoms()

        def _build_xyz_for_4coord_vecs(override_4coord):
            """Build a raw XYZ string using given vectors for 4-coord metals."""
            local_coords = [(0.0, 0.0, 0.0)] * n_atoms
            local_placed: set = set()

            # Offset each metal so multi-metal complexes don't overlap
            _metal_offset = 4.0  # Angstrom separation between metal centers
            for metal_rank, mi in enumerate(metal_indices):
                mx = _metal_offset * metal_rank
                local_coords[mi] = (mx, 0.0, 0.0)
                local_placed.add(mi)

                neighbors = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(mi).GetNeighbors()]
                n_coord = len(neighbors)
                if n_coord == 4:
                    vectors = override_4coord
                else:
                    vectors = _COORD_VECTORS_BASE.get(n_coord, _COORD_VECTORS_BASE.get(6, []))

                metal_sym = mol.GetAtomWithIdx(mi).GetSymbol()
                # For high-CN systems (>8), use Fibonacci sphere for even distribution
                if n_coord > 8:
                    golden_ratio = (1 + math.sqrt(5)) / 2
                    vectors = []
                    for fi in range(n_coord):
                        theta = math.acos(1 - 2 * (fi + 0.5) / n_coord)
                        phi = 2 * math.pi * fi / golden_ratio
                        vectors.append((math.sin(theta) * math.cos(phi),
                                        math.sin(theta) * math.sin(phi),
                                        math.cos(theta)))
                for i, nbr_idx in enumerate(neighbors):
                    donor_sym = mol.GetAtomWithIdx(nbr_idx).GetSymbol()
                    # TERMINAL M=E length (16.08.2026).  Default OFF -> byte-identical.
                    # ⚠ THIS site is the CN4 PATH (`_build_xyz_for_4coord_vecs`) and was
                    # the ONLY wiring in the first attempt -- which is why `me42` reported
                    # `REACH 0/24` on 16.08. and died with rc=3.  Oxo and nitrido systems are
                    # predominantly CN5-7.  The load-bearing wiring has since sat in the M-D SNAP
                    # (three sites, ~25180/25265/25471): that one snaps EVERY bonded
                    # M-D distance to its ideal length, independent of the coordination number,
                    # and does so BEFORE UFF.  This line stays as the CN4 special case.
                    # Without the switch or without a supplied band table, `kind="me"` yields
                    # the same value as before, because `_ml_me_band` then returns None.
                    _kind = "sigma"
                    if _delfin_env_int("DELFIN_FFFREE_ME_BOND_LEN", 0):
                        _kind = _ml_bond_kind(mol, mi, nbr_idx)
                    bond_len = _get_ml_bond_length(metal_sym, donor_sym, _kind)
                    if i < len(vectors):
                        vx, vy, vz = vectors[i]
                        mag = math.sqrt(vx**2 + vy**2 + vz**2)
                        if mag > 1e-8:
                            vx = vx/mag * bond_len
                            vy = vy/mag * bond_len
                            vz = vz/mag * bond_len
                    else:
                        angle = 2 * math.pi * i / n_coord
                        vx = bond_len * math.cos(angle)
                        vy = bond_len * math.sin(angle)
                        vz = 0.0
                    local_coords[nbr_idx] = (mx + vx, vy, vz)
                    local_placed.add(nbr_idx)

            # BFS to place remaining atoms — VSEPR-aware: distribute
            # unplaced neighbours according to hybridization (sp/sp2/sp3)
            # derived from total neighbour count.  Universal — no
            # SMILES/element-specific logic, only geometry.
            bond_len_default = 1.4
            _bfs_counter = [0]
            queue = list(local_placed)
            while queue:
                current = queue.pop(0)
                cx, cy, cz = local_coords[current]
                c_pos = (cx, cy, cz)
                atom = mol.GetAtomWithIdx(current)
                placed_nbrs = [
                    n.GetIdx() for n in atom.GetNeighbors()
                    if n.GetIdx() in local_placed
                ]
                unplaced_nbrs = [
                    n.GetIdx() for n in atom.GetNeighbors()
                    if n.GetIdx() not in local_placed
                ]
                if not unplaced_nbrs:
                    continue

                total_nbrs = len(placed_nbrs) + len(unplaced_nbrs)
                # Hybridization expectation by total neighbour count
                if total_nbrs >= 4:
                    expected_angle = 109.47
                elif total_nbrs == 3:
                    expected_angle = 120.0
                elif total_nbrs == 2:
                    expected_angle = 180.0
                else:
                    expected_angle = 120.0  # safe default
                cos_target = math.cos(math.radians(expected_angle))

                # Build "axis" = mean direction from placed nbrs back to centre.
                ax = ay = az = 0.0
                for oi in placed_nbrs:
                    ox, oy, oz = local_coords[oi]
                    ax += cx - ox
                    ay += cy - oy
                    az += cz - oz
                axis_mag = math.sqrt(ax * ax + ay * ay + az * az)
                if axis_mag < 1e-8 or not placed_nbrs:
                    # No reference axis — Fibonacci sphere fallback
                    _bfs_counter[0] += 1
                    golden = (1 + math.sqrt(5)) / 2
                    theta = math.acos(1 - 2 * (_bfs_counter[0] % 50 + 0.5) / 50)
                    phi = 2 * math.pi * _bfs_counter[0] / golden
                    ax = math.sin(theta) * math.cos(phi)
                    ay = math.sin(theta) * math.sin(phi)
                    az = math.cos(theta)
                    axis_mag = 1.0
                ax /= axis_mag; ay /= axis_mag; az /= axis_mag

                # Build orthonormal basis (u, v) perpendicular to axis
                if abs(az) < 0.9:
                    ux = -ay; uy = ax; uz = 0.0
                else:
                    ux = 1.0; uy = 0.0; uz = 0.0
                # Gram-Schmidt
                dot = ux * ax + uy * ay + uz * az
                ux -= dot * ax; uy -= dot * ay; uz -= dot * az
                u_mag = math.sqrt(ux * ux + uy * uy + uz * uz)
                if u_mag < 1e-8:
                    ux = 1.0; uy = 0.0; uz = 0.0; u_mag = 1.0
                ux /= u_mag; uy /= u_mag; uz /= u_mag
                vx = ay * uz - az * uy
                vy = az * ux - ax * uz
                vz = ax * uy - ay * ux

                # For each unplaced neighbour, place at angle expected_angle
                # from axis.  H-C-X angle = expected_angle gives:
                #   nbr_dir · axis = cos(expected_angle)
                # since axis points FROM placed-nbrs TO centre (i.e. away
                # from placed nbrs), this nbr_dir points OUT from centre
                # at expected_angle from the OUT-direction of axis-of-
                # placed.  The unplaced H's go OPPOSITE the placed nbr.
                # nbr_dir · axis = -cos(180-expected) = cos(expected)
                # Correct: H_dir = cos(180-expected)*axis + sin(180-expected)*perp
                # where 180-expected is the angle between nbr_dir and the
                # "away from centre toward placed nbr" direction.
                # Equivalently: H_dir = -cos(expected)*(-axis) + sin(expected)*perp
                # which simplifies to:
                cos_h = math.cos(math.pi - math.radians(expected_angle))  # =-cos_target
                sin_h = math.sin(math.pi - math.radians(expected_angle))
                n_unplaced = len(unplaced_nbrs)
                # Distribute unplaced nbrs evenly around axis (120° apart for sp3 3H,
                # 180° apart for sp2 2H, 360° for sp 1H).
                phase_step = 2 * math.pi / max(n_unplaced, 1)
                for k, nbr_idx in enumerate(unplaced_nbrs):
                    phase = k * phase_step
                    cp, sp = math.cos(phase), math.sin(phase)
                    perp_x = cp * ux + sp * vx
                    perp_y = cp * uy + sp * vy
                    perp_z = cp * uz + sp * vz
                    dx = cos_h * ax + sin_h * perp_x
                    dy = cos_h * ay + sin_h * perp_y
                    dz = cos_h * az + sin_h * perp_z
                    dx *= bond_len_default
                    dy *= bond_len_default
                    dz *= bond_len_default

                    # Clash check: rotate around axis (preserves angle)
                    nx, ny, nz = cx + dx, cy + dy, cz + dz
                    for _attempt in range(5):
                        too_close = False
                        for pi in local_placed:
                            px, py, pz = local_coords[pi]
                            dd = math.sqrt((nx - px) ** 2 + (ny - py) ** 2 + (nz - pz) ** 2)
                            if dd < 0.8:
                                too_close = True
                                break
                        if not too_close:
                            break
                        # Rotate around axis by 60° (preserves the
                        # cos_target angle, just shifts azimuth)
                        a60 = math.radians(60)
                        ca, sa = math.cos(a60), math.sin(a60)
                        # Rodrigues' formula for rotation around axis (ax,ay,az)
                        ddx = dx; ddy = dy; ddz = dz
                        dot_da = ddx * ax + ddy * ay + ddz * az
                        crx = ay * ddz - az * ddy
                        cry = az * ddx - ax * ddz
                        crz = ax * ddy - ay * ddx
                        dx = ddx * ca + crx * sa + ax * dot_da * (1 - ca)
                        dy = ddy * ca + cry * sa + ay * dot_da * (1 - ca)
                        dz = ddz * ca + crz * sa + az * dot_da * (1 - ca)
                        nx, ny, nz = cx + dx, cy + dy, cz + dz
                    local_coords[nbr_idx] = (nx, ny, nz)
                    local_placed.add(nbr_idx)
                    queue.append(nbr_idx)

            lines = []
            for i in range(n_atoms):
                atom = mol.GetAtomWithIdx(i)
                x, y, z = local_coords[i]
                lines.append(f"{atom.GetSymbol():4s} {x:12.6f} {y:12.6f} {z:12.6f}")
            return '\n'.join(lines) + '\n'

        # Check whether any metal is 4-coordinate (warrants dual geometry trial)
        has_4coord = any(
            len([nbr.GetIdx() for nbr in mol.GetAtomWithIdx(mi).GetNeighbors()]) == 4
            for mi in metal_indices
        )

        if has_4coord:
            # Try tetrahedral, square-planar, and (Iter-2 Subagent A) see-saw vectors.
            # Bit-exact HEAD when DELFIN_ALL_POLYHEDRA=0 (SS branch is gated off).
            tetra_vecs = [(1, 1, 1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1)]
            sq_vecs = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
            ss_vecs = [(0, 0, 1), (0, 0, -1), (1, 0, 0.2), (-1, 0, 0.2)]

            # Determine preferred geometry from metal identity
            metal_syms = [mol.GetAtomWithIdx(mi).GetSymbol() for mi in metal_indices]
            pref_geom = _PREFERRED_CN4_GEOMETRY.get(metal_syms[0], 'SQ')
            _all_poly_on = bool(_delfin_env_int("DELFIN_ALL_POLYHEDRA", 0))  # Iter-5: default 1→0

            xyz_tetra = _build_xyz_for_4coord_vecs(tetra_vecs)
            xyz_sq = _build_xyz_for_4coord_vecs(sq_vecs)
            xyz_ss = _build_xyz_for_4coord_vecs(ss_vecs) if _all_poly_on else None

            # OB UFF refinement for all candidates — through safe wrapper for
            # universal coordination constraints.
            xyz_tetra = _optimize_xyz_openbabel_safe(xyz_tetra, mol_template=mol)
            xyz_sq = _optimize_xyz_openbabel_safe(xyz_sq, mol_template=mol)
            if xyz_ss is not None:
                xyz_ss = _optimize_xyz_openbabel_safe(xyz_ss, mol_template=mol)

            # H-placement: trust OB UFF output above (environment-aware
            # energy minimisation).  Disabled `_fix_h_geometry_universal`
            # post-process — it produced 5x more H-clash violations by
            # snapping H to ideal VSEPR angles with arbitrary rotational
            # phase (methyl umbrella could point INTO nearby atoms).
            # See INSIGHTS_LOG 2026-04-29 ~12:00 UTC.

            # Score via temporary RDKit conformers
            try:
                mol_tmp = Chem.RWMol(mol)
                mol_tmp.RemoveAllConformers()
                conf_t = _xyz_to_rdkit_conformer(mol_tmp.GetMol(), xyz_tetra)
                conf_s = _xyz_to_rdkit_conformer(mol_tmp.GetMol(), xyz_sq)
                conf_ss = (
                    _xyz_to_rdkit_conformer(mol_tmp.GetMol(), xyz_ss)
                    if xyz_ss is not None else None
                )
                if conf_t is not None and conf_s is not None:
                    cid_t = mol_tmp.AddConformer(conf_t, assignId=True)
                    cid_s = mol_tmp.AddConformer(conf_s, assignId=True)
                    score_t = _geometry_quality_score(mol_tmp.GetMol(), cid_t)
                    score_s = _geometry_quality_score(mol_tmp.GetMol(), cid_s)
                    best_xyz = xyz_tetra if score_t <= score_s else xyz_sq
                    best_score = score_t if score_t <= score_s else score_s
                    # Iter-2 Subagent A: include SS in score-off when active.
                    if conf_ss is not None and xyz_ss is not None:
                        cid_ss = mol_tmp.AddConformer(conf_ss, assignId=True)
                        score_ss = _geometry_quality_score(mol_tmp.GetMol(), cid_ss)
                        if score_ss < best_score:
                            best_xyz = xyz_ss
                            best_score = score_ss
                    xyz_str = best_xyz
                else:
                    # Use metal-preferred geometry as fallback
                    xyz_str = xyz_tetra if pref_geom == 'TH' else xyz_sq
            except Exception:
                xyz_str = xyz_tetra if pref_geom == 'TH' else xyz_sq
        else:
            # Standard path for non-4-coord metals (use default vectors)
            xyz_str = _build_xyz_for_4coord_vecs(
                _COORD_VECTORS_BASE.get(4, [(1, 1, 1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1)])
            )
            xyz_str = _optimize_xyz_openbabel_safe(xyz_str, mol_template=mol)
            # H-fix disabled — see comment above.

        return xyz_str, None

    except Exception as e:
        return None, f"Manual metal embed error: {e}"


def _smiles_to_xyz_unsanitized_fallback(smiles: str) -> Tuple[Optional[str], Optional[str]]:
    """Last-resort fallback: embed without sanitization (handles valence/kekulize errors)."""
    if not RDKIT_AVAILABLE:
        return None, "RDKit is not installed"

    try:
        normalized = _normalize_metal_smiles(smiles)
        smiles_try = normalized or smiles
        try:
            params = Chem.SmilesParserParams()
            params.sanitize = False
            params.removeHs = False
            params.strictParsing = False
            mol = Chem.MolFromSmiles(smiles_try, params)
        except Exception:
            mol = Chem.MolFromSmiles(smiles_try, sanitize=False)
        if mol is None:
            return None, "Failed to parse SMILES (no sanitize)"

        # Save explicit H counts and original NoImplicit state, then zero out
        # H and set NoImplicit=True to prevent valence errors during embedding
        # (e.g. [CH+] with 3 bonds = valence 4 > permitted 3 for C+).
        saved_explicit_h = {}
        originally_no_implicit = set()
        for atom in mol.GetAtoms():
            if atom.GetNoImplicit():
                originally_no_implicit.add(atom.GetIdx())
            eh = atom.GetNumExplicitHs()
            if eh > 0:
                saved_explicit_h[atom.GetIdx()] = eh
                atom.SetNumExplicitHs(0)
            atom.SetNoImplicit(True)

        try:
            mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass

        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.useRandomCoords = True
        try:
            params.maxAttempts = 200
        except Exception:
            pass

        try:
            result = _embed_with_timeout(mol, params)
        except (TypeError, Exception):
            result = -1

        if result != 0:
            # Relaxed embedding for complex metal geometries (ferrocene,
            # multi-center Pd complexes) where ETKDG distance bounds fail.
            relaxed = AllChem.EmbedParameters()
            relaxed.randomSeed = 42
            relaxed.useRandomCoords = True
            relaxed.useExpTorsionAnglePrefs = False
            relaxed.useBasicKnowledge = False
            relaxed.enforceChirality = False
            try:
                relaxed.ignoreSmoothingFailures = True
            except Exception:
                pass
            result = _embed_with_timeout(mol, relaxed)

        if result != 0:
            return None, "Failed to generate 3D coordinates (unsanitized)"

        # Restore explicit H and add hydrogens with coordinates
        if not _prefer_no_sanitize(smiles):
            try:
                # Restore saved explicit H counts
                for idx, eh in saved_explicit_h.items():
                    mol.GetAtomWithIdx(idx).SetNumExplicitHs(eh)
                # Allow implicit H calculation only for organic subset atoms
                # (originally NoImplicit=False, e.g. bare C in phenyl rings).
                # Bracket atoms like [C], [CH+], [C@@] keep NoImplicit=True
                # to respect the SMILES-specified H count.
                for atom in mol.GetAtoms():
                    if (atom.GetIdx() not in originally_no_implicit
                            and atom.GetSymbol() not in _METAL_SET):
                        atom.SetNoImplicit(False)
                        atom.SetNumExplicitHs(0)
                mol.UpdatePropertyCache(strict=False)
                mol = Chem.AddHs(mol, addCoords=True)
                mol = _fix_zero_coord_hydrogens(mol)
            except Exception:
                pass

        xyz_content = _mol_to_xyz(mol)
        # H-fix disabled — see INSIGHTS_LOG 2026-04-29 ~12:00 UTC.
        # Stream-B Fix 2 (DELFIN_FIX_F19, default-OFF → byte-identical): repair
        # sp3-H tetrahedrality (AddHs places fallback CH3 H at 90/180°).
        xyz_content = _apply_f19_to_fallback_xyz(xyz_content, mol)
        return xyz_content, None
    except Exception as e:
        # Extra-permissive fallback for explicit valence errors
        if "Explicit valence" in str(e):
            try:
                normalized = _normalize_metal_smiles(smiles)
                smiles_try = normalized or smiles
                mol = Chem.MolFromSmiles(smiles_try, sanitize=False)
                if mol is None:
                    return None, "Failed to parse SMILES (no sanitize)"
                # Disable implicit H handling to avoid valence checks
                for atom in mol.GetAtoms():
                    atom.SetNoImplicit(True)
                    atom.SetNumExplicitHs(0)
                # Determinism: keep an explicit fixed seed on every branch.
                # The TypeError fallback must NOT drop randomSeed, otherwise
                # this explicit-valence-error path embeds non-reproducibly.
                _uns_seed = _deterministic_embed_seed(smiles)
                try:
                    result = AllChem.EmbedMolecule(
                        mol, useRandomCoords=True, randomSeed=_uns_seed
                    )
                except TypeError:
                    _uns_params = AllChem.EmbedParameters()
                    _uns_params.randomSeed = _uns_seed
                    _uns_params.useRandomCoords = True
                    result = AllChem.EmbedMolecule(mol, _uns_params)
                if result == 0:
                    xyz_content = _mol_to_xyz(mol)
                    # Stream-B Fix 2 (DELFIN_FIX_F19, default-OFF → byte-identical).
                    xyz_content = _apply_f19_to_fallback_xyz(xyz_content, mol)
                    return xyz_content, None
            except Exception:
                pass
        return None, f"Unsanitized fallback error: {e}"
