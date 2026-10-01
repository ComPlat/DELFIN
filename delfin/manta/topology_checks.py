"""Graph-signature topology checks, metal connectivity verification and the geometry predicates of the MANTA constructor's gates.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
import os
from collections import Counter
from typing import Dict, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.converter_flags import (
    DELFIN_COLLAPSE_BOND_MIN_SCALE,
    DELFIN_COLLAPSE_TETRA_VOL_MIN,
    DELFIN_SEVERE_DIST_MAX_ABS,
    DELFIN_SEVERE_DIST_MAX_SCALE,
    _delfin_env_int,
)
from delfin.manta.ml_tables import (
    Chem,
    OPENBABEL_AVAILABLE,
    RDKIT_AVAILABLE,
    _METAL_ATOMICNUMS,
    _METAL_SET,
    _get_ml_bond_length,
    _target_mc_dist,
    pybel,
)

logger = get_logger("delfin.smiles_converter")


def _count_xyz_clashes(xyz: str, threshold: float = 1.20) -> int:
    """Count non-bonded heavy/H atom-pair clashes in an XYZ block.

    Used by Iter-6 T3 (B4 rollback) and T1 (H-place clash recheck) as a
    cheap pure-Python proxy for "did this transformation introduce a
    physically nonsensical close contact?".  Counts pairs whose distance
    falls below ``threshold`` (Angstrom) but whose elements are NOT both
    metals (metal-metal short contacts are real bonds in dimers/clusters).
    """
    try:
        import numpy as _np
        positions: List[List[float]] = []
        syms: List[str] = []
        for line in xyz.strip().splitlines():
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                positions.append([
                    float(parts[1]), float(parts[2]), float(parts[3])
                ])
                syms.append(parts[0])
            except (ValueError, IndexError):
                continue
        if len(positions) < 2:
            return 0
        pos = _np.asarray(positions, dtype=float)
        n = pos.shape[0]
        diff = pos[:, None, :] - pos[None, :, :]
        dist = _np.linalg.norm(diff, axis=-1)
        cnt = 0
        for i in range(n):
            for j in range(i + 1, n):
                d = float(dist[i, j])
                if d < 1e-3 or d >= threshold:
                    continue
                if syms[i] in _METAL_SET and syms[j] in _METAL_SET:
                    continue
                cnt += 1
        return int(cnt)
    except Exception:
        return 0


def _xyz_to_canonical_smiles(xyz_delfin: str) -> Optional[str]:
    """Convert DELFIN-format XYZ to a canonical SMILES via OpenBabel.

    Metal–ligand bonds perceived by OB from interatomic distances are
    intentionally removed before canonicalisation so that the comparison
    reflects only the organic/ligand connectivity (bond perception for
    M–L is unreliable and distance-dependent).

    Returns the canonical SMILES string, or ``None`` on failure.
    """
    if not OPENBABEL_AVAILABLE:
        return None
    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        if not lines:
            return None
        std_xyz = f"{len(lines)}\n\n" + "\n".join(lines) + "\n"
        ob_mol = pybel.readstring('xyz', std_xyz)

        # Remove bonds involving metal atoms — OB perceives them from
        # distance alone which is geometry-dependent and unreliable.
        metal_ob_idxs = {a.OBAtom.GetIdx() for a in ob_mol.atoms
                         if a.atomicnum in _METAL_ATOMICNUMS}
        if metal_ob_idxs:
            raw = ob_mol.OBMol
            bonds_to_del = []
            for bond in pybel.ob.OBMolBondIter(raw):
                bi = bond.GetBeginAtomIdx()
                ei = bond.GetEndAtomIdx()
                if bi in metal_ob_idxs or ei in metal_ob_idxs:
                    bonds_to_del.append(bond)
            for bond in bonds_to_del:
                raw.DeleteBond(bond)

        # Canonical SMILES (fragment-aware — disconnected parts separated by '.')
        conv = pybel.ob.OBConversion()
        conv.SetOutFormat('can')
        raw_smi = conv.WriteString(ob_mol.OBMol).strip()
        if not raw_smi:
            return None
        # Normalise: sort dot-separated fragments for stable comparison
        parts = sorted(raw_smi.split('.'))
        return '.'.join(parts)
    except Exception:
        return None


def _roundtrip_ring_count_ok(xyz_delfin: str, original_smiles: str, tolerance: int = 2) -> bool:
    """Validate a 3D structure by round-tripping coordinates → SMILES via OpenBabel.

    Only organic rings (not containing metal atoms) are compared. Metal chelate
    and coordination rings are excluded because OB bond perception from XYZ is
    unreliable for metal-ligand bonds (e.g. Ir-N ~2.1 Å), which would cause
    false rejections for complexes like fac/mer-Ir(ppy)3.
    """
    if not OPENBABEL_AVAILABLE:
        return True
    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        if not lines:
            return True
        std_xyz = f"{len(lines)}\n\n" + "\n".join(lines) + "\n"
        rt_mol = pybel.readstring('xyz', std_xyz)
        # Count ORGANIC-ONLY rings from OB: exclude any ring that contains a
        # metal atom (OB may perceive M-L bonds from short atom-atom distances,
        # adding chelate rings to the SSSR that are absent in the SMILES graph).
        try:
            _metal_ob_idxs = {a.idx for a in rt_mol.atoms
                              if a.atomicnum in _METAL_ATOMICNUMS}
            rt_rings = sum(
                1 for ring in rt_mol.sssr
                if not any(i in _metal_ob_idxs for i in ring._path)
            )
        except Exception:
            # OB failed to kekulize or ring._path unavailable — be permissive.
            # Using total SSSR here would include chelate rings and cause false
            # rejections for aromatic metal complexes (e.g. Ir(ppy)3).
            return True
        # Organic-only ring count from original SMILES.
        # OB is used first (same engine as XYZ parsing → consistent ring perception).
        # RDKit SMILES parsing often fails for metal-cyclometallated SMILES where
        # all ring-closure bonds pass through the metal atom (e.g. Ir(ppy)3),
        # causing RDKit to report 0 organic rings and rejecting all conformers.
        orig_rings = None
        try:
            orig_mol_ob = pybel.readstring('smi', original_smiles)
            _metal_smi_idxs = {a.idx for a in orig_mol_ob.atoms
                               if a.atomicnum in _METAL_ATOMICNUMS}
            orig_rings = sum(
                1 for ring in orig_mol_ob.sssr
                if not any(i in _metal_smi_idxs for i in ring._path)
            )
        except Exception:
            pass
        if orig_rings is None and RDKIT_AVAILABLE:
            try:
                orig_mol_rd = Chem.MolFromSmiles(original_smiles, sanitize=False)
                if orig_mol_rd is not None:
                    try:
                        orig_mol_rd.UpdatePropertyCache(strict=False)
                    except Exception:
                        pass
                    metal_indices = {
                        a.GetIdx() for a in orig_mol_rd.GetAtoms()
                        if a.GetSymbol() in _METAL_SET
                    }
                    ring_info = orig_mol_rd.GetRingInfo()
                    orig_rings = sum(
                        1 for ring in ring_info.AtomRings()
                        if not any(idx in metal_indices for idx in ring)
                    )
            except Exception:
                pass
        if orig_rings is None:
            return True
        return abs(rt_rings - orig_rings) <= tolerance
    except Exception:
        return True


def _no_spurious_bonds(xyz_delfin: str, original_smiles: str) -> bool:
    """Return True if OB-perceived XYZ contains no clearly spurious homodiatomic bonds.

    Only checks a short whitelist of homodiatomic element pairs that are very
    unlikely to appear in typical coordination-chemistry SMILES but can be
    falsely perceived by OpenBabel from close atom distances in bad geometries:

        O-O  (peroxide artifact)
        F-F, Cl-Cl, Br-Br, I-I  (halogen-halogen artifact)

    N-N, P-P, S-S are intentionally NOT checked: these appear in legitimate
    ligands (hydrazine, phosphine dimers, disulfides) and also produce false
    positives for metal complexes where heteroatom donors from different
    ligands end up close in ETKDG-generated geometries.

    Returns True (permissive) on any error or if required libraries are absent.
    """
    # Homodiatomic pairs to check (frozenset of atomic number, atomic number)
    _CHECKED_HOMODIATOMIC = {
        frozenset([8, 8]),    # O-O
        frozenset([9, 9]),    # F-F
        frozenset([17, 17]),  # Cl-Cl
        frozenset([35, 35]),  # Br-Br
        frozenset([53, 53]),  # I-I
    }

    if not OPENBABEL_AVAILABLE or not RDKIT_AVAILABLE:
        return True
    try:
        # Homodiatomic pairs present in the original SMILES (e.g. real peroxides)
        orig_mol = Chem.MolFromSmiles(original_smiles, sanitize=False)
        if orig_mol is None:
            return True
        try:
            orig_mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass
        orig_homo: set = set()
        for bond in orig_mol.GetBonds():
            n1 = bond.GetBeginAtom().GetAtomicNum()
            n2 = bond.GetEndAtom().GetAtomicNum()
            pair = frozenset([n1, n2])
            if pair in _CHECKED_HOMODIATOMIC:
                orig_homo.add(pair)
        # For multi-metal complexes with bridging O atoms, O-O perception
        # by OB is expected (two oxo/hydroxo on neighbouring metals) — skip.
        metal_count = sum(1 for a in orig_mol.GetAtoms() if a.GetSymbol() in _METAL_SET)
        if metal_count >= 2:
            o_on_metal = 0
            for a in orig_mol.GetAtoms():
                if a.GetAtomicNum() == 8 and any(
                    n.GetSymbol() in _METAL_SET for n in a.GetNeighbors()
                ):
                    o_on_metal += 1
            if o_on_metal >= 2:
                orig_homo.add(frozenset([8, 8]))

        # Homodiatomic pairs OB perceives in XYZ
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        if not lines:
            return True
        std_xyz = f"{len(lines)}\n\n" + "\n".join(lines) + "\n"
        rt_mol = pybel.readstring('xyz', std_xyz)
        try:
            from openbabel import openbabel as _ob
            for bond in _ob.OBMolBondIter(rt_mol.OBMol):
                n1 = bond.GetBeginAtom().GetAtomicNum()
                n2 = bond.GetEndAtom().GetAtomicNum()
                pair = frozenset([n1, n2])
                if pair in _CHECKED_HOMODIATOMIC and pair not in orig_homo:
                    logger.debug("Spurious homodiatomic bond in XYZ: %s-%s", n1, n2)
                    return False
        except Exception:
            return True

        return True
    except Exception:
        return True


# ============================================================================
# Iter-8.1 multi-hapto safe-fallback (FPCFD per-class extras-filter)
# ============================================================================
#
# Per master_v3 full-pool (cb0ef52):
#   multi-hapto extras=26.61/fr (44% catastrophic) vs cf1d480 13.02 (14%)
#   hapto       extras= 8.84/fr ( 9% catastrophic) vs cf1d480  5.41 ( 4%)
# Forensics: results/iter8.1_multihapto_safefallback_design.md.  Two stacked
# failure modes — broken seed + OB-WRS rotor amplification.  This is a
# read-only post-filter (no geometry mutation): drop frames with extras > τ
# AND > 3×median; Best-of-K fallback (K = max(2, ⌈0.3·n⌉)) guarantees ≥2
# frames per SMILES.  Class-dispatched via _classify_complex_class (L3489).

# Per-class extra-heavy-bond thresholds (Iter-8.1 hybrid filter)
_ITER8_1_EXTRA_THRESHOLDS = {
    "multi_hapto": 25,   # full-pool baseline 26.61 → drop catastrophic tail
    "hapto":       20,   # full-pool baseline 8.84 → drop tail
    "multi_sigma": 9999, # already at <0.10 — no filter
    "sigma":       9999, # already at <0.16 — no filter
    "no_metal":    9999, # no metal, no extras — no filter
}


def _count_extra_heavy_bonds(xyz_delfin: str, original_smiles: str) -> int:
    """Fast detector-faithful extra-bond count for class-dispatched filtering.

    Mirrors quality_framework `find_topology_loss.compare_topology` extras
    computation: heavy-atom covalent-radius graph from XYZ minus expected
    SMILES bond multiset (organic + M-L combined).  Returns -1 if any error.

    Used by Iter-8.1 multi-hapto safe-fallback to identify and drop frames
    where ligand collapse produces ghost C-C / M-C / etc. bonds.
    """
    if not RDKIT_AVAILABLE:
        return -1
    try:
        import math
        from collections import Counter as _Counter
        mol = Chem.MolFromSmiles(original_smiles, sanitize=False)
        if mol is None:
            return -1
        try:
            mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass
        # Build SMILES heavy-pair multiset (sorted-tuple key)
        smi_pairs: _Counter = _Counter()
        for bond in mol.GetBonds():
            a = bond.GetBeginAtom()
            b = bond.GetEndAtom()
            if a.GetAtomicNum() <= 1 or b.GetAtomicNum() <= 1:
                continue
            smi_pairs[tuple(sorted([a.GetSymbol(), b.GetSymbol()]))] += 1
        # Build XYZ heavy-pair multiset via covalent-radius graph
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        atoms_xyz = []
        for line in lines:
            parts = line.split()
            if len(parts) < 4 or parts[0] == "H":
                continue
            try:
                atoms_xyz.append((parts[0], float(parts[1]),
                                  float(parts[2]), float(parts[3])))
            except ValueError:
                continue
        n = len(atoms_xyz)
        if n == 0:
            return -1
        xyz_pairs: _Counter = _Counter()
        # Covalent radii subset (Cordero 2008); fallback 1.50 for unknowns
        _COV = {"C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Cl": 1.02,
                "Br": 1.20, "I": 1.39, "P": 1.07, "S": 1.05, "B": 0.84,
                "Si": 1.11, "Se": 1.20, "As": 1.19, "Te": 1.38}
        # Organic vs metal-bond tolerance (matches detector defaults)
        _TOL_ORG = 0.45
        _TOL_METAL = 0.60
        for i in range(n):
            si, xi, yi, zi = atoms_xyz[i]
            ri = _COV.get(si, 1.50)
            mi = si in _METAL_SET
            for j in range(i + 1, n):
                sj, xj, yj, zj = atoms_xyz[j]
                rj = _COV.get(sj, 1.50)
                mj = sj in _METAL_SET
                tol = _TOL_METAL if (mi or mj) else _TOL_ORG
                d = math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2)
                if d <= (ri + rj + tol):
                    xyz_pairs[tuple(sorted([si, sj]))] += 1
        # n_extra = sum over keys of max(0, xyz - smi)
        n_extra = 0
        for key in set(xyz_pairs) | set(smi_pairs):
            d = xyz_pairs.get(key, 0) - smi_pairs.get(key, 0)
            if d > 0:
                n_extra += d
        return n_extra
    except Exception:
        return -1


def _organic_graph_signature(
    symbol_by_idx: Dict[int, str],
    adj: Dict[int, set],
    wl_rounds: int = 3,
) -> "frozenset | None":
    """Return a WL-like multiset signature for disconnected organic graphs.

    The signature is atom-order independent and captures connectivity pattern
    per component much stronger than plain element-count fragment signatures.
    Bond order is intentionally ignored (distance-perception ambiguity in XYZ);
    connectivity preservation is the primary target.
    """
    try:
        comp_mult: Dict[tuple, int] = {}
        visited: set = set()
        for start in sorted(adj.keys()):
            if start in visited:
                continue
            comp: List[int] = []
            stack = [start]
            while stack:
                node = stack.pop()
                if node in visited:
                    continue
                visited.add(node)
                comp.append(node)
                stack.extend(adj.get(node, set()) - visited)

            labels: Dict[int, str] = {
                i: str(symbol_by_idx.get(i, '?')) for i in comp
            }
            for _ in range(max(1, int(wl_rounds))):
                new_labels: Dict[int, str] = {}
                for i in comp:
                    neigh = sorted(
                        labels.get(j, '?') for j in adj.get(i, set()) if j in labels
                    )
                    new_labels[i] = labels[i] + "|" + ",".join(neigh)
                labels = new_labels

            comp_sig = (
                len(comp),
                tuple(sorted(symbol_by_idx.get(i, '?') for i in comp)),
                tuple(sorted(len(adj.get(i, set())) for i in comp)),
                tuple(sorted(labels[i] for i in comp)),
            )
            comp_mult[comp_sig] = comp_mult.get(comp_sig, 0) + 1

        return frozenset(comp_mult.items())
    except Exception:
        return None


def _component_stats_from_adj(adj: Dict[int, set]) -> Optional[Tuple[int, int, int]]:
    """Return ``(n_components, largest_component_size, n_nodes)`` for adjacency map."""
    try:
        if not adj:
            return 0, 0, 0
        visited: set = set()
        n_comp = 0
        largest = 0
        for start in adj.keys():
            if start in visited:
                continue
            n_comp += 1
            stack = [start]
            size = 0
            while stack:
                node = stack.pop()
                if node in visited:
                    continue
                visited.add(node)
                size += 1
                stack.extend(adj.get(node, set()) - visited)
            if size > largest:
                largest = size
        return n_comp, largest, len(adj)
    except Exception:
        return None


def _heavy_graph_edges_smiles(
    smiles: str,
) -> Optional[Tuple[Dict[int, str], set]]:
    """Return non-metal heavy symbols and bond-order-neutral edge set from SMILES."""
    if not RDKIT_AVAILABLE:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        try:
            mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass

        heavy_symbols: Dict[int, str] = {}
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() <= 1:
                continue
            if atom.GetSymbol() in _METAL_SET:
                continue
            heavy_symbols[atom.GetIdx()] = atom.GetSymbol()

        edges: set = set()
        for bond in mol.GetBonds():
            bi = bond.GetBeginAtomIdx()
            bj = bond.GetEndAtomIdx()
            if bi not in heavy_symbols or bj not in heavy_symbols:
                continue
            edges.add(tuple(sorted((int(bi), int(bj)))))
        return heavy_symbols, edges
    except Exception:
        return None


def _heavy_graph_edges_xyz(
    xyz_delfin: str,
    expected_heavy_symbols: Dict[int, str],
) -> Optional[set]:
    """Return bond-order-neutral non-metal heavy edge set perceived from XYZ."""
    if not OPENBABEL_AVAILABLE:
        return None
    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        if not lines:
            return None

        for idx, sym in expected_heavy_symbols.items():
            if idx >= len(lines):
                return None
            parts = lines[idx].split()
            if not parts or parts[0] != sym:
                return None

        std_xyz = f"{len(lines)}\n\n" + "\n".join(lines) + "\n"
        ob_mol = pybel.readstring('xyz', std_xyz).OBMol
        try:
            from openbabel import openbabel as _ob
        except ImportError:
            return None

        heavy_idx_set = set(int(i) for i in expected_heavy_symbols.keys())
        edges: set = set()
        for bond in _ob.OBMolBondIter(ob_mol):
            i1 = int(bond.GetBeginAtomIdx()) - 1
            i2 = int(bond.GetEndAtomIdx()) - 1
            if i1 not in heavy_idx_set or i2 not in heavy_idx_set:
                continue
            edges.add(tuple(sorted((i1, i2))))
        return edges
    except Exception:
        return None


def _heavy_graph_exact_match_ok(
    xyz_delfin: str,
    original_smiles: str,
) -> bool:
    """Return True if non-metal heavy connectivity matches exactly, ignoring bond order."""
    try:
        ref = _heavy_graph_edges_smiles(original_smiles)
        if ref is None:
            return True
        heavy_symbols, ref_edges = ref
        xyz_edges = _heavy_graph_edges_xyz(xyz_delfin, heavy_symbols)
        if xyz_edges is None:
            return True
        if xyz_edges == ref_edges:
            return True

        missing = sorted(ref_edges - xyz_edges)
        added = sorted(xyz_edges - ref_edges)
        if missing:
            logger.debug(
                "Heavy-graph mismatch: missing %d bond(s), first=%s",
                len(missing),
                missing[0],
            )
        if added:
            logger.debug(
                "Heavy-graph mismatch: added %d bond(s), first=%s",
                len(added),
                added[0],
            )
        return False
    except Exception:
        return True


def _shortest_path_length_excluding(
    adj: Dict[int, set],
    start: int,
    goal: int,
    blocked: int,
    max_depth: int = 6,
) -> Optional[int]:
    """Return shortest path length from ``start`` to ``goal`` while skipping ``blocked``."""
    from collections import deque
    if start == goal:
        return 0
    visited = {blocked, start}
    queue = deque([(start, 0)])
    while queue:
        node, depth = queue.popleft()
        if depth >= max_depth:
            continue
        for nbr in adj.get(node, set()):
            if nbr in visited:
                continue
            if nbr == goal:
                return depth + 1
            visited.add(nbr)
            queue.append((nbr, depth + 1))
    return None


def _cycle_size_signature_from_adj(
    adj: Dict[int, set],
    max_cycle_size: int = 8,
) -> Dict[int, Tuple[int, ...]]:
    """Return per-atom cycle-size signatures derived from the adjacency graph."""
    cycle_sizes: Dict[int, set] = {idx: set() for idx in adj}
    for node, neighbors in adj.items():
        nbrs = sorted(neighbors)
        for i in range(len(nbrs)):
            for j in range(i + 1, len(nbrs)):
                path_len = _shortest_path_length_excluding(
                    adj,
                    nbrs[i],
                    nbrs[j],
                    blocked=node,
                    max_depth=max(2, int(max_cycle_size) - 2),
                )
                if path_len is None:
                    continue
                cycle_size = path_len + 2
                if 3 <= cycle_size <= max_cycle_size:
                    cycle_sizes[node].add(int(cycle_size))
    return {
        idx: tuple(sorted(sizes))
        for idx, sizes in cycle_sizes.items()
    }


def _heavy_local_signature_multiset(
    symbol_by_idx: Dict[int, str],
    adj: Dict[int, set],
    wl_rounds: int = 2,
) -> Counter:
    """Return a multiset of local heavy-atom graph signatures."""
    cycle_sig = _cycle_size_signature_from_adj(adj)
    labels: Dict[int, str] = {
        idx: f"{symbol_by_idx.get(idx, '?')}|d{len(adj.get(idx, set()))}|c{','.join(map(str, cycle_sig.get(idx, ()) ))}"
        for idx in adj
    }
    for _ in range(max(1, int(wl_rounds))):
        new_labels: Dict[int, str] = {}
        for idx in adj:
            neigh = sorted(labels.get(nbr, '?') for nbr in adj.get(idx, set()))
            new_labels[idx] = labels[idx] + "|" + ";".join(neigh)
        labels = new_labels
    signatures = Counter()
    for idx in adj:
        signatures[
            (
                symbol_by_idx.get(idx, '?'),
                len(adj.get(idx, set())),
                cycle_sig.get(idx, ()),
                labels.get(idx, '?'),
            )
        ] += 1
    return signatures


def _heavy_local_signature_multiset_smiles(
    smiles: str,
) -> Optional[Counter]:
    """Return local heavy-atom environment multiset from SMILES."""
    if not RDKIT_AVAILABLE:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        try:
            mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass
        keep = {
            atom.GetIdx(): atom.GetSymbol()
            for atom in mol.GetAtoms()
            if atom.GetAtomicNum() > 1 and atom.GetSymbol() not in _METAL_SET
        }
        adj: Dict[int, set] = {idx: set() for idx in keep}
        for bond in mol.GetBonds():
            bi = bond.GetBeginAtomIdx()
            bj = bond.GetEndAtomIdx()
            if bi in adj and bj in adj:
                adj[bi].add(bj)
                adj[bj].add(bi)
        return _heavy_local_signature_multiset(keep, adj)
    except Exception:
        return None


def _heavy_local_signature_multiset_xyz(
    xyz_delfin: str,
    expected_heavy_symbols: Dict[int, str],
) -> Optional[Counter]:
    """Return local heavy-atom environment multiset from XYZ via OB perception."""
    xyz_edges = _heavy_graph_edges_xyz(xyz_delfin, expected_heavy_symbols)
    if xyz_edges is None:
        return None
    adj: Dict[int, set] = {idx: set() for idx in expected_heavy_symbols}
    for i1, i2 in xyz_edges:
        if i1 in adj and i2 in adj:
            adj[i1].add(i2)
            adj[i2].add(i1)
    return _heavy_local_signature_multiset(expected_heavy_symbols, adj)


def _heavy_local_signature_match_ok(
    xyz_delfin: str,
    original_smiles: str,
) -> bool:
    """Return True when local heavy-atom graph environments still match."""
    try:
        ref = _heavy_graph_edges_smiles(original_smiles)
        if ref is None:
            return True
        heavy_symbols, _ref_edges = ref
        sig_smiles = _heavy_local_signature_multiset_smiles(original_smiles)
        sig_xyz = _heavy_local_signature_multiset_xyz(xyz_delfin, heavy_symbols)
        if sig_smiles is None or sig_xyz is None:
            return True
        if sig_smiles == sig_xyz:
            return True
        logger.debug(
            "Heavy local-signature mismatch: smiles=%d xyz=%d",
            sum(sig_smiles.values()),
            sum(sig_xyz.values()),
        )
        return False
    except Exception:
        return True


def _nonmetal_fragment_ids(mol) -> Dict[int, int]:
    """Return connected-component ids for heavy non-metal atoms."""
    if not RDKIT_AVAILABLE or mol is None:
        return {}
    adj: Dict[int, set] = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
            continue
        adj[atom.GetIdx()] = set()
    for bond in mol.GetBonds():
        bi = bond.GetBeginAtomIdx()
        bj = bond.GetEndAtomIdx()
        if bi in adj and bj in adj:
            adj[bi].add(bj)
            adj[bj].add(bi)

    frag_id_by_atom: Dict[int, int] = {}
    next_frag_id = 0
    for start in sorted(adj):
        if start in frag_id_by_atom:
            continue
        stack = [start]
        while stack:
            node = stack.pop()
            if node in frag_id_by_atom:
                continue
            frag_id_by_atom[node] = next_frag_id
            stack.extend(adj.get(node, set()) - set(frag_id_by_atom))
        next_frag_id += 1
    return frag_id_by_atom


def _metal_aware_coordination_ok(
    mol,
    conf_id: int,
    hapto_groups: Optional[List[Tuple[int, List[int]]]] = None,
) -> bool:
    """Return True when metal-centered donor fragments match the intended model."""
    if not RDKIT_AVAILABLE or mol is None:
        return True
    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return True

    frag_id_by_atom = _nonmetal_fragment_ids(mol)
    by_metal: Dict[int, List[List[int]]] = {}
    hapto_atom_to_metal: Dict[int, int] = {}
    for metal_idx, group_atoms in (hapto_groups or []):
        by_metal.setdefault(int(metal_idx), []).append(list(group_atoms))
        for atom_idx in group_atoms:
            hapto_atom_to_metal[int(atom_idx)] = int(metal_idx)

    def _dist(i: int, j: int) -> float:
        pi = conf.GetAtomPosition(i)
        pj = conf.GetAtomPosition(j)
        return math.sqrt(
            (pi.x - pj.x) ** 2 + (pi.y - pj.y) ** 2 + (pi.z - pj.z) ** 2
        )

    for metal in mol.GetAtoms():
        if metal.GetSymbol() not in _METAL_SET:
            continue
        metal_idx = metal.GetIdx()
        metal_sym = metal.GetSymbol()
        metal_hapto_groups = by_metal.get(metal_idx, [])
        metal_hapto_atoms = {
            atom_idx for group_atoms in metal_hapto_groups for atom_idx in group_atoms
        }

        # Hapto consistency: each expected group should retain a plausible
        # metal-centroid distance after candidate selection.
        for group_atoms in metal_hapto_groups:
            if not group_atoms:
                continue
            pts = [conf.GetAtomPosition(atom_idx) for atom_idx in group_atoms]
            centroid = (
                sum(p.x for p in pts) / len(pts),
                sum(p.y for p in pts) / len(pts),
                sum(p.z for p in pts) / len(pts),
            )
            mpos = conf.GetAtomPosition(metal_idx)
            mc_dist = math.sqrt(
                (centroid[0] - mpos.x) ** 2
                + (centroid[1] - mpos.y) ** 2
                + (centroid[2] - mpos.z) ** 2
            )
            target_mc = _target_mc_dist(metal_sym, len(group_atoms))
            if abs(mc_dist - target_mc) > 0.70:
                return False

        # Fragment-aware donor consistency for non-hapto donors.
        expected_frag_donors: Dict[int, Counter] = {}
        expected_donor_atom_indices: set = set()
        for nbr in metal.GetNeighbors():
            donor_idx = nbr.GetIdx()
            if nbr.GetAtomicNum() <= 1 or nbr.GetSymbol() in _METAL_SET:
                continue
            if donor_idx in metal_hapto_atoms:
                continue
            expected_donor_atom_indices.add(donor_idx)
            frag_id = frag_id_by_atom.get(donor_idx)
            if frag_id is None:
                continue
            expected_frag_donors.setdefault(frag_id, Counter())[nbr.GetSymbol()] += 1

        if expected_frag_donors:
            observed_frag_donors: Dict[int, Counter] = {}
            for atom in mol.GetAtoms():
                atom_idx = atom.GetIdx()
                if atom_idx == metal_idx:
                    continue
                if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
                    continue
                if atom_idx in metal_hapto_atoms:
                    continue
                frag_id = frag_id_by_atom.get(atom_idx)
                if frag_id is None:
                    continue
                dist = _dist(metal_idx, atom_idx)
                target = _get_ml_bond_length(metal_sym, atom.GetSymbol())
                # Count atoms that genuinely occupy the coordination shell.
                if dist <= target + 0.45:
                    observed_frag_donors.setdefault(frag_id, Counter())[atom.GetSymbol()] += 1

            for frag_id, expected_counter in expected_frag_donors.items():
                observed_counter = observed_frag_donors.get(frag_id, Counter())
                for donor_sym, expected_count in expected_counter.items():
                    if observed_counter.get(donor_sym, 0) < expected_count:
                        return False

            # For secondary/non-hapto metals, reject foreign fragments that
            # enter the first coordination shell unexpectedly.
            if not metal_hapto_groups:
                expected_frag_ids = set(expected_frag_donors.keys())
                for atom in mol.GetAtoms():
                    atom_idx = atom.GetIdx()
                    if atom_idx == metal_idx:
                        continue
                    if atom.GetAtomicNum() <= 1 or atom.GetSymbol() in _METAL_SET:
                        continue
                    frag_id = frag_id_by_atom.get(atom_idx)
                    if frag_id is None or frag_id in expected_frag_ids:
                        continue
                    dist = _dist(metal_idx, atom_idx)
                    intrusion_limit = _get_ml_bond_length(metal_sym, atom.GetSymbol()) + 0.15
                    if dist <= intrusion_limit:
                        return False

    return True


def _organic_fragment_signature(smiles: str) -> "frozenset | None":
    """Return an order-independent connectivity signature for organic fragments."""
    if not RDKIT_AVAILABLE:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        try:
            mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass

        keep = [
            a.GetIdx() for a in mol.GetAtoms()
            if a.GetAtomicNum() > 1 and a.GetSymbol() not in _METAL_SET
        ]
        symbol_by_idx: Dict[int, str] = {
            idx: mol.GetAtomWithIdx(idx).GetSymbol() for idx in keep
        }
        adj: Dict[int, set] = {idx: set() for idx in keep}
        for bond in mol.GetBonds():
            bi = bond.GetBeginAtomIdx()
            bj = bond.GetEndAtomIdx()
            if bi in adj and bj in adj:
                adj[bi].add(bj)
                adj[bj].add(bi)
        return _organic_graph_signature(symbol_by_idx, adj)
    except Exception:
        return None


def _heavy_component_stats_smiles(smiles: str) -> Optional[Tuple[int, int, int]]:
    """Return heavy-atom component stats for SMILES, including metal atoms."""
    if not RDKIT_AVAILABLE:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return None
        try:
            mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass
        keep = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
        adj: Dict[int, set] = {idx: set() for idx in keep}
        for bond in mol.GetBonds():
            bi = bond.GetBeginAtomIdx()
            bj = bond.GetEndAtomIdx()
            if bi in adj and bj in adj:
                adj[bi].add(bj)
                adj[bj].add(bi)
        return _component_stats_from_adj(adj)
    except Exception:
        return None


def _organic_fragment_signature_xyz(xyz_delfin: str) -> "frozenset | None":
    """Return organic connectivity signature for DELFIN XYZ (via OB perception)."""
    if not OPENBABEL_AVAILABLE:
        return None
    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        if not lines:
            return None
        std_xyz = f"{len(lines)}\n\n" + "\n".join(lines) + "\n"
        ob_mol = pybel.readstring('xyz', std_xyz).OBMol
        try:
            from openbabel import openbabel as _ob
        except ImportError:
            return None

        metal_ob_idx = {a.GetIdx() for a in _ob.OBMolAtomIter(ob_mol)
                        if a.GetAtomicNum() in _METAL_ATOMICNUMS}

        # Build heavy-atom adjacency excluding metals.
        n_atoms = ob_mol.NumAtoms()
        adj: dict = {i: set() for i in range(1, n_atoms + 1)
                     if ob_mol.GetAtom(i).GetAtomicNum() not in _METAL_ATOMICNUMS
                     and ob_mol.GetAtom(i).GetAtomicNum() not in (0, 1)}
        symbol_by_idx: Dict[int, str] = {
            i: (
                Chem.GetPeriodicTable().GetElementSymbol(ob_mol.GetAtom(i).GetAtomicNum())
                if RDKIT_AVAILABLE else str(ob_mol.GetAtom(i).GetAtomicNum())
            )
            for i in adj
        }

        for bond in _ob.OBMolBondIter(ob_mol):
            i1 = bond.GetBeginAtomIdx()
            i2 = bond.GetEndAtomIdx()
            if i1 in adj and i2 in adj:
                adj[i1].add(i2)
                adj[i2].add(i1)

        return _organic_graph_signature(symbol_by_idx, adj)
    except Exception:
        return None


def _heavy_component_stats_xyz(xyz_delfin: str) -> Optional[Tuple[int, int, int]]:
    """Return heavy-atom component stats for XYZ via Open Babel, incl. metals."""
    if not OPENBABEL_AVAILABLE:
        return None
    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        if not lines:
            return None
        std_xyz = f"{len(lines)}\n\n" + "\n".join(lines) + "\n"
        ob_mol = pybel.readstring('xyz', std_xyz).OBMol
        try:
            from openbabel import openbabel as _ob
        except ImportError:
            return None

        n_atoms = ob_mol.NumAtoms()
        adj: Dict[int, set] = {
            i: set()
            for i in range(1, n_atoms + 1)
            if ob_mol.GetAtom(i).GetAtomicNum() not in (0, 1)
        }
        for bond in _ob.OBMolBondIter(ob_mol):
            i1 = bond.GetBeginAtomIdx()
            i2 = bond.GetEndAtomIdx()
            if i1 in adj and i2 in adj:
                adj[i1].add(i2)
                adj[i2].add(i1)
        return _component_stats_from_adj(adj)
    except Exception:
        return None


def _global_heavy_connectivity_ok(
    xyz_delfin: str,
    original_smiles: str,
    max_extra_components: int = 1,
    min_largest_frac: float = 0.70,
) -> bool:
    """Reject severely fragmented heavy-atom graphs vs original SMILES.

    Complements ``_fragment_topology_ok``: the organic-only signature can miss
    cases where metal-ligand connectivity collapses while organic fragments stay
    internally intact.
    """
    orig_stats = _heavy_component_stats_smiles(original_smiles)
    xyz_stats = _heavy_component_stats_xyz(xyz_delfin)
    if orig_stats is None or xyz_stats is None:
        return True

    o_comp, o_largest, o_total = orig_stats
    x_comp, x_largest, _x_total = xyz_stats

    if x_comp > o_comp + int(max_extra_components):
        logger.debug(
            "Heavy-graph fragmentation: SMILES has %d component(s), XYZ has %d",
            o_comp, x_comp,
        )
        return False

    # If original is mostly one connected heavy graph, demand that XYZ keeps
    # a comparably large main component.
    if o_total > 0 and o_largest >= max(3, int(math.ceil(0.80 * o_total))):
        min_keep = max(2, int(math.floor(float(min_largest_frac) * float(o_largest))))
        if x_largest < min_keep:
            logger.debug(
                "Heavy-graph largest component collapsed: SMILES=%d, XYZ=%d (min=%d)",
                o_largest, x_largest, min_keep,
            )
            return False

    return True


def _fragment_topology_relaxed_fallback_ok(
    xyz_delfin: str,
    original_smiles: str,
) -> bool:
    """Allow only mildly fragmented fallback candidates.

    Used when strict fragment-topology checks reject all sampling conformers:
    keep only structures where the organic fragment signature is unchanged and
    heavy-atom connectivity degradation is limited.
    """
    try:
        orig_sig = _organic_fragment_signature(original_smiles)
        xyz_sig = _organic_fragment_signature_xyz(xyz_delfin)
        if orig_sig is None or xyz_sig is None or orig_sig != xyz_sig:
            return False
        return _global_heavy_connectivity_ok(
            xyz_delfin,
            original_smiles,
            max_extra_components=3,
            min_largest_frac=0.85,
        )
    except Exception:
        return False


def _metal_donor_distances_realistic(
    xyz_delfin: str,
    mol_template=None,
    min_abs_ml: float = 1.70,
    max_frac: float = 1.50,
) -> bool:
    """Reject structures with unphysical metal-donor distances.

    Two checks per metal atom in the XYZ:

    1. No heavy atom (any element, coordinating or not) closer than
       ``min_abs_ml`` to any metal — catches collapsed alt-binding
       artifacts where OB perception created false short bonds.
    2. When ``mol_template`` is available, each bonded non-metal donor
       must stay below ``max_frac`` × ``_get_ml_bond_length`` to flag
       broken connectivity.
    """
    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        atoms: List[Tuple[str, Tuple[float, float, float]]] = []
        for line in lines:
            parts = line.split()
            if len(parts) < 4:
                return True
            atoms.append(
                (parts[0], (float(parts[1]), float(parts[2]), float(parts[3])))
            )
        metals = [
            (sym, pos) for sym, pos in atoms if sym in _METAL_SET
        ]
        if not metals:
            return True

        min_abs_sq = float(min_abs_ml) * float(min_abs_ml)
        for m_sym, m_pos in metals:
            for a_sym, a_pos in atoms:
                if a_sym == m_sym and a_pos == m_pos:
                    continue
                if a_sym == 'H':
                    continue
                dx = m_pos[0] - a_pos[0]
                dy = m_pos[1] - a_pos[1]
                dz = m_pos[2] - a_pos[2]
                dsq = dx * dx + dy * dy + dz * dz
                if 0.0 < dsq < min_abs_sq:
                    return False

        if mol_template is not None and RDKIT_AVAILABLE:
            if len(atoms) != mol_template.GetNumAtoms():
                return True
            coords = [pos for _sym, pos in atoms]
            for atom in mol_template.GetAtoms():
                if atom.GetSymbol() not in _METAL_SET:
                    continue
                m_idx = atom.GetIdx()
                m_sym_t = atom.GetSymbol()
                for nbr in atom.GetNeighbors():
                    if nbr.GetSymbol() in _METAL_SET:
                        continue
                    if nbr.GetAtomicNum() <= 1:
                        continue
                    n_idx = nbr.GetIdx()
                    dx = coords[m_idx][0] - coords[n_idx][0]
                    dy = coords[m_idx][1] - coords[n_idx][1]
                    dz = coords[m_idx][2] - coords[n_idx][2]
                    d = math.sqrt(dx * dx + dy * dy + dz * dz)
                    ideal = float(_get_ml_bond_length(m_sym_t, nbr.GetSymbol()))
                    if ideal <= 0:
                        continue
                    if d > float(max_frac) * ideal:
                        return False
        return True
    except Exception:
        return False


def _verify_metal_connectivity(
    xyz_delfin: str,
    mol_template,
    max_donor_frac: float = 1.60,
    min_donor_frac: float = 0.75,
    min_nonddonor_frac: float = 0.80,
) -> bool:
    """Verify that the XYZ preserves the metal-donor connectivity from the template.

    Fundamental principles:
    1. Every defined donor must be within ``max_donor_frac × ideal`` of
       its metal (not drifted away).
    2. Every defined donor must be farther than ``min_donor_frac × ideal``
       from its metal (not collapsed onto it).
    3. Bridging donors (bonded to 2+ metals) only need to satisfy the
       distance window for *at least one* of their metals — they cannot
       be at ideal distance from all metals simultaneously.
    4. Non-bonded atoms must not have collapsed onto a metal
       (< min_nondonor_frac × ideal).
    """
    if not RDKIT_AVAILABLE or mol_template is None:
        return True
    try:
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        if len(lines) != mol_template.GetNumAtoms():
            return True
        coords: List[Tuple[float, float, float]] = []
        for line in lines:
            parts = line.split()
            if len(parts) < 4:
                return True
            coords.append((float(parts[1]), float(parts[2]), float(parts[3])))

        # Pre-compute which donor atoms are bridging (bonded to 2+ metals)
        # and collect all (metal_idx, donor_idx) pairs for deferred bridging check.
        bridging_pairs: Dict[int, List[Tuple[int, str, float]]] = {}
        for atom in mol_template.GetAtoms():
            if atom.GetSymbol() not in _METAL_SET:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() <= 1 or nbr.GetSymbol() in _METAL_SET:
                    continue
                n_metal_nbrs = sum(
                    1 for nn in nbr.GetNeighbors()
                    if nn.GetSymbol() in _METAL_SET
                )
                if n_metal_nbrs >= 2:
                    d_idx = nbr.GetIdx()
                    m_idx = atom.GetIdx()
                    m_sym = atom.GetSymbol()
                    d_sym = nbr.GetSymbol()
                    mx, my, mz = coords[m_idx]
                    dx, dy, dz = coords[d_idx]
                    d = math.sqrt((mx - dx) ** 2 + (my - dy) ** 2 + (mz - dz) ** 2)
                    ideal = float(_get_ml_bond_length(m_sym, d_sym))
                    if ideal > 0:
                        bridging_pairs.setdefault(d_idx, []).append((m_idx, m_sym, d / ideal))

        for atom in mol_template.GetAtoms():
            if atom.GetSymbol() not in _METAL_SET:
                continue
            m_idx = atom.GetIdx()
            m_sym = atom.GetSymbol()
            mx, my, mz = coords[m_idx]
            bonded_set = {nbr.GetIdx() for nbr in atom.GetNeighbors()}

            for d_idx in bonded_set:
                d_atom = mol_template.GetAtomWithIdx(d_idx)
                if d_atom.GetAtomicNum() <= 1:
                    continue
                d_sym = d_atom.GetSymbol()
                if d_sym in _METAL_SET:
                    continue
                dx, dy, dz = coords[d_idx]
                d = math.sqrt((mx - dx) ** 2 + (my - dy) ** 2 + (mz - dz) ** 2)
                ideal = float(_get_ml_bond_length(m_sym, d_sym))
                if ideal <= 0:
                    continue
                ratio = d / ideal

                # Bridging donors: defer to collective check below
                if d_idx in bridging_pairs:
                    continue

                # Terminal donor: must be within [min_donor_frac, max_donor_frac] × ideal
                if ratio > max_donor_frac or ratio < min_donor_frac:
                    return False

            # Non-donor collapse check — only for atoms ≥3 bonds away from
            # ANY metal.  Atoms within 2 bonds of any metal center are
            # naturally close in macrocyclic/chelating/bimetallic systems
            # (e.g. porphyrin C_alpha, bridging-ligand fragments).
            near_any_metal: set = set()
            for m_atom2 in mol_template.GetAtoms():
                if m_atom2.GetSymbol() not in _METAL_SET:
                    continue
                near_any_metal.add(m_atom2.GetIdx())
                for nbr1 in m_atom2.GetNeighbors():
                    near_any_metal.add(nbr1.GetIdx())
                    for nbr2 in nbr1.GetNeighbors():
                        near_any_metal.add(nbr2.GetIdx())
            for other in mol_template.GetAtoms():
                o_idx = other.GetIdx()
                if o_idx == m_idx or o_idx in near_any_metal:
                    continue
                if other.GetAtomicNum() <= 1:
                    continue
                if other.GetSymbol() in _METAL_SET:
                    continue
                ox, oy, oz = coords[o_idx]
                d = math.sqrt((mx - ox) ** 2 + (my - oy) ** 2 + (mz - oz) ** 2)
                ideal = float(_get_ml_bond_length(m_sym, other.GetSymbol()))
                if ideal <= 0:
                    continue
                if d < min_nonddonor_frac * ideal:
                    return False

        # Bridging donor check: each bridging donor must satisfy the distance
        # window for AT LEAST ONE of its metal partners.
        for d_idx, metal_entries in bridging_pairs.items():
            any_ok = False
            for _m_idx, _m_sym, ratio in metal_entries:
                if min_donor_frac <= ratio <= max_donor_frac:
                    any_ok = True
                    break
            if not any_ok:
                return False

        return True
    except Exception:
        return False


def _fragment_topology_ok(xyz_delfin: str, original_smiles: str) -> bool:
    """Return True if topology is consistent with the original SMILES."""
    orig_sig = _organic_fragment_signature(original_smiles)
    xyz_sig = _organic_fragment_signature_xyz(xyz_delfin)
    if orig_sig is not None and xyz_sig is not None and orig_sig != xyz_sig:
        def _total_frags(sig):
            return sum(mult for _, mult in sig)

        logger.debug(
            "Organic graph mismatch: SMILES has %d fragment(s), XYZ has %d",
            _total_frags(orig_sig), _total_frags(xyz_sig),
        )
        return False

    if not _global_heavy_connectivity_ok(xyz_delfin, original_smiles):
        return False

    return True


def _has_atom_clash(mol, conf_id: int, min_dist: float = 0.5) -> bool:
    """Return True if any pair of non-bonded atoms is closer than *min_dist* Å.

    Checks only heavy atoms (non-H) for efficiency.  Bonded atom pairs
    are excluded since their distance is governed by bond length.
    """
    conf = mol.GetConformer(conf_id)
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    bonded = set()
    for bond in mol.GetBonds():
        bonded.add((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        bonded.add((bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()))
    for i in range(len(heavy)):
        pi = conf.GetAtomPosition(heavy[i])
        for j in range(i + 1, len(heavy)):
            if (heavy[i], heavy[j]) in bonded:
                continue
            pj = conf.GetAtomPosition(heavy[j])
            dx = pi.x - pj.x
            dy = pi.y - pj.y
            dz = pi.z - pj.z
            if dx*dx + dy*dy + dz*dz < min_dist * min_dist:
                return True
    return False


def _has_unphysical_metal_nonbonded_contact(
    mol,
    conf_id: int,
    min_abs_dist: float = 1.10,
    rel_ml_scale: float = 0.55,
) -> bool:
    """Return True for unrealistically short non-bonded metal-heavy contacts.

    Some generated conformers place a non-coordinating heavy atom extremely
    close to a metal center (e.g. <1.0 A) while no bond exists in the graph.
    These geometries are physically implausible and should be rejected.
    """
    if not RDKIT_AVAILABLE:
        return False
    try:
        conf = mol.GetConformer(conf_id)
        for m_atom in mol.GetAtoms():
            if m_atom.GetSymbol() not in _METAL_SET:
                continue
            m_idx = m_atom.GetIdx()
            m_pos = conf.GetAtomPosition(m_idx)
            m_sym = m_atom.GetSymbol()
            for atom in mol.GetAtoms():
                a_idx = atom.GetIdx()
                if a_idx == m_idx:
                    continue
                if atom.GetAtomicNum() <= 1:
                    continue
                if atom.GetSymbol() in _METAL_SET:
                    continue
                if mol.GetBondBetweenAtoms(m_idx, a_idx) is not None:
                    continue

                a_pos = conf.GetAtomPosition(a_idx)
                d = math.sqrt(
                    (m_pos.x - a_pos.x) ** 2
                    + (m_pos.y - a_pos.y) ** 2
                    + (m_pos.z - a_pos.z) ** 2
                )
                expected_ml = _get_ml_bond_length(m_sym, atom.GetSymbol())
                min_d = max(float(min_abs_dist), float(rel_ml_scale) * float(expected_ml))
                if d < min_d:
                    return True
    except Exception:
        return False
    return False


def _has_unphysical_oco_geometry(
    mol,
    conf_id: int,
    min_oco_angle: float = 100.0,
    max_oco_angle: float = 145.0,
    min_co_dist: float = 1.05,
    max_co_dist: float = 1.45,
) -> bool:
    """Return True if a carboxyl-like O-C-O unit is unrealistically distorted.

    For carbon atoms bound to exactly two oxygen atoms (typical carboxyl/carbonyl
    motif), the O-C-O angle should be trigonal-planar-like (~120 deg), not
    linear like free CO2. Distances are also checked against broad covalent
    ranges to catch collapsed or stretched C-O bonds.
    """
    if not RDKIT_AVAILABLE:
        return False
    try:
        conf = mol.GetConformer(conf_id)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 6:
                continue
            if atom.GetSymbol() in _METAL_SET:
                continue

            o_neighbors = [
                n for n in atom.GetNeighbors()
                if n.GetAtomicNum() == 8 and n.GetSymbol() not in _METAL_SET
            ]
            if len(o_neighbors) != 2:
                continue

            c_idx = atom.GetIdx()
            c_pos = conf.GetAtomPosition(c_idx)
            o1_idx = o_neighbors[0].GetIdx()
            o2_idx = o_neighbors[1].GetIdx()
            o1_pos = conf.GetAtomPosition(o1_idx)
            o2_pos = conf.GetAtomPosition(o2_idx)

            d1 = math.sqrt(
                (c_pos.x - o1_pos.x) ** 2
                + (c_pos.y - o1_pos.y) ** 2
                + (c_pos.z - o1_pos.z) ** 2
            )
            d2 = math.sqrt(
                (c_pos.x - o2_pos.x) ** 2
                + (c_pos.y - o2_pos.y) ** 2
                + (c_pos.z - o2_pos.z) ** 2
            )
            if d1 < min_co_dist or d1 > max_co_dist or d2 < min_co_dist or d2 > max_co_dist:
                return True

            v1 = (o1_pos.x - c_pos.x, o1_pos.y - c_pos.y, o1_pos.z - c_pos.z)
            v2 = (o2_pos.x - c_pos.x, o2_pos.y - c_pos.y, o2_pos.z - c_pos.z)
            m1 = math.sqrt(v1[0] ** 2 + v1[1] ** 2 + v1[2] ** 2)
            m2 = math.sqrt(v2[0] ** 2 + v2[1] ** 2 + v2[2] ** 2)
            if m1 < 1e-10 or m2 < 1e-10:
                return True
            cos_a = max(-1.0, min(1.0, (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (m1 * m2)))
            angle = math.degrees(math.acos(cos_a))
            if angle < min_oco_angle or angle > max_oco_angle:
                return True
    except Exception:
        return False
    return False


def _flatten_sp2_atoms_xyz(xyz_delfin: str, mol_template) -> str:
    """Project every non-metal sp2 3-coordinate atom onto its neighbours' plane.

    UFF torsion constraints keep ring-junction atoms close to planar but do
    not fully enforce planarity.  This post-UFF polish moves any such atom
    exactly into the plane defined by its three heavy (non-metal) neighbours.
    The shift is purely geometric and acts on each atom independently, so it
    cannot introduce new topology, change ligand connectivity, or collapse
    donors onto the metal — the largest possible move for a 30° pyramidal
    deviation is below 0.25 Å for a typical 1.4 Å bond, well inside the
    graph gate's own bond-distance window.
    """
    if not RDKIT_AVAILABLE or mol_template is None or not xyz_delfin:
        return xyz_delfin
    try:
        import numpy as np
        lines = [l for l in xyz_delfin.splitlines() if l.strip()]
        if len(lines) != mol_template.GetNumAtoms():
            return xyz_delfin
        symbols: List[str] = []
        coords_list: List[List[float]] = []
        for line in lines:
            parts = line.split()
            if len(parts) < 4:
                return xyz_delfin
            symbols.append(parts[0])
            coords_list.append([float(parts[1]), float(parts[2]), float(parts[3])])
        coords = np.array(coords_list, dtype=float)

        def _is_sp2_graph(_atom):
            for _b in _atom.GetBonds():
                if (
                    _b.GetBondType() == Chem.BondType.AROMATIC
                    or _b.GetBondType() == Chem.BondType.DOUBLE
                    or _b.GetIsAromatic()
                    or _b.GetBondTypeAsDouble() >= 1.5
                ):
                    return True
            return False

        # ERDBEBEN E1' (default-off DELFIN_FFFREE_DONOR_PYRAMIDAL -> byte-identical): the bond-graph sp2
        # test above fires on HYPERVALENT PYRAMIDAL lone-pair donors that merely CARRY a double bond
        # (sulfoxide S=O, phosphine-oxide P=O, sulfone S) -> the flattener wrongly projects the correctly
        # built PYRAMIDAL donor (substituent-angle-sum ~298-318) onto its neighbour plane (->360 planar),
        # so the lone pair no longer points at the metal and M-D can't close (Ru-S 2.8-3.0 vs ~2.28).
        # Root fix, first-principles + universal: RDKit's CONJUGATION-AWARE hybridisation is authoritative --
        # an atom RDKit calls SP3 is PYRAMIDAL and must NOT be flattened, whatever its bond orders (RDKit
        # correctly gives sulfoxide/sulfone/thioether S, phosphine P, aqua/ether O, amine N -> SP3; amide/
        # aniline/pyridine/carbonyl -> SP2, still flattened).  Verified on TEXMIT: the template BUILDS the
        # sulfoxide-S pyramidal (298); this flatten step was breaking it.  Matches the landed donor-hyb
        # detector's expected-pyramidal rule.  Never element- or refcode-specific.
        _protect_pyramidal = _delfin_env_int("DELFIN_FFFREE_DONOR_PYRAMIDAL", 0)
        for atom in mol_template.GetAtoms():
            if atom.GetSymbol() in _METAL_SET:
                continue
            if atom.GetAtomicNum() <= 1:
                continue
            # sp2 from bond graph rather than flags so flatten acts
            # consistently on sanitised and unsanitised mols alike.
            if not _is_sp2_graph(atom):
                continue
            if _protect_pyramidal:
                # Scope to COORDINATING pyramidal donors (SP3 AND bonded to a metal).  RDKit's SP3 is the
                # pyramidal signal, but only a donor that COORDINATES the metal must keep its lone pair
                # pointing at it -- so only THOSE must not be flattened.  A NON-coordinating hypervalent S
                # (CITMUR's O-bound-DMSO sulfoxide-S: SP3 but bonded_metal=False) does not help the crystal
                # match; protecting it shifted CITMUR's best-valid frame worse and lost its CCDC isomer
                # (the E1' full:1000 A/B caught it).  Restricting to metal-bonded donors keeps TEXMIT's
                # Ru-bound sulfoxide-S win while dropping the CITMUR regression -- "fix the coordinating
                # donor without touching the rest".
                try:
                    if (str(atom.GetHybridization()) == "SP3"
                            and any(nb.GetSymbol() in _METAL_SET for nb in atom.GetNeighbors())):
                        continue
                except Exception:
                    pass
            # THE METAL IS A sigma PARTNER (DELFIN_FFFREE_DONOR_SIGMA_COUNT, default OFF).
            #
            # heavy_nbrs explicitly throws the metal away.  A donor with THREE heavy
            # neighbours PLUS metal thus really has FOUR sigma partners -- it is tetrahedral --
            # and is nevertheless projected here onto the plane of its three neighbours, i.e.
            # flattened.  That is the documented root in its most visible form:
            # hybridisation is decided on a METAL-FREE graph, and the donor counts one
            # partner too few.
            #
            # The protection above only kicks in when RDKit calls the donor SP3 -- and exactly
            # that it does not do when the metal was removed before typing: with three
            # neighbours the atom reads as SP2, and then flattening is by definition
            # right.  Hence here not the hybridisation OPINION, but COUNTING the
            # sigma partners, the metal counted in and the H too.
            #
            # The user saw it on VURMIE: "it is sp2 trigonal planar but has to coordinate
            # as sp3".  From now on it is measured by find_donor_hybridisation.
            _n_metal_nb = 0
            _n_h_nb = 0
            for _nb in atom.GetNeighbors():
                if _nb.GetSymbol() in _METAL_SET:
                    _n_metal_nb += 1
                elif _nb.GetAtomicNum() == 1:
                    _n_h_nb += 1
            heavy_nbrs = [
                n.GetIdx()
                for n in atom.GetNeighbors()
                if n.GetAtomicNum() > 1 and n.GetSymbol() not in _METAL_SET
            ]
            # THREE PARTNERS FROM PERIOD 3 ON ARE PYRAMIDAL -- that was the error of the first version.
            #
            # It protected only the case "four sigma partners -> tetrahedral" and measured affected=0.
            # The new axis find_donor_hybridisation then said where it really happens,
            # instead of me guessing it -- 14 of 24 conspicuous systems are pyramidal_flattened,
            # and the elements speak for themselves:
            #     NEYCUP Ti-Se(3sig) 37,42   QADFAC Pd-Se(3sig) 37,42   QADFEG Pt-Se(3sig) 37,42
            #     HIXBIA Pt-Te(3sig) 22,16   HIXBOG Pd-Te(3sig) 22,16
            # Three times resp. twice EXACTLY the same deviation across different metals -- a
            # construction rule, not noise.  Se and Te, with two heavy neighbours plus
            # metal, have exactly THREE partners and thus fell through the >=4 grid.
            #
            # The documented rule, now completely captured:
            #     >= 4 partners (metal counted in)       -> tetrahedral, not flat
            #     == 3 partners, period 2 (B,C,N,O,F)    -> planar, metal IN the plane: flat OK
            #     == 3 partners, period 3+ (P,S,Se,As..) -> PYRAMIDAL, not flat
            # Confirmed by the crystals (58770 of them): S|3 has p50 97.77 degrees -- clearly
            # below the 120 degrees of a plane -- while N|3 sits at 119.94.
            _sig_tot = len(heavy_nbrs) + _n_metal_nb + _n_h_nb
            _period3 = atom.GetSymbol() not in ("B", "C", "N", "O", "F", "H")
            if (_n_metal_nb
                    and os.environ.get("DELFIN_FFFREE_DONOR_SIGMA_COUNT", "0") == "1"
                    and (_sig_tot >= 4 or (_sig_tot == 3 and _period3))):
                continue          # tetrahedral or pyramidal -- in both cases NOT flat
            if len(heavy_nbrs) != 3:
                continue
            idx = atom.GetIdx()
            a, b, c = heavy_nbrs
            pa, pb, pc = coords[a], coords[b], coords[c]
            normal = np.cross(pb - pa, pc - pa)
            n_norm = float(np.linalg.norm(normal))
            if n_norm < 1e-9:
                continue
            normal = normal / n_norm
            centroid = (pa + pb + pc) / 3.0
            # Signed distance from X to plane(A,B,C)
            d = float(np.dot(coords[idx] - centroid, normal))
            coords[idx] = coords[idx] - d * normal

        out_lines: List[str] = []
        for sym, (x, y, z) in zip(symbols, coords):
            out_lines.append(f"{sym:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
        return "\n".join(out_lines) + "\n"
    except Exception as exc:
        logger.debug("_flatten_sp2_atoms_xyz failed: %s", exc)
        return xyz_delfin


def _has_pi_ring_nonplanarity(
    mol,
    conf_id: int,
    max_ring_rms: Optional[float] = None,
) -> bool:
    """Return True if unsaturated 5-7 membered rings are strongly non-planar.

    When ``max_ring_rms`` is ``None`` (the default going forward), the tolerance
    scales with ring size as ``0.25 * mean_in-ring_bond_length``.  A literal
    value can still be passed for call sites that want the historical absolute
    threshold.
    """
    if not RDKIT_AVAILABLE:
        return False
    try:
        import numpy as np
    except Exception:
        return False
    try:
        conf = mol.GetConformer(conf_id)
        try:
            Chem.GetSymmSSSR(mol)
        except Exception:
            pass
        ri = mol.GetRingInfo()
        if ri is None:
            return False

        for ring in ri.AtomRings():
            if len(ring) < 5 or len(ring) > 7:
                continue
            if any(
                mol.GetAtomWithIdx(i).GetAtomicNum() <= 1
                or mol.GetAtomWithIdx(i).GetSymbol() in _METAL_SET
                for i in ring
            ):
                continue

            unsat = 0
            valid_edges = 0
            for i in range(len(ring)):
                a = ring[i]
                b = ring[(i + 1) % len(ring)]
                bond = mol.GetBondBetweenAtoms(a, b)
                if bond is None:
                    continue
                valid_edges += 1
                if bond.GetBondType() != Chem.BondType.SINGLE or bond.GetIsAromatic():
                    unsat += 1
            if valid_edges < max(3, len(ring) - 1):
                continue
            if unsat < 2:
                continue

            pts = np.array([
                [
                    conf.GetAtomPosition(a).x,
                    conf.GetAtomPosition(a).y,
                    conf.GetAtomPosition(a).z,
                ]
                for a in ring
            ], dtype=float)
            ctr = pts.mean(axis=0)
            q = pts - ctr
            _u, _s, vh = np.linalg.svd(q, full_matrices=False)
            n = vh[-1]
            n_norm = float(np.linalg.norm(n))
            if n_norm < 1e-12:
                continue
            n = n / n_norm
            d = np.abs(q @ n)
            rms = float(np.sqrt(np.mean(d * d)))
            if max_ring_rms is None:
                edges = np.linalg.norm(
                    np.diff(np.vstack([pts, pts[:1]]), axis=0), axis=1
                )
                mean_bond = float(edges.mean()) if edges.size else 1.4
                tol = 0.25 * mean_bond
            else:
                tol = max_ring_rms
            if rms > tol:
                return True
    except Exception:
        return False
    return False


def _has_collapsed_sp3_centre(mol, conf, vol_min: float) -> bool:
    """Return True if any sp3 centre has PLANAR-COLLAPSED: a non-aromatic sp3 atom with >=2 heavy
    neighbours and >=4 bonded neighbours whose four nearest neighbours have a tetrahedron volume below
    ``vol_min`` A^3.  A real sp3 centre is tetrahedral (~1.5-2.5 A^3) even inside a flat ring; only a
    genuine collapse flattens the centre itself -- so this is silent-on-clean by construction (clean-CCDC
    floor 0.81).  Mirrors signal (A) of the WEDDELL ``find_planar_collapse`` eye axis; uses the RDKit
    bond graph (robust -- not distance-based adjacency, which a collapse would confuse)."""
    try:
        for a in mol.GetAtoms():
            if a.GetSymbol() in _METAL_SET or a.GetIsAromatic():
                continue
            if a.GetHybridization() != Chem.HybridizationType.SP3:
                continue
            nbrs = list(a.GetNeighbors())
            if len(nbrs) < 4 or sum(1 for nb in nbrs if nb.GetSymbol() != "H") < 2:
                continue
            pa = conf.GetAtomPosition(a.GetIdx())

            def _d2(nb):
                q = conf.GetAtomPosition(nb.GetIdx())
                return (q.x - pa.x) ** 2 + (q.y - pa.y) ** 2 + (q.z - pa.z) ** 2

            near = sorted(nbrs, key=_d2)[:4]
            p = [conf.GetAtomPosition(nb.GetIdx()) for nb in near]
            v1 = (p[1].x - p[0].x, p[1].y - p[0].y, p[1].z - p[0].z)
            v2 = (p[2].x - p[0].x, p[2].y - p[0].y, p[2].z - p[0].z)
            v3 = (p[3].x - p[0].x, p[3].y - p[0].y, p[3].z - p[0].z)
            cx = v1[1] * v2[2] - v1[2] * v2[1]
            cy = v1[2] * v2[0] - v1[0] * v2[2]
            cz = v1[0] * v2[1] - v1[1] * v2[0]
            vol = abs(cx * v3[0] + cy * v3[1] + cz * v3[2]) / 6.0
            if vol < vol_min:
                return True
    except Exception:
        return False
    return False


def _has_severe_covalent_distortion(
    mol,
    conf_id: int,
    max_abs_bond: Optional[float] = None,
    max_covalent_scale: Optional[float] = None,
) -> bool:
    """Return True if non-metal covalent bonds are unrealistically stretched.

    Thresholds default to the module-level ``DELFIN_SEVERE_DIST_MAX_ABS``
    and ``DELFIN_SEVERE_DIST_MAX_SCALE`` env-tunable constants so the
    gate can be relaxed or tightened without touching code.  Explicit
    overrides still win.

    When ``DELFIN_FFFREE_COLLAPSE_REJECT`` is set (default off ->
    byte-identical), the gate also rejects a PLANAR-COLLAPSED conformer via
    the symmetric lower bound it otherwise lacks: a covalent bond compressed
    below ``DELFIN_COLLAPSE_BOND_MIN_SCALE`` x the covalent-radius sum, or an
    sp3 centre flattened into a plane (see ``_has_collapsed_sp3_centre``).
    """
    if not RDKIT_AVAILABLE:
        return False
    if max_abs_bond is None:
        max_abs_bond = DELFIN_SEVERE_DIST_MAX_ABS
    if max_covalent_scale is None:
        max_covalent_scale = DELFIN_SEVERE_DIST_MAX_SCALE
    # read the on/off gate at CALL time (matches E1' DONOR_PYRAMIDAL) so the loop's --on, set after this
    # module is imported, actually activates it; the thresholds stay module-level (never toggled per side).
    _collapse_reject = _delfin_env_int("DELFIN_FFFREE_COLLAPSE_REJECT", 0) == 1
    try:
        conf = mol.GetConformer(conf_id)
        pt = Chem.GetPeriodicTable()
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetSymbol() in _METAL_SET or a2.GetSymbol() in _METAL_SET:
                continue

            p1 = conf.GetAtomPosition(a1.GetIdx())
            p2 = conf.GetAtomPosition(a2.GetIdx())
            d = math.sqrt(
                (p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2 + (p1.z - p2.z) ** 2
            )
            if d > max_abs_bond:
                return True

            try:
                rc1 = float(pt.GetRcovalent(a1.GetAtomicNum()))
                rc2 = float(pt.GetRcovalent(a2.GetAtomicNum()))
                if rc1 > 0 and rc2 > 0 and d > max_covalent_scale * (rc1 + rc2):
                    return True
                # collapse-reject (gated): a bond crushed below every real bond order = squashed frame
                if (_collapse_reject and rc1 > 0 and rc2 > 0
                        and d < DELFIN_COLLAPSE_BOND_MIN_SCALE * (rc1 + rc2)):
                    return True
            except Exception:
                pass
        # collapse-reject (gated): an sp3 centre flattened into a plane (distance-geometry degeneracy)
        if _collapse_reject and _has_collapsed_sp3_centre(mol, conf, DELFIN_COLLAPSE_TETRA_VOL_MIN):
            return True
    except Exception:
        return False
    return False
