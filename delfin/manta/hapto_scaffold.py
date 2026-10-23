"""Sphere-based constructive hapto scaffold: ansa bridges, shared eta rings, centroid bias, ring planarity, hapto geometry correction and BFS propagation in the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.converter_flags import (
    _class_conditional_flag,
    _delfin_env_int,
)
from delfin.manta.ml_tables import (
    Chem,
    Point3D,
    RDKIT_AVAILABLE,
    _COVALENT_RADII,
    _METAL_SET,
    _get_ml_bond_length,
    _target_mc_dist,
)

logger = get_logger("delfin.smiles_converter")


def _find_ansa_bridges(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
) -> List[Tuple[int, int, int]]:
    """Find bridging atoms (Si, C, Ge, etc.) connecting two hapto groups on the same metal.

    Returns list of (bridge_atom_idx, group_index_1, group_index_2).
    """
    if not RDKIT_AVAILABLE or mol is None or len(hapto_groups) < 2:
        return []

    # Build map: atom_idx -> list of (group_index, metal_idx)
    atom_to_groups: Dict[int, List[Tuple[int, int]]] = {}
    for gi, (metal_idx, grp) in enumerate(hapto_groups):
        for a_idx in grp:
            atom_to_groups.setdefault(a_idx, []).append((gi, metal_idx))

    bridges: List[Tuple[int, int, int]] = []
    seen_bridges: set = set()
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if sym in _METAL_SET or sym == 'H':
            continue
        a_idx = atom.GetIdx()
        if a_idx in atom_to_groups:
            continue  # skip hapto ring atoms themselves

        # Check if this atom connects to atoms from 2+ different groups on the same metal
        connected_groups: Dict[int, set] = {}  # metal_idx -> set of group indices
        for nbr in atom.GetNeighbors():
            ni = nbr.GetIdx()
            if ni in atom_to_groups:
                for gi, mi in atom_to_groups[ni]:
                    connected_groups.setdefault(mi, set()).add(gi)

        for _mi, group_set in connected_groups.items():
            if len(group_set) >= 2:
                glist = sorted(group_set)
                for i in range(len(glist)):
                    for j in range(i + 1, len(glist)):
                        key = (a_idx, glist[i], glist[j])
                        if key not in seen_bridges:
                            seen_bridges.add(key)
                            bridges.append(key)

    return bridges


def _detect_shared_ring_eta_groups(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
) -> List[Dict[str, object]]:
    """P4 Layer 2 — detect clusters of eta-groups fused into one macrocycle.

    Generalises :func:`_find_ansa_bridges` (which only sees SINGLE-ATOM
    bridges) to multi-atom bridge PATHS (e.g. a benzene ring connecting two
    eta(C=C) groups, as in MEWCIA / MIRSUE).  For each metal, BFS between
    every pair of its eta-groups through non-eta / non-metal atoms; a pair
    is *bridge-connected* when such a path exists (length <= ``max_path``).
    Pairwise bridge-connected eta-groups are unioned into clusters; a
    cluster of >= 2 groups that share one macrocyclic ring is returned.

    Returns a list of dicts, one per metal-cluster::

        {"metal": metal_idx,
         "groups": [g_idx, ...],            # sorted local group indices
         "bridge_atoms": frozenset(...),    # all path atoms between groups
         "n_groups": int}

    Graph/group-theory only — no SMILES regex, no refcode keying.
    Deterministic: groups & atoms sorted by index throughout.
    """
    if not RDKIT_AVAILABLE or mol is None or len(hapto_groups) < 2:
        return []
    from collections import deque

    # Map atom -> set of group indices it belongs to (per metal).
    atom_to_group: Dict[int, int] = {}
    group_metal: Dict[int, int] = {}
    group_atoms: Dict[int, set] = {}
    for gi, (metal_idx, grp) in enumerate(hapto_groups):
        group_metal[gi] = metal_idx
        group_atoms[gi] = set(grp)
        for a in grp:
            atom_to_group[a] = gi

    all_eta_atoms = set(atom_to_group.keys())
    max_path = 6  # benzene bridge = path length ~3-4; allow some slack

    # Per metal, find bridge-connected eta-group pairs + their path atoms.
    groups_by_metal: Dict[int, List[int]] = {}
    for gi, mi in group_metal.items():
        groups_by_metal.setdefault(mi, []).append(gi)

    def _bridge_path(src_grp: int, dst_grp: int) -> Optional[List[int]]:
        """Shortest BFS path of non-eta / non-metal atoms linking any atom
        of src_grp to any atom of dst_grp.  Returns the interior path atoms
        (excluding the eta endpoints), or None if no short bridge exists."""
        dst_set = group_atoms[dst_grp]
        # Seed from atoms directly bonded to a src-group atom (interior only).
        start_atoms: List[Tuple[int, List[int]]] = []
        seeded: set = set()
        for a in sorted(group_atoms[src_grp]):
            for nb in mol.GetAtomWithIdx(a).GetNeighbors():
                ni = nb.GetIdx()
                if (ni in all_eta_atoms or nb.GetSymbol() in _METAL_SET):
                    continue
                if ni in seeded:
                    continue
                seeded.add(ni)
                start_atoms.append((ni, [ni]))
        q = deque(start_atoms)
        visited = set(seeded)
        while q:
            cur, path = q.popleft()
            # Reached an atom bonded to the destination group?
            for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
                ni = nb.GetIdx()
                if ni in dst_set:
                    return path
            if len(path) >= max_path:
                continue
            for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
                ni = nb.GetIdx()
                if (ni in all_eta_atoms or nb.GetSymbol() in _METAL_SET
                        or ni in visited):
                    continue
                visited.add(ni)
                q.append((ni, path + [ni]))
        return None

    results: List[Dict[str, object]] = []
    for mi in sorted(groups_by_metal.keys()):
        glist = sorted(groups_by_metal[mi])
        if len(glist) < 2:
            continue
        # Build bridge graph among this metal's eta-groups.
        adj: Dict[int, set] = {g: set() for g in glist}
        pair_paths: Dict[Tuple[int, int], List[int]] = {}
        for ii in range(len(glist)):
            for jj in range(ii + 1, len(glist)):
                gi, gj = glist[ii], glist[jj]
                path = _bridge_path(gi, gj)
                if path is not None:
                    adj[gi].add(gj)
                    adj[gj].add(gi)
                    pair_paths[(gi, gj)] = path
        # Union-find / connected components over the bridge graph.
        seen_g: set = set()
        for g0 in glist:
            if g0 in seen_g:
                continue
            comp: List[int] = []
            stack = [g0]
            while stack:
                cur = stack.pop()
                if cur in seen_g:
                    continue
                seen_g.add(cur)
                comp.append(cur)
                for nb in sorted(adj[cur]):
                    if nb not in seen_g:
                        stack.append(nb)
            comp.sort()
            if len(comp) < 2:
                continue
            bridge_atoms: set = set()
            comp_set = set(comp)
            for (gi, gj), path in pair_paths.items():
                if gi in comp_set and gj in comp_set:
                    bridge_atoms.update(path)
            results.append({
                "metal": mi,
                "groups": comp,
                "bridge_atoms": frozenset(bridge_atoms),
                "n_groups": len(comp),
            })
    return results


def _apply_hapto_centroid_bias(
    mol,
    conf_id: int,
    hapto_groups: List[Tuple[int, List[int]]],
    min_centroid_sep: float = 2.5,
) -> bool:
    """Conservative eta-group regularization without aggressive expansion.

    Applies only rigid group translations (no internal ring distortion):
    1) enforce a minimum metal-centroid distance for each eta-group,
    2) apply a light anti-collapse repulsion only when groups are too close.
    """
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return False
    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return False

    def _centroid(indices: List[int]) -> Tuple[float, float, float]:
        pts = [conf.GetAtomPosition(i) for i in indices]
        return (
            sum(p.x for p in pts) / len(pts),
            sum(p.y for p in pts) / len(pts),
            sum(p.z for p in pts) / len(pts),
        )

    def _norm(vx: float, vy: float, vz: float) -> float:
        return math.sqrt(vx * vx + vy * vy + vz * vz)

    def _translate(indices: List[int], dx: float, dy: float, dz: float) -> None:
        for idx in indices:
            p = conf.GetAtomPosition(idx)
            conf.SetAtomPosition(idx, Point3D(p.x + dx, p.y + dy, p.z + dz))

    def _group_min_distance(a: List[int], b: List[int]) -> float:
        out = float("inf")
        for ai in a:
            pa = conf.GetAtomPosition(ai)
            for bi in b:
                pb = conf.GetAtomPosition(bi)
                d = _norm(pa.x - pb.x, pa.y - pb.y, pa.z - pb.z)
                if d < out:
                    out = d
        return out if math.isfinite(out) else 999.0

    # Normalize/validate groups once.
    group_entries: List[Dict[str, object]] = []
    for metal_idx, grp in hapto_groups:
        if metal_idx < 0 or metal_idx >= mol.GetNumAtoms():
            continue
        grp_valid = sorted(set(i for i in grp if 0 <= i < mol.GetNumAtoms()))
        if len(grp_valid) < 2:
            continue
        group_entries.append({
            "metal": metal_idx,
            "atoms": grp_valid,
        })
    if not group_entries:
        return False

    by_metal: Dict[int, List[int]] = {}
    for gi, ge in enumerate(group_entries):
        by_metal.setdefault(int(ge["metal"]), []).append(gi)

    moved = False

    # Step 1: enforce minimal metal-centroid distance by rigidly pushing groups.
    for metal_idx, gidxs in by_metal.items():
        mpos = conf.GetAtomPosition(metal_idx)
        msym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
        for gi in gidxs:
            atoms = list(group_entries[gi]["atoms"])
            cx, cy, cz = _centroid(atoms)
            vx, vy, vz = (cx - mpos.x, cy - mpos.y, cz - mpos.z)
            dist = _norm(vx, vy, vz)
            min_mc = _target_mc_dist(msym, len(atoms))

            if dist < min_mc:
                if dist < 1e-8:
                    ux, uy, uz = (1.0, 0.0, 0.0)
                else:
                    ux, uy, uz = (vx / dist, vy / dist, vz / dist)
                shift = float(min_mc - dist)
                _translate(atoms, ux * shift, uy * shift, uz * shift)
                moved = True

    # Detect ansa bridges for separation limits
    ansa_bridged_pairs: set = set()
    try:
        bridges = _find_ansa_bridges(mol, hapto_groups)
        for _bridge_atom, g1, g2 in bridges:
            ansa_bridged_pairs.add((min(g1, g2), max(g1, g2)))
    except Exception:
        pass

    # Step 2: light anti-collapse guard: only separate groups that are too close.
    safe_min_cc = max(1.85, min(2.20, float(min_centroid_sep)))
    for metal_idx, gidxs in by_metal.items():
        if len(gidxs) < 2:
            continue
        for _iter in range(4):
            changed = False
            centroids: Dict[int, Tuple[float, float, float]] = {}
            for gi in gidxs:
                atoms = list(group_entries[gi]["atoms"])
                centroids[gi] = _centroid(atoms)

            for ii in range(len(gidxs)):
                gi = gidxs[ii]
                ai = set(group_entries[gi]["atoms"])
                ci = centroids[gi]
                for jj in range(ii + 1, len(gidxs)):
                    gj = gidxs[jj]
                    aj = set(group_entries[gj]["atoms"])
                    if ai.intersection(aj):
                        continue

                    cj = centroids[gj]
                    dx, dy, dz = (cj[0] - ci[0], cj[1] - ci[1], cj[2] - ci[2])
                    dcc = _norm(dx, dy, dz)
                    eta_i = len(group_entries[gi]["atoms"])
                    eta_j = len(group_entries[gj]["atoms"])
                    pair_min_cc = safe_min_cc
                    if eta_i >= 4 and eta_j >= 4:
                        pair_min_cc = max(pair_min_cc, 3.00)
                    elif eta_i >= 4 or eta_j >= 4:
                        pair_min_cc = max(pair_min_cc, 2.60)
                    else:  # eta3/eta3
                        pair_min_cc = max(pair_min_cc, 2.30)

                    # Ansa-bridged pairs: limit max centroid separation
                    # Si-C bond ~1.88 A -> max separation ~ 2*1.88 + ring_radii
                    pair_key = (min(gi, gj), max(gi, gj))
                    max_sep = float('inf')
                    if pair_key in ansa_bridged_pairs:
                        max_sep = 5.0  # conservative limit for Si-bridged systems

                    min_cc = _group_min_distance(
                        list(group_entries[gi]["atoms"]),
                        list(group_entries[gj]["atoms"]),
                    )
                    if min_cc >= pair_min_cc:
                        continue

                    if dcc < 1e-8:
                        ux, uy, uz = (1.0, 0.0, 0.0)
                    else:
                        ux, uy, uz = (dx / dcc, dy / dcc, dz / dcc)

                    shift = min(0.45, 0.5 * float(pair_min_cc - min_cc))
                    # Limit shift for ansa-bridged pairs
                    if pair_key in ansa_bridged_pairs and dcc + 2 * shift > max_sep:
                        shift = max(0.0, 0.5 * (max_sep - dcc))
                    _translate(list(group_entries[gi]["atoms"]), -ux * shift, -uy * shift, -uz * shift)
                    _translate(list(group_entries[gj]["atoms"]), ux * shift, uy * shift, uz * shift)
                    changed = True
                    moved = True

            if not changed:
                break

    return moved


def _enforce_hapto_ring_planarity(
    mol,
    conf_id: int,
    hapto_groups: List[Tuple[int, List[int]]],
    target_cc: float = 1.40,
    n_iter: int = 5,
) -> bool:
    """Place hapto ring atoms as regular polygon with correct C-C distances.

    For cyclic groups (eta>=3, all atoms have >=2 intra-group neighbors):
      constructs a regular polygon in the plane perpendicular to the
      metal-centroid axis with C-C = target_cc.
    For chain groups: iteratively regularizes C-C distances.
    Substituent atoms (H + direct non-ring neighbors) are rigidly
    translated along with their parent ring atom.
    """
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return False
    try:
        import numpy as np
    except ImportError:
        return False
    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return False

    moved_substituents: set = set()
    changed = False
    for metal_idx, grp in hapto_groups:
        atoms = sorted(set(i for i in grp if 0 <= i < mol.GetNumAtoms()))
        eta = len(atoms)
        if eta < 3:
            continue
        if metal_idx < 0 or metal_idx >= mol.GetNumAtoms():
            continue

        # Current positions
        pts = np.array([
            [conf.GetAtomPosition(i).x,
             conf.GetAtomPosition(i).y,
             conf.GetAtomPosition(i).z]
            for i in atoms
        ], dtype=float)
        mpos = np.array([
            conf.GetAtomPosition(metal_idx).x,
            conf.GetAtomPosition(metal_idx).y,
            conf.GetAtomPosition(metal_idx).z,
        ])
        centroid = pts.mean(axis=0)

        # Normal = metal -> centroid direction
        mc_vec = centroid - mpos
        mc_dist = float(np.linalg.norm(mc_vec))
        if mc_dist < 1e-8:
            q = pts - centroid
            try:
                _u, _s, vh = np.linalg.svd(q, full_matrices=False)
            except np.linalg.LinAlgError:
                continue
            normal = vh[-1].copy()
        else:
            normal = mc_vec / mc_dist
        n_norm = float(np.linalg.norm(normal))
        if n_norm < 1e-12:
            continue
        normal = normal / n_norm

        # In-plane orthonormal basis
        if abs(normal[0]) < 0.9:
            ref = np.array([1.0, 0.0, 0.0])
        else:
            ref = np.array([0.0, 1.0, 0.0])
        u = np.cross(normal, ref)
        u /= np.linalg.norm(u)
        v = np.cross(normal, u)
        v /= np.linalg.norm(v)

        # Build intra-group adjacency
        atom_set = set(atoms)
        idx_map = {a: i for i, a in enumerate(atoms)}
        adj = [[] for _ in range(eta)]
        for i_local, atom_idx in enumerate(atoms):
            a = mol.GetAtomWithIdx(atom_idx)
            for nbr in a.GetNeighbors():
                ni = nbr.GetIdx()
                if ni in atom_set and ni != atom_idx:
                    j_local = idx_map[ni]
                    if j_local not in adj[i_local]:
                        adj[i_local].append(j_local)

        is_cyclic = all(len(adj[i]) >= 2 for i in range(eta))

        metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
        target_mc = _target_mc_dist(metal_sym, eta)

        if is_cyclic:
            # --- Regular polygon placement for cyclic rings ---
            # Find ring traversal order via graph walk
            order = []
            visited_r = [False] * eta
            cur, prev = 0, -1
            for _ in range(eta):
                order.append(cur)
                visited_r[cur] = True
                nxt = -1
                for nb in adj[cur]:
                    if nb != prev and not visited_r[nb]:
                        nxt = nb
                        break
                if nxt == -1:
                    break
                prev, cur = cur, nxt

            if len(order) != eta:
                # Fallback: angular ordering
                rel = pts - centroid
                angles = np.arctan2(rel @ v, rel @ u)
                order = list(np.argsort(angles))

            # Circumradius for regular polygon with side = target_cc
            R = target_cc / (2.0 * np.sin(np.pi / eta))

            # Starting angle: preserve first atom's angular position
            rel0 = pts[order[0]] - centroid
            start_angle = float(np.arctan2(rel0 @ v, rel0 @ u))

            new_centroid = mpos + normal * target_mc

            pts_new = np.zeros_like(pts)
            for k, idx_in_order in enumerate(order):
                angle = start_angle + 2.0 * np.pi * k / eta
                pts_new[idx_in_order] = (
                    new_centroid
                    + R * np.cos(angle) * u
                    + R * np.sin(angle) * v
                )
        else:
            # --- Chain/non-cyclic: project + iterative regularization ---
            q = pts - centroid
            proj = q - np.outer(q @ normal, normal)
            pts_new = proj + centroid

            edges = []
            for i_local in range(eta):
                for j_local in adj[i_local]:
                    if j_local > i_local:
                        edges.append((i_local, j_local))

            if edges:
                for _it in range(20):
                    for i_l, j_l in edges:
                        pi = pts_new[i_l]
                        pj = pts_new[j_l]
                        vec = pj - pi
                        dist = float(np.linalg.norm(vec))
                        if dist < 1e-8:
                            vec = np.random.default_rng(42).standard_normal(3)
                            dist = float(np.linalg.norm(vec))
                        correction = 0.5 * (dist - target_cc) / dist
                        pts_new[i_l] = pi + correction * vec
                        pts_new[j_l] = pj - correction * vec

                ctr2 = pts_new.mean(axis=0)
                q2 = pts_new - ctr2
                proj2 = q2 - np.outer(q2 @ normal, normal)
                pts_new = proj2 + ctr2

            # Adjust M-centroid distance
            centroid_new = pts_new.mean(axis=0)
            mc_v = centroid_new - mpos
            mc_d = float(np.linalg.norm(mc_v))
            if mc_d > 1e-8:
                shift = (target_mc / mc_d - 1.0) * mc_v
                pts_new = pts_new + shift

        # Write back ring atom coordinates
        new_centroid_final = pts_new.mean(axis=0)
        for i_local, atom_idx in enumerate(atoms):
            conf.SetAtomPosition(
                atom_idx,
                Point3D(float(pts_new[i_local, 0]),
                        float(pts_new[i_local, 1]),
                        float(pts_new[i_local, 2])),
            )

        # Place substituents: correct bond distance + radial outward direction
        all_hapto_atoms = set()
        for _mi, _gi in hapto_groups:
            all_hapto_atoms.update(_gi)

        # Detect bridge atoms (bonded to ring atoms in 2+ different groups)
        bridge_atoms: set = set()
        for pot in range(mol.GetNumAtoms()):
            if pot in all_hapto_atoms:
                continue
            pa = mol.GetAtomWithIdx(pot)
            grps_connected: set = set()
            for pnbr in pa.GetNeighbors():
                pni = pnbr.GetIdx()
                for g_idx, (_gm, g_atoms) in enumerate(hapto_groups):
                    if pni in set(g_atoms):
                        grps_connected.add(g_idx)
            if len(grps_connected) >= 2:
                bridge_atoms.add(pot)

        for i_local, atom_idx in enumerate(atoms):
            ring_pos = pts_new[i_local]
            outward = ring_pos - new_centroid_final
            outward_len = float(np.linalg.norm(outward))
            if outward_len > 1e-8:
                outward_unit = outward / outward_len
            else:
                outward_unit = u  # fallback

            a = mol.GetAtomWithIdx(atom_idx)
            subs = []
            for nbr in a.GetNeighbors():
                ni = nbr.GetIdx()
                if ni in atom_set or ni == metal_idx or ni in moved_substituents:
                    continue
                if ni in all_hapto_atoms or ni in bridge_atoms:
                    continue
                subs.append(nbr)

            if not subs:
                continue

            # For single substituent: place radially outward
            # For multiple: spread in a fan around the outward direction
            for s_k, nbr in enumerate(subs):
                ni = nbr.GetIdx()
                nbr_sym = nbr.GetSymbol()
                if nbr_sym == 'H':
                    bond_d = 1.08
                elif nbr_sym == 'Si':
                    bond_d = 1.87
                else:
                    bond_d = 1.50

                if len(subs) == 1:
                    # Single sub: radially outward, tilted away from metal
                    direction = outward_unit + 0.3 * normal
                    d_norm = float(np.linalg.norm(direction))
                    if d_norm > 1e-8:
                        direction = direction / d_norm
                    else:
                        direction = outward_unit
                else:
                    # Multiple subs: spread around outward direction
                    tilt_angle = (s_k - (len(subs) - 1) / 2.0) * 1.2
                    direction = (
                        outward_unit * np.cos(tilt_angle)
                        + normal * np.sin(tilt_angle)
                    )
                    d_norm = float(np.linalg.norm(direction))
                    if d_norm > 1e-8:
                        direction = direction / d_norm

                sub_pos = ring_pos + bond_d * direction
                conf.SetAtomPosition(
                    ni,
                    Point3D(float(sub_pos[0]),
                            float(sub_pos[1]),
                            float(sub_pos[2])),
                )
                moved_substituents.add(ni)
        changed = True

    # --- Post-processing: ansa bridge constraint + bridge atom placement ---
    if changed:
        try:
            import numpy as np
            conf = mol.GetConformer(conf_id)
        except Exception:
            return changed

        all_hapto_set = set()
        group_of_atom: dict = {}
        for g_idx, (_mi, _gi) in enumerate(hapto_groups):
            all_hapto_set.update(_gi)
            for ai in _gi:
                group_of_atom[ai] = g_idx

        # Find bridge atoms connecting different hapto groups
        bridges = []
        bridge_atom_set: set = set()
        for pot in range(mol.GetNumAtoms()):
            if pot in all_hapto_set:
                continue
            pa = mol.GetAtomWithIdx(pot)
            connections = []
            for pnbr in pa.GetNeighbors():
                pni = pnbr.GetIdx()
                if pni in group_of_atom:
                    connections.append((pni, group_of_atom[pni]))
            groups_seen = set(gi for _, gi in connections)
            if len(groups_seen) >= 2:
                bridges.append((pot, connections))
                bridge_atom_set.add(pot)

        # Collect movable atoms per group (ring + substituents, excluding bridges)
        group_atom_sets: dict = {}
        for g_idx, (_mi, _gi) in enumerate(hapto_groups):
            g_set = set(_gi)
            for ai in _gi:
                a = mol.GetAtomWithIdx(ai)
                for nbr in a.GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni not in all_hapto_set and ni != _mi and ni not in bridge_atom_set:
                        g_set.add(ni)
            group_atom_sets[g_idx] = g_set

        def _get_pos(idx):
            p = conf.GetAtomPosition(idx)
            return np.array([p.x, p.y, p.z])

        def _set_pos(idx, arr):
            conf.SetAtomPosition(idx, Point3D(float(arr[0]),
                                               float(arr[1]),
                                               float(arr[2])))

        def _rodrigues_rotate(points, center, axis, angle):
            k = axis / np.linalg.norm(axis)
            cos_a, sin_a = np.cos(angle), np.sin(angle)
            result = []
            for p in points:
                v = p - center
                v_rot = (v * cos_a + np.cross(k, v) * sin_a
                         + k * np.dot(k, v) * (1.0 - cos_a))
                result.append(center + v_rot)
            return result

        def _orient_bridge_atom_inward(ring_atoms, all_atoms, c_bridge,
                                        centroid_self, centroid_other, normal):
            """Rotate ring around its normal so c_bridge faces toward other ring."""
            toward = centroid_other - centroid_self
            toward_proj = toward - np.dot(toward, normal) * normal
            tp_len = float(np.linalg.norm(toward_proj))
            if tp_len < 0.1:
                return  # can't determine direction
            toward_proj /= tp_len

            v_bridge = _get_pos(c_bridge) - centroid_self
            v_proj = v_bridge - np.dot(v_bridge, normal) * normal
            vp_len = float(np.linalg.norm(v_proj))
            if vp_len < 1e-8:
                return
            v_proj /= vp_len

            cos_r = float(np.clip(np.dot(v_proj, toward_proj), -1.0, 1.0))
            sin_r = float(np.dot(np.cross(v_proj, toward_proj), normal))
            rot_angle = np.arctan2(sin_r, cos_r)
            if abs(rot_angle) < 0.01:
                return

            pts = [_get_pos(a) for a in all_atoms]
            pts_rot = _rodrigues_rotate(pts, centroid_self, normal, rot_angle)
            for ai, new_p in zip(all_atoms, pts_rot):
                _set_pos(ai, new_p)

        # Process each ansa bridge using analytical tilt angle
        for bridge_idx, connections in bridges:
            bridge_sym = mol.GetAtomWithIdx(bridge_idx).GetSymbol()
            if bridge_sym == 'Si':
                bridge_bond = 1.87
                bridge_angle_deg = 93.0
            elif bridge_sym == 'C':
                bridge_bond = 1.54
                bridge_angle_deg = 109.5
            else:
                bridge_bond = 1.80
                bridge_angle_deg = 100.0
            d_target = 2.0 * bridge_bond * np.sin(np.radians(bridge_angle_deg) / 2.0)

            by_group: dict = {}
            for ring_atom, g_idx in connections:
                by_group.setdefault(g_idx, []).append(ring_atom)
            group_indices = list(by_group.keys())
            if len(group_indices) < 2:
                continue
            gi, gj = group_indices[0], group_indices[1]
            c_a = by_group[gi][0]
            c_b = by_group[gj][0]
            metal_idx = hapto_groups[gi][0]
            metal_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
            grp_i = hapto_groups[gi][1]
            grp_j = hapto_groups[gj][1]
            atoms_i = list(group_atom_sets.get(gi, set(grp_i)))
            atoms_j = list(group_atom_sets.get(gj, set(grp_j)))

            eta_i = len(grp_i)
            eta_j = len(grp_j)
            r_i = _target_mc_dist(metal_sym, eta_i)
            r_j = _target_mc_dist(metal_sym, eta_j)
            R_i = target_cc / (2.0 * np.sin(np.pi / max(eta_i, 3)))
            R_j = target_cc / (2.0 * np.sin(np.pi / max(eta_j, 3)))
            r = (r_i + r_j) / 2.0
            R = (R_i + R_j) / 2.0

            # Analytical tilt: psi = arctan(R/r) + arcsin(d_target/(2*sqrt(r²+R²)))
            hyp = np.sqrt(r ** 2 + R ** 2)
            ratio = d_target / (2.0 * hyp)
            if abs(ratio) > 1.0:
                continue  # bridge too long for ring geometry
            psi_target = np.arctan2(R, r) + np.arcsin(ratio)
            alpha_target = 2.0 * psi_target

            pos_m = _get_pos(metal_idx)
            cent_i = np.mean([_get_pos(a) for a in grp_i], axis=0)
            cent_j = np.mean([_get_pos(a) for a in grp_j], axis=0)
            n_i = cent_i - pos_m
            n_j = cent_j - pos_m
            ni_len = float(np.linalg.norm(n_i))
            nj_len = float(np.linalg.norm(n_j))
            if ni_len < 1e-8 or nj_len < 1e-8:
                continue
            n_i /= ni_len
            n_j /= nj_len

            cos_alpha_cur = float(np.dot(n_i, n_j))
            alpha_current = np.arccos(np.clip(cos_alpha_cur, -1.0, 1.0))

            delta_tilt = (alpha_current - alpha_target) / 2.0
            if delta_tilt < 0.01:
                pass  # skip tilt, go straight to orient + place bridge
            else:
                # Tilt axis: perpendicular to plane of M, cent_i, cent_j
                rot_axis = np.cross(n_i, n_j)
                ra_len = float(np.linalg.norm(rot_axis))
                if ra_len < 1e-8:
                    if abs(n_i[0]) < 0.9:
                        rot_axis = np.cross(n_i, np.array([1.0, 0.0, 0.0]))
                    else:
                        rot_axis = np.cross(n_i, np.array([0.0, 1.0, 0.0]))
                    ra_len = float(np.linalg.norm(rot_axis))
                rot_axis /= ra_len

                # Check direction: +delta should increase n_i·n_j (reduce angle)
                test_ci = _rodrigues_rotate([cent_i], pos_m, rot_axis, 0.01)[0]
                test_ni = test_ci - pos_m
                test_ni /= np.linalg.norm(test_ni)
                if float(np.dot(test_ni, n_j)) < cos_alpha_cur:
                    rot_axis = -rot_axis

                # Apply single tilt to ring_i (+delta_tilt)
                pts_i = [_get_pos(a) for a in atoms_i]
                pts_i_rot = _rodrigues_rotate(pts_i, pos_m, rot_axis, +delta_tilt)
                for ai, new_p in zip(atoms_i, pts_i_rot):
                    _set_pos(ai, new_p)

                # Apply single tilt to ring_j (-delta_tilt)
                pts_j = [_get_pos(a) for a in atoms_j]
                pts_j_rot = _rodrigues_rotate(pts_j, pos_m, rot_axis, -delta_tilt)
                for aj, new_p in zip(atoms_j, pts_j_rot):
                    _set_pos(aj, new_p)

            # Orient bridge atoms to face each other (now that rings are tilted)
            cent_i = np.mean([_get_pos(a) for a in grp_i], axis=0)
            cent_j = np.mean([_get_pos(a) for a in grp_j], axis=0)
            n_i = cent_i - pos_m
            ni_len = float(np.linalg.norm(n_i))
            n_j = cent_j - pos_m
            nj_len = float(np.linalg.norm(n_j))
            if ni_len > 1e-8:
                n_i /= ni_len
            if nj_len > 1e-8:
                n_j /= nj_len
            _orient_bridge_atom_inward(
                grp_i, atoms_i, c_a, cent_i, cent_j, n_i)
            _orient_bridge_atom_inward(
                grp_j, atoms_j, c_b, cent_j, cent_i, n_j)

            # --- Place bridge atom with correct geometry ---
            pos_a = _get_pos(c_a)
            pos_b = _get_pos(c_b)
            mid = (pos_a + pos_b) / 2.0
            d_ab = float(np.linalg.norm(pos_a - pos_b))
            half_d = d_ab / 2.0
            h_sq = bridge_bond ** 2 - half_d ** 2
            h = float(np.sqrt(max(h_sq, 0.01)))

            away = mid - _get_pos(metal_idx)
            a_norm = float(np.linalg.norm(away))
            if a_norm > 1e-8:
                away = away / a_norm
            else:
                away = np.array([0.0, 1.0, 0.0])
            bridge_pos = mid + h * away
            _set_pos(bridge_idx, bridge_pos)

            # Place bridge substituents in local frame
            pa = mol.GetAtomWithIdx(bridge_idx)
            ca_cb = pos_b - pos_a
            x_ax = ca_cb / max(d_ab, 1e-8)
            z_ax = away
            y_ax = np.cross(z_ax, x_ax)
            y_norm = float(np.linalg.norm(y_ax))
            if y_norm > 1e-8:
                y_ax = y_ax / y_norm
            else:
                y_ax = np.array([0.0, 0.0, 1.0])

            sub_count = 0
            for pnbr in pa.GetNeighbors():
                pni = pnbr.GetIdx()
                if pni in all_hapto_set:
                    continue
                if pnbr.GetSymbol() == 'C':
                    sub_bond = 1.87
                elif pnbr.GetSymbol() == 'H':
                    sub_bond = 1.48
                else:
                    sub_bond = 1.50
                sub_angle = (sub_count - 0.5) * 2.1
                direction = (np.cos(sub_angle) * y_ax
                             + np.sin(sub_angle) * z_ax)
                d_norm = float(np.linalg.norm(direction))
                if d_norm > 1e-8:
                    direction /= d_norm
                sub_pos = bridge_pos + sub_bond * direction
                _set_pos(pni, sub_pos)
                sub_count += 1

    return changed


def _correct_hapto_geometry(
    mol,
    conf_id: int,
    hapto_groups: List[Tuple[int, List[int]]],
    target_cc: float = 1.40,
) -> bool:
    """Topology-preserving hapto geometry correction for ETKDG embeddings.

    Uses rigid-body transforms on connected fragments to fix M-centroid
    distances and ring orientations while preserving all non-metal bond
    lengths and angles from the ETKDG embedding.

    Steps:
      1. Identify fragment for each hapto group (BFS from ring atoms,
         excluding metal and other groups' ring atoms).
      2. For groups on the same metal, enforce angular separation
         (opposite sides for 2 groups, ~120° for 3, etc.).
      3. Rigid-body translate each fragment so M-centroid = target distance.
      4. Flatten ring atoms onto plane perpendicular to M-centroid axis (SVD).
      5. Gentle C-C regularization via spring iterations.
      6. For ansa-bridged systems, adjust tilt and place bridge atoms.
    """
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return False
    try:
        import numpy as np
    except ImportError:
        return False
    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return False

    # PERF (BYTE-IDENTICAL): mirror the conformer into a local float64 array and
    # route every _gp/_sp through it (eliminating the per-call RDKit boundary,
    # the dominant Python cost in the rotation / clash passes), flushing back to
    # the conformer once at the end.  _gp returned ``np.array([p.x, p.y, p.z])``
    # (float64) and _sp wrote ``Point3D(float(...))``; the conformer stores
    # doubles so the round-trip is bit-identical (verified array_equal).  The
    # only main-body returns are the early-exit guards above (before any _gp/_sp)
    # and the final ``return changed`` (flush precedes it), so no write is lost.
    _n_atoms_chg = mol.GetNumAtoms()
    coords = np.empty((_n_atoms_chg, 3), dtype=float)
    for _i in range(_n_atoms_chg):
        _p = conf.GetAtomPosition(_i)
        coords[_i, 0] = _p.x
        coords[_i, 1] = _p.y
        coords[_i, 2] = _p.z

    def _gp(i):
        return coords[i].copy()

    def _sp(i, arr):
        coords[i, 0] = float(arr[0])
        coords[i, 1] = float(arr[1])
        coords[i, 2] = float(arr[2])

    # -- Collect all hapto atoms and bridge atoms --
    all_hapto = set()
    group_of = {}
    for gi, (mi, catoms) in enumerate(hapto_groups):
        for a in catoms:
            all_hapto.add(a)
            group_of[a] = gi

    metals = set(mi for mi, _ in hapto_groups)
    # Also include ALL metal atoms (not just hapto metals) so that
    # BFS fragments don't cross into other metals' coordination spheres.
    for ai in range(mol.GetNumAtoms()):
        if mol.GetAtomWithIdx(ai).GetSymbol() in _METAL_SET:
            metals.add(ai)

    bridge_atoms: set = set()
    for ai in range(mol.GetNumAtoms()):
        if ai in all_hapto or ai in metals:
            continue
        a = mol.GetAtomWithIdx(ai)
        groups_seen = set()
        for nbr in a.GetNeighbors():
            ni = nbr.GetIdx()
            if ni in group_of:
                groups_seen.add(group_of[ni])
        if len(groups_seen) >= 2:
            bridge_atoms.add(ai)

    # -- Build exclusion set: atoms in rings containing sigma donors to
    #    other metals.  These must NOT be pulled along when rotating
    #    a hapto fragment, as they belong to another metal's coordination
    #    sphere (e.g. pyridine ring bonded to Pt in an Fe/Pt complex). --
    sigma_ring_exclude: set = set()
    ri = mol.GetRingInfo()
    try:
        all_rings = ri.AtomRings()
    except Exception:
        all_rings = []
    for m_idx in metals:
        for nbr in mol.GetAtomWithIdx(m_idx).GetNeighbors():
            ni = nbr.GetIdx()
            if ni in all_hapto or ni in metals:
                continue
            # ni is a sigma donor to metal m_idx
            # Exclude all ring atoms of rings containing ni
            for ring in all_rings:
                if ni in ring:
                    sigma_ring_exclude.update(ring)

    # Phase 6B (2026-05-12): per Wave-3A forensics, the BFS below over-extends
    # in single-metal-multi-ligand complexes.  Sigma-donors of the SAME metal
    # (and atoms one-hop from them) are pulled into the η-fragment by the
    # existing "if nni in metals and nni != mi" check (only OTHER metals
    # excluded).  Step 2 then rigid-translates the swollen fragment to the
    # η-target distance, dragging σ-coord backbone atoms (CH2/Ph) to ~2.5 Å
    # from M → spurious M-L extras and FICNAG/MIPSOW/Ir-η2-class failures.
    # Fix (env-gated DELFIN_HAPTO_FRAG_STRICT default 0): pre-compute σ-donor
    # set per metal INCLUDING same-metal σ-bonds, exclude both σ-donors and
    # their immediate neighbors from BFS fragment growth.
    same_metal_sigma_exclude: set = set()
    if hapto_groups and _class_conditional_flag(
        "DELFIN_HAPTO_FRAG_STRICT", mol, default=0,
        default_classes=("hapto", "multi_hapto"),
    ):
        try:
            for atom in mol.GetAtoms():
                if atom.GetIdx() not in metals:
                    continue
                m_idx_local = atom.GetIdx()
                for nbr in atom.GetNeighbors():
                    ni_local = nbr.GetIdx()
                    if ni_local in all_hapto:
                        continue  # this neighbor is an η-atom, not σ-donor
                    # ni_local is a σ-donor to metal m_idx_local
                    same_metal_sigma_exclude.add(ni_local)
                    # also exclude one-hop neighbors of the σ-donor
                    for nn in mol.GetAtomWithIdx(ni_local).GetNeighbors():
                        nni_local = nn.GetIdx()
                        if nni_local in metals or nni_local in all_hapto:
                            continue
                        same_metal_sigma_exclude.add(nni_local)
        except Exception:
            same_metal_sigma_exclude.clear()

    # -- Find fragment for each group (BFS, excluding metal/other groups/bridges) --
    fragments: List[set] = []
    for gi, (mi, catoms) in enumerate(hapto_groups):
        frag = set(catoms)
        queue = list(catoms)
        while queue:
            cur = queue.pop(0)
            for nbr in mol.GetAtomWithIdx(cur).GetNeighbors():
                ni = nbr.GetIdx()
                if ni in frag or ni in metals or ni in bridge_atoms:
                    continue
                if ni in all_hapto and group_of.get(ni) != gi:
                    continue
                # Exclude sigma donors to other metals and their ring atoms
                if ni in sigma_ring_exclude:
                    continue
                # Phase 6B: also exclude same-metal σ-donors + 1-hop nbrs
                if ni in same_metal_sigma_exclude:
                    continue
                bonded_to_other = False
                for nn in mol.GetAtomWithIdx(ni).GetNeighbors():
                    nni = nn.GetIdx()
                    if nni in metals and nni != mi:
                        bonded_to_other = True
                        break
                if bonded_to_other:
                    continue
                frag.add(ni)
                queue.append(ni)
        fragments.append(frag)

    # -- Group hapto entries by metal --
    by_metal: Dict[int, List[int]] = {}
    for gi, (mi, _) in enumerate(hapto_groups):
        by_metal.setdefault(mi, []).append(gi)

    changed = False

    # =====================================================================
    # STEP 1: Place hapto groups at analytically ideal positions on the
    # coordination sphere.  Uses rigid-body rotation around the metal
    # centre (preserves ALL internal distances within each fragment).
    #
    # 2 groups without ansa bridge → sandwich (~175°)
    # 2 groups with ansa bridge   → bent (~130°)
    # 3+ groups                   → evenly distributed in a plane
    # =====================================================================

    def _rodrigues(v, k, angle):
        """Rodrigues rotation of vector *v* around unit axis *k*."""
        ca, sa = np.cos(angle), np.sin(angle)
        return v * ca + np.cross(k, v) * sa + k * np.dot(k, v) * (1.0 - ca)

    def _rotate_fragment(frag_indices, centre, cur_dir, tgt_dir):
        """Rigid-body rotation of fragment so *cur_dir* maps to *tgt_dir*."""
        cd = np.asarray(cur_dir, dtype=float)
        td = np.asarray(tgt_dir, dtype=float)
        cos_a = float(np.clip(np.dot(cd, td), -1.0, 1.0))
        if cos_a > 0.9999:
            return  # already aligned
        cross = np.cross(cd, td)
        cl = float(np.linalg.norm(cross))
        if cl < 1e-8:
            # Anti-parallel – pick any perpendicular axis
            perp = np.array([1, 0, 0]) if abs(cd[0]) < 0.9 else np.array([0, 1, 0])
            cross = np.cross(cd, perp)
            cl = float(np.linalg.norm(cross))
        k = cross / cl
        rot_angle = np.arccos(cos_a)
        for ai in frag_indices:
            v = _gp(ai) - centre
            _sp(ai, centre + _rodrigues(v, k, rot_angle))

    def _groups_bridged(gi, gj):
        """True if groups *gi* and *gj* share an ansa bridge atom."""
        for bi in bridge_atoms:
            gs = set()
            for nbr in mol.GetAtomWithIdx(bi).GetNeighbors():
                ni = nbr.GetIdx()
                if ni in group_of:
                    gs.add(group_of[ni])
            if gi in gs and gj in gs:
                return True
        return False

    for mi, gidxs in by_metal.items():
        if len(gidxs) < 2:
            continue
        mpos = _gp(mi)
        n_grp = len(gidxs)

        # Current centroid unit-directions from metal
        cur_dirs: Dict[int, np.ndarray] = {}
        for gi in gidxs:
            c = np.mean([_gp(a) for a in hapto_groups[gi][1]], axis=0)
            v = c - mpos
            d = float(np.linalg.norm(v))
            cur_dirs[gi] = v / d if d > 1e-8 else np.array([1.0, 0, 0])

        # --- Compute ideal directions on sphere ---
        ideal_dirs: Dict[int, np.ndarray] = {}

        if n_grp == 2:
            gi, gj = gidxs
            bridged = _groups_bridged(gi, gj)
            target_angle = np.radians(130.0 if bridged else 175.0)

            d_i = cur_dirs[gi]
            ideal_dirs[gi] = d_i  # keep group i fixed

            # Rotation axis: perpendicular to d_i in the (d_i, d_j) plane
            cross = np.cross(d_i, cur_dirs[gj])
            cl = float(np.linalg.norm(cross))
            if cl < 1e-8:
                perp = np.array([1, 0, 0]) if abs(d_i[0]) < 0.9 else np.array([0, 1, 0])
                cross = np.cross(d_i, perp)
                cl = float(np.linalg.norm(cross))
            rot_ax = cross / cl
            ideal_dirs[gj] = _rodrigues(d_i, rot_ax, target_angle)

        elif n_grp >= 3:
            # Best-fit plane through current centroid directions
            pts = np.array([cur_dirs[gi] for gi in gidxs])
            try:
                _, _, vh = np.linalg.svd(
                    pts - pts.mean(axis=0), full_matrices=False)
                px = vh[0] / np.linalg.norm(vh[0])
                py = vh[1] / np.linalg.norm(vh[1])
            except np.linalg.LinAlgError:
                px = np.array([1, 0, 0])
                py = np.array([0, 1, 0])

            # Starting angle from first group's projection
            start = np.arctan2(
                float(np.dot(cur_dirs[gidxs[0]], py)),
                float(np.dot(cur_dirs[gidxs[0]], px)),
            )
            for k, gi in enumerate(gidxs):
                a = start + 2.0 * np.pi * k / n_grp
                d = np.cos(a) * px + np.sin(a) * py
                dl = float(np.linalg.norm(d))
                ideal_dirs[gi] = d / dl if dl > 1e-8 else px

        # --- Rotate each fragment to its ideal direction ---
        for gi in gidxs:
            if gi not in ideal_dirs:
                continue
            _rotate_fragment(fragments[gi], mpos, cur_dirs[gi], ideal_dirs[gi])
            changed = True

    # =====================================================================
    # STEP 2: Rigid-body translate fragments to target M-centroid distance
    # =====================================================================
    for gi, (mi, catoms) in enumerate(hapto_groups):
        mpos = _gp(mi)
        centroid = np.mean([_gp(a) for a in catoms], axis=0)
        mc_vec = centroid - mpos
        mc_dist = float(np.linalg.norm(mc_vec))
        if mc_dist < 1e-8:
            mc_dir = np.array([1.0, 0.0, 0.0])
            mc_dist = 1e-8
        else:
            mc_dir = mc_vec / mc_dist

        metal_sym = mol.GetAtomWithIdx(mi).GetSymbol()
        target_mc = _target_mc_dist(metal_sym, len(catoms))

        shift = mc_dir * (target_mc - mc_dist)
        if abs(target_mc - mc_dist) > 0.01:
            for ai in fragments[gi]:
                _sp(ai, _gp(ai) + shift)
            changed = True

    # =====================================================================
    # STEP 3: Rotate ring plane to be perpendicular to M-centroid axis
    # Uses rigid-body rotation (preserves ALL internal distances).
    # =====================================================================
    for gi, (mi, catoms) in enumerate(hapto_groups):
        atoms = sorted(set(catoms))
        eta = len(atoms)
        if eta < 3:
            continue
        mpos = _gp(mi)
        pts = np.array([_gp(a) for a in atoms])
        centroid = pts.mean(axis=0)
        mc_vec = centroid - mpos
        mc_dist = float(np.linalg.norm(mc_vec))
        if mc_dist < 1e-8:
            continue
        target_normal = mc_vec / mc_dist

        # Find current ring plane normal via SVD
        q = pts - centroid
        try:
            _u, _s, vh = np.linalg.svd(q, full_matrices=False)
        except np.linalg.LinAlgError:
            continue
        ring_normal = vh[-1].copy()
        # Ensure ring_normal points same way as target_normal
        if np.dot(ring_normal, target_normal) < 0:
            ring_normal = -ring_normal

        # Compute rotation from ring_normal to target_normal
        cos_a = float(np.clip(np.dot(ring_normal, target_normal), -1, 1))
        if cos_a > 0.9999:
            continue  # already aligned
        rot_axis = np.cross(ring_normal, target_normal)
        ra_len = float(np.linalg.norm(rot_axis))
        if ra_len < 1e-10:
            continue
        rot_axis /= ra_len
        angle = np.arccos(cos_a)
        cos_r = np.cos(angle)
        sin_r = np.sin(angle)
        k = rot_axis

        # Apply Rodrigues rotation to entire fragment (centered at centroid)
        for ai in fragments[gi]:
            v = _gp(ai) - centroid
            v_rot = (v * cos_r + np.cross(k, v) * sin_r
                     + k * np.dot(k, v) * (1.0 - cos_r))
            _sp(ai, centroid + v_rot)
        changed = True

    # =====================================================================
    # STEP 4: Ring shape correction.
    # For cyclic rings: reconstruct as regular polygon (optimal geometry
    # for Cp/arene hapto groups).  For non-cyclic: gentle spring correction.
    # Only ring atoms + their direct H substituents are moved.
    # Polygon starting angle optimized to minimize clashes with environment.
    # =====================================================================
    # Determine max number of groups on a single metal (for adaptive behavior)
    max_groups_per_metal = max(
        (len(gidxs) for gidxs in by_metal.values()), default=1
    )

    for gi, (mi, catoms) in enumerate(hapto_groups):
        atoms = sorted(set(catoms))
        eta = len(atoms)
        if eta < 3:
            continue
        atom_set = set(atoms)

        # Build edges (intra-ring bonds)
        edges = []
        for ai in atoms:
            for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
                ni = nbr.GetIdx()
                if ni in atom_set and ni > ai:
                    edges.append((ai, ni))

        if not edges:
            continue

        # Collect H substituents for each ring atom
        h_subs: Dict[int, List[int]] = {}
        for ai in atoms:
            hs = []
            for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
                ni = nbr.GetIdx()
                if nbr.GetSymbol() == 'H' and ni not in all_hapto:
                    hs.append(ni)
            h_subs[ai] = hs

        # Reconstruct cyclic rings as regular polygons.
        # Choose starting angle to minimize clashes with non-ring atoms.
        mpos = _gp(mi)
        ring_pts = np.array([_gp(a) for a in atoms])
        ring_centroid = ring_pts.mean(axis=0)
        mc_vec = ring_centroid - mpos
        mc_norm = float(np.linalg.norm(mc_vec))
        if mc_norm < 1e-8:
            continue
        plane_normal = mc_vec / mc_norm

        # In-plane basis
        if abs(plane_normal[0]) < 0.9:
            ref = np.array([1.0, 0.0, 0.0])
        else:
            ref = np.array([0.0, 1.0, 0.0])
        u = np.cross(plane_normal, ref)
        u /= np.linalg.norm(u)
        v = np.cross(plane_normal, u)
        v /= np.linalg.norm(v)

        # Build adjacency for ring traversal
        adj_r = [[] for _ in range(eta)]
        idx_map = {a: i for i, a in enumerate(atoms)}
        for ai in atoms:
            for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
                ni = nbr.GetIdx()
                if ni in atom_set and ni != ai:
                    jl = idx_map[ni]
                    il = idx_map[ai]
                    if jl not in adj_r[il]:
                        adj_r[il].append(jl)

        is_cyclic = all(len(adj_r[i]) >= 2 for i in range(eta))

        if not is_cyclic:
            # Non-cyclic groups (rare): gentle spring-based C-C correction
            for _it in range(30):
                max_err = 0.0
                for ai, aj in edges:
                    pi = _gp(ai)
                    pj = _gp(aj)
                    vec = pj - pi
                    d = float(np.linalg.norm(vec))
                    if d < 1e-8:
                        continue
                    err = d - target_cc
                    max_err = max(max_err, abs(err))
                    c = 0.15 * err / d * vec
                    c = c - np.dot(c, plane_normal) * plane_normal
                    _sp(ai, pi + c)
                    _sp(aj, pj - c)
                if max_err < 0.1:
                    break
            changed = True
            continue

        # -- Cyclic rings: always reconstruct as regular polygon --
        # For each ring atom, find its "branch" — all fragment atoms
        # reachable exclusively through that ring atom.  These move
        # rigidly with the ring atom (preserves substituent bond lengths).

        # Build branch for each ring atom (BFS into fragment, not
        # crossing other ring atoms or metal or bridges)
        branches: Dict[int, List[int]] = {a: [] for a in atoms}
        for ai in atoms:
            visited_br: set = set(atoms) | metals | bridge_atoms | sigma_ring_exclude
            queue_br = []
            for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
                ni = nbr.GetIdx()
                if ni not in visited_br:
                    queue_br.append(ni)
                    visited_br.add(ni)
            while queue_br:
                cur_br = queue_br.pop(0)
                branches[ai].append(cur_br)
                for nbr in mol.GetAtomWithIdx(cur_br).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni not in visited_br:
                        visited_br.add(ni)
                        queue_br.append(ni)

        # Ring traversal order
        order = []
        visited_r = [False] * eta
        cur, prev = 0, -1
        for _ in range(eta):
            order.append(cur)
            visited_r[cur] = True
            nxt = -1
            for nb in adj_r[cur]:
                if nb != prev and not visited_r[nb]:
                    nxt = nb
                    break
            if nxt == -1:
                break
            prev, cur = cur, nxt

        if len(order) != eta:
            rel = ring_pts - ring_centroid
            angles_f = np.arctan2(rel @ v, rel @ u)
            order = list(np.argsort(angles_f))

        # Circumradius for regular polygon
        R = target_cc / (2.0 * np.sin(np.pi / eta))

        # Collect atom indices that belong to THIS group's branches
        # (these move with the ring, so exclude from clash scoring)
        own_branch_atoms = set(atoms)
        for ai in atoms:
            own_branch_atoms.update(branches[ai])

        # Collect non-ring, non-metal, non-own-branch heavy atom positions
        # for clash scoring
        env_atoms = []
        for eidx in range(mol.GetNumAtoms()):
            if eidx in own_branch_atoms or eidx in metals or eidx in bridge_atoms:
                continue
            if mol.GetAtomWithIdx(eidx).GetSymbol() == 'H':
                continue
            env_atoms.append(_gp(eidx))

        # PERF (BYTE-IDENTICAL): stack the (invariant) environment positions so
        # the per-(env × ring-point) distance — the profiled norm-storm here —
        # can be evaluated with one matmul-ddot row-norm per env atom instead of
        # eta scalar np.linalg.norm calls.  np.matmul(d[:,None,:], d[:,:,None])
        # dispatches the SAME per-row BLAS ddot kernel as scalar
        # np.linalg.norm(v)=sqrt(v.dot(v)) → each distance is bit-identical
        # (FMA/summation order preserved; axis-norm / (x**2).sum drift ~1 ULP
        # and are rejected).  The score ``+=`` is kept SCALAR in the exact
        # original ``for ep: for pi:`` order so the float summation — and the
        # selected best_pts — is unchanged.
        env_arr = (np.asarray(env_atoms, dtype=float)
                   if env_atoms else np.zeros((0, 3), dtype=float))
        n_env = env_arr.shape[0]
        # Try multiple starting angles, pick the one with fewest clashes
        best_pts = None
        best_clash_score = float('inf')
        n_angles = 12
        for a_try in range(n_angles):
            start_angle = 2.0 * np.pi * a_try / n_angles
            pts_try = np.zeros((eta, 3))
            for k, idx_in_order in enumerate(order):
                angle = start_angle + 2.0 * np.pi * k / eta
                pts_try[idx_in_order] = (
                    ring_centroid
                    + R * np.cos(angle) * u
                    + R * np.sin(angle) * v
                )
            # Score: sum of 1/d for close contacts with environment
            clash_score = 0.0
            if n_env:
                for e_i in range(n_env):
                    deltas = pts_try - env_arr[e_i]
                    dists = np.sqrt(
                        np.matmul(deltas[:, None, :],
                                  deltas[:, :, None]).reshape(-1)
                    )
                    for p_i in range(eta):
                        d = float(dists[p_i])
                        if d < 1.5:
                            clash_score += (1.5 - d) ** 2
            if clash_score < best_clash_score:
                best_clash_score = clash_score
                best_pts = pts_try

        if best_pts is None:
            continue

        # Apply polygon positions, moving entire branches with each ring atom
        for i_loc, ai in enumerate(atoms):
            delta = best_pts[i_loc] - _gp(ai)
            if float(np.linalg.norm(delta)) < 0.001:
                continue
            _sp(ai, best_pts[i_loc])
            for bi in branches[ai]:
                _sp(bi, _gp(bi) + delta)

        # Orient branches outward: if a branch root atom is closer to the
        # ring centroid than its ring parent, reflect it across the parent
        # in the outward direction (away from centroid).
        for ai in atoms:
            if not branches[ai]:
                continue
            ring_pos = _gp(ai)
            outward = ring_pos - ring_centroid
            outward_len = float(np.linalg.norm(outward))
            if outward_len < 1e-8:
                continue
            outward_unit = outward / outward_len
            # Check first branch atom (root of substituent)
            root = branches[ai][0]
            root_pos = _gp(root)
            branch_dir = root_pos - ring_pos
            # If branch points inward (toward centroid), reflect it
            if float(np.dot(branch_dir, outward_unit)) < 0:
                # Reflect all branch atoms across the plane perpendicular
                # to outward through the ring atom
                for bi in branches[ai]:
                    bp = _gp(bi)
                    v = bp - ring_pos
                    proj = float(np.dot(v, outward_unit))
                    # Reflect across plane perpendicular to outward
                    _sp(bi, bp - 2.0 * proj * outward_unit)

        changed = True

    # =====================================================================
    # STEP 5: Ansa bridge handling
    # =====================================================================
    if bridge_atoms and changed:
        for bridge_idx in bridge_atoms:
            ba = mol.GetAtomWithIdx(bridge_idx)
            bridge_sym = ba.GetSymbol()
            if bridge_sym == 'Si':
                bridge_bond = 1.87
            elif bridge_sym == 'C':
                bridge_bond = 1.54
            else:
                bridge_bond = 1.80

            # Find ring atoms this bridge connects
            ring_nbrs = []
            for nbr in ba.GetNeighbors():
                ni = nbr.GetIdx()
                if ni in all_hapto:
                    ring_nbrs.append((ni, group_of[ni]))
            if len(ring_nbrs) < 2:
                continue

            # Place bridge at midpoint of connected ring atoms, pushed outward
            ring_pos = [_gp(rn) for rn, _ in ring_nbrs]
            mid = np.mean(ring_pos, axis=0)

            # Find the metal for outward direction
            mi_bridge = hapto_groups[ring_nbrs[0][1]][0]
            mpos = _gp(mi_bridge)
            away = mid - mpos
            aw_len = float(np.linalg.norm(away))
            if aw_len > 1e-8:
                away /= aw_len
            else:
                away = np.array([0.0, 1.0, 0.0])

            # Distance between the two ring attachment points
            if len(ring_pos) >= 2:
                d_ab = float(np.linalg.norm(ring_pos[0] - ring_pos[1]))
                half_d = d_ab / 2.0
                h_sq = bridge_bond ** 2 - half_d ** 2
                h = float(np.sqrt(max(h_sq, 0.01)))
            else:
                h = bridge_bond * 0.7

            bridge_pos = mid + h * away

            # Ensure bridge atom is at least 2.5A from metal
            bm_vec = bridge_pos - mpos
            bm_dist = float(np.linalg.norm(bm_vec))
            min_bridge_metal = 2.5
            if bm_dist < min_bridge_metal and bm_dist > 1e-8:
                bridge_pos = mpos + bm_vec / bm_dist * min_bridge_metal
            elif bm_dist < 1e-8:
                bridge_pos = mpos + away * min_bridge_metal

            _sp(bridge_idx, bridge_pos)
            changed = True

            # Place non-hapto substituents on the bridge atom
            sub_count = 0
            ca_cb = ring_pos[1] - ring_pos[0] if len(ring_pos) >= 2 else np.array([1, 0, 0])
            ca_cb_len = float(np.linalg.norm(ca_cb))
            if ca_cb_len > 1e-8:
                x_ax = ca_cb / ca_cb_len
            else:
                x_ax = np.array([1.0, 0.0, 0.0])
            z_ax = away
            y_ax = np.cross(z_ax, x_ax)
            y_norm = float(np.linalg.norm(y_ax))
            if y_norm > 1e-8:
                y_ax /= y_norm
            else:
                y_ax = np.array([0.0, 0.0, 1.0])

            for nbr in ba.GetNeighbors():
                ni = nbr.GetIdx()
                if ni in all_hapto:
                    continue
                nsym = nbr.GetSymbol()
                if nsym == 'C':
                    sub_bond = 1.87
                elif nsym == 'H':
                    sub_bond = 1.48
                else:
                    sub_bond = 1.50
                sub_angle = (sub_count - 0.5) * 2.1
                direction = (np.cos(sub_angle) * y_ax
                             + np.sin(sub_angle) * z_ax)
                d_n = float(np.linalg.norm(direction))
                if d_n > 1e-8:
                    direction /= d_n
                _sp(ni, bridge_pos + sub_bond * direction)
                sub_count += 1

    # =====================================================================
    # STEP 6: Fix H orientations on ring carbons (radially outward + tilt)
    # =====================================================================
    for gi, (mi, catoms) in enumerate(hapto_groups):
        atoms = sorted(set(catoms))
        mpos = _gp(mi)
        centroid = np.mean([_gp(a) for a in atoms], axis=0)
        mc_vec = centroid - mpos
        mc_dist = float(np.linalg.norm(mc_vec))
        if mc_dist < 1e-8:
            continue
        normal = mc_vec / mc_dist

        for ai in atoms:
            ring_pos = _gp(ai)
            outward = ring_pos - centroid
            out_len = float(np.linalg.norm(outward))
            if out_len > 1e-8:
                outward_unit = outward / out_len
            else:
                continue

            a = mol.GetAtomWithIdx(ai)
            h_nbrs = []
            for nbr in a.GetNeighbors():
                ni = nbr.GetIdx()
                if nbr.GetSymbol() == 'H' and ni not in all_hapto:
                    h_nbrs.append(ni)

            if not h_nbrs:
                continue

            # Place H atoms radially outward, tilted away from metal
            direction = outward_unit + 0.3 * normal
            d_n = float(np.linalg.norm(direction))
            if d_n > 1e-8:
                direction /= d_n
            for hi in h_nbrs:
                _sp(hi, ring_pos + 1.08 * direction)
                changed = True

    # =====================================================================
    # STEP 7: Clash resolution — push apart any non-bonded atoms that
    # are too close.  Hapto ring atoms are frozen; only non-ring atoms
    # are pushed outward from the metal or from each other.
    # Also enforce minimum distance from metals for non-bonded atoms.
    # =====================================================================
    if changed:
        # Determine minimum allowed distances
        frozen = set()
        for _mi, catoms in hapto_groups:
            frozen.update(catoms)
            frozen.add(_mi)

        # Build set of atom pairs that are directly bonded (should not be
        # pushed apart — bond length is the embedding's responsibility)
        bonded_pairs: set = set()
        for bond in mol.GetBonds():
            ai, aj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bonded_pairs.add((min(ai, aj), max(ai, aj)))

        # Build "rigid group" for each non-frozen heavy atom: the atom
        # itself plus all its directly bonded H atoms that are also not
        # frozen.  When a heavy atom is pushed, its H group moves with it.
        rigid_group: Dict[int, List[int]] = {}
        h_parent: Dict[int, int] = {}  # H atom → parent heavy atom
        for ai in range(mol.GetNumAtoms()):
            if ai in frozen:
                continue
            a = mol.GetAtomWithIdx(ai)
            if a.GetSymbol() == 'H':
                continue
            hs = []
            for nbr in a.GetNeighbors():
                ni = nbr.GetIdx()
                if nbr.GetSymbol() == 'H' and ni not in frozen:
                    hs.append(ni)
                    h_parent[ni] = ai
            rigid_group[ai] = hs

        def _push_atom(idx, displacement):
            """Push atom + its bonded H group by displacement vector.
            If idx is an H atom, redirect to its parent heavy atom
            (unless the parent is frozen)."""
            if idx in h_parent:
                parent = h_parent[idx]
                if parent not in frozen:
                    idx = parent
                # else: push H directly (parent frozen, can't move it)
            if idx in frozen:
                return  # safety: never move frozen atoms
            _sp(idx, _gp(idx) + displacement)
            if idx in rigid_group:
                for hi in rigid_group[idx]:
                    _sp(hi, _gp(hi) + displacement)

        n_atoms = mol.GetNumAtoms()
        # PERF (BYTE-IDENTICAL): the candidate (i, j) pairs, their frozen flags
        # and the static min_d threshold depend only on frozen / bonded_pairs /
        # symbols — none change across the 40 passes — so hoist them once in the
        # exact original (outer i, inner j) order; the per-pass mutation
        # sequence is unchanged.  The distance stays scalar
        # ``sqrt(diff.dot(diff))`` (= np.linalg.norm on a 1-D vector, same BLAS
        # ddot → bit-identical) because the loop pushes atoms mid-pass and is
        # therefore inherently sequential.
        _sym_chg = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(n_atoms)]
        step6_pairs = []
        for i in range(n_atoms):
            i_frozen = i in frozen
            i_is_metal = i in metals
            si = _sym_chg[i]
            for j in range(i + 1, n_atoms):
                j_frozen = j in frozen
                if i_frozen and j_frozen:
                    continue
                if (i, j) in bonded_pairs:
                    continue
                sj = _sym_chg[j]
                j_is_metal = j in metals
                if i_is_metal or j_is_metal:
                    other_sym = sj if i_is_metal else si
                    min_d = 2.5 if other_sym == 'H' else 2.2
                elif si == 'H' and sj == 'H':
                    min_d = 1.5
                elif si == 'H' or sj == 'H':
                    min_d = 1.0
                else:
                    min_d = 1.2
                step6_pairs.append((i, j, i_frozen, j_frozen, min_d))

        for _pass in range(40):
            any_push = False
            for i, j, i_frozen, j_frozen, min_d in step6_pairs:
                pi = _gp(i)
                pj = _gp(j)
                diff = pi - pj
                d = float(np.sqrt(diff.dot(diff)))

                if d >= min_d:
                    continue

                if d < 1e-8:
                    direction = np.random.default_rng(
                        42 + _pass).standard_normal(3)
                    direction /= np.linalg.norm(direction)
                else:
                    direction = (pj - pi) / d

                gap = min_d - d
                if i_frozen:
                    _push_atom(j, gap * direction)
                elif j_frozen:
                    _push_atom(i, -gap * direction)
                else:
                    _push_atom(i, -0.5 * gap * direction)
                    _push_atom(j, 0.5 * gap * direction)
                any_push = True

            if not any_push:
                break

        # Final bond repair: fix any bonds distorted by the clash pushes.
        # For each bond, if the distance deviates >20% from the expected
        # bond length, move the lighter/non-frozen atom to correct it.
        for bond in mol.GetBonds():
            ai = bond.GetBeginAtomIdx()
            aj = bond.GetEndAtomIdx()
            si = mol.GetAtomWithIdx(ai).GetSymbol()
            sj = mol.GetAtomWithIdx(aj).GetSymbol()
            if si in _METAL_SET or sj in _METAL_SET:
                continue  # skip metal bonds
            pi = _gp(ai)
            pj = _gp(aj)
            d = float(np.linalg.norm(pi - pj))
            # Expected bond lengths
            if si == 'H' or sj == 'H':
                expected = 1.08
            elif si == 'Si' or sj == 'Si':
                expected = 1.87
            else:
                expected = 1.50
            if d < 1e-8 or abs(d - expected) / expected < 0.20:
                continue
            # Move the non-frozen (or lighter) atom
            direction = (pj - pi) / d
            correction = expected - d
            ai_frozen = ai in frozen
            aj_frozen = aj in frozen
            if ai_frozen and not aj_frozen:
                _sp(aj, pj + correction * direction)
            elif aj_frozen and not ai_frozen:
                _sp(ai, pi - correction * direction)
            elif not ai_frozen and not aj_frozen:
                # Move H toward its parent, or split equally
                if si == 'H':
                    _sp(ai, pi - correction * direction)
                elif sj == 'H':
                    _sp(aj, pj + correction * direction)
                else:
                    _sp(ai, pi - 0.5 * correction * direction)
                    _sp(aj, pj + 0.5 * correction * direction)

    # Flush the local array back to the conformer (bit-identical write-back).
    for _i in range(_n_atoms_chg):
        conf.SetAtomPosition(_i, Point3D(
            float(coords[_i, 0]), float(coords[_i, 1]), float(coords[_i, 2])))

    return changed


# ---------------------------------------------------------------------------
# BFS propagation for non-hapto atoms after Phase 1 correction
# ---------------------------------------------------------------------------

def _propagate_non_hapto_atoms(
    mol,
    conf_id,
    hapto_groups,
    extra_fixed_indices: Optional[set] = None,
):
    """Re-place non-hapto atoms using BFS from the fixed hapto scaffold.

    After Phase 1 corrects hapto ring geometry, atoms NOT reachable from
    hapto rings (without crossing metals) may still be in bad ETKDG positions.
    These are typically sigma ligands and their substituents.  This function
    identifies them and re-places them using VSEPR-based local geometry,
    followed by bond-length relaxation and clash resolution.
    """
    from collections import deque
    import numpy as np

    conf = mol.GetConformer(conf_id)
    n_atoms = mol.GetNumAtoms()

    # PERF (BYTE-IDENTICAL): mirror the conformer into a local float64 array.
    # _gp returned ``np.array([p.x, p.y, p.z])`` (float64) and _sp wrote it
    # back via ``Point3D(float(...))``; the conformer stores doubles, so the
    # read→array→write round-trip is bit-identical (verified array_equal).  The
    # placement BFS and the two relaxation/clash passes mutate positions in
    # place and read them back sequentially, so every _gp/_sp is routed through
    # the local array (eliminating the per-call RDKit boundary, the dominant
    # Python cost) and the conformer is flushed once at the end.
    coords = np.empty((n_atoms, 3), dtype=float)
    for _i in range(n_atoms):
        _p = conf.GetAtomPosition(_i)
        coords[_i, 0] = _p.x
        coords[_i, 1] = _p.y
        coords[_i, 2] = _p.z

    def _gp(i):
        return coords[i].copy()

    def _sp(i, pos):
        coords[i, 0] = float(pos[0])
        coords[i, 1] = float(pos[1])
        coords[i, 2] = float(pos[2])

    # Identify metals and hapto ring atoms
    metals = set()
    hapto_ring_atoms = set()
    for mi, catoms in hapto_groups:
        metals.add(mi)
        for a in catoms:
            hapto_ring_atoms.add(a)

    # Compute hapto-reachable set: BFS from ring atoms, not crossing metals.
    # These atoms were already corrected by branch-aware rotation in Phase 1.
    hapto_reachable = set(hapto_ring_atoms)
    bfs_q = deque(list(hapto_ring_atoms))
    while bfs_q:
        cur = bfs_q.popleft()
        for nbr in mol.GetAtomWithIdx(cur).GetNeighbors():
            ni = nbr.GetIdx()
            if ni not in hapto_reachable and ni not in metals:
                hapto_reachable.add(ni)
                bfs_q.append(ni)

    # Fixed = metals + hapto_reachable (already in correct positions)
    fixed = metals | hapto_reachable
    if extra_fixed_indices:
        fixed = fixed | set(extra_fixed_indices)

    # Find atoms that need re-placement
    needs_placement = set()
    for ai in range(n_atoms):
        if ai not in fixed:
            needs_placement.add(ai)

    if not needs_placement:
        return  # all atoms already in good positions

    # Bond length helper
    def _bl(sym_a, sym_b):
        if sym_a == 'H' or sym_b == 'H':
            other = sym_b if sym_a == 'H' else sym_a
            return {'C': 1.08, 'N': 1.01, 'O': 0.96, 'Si': 1.48,
                    'S': 1.34, 'P': 1.42}.get(other, 1.08)
        pair = frozenset([sym_a, sym_b])
        _bl_map = {
            frozenset(['C', 'C']): 1.50, frozenset(['C', 'N']): 1.47,
            frozenset(['C', 'O']): 1.43, frozenset(['C', 'Si']): 1.87,
            frozenset(['C', 'S']): 1.82, frozenset(['C', 'P']): 1.84,
            frozenset(['N', 'N']): 1.45, frozenset(['N', 'O']): 1.40,
            frozenset(['Si', 'Si']): 2.34,
        }
        if pair in _bl_map:
            return _bl_map[pair]
        r1 = _COVALENT_RADII.get(sym_a, 0.76)
        r2 = _COVALENT_RADII.get(sym_b, 0.76)
        return r1 + r2

    # BFS from fixed atoms to place non-hapto atoms
    placed = set(fixed)
    queue = deque()
    in_queue = set()

    for ai in fixed:
        for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
            ni = nbr.GetIdx()
            if ni not in placed and ni not in in_queue:
                queue.append((ni, ai))
                in_queue.add(ni)

    while queue:
        ai, parent_idx = queue.popleft()
        if ai in placed:
            continue

        parent_pos = _gp(parent_idx)
        child_sym = mol.GetAtomWithIdx(ai).GetSymbol()
        parent_sym = mol.GetAtomWithIdx(parent_idx).GetSymbol()

        if parent_sym in _METAL_SET:
            bond_len = float(_get_ml_bond_length(parent_sym, child_sym))
        elif child_sym in _METAL_SET:
            bond_len = float(_get_ml_bond_length(child_sym, parent_sym))
        else:
            bond_len = _bl(parent_sym, child_sym)

        # VSEPR: direction away from already-placed neighbors of parent
        parent_atom = mol.GetAtomWithIdx(parent_idx)
        used_dirs = []
        for pnbr in parent_atom.GetNeighbors():
            pni = pnbr.GetIdx()
            if pni in placed and pni != ai:
                d = _gp(pni) - parent_pos
                d_len = float(np.linalg.norm(d))
                if d_len > 1e-8:
                    used_dirs.append(d / d_len)

        if len(used_dirs) == 0:
            direction = np.array([0.0, 0.0, 1.0])
        elif len(used_dirs) == 1:
            d0 = used_dirs[0]
            ref = (np.array([1.0, 0.0, 0.0]) if abs(d0[0]) < 0.9
                   else np.array([0.0, 1.0, 0.0]))
            perp = np.cross(d0, ref)
            perp = perp / max(float(np.linalg.norm(perp)), 1e-12)
            # ~109.5° from existing bond (tetrahedral)
            direction = (-d0 * np.cos(np.radians(70.5))
                         + perp * np.sin(np.radians(70.5)))
        elif len(used_dirs) == 2:
            avg = (used_dirs[0] + used_dirs[1]) / 2.0
            avg_len = float(np.linalg.norm(avg))
            if avg_len > 1e-8:
                direction = -avg / avg_len
            else:
                perp = np.cross(used_dirs[0], used_dirs[1])
                p_len = float(np.linalg.norm(perp))
                direction = (perp / p_len if p_len > 1e-8
                             else np.array([0.0, 0.0, 1.0]))
        else:
            avg = sum(used_dirs) / len(used_dirs)
            avg_len = float(np.linalg.norm(avg))
            direction = (-avg / avg_len if avg_len > 1e-8
                         else np.array([0.0, 0.0, 1.0]))

        d_norm = float(np.linalg.norm(direction))
        if d_norm > 1e-8:
            direction = direction / d_norm

        _sp(ai, parent_pos + bond_len * direction)
        placed.add(ai)

        for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
            ni = nbr.GetIdx()
            if ni not in placed and ni not in in_queue:
                queue.append((ni, ai))
                in_queue.add(ni)

    # Bond-length relaxation for non-metal bonds
    bond_targets = []
    for bond in mol.GetBonds():
        bi = bond.GetBeginAtomIdx()
        bj = bond.GetEndAtomIdx()
        si = mol.GetAtomWithIdx(bi).GetSymbol()
        sj = mol.GetAtomWithIdx(bj).GetSymbol()
        if si in _METAL_SET or sj in _METAL_SET:
            continue
        # Only relax bonds involving at least one re-placed atom
        if bi not in needs_placement and bj not in needs_placement:
            continue
        bond_targets.append((bi, bj, _bl(si, sj)))

    rng = np.random.default_rng(42)
    for _relax_pass in range(60):
        max_err = 0.0
        forces = {}
        for bi, bj, target in bond_targets:
            diff = coords[bj] - coords[bi]
            # bit-identical to float(np.linalg.norm(diff)) (1-D norm = sqrt of
            # the BLAS ddot); .dot is the fastest scalar form (verified).
            d = float(np.sqrt(diff.dot(diff)))
            err = abs(d - target)
            if err < 0.01:
                continue
            if d < 1e-8:
                diff = rng.standard_normal(3)
                d = float(np.linalg.norm(diff))
            unit = diff / d
            force = 0.3 * (d - target) * unit
            forces.setdefault(bi, np.zeros(3))
            forces.setdefault(bj, np.zeros(3))
            forces[bi] = forces[bi] + force
            forces[bj] = forces[bj] - force
            if err > max_err:
                max_err = err
        if max_err < 0.05:
            break
        for ai, f in forces.items():
            if ai in metals:
                continue
            w = 0.1 if ai in hapto_ring_atoms else 1.0
            _sp(ai, _gp(ai) + f * w)

    # Clash resolution for re-placed atoms
    bonded_pairs = set()
    for bond in mol.GetBonds():
        bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonded_pairs.add((bi, bj))
        bonded_pairs.add((bj, bi))

    # PERF (BYTE-IDENTICAL): the active (i, j) pair set and each pair's static
    # min_d threshold depend only on symbols / bonded / needs_placement, none
    # of which change across the 20 passes — hoist them once in the exact
    # original (outer i, inner j) iteration order so the per-pass mutation
    # sequence is unchanged.  The distance stays scalar (the loop mutates
    # positions mid-pass, so it is inherently sequential).
    _sym_p = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(n_atoms)]
    clash_pairs = []
    for i in range(n_atoms):
        si = _sym_p[i]
        si_metal = si in _METAL_SET
        si_h = si == 'H'
        i_in_np = i in needs_placement
        for j in range(i + 1, n_atoms):
            if (i, j) in bonded_pairs:
                continue
            if not i_in_np and j not in needs_placement:
                continue
            sj = _sym_p[j]
            if si_metal or sj in _METAL_SET:
                min_d = 2.0
            elif si_h and sj == 'H':
                min_d = 1.5
            elif si_h or sj == 'H':
                min_d = 1.0
            else:
                min_d = 1.2
            clash_pairs.append((i, j, min_d))

    for _pass in range(20):
        any_push = False
        for i, j, min_d in clash_pairs:
            pi = _gp(i)
            pj = _gp(j)
            diff = pi - pj
            d = float(np.sqrt(diff.dot(diff)))

            if d >= min_d:
                continue
            if d < 1e-8:
                direction = rng.standard_normal(3)
                direction /= max(float(np.linalg.norm(direction)), 1e-12)
            else:
                direction = (pj - pi) / d

            gap = min_d - d
            i_fixed = i in fixed
            j_fixed = j in fixed
            if i_fixed:
                _sp(j, pj + gap * direction)
            elif j_fixed:
                _sp(i, pi - gap * direction)
            else:
                _sp(i, pi - 0.5 * gap * direction)
                _sp(j, pj + 0.5 * gap * direction)
            any_push = True
        if not any_push:
            break

    # Flush the local array back to the conformer (bit-identical write-back).
    for _i in range(n_atoms):
        conf.SetAtomPosition(_i, Point3D(
            float(coords[_i, 0]), float(coords[_i, 1]), float(coords[_i, 2])))


# ---------------------------------------------------------------------------
# Sphere-based constructive hapto scaffold builder
# ---------------------------------------------------------------------------

def _build_hapto_scaffold(
    mol,
    hapto_groups: List[Tuple[int, List[int]]],
) -> bool:
    """Construct hapto complex geometry using sphere-based ligand placement.

    Instead of relying on RDKit's distance geometry (ETKDG), which has no
    knowledge of eta-coordination, this builds the geometry analytically:

    1. Place each metal at a fixed position
    2. Distribute hapto-group centroids on a sphere (r = M-centroid distance)
    3. Construct each hapto ring as a regular polygon perpendicular to M->centroid
    4. Place sigma-bonded ligands in remaining coordination directions
    5. Place substituents with correct bond lengths and angles
    6. BFS-propagate remaining atoms with local VSEPR geometry rules

    The mol gets a new conformer with the constructed coordinates.
    Returns True on success.
    """
    if not RDKIT_AVAILABLE or mol is None or not hapto_groups:
        return False
    try:
        import numpy as np
    except ImportError:
        return False

    n_atoms = mol.GetNumAtoms()
    coords = np.full((n_atoms, 3), np.nan)
    placed: set = set()

    # ---- Classify atoms ----
    by_metal: Dict[int, List[List[int]]] = {}
    all_hapto_atoms: set = set()
    atom_to_hapto_group: Dict[int, int] = {}  # atom_idx -> group index
    for g_idx, (metal_idx, grp) in enumerate(hapto_groups):
        by_metal.setdefault(metal_idx, []).append(grp)
        all_hapto_atoms.update(grp)
        for a in grp:
            atom_to_hapto_group[a] = g_idx
    metal_indices = sorted(by_metal.keys())

    sigma_donors: Dict[int, List[int]] = {}
    for mi in metal_indices:
        donors = []
        for nbr in mol.GetAtomWithIdx(mi).GetNeighbors():
            ni = nbr.GetIdx()
            if ni not in all_hapto_atoms and nbr.GetSymbol() not in _METAL_SET:
                donors.append(ni)
        sigma_donors[mi] = donors

    # Detect ansa bridges (multiple bridges per group pair supported)
    ansa_bridge_map: Dict[Tuple[int, int], List[int]] = {}
    all_bridge_atoms: set = set()
    try:
        bridges = _find_ansa_bridges(mol, hapto_groups)
        for bridge_atom, gi, gj in bridges:
            key = (min(gi, gj), max(gi, gj))
            ansa_bridge_map.setdefault(key, []).append(bridge_atom)
            all_bridge_atoms.add(bridge_atom)
    except Exception:
        pass

    # ---- Helper functions ----
    def _ortho_basis(normal):
        n = normal / max(float(np.linalg.norm(normal)), 1e-12)
        ref = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        u = np.cross(n, ref)
        u = u / max(float(np.linalg.norm(u)), 1e-12)
        v = np.cross(n, u)
        v = v / max(float(np.linalg.norm(v)), 1e-12)
        return u, v

    def _sphere_dirs(n):
        if n <= 0:
            return []
        if n == 1:
            return [np.array([0.0, 0.0, 1.0])]
        if n == 2:
            return [np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, -1.0])]
        dirs = []
        golden = (1.0 + np.sqrt(5.0)) / 2.0
        for k in range(n):
            theta = np.arccos(1.0 - 2.0 * (k + 0.5) / n)
            phi = 2.0 * np.pi * k / golden
            dirs.append(np.array([
                np.sin(theta) * np.cos(phi),
                np.sin(theta) * np.sin(phi),
                np.cos(theta),
            ]))
        return dirs

    def _ring_traversal(grp):
        atom_set = set(grp)
        idx_map = {a: i for i, a in enumerate(grp)}
        eta = len(grp)
        adj = [[] for _ in range(eta)]
        for i_loc, atom_idx in enumerate(grp):
            for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                ni = nbr.GetIdx()
                if ni in atom_set and ni != atom_idx:
                    j_loc = idx_map[ni]
                    if j_loc not in adj[i_loc]:
                        adj[i_loc].append(j_loc)
        is_cyclic = all(len(adj[i]) >= 2 for i in range(eta))
        if is_cyclic and eta >= 3:
            order = []
            visited = [False] * eta
            cur, prev = 0, -1
            for _ in range(eta):
                order.append(cur)
                visited[cur] = True
                nxt = -1
                for nb in adj[cur]:
                    if nb != prev and not visited[nb]:
                        nxt = nb
                        break
                if nxt == -1:
                    break
                prev, cur = cur, nxt
            if len(order) == eta:
                return [grp[i] for i in order]
        return grp

    def _build_chain_order(grp, mol_ref):
        """Order a non-cyclic hapto group as a chain following molecular topology.

        For disconnected hapto groups (e.g., two eta2 pairs bridged by
        non-hapto atoms), determines ordering by shortest path through
        the full molecular graph.
        """
        grp_set = set(grp)
        adj_intra: Dict[int, List[int]] = {a: [] for a in grp}
        for a in grp:
            for nb in mol_ref.GetAtomWithIdx(a).GetNeighbors():
                if nb.GetIdx() in grp_set:
                    adj_intra[a].append(nb.GetIdx())

        # Find connected components within the group
        comp_id: Dict[int, int] = {}
        components: List[List[int]] = []
        for a in grp:
            if a in comp_id:
                continue
            comp: List[int] = []
            stack = [a]
            while stack:
                cur = stack.pop()
                if cur in comp_id:
                    continue
                comp_id[cur] = len(components)
                comp.append(cur)
                for nb in adj_intra[cur]:
                    if nb not in comp_id:
                        stack.append(nb)
            components.append(comp)

        if len(components) == 1:
            # Single connected component: simple chain walk
            start = grp[0]
            for a in grp:
                if len(adj_intra[a]) <= 1:
                    start = a
                    break
            chain = [start]
            visited = {start}
            while len(chain) < len(grp):
                cur = chain[-1]
                nxt = None
                for nb in adj_intra[cur]:
                    if nb not in visited:
                        nxt = nb
                        break
                if nxt is None:
                    break
                chain.append(nxt)
                visited.add(nxt)
            return chain

        # Multiple disconnected components: order them by BFS shortest
        # path through non-hapto atoms in the molecular graph
        def _bfs_dist(src, tgt_set):
            """BFS from src to any atom in tgt_set, through non-metal atoms."""
            from collections import deque
            q = deque([(src, 0)])
            seen = {src}
            while q:
                cur, d = q.popleft()
                if cur in tgt_set:
                    return d
                for nb in mol_ref.GetAtomWithIdx(cur).GetNeighbors():
                    ni = nb.GetIdx()
                    if ni not in seen and nb.GetSymbol() not in _METAL_SET:
                        seen.add(ni)
                        q.append((ni, d + 1))
            return 999

        # Order components by shortest path between consecutive pairs
        ordered_comps = [components[0]]
        remaining = list(range(1, len(components)))
        while remaining:
            last_comp = ordered_comps[-1]
            best_ci = remaining[0]
            best_d = 999
            for ci in remaining:
                # Find min BFS distance from any atom in last_comp to comp[ci]
                comp_set = set(components[ci])
                for a in last_comp:
                    d = _bfs_dist(a, comp_set)
                    if d < best_d:
                        best_d = d
                        best_ci = ci
            ordered_comps.append(components[best_ci])
            remaining.remove(best_ci)

        # Within each component, order as chain from endpoint
        chain = []
        for comp in ordered_comps:
            if len(comp) == 1:
                chain.append(comp[0])
                continue
            start = comp[0]
            for a in comp:
                if len(adj_intra[a]) <= 1:
                    start = a
                    break
            sub = [start]
            visited = {start}
            while len(sub) < len(comp):
                cur = sub[-1]
                nxt = None
                for nb in adj_intra[cur]:
                    if nb not in visited:
                        nxt = nb
                        break
                if nxt is None:
                    break
                sub.append(nxt)
                visited.add(nxt)
            # Orient: if previous chain end is closer to sub[-1], reverse
            if chain:
                prev_end = chain[-1]
                d_fwd = _bfs_dist(prev_end, {sub[0]})
                d_rev = _bfs_dist(prev_end, {sub[-1]})
                if d_rev < d_fwd:
                    sub = sub[::-1]
            chain.extend(sub)
        return chain

    def _bond_len(sym1, sym2):
        pair = frozenset([sym1, sym2])
        if 'H' in pair:
            other = sym2 if sym1 == 'H' else sym1
            return {'C': 1.08, 'N': 1.01, 'O': 0.96, 'Si': 1.48, 'B': 1.19}.get(other, 1.08)
        if pair == frozenset(['Si', 'C']):
            return 1.87
        if pair == frozenset(['C', 'C']):
            return 1.50
        if pair == frozenset(['C', 'N']):
            return 1.47
        if pair == frozenset(['C', 'O']):
            return 1.43
        r1 = _COVALENT_RADII.get(sym1, 0.76)
        r2 = _COVALENT_RADII.get(sym2, 0.76)
        return r1 + r2

    def _rodrigues(points, center, axis, angle):
        k = axis / max(float(np.linalg.norm(axis)), 1e-12)
        c, s = np.cos(angle), np.sin(angle)
        out = []
        for p in points:
            v = p - center
            v_rot = v * c + np.cross(k, v) * s + k * float(np.dot(k, v)) * (1.0 - c)
            out.append(center + v_rot)
        return out

    # ---- Place metals ----
    if len(metal_indices) == 1:
        coords[metal_indices[0]] = [0.0, 0.0, 0.0]
        placed.add(metal_indices[0])
    else:
        for i, mi in enumerate(metal_indices):
            if i == 0:
                coords[mi] = [0.0, 0.0, 0.0]
            else:
                prev_mi = metal_indices[i - 1]
                si = mol.GetAtomWithIdx(mi).GetSymbol()
                sj = mol.GetAtomWithIdx(prev_mi).GetSymbol()
                bond = mol.GetBondBetweenAtoms(mi, prev_mi)
                d = (_COVALENT_RADII.get(si, 1.5) + _COVALENT_RADII.get(sj, 1.5) + 0.1
                     if bond is not None else 3.5)
                coords[mi] = coords[prev_mi] + np.array([d, 0.0, 0.0])
            placed.add(mi)

    # ---- Build global group list ----
    group_list: List[Tuple[int, List[int]]] = []
    for mi in metal_indices:
        for grp in by_metal[mi]:
            group_list.append((mi, grp))

    # ---- P4 Layer 2: DELFIN_HAPTO_SHARED_RING_FIX detection -------------
    # Detect clusters of eta-groups fused into one macrocycle (multi-atom
    # bridge paths, invisible to _find_ansa_bridges).  The legacy linked-
    # pair branch assigns every such pair the SAME [0,0,1] direction +
    # combined_centroid → exact overlap (MEWCIA: 4 pairs at d=0).  When the
    # flag is set, the cluster's groups are pre-assigned DISTINCT azimuthal
    # directions on a circle around the metal so the existing ring-builder
    # places each at a distinct centroid.  Default OFF → empty map →
    # byte-identical (legacy collapse preserved).
    _shared_ring_dirs: Dict[Tuple[int, int], "np.ndarray"] = {}
    _shared_ring_cluster_members: Dict[int, set] = {}
    if _delfin_env_int("DELFIN_HAPTO_SHARED_RING_FIX", 0):
        try:
            _clusters = _detect_shared_ring_eta_groups(mol, hapto_groups)
        except Exception as _sr_exc:
            logger.debug("shared-ring detection failed: %s", _sr_exc)
            _clusters = []
        # Map global hapto-group index -> (metal, local index within metal).
        _gidx_to_local: Dict[int, Tuple[int, int]] = {}
        _per_metal_counter: Dict[int, int] = {}
        for _g_idx, (_m_idx, _grp) in enumerate(hapto_groups):
            _loc = _per_metal_counter.get(_m_idx, 0)
            _gidx_to_local[_g_idx] = (_m_idx, _loc)
            _per_metal_counter[_m_idx] = _loc + 1
        for _cl in _clusters:
            _cl_metal = int(_cl["metal"])
            _cl_groups = list(_cl["groups"])  # global g_idx, sorted
            _n_cl = len(_cl_groups)
            if _n_cl < 2:
                continue
            _members = _shared_ring_cluster_members.setdefault(_cl_metal, set())
            for _k, _g_idx in enumerate(_cl_groups):
                _mloc = _gidx_to_local.get(_g_idx)
                if _mloc is None or _mloc[0] != _cl_metal:
                    continue
                _local = _mloc[1]
                _members.add(_local)
                # Distinct azimuth on a circle (2*pi*k/N) in the xy-plane,
                # deterministic by sorted group order.
                _az = 2.0 * np.pi * _k / _n_cl
                _shared_ring_dirs[(_cl_metal, _local)] = np.array([
                    np.cos(_az), np.sin(_az), 0.0,
                ])

    # ---- For each metal: build coordination sphere ----
    for mi in metal_indices:
        m_pos = coords[mi].copy()
        m_sym = mol.GetAtomWithIdx(mi).GetSymbol()
        groups = by_metal[mi]
        n_hapto = len(groups)
        n_sigma = len(sigma_donors.get(mi, []))
        n_total = n_hapto + n_sigma
        if n_total == 0:
            continue

        # Map local group indices -> global
        local_to_global = []
        for g_global, (gm, _) in enumerate(group_list):
            if gm == mi:
                local_to_global.append(g_global)

        # Detect ansa constraints for this metal
        ansa_angle = None
        ansa_local_i = None
        ansa_local_j = None
        ansa_bridge_indices: List[int] = []
        for (gi_g, gj_g), br_list in ansa_bridge_map.items():
            try:
                li = local_to_global.index(gi_g)
                lj = local_to_global.index(gj_g)
            except ValueError:
                continue
            ansa_local_i, ansa_local_j = li, lj
            ansa_bridge_indices = list(br_list)
            g_i, g_j = groups[li], groups[lj]
            eta_i, eta_j = len(g_i), len(g_j)
            r_i = _target_mc_dist(m_sym, eta_i)
            r_j = _target_mc_dist(m_sym, eta_j)
            R_i = 1.40 / (2.0 * np.sin(np.pi / max(eta_i, 3)))
            R_j = 1.40 / (2.0 * np.sin(np.pi / max(eta_j, 3)))
            r_avg = (r_i + r_j) / 2.0
            R_avg = (R_i + R_j) / 2.0
            # Use first bridge for angle computation
            br_sym = mol.GetAtomWithIdx(br_list[0]).GetSymbol()
            br_bond = {'Si': 1.87, 'C': 1.54, 'Ge': 1.94}.get(br_sym, 1.80)
            br_angle_deg = {'Si': 93.0, 'C': 109.5, 'Ge': 95.0}.get(br_sym, 100.0)
            d_target = 2.0 * br_bond * np.sin(np.radians(br_angle_deg) / 2.0)
            hyp = np.sqrt(r_avg ** 2 + R_avg ** 2)
            ratio = d_target / (2.0 * hyp)
            if abs(ratio) <= 1.0:
                psi = np.arctan2(R_avg, r_avg) + np.arcsin(ratio)
                ansa_angle = 2.0 * psi
            else:
                ansa_angle = np.radians(140.0)
            break

        # Generate direction vectors
        hapto_dirs: List[Optional[np.ndarray]] = [None] * n_hapto

        if ansa_angle is not None and ansa_local_i is not None:
            half = ansa_angle / 2.0
            d1 = np.array([0.0, np.sin(half), np.cos(half)])
            d2 = np.array([0.0, -np.sin(half), np.cos(half)])
            hapto_dirs[ansa_local_i] = d1 / float(np.linalg.norm(d1))
            hapto_dirs[ansa_local_j] = d2 / float(np.linalg.norm(d2))

        # P4 Layer 2: pre-assign DISTINCT azimuthal directions to shared-ring
        # cluster members (DELFIN_HAPTO_SHARED_RING_FIX).  This claims the
        # cluster groups before the linked-pair detection below (which only
        # fires on groups with hapto_dirs[li] is None), so the existing
        # ring-builder places each group at a distinct centroid on the
        # circle — eliminating the exact-overlap collapse.  Default OFF →
        # _shared_ring_dirs empty → no-op → byte-identical.
        _cl_members = _shared_ring_cluster_members.get(mi)
        if _cl_members:
            for _local in sorted(_cl_members):
                if 0 <= _local < n_hapto and hapto_dirs[_local] is None:
                    _d = _shared_ring_dirs.get((mi, _local))
                    if _d is not None:
                        _nrm = float(np.linalg.norm(_d))
                        if _nrm > 1e-9:
                            hapto_dirs[_local] = _d / _nrm

        # Linked groups: hapto pairs connected through short non-hapto paths
        # (e.g., COD eta2+eta2 bridged by CH2-CH2).
        # Detect and mark for combined rectangular placement.
        linked_pairs: List[Tuple[int, int]] = []  # (local_i, local_j)
        if n_hapto >= 2:
            from collections import deque
            linked_used: set = set()
            for li in range(n_hapto):
                if li in linked_used or hapto_dirs[li] is not None:
                    continue
                for lj in range(li + 1, n_hapto):
                    if lj in linked_used or hapto_dirs[lj] is not None:
                        continue
                    grp_i_set = set(groups[li])
                    grp_j_set = set(groups[lj])
                    linked = False
                    for start in groups[li]:
                        q = deque([(start, 0)])
                        seen = {start}
                        while q:
                            cur, d_bfs = q.popleft()
                            if cur in grp_j_set:
                                linked = True
                                break
                            if d_bfs >= 4:
                                continue
                            for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
                                ni = nb.GetIdx()
                                if ni not in seen and nb.GetSymbol() not in _METAL_SET:
                                    seen.add(ni)
                                    q.append((ni, d_bfs + 1))
                        if linked:
                            break
                    if linked:
                        linked_pairs.append((li, lj))
                        linked_used.add(li)
                        linked_used.add(lj)
                        # Assign same direction (will be placed as rectangle)
                        hapto_dirs[li] = np.array([0.0, 0.0, 1.0])
                        hapto_dirs[lj] = np.array([0.0, 0.0, 1.0])
                        break

        # ---- Assign ideal directions for unassigned hapto groups ----
        unassigned_hapto = [i for i in range(n_hapto) if hapto_dirs[i] is None]
        assigned_dirs = [hapto_dirs[i] for i in range(n_hapto)
                         if hapto_dirs[i] is not None]

        # In multi-metal systems, compute anti-M2 direction so hapto
        # groups point away from the other metal(s).
        _anti_other_metal = None
        if len(metal_indices) > 1:
            other_positions = [coords[mj] for mj in metal_indices if mj != mi
                               and mj in placed]
            if other_positions:
                avg_other = np.mean(other_positions, axis=0)
                away = m_pos - avg_other
                away_len = float(np.linalg.norm(away))
                if away_len > 1e-8:
                    _anti_other_metal = away / away_len

        if unassigned_hapto:
            if not assigned_dirs:
                # All hapto groups unassigned: use analytical placement
                if len(unassigned_hapto) == 1:
                    hapto_dirs[unassigned_hapto[0]] = (_anti_other_metal
                        if _anti_other_metal is not None
                        else np.array([0.0, 0.0, 1.0]))
                elif len(unassigned_hapto) == 2:
                    # Sandwich: ~175° apart (nearly anti-parallel).
                    # In multi-metal: orient perpendicular to M-M axis.
                    if _anti_other_metal is not None:
                        u_s, _ = _ortho_basis(_anti_other_metal)
                        hapto_dirs[unassigned_hapto[0]] = u_s
                        a175 = np.radians(175.0)
                        hapto_dirs[unassigned_hapto[1]] = (
                            np.cos(a175) * u_s + np.sin(a175)
                            * _anti_other_metal)
                    else:
                        hapto_dirs[unassigned_hapto[0]] = np.array([0.0, 0.0, 1.0])
                        a175 = np.radians(175.0)
                        hapto_dirs[unassigned_hapto[1]] = np.array(
                            [0.0, np.sin(a175), np.cos(a175)])
                else:
                    # Evenly distributed in a plane (120° for 3, 90° for 4, …)
                    for k, ui in enumerate(unassigned_hapto):
                        a = 2.0 * np.pi * k / len(unassigned_hapto)
                        hapto_dirs[ui] = np.array(
                            [np.cos(a), np.sin(a), 0.0])
            else:
                # Some already assigned (ansa/linked): place rest maximally
                # separated from all existing directions.
                used = list(assigned_dirs)
                candidates = _sphere_dirs(60)
                for ui in unassigned_hapto:
                    best_dir = None
                    best_max_dot = 1.0
                    for cd in candidates:
                        max_dot = max(float(np.dot(cd, ud)) for ud in used)
                        if max_dot < best_max_dot:
                            best_max_dot = max_dot
                            best_dir = cd
                    hapto_dirs[ui] = (best_dir if best_dir is not None
                                      else np.array([0.0, 0.0, 1.0]))
                    used.append(hapto_dirs[ui])

        # Fill sigma donor directions, maximally separated from hapto dirs.
        # Piano-stool special case: 1 hapto + n sigma → place sigma donors
        # analytically on a cone at ~125° from the hapto axis.
        sigma_dir_list = []
        if n_sigma > 0:
            all_used = [hd for hd in hapto_dirs if hd is not None]
            if n_hapto == 1 and len(all_used) == 1:
                hapto_ax = all_used[0] / max(float(np.linalg.norm(all_used[0])), 1e-12)
                cone_angle = np.radians(125.0)
                u_cone, v_cone = _ortho_basis(hapto_ax)
                for k in range(n_sigma):
                    phi = 2.0 * np.pi * k / n_sigma
                    d = (np.cos(cone_angle) * hapto_ax
                         + np.sin(cone_angle) * (np.cos(phi) * u_cone + np.sin(phi) * v_cone))
                    sigma_dir_list.append(d / max(float(np.linalg.norm(d)), 1e-12))
            else:
                candidates = _sphere_dirs(max(n_sigma + len(all_used) + 8, 60))
                candidates.sort(key=lambda cd: max(
                    float(np.dot(cd, ud)) for ud in all_used) if all_used else 0.0)
                for cd in candidates:
                    if len(sigma_dir_list) >= n_sigma:
                        break
                    ok = True
                    for ud in all_used + sigma_dir_list:
                        if float(np.dot(cd, ud)) > 0.50:
                            ok = False
                            break
                    if ok:
                        sigma_dir_list.append(cd)
            while len(sigma_dir_list) < n_sigma:
                sigma_dir_list.append(np.array([1.0, 0.0, 0.0]))

        for i in range(n_hapto):
            if hapto_dirs[i] is None:
                hapto_dirs[i] = np.array([0.0, 0.0, 1.0])

        # ---- Place linked pairs as rectangles first ----
        linked_placed_locals: set = set()
        for li, lj in linked_pairs:
            direction = hapto_dirs[li]
            grp_i, grp_j = groups[li], groups[lj]
            mc_i = _target_mc_dist(m_sym, len(grp_i))
            mc_j = _target_mc_dist(m_sym, len(grp_j))
            mc_avg = (mc_i + mc_j) / 2.0
            combined_centroid = m_pos + direction * mc_avg
            u_lnk, v_lnk = _ortho_basis(direction)

            # Find which atoms bridge the two groups, to orient the rectangle
            # Chain order: build chain for combined atoms through mol graph
            all_linked = list(grp_i) + list(grp_j)
            chain_i = _build_chain_order(grp_i, mol)
            chain_j = _build_chain_order(grp_j, mol)

            # Orient chains: find which ends face each other (BFS shortest)
            from collections import deque
            def _bfs_d(src, tgt):
                q = deque([(src, 0)])
                seen = {src}
                while q:
                    cur, dd = q.popleft()
                    if cur == tgt:
                        return dd
                    if dd >= 6:
                        return 99
                    for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
                        ni = nb.GetIdx()
                        if ni not in seen and nb.GetSymbol() not in _METAL_SET:
                            seen.add(ni)
                            q.append((ni, dd + 1))
                return 99

            # Find the pair of endpoints with shortest BFS distance
            best_di, best_dj = 0, 0
            best_bfs = 99
            for di_end in [0, -1]:
                for dj_end in [0, -1]:
                    bd = _bfs_d(chain_i[di_end], chain_j[dj_end])
                    if bd < best_bfs:
                        best_bfs = bd
                        best_di, best_dj = di_end, dj_end
            # Orient chains so facing ends are the last atoms
            if best_di == 0:
                chain_i = chain_i[::-1]
            if best_dj == 0:
                chain_j = chain_j[::-1]

            # Place as rectangle: two parallel rows separated by bridge gap
            # Bridge gap ≈ n_bridge_atoms * 1.54A with 109.5° angles
            bridge_gap = max(best_bfs * 1.00, 2.2)  # empirical
            half_gap = bridge_gap / 2.0
            cc = 1.40  # C-C distance within each pair

            for k, atom_idx in enumerate(chain_i):
                offset_v = (k - (len(chain_i) - 1) / 2.0) * cc
                coords[atom_idx] = combined_centroid + half_gap * u_lnk + offset_v * v_lnk
                placed.add(atom_idx)
            for k, atom_idx in enumerate(chain_j):
                offset_v = (k - (len(chain_j) - 1) / 2.0) * cc
                coords[atom_idx] = combined_centroid - half_gap * u_lnk + offset_v * v_lnk
                placed.add(atom_idx)
            linked_placed_locals.add(li)
            linked_placed_locals.add(lj)

        # ---- Build hapto rings/chains on sphere ----
        for g_local, grp in enumerate(groups):
            if g_local in linked_placed_locals:
                continue  # already placed as linked rectangle
            direction = hapto_dirs[g_local]
            eta = len(grp)
            mc_dist = _target_mc_dist(m_sym, eta)
            centroid = m_pos + direction * mc_dist

            if eta >= 3:
                R_ring = 1.40 / (2.0 * np.sin(np.pi / eta))
            else:
                R_ring = 0.69

            u, v = _ortho_basis(direction)
            ordered = _ring_traversal(grp)

            # Check if group is cyclic (all atoms have >=2 intra-group bonds)
            grp_set = set(grp)
            is_cyclic = True
            for a_idx in grp:
                intra = sum(1 for nb in mol.GetAtomWithIdx(a_idx).GetNeighbors()
                            if nb.GetIdx() in grp_set)
                if intra < 2:
                    is_cyclic = False
                    break

            if eta == 2:
                coords[ordered[0]] = centroid + R_ring * u
                coords[ordered[1]] = centroid - R_ring * u
            elif is_cyclic:
                # Regular polygon (closed ring)
                for k, atom_idx in enumerate(ordered):
                    angle = 2.0 * np.pi * k / eta
                    coords[atom_idx] = centroid + R_ring * (
                        np.cos(angle) * u + np.sin(angle) * v
                    )
            else:
                # Non-cyclic group: check for disconnected components
                # (e.g., COD eta4 = two eta2 pairs bridged by CH2-CH2)
                chain = _build_chain_order(grp, mol)
                # Find connected components within the hapto group
                comp_list = []
                cur_comp = [chain[0]]
                for k in range(1, len(chain)):
                    # Check if chain[k] bonds to chain[k-1] within grp
                    bonded = False
                    for nb in mol.GetAtomWithIdx(chain[k]).GetNeighbors():
                        if nb.GetIdx() == chain[k - 1]:
                            bonded = True
                            break
                    if bonded:
                        cur_comp.append(chain[k])
                    else:
                        comp_list.append(cur_comp)
                        cur_comp = [chain[k]]
                comp_list.append(cur_comp)

                if len(comp_list) >= 2:
                    # Disconnected components: place as parallel pairs
                    # on a rectangle to keep both gaps bridgeable
                    n_comp = len(comp_list)
                    for ci, comp in enumerate(comp_list):
                        # Angle offset for this component on the circle
                        comp_center_angle = 2.0 * np.pi * ci / n_comp
                        cc_half = 1.40 * (len(comp) - 1) / 2.0
                        for k, atom_idx in enumerate(comp):
                            # Spread atoms within component along v axis
                            offset = (k - (len(comp) - 1) / 2.0) * 1.40
                            pos = centroid + (
                                R_ring * np.cos(comp_center_angle) * u
                                + offset * v
                            )
                            coords[atom_idx] = pos
                else:
                    # Single connected chain: place along arc
                    step_angle = 2.0 * np.arcsin(
                        min(1.40 / (2.0 * R_ring), 1.0))
                    total_span = step_angle * (eta - 1)
                    start_angle = -total_span / 2.0
                    for k, atom_idx in enumerate(chain):
                        angle = start_angle + k * step_angle
                        coords[atom_idx] = centroid + R_ring * (
                            np.cos(angle) * u + np.sin(angle) * v
                        )
            for atom_idx in grp:
                placed.add(atom_idx)

        # ---- Orient ansa-bridged rings so bridge carbons face each other ----
        if ansa_bridge_indices and ansa_local_i is not None:
            grp_i = groups[ansa_local_i]
            grp_j = groups[ansa_local_j]
            # Collect all neighbors of ALL bridge atoms
            br_nbrs: set = set()
            for _bi in ansa_bridge_indices:
                for _bn in mol.GetAtomWithIdx(_bi).GetNeighbors():
                    br_nbrs.add(_bn.GetIdx())
            for grp_src, dir_src, grp_tgt_local in [
                (grp_i, hapto_dirs[ansa_local_i], ansa_local_j),
                (grp_j, hapto_dirs[ansa_local_j], ansa_local_i),
            ]:
                c_src = next((a for a in grp_src if a in br_nbrs), None)
                if c_src is None:
                    continue
                cent_src = m_pos + dir_src * _target_mc_dist(m_sym, len(grp_src))
                dir_tgt = hapto_dirs[grp_tgt_local]
                cent_tgt = m_pos + dir_tgt * _target_mc_dist(
                    m_sym, len(groups[grp_tgt_local])
                )
                toward = cent_tgt - cent_src
                toward_proj = toward - float(np.dot(toward, dir_src)) * dir_src
                tp_len = float(np.linalg.norm(toward_proj))
                if tp_len < 1e-8:
                    continue
                toward_proj = toward_proj / tp_len
                vc = coords[c_src] - cent_src
                vc_proj = vc - float(np.dot(vc, dir_src)) * dir_src
                vp_len = float(np.linalg.norm(vc_proj))
                if vp_len < 1e-8:
                    continue
                vc_proj = vc_proj / vp_len
                cos_r = float(np.clip(np.dot(vc_proj, toward_proj), -1.0, 1.0))
                sin_r = float(np.dot(np.cross(vc_proj, toward_proj), dir_src))
                rot_angle = float(np.arctan2(sin_r, cos_r))
                if abs(rot_angle) > 0.01:
                    pts = [coords[a].copy() for a in grp_src]
                    pts_rot = _rodrigues(pts, cent_src, dir_src, rot_angle)
                    for a_idx, new_p in zip(grp_src, pts_rot):
                        coords[a_idx] = new_p

        # ---- Place sigma donors ----
        for s_idx, donor_idx in enumerate(sigma_donors.get(mi, [])):
            if donor_idx in placed:
                continue
            direction = (sigma_dir_list[s_idx] if s_idx < len(sigma_dir_list)
                         else np.array([1.0, 0.0, 0.0]))
            donor_sym = mol.GetAtomWithIdx(donor_idx).GetSymbol()
            bl = _get_ml_bond_length(m_sym, donor_sym)
            coords[donor_idx] = m_pos + direction * bl
            placed.add(donor_idx)

    # ---- Place substituents on hapto ring atoms ----
    # Welle-3 T7.1 (2026-05-15): multi-substituent tilt direction fix.
    # Legacy formula (len(subs)>=2): tilt=(s_k-(len-1)/2)*1.2 rad splays
    # alternating substituents ABOVE and BELOW the ring plane.  For sp3
    # endpoints of η4-dienes this is geometrically correct (axial/equa-
    # torial H), but for ring-fused sp2 carbons (indenyl/naphthyl-type
    # η-ligands, η3-arene-derived allyls) the s_k=0 substituent gets
    # tilted -34° TOWARD the metal, dragging the pendant ring into the
    # coordination sphere.  Audit on 100 hapto SMILES (Welle-3 T7.1):
    # η3: 100%, η4: 64%, η5: 17%, η6: 12% high-OOP (>30°) groups, with
    # 20/44 cases driven by this multi-sub formula pulling the second
    # sub into the M-cone.  User pattern: "mostly fails at the substituents
    # on the hapto directions and the hapto-ring connections".
    #
    # Fix (DELFIN_HAPTO_SUB_FIX=1, default 0):
    #   * Replace signed tilt by ABS tilt (cos(|t|), sin(|t|)) so both
    #     subs lift AWAY from metal (no negative-z lobe).
    #   * Add azimuthal spread via tangent vector (perpendicular to
    #     outward and normal) so the subs are not coplanar in the
    #     M-outward plane (would clash).
    #   * Single-sub case unchanged (legacy +0.3*normal is fine).
    _sub_fix = _delfin_env_int("DELFIN_HAPTO_SUB_FIX", 0) == 1
    for metal_idx, grp in hapto_groups:
        ring_set = set(grp)
        m_pos = coords[metal_idx]
        centroid = np.mean([coords[a] for a in grp], axis=0)
        mc_dir = centroid - m_pos
        mc_len = float(np.linalg.norm(mc_dir))
        normal = mc_dir / mc_len if mc_len > 1e-8 else np.array([0.0, 0.0, 1.0])

        for atom_idx in grp:
            ring_pos = coords[atom_idx]
            outward = ring_pos - centroid
            outward_len = float(np.linalg.norm(outward))
            outward_unit = (outward / outward_len if outward_len > 1e-8
                            else np.array([1.0, 0.0, 0.0]))
            subs = []
            for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                ni = nbr.GetIdx()
                if ni in ring_set or ni == metal_idx or ni in placed:
                    continue
                if ni in all_bridge_atoms:
                    continue
                subs.append(nbr)

            # Tangent: perp to outward + normal (in ring plane, perp to outward)
            _tangent = np.cross(normal, outward_unit)
            _tl = float(np.linalg.norm(_tangent))
            _tangent = (_tangent / _tl if _tl > 1e-8
                        else np.array([0.0, 1.0, 0.0]))
            for s_k, nbr in enumerate(subs):
                ni = nbr.GetIdx()
                bl = _bond_len(
                    mol.GetAtomWithIdx(atom_idx).GetSymbol(), nbr.GetSymbol()
                )
                if len(subs) == 1:
                    direction = outward_unit + 0.3 * normal
                elif _sub_fix:
                    # T7.1: all subs lift AWAY from metal (abs tilt),
                    # spread tangentially in ring plane to avoid clash.
                    tilt = 0.55  # ~31.5° away-from-M for all subs
                    azi = (s_k - (len(subs) - 1) / 2.0) * np.radians(25.0)
                    direction = (outward_unit * np.cos(tilt) * np.cos(azi)
                                 + _tangent * np.cos(tilt) * np.sin(azi)
                                 + normal * np.sin(tilt))
                else:
                    tilt = (s_k - (len(subs) - 1) / 2.0) * 1.2
                    direction = outward_unit * np.cos(tilt) + normal * np.sin(tilt)
                d_norm = float(np.linalg.norm(direction))
                if d_norm > 1e-8:
                    direction = direction / d_norm
                coords[ni] = ring_pos + bl * direction
                placed.add(ni)

    # ---- Place ALL ansa bridge atoms ----
    for (gi, gj), bridge_list in ansa_bridge_map.items():
        grp_i = hapto_groups[gi][1]
        grp_j = hapto_groups[gj][1]
        metal_idx = hapto_groups[gi][0]

        for bridge_idx in bridge_list:
            if bridge_idx in placed:
                continue
            bridge_atom = mol.GetAtomWithIdx(bridge_idx)
            bridge_sym = bridge_atom.GetSymbol()
            br_nbrs = set(n.GetIdx() for n in bridge_atom.GetNeighbors())

            c_i = next((a for a in grp_i if a in br_nbrs), None)
            c_j = next((a for a in grp_j if a in br_nbrs), None)
            if c_i is not None and c_j is not None:
                pos_a, pos_b = coords[c_i], coords[c_j]
                mid = (pos_a + pos_b) / 2.0
                bl = _bond_len(bridge_sym, 'C')
                d_ab = float(np.linalg.norm(pos_a - pos_b))
                h_sq = bl ** 2 - (d_ab / 2.0) ** 2
                h = float(np.sqrt(max(h_sq, 0.01)))
                away = mid - coords[metal_idx]
                a_norm = float(np.linalg.norm(away))
                away = away / a_norm if a_norm > 1e-8 else np.array([0.0, 1.0, 0.0])
                coords[bridge_idx] = mid + h * away
                placed.add(bridge_idx)

                # Place bridge substituents using robust local frame
                ca_cb = pos_b - pos_a
                x_ax = ca_cb / max(d_ab, 1e-8)
                z_ax = away
                y_ax = np.cross(z_ax, x_ax)
                y_norm = float(np.linalg.norm(y_ax))
                if y_norm < 1e-8:
                    # Degenerate: z_ax parallel to x_ax, pick perpendicular
                    ref = (np.array([1.0, 0.0, 0.0]) if abs(z_ax[0]) < 0.9
                           else np.array([0.0, 1.0, 0.0]))
                    y_ax = np.cross(z_ax, ref)
                    y_norm = float(np.linalg.norm(y_ax))
                y_ax = y_ax / max(y_norm, 1e-12)

                sub_count = 0
                for pnbr in bridge_atom.GetNeighbors():
                    pni = pnbr.GetIdx()
                    if pni in placed or pni in all_hapto_atoms:
                        continue
                    sub_bl = _bond_len(bridge_sym, pnbr.GetSymbol())
                    sub_angle = (sub_count - 0.5) * 2.1
                    direction = np.cos(sub_angle) * y_ax + np.sin(sub_angle) * z_ax
                    d_norm = float(np.linalg.norm(direction))
                    if d_norm > 1e-8:
                        direction = direction / d_norm
                    coords[pni] = coords[bridge_idx] + sub_bl * direction
                    placed.add(pni)
                    sub_count += 1
                sub_count += 1

    # ---- Pre-compute ring info for inline ring placement during BFS ----
    try:
        Chem.FastFindRings(mol)
        _ri = mol.GetRingInfo()
        _all_rings = list(_ri.AtomRings())
    except Exception:
        _all_rings = []
    _organic_rings = []
    for ring in _all_rings:
        ring_set = set(ring)
        # Iter-8.10 (2026-05-11): skip rings that have ANY overlap with
        # hapto atoms (not just rings ENTIRELY inside hapto set). Previous
        # condition `<= all_hapto_atoms` left mixed rings (hapto + non-hapto
        # members) to be polygon-placed by _place_ring_inline, which
        # OVERWRITES the already-placed hapto coordinates with a generic
        # n-gon. Caused GIQRAA Fe CN10 (2× η4-diene): only 5/8 Fe-C bonds
        # at correct 2.011Å, 3 hapto-Cs scattered at 2.78-3.0Å. Same
        # ZOQTEE Ru CN12: 9/12 correct.  Fix (Subagent-validated):
        #   GIQRAA: 5/8 → 8/8 hapto Fe-C correct
        #   ZOQTEE: 9/12 → 12/12 hapto Ru-C correct
        # Bug pre-dates HEAD (81f8a1f also has it). No env-flag needed —
        # the previous behavior is objectively broken for multi-fragment
        # hapto. Non-hapto ring atoms fall through to BFS-VSEPR placement.
        if ring_set & all_hapto_atoms:
            continue
        if any(mol.GetAtomWithIdx(a).GetSymbol() in _METAL_SET for a in ring):
            continue
        _organic_rings.append(ring)
    _atom_to_rings: Dict[int, List[int]] = {}
    for ri_idx, ring in enumerate(_organic_rings):
        for a in ring:
            _atom_to_rings.setdefault(a, []).append(ri_idx)
    _ring_done: set = set()

    def _place_ring_inline(ring_tuple, anchor_idx):
        ring_set = set(ring_tuple)
        rsize = len(ring_tuple)
        cc_ring = 1.40 if rsize <= 6 else 1.50
        radius = cc_ring / (2.0 * np.sin(np.pi / max(rsize, 3)))

        anchor_pos = coords[anchor_idx]
        parent_of_anchor = None
        for nbr in mol.GetAtomWithIdx(anchor_idx).GetNeighbors():
            ni = nbr.GetIdx()
            if ni in placed and ni not in ring_set:
                parent_of_anchor = ni
                break
        if parent_of_anchor is not None:
            bond_dir = anchor_pos - coords[parent_of_anchor]
            bd_len = float(np.linalg.norm(bond_dir))
            if bd_len > 1e-8:
                bond_dir = bond_dir / bd_len
            else:
                bond_dir = np.array([1.0, 0.0, 0.0])
        else:
            bond_dir = np.array([1.0, 0.0, 0.0])

        # Iter-8.8 (2026-05-11): Per-ring chemistry-correct dispatch.
        #
        # Two geometric regimes, decided by ring membership in
        # all_hapto_atoms:
        #
        # (A) Hapto-coordinated ring (Cp / arene / diene η-bonded to
        #     metal): bond_dir is the metal→ring axis = ring NORMAL.
        #     Ring plane is perpendicular to bond_dir.  81f8a1f formula.
        #     (User pattern: hapto fragments destroyed in HEAD post-Iter-5)
        #
        # (B) Pendant organic ring (phenyl substituent on σ-donor,
        #     fused-aromatic side-arm): bond_dir lies IN the ring plane
        #     (the σ-bond is one ring edge).  HEAD/Iter-5 formula.
        #
        # env-flag DELFIN_HAPTO_LEGACY_RING_PLANE: default 0 (HEAD/Iter-5
        # formula active across all rings).  Default REVERTED 2026-05-11
        # after smoke500 with default=1 showed sigma -2.8pp topo
        # regression with no measurable hapto/multi-hapto compensation;
        # leaves dispatch as opt-in env=1 for per-SMILES experiments.
        _legacy_ring_plane = bool(
            _delfin_env_int("DELFIN_HAPTO_LEGACY_RING_PLANE", 0)
        )
        _ring_is_hapto = (
            _legacy_ring_plane
            and any(a in all_hapto_atoms for a in ring_tuple)
        )

        if _ring_is_hapto:
            # 81f8a1f hapto-ring: bond_dir is ring NORMAL.  Polygon centre
            # at radius along bond_dir from anchor.  Anchor preserved at
            # angle = anchor_offset.
            u_r, v_r = _ortho_basis(bond_dir)
            ring_center = anchor_pos + bond_dir * radius
        else:
            # HEAD/Iter-5 pendant ring: bond_dir IN ring plane.  Polygon
            # centre at radius along bond_dir from anchor.  Anchor at
            # angle=0 in (-bond_dir, perp_in_plane) basis.
            ring_center = anchor_pos + bond_dir * radius
            perp_in_plane, _ = _ortho_basis(bond_dir)

        ring_adj = {a: [] for a in ring_tuple}
        for a in ring_tuple:
            for nbr in mol.GetAtomWithIdx(a).GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    ring_adj[a].append(nbr.GetIdx())
        chain = [anchor_idx]
        visited_r = {anchor_idx}
        while len(chain) < rsize:
            cur = chain[-1]
            nxt = None
            for nb in ring_adj[cur]:
                if nb not in visited_r:
                    nxt = nb
                    break
            if nxt is None:
                break
            chain.append(nxt)
            visited_r.add(nxt)
        if len(chain) != rsize:
            return False

        if _ring_is_hapto:
            # 81f8a1f formula: anchor preserved at its angle in (u_r, v_r)
            anchor_offset = float(np.arctan2(
                float(np.dot(anchor_pos - ring_center, v_r)),
                float(np.dot(anchor_pos - ring_center, u_r)),
            ))
            for k, atom_idx in enumerate(chain):
                angle = anchor_offset + 2.0 * np.pi * k / rsize
                coords[atom_idx] = ring_center + radius * (
                    np.cos(angle) * u_r + np.sin(angle) * v_r
                )
                placed.add(atom_idx)
        else:
            for k, atom_idx in enumerate(chain):
                angle = 2.0 * np.pi * k / rsize
                coords[atom_idx] = ring_center + radius * (
                    np.cos(angle) * (-bond_dir) + np.sin(angle) * perp_in_plane
                )
                placed.add(atom_idx)

        # Iter-22 (2026-05-20): fused-aromatic coplanar enforcement.
        # If this ring is ortho-fused (shares >=2 atoms = a bond-edge) to an
        # already-placed organic ring, project its non-shared atoms onto the
        # partner ring's plane so the fused system stays coplanar.  Prevents
        # the post-UFF ring-twist that derive_xyz_bonds (organic_tol 0.40)
        # mis-reads as spurious C-C bonds — the hapto %topo bottleneck per
        # forensik_hapto_class_pre_uff_2026_05_18 (74.1% broken via secondary
        # fused-aromatic ligand).  In-plane n-gon arrangement is preserved;
        # only the out-of-plane component is removed (UFF relaxes the rest).
        # Graph-only, default-ON hapto+multi_hapto, bit-exact for sigma.
        if _class_conditional_flag(
            "DELFIN_5F_FUSED_AROMATIC_COPLANAR", mol, default=0,
            default_classes=["hapto", "multi_hapto"],
        ):
            try:
                _best = None
                for _other in _organic_rings:
                    _os = set(_other)
                    if _os == ring_set or not _os.issubset(placed):
                        continue
                    _sh = ring_set & _os
                    if len(_sh) < 2:
                        continue
                    if _best is None or len(_sh) > len(_best[0]):
                        _pts = np.array([coords[a] for a in _other], dtype=float)
                        _cen = _pts.mean(axis=0)
                        _, _, _vh = np.linalg.svd(_pts - _cen)
                        _nl = float(np.linalg.norm(_vh[2]))
                        if _nl > 1e-9:
                            _best = (sorted(_sh), _vh[2] / _nl)
                if _best is not None:
                    _shared_atoms, _normal = _best
                    _plane_pt = coords[_shared_atoms[0]].astype(float)
                    for _ai in chain:
                        if _ai in _shared_atoms:
                            continue
                        _d = float(np.dot(coords[_ai] - _plane_pt, _normal))
                        coords[_ai] = coords[_ai] - _d * _normal
            except Exception:
                pass
        return True

    # ---- BFS propagate remaining atoms with VSEPR local geometry ----
    for _iteration in range(n_atoms * 3):
        progress = False
        for atom in mol.GetAtoms():
            ai = atom.GetIdx()
            if ai in placed:
                continue
            parent = None
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in placed:
                    parent = nbr
                    break
            if parent is None:
                continue

            # If this atom belongs to an unplaced organic ring, place
            # the entire ring analytically instead of via BFS.
            ring_handled = False
            for _ri_idx in _atom_to_rings.get(ai, []):
                if _ri_idx in _ring_done:
                    continue
                _ring_t = _organic_rings[_ri_idx]
                _ring_s = set(_ring_t)
                _ring_anchors = [a for a in _ring_t if a in placed]
                if _ring_anchors:
                    if _place_ring_inline(_ring_t, _ring_anchors[0]):
                        _ring_done.add(_ri_idx)
                        ring_handled = True
                        progress = True
                        break
            if ring_handled:
                continue

            pi = parent.GetIdx()
            parent_pos = coords[pi]
            child_sym = atom.GetSymbol()
            parent_sym = parent.GetSymbol()

            if parent_sym in _METAL_SET:
                bl = _get_ml_bond_length(parent_sym, child_sym)
            elif child_sym in _METAL_SET:
                bl = _get_ml_bond_length(child_sym, parent_sym)
            else:
                bl = _bond_len(parent_sym, child_sym)

            # VSEPR: direction away from already-placed neighbors of parent
            used_dirs = []
            for pnbr in parent.GetNeighbors():
                if pnbr.GetIdx() in placed and pnbr.GetIdx() != ai:
                    d = coords[pnbr.GetIdx()] - parent_pos
                    d_len = float(np.linalg.norm(d))
                    if d_len > 1e-8:
                        used_dirs.append(d / d_len)

            # Metal avoidance: if child is NOT bonded to any metal,
            # add metal directions as repulsive so BFS pushes away.
            _child_bonded_metals = set()
            for nbr in atom.GetNeighbors():
                if nbr.GetSymbol() in _METAL_SET:
                    _child_bonded_metals.add(nbr.GetIdx())
            for mi in metal_indices:
                if mi in _child_bonded_metals:
                    continue
                m_to_p = parent_pos - coords[mi]
                m_to_p_len = float(np.linalg.norm(m_to_p))
                if m_to_p_len < 4.0 and m_to_p_len > 1e-8:
                    used_dirs.append(-m_to_p / m_to_p_len)

            if len(used_dirs) == 0:
                direction = np.array([0.0, 0.0, 1.0])
            elif len(used_dirs) == 1:
                d0 = used_dirs[0]
                ref = (np.array([1.0, 0.0, 0.0]) if abs(d0[0]) < 0.9
                       else np.array([0.0, 1.0, 0.0]))
                perp = np.cross(d0, ref)
                perp = perp / max(float(np.linalg.norm(perp)), 1e-12)
                direction = (-d0 * np.cos(np.radians(70.5))
                             + perp * np.sin(np.radians(70.5)))
            elif len(used_dirs) == 2:
                avg = (used_dirs[0] + used_dirs[1]) / 2.0
                avg_len = float(np.linalg.norm(avg))
                if avg_len > 1e-8:
                    direction = -avg / avg_len
                else:
                    perp = np.cross(used_dirs[0], used_dirs[1])
                    p_len = float(np.linalg.norm(perp))
                    direction = (perp / p_len if p_len > 1e-8
                                 else np.array([0.0, 0.0, 1.0]))
            else:
                avg = sum(used_dirs) / len(used_dirs)
                avg_len = float(np.linalg.norm(avg))
                direction = (-avg / avg_len if avg_len > 1e-8
                             else np.array([0.0, 0.0, 1.0]))

            d_norm = float(np.linalg.norm(direction))
            if d_norm > 1e-8:
                direction = direction / d_norm

            coords[ai] = parent_pos + bl * direction
            placed.add(ai)
            progress = True

        if not progress:
            break

    # Handle orphan atoms
    rng = np.random.default_rng(42)
    for ai in range(n_atoms):
        if ai not in placed:
            coords[ai] = rng.standard_normal(3) * 3.0

    # ---- Bond-length relaxation (fix topology) ----
    # BFS only uses one parent per atom; bonds to other placed atoms may be
    # stretched.  Iteratively correct all bond lengths toward their targets.
    bond_targets: List[Tuple[int, int, float]] = []
    for bond in mol.GetBonds():
        bi = bond.GetBeginAtomIdx()
        bj = bond.GetEndAtomIdx()
        si = mol.GetAtomWithIdx(bi).GetSymbol()
        sj = mol.GetAtomWithIdx(bj).GetSymbol()
        if si in _METAL_SET or sj in _METAL_SET:
            continue  # skip metal bonds - already placed correctly
        target_bl = _bond_len(si, sj)
        bond_targets.append((bi, bj, target_bl))

    # Direct neighbours of any metal (σ-donors) — placed at the correct
    # M-D bond length by the σ-donor placement block.  Their downstream
    # ligand bonds (e.g. Si-Cl on a SiCl3 σ-donor, or N-C inside an
    # imidazole ring) are spring-constrained while the M-D bond itself
    # is not (skipped at L12554), so the spring forces drag the donor
    # away from the metal toward its other ligand atoms.  Freezing
    # σ-donors here keeps every M-D distance pinned at the cone target,
    # generalising the hapto-atom freeze to all coordination donors.
    metal_donors_frozen: set = set()
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        for nbr in atom.GetNeighbors():
            ni = nbr.GetIdx()
            if nbr.GetSymbol() in _METAL_SET:
                continue
            if ni in all_hapto_atoms:
                continue
            metal_donors_frozen.add(ni)

    # Welle-5l Track-4 (2026-05-18): conditional revert of the d0c345c
    # σ-donor HARD-freeze.  When the molecule carries an sp3-C donor or
    # a terminal CH₃ on a σ-donor (universal RDKit graph features), the
    # champion e6761e4 mechanism applies a *partial* relaxation
    # ``w = 0.3`` instead of fully freezing — recovering CH₃ umbrella
    # realism + sp3-C tetrahedral geometry on 29-Ni / WUXQAK / ACPICF
    # patterns while preserving the topology + M-D win on the rest of
    # the pool.  Default-OFF: env-flag DELFIN_5L_T4_SOFT_DONOR_FEATURE_GATE
    # must be set to ``1`` for any behavioural change.  See
    # ``delfin/_e6761e4_soft_donor.py`` for the feature detectors and
    # ``iters/welle5l_RETRO_3_hidden_wins_2026_05_18.md`` Section 2 for
    # the master-rank analysis (29/136 metric wins for e6761e4).
    try:
        from delfin.manta._e6761e4_soft_donor import (
            apply_soft_relaxation as _t4_apply_soft,
            donor_weight_for_atom as _t4_donor_weight,
        )
        _t4_soft_active = bool(_t4_apply_soft(mol))
    except Exception:
        _t4_soft_active = False
        _t4_donor_weight = None  # type: ignore[assignment]

    # Spring-based relaxation: accumulate all bond forces, then apply
    for _relax_pass in range(120):
        forces = np.zeros_like(coords)
        max_err = 0.0
        for bi, bj, target_bl in bond_targets:
            diff = coords[bj] - coords[bi]
            d = float(np.linalg.norm(diff))
            err = abs(d - target_bl)
            if err < 0.01:
                continue
            if d < 1e-8:
                diff = rng.standard_normal(3)
                d = float(np.linalg.norm(diff))
            unit = diff / d
            # Spring force proportional to displacement
            force = 0.3 * (d - target_bl) * unit
            forces[bi] += force
            forces[bj] -= force
            if err > max_err:
                max_err = err
        if max_err < 0.05:
            break
        # Apply forces — metals and hapto-ring atoms are frozen.  σ-donors
        # are frozen by default (post-d0c345c) UNLESS the Welle-5l T4
        # feature gate is active for this molecule, in which case they
        # receive a partial weight (default w = 0.3, env-tunable via
        # DELFIN_5L_T4_SOFT_DONOR_WEIGHT).  Default-OFF behaviour =
        # bit-exact identical to pre-T4.
        for ai in range(n_atoms):
            if mol.GetAtomWithIdx(ai).GetSymbol() in _METAL_SET:
                continue
            if ai in all_hapto_atoms:
                continue
            if ai in metal_donors_frozen:
                if not _t4_soft_active or _t4_donor_weight is None:
                    continue
                _w_ai = _t4_donor_weight(mol, ai)
                if _w_ai <= 0.0:
                    continue
                coords[ai] = coords[ai] + forces[ai] * _w_ai
                continue
            coords[ai] = coords[ai] + forces[ai]

    # ---- Clash resolution (push apart, preserve intra-ring geometry) ----
    # Build metal bond set for metal proximity check
    metal_bonded: set = set()
    for mi in metal_indices:
        for nbr in mol.GetAtomWithIdx(mi).GetNeighbors():
            metal_bonded.add((mi, nbr.GetIdx()))
            metal_bonded.add((nbr.GetIdx(), mi))

    rng_clash = np.random.default_rng(123)
    for _pass in range(30):
        moved = False
        for i in range(n_atoms):
            is_metal_i = mol.GetAtomWithIdx(i).GetSymbol() in _METAL_SET
            for j in range(i + 1, n_atoms):
                is_metal_j = mol.GetAtomWithIdx(j).GetSymbol() in _METAL_SET
                if is_metal_i and is_metal_j:
                    continue
                diff = coords[j] - coords[i]
                d = float(np.linalg.norm(diff))
                # Metal proximity: non-bonded atoms too close to metal
                if (is_metal_i or is_metal_j) and (i, j) not in metal_bonded:
                    min_ml = 2.8  # minimum non-bonded M-X distance
                    if d < min_ml:
                        if d < 1e-8:
                            direction = rng_clash.standard_normal(3)
                            direction /= max(np.linalg.norm(direction), 1e-12)
                            push = 0.5 * direction
                        else:
                            push = 0.8 * (min_ml - d) / d * diff
                        non_metal = j if is_metal_i else i
                        metal_atom = i if is_metal_i else j
                        sign = 1.0 if non_metal == j else -1.0
                        # Allow push if atom is not hapto, OR if it's hapto
                        # but not bonded to THIS metal (bimetallic case).
                        is_own_hapto = (non_metal in all_hapto_atoms
                                        and (metal_atom, non_metal) in metal_bonded)
                        if not is_own_hapto:
                            coords[non_metal] = coords[non_metal] + sign * push
                            moved = True
                    continue
                if is_metal_i or is_metal_j:
                    continue
                if d < 0.8:
                    if d < 1e-8:
                        direction = rng_clash.standard_normal(3)
                        direction /= max(np.linalg.norm(direction), 1e-12)
                        push = 0.4 * direction
                    else:
                        push = 0.5 * (0.8 - d) / d * diff
                    gi = atom_to_hapto_group.get(i, -1)
                    gj = atom_to_hapto_group.get(j, -2)
                    same_group = (gi == gj and gi >= 0)
                    if not same_group:
                        can_move_i = i not in all_hapto_atoms
                        can_move_j = j not in all_hapto_atoms
                        if can_move_i:
                            coords[i] = coords[i] - push * 0.5
                            moved = True
                        if can_move_j:
                            coords[j] = coords[j] + push * 0.5
                            moved = True
        if not moved:
            break

    # ---- Rigid-body separation for multi-metal systems ----
    # If hapto atoms of one metal intrude into another metal's coordination
    # sphere, translate the entire hapto fragment (metal + all its bonded
    # atoms) as a rigid body to increase M-M distance.
    if len(metal_indices) > 1:
        min_nonbonded_ml = 2.8
        for _sep_pass in range(5):
            any_moved = False
            for mi in metal_indices:
                mi_hapto = set()
                for grp in by_metal.get(mi, []):
                    mi_hapto.update(grp)
                if not mi_hapto:
                    continue
                for mj in metal_indices:
                    if mj == mi:
                        continue
                    mj_pos = coords[mj]
                    worst_intrusion = 0.0
                    push_dir = np.zeros(3)
                    for ha in mi_hapto:
                        d = float(np.linalg.norm(coords[ha] - mj_pos))
                        if d < min_nonbonded_ml and (mj, ha) not in metal_bonded:
                            intrusion = min_nonbonded_ml - d
                            if intrusion > worst_intrusion:
                                worst_intrusion = intrusion
                                v = coords[ha] - mj_pos
                                vl = float(np.linalg.norm(v))
                                push_dir = v / vl if vl > 1e-8 else np.array([1, 0, 0])
                    if worst_intrusion > 0.05:
                        shift = push_dir * worst_intrusion * 1.2
                        mi_frag = {mi} | mi_hapto
                        for nbr in mol.GetAtomWithIdx(mi).GetNeighbors():
                            mi_frag.add(nbr.GetIdx())
                        for ai in mi_frag:
                            coords[ai] = coords[ai] + shift
                        any_moved = True
            if not any_moved:
                break

    # ---- Write conformer to mol ----
    conf = Chem.Conformer(n_atoms)
    for i in range(n_atoms):
        conf.SetAtomPosition(i, Point3D(
            float(coords[i, 0]), float(coords[i, 1]), float(coords[i, 2])))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return True
