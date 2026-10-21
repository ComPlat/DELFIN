"""Multi-metal scaffold, aromatic ring snapping and scaling, bridging donors, lone-pair tilt and ligand alignment of the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.converter_flags import (
    _delfin_env_int,
)
from delfin.manta.ml_tables import (
    Chem,
    RDKIT_AVAILABLE,
    _COVALENT_RADII,
    _METAL_SET,
    _get_ml_bond_length,
)

logger = get_logger("delfin.smiles_converter")


def _build_multimetal_scaffold(
    mol,
    metal_indices: List[int],
    bridging_donors: List[Tuple[int, List[int]]],
) -> Optional[Dict[int, Tuple[float, float, float]]]:
    """Place metals + bridging donors as a rigid scaffold.

    Strategy:
    1. Metal1 at origin
    2. Each bridging donor at ideal M1-D distance along geometry vectors
    3. Metal2 placed so M2-bridge distances match ideals
    4. Returns {atom_idx: (x,y,z)} for all scaffold atoms

    This ensures the M-bridge-M core has correct geometry BEFORE
    peripheral donors are placed.
    """
    if len(metal_indices) < 2 or not bridging_donors:
        return None

    try:
        import numpy as np

        coords: Dict[int, Tuple[float, float, float]] = {}
        m1 = metal_indices[0]
        m2 = metal_indices[1]
        m1_sym = mol.GetAtomWithIdx(m1).GetSymbol()
        m2_sym = mol.GetAtomWithIdx(m2).GetSymbol()

        # Place M1 at origin.
        coords[m1] = (0.0, 0.0, 0.0)

        # Identify bridging donor atoms shared between M1 and M2.
        shared_bridges = []
        for d_idx, m_list in bridging_donors:
            if m1 in m_list and m2 in m_list:
                shared_bridges.append(d_idx)

        if not shared_bridges:
            return None

        # Place first bridging donor on +x axis at ideal M1-D distance.
        d0 = shared_bridges[0]
        d0_sym = mol.GetAtomWithIdx(d0).GetSymbol()
        d_m1 = float(_get_ml_bond_length(m1_sym, d0_sym))
        d_m2 = float(_get_ml_bond_length(m2_sym, d0_sym))
        coords[d0] = (d_m1, 0.0, 0.0)

        # Place M2: on the M1-D0-M2 plane.
        # M2 is at distance d_m2 from D0.
        # M1-D0 = d_m1 along x. M2 must be at d_m2 from D0.
        # Choose M1-D0-M2 angle ~130° (typical for μ-O bridge).
        bridge_angle_rad = math.radians(130.0)
        m2_x = d_m1 + d_m2 * math.cos(math.pi - bridge_angle_rad)
        m2_y = d_m2 * math.sin(math.pi - bridge_angle_rad)
        coords[m2] = (m2_x, m2_y, 0.0)

        # Place additional bridging donors (if any) between the two metals.
        for k, dk in enumerate(shared_bridges[1:], 1):
            dk_sym = mol.GetAtomWithIdx(dk).GetSymbol()
            dk_m1 = float(_get_ml_bond_length(m1_sym, dk_sym))
            dk_m2 = float(_get_ml_bond_length(m2_sym, dk_sym))
            # Place below the M1-M2 plane (alternating above/below)
            angle_offset = math.radians(130.0)
            sign = -1.0 if k % 2 == 1 else 1.0
            mid_x = (coords[m1][0] + coords[m2][0]) / 2.0
            mid_y = (coords[m1][1] + coords[m2][1]) / 2.0
            coords[dk] = (mid_x, mid_y, sign * dk_m1 * 0.8)

        return coords
    except Exception as exc:
        logger.debug("_build_multimetal_scaffold failed: %s", exc)
        return None


def _snap_fused_aromatic_groups_to_plane(
    mol,
    conf_id: int = 0,
    rms_threshold: float = 0.03,
) -> bool:
    """Project *fused* aromatic ring systems onto a common best-fit plane.

    A fused aromatic group = set of atoms that belong to two or more
    rings which share at least one edge, where each ring is aromatic
    (or has >= 2 unsaturated bonds).  Examples: carbazole (5+6+6),
    fluorene (5+6+6), naphthalene (6+6), indole (5+6).  Individually
    every member ring is planar but UFF can bend the ring-ring
    dihedral, leaving the collective pi-system "knickt".

    Approach per fused group:
      1. Collect all heavy atoms of all member rings.
      2. SVD best-fit plane.
      3. Minimal-movement projection (each atom shifted only by its
         component perpendicular to the plane).

    Biaryls (aryl-aryl linked by a single bond but NOT sharing an
    edge) stay unaffected: their rings are NOT fused and remain free
    to rotate.

    Returns True if any group was modified.
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

        # Qualifying rings: 5-7 members, >=2 unsat bonds, no metal, no H.
        rings = []
        for r in ri.AtomRings():
            if len(r) < 5 or len(r) > 7:
                continue
            if any(
                mol.GetAtomWithIdx(i).GetAtomicNum() <= 1
                or mol.GetAtomWithIdx(i).GetSymbol() in _METAL_SET
                for i in r
            ):
                continue
            unsat = 0
            for i in range(len(r)):
                a, b = r[i], r[(i + 1) % len(r)]
                bond = mol.GetBondBetweenAtoms(a, b)
                if bond is None:
                    continue
                if bond.GetBondType() != Chem.BondType.SINGLE or bond.GetIsAromatic():
                    unsat += 1
            if unsat < 2:
                continue
            rings.append(set(r))

        if not rings:
            return False

        # Union-find: rings that share at least one atom belong to one group.
        parent = list(range(len(rings)))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        for i in range(len(rings)):
            for j in range(i + 1, len(rings)):
                if rings[i] & rings[j]:
                    union(i, j)

        groups: Dict[int, set] = {}
        for i in range(len(rings)):
            key = find(i)
            groups.setdefault(key, set()).update(rings[i])

        changed = False
        for atoms in groups.values():
            if len(atoms) < 8:
                # Single ring - handled by per-ring snap.
                continue
            idxs = sorted(atoms)
            pts = np.array([
                [conf.GetAtomPosition(i).x,
                 conf.GetAtomPosition(i).y,
                 conf.GetAtomPosition(i).z]
                for i in idxs
            ])
            centroid = pts.mean(axis=0)
            cp = pts - centroid
            try:
                _, _, vh = np.linalg.svd(cp, full_matrices=False)
            except np.linalg.LinAlgError:
                continue
            normal = vh[-1]
            nn = np.linalg.norm(normal)
            if nn < 1e-9:
                continue
            normal /= nn
            offsets = cp @ normal
            rms = float(np.sqrt(float((offsets * offsets).mean())))
            if rms < rms_threshold:
                continue
            # Project atoms onto plane (minimal movement).
            for idx, aidx in enumerate(idxs):
                newp = pts[idx] - offsets[idx] * normal
                conf.SetAtomPosition(
                    int(aidx),
                    (float(newp[0]), float(newp[1]), float(newp[2])),
                )
            changed = True
        return changed
    except Exception as exc:
        logger.debug("_snap_fused_aromatic_groups_to_plane failed: %s", exc)
        return False


def _scale_aromatic_rings_to_ideal_cc(
    mol,
    conf_id: int = 0,
    target_cc: float = 1.39,
    lo: float = 1.34,
    hi: float = 1.44,
) -> bool:
    """Rescale squashed/stretched aromatic 5/6 rings to ideal aromatic C-C.

    The rigid/legacy assembly paths sometimes emit aromatic rings with
    in-ring C-C bond lengths well below the physical ~1.39 Å (observed
    1.20-1.30 Å on ATABOM/AFALOK), which crushes the ring-attached H toward
    each other (H-H < 0.9 Å) and trips the declash gate.  Declash cannot fix
    this — the ring bonds are rigid, zero DOF — because the *size* is wrong,
    not the rotation.

    This pass uniformly rescales each affected ring about its own centroid so
    the mean in-ring bond length becomes ``target_cc``.  Each ring atom and
    its entire pendant substituent *branch* (everything reachable through that
    ring atom, not crossing another ring atom or a metal) move rigidly by the
    same per-atom displacement, so substituent bond lengths/angles and ring-H
    are preserved exactly — only the ring breathes.

    Universal (RDKit aromaticity / ring perception, never SMILES-specific):

    * Rings of size 5 or 6 only.
    * >= 2 unsaturated (aromatic/double/triple) in-ring bonds, matching the
      planarity gate's definition (skips chair/boat saturated rings).
    * Rings containing a metal or an H ring-member are skipped (metallacycles
      and hapto faces have their own handling).
    * **M-D invariant:** any ring with a branch atom bonded to a metal is
      skipped entirely so the metal-coordination geometry is never disturbed.
    * **NEVER-WORSE:** a ring whose mean in-ring C-C already lies in
      ``[lo, hi]`` is left untouched.

    Returns True if any ring was modified.
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

        n_atoms = mol.GetNumAtoms()
        metal_idx = {
            a.GetIdx()
            for a in mol.GetAtoms()
            if a.GetSymbol() in _METAL_SET
        }

        def _pos(i):
            p = conf.GetAtomPosition(int(i))
            return np.array([p.x, p.y, p.z], dtype=float)

        # --- FUSE-AWARE SCOPING (root fix, never-worse-by-construction) -------
        # Reconstructing each aromatic ring about ITS OWN centroid is only sound
        # for an ISOLATED monocyclic aryl.  On a FUSED poly-aromatic backbone
        # (naphthalene / acenaphthene / quinoline / phenanthroline / indole /
        # carbazole ...) two rings SHARE atoms: reconstructing ring A drags the
        # fused ring B (its shared atom's rigid branch carries B along), then B
        # is rebuilt about a DIFFERENT centroid dragging A back — the shared /
        # adjacent atoms get inconsistent targets and SUPERIMPOSE, tearing the
        # backbone (measured: SIBXOS C30-C32 @ 0.07 Å; BIBMEH torn C-N / C-C).
        # So any ring that shares an atom or edge with ANOTHER aromatic ring is
        # part of a fused π-system and must NOT be reconstructed in isolation.
        # We reuse the fused-component union logic from ``_arom_planarize`` and
        # skip every fused ring, restoring the pass to its intended
        # isolated-monocyclic-aryl case (crushed pendant phenyls, e.g.
        # ATABOM / AFALOK — they share no ring atom with any other ring and are
        # untouched).  The hardened guard below is the belt-and-suspenders for
        # any fusion this skip cannot see.
        def _is_arom_candidate(_r):
            _e = len(_r)
            if _e < 5 or _e > 6:
                return False
            if any(mol.GetAtomWithIdx(i).GetAtomicNum() <= 1
                   or mol.GetAtomWithIdx(i).GetSymbol() in _METAL_SET for i in _r):
                return False
            _u = 0
            for _i in range(_e):
                _a, _b = _r[_i], _r[(_i + 1) % _e]
                _bd = mol.GetBondBetweenAtoms(_a, _b)
                if _bd is None:
                    continue
                if _bd.GetBondType() != Chem.BondType.SINGLE or _bd.GetIsAromatic():
                    _u += 1
            return _u >= 2
        _fused_atoms: set = set()
        try:
            from delfin.manta._arom_planarize import _fuse_components as _arom_fuse
            _cand_rings = [tuple(r) for r in ri.AtomRings() if _is_arom_candidate(r)]
            if len(_cand_rings) >= 2:
                for _comp in _arom_fuse(_cand_rings):
                    _cs = set(_comp)
                    if sum(1 for _r in _cand_rings if set(_r) <= _cs) >= 2:
                        _fused_atoms |= _cs
        except Exception:
            _fused_atoms = set()

        changed = False
        for ring in ri.AtomRings():
            eta = len(ring)
            if eta < 5 or eta > 6:
                continue
            ring_set = set(ring)
            if ring_set & _fused_atoms:
                continue  # fused π-system ring — skip isolated rebuild (root fix)
            if any(
                mol.GetAtomWithIdx(i).GetAtomicNum() <= 1
                or mol.GetAtomWithIdx(i).GetSymbol() in _METAL_SET
                for i in ring
            ):
                continue
            # count unsaturated in-ring bonds (same definition as the gate)
            unsat = 0
            for i in range(eta):
                a, b = ring[i], ring[(i + 1) % eta]
                bond = mol.GetBondBetweenAtoms(a, b)
                if bond is None:
                    continue
                if bond.GetBondType() != Chem.BondType.SINGLE or bond.GetIsAromatic():
                    unsat += 1
            if unsat < 2:
                continue

            # in-ring bond lengths
            lens = []
            for i in range(eta):
                a, b = ring[i], ring[(i + 1) % eta]
                if mol.GetBondBetweenAtoms(a, b) is None:
                    continue
                lens.append(float(np.linalg.norm(_pos(a) - _pos(b))))
            if not lens:
                continue
            mean_cc = float(np.mean(lens))
            min_cc = float(np.min(lens))
            max_cc = float(np.max(lens))
            if mean_cc < 1e-6:
                continue
            # NEVER-WORSE: leave physically-sized AND undistorted rings alone.
            # A ring is "good" only when its mean sits in [lo, hi] AND every
            # individual bond is within a tolerance band of the ideal (so a
            # ring whose mean is coincidentally ~1.39 but is badly distorted —
            # min 0.5 / max 2.7 — is still reconstructed).
            good = (lo <= mean_cc <= hi) and (min_cc >= lo - 0.05) \
                and (max_cc <= hi + 0.10)
            if good:
                continue

            # Build each ring atom's rigid branch (BFS into the fragment, not
            # crossing other ring atoms; STOP before stepping onto a metal so
            # we never carry the metal in a branch).  A ring atom whose branch
            # reaches a metal is "anchored": it bears the M-coordination and
            # must NOT move (M-D invariant).  Pendant aryls on P/donors are the
            # common case — exactly one anchored (ipso) atom.
            branches = {}
            anchored = set()
            for ai in ring:
                visited = set(ring_set)
                br = []
                queue = []
                for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni in metal_idx:
                        anchored.add(ai)
                        continue  # do not traverse into / move the metal
                    if ni not in visited:
                        visited.add(ni)
                        queue.append(ni)
                while queue:
                    cur = queue.pop(0)
                    br.append(cur)
                    for nbr in mol.GetAtomWithIdx(cur).GetNeighbors():
                        ni = nbr.GetIdx()
                        if ni in metal_idx:
                            anchored.add(ai)
                            continue  # stop at the metal boundary
                        if ni not in visited:
                            visited.add(ni)
                            queue.append(ni)
                branches[ai] = br
            # >=2 anchored atoms: reshaping would move an M-coordinating bond;
            # too risky, skip (chelating diaryl, fused metallacycle, …).
            if len(anchored) >= 2:
                continue

            # --- Reconstruct ring as a regular polygon at ideal C-C ---
            # Uniform scaling cannot repair a ring whose bonds range from 0.5
            # to 2.7 Å (build-distorted).  Instead lay out a regular eta-gon of
            # side ``target_cc`` in the ring's own best-fit plane, centred on
            # the current centroid, and assign each ring atom to consecutive
            # polygon vertices following the ring's bond connectivity so the
            # in-ring topology (and substituent attachment order) is preserved.
            pts = np.array([_pos(a) for a in ring], dtype=float)
            centroid = pts.mean(axis=0)
            centered = pts - centroid
            try:
                _, _, vh = np.linalg.svd(centered, full_matrices=False)
            except np.linalg.LinAlgError:
                continue
            normal = vh[-1]
            nn = float(np.linalg.norm(normal))
            if nn < 1e-9:
                continue
            normal = normal / nn
            # in-plane orthonormal basis
            ref = np.array([1.0, 0.0, 0.0]) if abs(normal[0]) < 0.9 \
                else np.array([0.0, 1.0, 0.0])
            u = np.cross(normal, ref)
            un = float(np.linalg.norm(u))
            if un < 1e-9:
                continue
            u = u / un
            v = np.cross(normal, u)

            # connectivity traversal order around the ring
            idx_map = {a: k for k, a in enumerate(ring)}
            adj = [[] for _ in range(eta)]
            for ai in ring:
                for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni in ring_set and ni != ai:
                        adj[idx_map[ai]].append(idx_map[ni])
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
                # connectivity incomplete -> fall back to angular order so we
                # still emit a sane polygon rather than skipping
                rel = centered
                ang = np.arctan2(rel @ v, rel @ u)
                order = list(np.argsort(ang))

            # circumradius for regular polygon with side = target_cc
            R = target_cc / (2.0 * np.sin(np.pi / eta))
            # anchor the polygon's starting angle to the current angular
            # position of the first traversed atom -> minimal in-plane rotation
            rel0 = centered[order[0]]
            start_angle = float(np.arctan2(np.dot(rel0, v), np.dot(rel0, u)))

            # Lay out the ideal polygon (per ring atom).
            poly = {}
            ok = True
            for k, loc in enumerate(order):
                ai = ring[loc]
                angle = start_angle + 2.0 * np.pi * k / eta
                p = centroid + R * np.cos(angle) * u + R * np.sin(angle) * v
                if not np.all(np.isfinite(p)):
                    ok = False
                    break
                poly[ai] = p
            if not ok:
                continue

            # M-D invariant: if one ring atom is anchored (its branch carries an
            # M-coordinating bond), translate the whole polygon in-plane so that
            # anchored atom lands exactly on its current position.  Then it (and
            # its metal-bearing branch) stay put while the rest of the ring +
            # their non-metal branches move.
            if anchored:
                anc = next(iter(anchored))
                shift = _pos(anc) - poly[anc]
                for ai in ring:
                    poly[ai] = poly[ai] + shift

            new_positions = {}
            for ai in ring:
                if ai in anchored:
                    continue  # frozen: keeps the M-D bond intact
                new_p = poly[ai]
                delta = new_p - _pos(ai)
                new_positions[ai] = new_p
                for bi in branches.get(ai, ()):  # rigid branch translation
                    bp = _pos(bi) + delta
                    if not np.all(np.isfinite(bp)):
                        ok = False
                        break
                    new_positions[bi] = bp
                if not ok:
                    break
            if not ok or not new_positions:
                continue

            # --- per-ring NEVER-WORSE guard ---
            # Reconstructing one ring must not shove its atoms into the rest of
            # the structure: compute the closest contact between any MOVED atom
            # and any non-moved atom before and after.  Keep the move only if it
            # does not introduce a new sub-0.9 Å clash that is worse than what
            # was there before.  Pure ring-size repair always passes; a move
            # that would create a fresh inter-fragment crash is reverted.
            moved = list(new_positions.keys())
            moved_set = set(moved)
            # --- HARD crush floors (fused-ring safety, absolute) --------------
            # The moved-vs-NON-moved check below is BLIND to a superposition
            # WITHIN the moved set itself (the fused-ring tug-of-war signature:
            # ring A places the shared atoms, ring B's later rebuild drags them
            # to an inconsistent target -> two moved atoms land on top of each
            # other).  Add (1) a moved-vs-MOVED heavy crush floor and (2) a hard
            # bond-length-collapse floor over every bond touching a moved atom.
            # A legitimate isolated-ring repair yields a clean 1.39 Å polygon
            # (moved-moved min ~1.0 Å, no collapsed bond) so it always passes;
            # only a crush/collapse — which is never physical — is reverted.
            _FLOOR = 0.90
            _mh = [j for j in moved
                   if mol.GetAtomWithIdx(j).GetAtomicNum() > 1]
            _crush = False
            if len(_mh) >= 2:
                _MP = np.array([new_positions[j] for j in _mh])
                for _k in range(len(_mh) - 1):
                    _dd = np.linalg.norm(_MP[_k + 1:] - _MP[_k], axis=1)
                    if _dd.size and float(_dd.min()) < _FLOOR:
                        _crush = True
                        break
            if not _crush:
                for _b in mol.GetBonds():
                    _a1, _a2 = _b.GetBeginAtomIdx(), _b.GetEndAtomIdx()
                    if _a1 not in moved_set and _a2 not in moved_set:
                        continue
                    if (mol.GetAtomWithIdx(_a1).GetAtomicNum() <= 1
                            or mol.GetAtomWithIdx(_a2).GetAtomicNum() <= 1):
                        continue
                    _p1 = new_positions.get(_a1, _pos(_a1))
                    _p2 = new_positions.get(_a2, _pos(_a2))
                    if float(np.linalg.norm(_p1 - _p2)) < _FLOOR:
                        _crush = True
                        break
            if _crush:
                continue  # revert: reconstruction crushed atoms / collapsed a bond
            others = [j for j in range(n_atoms) if j not in moved_set]
            if others:
                old_moved = np.array([_pos(j) for j in moved])
                new_moved = np.array([new_positions[j] for j in moved])
                oth = np.array([_pos(j) for j in others])

                def _min_contact(mpos):
                    mn = np.inf
                    for k in range(mpos.shape[0]):
                        dd = np.linalg.norm(oth - mpos[k], axis=1)
                        if dd.size:
                            mn = min(mn, float(dd.min()))
                    return mn
                before = _min_contact(old_moved)
                after = _min_contact(new_moved)
                CLASH = 0.90
                if after < before and after < CLASH:
                    continue  # revert: would create a worse clash
            for idx, p in new_positions.items():
                conf.SetAtomPosition(
                    int(idx), (float(p[0]), float(p[1]), float(p[2]))
                )
            changed = True
        return changed
    except Exception as exc:
        logger.debug("_scale_aromatic_rings_to_ideal_cc failed: %s", exc)
        return False


def _snap_aromatic_rings_to_plane(
    mol,
    conf_id: int = 0,
    rms_threshold: float = 0.05,
) -> bool:
    """Project aromatic / unsaturated ring atoms onto their best-fit plane.

    UFF relaxation can buckle 5-7 membered aromatic rings just enough that
    the downstream planarity gate (``_has_pi_ring_nonplanarity``) rejects
    otherwise-valid isomers.  Rather than relaxing the gate (which would
    admit genuinely bad geometries), we project each such ring onto its
    SVD best-fit plane.  The projection is *minimal-movement* (each atom
    moves only the component perpendicular to the plane), so bond
    lengths, bond angles in the plane, and the surrounding heavy-atom
    positions are preserved to first order.

    Only rings with >= 2 unsaturated (aromatic / double / triple) bonds
    are touched, to avoid flattening chair / boat conformations of
    saturated rings.  Rings containing a metal are skipped (metallacycles
    have their own planarity handling).

    Returns True if any ring was modified.
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

        changed = False
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
            for i in range(len(ring)):
                a, b = ring[i], ring[(i + 1) % len(ring)]
                bond = mol.GetBondBetweenAtoms(a, b)
                if bond is None:
                    continue
                if bond.GetBondType() != Chem.BondType.SINGLE or bond.GetIsAromatic():
                    unsat += 1
            if unsat < 2:
                continue

            pts = np.array([
                [conf.GetAtomPosition(i).x,
                 conf.GetAtomPosition(i).y,
                 conf.GetAtomPosition(i).z]
                for i in ring
            ])
            centroid = pts.mean(axis=0)
            centered = pts - centroid
            try:
                _, _, vh = np.linalg.svd(centered, full_matrices=False)
            except np.linalg.LinAlgError:
                continue
            normal = vh[-1]
            nn = np.linalg.norm(normal)
            if nn < 1e-9:
                continue
            normal /= nn
            # Out-of-plane deviations
            offsets = centered @ normal
            rms = float(np.sqrt(float((offsets * offsets).mean())))
            #
            # Aromatic-H plane invariant: ring-attached H must follow the
            # π-system as a rigid body.  RDKit ETKDG on kekulized SMILES
            # (especially with cations like [N+]) produces a planar
            # heavy-atom ring but places ring-H 0.5-1.4 Å out of plane.
            # Project each ring-H onto the ring plane via outward-radial
            # placement at the element-specific ideal X-H bond length.
            #
            # H-projection runs INDEPENDENTLY of heavy-atom snap: even when
            # ring is already planar (rms < threshold) the H atoms can still
            # be out of plane.  Ring-atom snap only fires when rms >= threshold.
            _ideal_xh = {"C": 1.08, "N": 1.01, "O": 0.96, "S": 1.34}
            ring_did_snap = False
            if rms >= rms_threshold:
                # Minimal-movement projection: subtract normal-component
                # from each ring atom.  Substituent heavy atoms unchanged.
                for idx, atom_idx in enumerate(ring):
                    new_p = pts[idx] - offsets[idx] * normal
                    conf.SetAtomPosition(
                        int(atom_idx),
                        (float(new_p[0]), float(new_p[1]), float(new_p[2])),
                    )
                ring_did_snap = True
                changed = True

            # Step 2: project ring-attached H atoms onto current ring plane
            # (uses CURRENT parent position regardless of whether ring was
            # snapped this pass — handles the cf1d480 case where ring was
            # already planar but H was placed off-plane by ETKDG).
            for atom_idx in ring:
                try:
                    parent_atom = mol.GetAtomWithIdx(int(atom_idx))
                    h_nbrs = [
                        n for n in parent_atom.GetNeighbors()
                        if n.GetAtomicNum() == 1
                    ]
                    if len(h_nbrs) != 1:
                        continue  # CH2/NH2 etc.: not a single-H aromatic
                    h_idx = h_nbrs[0].GetIdx()
                    parent_pos_obj = conf.GetAtomPosition(int(atom_idx))
                    parent_pos = np.array([
                        parent_pos_obj.x, parent_pos_obj.y, parent_pos_obj.z
                    ])
                    h_pos = conf.GetAtomPosition(h_idx)
                    h_arr = np.array([h_pos.x, h_pos.y, h_pos.z])
                    h_oop = float(abs(np.dot(h_arr - centroid, normal)))
                    if h_oop <= 0.10:
                        continue  # already in-plane, leave alone
                    outward = parent_pos - centroid  # in-plane vector
                    n_out = float(np.linalg.norm(outward))
                    if n_out < 1e-6:
                        continue
                    outward_unit = outward / n_out
                    sym = parent_atom.GetSymbol()
                    d_xh = _ideal_xh.get(sym, 1.08)
                    h_new = parent_pos + d_xh * outward_unit
                    conf.SetAtomPosition(
                        int(h_idx),
                        (float(h_new[0]), float(h_new[1]), float(h_new[2])),
                    )
                    changed = True
                except Exception:
                    continue
        # Phase 5A (2026-05-12): universal H-bond-length normalization for
        # non-ring/sp3/methyl/orphan H atoms.  Per Wave-2A forensics, η-isomer
        # emit paths (HD-TA, OB-WRS) bypass the normal UFF+snap pipeline and
        # produce collapsed C-H bonds (D-JESNAA01 H min 0.54 Å observed).
        # Aromatic-H projection above only fixes ring-attached H; the
        # additional pass below normalizes any other H whose parent-distance
        # is outside [0.92, 1.20] Å to ideal X-H along the existing
        # parent→H direction.  Env-gated default OFF for safety.
        if _delfin_env_int("DELFIN_H_UNIVERSAL_NORMALIZE", 0):
            try:
                for atom in mol.GetAtoms():
                    if atom.GetAtomicNum() != 1:
                        continue
                    nbrs = atom.GetNeighbors()
                    if len(nbrs) != 1:
                        continue
                    parent = nbrs[0]
                    if parent.GetIsAromatic() and parent.IsInRing():
                        continue
                    h_idx = atom.GetIdx()
                    p_idx = parent.GetIdx()
                    p_pos_o = conf.GetAtomPosition(p_idx)
                    h_pos_o = conf.GetAtomPosition(h_idx)
                    p_pos = np.array([p_pos_o.x, p_pos_o.y, p_pos_o.z])
                    h_pos = np.array([h_pos_o.x, h_pos_o.y, h_pos_o.z])
                    vec = h_pos - p_pos
                    d = float(np.linalg.norm(vec))
                    if d < 0.1:
                        continue
                    if 0.92 <= d <= 1.20:
                        continue
                    sym = parent.GetSymbol()
                    d_xh = _ideal_xh.get(sym, 1.08)
                    h_new = p_pos + vec / d * d_xh
                    conf.SetAtomPosition(
                        int(h_idx),
                        (float(h_new[0]), float(h_new[1]), float(h_new[2])),
                    )
                    changed = True
            except Exception as exc:
                logger.debug("H universal-normalize skipped: %s", exc)
        return changed
    except Exception as exc:
        logger.debug("_snap_aromatic_rings_to_plane failed: %s", exc)
        return False


def _snap_aromatic_rings_in_xyz(
    xyz_str: str,
    mol_template,
    rms_threshold: float = 0.05,
) -> str:
    """Parse XYZ, call :func:`_snap_aromatic_rings_to_plane`, write XYZ back.

    Used post-UFF so downstream planarity gates see rings that carry at
    most ``rms_threshold`` of residual buckle per atom.  Hydrogens are
    added to the template when needed so the conformer atom count
    matches the XYZ.
    """
    if not RDKIT_AVAILABLE or mol_template is None:
        return xyz_str
    try:
        lines = [l for l in xyz_str.strip().splitlines() if l.strip()]
        if not lines:
            return xyz_str
        mol = Chem.RWMol(mol_template)
        mol.RemoveAllConformers()
        n = mol.GetNumAtoms()
        if len(lines) == n:
            mol_use = mol.GetMol()
        else:
            mol_h = Chem.AddHs(mol)
            if len(lines) == mol_h.GetNumAtoms():
                mol_use = mol_h
            else:
                return xyz_str
        conf = Chem.Conformer(mol_use.GetNumAtoms())
        for i, l in enumerate(lines):
            parts = l.split()
            if len(parts) < 4:
                return xyz_str
            conf.SetAtomPosition(
                i, (float(parts[1]), float(parts[2]), float(parts[3]))
            )
        cid = mol_use.AddConformer(conf, assignId=True)
        # Aryl-ring-size pass FIRST (env-gated, default OFF, byte-identical
        # when unset): the rigid/legacy build paths can emit aromatic rings
        # with crushed in-ring C-C (~1.20-1.30 Å) that squeeze ring-H into a
        # < 0.9 Å clash the declash gate cannot resolve (ring bonds are
        # rigid).  Rescaling each ring to the ideal aromatic C-C before the
        # plane snap fixes the root cause (ring SIZE, not rotation).
        changed_size = False
        if _delfin_env_int("DELFIN_FFFREE_ARYL_RING_SIZE", 0):
            changed_size = _scale_aromatic_rings_to_ideal_cc(mol_use, cid)
        # Fused-group snap (collective plane for carbazole-style fused
        # aromatics), then per-ring snap as a second pass so any residual
        # single-ring buckle is also removed.
        changed_fused = _snap_fused_aromatic_groups_to_plane(
            mol_use, cid, rms_threshold=rms_threshold,
        )
        changed_rings = _snap_aromatic_rings_to_plane(
            mol_use, cid, rms_threshold=rms_threshold
        )
        changed = bool(changed_size or changed_fused or changed_rings)
        if not changed:
            return xyz_str
        new_conf = mol_use.GetConformer(cid)
        out = []
        for i in range(mol_use.GetNumAtoms()):
            atom = mol_use.GetAtomWithIdx(i)
            p = new_conf.GetAtomPosition(i)
            out.append(
                f"{atom.GetSymbol():4s} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}"
            )
        return '\n'.join(out) + '\n'
    except Exception as exc:
        logger.debug("_snap_aromatic_rings_in_xyz failed: %s", exc)
        return xyz_str


def _scale_aromatic_rings_in_xyz_from_smiles(xyz_str: str, smiles: str) -> str:
    """Universal aryl-ring-size finalizer for a finished XYZ.

    Delegates to :func:`_scale_aromatic_rings_in_xyz_geom`, which perceives
    rings purely from the geometry.  The metal / hapto build path mutates the
    molecule (dative-bond rewrites, charge/H adjustments) so a SMILES-rebuilt
    mol no longer maps 1:1 onto the emitted geometry — geometry-only perception
    is the only path-independent way to reach those squashed rings.  The
    ``smiles`` argument is retained for call-site compatibility but unused.

    Env-gated default OFF (``DELFIN_FFFREE_ARYL_RING_SIZE``).
    """
    return _scale_aromatic_rings_in_xyz_geom(xyz_str)


def _scale_aromatic_rings_in_xyz_geom(xyz_str: str) -> str:
    """Geometry-only aromatic 5/6-ring size normalization on a finished XYZ.

    This is the universal, path-independent finalizer.  The metal / hapto build
    path mutates the molecule (dative-bond rewrites, charge/H adjustments) so a
    SMILES-rebuilt mol no longer maps 1:1 onto the emitted geometry — the only
    thing we can trust is the XYZ itself.  This pass therefore perceives bonds
    purely from interatomic distances (covalent radii + tolerance), finds
    near-planar 5/6 carbon rings, and reconstructs any whose in-ring C-C is
    crushed/stretched to a regular polygon at the ideal aromatic C-C (1.39 Å),
    rigidly dragging each ring atom's substituent branch.

    Universal (never SMILES-specific).  M-D safe: a ring atom whose branch
    reaches a metal is "anchored" (frozen) so the M-coordination bond is never
    moved; rings with >= 2 anchors are skipped.  NEVER-WORSE: physically-sized,
    undistorted rings are left untouched.  Env-gated (default OFF); returns the
    input unchanged on any failure or non-finite result.
    """
    if not xyz_str:
        return xyz_str
    if not _delfin_env_int("DELFIN_FFFREE_ARYL_RING_SIZE", 0):
        return xyz_str
    try:
        import numpy as np
        target_cc = 1.39
        lo, hi = 1.34, 1.44

        raw = [l for l in xyz_str.splitlines() if l.strip()]
        # tolerate an optional XYZ header (count + comment)
        coord_lines = []
        for l in raw:
            p = l.split()
            if len(p) >= 4:
                try:
                    float(p[1]); float(p[2]); float(p[3])
                except ValueError:
                    continue
                coord_lines.append(l)
        n = len(coord_lines)
        if n < 5:
            return xyz_str
        syms = []
        coords = np.zeros((n, 3), dtype=float)
        for i, l in enumerate(coord_lines):
            p = l.split()
            syms.append(p[0])
            coords[i] = (float(p[1]), float(p[2]), float(p[3]))

        metal_set = _METAL_SET
        metal_idx = {i for i, s in enumerate(syms) if s in metal_set}

        # --- bond graph from covalent radii (heavy + H), exclude metals ---
        def _rad(s):
            return _COVALENT_RADII.get(s, 0.76)

        adj = [[] for _ in range(n)]
        heavy = [i for i in range(n) if syms[i] != 'H' and i not in metal_idx]
        # pairwise within a generous cutoff
        for ii in range(len(heavy)):
            i = heavy[ii]
            for jj in range(ii + 1, len(heavy)):
                j = heavy[jj]
                d = float(np.linalg.norm(coords[i] - coords[j]))
                cut = _rad(syms[i]) + _rad(syms[j]) + 0.40
                if 0.4 < d <= cut:
                    adj[i].append(j)
                    adj[j].append(i)
        # H bonds (for branch dragging only; H not part of ring search)
        for i in range(n):
            if syms[i] != 'H':
                continue
            best, bestd = -1, 1e9
            for j in heavy:
                d = float(np.linalg.norm(coords[i] - coords[j]))
                cut = _rad('H') + _rad(syms[j]) + 0.40
                if d <= cut and d < bestd:
                    best, bestd = j, d
            if best >= 0:
                adj[i].append(best)
                adj[best].append(i)
        # metal-coordination edges (so branch BFS can detect "reaches metal")
        for m in metal_idx:
            for j in heavy:
                d = float(np.linalg.norm(coords[m] - coords[j]))
                cut = _COVALENT_RADII.get(syms[m], 1.35) + _rad(syms[j]) + 0.55
                if d <= cut:
                    adj[m].append(j)
                    adj[j].append(m)

        # --- find carbon 5/6 rings via bounded DFS over the carbon subgraph ---
        carbons = [i for i in heavy if syms[i] == 'C']
        cset = set(carbons)
        cadj = {i: [j for j in adj[i] if j in cset] for i in carbons}
        rings = set()

        def _dfs(start, cur, path, depth):
            if depth > 6:
                return
            for nb in cadj[cur]:
                if nb == start and len(path) in (5, 6):
                    rings.add(frozenset(path))
                elif nb not in path and len(path) < 6:
                    _dfs(start, nb, path + [nb], depth + 1)

        for s in carbons:
            _dfs(s, s, [s], 1)

        # keep rings where every member has exactly 2 in-ring C neighbours
        clean_rings = []
        for r in rings:
            rl = list(r)
            if all(len([x for x in cadj[m] if x in r]) == 2 for m in rl):
                clean_rings.append(rl)

        # FUSE-AWARE SCOPING (root fix): reconstructing a ring of a FUSED
        # π-system about its own centroid tears the backbone (see the RDKit
        # twin _scale_aromatic_rings_to_ideal_cc — SIBXOS acenaphthene is an
        # all-carbon fused system this geometry pass also perceives).  Skip any
        # ring that shares an atom with another perceived ring; only isolated
        # monocyclic aryls are reconstructed.  Reuses the _arom_planarize union.
        _fused_atoms = set()
        try:
            from delfin.manta._arom_planarize import _fuse_components as _arom_fuse
            _cr_t = [tuple(r) for r in clean_rings]
            if len(_cr_t) >= 2:
                for _comp in _arom_fuse(_cr_t):
                    _cs = set(_comp)
                    if sum(1 for _r in _cr_t if set(_r) <= _cs) >= 2:
                        _fused_atoms |= _cs
        except Exception:
            _fused_atoms = set()

        def _pos(i):
            return coords[i]

        changed = False
        for ring in clean_rings:
            eta = len(ring)
            ring_set = set(ring)
            if ring_set & _fused_atoms:
                continue  # fused π-system ring — skip isolated rebuild (root fix)
            # ordered ring traversal via connectivity
            radj = {m: [x for x in cadj[m] if x in ring_set] for m in ring}
            order = []
            visited = set()
            cur, prev = ring[0], -1
            for _ in range(eta):
                order.append(cur)
                visited.add(cur)
                nxt = -1
                for nb in radj[cur]:
                    if nb != prev and nb not in visited:
                        nxt = nb
                        break
                if nxt == -1:
                    break
                prev, cur = cur, nxt
            if len(order) != eta:
                continue

            lens = [float(np.linalg.norm(_pos(order[k]) - _pos(order[(k + 1) % eta])))
                    for k in range(eta)]
            mean_cc = float(np.mean(lens))
            min_cc = float(np.min(lens))
            max_cc = float(np.max(lens))
            if mean_cc < 1e-6:
                continue

            # Aromatic (sp2) vs saturated (sp3) discriminator that survives a
            # collapsed geometry (where loose distance cutoffs would create
            # spurious "bonds").  A benzene-type ring carbon carries at most ONE
            # substituent — either a single H (aromatic C-H) or a single heavy
            # group — while a cyclohexane CH2 carries TWO hydrogens.  Counting
            # substituents with a STRICT covalent cutoff (immune to the crush
            # contacts) and requiring <= 1 H and <= 1 heavy substituent per ring
            # member excludes chair/boat saturated rings (which must never be
            # flattened) without a planarity test that a distorted aromatic ring
            # would itself fail.  The reconstruction re-planarizes the ring.
            def _strict_subs(m):
                nh = 0
                nheavy = 0
                pm = _pos(m)
                for nb in adj[m]:
                    if nb in ring_set or nb in metal_idx:
                        continue
                    d = float(np.linalg.norm(pm - _pos(nb)))
                    cut = _rad(syms[m]) + _rad(syms[nb]) + 0.18
                    if d > cut:
                        continue  # spurious crush contact, not a real bond
                    if syms[nb] == 'H':
                        nh += 1
                    else:
                        nheavy += 1
                return nh, nheavy
            saturated = False
            for m in ring:
                nh, nheavy = _strict_subs(m)
                if nh >= 2 or nheavy >= 2 or (nh + nheavy) >= 2:
                    saturated = True
                    break
            if saturated:
                continue

            rp = np.array([_pos(m) for m in ring])
            rc = rp.mean(axis=0)
            rcent = rp - rc
            try:
                _, _, rvh = np.linalg.svd(rcent, full_matrices=False)
            except np.linalg.LinAlgError:
                continue
            rnormal = rvh[-1]
            rnn = float(np.linalg.norm(rnormal))
            if rnn < 1e-9:
                continue
            rnormal = rnormal / rnn

            # NEVER-WORSE: physically-sized AND undistorted -> leave alone
            good = (lo <= mean_cc <= hi) and (min_cc >= lo - 0.05) \
                and (max_cc <= hi + 0.10)
            if good:
                continue

            # branches (BFS into fragment, stop before metals); anchored atoms
            branches = {}
            anchored = set()
            for ai in ring:
                seen = set(ring_set)
                br = []
                q = []
                for nb in adj[ai]:
                    if nb in metal_idx:
                        anchored.add(ai)
                        continue
                    if nb not in seen:
                        seen.add(nb)
                        q.append(nb)
                while q:
                    c = q.pop(0)
                    br.append(c)
                    for nb in adj[c]:
                        if nb in metal_idx:
                            anchored.add(ai)
                            continue
                        if nb not in seen:
                            seen.add(nb)
                            q.append(nb)
                branches[ai] = br
            if len(anchored) >= 2:
                continue

            # in-plane basis + ideal polygon at ideal C-C
            ref = np.array([1.0, 0.0, 0.0]) if abs(rnormal[0]) < 0.9 \
                else np.array([0.0, 1.0, 0.0])
            u = np.cross(rnormal, ref)
            un = float(np.linalg.norm(u))
            if un < 1e-9:
                continue
            u = u / un
            v = np.cross(rnormal, u)
            R = target_cc / (2.0 * np.sin(np.pi / eta))
            rel0 = rcent[ring.index(order[0])]
            start_angle = float(np.arctan2(np.dot(rel0, v), np.dot(rel0, u)))

            poly = {}
            ok = True
            for k, ai in enumerate(order):
                angle = start_angle + 2.0 * np.pi * k / eta
                p = rc + R * np.cos(angle) * u + R * np.sin(angle) * v
                if not np.all(np.isfinite(p)):
                    ok = False
                    break
                poly[ai] = p
            if not ok:
                continue
            if anchored:
                anc = next(iter(anchored))
                shift = _pos(anc) - poly[anc]
                for ai in ring:
                    poly[ai] = poly[ai] + shift

            updates = {}
            for ai in ring:
                if ai in anchored:
                    continue
                new_p = poly[ai]
                delta = new_p - _pos(ai)
                updates[ai] = new_p
                for bi in branches.get(ai, ()):
                    bp = _pos(bi) + delta
                    if not np.all(np.isfinite(bp)):
                        ok = False
                        break
                    updates[bi] = bp
                if not ok:
                    break
            if not ok or not updates:
                continue
            # --- HARD NEVER-WORSE guard (this pass had NONE) ------------------
            # Never commit a reconstruction that crushes two heavy atoms of the
            # moved set together (fused tug-of-war), drives a moved atom into a
            # fixed atom, or collapses a bond.  A clean isolated-ring repair
            # (1.39 Å polygon, rigid substituent drag) never trips these.
            _FLOOR = 0.90
            _mset = set(updates.keys())
            _mh = [j for j in _mset if syms[j] != 'H']
            _crush = False
            if len(_mh) >= 2:
                _MP = np.array([updates[j] for j in _mh])
                for _k in range(len(_mh) - 1):
                    _dd = np.linalg.norm(_MP[_k + 1:] - _MP[_k], axis=1)
                    if _dd.size and float(_dd.min()) < _FLOOR:
                        _crush = True
                        break
            if not _crush:
                _oh = [j for j in heavy if j not in _mset]
                if _oh and _mh:
                    _OP = np.array([_pos(j) for j in _oh])
                    for j in _mh:
                        _dd = np.linalg.norm(_OP - updates[j], axis=1)
                        if _dd.size and float(_dd.min()) < _FLOOR:
                            _crush = True
                            break
            if not _crush:
                for j in _mset:
                    if syms[j] == 'H':
                        continue
                    _pj = updates[j]
                    for _nb in adj[j]:
                        if syms[_nb] == 'H' or _nb in metal_idx:
                            continue
                        _pn = updates.get(_nb, _pos(_nb))
                        if float(np.linalg.norm(_pj - _pn)) < _FLOOR:
                            _crush = True
                            break
                    if _crush:
                        break
            if _crush:
                continue  # revert this ring: reconstruction would crush / collapse
            for idx, p in updates.items():
                coords[idx] = p
            changed = True

        if not changed:
            return xyz_str
        out = []
        for i in range(n):
            out.append(f"{syms[i]:4s} {coords[i,0]:12.6f} "
                       f"{coords[i,1]:12.6f} {coords[i,2]:12.6f}")
        return '\n'.join(out) + '\n'
    except Exception as exc:
        logger.debug("_scale_aromatic_rings_in_xyz_geom failed: %s", exc)
        return xyz_str


def _snap_bridging_donors_to_compromise(
    xyz_str: str,
    mol_template,
    metal_indices: List[int],
    bridging: List[Tuple[int, List[int]]],
) -> str:
    """Move every bridging donor onto the metal-metal line at a distance
    satisfying both metal--donor bond-length ideals simultaneously.

    The per-metal builds in the multinuclear coupled-enumeration loop
    each overwrite the bridging donor with that metal's polyhedron
    vertex, so after the two sequential builds the bridge is pinned at
    metal 2's ideal but arbitrarily far from metal 1 (failing Rule 1).
    This helper analytically solves for the bridge position on the
    ``metal_1 - metal_2`` line that sits at exactly ``d_M1^ideal`` from
    metal 1 (which automatically yields the correct ``d_M2^ideal``
    when the two metals are positioned by ``_build_multimetal_scaffold``
    on the M1-D-M2 line).  If the metals have drifted off the line
    (e.g. after template ETKDG plus two rescales), we snap the bridge
    to the fractional position ``t = d_M1^ideal / (d_M1^ideal + d_M2^ideal)``
    along the metal-metal segment; that keeps both distances within
    ``~[0.9, 1.1] * ideal`` even when the M-M separation is off by a
    few tenths of an Å, which is well inside the graph-gate bridging
    window [0.55, 3.00].
    """
    try:
        lines = [l for l in xyz_str.splitlines() if l.strip()]
        if len(lines) != mol_template.GetNumAtoms():
            return xyz_str
        coords: List[List[float]] = []
        for ln in lines:
            p = ln.split()
            coords.append([float(p[1]), float(p[2]), float(p[3])])
        if len(metal_indices) < 2:
            return xyz_str
        m1, m2 = int(metal_indices[0]), int(metal_indices[1])
        m1_sym = mol_template.GetAtomWithIdx(m1).GetSymbol()
        m2_sym = mol_template.GetAtomWithIdx(m2).GetSymbol()
        p1 = coords[m1]
        p2 = coords[m2]
        vec = [p2[k] - p1[k] for k in range(3)]
        vec_norm = math.sqrt(sum(v * v for v in vec))
        if vec_norm < 1e-6:
            return xyz_str
        ux = [v / vec_norm for v in vec]
        changed = False
        for d_idx, m_list in bridging:
            if m1 not in m_list or m2 not in m_list:
                continue
            d_sym = mol_template.GetAtomWithIdx(d_idx).GetSymbol()
            d_m1 = float(_get_ml_bond_length(m1_sym, d_sym))
            d_m2 = float(_get_ml_bond_length(m2_sym, d_sym))
            if d_m1 <= 0 or d_m2 <= 0:
                continue
            t = d_m1 / max(d_m1 + d_m2, 1e-6)
            target = [p1[k] + t * vec_norm * ux[k] for k in range(3)]
            coords[d_idx] = target
            changed = True
        if not changed:
            return xyz_str
        out = []
        for i, ln in enumerate(lines):
            sym = ln.split()[0]
            x, y, z = coords[i]
            out.append(f"{sym:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
        return "\n".join(out) + "\n"
    except Exception as exc:
        logger.debug("_snap_bridging_donors_to_compromise failed: %s", exc)
        return xyz_str


def _compute_lp_tilt_rotations(
    donors: List[int],
    pts,
    m_pos,
    mol,
    levels=(0.5, 0.9),
):
    """Return a list of (R_tilt, pivot) for partial LP-alignment corrections.

    Multiple correction levels let the optimiser find a compromise between
    LP-alignment (which demands rotation toward M) and inter-fragment
    clash-avoidance (which may prefer a slightly misaligned LP over a
    colliding backbone).  Each level k in (0,1] rotates by
    ``k * (theta_full - residual_10deg)``; 0 is equivalent to no tilt
    (omitted from the returned list, handled by the caller).

    LP-direction uses the anti-bisector of D -> non-metal-neighbour
    vectors: same universal rule for sp/sp2/sp3.  Only applicable to
    monodentate non-bridging donors that have ring/atom neighbours.
    """
    try:
        import numpy as _np
    except Exception:
        return []
    if len(donors) != 1:
        return []
    d_idx = donors[0]
    d_atom = mol.GetAtomWithIdx(d_idx)
    if sum(1 for nb in d_atom.GetNeighbors() if nb.GetSymbol() in _METAL_SET) >= 2:
        return []
    non_metal_nbrs = [
        nb.GetIdx() for nb in d_atom.GetNeighbors()
        if nb.GetSymbol() not in _METAL_SET
    ]
    if not non_metal_nbrs:
        return []
    d_pos = pts[d_idx]
    lp = _np.zeros(3)
    for nb_idx in non_metal_nbrs:
        v = pts[nb_idx] - d_pos
        vn = _np.linalg.norm(v)
        if vn > 1e-8:
            lp += v / vn
    lp_norm = _np.linalg.norm(lp)
    if lp_norm < 1e-6:
        return []
    lp_unit = -lp / lp_norm
    dm = m_pos - d_pos
    dm_norm = _np.linalg.norm(dm)
    if dm_norm < 1e-6:
        return []
    dm_unit = dm / dm_norm
    cos_t = max(-1.0, min(1.0, float(_np.dot(lp_unit, dm_unit))))
    theta_full = float(math.acos(cos_t))
    if theta_full <= math.radians(15.0):
        return []
    axis = _np.cross(lp_unit, dm_unit)
    axis_norm = _np.linalg.norm(axis)
    if axis_norm < 1e-6:
        tmp = _np.array([1.0, 0.0, 0.0])
        if abs(float(_np.dot(lp_unit, tmp))) > 0.9:
            tmp = _np.array([0.0, 1.0, 0.0])
        axis = _np.cross(lp_unit, tmp)
        axis_norm = _np.linalg.norm(axis)
        if axis_norm < 1e-6:
            return []
    axis /= axis_norm
    theta_target = max(0.0, theta_full - math.radians(10.0))
    out = []
    for k in levels:
        theta = theta_target * float(k)
        if theta < 1e-4:
            continue
        c = math.cos(theta); s = math.sin(theta)
        K = _np.array([
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ])
        R_tilt = _np.eye(3) + s * K + (1.0 - c) * K @ K
        out.append((R_tilt, d_pos))
    return out


def _align_and_orient_ligands(
    coords,
    mol,
    metal_idx: int,
    donor_atom_indices: List[int],
    fragments: Optional[List[set]] = None,
    n_rot_mono: int = 12,
    n_rot_bi: int = 12,
    passes: int = 5,
    lp_weight: float = 0.0,
    sym_weight: float = 0.0,
):
    """Rotate each ligand fragment so donor stays at its polyhedron vertex
    and the rest of the fragment is oriented to jointly minimise:
      (a) inter-fragment clash (existing covalent-radius term + M-intrusion),
      (b) lone-pair / M-D misalignment for monodentate sp/sp2/sp3 donors.

    Trials per fragment = {no-tilt, tilt@50%, tilt@90%} x {N spin angles}.
    Passes iterate in alternating fragment order (approximate simultaneous
    optimisation -- breaks the greedy "first fragment wins" bias).

    Modifies ``coords`` in place.  Accepts ``list`` of tuples or numpy array.
    """
    if not RDKIT_AVAILABLE or mol is None:
        return coords
    try:
        import numpy as _np
    except Exception:
        return coords
    try:
        n_atoms = mol.GetNumAtoms()
        if len(coords) != n_atoms:
            return coords

        # Normalise to a numpy working copy.
        is_list = not isinstance(coords, _np.ndarray)
        pts = _np.array([list(coords[i]) for i in range(n_atoms)], dtype=float)

        # Fragment decomposition (non-metal connected components).
        if fragments is None:
            non_metal = {
                a.GetIdx() for a in mol.GetAtoms()
                if a.GetSymbol() not in _METAL_SET and a.GetAtomicNum() > 1
            }
            adj: Dict[int, set] = {i: set() for i in non_metal}
            for bond in mol.GetBonds():
                bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if bi in non_metal and bj in non_metal:
                    adj[bi].add(bj)
                    adj[bj].add(bi)
            visited: set = set()
            fragments = []
            for start in sorted(non_metal):
                if start in visited:
                    continue
                stack = [start]
                frag: set = set()
                while stack:
                    node = stack.pop()
                    if node in visited:
                        continue
                    visited.add(node)
                    frag.add(node)
                    for nb in adj.get(node, ()):
                        if nb not in visited:
                            stack.append(nb)
                fragments.append(frag)

        donor_set = set(donor_atom_indices)
        # Per-fragment heavy atom list + donor list.
        frag_info: List[dict] = []
        for frag in fragments:
            frag_heavy = sorted(
                i for i in frag
                if mol.GetAtomWithIdx(i).GetAtomicNum() > 1
            )
            frag_donors = [i for i in frag_heavy if i in donor_set]
            if len(frag_heavy) <= len(frag_donors):
                continue  # nothing to rotate
            frag_info.append({
                "heavy": frag_heavy,
                "donors": frag_donors,
            })

        m_pos = pts[metal_idx]
        r_cov = {
            i: _COVALENT_RADII.get(
                mol.GetAtomWithIdx(i).GetSymbol(), 0.75
            )
            for i in range(n_atoms)
        }

        def _clash_for_frag(frag_heavy, frag_pts):
            """Penalty between this fragment's non-donor atoms and every
            heavy atom of all OTHER fragments, plus the metal itself.

            PERF (BYTE-IDENTICAL): the per-(point x other_heavy) distance was
            the profiled norm-storm (~8.1M scalar np.linalg.norm calls on
            ligand-rich complexes).  It is replaced by ONE matmul-ddot row
            norm per fragment point: np.matmul(d[:,None,:], d[:,:,None]) ->
            sqrt dispatches the SAME BLAS ddot kernel as scalar
            np.linalg.norm(v) = sqrt(v.dot(v)), so each distance is
            bit-identical (FMA/summation order preserved).  The score
            accumulation loop is kept SCALAR in the exact original order
            (frag_pts.items() order, metal-term first, then j in other_heavy
            order) so the float += summation order -- and hence the emitted
            geometry -- is unchanged.  Invariant arrays (other-heavy index
            list, their positions, their radii) are hoisted out of the
            per-point loop (recomputed identically before).
            """
            score = 0.0
            other_heavy = [
                j for info in frag_info
                if set(info["heavy"]) != set(frag_heavy)
                for j in info["heavy"]
            ]
            m_sym = mol.GetAtomWithIdx(metal_idx).GetSymbol()
            # Hoisted per-call invariants for the other-heavy inner loop.
            n_other = len(other_heavy)
            if n_other:
                other_pts = _np.asarray([pts[j] for j in other_heavy], dtype=float)
                other_rcov = [r_cov[j] for j in other_heavy]
            for i, p in frag_pts.items():
                ri = r_cov[i]
                sym_i = mol.GetAtomWithIdx(i).GetSymbol()
                p_arr = _np.asarray(p)
                d = float(_np.linalg.norm(p_arr - m_pos))
                try:
                    ml_ref = float(_get_ml_bond_length(m_sym, sym_i))
                except Exception:
                    ml_ref = 2.2
                thr_m = max(1.20, 0.80 * ml_ref)
                if d < thr_m:
                    score += 5.0 * (thr_m - d) ** 2
                if n_other:
                    # Per-row Euclidean norm, BIT-IDENTICAL to a loop of
                    # float(np.linalg.norm(p - pts[j])): batched matmul uses
                    # the same per-row ddot kernel (verified array_equal).
                    deltas = p_arr - other_pts
                    dists = _np.sqrt(
                        _np.matmul(deltas[:, None, :], deltas[:, :, None]).reshape(-1)
                    )
                    # thr kept in the original scalar 1.3*(ri+rj) form:
                    # distributing the 1.3 would change the FMA/rounding and
                    # break bit-identity, so only the norm is vectorized.
                    for k in range(n_other):
                        d = float(dists[k])
                        thr = 1.3 * (ri + other_rcov[k])
                        if d < thr:
                            score += (thr - d) ** 2
            return score

        def _lp_penalty_for_frag(donors, frag_pts):
            """LP-M misalignment penalty for a monodentate donor.

            Reads the donor's non-metal ring/chain neighbour positions from
            the candidate ``frag_pts`` (the rest of the fragment), builds
            the anti-bisector lone-pair vector, and penalises large angular
            deviation from the D->M direction.  15 deg deadband; quadratic
            in excess so small tilts don't dominate the clash landscape.
            """
            if len(donors) != 1:
                return 0.0
            d_idx = donors[0]
            d_atom = mol.GetAtomWithIdx(d_idx)
            if sum(1 for nb in d_atom.GetNeighbors()
                   if nb.GetSymbol() in _METAL_SET) >= 2:
                return 0.0
            nbrs = [nb.GetIdx() for nb in d_atom.GetNeighbors()
                    if nb.GetSymbol() not in _METAL_SET]
            if not nbrs:
                return 0.0
            d_pos = pts[d_idx]
            lp = _np.zeros(3)
            for nb_idx in nbrs:
                if nb_idx in frag_pts:
                    p = _np.asarray(frag_pts[nb_idx])
                else:
                    p = pts[nb_idx]
                v = p - d_pos
                vn = _np.linalg.norm(v)
                if vn > 1e-8:
                    lp += v / vn
            lp_norm = _np.linalg.norm(lp)
            if lp_norm < 1e-6:
                return 0.0
            lp_unit = -lp / lp_norm
            dm = m_pos - d_pos
            dm_norm = _np.linalg.norm(dm)
            if dm_norm < 1e-6:
                return 0.0
            dm_unit = dm / dm_norm
            cos_t = max(-1.0, min(1.0, float(_np.dot(lp_unit, dm_unit))))
            ang = math.acos(cos_t)
            excess = max(0.0, ang - math.radians(15.0))
            return lp_weight * excess * excess

        def _frag_score(heavy, donors, frag_pts):
            return (
                _clash_for_frag(heavy, frag_pts)
                + _lp_penalty_for_frag(donors, frag_pts)
            )

        def _rotate_about(axis, pivot, angle_deg, atom_pts):
            a = _np.asarray(axis, dtype=float)
            an = float(_np.linalg.norm(a))
            if an < 1e-10:
                return atom_pts
            u = a / an
            c = math.cos(math.radians(angle_deg))
            s = math.sin(math.radians(angle_deg))
            ux, uy, uz = float(u[0]), float(u[1]), float(u[2])
            R = _np.array([
                [c + ux * ux * (1 - c),
                 ux * uy * (1 - c) - uz * s,
                 ux * uz * (1 - c) + uy * s],
                [uy * ux * (1 - c) + uz * s,
                 c + uy * uy * (1 - c),
                 uy * uz * (1 - c) - ux * s],
                [uz * ux * (1 - c) - uy * s,
                 uz * uy * (1 - c) + ux * s,
                 c + uz * uz * (1 - c)],
            ])
            pivot = _np.asarray(pivot, dtype=float)
            out = {}
            for idx, p in atom_pts.items():
                out[idx] = tuple((_np.asarray(p) - pivot) @ R.T + pivot)
            return out

        frag_count = len(frag_info)
        for _pass in range(passes):
            changed = False
            # Alternate fragment order each pass so no single fragment
            # dominates the greedy cascade.  On even passes, forward;
            # odd passes, reverse.  Deterministic, no RNG.
            order = range(frag_count) if _pass % 2 == 0 else range(frag_count - 1, -1, -1)
            for idx in order:
                info = frag_info[idx]
                heavy = info["heavy"]
                donors = info["donors"]
                non_donor = [i for i in heavy if i not in donor_set]
                if not non_donor:
                    continue
                if len(donors) == 1:
                    d = donors[0]
                    axis = pts[d] - m_pos
                    pivot = m_pos
                    n_rot = n_rot_mono
                elif len(donors) == 2:
                    axis = pts[donors[1]] - pts[donors[0]]
                    pivot = 0.5 * (pts[donors[0]] + pts[donors[1]])
                    n_rot = n_rot_bi
                else:
                    continue

                atom_pts = {i: tuple(pts[i]) for i in non_donor}
                best_score = _frag_score(heavy, donors, atom_pts)
                best_pts = atom_pts
                # Multi-level tilt: 50% + 90% of the ideal LP->M correction.
                # Combined with the no-tilt base gives three starting points
                # for the spin-search, so the optimiser can trade clash
                # against LP-alignment in finer steps than binary on/off.
                tilts = _compute_lp_tilt_rotations(donors, pts, m_pos, mol)
                base_sets = [atom_pts]
                for (R_tilt, tilt_pivot) in tilts:
                    tilted = {
                        i: tuple(
                            R_tilt @ (_np.asarray(atom_pts[i]) - tilt_pivot)
                            + tilt_pivot
                        )
                        for i in non_donor
                    }
                    base_sets.append(tilted)
                step = 360.0 / max(1, n_rot)
                for base in base_sets:
                    s0 = _frag_score(heavy, donors, base)
                    if s0 < best_score - 1e-6:
                        best_score = s0
                        best_pts = base
                    for k in range(1, n_rot):
                        angle = step * k
                        candidate = _rotate_about(axis, pivot, angle, base)
                        s = _frag_score(heavy, donors, candidate)
                        if s < best_score - 1e-6:
                            best_score = s
                            best_pts = candidate
                if best_pts is not atom_pts:
                    for i, p in best_pts.items():
                        pts[i] = _np.array(p, dtype=float)
                    changed = True
            if not changed:
                break

        # Write back.
        if is_list:
            for i in range(n_atoms):
                coords[i] = (float(pts[i, 0]), float(pts[i, 1]), float(pts[i, 2]))
        else:
            coords[:] = pts
        return coords
    except Exception as exc:
        logger.debug("_align_and_orient_ligands failed: %s", exc)
        return coords
