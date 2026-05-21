"""HAPTO_RIGID_V2 (2026-05-21) — rigid ring-substituent shell rebuilder.

User-diagnosed + metric-confirmed (hapto_diagnose.py, 0-FP on CCDC): 91.1% of
hapto structures are topologically broken via TWO universal failure modes that
originate in ``_build_hapto_scaffold``:

  (1) Ring-substituent collapse — substituents on RING carbons (both the η-ring
      and pendant organic rings such as the P3N3 cyclophosphazene in VIXXEF) are
      placed COPLANAR with the ring instead of lifted out of the ring plane
      (sp3/sp2 VSEPR).  Concretely ``_place_ring_inline`` lays a ring in one
      plane and the BFS-VSEPR substituent walk leaves each substituent + its H
      in that same plane → a P with two F at the same z as the ring N collapses
      (F-N at 0.84 Å, F bonded to 4), and whole substituents flatten into one
      plane.  Their H collapse (C-H < 0.8 Å) for the same reason.
  (2) σ-donor backbone clash into the η-cone.

The fix is universal and geometric: for EVERY ring (η-coordinated or pendant
organic) we recompute the ring-plane normal from the already-placed ring atoms
and re-grow every NON-ring substituent of every ring atom out of that plane with
correct sp2 (in-plane radial) / sp3 (alternating above/below) geometry; deeper
sub-trees then grow by a VSEPR walk that is biased away from the ring plane and
away from every metal.  A hard non-overlap guard + bond-length spring finish.

Ring atoms themselves, metals, σ-donors (metal-pinned) and ansa bridges stay
frozen — the rigid η-construction the legacy scaffold gets right is preserved.

Gated behind ``DELFIN_HAPTO_RIGID_V2`` (class-conditional hapto/multi_hapto,
default-OFF); bit-exact when off.  Operates purely on the numpy ``coords`` array
+ the RDKit ``mol`` graph already in ``_build_hapto_scaffold`` scope.
"""
from __future__ import annotations

from collections import deque
from typing import Dict, List, Set, Tuple

import numpy as np

_COV = {
    'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'P': 1.07,
    'S': 1.05, 'Cl': 1.02, 'Br': 1.20, 'I': 1.39, 'B': 0.84, 'Si': 1.11,
    'Se': 1.20, 'As': 1.19, 'Te': 1.38, 'Ge': 1.20,
}


def _bond_len(s1: str, s2: str) -> float:
    pair = frozenset([s1, s2])
    if 'H' in pair:
        other = s2 if s1 == 'H' else s1
        return {'C': 1.08, 'N': 1.01, 'O': 0.96, 'Si': 1.48, 'B': 1.19,
                'P': 1.42, 'S': 1.34}.get(other, 1.08)
    if pair == frozenset(['Si', 'C']):
        return 1.87
    if pair == frozenset(['C', 'C']):
        return 1.50
    if pair == frozenset(['C', 'N']):
        return 1.47
    if pair == frozenset(['C', 'O']):
        return 1.43
    if pair == frozenset(['P', 'N']):
        return 1.58
    if pair == frozenset(['P', 'F']):
        return 1.55
    if pair == frozenset(['P', 'C']):
        return 1.84
    return _COV.get(s1, 0.76) + _COV.get(s2, 0.76)


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    return v / n if n > 1e-9 else v


def _ortho(n: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    n = _unit(n)
    ref = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = _unit(np.cross(n, ref))
    v = _unit(np.cross(n, u))
    return u, v


_MAX_DEG = {"C": 4, "N": 4, "O": 3, "S": 6, "P": 6, "F": 1, "Cl": 1, "Br": 1,
            "I": 1, "B": 4, "Si": 4, "Se": 4, "As": 5, "Te": 6, "H": 1}
_METALS = set(
    "Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os "
    "Ir Pt Au Hg La Ce Nd Sm Gd Lu U Th".split()
)


def count_topology_breaks(mol, conf_id: int = 0) -> int:
    """Count the headline topology breaks (mirrors hapto_diagnose.py):
    per-atom heavy over-coordination + heavy-heavy collapse + X-H collapse +
    superposition.  Used to decide whether the rigid scaffold is strictly
    better than the alternative candidate, so prefer-scaffold only fires when
    it does not regress.  Geometric (no SMILES needed)."""
    try:
        conf = mol.GetConformer(conf_id)
    except Exception:
        return 10 ** 6
    n = mol.GetNumAtoms()
    syms = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(n)]
    P = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                   conf.GetAtomPosition(i).z] for i in range(n)], dtype=float)
    # geometric heavy-bond graph
    hadj: Dict[int, List[int]] = {i: [] for i in range(n)}
    breaks = 0
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(P[i] - P[j]))
            si, sj = syms[i], syms[j]
            if si in _METALS or sj in _METALS:
                continue
            ideal = _bond_len(si, sj)
            if d < 1.30 * ideal:  # geometric bond
                if si != "H" and sj != "H":
                    hadj[i].append(j); hadj[j].append(i)
                    if d < 0.70 * ideal:
                        breaks += 1  # heavy_collapse
                else:
                    if d < 0.80:
                        breaks += 1  # xh_collapse
            if d < 0.5:
                breaks += 1  # superposition
    for c in range(n):
        if syms[c] in _METALS or syms[c] == "H":
            continue
        mx = _MAX_DEG.get(syms[c], 4)
        if len(hadj[c]) > mx:
            breaks += 1  # topology_overcoord
    return breaks


def count_topology_breaks_xyz(xyz_block: str) -> int:
    """Headline topology breaks (over-coord + heavy/X-H collapse + superposition)
    on a DELFIN coordinate block (header-less ``sym x y z`` lines).  Used to
    filter σ-isomers that the OB-UFF leg folded worse than the primary frame."""
    syms: List[str] = []
    P: List[List[float]] = []
    for ln in str(xyz_block).splitlines():
        p = ln.split()
        if len(p) >= 4:
            try:
                P.append([float(p[1]), float(p[2]), float(p[3])]); syms.append(p[0])
            except ValueError:
                pass
    n = len(syms)
    if n < 3:
        return 0
    Pa = np.array(P, dtype=float)
    hadj: Dict[int, List[int]] = {i: [] for i in range(n)}
    breaks = 0
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(Pa[i] - Pa[j]))
            si, sj = syms[i], syms[j]
            if si in _METALS or sj in _METALS:
                continue
            ideal = _bond_len(si, sj)
            if d < 1.30 * ideal:
                if si != "H" and sj != "H":
                    hadj[i].append(j); hadj[j].append(i)
                    if d < 0.70 * ideal:
                        breaks += 1
                elif d < 0.80:
                    breaks += 1
            if d < 0.5:
                breaks += 1
    for c in range(n):
        if syms[c] in _METALS or syms[c] == "H":
            continue
        if len(hadj[c]) > _MAX_DEG.get(syms[c], 4):
            breaks += 1
    return breaks


def _vsepr_open_dirs(used: List[np.ndarray], n_open: int, degree: int,
                     anti_metal: np.ndarray) -> List[np.ndarray]:
    """Return ``n_open`` unit directions that fill the VSEPR slots of a centre
    whose total coordination is ``degree``, given the already-occupied unit
    directions ``used``.  Siblings are spread apart (no collapse); a mild
    anti-metal bias orients otherwise-free frames away from the metal."""
    # ideal X-center-X angle from the steric number (= degree)
    if degree <= 2:
        ang = np.radians(180.0)
    elif degree == 3:
        ang = np.radians(120.0)
    else:
        ang = np.radians(109.5)

    out: List[np.ndarray] = []
    cur = [u.copy() for u in used]

    def _pick_axis_for_one(d0: np.ndarray) -> np.ndarray:
        ref = (np.array([1.0, 0.0, 0.0]) if abs(d0[0]) < 0.9
               else np.array([0.0, 1.0, 0.0]))
        perp = _unit(np.cross(d0, ref))
        if float(np.linalg.norm(anti_metal)) > 1e-8:
            # rotate perp toward the anti-metal half-space
            if float(np.dot(perp, anti_metal)) < 0:
                perp = -perp
        return -d0 * np.cos(np.pi - ang) + perp * np.sin(np.pi - ang)

    for _ in range(n_open):
        if not cur:
            # no occupied dirs yet: first one points along anti-metal (or z)
            d = anti_metal if float(np.linalg.norm(anti_metal)) > 1e-8 \
                else np.array([0.0, 0.0, 1.0])
            d = _unit(d)
        elif len(cur) == 1:
            d = _unit(_pick_axis_for_one(cur[0]))
        elif len(cur) == 2:
            # third slot: perpendicular to the plane of the two, on the side
            # opposite their sum (and biased away from metal)
            s = cur[0] + cur[1]
            sl = float(np.linalg.norm(s))
            if sl > 1e-3:
                bis = s / sl
                perp = _unit(np.cross(cur[0], cur[1]))
                # tetrahedral: tilt the new bond away from the bisector
                d = _unit(-bis * np.cos(np.radians(54.75))
                          + perp * np.sin(np.radians(54.75)))
                if float(np.linalg.norm(anti_metal)) > 1e-8 \
                        and float(np.dot(d, anti_metal)) < 0:
                    d = _unit(-bis * np.cos(np.radians(54.75))
                              - perp * np.sin(np.radians(54.75)))
            else:
                d = _unit(np.cross(cur[0], cur[1]))
        else:
            avg = sum(cur) / len(cur)
            al = float(np.linalg.norm(avg))
            d = -avg / al if al > 1e-8 else _unit(anti_metal
                                                  if float(np.linalg.norm(anti_metal)) > 1e-8
                                                  else np.array([0.0, 0.0, 1.0]))
        out.append(d)
        cur.append(d)
    return out


def rebuild_substituent_shells(
    mol,
    coords: np.ndarray,
    hapto_groups: List[Tuple[int, List[int]]],
    metal_indices: List[int],
    all_hapto_atoms: Set[int],
    all_bridge_atoms: Set[int],
    metal_set,
) -> int:
    """Re-grow non-ring substituents of every ring atom out of the ring plane and
    propagate deeper sub-trees away from ring/metal.  Mutates ``coords`` in
    place; returns number of atoms re-placed."""
    try:
        from rdkit import Chem
        Chem.FastFindRings(mol)
        ri = mol.GetRingInfo()
        atom_rings = [tuple(r) for r in ri.AtomRings()]
    except Exception:
        atom_rings = []

    n_atoms = mol.GetNumAtoms()
    metal_pos = {mi: coords[mi].copy() for mi in metal_indices}

    # FROZEN anchors: metals, η-ring atoms, σ-donors (metal-pinned), ansa bridges
    sigma_donors: Set[int] = set()
    for mi in metal_indices:
        for nbr in mol.GetAtomWithIdx(mi).GetNeighbors():
            ni = nbr.GetIdx()
            if ni not in all_hapto_atoms and nbr.GetSymbol() not in metal_set:
                sigma_donors.add(ni)
    frozen: Set[int] = (set(metal_indices) | set(all_hapto_atoms)
                        | set(all_bridge_atoms) | set(sigma_donors))

    # Ring set the substituents hang off.  Two sources:
    #   (a) The η-groups themselves — RDKit's GetRingInfo FRAGMENTS a metal-
    #       coordinated η-ring into many tiny metal-fused pseudo-rings (Ru bonded
    #       to every ring C), so it cannot be recovered from atom_rings.  We add
    #       each hapto group explicitly with normal = metal→centroid axis; these
    #       atoms stay FROZEN but their substituents are re-grown out of the
    #       ring plane (the #1 user failure mode).
    #   (b) Pendant organic rings (no metal) from GetRingInfo.
    ring_of: Dict[int, List[int]] = {}
    rings: List[Tuple[Tuple[int, ...], np.ndarray, np.ndarray]] = []
    _seen_ring_keys: Set[frozenset] = set()

    for _mi, grp in hapto_groups:
        if len(grp) < 2:
            continue
        key = frozenset(grp)
        if key in _seen_ring_keys:
            continue
        _seen_ring_keys.add(key)
        pts = np.array([coords[a] for a in grp], dtype=float)
        cen = pts.mean(axis=0)
        ax = cen - coords[_mi]
        normal = _unit(ax) if float(np.linalg.norm(ax)) > 1e-8 \
            else np.array([0.0, 0.0, 1.0])
        rings.append((tuple(grp), cen, normal))
        for a in grp:
            ring_of.setdefault(a, []).append(len(rings) - 1)

    for r in atom_rings:
        if any(mol.GetAtomWithIdx(a).GetSymbol() in metal_set for a in r):
            continue
        key = frozenset(r)
        if key in _seen_ring_keys:
            continue
        _seen_ring_keys.add(key)
        pts = np.array([coords[a] for a in r], dtype=float)
        cen = pts.mean(axis=0)
        try:
            _, _, vh = np.linalg.svd(pts - cen)
            normal = _unit(vh[2])
        except Exception:
            normal = np.array([0.0, 0.0, 1.0])
        rings.append((r, cen, normal))
        for a in r:
            ring_of.setdefault(a, []).append(len(rings) - 1)

    placed: Set[int] = set(frozen) | set(a for r, _, _ in rings for a in r)
    replaced: List[int] = []
    queue: deque = deque()

    def _anti_metal(pos: np.ndarray) -> np.ndarray:
        b = np.zeros(3)
        for mi in metal_indices:
            v = pos - metal_pos[mi]
            d = float(np.linalg.norm(v))
            if 1e-8 < d < 4.5:
                b += (v / d) * (4.5 - d)
        return _unit(b) if float(np.linalg.norm(b)) > 1e-8 else np.zeros(3)

    # ---- 1. re-place non-ring substituents of every ring atom, out-of-plane --
    # Universal VSEPR: a ring atom's substituents fill the directions NOT taken
    # by its already-placed neighbours (the two ring members AND any frozen
    # exocyclic attachment such as a Cp-ring carbon or σ-donor).  For a single
    # substituent on an sp2 ring atom this is the in-plane radial bisector; for
    # sp3 ring atoms (e.g. a P bearing 2 F) the substituents split above/below
    # the ring plane.  Crucially this counts the exocyclic attachment so the
    # substituent never points back toward the metal/ring it is fused to.
    for ridx, (r, cen, normal) in enumerate(rings):
        rset = set(r)
        for a in r:
            atom = mol.GetAtomWithIdx(a)
            apos = coords[a]
            sym_a = atom.GetSymbol()
            arom = atom.GetIsAromatic()
            # in-plane radial (outward bisector of the two ring neighbours)
            ring_nbr_pos = [coords[nb.GetIdx()] for nb in atom.GetNeighbors()
                            if nb.GetIdx() in rset]
            if len(ring_nbr_pos) >= 2:
                inn = sum(_unit(rp - apos) for rp in ring_nbr_pos)
                radial = _unit(-inn)
            else:
                radial = _unit(apos - cen)
            if float(np.linalg.norm(radial)) < 1e-8:
                radial = _unit(apos - cen)
            # directions already occupied by PLACED neighbours (ring + frozen
            # exocyclic): substituents must avoid these.
            occ = []
            for nb in atom.GetNeighbors():
                nj = nb.GetIdx()
                if nj in placed or nj in rset:
                    occ.append(_unit(coords[nj] - apos))
            subs = [nb.GetIdx() for nb in atom.GetNeighbors()
                    if nb.GetIdx() not in rset
                    and nb.GetIdx() not in metal_indices
                    and nb.GetIdx() not in frozen
                    and nb.GetIdx() not in placed]
            ns = len(subs)
            for k, ni in enumerate(subs):
                sym_s = mol.GetAtomWithIdx(ni).GetSymbol()
                bl = _bond_len(sym_a, sym_s)
                if ns == 1 and arom:
                    direction = radial
                elif ns == 1 and occ:
                    # one substituent: opposite the mean of occupied dirs.
                    avg = sum(occ) / len(occ)
                    al = float(np.linalg.norm(avg))
                    base = -avg / al if al > 1e-8 else radial
                    # If the centre is (near-)sp3 with >=3 coplanar occupied
                    # neighbours, the 4th substituent MUST leave the ring plane
                    # (else it lands in the ring centre).  Force an out-of-plane
                    # lift in that case; otherwise nudge toward the radial.
                    if len(occ) >= 3 and abs(float(np.dot(base, normal))) < 0.30:
                        s_lift = 1.0 if float(np.dot(radial, normal)) >= 0 else -1.0
                        direction = _unit(base + 0.9 * s_lift * normal)
                    else:
                        direction = _unit(base + 0.3 * radial)
                elif ns == 1:
                    direction = _unit(radial + 0.5 * normal)
                else:
                    # multiple substituents (PF2, CMe2, …): split above/below
                    # the ring plane, splayed radially outward (tetrahedral).
                    sign = 1.0 if (k % 2 == 0) else -1.0
                    lift = np.radians(54.0)
                    direction = _unit(radial * np.cos(lift)
                                      + sign * normal * np.sin(lift))
                coords[ni] = apos + bl * _unit(direction)
                placed.add(ni)
                replaced.append(ni)
                # enqueue ni's children (ni itself is now placed); the walk
                # grows the rest of the substituent sub-tree from here.
                for cn in mol.GetAtomWithIdx(ni).GetNeighbors():
                    cj = cn.GetIdx()
                    if cj not in placed and cn.GetSymbol() not in metal_set:
                        queue.append((cj, ni))

    # ---- 2. seed σ-donor backbones away from cones ----
    for d_idx in sigma_donors:
        for nbr in mol.GetAtomWithIdx(d_idx).GetNeighbors():
            ni = nbr.GetIdx()
            if ni in placed or ni in metal_indices:
                continue
            queue.append((ni, d_idx))

    # ---- 3. VSEPR walk for deeper non-ring sub-trees ----
    pending: deque = deque(queue)
    guard = 0
    while pending and guard < n_atoms * 6:
        guard += 1
        ai, pi = pending.popleft()
        if ai in placed:
            continue
        # if this atom is part of a ring, it was already placed; skip
        if ai in ring_of:
            placed.add(ai)
            for nbr in mol.GetAtomWithIdx(ai).GetNeighbors():
                nj = nbr.GetIdx()
                if nj not in placed and nbr.GetSymbol() not in metal_set:
                    pending.append((nj, ai))
            continue
        if pi not in placed:
            pending.append((ai, pi))
            continue
        ppos = coords[pi]
        sym_p = mol.GetAtomWithIdx(pi).GetSymbol()
        # Place ALL currently-unplaced substituents of this parent JOINTLY so
        # siblings (e.g. the three H of a methyl, gem-dimethyl) get distinct
        # tetrahedral/trigonal directions instead of greedily collapsing onto
        # one slot.
        siblings = [ai]
        for pn in mol.GetAtomWithIdx(pi).GetNeighbors():
            pj = pn.GetIdx()
            if pj == ai or pj in placed:
                continue
            if pn.GetSymbol() in metal_set:
                continue
            siblings.append(pj)
        used = []
        for pn in mol.GetAtomWithIdx(pi).GetNeighbors():
            pj = pn.GetIdx()
            if pj in placed and pj not in siblings:
                used.append(_unit(coords[pj] - ppos))
        for mi in metal_indices:
            mv = ppos - metal_pos[mi]
            ml = float(np.linalg.norm(mv))
            if 1e-8 < ml < 4.0:
                used.append(_unit(-mv))
        deg = mol.GetAtomWithIdx(pi).GetDegree()
        n_open = len(siblings)
        # Build n_open new directions filling the VSEPR slots left by `used`.
        new_dirs = _vsepr_open_dirs(used, n_open, deg, _anti_metal(ppos))
        for sib, direction in zip(siblings, new_dirs):
            sym_c = mol.GetAtomWithIdx(sib).GetSymbol()
            bl = _bond_len(sym_p, sym_c)
            coords[sib] = ppos + bl * _unit(direction)
            placed.add(sib)
            replaced.append(sib)
            for nbr in mol.GetAtomWithIdx(sib).GetNeighbors():
                nj = nbr.GetIdx()
                if nj not in placed and nbr.GetSymbol() not in metal_set:
                    pending.append((nj, sib))

    if not replaced:
        return 0
    replaced_set = set(replaced)

    # ---- 4. bond-length spring on newly placed sub-trees ----
    bond_targets: List[Tuple[int, int, float]] = []
    for bond in mol.GetBonds():
        bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        si = mol.GetAtomWithIdx(bi).GetSymbol()
        sj = mol.GetAtomWithIdx(bj).GetSymbol()
        if si in metal_set or sj in metal_set:
            continue
        if bi not in replaced_set and bj not in replaced_set:
            continue
        bond_targets.append((bi, bj, _bond_len(si, sj)))
    for _ in range(80):
        forces = np.zeros_like(coords)
        max_err = 0.0
        for bi, bj, tb in bond_targets:
            diff = coords[bj] - coords[bi]
            d = float(np.linalg.norm(diff))
            if d < 1e-8:
                diff = np.array([0.1, 0.2, 0.3]); d = float(np.linalg.norm(diff))
            err = d - tb
            max_err = max(max_err, abs(err))
            f = 0.3 * err * diff / d
            forces[bi] += f
            forces[bj] -= f
        if max_err < 0.03:
            break
        for ai in replaced:
            coords[ai] = coords[ai] + forces[ai]

    # ---- 5. hard non-overlap guard — runs LAST so the bond-spring cannot
    # re-overlap sibling atoms (e.g. two methyl H pulled back together).  Floor
    # is generous for H-H so sibling hydrogens never count as a superposition.
    #
    # Bond-preservation contract (Iter-26 H-fix): a terminal hydrogen (single
    # heavy neighbour) MUST stay on its parent's bond sphere — a free
    # translation that "fixes" an H-vs-C clash by pushing the H back through its
    # own parent SHORTENS the C-H bond, and the downstream OB-UFF then re-extends
    # it in a now-wrong direction straight into the clashing atom (the observed
    # xh_collapse: C-H pulled to ~0.9 Å here, OB-UFF re-grows it to 1.08 Å onto
    # the neighbour at 0.42 Å).  For such hydrogens we therefore SLIDE the H
    # tangentially on the bond sphere (rotate about the parent) instead of
    # translating, and re-project to the exact bond length every step.  Heavy
    # atoms keep the original free-push behaviour.
    parent_of: Dict[int, int] = {}
    for ai in replaced_set:
        a = mol.GetAtomWithIdx(ai)
        heavy_nbrs = [nb.GetIdx() for nb in a.GetNeighbors()
                      if nb.GetSymbol() != "H" and nb.GetSymbol() not in metal_set]
        if a.GetSymbol() == "H" and len(heavy_nbrs) == 1:
            parent_of[ai] = heavy_nbrs[0]

    def _slide_terminal_h(h_idx: int, p_idx: int, away: np.ndarray,
                          bl_h: float) -> bool:
        """Rotate terminal H around its parent so it moves toward ``away``
        (a unit push direction) while preserving the parent–H bond length.
        Returns True if it actually moved."""
        bond = coords[h_idx] - coords[p_idx]
        r = float(np.linalg.norm(bond))
        if r < 1e-6:
            # degenerate: re-seed along the away direction
            coords[h_idx] = coords[p_idx] + bl_h * _unit(away if
                float(np.linalg.norm(away)) > 1e-8 else np.array([0.0, 0.0, 1.0]))
            return True
        bdir = bond / r
        # tangential component of the desired push (perp to the bond)
        tang = away - float(np.dot(away, bdir)) * bdir
        tl = float(np.linalg.norm(tang))
        if tl < 1e-6:
            # push is collinear with the bond: pick an arbitrary perpendicular
            ref = (np.array([1.0, 0.0, 0.0]) if abs(bdir[0]) < 0.9
                   else np.array([0.0, 1.0, 0.0]))
            tang = _unit(np.cross(bdir, ref))
        else:
            tang = tang / tl
        # rotate bond direction toward the tangent by a fixed small step
        step = np.radians(18.0)
        new_dir = _unit(bdir * np.cos(step) + tang * np.sin(step))
        coords[h_idx] = coords[p_idx] + bl_h * new_dir
        return True

    for _ in range(120):
        moved = False
        for ai in replaced:
            sa = mol.GetAtomWithIdx(ai).GetSymbol()
            for bj in range(n_atoms):
                if bj == ai:
                    continue
                sb = mol.GetAtomWithIdx(bj).GetSymbol()
                if sa in metal_set or sb in metal_set:
                    continue
                if mol.GetBondBetweenAtoms(ai, bj) is not None:
                    continue
                diff = coords[ai] - coords[bj]
                d = float(np.linalg.norm(diff))
                # H-H needs a wider floor than 0.80·Σcov(=0.50) to clear the
                # superposition (0.5 Å) + xh-clash thresholds.
                if sa == "H" and sb == "H":
                    floor = 1.40
                elif sa == "H" or sb == "H":
                    floor = 0.90 * (_COV.get(sa, 0.76) + _COV.get(sb, 0.76))
                else:
                    floor = 0.80 * (_COV.get(sa, 0.76) + _COV.get(sb, 0.76))
                if d < floor:
                    away = (diff / d if d > 1e-8
                            else np.array([0.31, 0.11, 0.21]))
                    if ai in parent_of:
                        # bond-preserving tangential slide
                        p_idx = parent_of[ai]
                        bl_h = _bond_len(mol.GetAtomWithIdx(p_idx).GetSymbol(),
                                         sa)
                        _slide_terminal_h(ai, p_idx, away, bl_h)
                    else:
                        coords[ai] = coords[ai] + (floor - d) * 0.6 * away
                    moved = True
        if not moved:
            break

    # ---- 6. final terminal-H bond snap — guarantees every re-placed terminal
    # hydrogen sits at the exact element bond length from its parent (no residual
    # under/over-stretch the downstream OB-UFF could amplify into a collapse).
    for h_idx, p_idx in parent_of.items():
        bl_h = _bond_len(mol.GetAtomWithIdx(p_idx).GetSymbol(),
                         mol.GetAtomWithIdx(h_idx).GetSymbol())
        bond = coords[h_idx] - coords[p_idx]
        r = float(np.linalg.norm(bond))
        if r > 1e-6:
            coords[h_idx] = coords[p_idx] + bl_h * (bond / r)

    return len(replaced_set)


def constrained_uff_relax(
    mol,
    coords: np.ndarray,
    hapto_groups: List[Tuple[int, List[int]]],
    metal_indices: List[int],
    all_hapto_atoms: Set[int],
    all_bridge_atoms: Set[int],
    metal_set,
    max_iters: int = 400,
) -> int:
    """Constrained UFF relaxation of the rigid-V2 build (user-proposed 2026-05-21).

    The rigid placement gets TOPOLOGY (connectivity + H) right but leaves residual
    inter-/intra-ligand vdW CLASHES and bond-length DISTORTION.  This pass resolves
    them with a UFF minimisation that holds the COORDINATION CORE rigid while only
    the ligand shell relaxes.  It is fully universal / graph-only (no per-refcode
    fitting) and deterministic.

    Mechanism (the three fundamental rules the user identified):
      (1) EXCLUDE the metal from the force field.  UFF cannot type unparametrised
          4d/5d transition metals — leaving them in injects garbage gradients.  We
          build a metal-FREE editable copy of ``mol`` (delete every metal atom),
          so the FF only sees the organic ligand graph.
      (2) FREEZE the coordination core: metal (already removed) + every η-ring atom
          + every σ-donor (the atoms directly bonded to a metal) + ansa-bridge
          atoms.  These are added as UFF fixed points so the rigid η-cone and the
          metal-pinned donor positions cannot drift.  Because the donors are frozen
          and the metal is never moved, every M–donor distance is preserved exactly
          (the user's "pin M-D" requirement is met implicitly + exactly, with no
          spurious metal force-field terms).
      (3) RELAX the rest.  The free ligand atoms minimise against UFF vdW + bonded
          (bond/angle/torsion) terms, so inter-ligand and intra-ligand clashes
          resolve and stretched/compressed bonds normalise — WITHOUT moving the
          frozen core.

    Mutates ``coords`` in place; returns the number of free (relaxed) atoms, or 0
    if the relaxation could not run (the rigid build is then left untouched).
    Deterministic: UFF minimisation from a fixed start with fixed points is a
    deterministic descent; no randomness is introduced.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception:
        return 0

    n_atoms = mol.GetNumAtoms()
    if n_atoms < 4:
        return 0

    # --- scope guard: SAFE only for a single coordination centre --------------
    # Deleting the metal to build the FF cuts every metal bond.  For systems with
    # more than one metal — or a main-group bridge metal (Sn/Pb/Ge/Sb/Bi acting
    # as a metalloligand, e.g. a Cl3Sn–Ru bridge) — that fragmentation can
    # DISCONNECT the ligand graph into pieces only held together through the
    # deleted metal(s).  UFF then has no term coupling those pieces and the
    # downstream per-isomer OB-UFF amplifies the perturbation into collapses
    # (validated failure: HUCPIH, a Cl3Sn–Ru η5).  Restrict the relaxation to the
    # single-metal case where the organic ligand graph stays connected.
    _MAIN_GROUP_BRIDGE = {"Sn", "Pb", "Ge", "Sb", "Bi"}
    if len(metal_indices) != 1:
        return 0
    for a in range(n_atoms):
        if a in set(metal_indices):
            continue
        if mol.GetAtomWithIdx(a).GetSymbol() in _MAIN_GROUP_BRIDGE:
            # a main-group atom directly bonded to the metal = metalloligand
            # bridge → unsafe fragmentation.
            for nbr in mol.GetAtomWithIdx(a).GetNeighbors():
                if nbr.GetIdx() in set(metal_indices):
                    return 0

    # --- frozen coordination core (metal handled separately via deletion) ------
    sigma_donors: Set[int] = set()
    metal_idx_set = set(metal_indices)
    for mi in metal_indices:
        for nbr in mol.GetAtomWithIdx(mi).GetNeighbors():
            ni = nbr.GetIdx()
            if nbr.GetSymbol() not in metal_set:
                sigma_donors.add(ni)
    frozen_core: Set[int] = (set(all_hapto_atoms) | set(all_bridge_atoms)
                             | sigma_donors)
    # also freeze every atom the converter flagged as hapto via the groups
    for _mi, grp in hapto_groups:
        for a in grp:
            frozen_core.add(a)
    frozen_core -= metal_idx_set  # metals are removed, not frozen

    # --- build a metal-FREE editable copy + index map --------------------------
    rw = Chem.RWMol(mol)
    # remove metals in descending index order so earlier indices stay valid
    for mi in sorted(metal_idx_set, reverse=True):
        rw.RemoveAtom(mi)
    frag = rw.GetMol()
    # old index -> new index (metals dropped)
    old_to_new: Dict[int, int] = {}
    new_to_old: Dict[int, int] = {}
    cur = 0
    for old in range(n_atoms):
        if old in metal_idx_set:
            continue
        old_to_new[old] = cur
        new_to_old[cur] = old
        cur += 1
    n_frag = frag.GetNumAtoms()
    if n_frag < 3:
        return 0

    # The fragment may have radicals/odd valences where metal bonds were cut.
    # UFF does not need a clean sanitise; it only needs ring perception + a
    # conformer.  Use a permissive partial sanitise so typing succeeds.
    try:
        Chem.SanitizeMol(
            frag,
            sanitizeOps=(Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                         | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                         | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                         | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
                         | Chem.SanitizeFlags.SANITIZE_SYMMRINGS),
            catchErrors=True,
        )
    except Exception:
        pass

    # --- conformer from current coords (only non-metal atoms) ------------------
    from rdkit.Geometry import Point3D
    conf = Chem.Conformer(n_frag)
    for new in range(n_frag):
        old = new_to_old[new]
        conf.SetAtomPosition(new, Point3D(float(coords[old, 0]),
                                          float(coords[old, 1]),
                                          float(coords[old, 2])))
    frag.RemoveAllConformers()
    frag.AddConformer(conf, assignId=True)

    if not AllChem.UFFHasAllMoleculeParams(frag):
        # a non-organic ligand atom (e.g. a main-group donor) is untypable;
        # bail rather than relax against garbage.
        return 0

    try:
        ff = AllChem.UFFGetMoleculeForceField(frag)
    except Exception:
        return 0
    if ff is None:
        return 0
    ff.Initialize()

    # --- freeze the coordination core (fixed points) ---------------------------
    n_free = 0
    for new in range(n_frag):
        old = new_to_old[new]
        if old in frozen_core:
            ff.AddFixedPoint(new)
        else:
            n_free += 1
    if n_free == 0:
        return 0

    # snapshot the pre-relaxation coordinates so we can REVERT if the relaxation
    # makes the structure worse (break-aware accept/reject, below).
    syms_all = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(n_atoms)]
    coords_before = coords.copy()
    breaks_before = _diagnose_breaks(syms_all, coords)

    # deterministic descent; fixed points hold the core, free atoms relax.
    try:
        ff.Minimize(maxIts=int(max_iters), forceTol=1e-4, energyTol=1e-6)
    except Exception:
        return 0

    # --- write relaxed positions back (free atoms only; core is unchanged) -----
    relaxed_conf = frag.GetConformer()
    for new in range(n_frag):
        old = new_to_old[new]
        if old in frozen_core:
            continue
        p = relaxed_conf.GetAtomPosition(new)
        coords[old, 0] = p.x
        coords[old, 1] = p.y
        coords[old, 2] = p.z

    # --- break-aware accept/reject (safe-by-construction) ----------------------
    # The metal-free fragmentation cuts every metal bond; for bi-/poly-metallic
    # or main-group-donor systems (e.g. a Sn–Ru bridge) the fragment graph can be
    # disconnected, and the UFF minimum then collapses ligand pieces into each
    # other.  Universal guard: keep the relaxation ONLY if it did not INCREASE the
    # headline topology breaks (over-coord + heavy/X-H collapse + superposition).
    # Otherwise revert to the rigid-V2 build — the relaxation is never a
    # net-topology regression, by construction.
    breaks_after = _diagnose_breaks(syms_all, coords)
    if breaks_after > breaks_before:
        coords[:] = coords_before
        return 0
    return n_free


def _diagnose_breaks(syms: List[str], P: np.ndarray) -> int:
    """Authoritative structural-break count, using the SAME geometric bond model
    as the calibration metric (``delfin._bond_decollapse``) so the accept/reject
    guard agrees exactly with the 0-FP-on-CCDC diagnostic.  Sums the topology +
    distortion criteria the relaxation could plausibly affect: heavy
    over-coordination, heavy-heavy collapse, X-H collapse, bond-length
    distortion, H over-coordination, and superposition.  Metal bonds excluded
    (consistent with the metric)."""
    try:
        import delfin._bond_decollapse as _bd
    except Exception:
        return 0
    n = len(syms)
    if n < 3:
        return 0
    b = _bd._geometric_bonds(syms, P)
    hadj: Dict[int, List[int]] = {i: [] for i in range(n)}
    for i, j in b:
        if syms[i] != "H" and syms[j] != "H":
            hadj[i].append(j); hadj[j].append(i)
    brk = 0
    # heavy over-coordination
    for c in range(n):
        if _bd._is_metal(syms[c]) or syms[c] == "H":
            continue
        if len(hadj[c]) > _MAX_DEG.get(syms[c], 4):
            brk += 1
    # bonded collapse + bond-length distortion
    for i, j in b:
        si, sj = syms[i], syms[j]
        if _bd._is_metal(si) or _bd._is_metal(sj):
            continue
        d = float(np.linalg.norm(P[i] - P[j]))
        if si == "H" or sj == "H":
            if d < 0.80:
                brk += 1                       # xh_collapse
            continue
        ideal = _bd._ideal_bond(si, sj)
        if d < 0.70 * ideal:
            brk += 1                           # heavy_collapse
        dev = (d - ideal) / ideal
        if dev < -0.19 or dev > 0.03:
            brk += 1                           # bond_distort
    # H over-coordination (H within bonding range of >=2 heavy non-metal atoms)
    for h in range(n):
        if syms[h] != "H":
            continue
        nb = sum(1 for j in range(n) if j != h and syms[j] != "H"
                 and not _bd._is_metal(syms[j])
                 and float(np.linalg.norm(P[h] - P[j])) < 1.3 * _bd._ideal_bond(syms[h], syms[j]))
        if nb >= 2:
            brk += 1
    # superposition (any pair < 0.5 A)
    for i in range(n):
        for j in range(i + 1, n):
            if float(np.linalg.norm(P[i] - P[j])) < 0.5:
                brk += 1
    return brk
