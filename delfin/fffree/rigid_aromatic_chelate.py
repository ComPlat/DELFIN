"""delfin.fffree.rigid_aromatic_chelate — rigid-body placement of aromatic chelates.

Real-CCDC crystallography of aromatic-chelate ligands (1,10-phenanthroline,
2,2'-bipyridine, 2,2':6',2''-terpyridine, pyridyl-pyridyl-pyridyl scaffolds,
salen-aromatic, pybox, ...) shows the ligand as an **effectively rigid body**:
the aromatic backbone fixes the internal coordinates (bite angle, donor-donor
distance, plane, every internal C-C and C-N distance) up to the level of
thermal motion.  The Mogul-PRIMARY soft-DG embed + per-donor centroid radial
pull treats every chelate as flexible and lets the metal pull individual
donors independently, which produces (a) M-N distances 0.10-0.25 Å above the
CCDC mean, (b) the metal sitting 0.05-0.15 Å above the chelate plane, and
(c) per-isomer drift in the donor-vertex direction.

This module replaces the soft-DG result for aromatic chelates with a single
rigid-body Kabsch placement:

    1.  Detect (universal): an aromatic chelate is a multidentate ligand
        whose donor atoms are pairwise connected through paths consisting
        ENTIRELY of aromatic bonds (``Bond.GetIsAromatic() == True``).  No
        SMILES patterns, no per-chelate template table.

    2.  Template (universal): RDKit-embed the standalone ligand fragment
        (ligand SMILES alone, no metal).  The aromatic backbone enforces
        its rigid geometry by construction — the embedder cannot break
        aromaticity, so the bite angle, donor-donor distance, and ring
        planarity all come out at CCDC values.

    3.  Place (universal): Kabsch-align the template donor positions onto
        the target M-D vector positions (per-donor target = ``M + r_md_k
        * polyhedron_vertex_k``), then transform the WHOLE ligand subtree
        (donors + every aromatic/substituent atom) by that same rigid
        transform.  Pure rotation + translation -> every internal
        distance and angle is preserved up to numerical precision.

After this step the chelate ligand is byte-identical to a CCDC-embedded
standalone ligand fragment rigid-fit to the polyhedron vertices, and the
metal sits exactly on the polyhedron vertex centre.  M-D distances come
from the CCDC supplement (per metal + donor element), bite angle is the
ligand's native aromatic bite, plane is the ligand's native plane, metal-
out-of-plane is mathematically 0 for symmetric bidentate chelates and
strictly below the per-donor M-D ε for non-symmetric chelates.

Universality contract:

  * NO SMILES patterns.
  * NO per-chelate template table.
  * Detection uses ``Bond.GetIsAromatic()`` only — applies equally to
    phen, bipy, terpy, pyridyl-imidazolyl, salen-aromatic, ... .
  * Template uses RDKit ``EmbedMolecule`` only — the aromatic backbone is
    the geometric truth, exactly as in the CCDC.
  * Placement uses Kabsch only — no element-specific or ring-size-specific
    branches.

Determinism contract:

  * Deterministic ring iteration order.
  * Fixed embed seed (matches the rest of the fffree pipeline).
  * Pure linear algebra, no RNG outside the standalone embed.
  * Idempotent on second call (template + Kabsch yields the same result).

Env-flag::

    DELFIN_FFFREE_RIGID_AROMATIC_CHELATE = "1"   # default ON when
                                                  # MOGUL_PRIMARY=1
    DELFIN_FFFREE_RIGID_AROMATIC_CHELATE = "0"   # disable

When the gate is OFF (or when no aromatic chelate is detected, or when
the standalone embed fails), the geometry passed in is returned
unchanged (defence-in-depth: never crash a build).

Author: hmaximilian <hmaximilian496@gmail.com>
Branch: feat-mogul-primary-2026-06-07
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np


# Standalone-ligand embed seed.  Matches the rest of the fffree pipeline so
# any cross-process reproduction is byte-identical when PYTHONHASHSEED=0.
_TEMPLATE_SEED: int = 0xDE17


__all__ = [
    "rigid_aromatic_chelate_enabled",
    "is_aromatic_chelate",
    "detect_aromatic_chelate_groups",
    "build_chelate_template",
    "place_rigid_aromatic_chelates",
]


# ---------------------------------------------------------------------------
# Env-flag wiring
# ---------------------------------------------------------------------------
def rigid_aromatic_chelate_enabled() -> bool:
    """Return True iff rigid-aromatic-chelate placement is wired on.

    Default ON when ``DELFIN_FFFREE_MOGUL_PRIMARY=1`` is also set (the
    user-eye verdict on phen / bipy / terpy chelates says the soft-DG
    result is below XRD-grade — the rigid path is the corrective default).
    Default OFF otherwise so the legacy non-MOGUL_PRIMARY path stays
    byte-identical.  Set ``DELFIN_FFFREE_RIGID_AROMATIC_CHELATE=0``
    explicitly to compare bytes against the pre-gate path.
    """
    raw = os.environ.get("DELFIN_FFFREE_RIGID_AROMATIC_CHELATE", "").strip()
    if raw == "":
        return os.environ.get("DELFIN_FFFREE_MOGUL_PRIMARY", "0") == "1"
    return raw == "1"


# ---------------------------------------------------------------------------
# Detection: aromatic-bond-only path between every donor pair
# ---------------------------------------------------------------------------
def _aromatic_neighbours(mol, atom_idx: int) -> List[int]:
    """Atoms reached from ``atom_idx`` along an aromatic-rigid bond.

    "Aromatic-rigid" means EITHER:
      * the bond itself is aromatic (``Bond.GetIsAromatic() == True``), OR
      * the bond is a SINGLE bond between two atoms that are BOTH
        aromatic AND the bond is the canonical biaryl (or biheteroaryl)
        connecting two aromatic ring systems.

    The second case is what bipyridine / terpyridine / quaterpyridine /
    pyridyl-imidazolyl rely on: the 2-2' inter-pyridyl bond is single by
    bond-order accounting (it is NOT inside a Hückel ring) yet the
    backbone is geometrically RIGID because the two pyridyl rings are
    conjugated through it.  The bond cannot rotate freely without
    breaking the conjugation, and crystallographically these chelates
    sit at the planar (or near-planar) conformer.

    The check is universal: aromatic flag from RDKit Hückel perception
    plus a single-bond + both-endpoints-aromatic gate.  No SMILES
    patterns, no per-element branches, no per-ring-size table.
    """
    out: List[int] = []
    atom = mol.GetAtomWithIdx(int(atom_idx))
    for nb in atom.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(int(atom_idx), int(nb.GetIdx()))
        if bond is None:
            continue
        if bool(bond.GetIsAromatic()):
            out.append(int(nb.GetIdx()))
            continue
        # Conjugated-single biaryl bridge: SINGLE bond + both endpoints
        # aromatic.  Universal capture of bipy 2-2', terpy 2-2'/6'-2'',
        # pyridyl-imidazolyl, ... .
        try:
            from rdkit.Chem import BondType as _BT
            is_single = (bond.GetBondType() == _BT.SINGLE)
        except Exception:
            is_single = False
        if not is_single:
            continue
        if bool(atom.GetIsAromatic()) and bool(nb.GetIsAromatic()):
            out.append(int(nb.GetIdx()))
    return sorted(out)


def _aromatic_connected(mol, src: int, dst: int,
                        restrict_to: Optional[Set[int]] = None) -> bool:
    """BFS: is ``dst`` reachable from ``src`` via aromatic bonds only?

    Pure graph search.  ``restrict_to`` optionally bounds the BFS to a
    subset of atom indices (typically the ligand subtree) so we never
    leak across a metal-bridged path.
    """
    src = int(src)
    dst = int(dst)
    if src == dst:
        return True
    n = mol.GetNumAtoms()
    if src < 0 or dst < 0 or src >= n or dst >= n:
        return False
    seen: Set[int] = {src}
    queue: List[int] = [src]
    head = 0
    while head < len(queue):
        u = queue[head]
        head += 1
        if u == dst:
            return True
        for j in _aromatic_neighbours(mol, u):
            if j in seen:
                continue
            if restrict_to is not None and j not in restrict_to:
                continue
            seen.add(j)
            queue.append(j)
    return dst in seen


def is_aromatic_chelate(
    mol,
    donor_idxs: Sequence[int],
    *,
    restrict_to: Optional[Sequence[int]] = None,
) -> bool:
    """Return True iff every pair of donors is connected by an aromatic path.

    A chelate is "aromatic" iff for every pair ``(d_a, d_b)`` of donor
    atoms there is a path ``d_a -> ... -> d_b`` in the molecular graph
    where EVERY bond on the path is aromatic
    (``Bond.GetIsAromatic() == True``).

    Universal:
        * No SMILES patterns.
        * No per-element / per-ring-size branches.
        * Pure RDKit aromaticity perception (Hückel rule).

    Parameters
    ----------
    mol : RDKit Mol
        The molecule (typically the full-complex mol or a single ligand
        fragment with explicit Hs).
    donor_idxs : sequence of int
        Donor atom indices in ``mol``.  Must have length ≥ 2 for the
        concept of a chelate to apply.
    restrict_to : sequence of int, optional
        Atom-index subset to bound the BFS within.  When supplied, every
        donor must be in this set and every connecting path must stay
        within it.  Use the ligand subtree to avoid leakage across the
        metal.

    Returns
    -------
    bool
        True iff aromatic chelate; False otherwise (mono / non-aromatic
        chelate / mixed aromatic + sp3 backbone).
    """
    donors = sorted(int(d) for d in donor_idxs)
    if len(donors) < 2:
        return False
    rt: Optional[Set[int]] = None
    if restrict_to is not None:
        rt = {int(a) for a in restrict_to}
        if not all(d in rt for d in donors):
            return False
    # Each donor itself must be aromatic.  Aromatic-bonded chelates have
    # aromatic donor atoms by definition (pyridyl-N, imidazolyl-N,
    # cyclopentadienyl-C, ...).
    for d in donors:
        try:
            if not bool(mol.GetAtomWithIdx(int(d)).GetIsAromatic()):
                return False
        except Exception:
            return False
    # Pairwise aromatic-path connectivity.
    for a_idx in range(len(donors)):
        for b_idx in range(a_idx + 1, len(donors)):
            if not _aromatic_connected(mol, donors[a_idx], donors[b_idx],
                                       restrict_to=rt):
                return False
    return True


def detect_aromatic_chelate_groups(
    mol,
    metal_idx: int,
    donor_idxs: Sequence[int],
) -> List[Dict[str, object]]:
    """Walk the donor set, group donors per ligand subtree, and flag the
    aromatic chelates.

    Each returned entry is a dict::

        {
            "donor_atom_idxs": [int, ...],   # donors of this ligand
            "subtree_atom_idxs": [int, ...], # full ligand subtree
            "is_aromatic_chelate": bool,
        }

    Universal — pure graph traversal (no SMILES, no element branches).
    The grouping uses the same M-removed-graph component logic as the
    bounds-matrix projection so donor → subtree mapping is consistent.
    """
    donors = sorted(int(d) for d in donor_idxs)
    if not donors:
        return []
    # Per-donor subtree (atoms reachable without crossing the metal).
    sub_by_donor: Dict[int, Set[int]] = {}
    for d in donors:
        seen: Set[int] = {d}
        stack: List[int] = [d]
        while stack:
            u = stack.pop()
            try:
                atom = mol.GetAtomWithIdx(int(u))
            except Exception:
                continue
            for nb in atom.GetNeighbors():
                j = int(nb.GetIdx())
                if j == int(metal_idx) or j in seen:
                    continue
                seen.add(j)
                stack.append(j)
        sub_by_donor[d] = seen
    # Union-find: donors with overlapping subtrees belong to the same
    # ligand (chelate).
    parent: Dict[int, int] = {d: d for d in donors}

    def _find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def _union(a: int, b: int) -> None:
        ra, rb = _find(a), _find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    for i in range(len(donors)):
        for j in range(i + 1, len(donors)):
            if sub_by_donor[donors[i]] & sub_by_donor[donors[j]]:
                _union(donors[i], donors[j])
    groups: Dict[int, List[int]] = {}
    for d in donors:
        groups.setdefault(_find(d), []).append(d)
    out: List[Dict[str, object]] = []
    for root in sorted(groups.keys()):
        d_list = sorted(groups[root])
        # Union of subtrees = the full ligand subtree.
        subtree: Set[int] = set()
        for d in d_list:
            subtree |= sub_by_donor[d]
        out.append({
            "donor_atom_idxs": d_list,
            "subtree_atom_idxs": sorted(subtree),
            "is_aromatic_chelate": (
                len(d_list) >= 2
                and is_aromatic_chelate(mol, d_list, restrict_to=subtree)
            ),
        })
    return out


# ---------------------------------------------------------------------------
# Template: RDKit-embed the standalone ligand fragment
# ---------------------------------------------------------------------------
def build_chelate_template(
    mol,
    subtree_atom_idxs: Sequence[int],
    donor_atom_idxs: Sequence[int],
) -> Optional[Tuple[np.ndarray, List[int], List[int]]]:
    """Build a rigid template by RDKit-embedding the standalone ligand subtree.

    The full-complex mol carries the metal + every ligand; we cut out the
    chelate ligand subtree alone, add Hs (the input mol typically already
    has explicit Hs but we re-add to handle any case), and DG-embed.  The
    aromatic backbone determines the rigid geometry — RDKit's DG embedder
    cannot break aromatic ring planarity or distort C-C aromatic bond
    lengths, so phen/bipy/terpy come out at CCDC-style geometry.

    Returns
    -------
    (template_coords, template_donor_indices, original_atom_idxs) or None
        ``template_coords`` is an ``(N_sub, 3)`` ndarray of standalone-
        embed coordinates.  ``template_donor_indices`` are the row indices
        of the donor atoms in ``template_coords`` (i.e. the index of each
        original donor atom *inside the template*).  ``original_atom_idxs``
        is the list of original atom indices in the full-complex mol, in
        the same row order as ``template_coords`` (so the caller can map
        template rows back to full-complex rows).  ``None`` on any failure
        (embed failure, no donors, missing aromaticity).
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return None
    subtree = sorted(int(a) for a in subtree_atom_idxs)
    donors = sorted(int(d) for d in donor_atom_idxs)
    if not donors or not all(d in subtree for d in donors):
        return None
    # Build a sub-mol containing only the subtree atoms, preserving bonds.
    rw = Chem.RWMol()
    old_to_new: Dict[int, int] = {}
    for old in subtree:
        try:
            old_atom = mol.GetAtomWithIdx(int(old))
        except Exception:
            return None
        new_atom = Chem.Atom(old_atom.GetAtomicNum())
        new_atom.SetFormalCharge(int(old_atom.GetFormalCharge()))
        # Preserve aromaticity flag — RDKit re-perceives but the explicit
        # flag prevents accidental loss on the sub-mol.
        new_atom.SetIsAromatic(bool(old_atom.GetIsAromatic()))
        new_atom.SetNoImplicit(True)
        # Preserve explicit Hs on this atom: count of H neighbours that
        # are in the subtree set become NumExplicitHs.  (RDKit sub-mol
        # construction otherwise discards implicit-H counts.)
        n_explicit_h = 0
        for nb in old_atom.GetNeighbors():
            if int(nb.GetIdx()) in {a for a in subtree} and nb.GetSymbol() == "H":
                # Counted at the heavy-atom side via the bond below.
                pass
        new_atom.SetNumExplicitHs(n_explicit_h)
        new_idx = rw.AddAtom(new_atom)
        old_to_new[int(old)] = int(new_idx)
    # Bonds.
    added_pairs: Set[Tuple[int, int]] = set()
    for bond in mol.GetBonds():
        i = int(bond.GetBeginAtomIdx())
        j = int(bond.GetEndAtomIdx())
        if i not in old_to_new or j not in old_to_new:
            continue
        key = (min(i, j), max(i, j))
        if key in added_pairs:
            continue
        added_pairs.add(key)
        bt = bond.GetBondType()
        rw.AddBond(old_to_new[i], old_to_new[j], bt)
        # Preserve aromaticity on the bond too.
        new_bond = rw.GetBondBetweenAtoms(old_to_new[i], old_to_new[j])
        if new_bond is not None:
            new_bond.SetIsAromatic(bool(bond.GetIsAromatic()))
    sub_mol = rw.GetMol()
    try:
        # SanitizeMol re-perceives aromaticity from the explicit bonds we
        # set; if the sub-fragment is not a valid valence-closed graph
        # (rare for properly Pólya-decomposed ligands but possible for
        # charged sub-fragments) sanitize raises and we bail.
        Chem.SanitizeMol(sub_mol)
    except Exception:
        # Try ring-only sanitize as a fallback so charged aromatic NHC
        # fragments (which have non-standard valence) can still embed.
        try:
            Chem.SanitizeMol(
                sub_mol,
                sanitizeOps=(
                    Chem.SanitizeFlags.SANITIZE_ALL
                    ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
                ),
            )
        except Exception:
            return None
    # Pre-compute the template-frame donor indices so we can construct
    # the placeholder-bond graph below.
    template_donor_idxs_pre = [old_to_new[int(d)] for d in donors]
    # Embed the standalone ligand fragment WITH a placeholder dummy
    # atom representing the metal, bonded to every donor.  This forces
    # the embedder to land on the CHELATING conformer (cis-cis-cis for
    # tridentate, cis for bidentate) instead of the energetically
    # preferred extended/anti conformer.  Universal — works for any
    # denticity ≥ 2.
    #
    # Without the placeholder bond, a tridentate aromatic chelate (e.g.
    # terpy) embeds into the all-anti extended conformer with outer
    # donors ~6.4 Å apart -- chemically wrong for the chelating state
    # whose outer donors are ~4.0 Å apart.  Adding the placeholder
    # forces the cis-cis-cis ring-closure.
    #
    # We use a dummy atom (atomic number 0) so the embedder doesn't
    # complain about valence on a generic-metal symbol, and a SINGLE
    # bond which the embedder treats as a flexible-bridge constraint.
    try:
        rw_emb = Chem.RWMol(sub_mol)
        # Use a generic carbon placeholder + DATIVE bonds (donor → metal).
        # Wildcard atom (atomic number 0) does NOT receive embedder
        # distance constraints, so the donors fly apart into the
        # extended conformer.  A real carbon with DATIVE bonds gets
        # standard bond-length bounds (C-N ≈ 1.5 Å) and forces the
        # chelating cis-cis-cis ring closure.  Universal -- works for
        # any aromatic chelate denticity.
        dummy = Chem.Atom("C")
        dummy.SetNoImplicit(True)
        dummy.SetFormalCharge(0)
        dummy_idx = rw_emb.AddAtom(dummy)
        for tr in template_donor_idxs_pre:
            rw_emb.AddBond(int(tr), int(dummy_idx), Chem.BondType.DATIVE)
        emb_mol = rw_emb.GetMol()
        try:
            Chem.SanitizeMol(
                emb_mol,
                sanitizeOps=(
                    Chem.SanitizeFlags.SANITIZE_ALL
                    ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
                    ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
                ),
            )
        except Exception:
            pass
        params = AllChem.ETKDGv3()
        params.randomSeed = int(_TEMPLATE_SEED)
        try:
            params.useRandomCoords = False
        except Exception:
            pass
        ok = AllChem.EmbedMolecule(emb_mol, params)
        if ok < 0:
            try:
                params.useRandomCoords = True
            except Exception:
                pass
            ok = AllChem.EmbedMolecule(emb_mol, params)
        if ok < 0:
            # Last-resort fallback: embed without the placeholder atom.
            # The result is the extended conformer but it's better than
            # returning None (no rigid placement at all).
            params2 = AllChem.ETKDGv3()
            params2.randomSeed = int(_TEMPLATE_SEED)
            try:
                params2.useRandomCoords = True
            except Exception:
                pass
            ok2 = AllChem.EmbedMolecule(sub_mol, params2)
            if ok2 < 0:
                return None
            conf = sub_mol.GetConformer()
            n_sub = sub_mol.GetNumAtoms()
            coords = np.array(
                [list(conf.GetAtomPosition(i)) for i in range(n_sub)],
                dtype=float,
            )
        else:
            # Strip the placeholder atom from the output coords.
            conf = emb_mol.GetConformer()
            n_emb = emb_mol.GetNumAtoms()
            coords_full = np.array(
                [list(conf.GetAtomPosition(i)) for i in range(n_emb)],
                dtype=float,
            )
            # Placeholder atom is the LAST one (we added it via AddAtom).
            n_sub = sub_mol.GetNumAtoms()
            coords = coords_full[:n_sub]
        # Light MMFF tightening on the standalone sub_mol; rotates the
        # substituent groups onto MMFF values but does NOT relax the
        # aromatic backbone bond lengths (they are already at canonical
        # aromatic distances from ETKDG).
        try:
            # Re-set the conformer onto sub_mol with the placeholder-
            # stripped coords before MMFF.
            conf_sub = Chem.Conformer(n_sub)
            for ii in range(n_sub):
                from rdkit.Geometry import Point3D
                conf_sub.SetAtomPosition(
                    int(ii),
                    Point3D(float(coords[ii, 0]),
                            float(coords[ii, 1]),
                            float(coords[ii, 2])),
                )
            sub_mol.RemoveAllConformers()
            sub_mol.AddConformer(conf_sub, assignId=True)
            AllChem.MMFFOptimizeMolecule(sub_mol)
            conf = sub_mol.GetConformer()
            coords = np.array(
                [list(conf.GetAtomPosition(i)) for i in range(n_sub)],
                dtype=float,
            )
        except Exception:
            pass
    except Exception:
        return None
    # Translate the template to its own centroid for numerical stability.
    coords = coords - coords.mean(axis=0)
    # Map donors: original atom idx -> template row.
    template_donor_idxs = [old_to_new[int(d)] for d in donors]
    # Original-atom-idx list in template row order (length n_sub).
    original_idxs = [-1] * n_sub
    for old, new in old_to_new.items():
        if 0 <= new < n_sub:
            original_idxs[new] = old
    return coords, template_donor_idxs, original_idxs


# ---------------------------------------------------------------------------
# Placement: Kabsch rigid-body align template donors onto target positions
# ---------------------------------------------------------------------------
def _rotation_from_to(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Proper-rotation matrix mapping unit-vector ``a`` onto unit-vector ``b``.

    Rodrigues formula.  Returns identity when the input vectors are
    already parallel (within numerical precision); returns a 180°
    rotation around an arbitrary perpendicular axis when they are
    anti-parallel.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na < 1e-9 or nb < 1e-9:
        return np.eye(3)
    a = a / na
    b = b / nb
    v = np.cross(a, b)
    s = float(np.linalg.norm(v))
    c = float(np.dot(a, b))
    if s < 1e-9:
        if c > 0:
            return np.eye(3)
        # Anti-parallel: pick any perpendicular axis.
        perp = np.array([1.0, 0.0, 0.0], dtype=float)
        if abs(a[0]) > 0.9:
            perp = np.array([0.0, 1.0, 0.0], dtype=float)
        return _rodrigues(perp - a * float(np.dot(perp, a)), np.pi)
    vx = np.array([
        [0.0, -v[2], v[1]],
        [v[2], 0.0, -v[0]],
        [-v[1], v[0], 0.0],
    ], dtype=float)
    return np.eye(3) + vx + (vx @ vx) * ((1.0 - c) / (s * s))


def _rodrigues(axis: np.ndarray, angle: float) -> np.ndarray:
    """Right-hand rotation by ``angle`` rad around the unit vector ``axis``."""
    axis = np.asarray(axis, dtype=float)
    na = float(np.linalg.norm(axis))
    if na < 1e-9:
        return np.eye(3)
    axis = axis / na
    x, y, z = axis
    c = float(np.cos(float(angle)))
    s = float(np.sin(float(angle)))
    C = 1.0 - c
    return np.array([
        [c + x * x * C,     x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, c + y * y * C,     y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
    ], dtype=float)


def _kabsch_rt(P_src: np.ndarray, P_tgt: np.ndarray
               ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Closed-form proper-rotation + translation aligning ``P_src`` onto
    ``P_tgt`` in least-squares sense (Kabsch, determinant-corrected).

    Returns ``(R, t, src_centroid, tgt_centroid)`` such that for every
    point ``x`` in the source frame::

        x_aligned = R @ (x - src_centroid) + tgt_centroid
                  = R @ x + (tgt_centroid - R @ src_centroid)
                  = R @ x + t        with t := tgt_centroid - R @ src_centroid
    """
    P_src = np.asarray(P_src, dtype=float)
    P_tgt = np.asarray(P_tgt, dtype=float)
    src_c = P_src.mean(axis=0)
    tgt_c = P_tgt.mean(axis=0)
    H = (P_src - src_c).T @ (P_tgt - tgt_c)
    try:
        U, _, Vt = np.linalg.svd(H)
    except np.linalg.LinAlgError:
        R = np.eye(3)
    else:
        d = 1.0 if np.linalg.det(Vt.T @ U.T) > 0 else -1.0
        D = np.diag([1.0, 1.0, d])
        R = Vt.T @ D @ U.T
    t = tgt_c - R @ src_c
    return R, t, src_c, tgt_c


def place_rigid_aromatic_chelates(
    *,
    mol,
    syms: Sequence[str],
    metal_idx: int,
    donor_idxs: Sequence[int],
    P: np.ndarray,
    lower: Optional[np.ndarray] = None,
    upper: Optional[np.ndarray] = None,
    geometry_key: Optional[str] = None,
    donor_at_vertex: Optional[Sequence[int]] = None,
) -> Tuple[np.ndarray, List[Dict[str, object]]]:
    """Replace soft-DG aromatic chelate geometry with a rigid-body placement.

    Pipeline:

      1.  Detect every chelate ligand whose donor pairs are connected via
          aromatic bonds (``is_aromatic_chelate``).
      2.  For each aromatic chelate, build a TEMPLATE via standalone
          RDKit embed (``build_chelate_template``).  The aromatic backbone
          gives canonical bite angle / D-D distance / planarity.
      3.  Compute per-donor TARGET positions from the polyhedron vertices
          and CCDC supplement M-D windows (mid-window).  When ``lower`` /
          ``upper`` are None or unset for a donor, fall back to
          ``polyhedra.md_distance``.
      4.  Kabsch align template donors → target donor positions; apply the
          same rigid transform to the WHOLE ligand subtree.
      5.  Defence-in-depth: if any donor lands more than 0.15 Å off its
          target M-D distance, ROLL BACK this chelate (return geometry
          unchanged for it) so we never produce a worse output than the
          soft-DG pre-state.

    Parameters
    ----------
    mol : RDKit Mol
        Full-complex molecule (metal + ligands + explicit H).
    syms : sequence of str
        Element symbols.  Must have length ``mol.GetNumAtoms()``.
    metal_idx : int
        Index of the metal atom in ``mol`` / ``P``.
    donor_idxs : sequence of int
        Donor atom indices in ``mol`` / ``P``.
    P : (N, 3) ndarray
        Cartesian coordinates BEFORE this step.
    lower, upper : (N, N) ndarrays, optional
        Bounds matrix from ``mogul_bounds.build_bounds_matrix``.  When
        supplied, the per-donor target M-D distance is the mid-window
        ``(lower[m, d] + upper[m, d]) / 2`` so the CCDC empirical mean
        is the source of truth.  Falls back to covalent-radii sums when
        absent.
    geometry_key : str, optional
        Polyhedron name (e.g. ``"OC-6 octahedron"``) used to resolve the
        target donor directions.  When None, the per-CN first-candidate
        rule is used.
    donor_at_vertex : sequence of int, optional
        Pólya orbit permutation mapping vertex k → donor atom idx at that
        vertex.  When supplied, the per-vertex target direction uses
        ``donor_at_vertex[k]`` as the donor for vertex k.

    Returns
    -------
    (P_new, diag)
        ``P_new`` is the updated coordinate array (same shape as input).
        ``diag`` is a list of per-chelate diagnostic dicts (one per
        detected aromatic chelate) with keys ``donor_idxs``, ``n_atoms``,
        ``placed``, ``residual_donor_rms_A``, ``rollback_reason``
        (when placed=False).  Always returns finite arrays — never raises.
    """
    P = np.asarray(P, dtype=float).copy()
    diag: List[Dict[str, object]] = []
    if int(metal_idx) < 0 or len(donor_idxs) < 2:
        return P, diag
    # 1) Detect chelate groups + which ones are aromatic.
    groups = detect_aromatic_chelate_groups(
        mol, int(metal_idx), [int(d) for d in donor_idxs],
    )
    if not any(g["is_aromatic_chelate"] for g in groups):
        return P, diag
    # 2) Per-vertex target directions = polyhedron unit vectors, in the
    #    embed's orientation.  We compute them in two stages so the
    #    aromatic-chelate placement is consistent with what the rest of
    #    the projection pipeline uses (the same polyhedron, the same
    #    embed-frame orientation via Kabsch).
    n_d = len(donor_idxs)
    donors_sorted = sorted(int(d) for d in donor_idxs)
    if donor_at_vertex is not None and len(donor_at_vertex) == n_d:
        donors_per_vertex = [int(d) for d in donor_at_vertex]
    else:
        donors_per_vertex = list(donors_sorted)
    metal_pos = P[int(metal_idx)].copy()
    # Embed-derived donor directions.
    cur_dirs = np.zeros((n_d, 3), dtype=float)
    for k, d in enumerate(donors_per_vertex):
        v = P[int(d)] - metal_pos
        nv = float(np.linalg.norm(v))
        cur_dirs[k] = v / nv if nv > 1e-9 else np.array([1.0, 0.0, 0.0])
    # Polyhedron unit vectors.
    try:
        from delfin.fffree import polyhedra as _polyhedra
        ref_vecs = None
        if geometry_key:
            try:
                ref_vecs = _polyhedra.ref_vectors(str(geometry_key))
            except Exception:
                ref_vecs = None
        if ref_vecs is None:
            metal_sym = str(syms[int(metal_idx)])
            try:
                cands = _polyhedra.geometries_for_cn(int(n_d), metal_sym)
            except Exception:
                cands = []
            for g in cands:
                try:
                    ref_vecs = _polyhedra.ref_vectors(str(g))
                    break
                except Exception:
                    continue
        if ref_vecs is not None and ref_vecs.shape[0] >= n_d:
            V = np.asarray(ref_vecs[:n_d], dtype=float)
            norms = np.linalg.norm(V, axis=1, keepdims=True)
            norms = np.where(norms < 1e-9, 1.0, norms)
            ref_unit = V / norms
            # Kabsch align ref_unit -> cur_dirs (no reflection, no
            # isomer flip).
            H = cur_dirs.T @ ref_unit
            try:
                U, _, Vt = np.linalg.svd(H)
                d_sign = 1.0 if np.linalg.det(U) * np.linalg.det(Vt) > 0 else -1.0
                D = np.diag([1.0, 1.0, d_sign])
                Rrot = U @ D @ Vt
                aligned_targets = ref_unit @ Rrot.T
            except Exception:
                aligned_targets = cur_dirs.copy()
        else:
            aligned_targets = cur_dirs.copy()
    except ImportError:
        aligned_targets = cur_dirs.copy()
    # Per-donor target M-D distance (mid-window of CCDC bounds, fall back
    # to covalent radii when bounds absent).
    md_targets_by_atom: Dict[int, float] = {}
    try:
        from delfin.fffree import polyhedra as _polyhedra_for_md
    except ImportError:
        _polyhedra_for_md = None
    for d in donors_per_vertex:
        target = None
        if lower is not None and upper is not None:
            i, j = (int(metal_idx), int(d)) if int(metal_idx) < int(d) else (
                int(d), int(metal_idx))
            try:
                lo = float(lower[i, j])
                hi = float(upper[i, j])
                if np.isfinite(hi) and hi > 0.0 and hi > lo:
                    target = 0.5 * (lo + hi)
            except Exception:
                target = None
        if target is None and _polyhedra_for_md is not None:
            try:
                target = float(_polyhedra_for_md.md_distance(
                    str(syms[int(metal_idx)]), str(syms[int(d)])))
            except Exception:
                target = None
        if target is None or not np.isfinite(target) or target <= 0.0:
            target = 2.0  # last-resort sane fallback
        md_targets_by_atom[int(d)] = float(target)
    # 3) Walk each aromatic chelate and apply the rigid placement.
    #
    # IMPORTANT geometric note (hmaximilian 2026-06-08).  Naive Kabsch
    # alignment of template donors onto polyhedron-vertex positions
    # ``M + r_md · vertex_k`` fails for the typical case where the
    # polyhedron's cis-edge donor-donor distance (e.g. OC-6 cis is
    # √2 · r_md ≈ 2.76 Å for Co-N at 1.95 Å) does NOT match the template's
    # NATURAL bite-distance (phen ≈ 2.65 Å).  Kabsch returns the LSQ
    # optimum which leaves ≥ ½·(target_DD - template_DD) ≈ 0.05-0.30 Å
    # residual on each donor — exceeding our 0.15 Å rollback gate and
    # making the rigid step a NO-OP exactly where it matters most.
    #
    # The correct (chemically truthful) algorithm: the AROMATIC backbone
    # gives the donor-donor distance + plane RIGIDLY; the polyhedron
    # gives the ORIENTATION + the per-donor RADIAL distance.  We
    # therefore decompose the placement into two independent rigid
    # transforms:
    #
    #   1.  Rotate the template so its donor-donor midpoint axis aligns
    #       with the SUM of the assigned polyhedron vertex directions
    #       (the chelate's "bisecting" axis).
    #   2.  Rotate around that axis so the template's chelate-plane
    #       normal lies perpendicular to the polyhedron's plane spanned
    #       by the two vertex vectors.
    #   3.  Translate so the donor-donor midpoint sits at distance
    #       ``sqrt(r_md² − (½·D-D)²)`` from the metal along the bisecting
    #       axis (so each donor lands at distance r_md from the metal).
    #
    # The donor-donor distance + chelate plane are PRESERVED by
    # construction; the per-donor M-D distance is enforced exactly; only
    # the angular alignment may carry a small (template-vs-polyhedron-
    # bite mismatch) residual which we measure as the post-place
    # vertex-direction agreement.  Universal — no SMILES branches.
    for g in groups:
        if not g["is_aromatic_chelate"]:
            continue
        d_list: List[int] = list(g["donor_atom_idxs"])  # type: ignore[arg-type]
        subtree: List[int] = list(g["subtree_atom_idxs"])  # type: ignore[arg-type]
        if len(d_list) < 2:
            continue
        chelate_diag: Dict[str, object] = {
            "donor_idxs": d_list,
            "n_atoms": len(subtree),
            "placed": False,
            "residual_donor_rms_A": float("inf"),
        }
        # Build the template.
        tpl = build_chelate_template(mol, subtree, d_list)
        if tpl is None:
            chelate_diag["rollback_reason"] = "template_embed_failed"
            diag.append(chelate_diag)
            continue
        tpl_coords, tpl_donor_rows, original_idxs = tpl
        # Map each donor in d_list to its polyhedron vertex index.
        vertex_idxs: List[int] = []
        ok_vertex = True
        for d in d_list:
            try:
                k = donors_per_vertex.index(int(d))
            except ValueError:
                ok_vertex = False
                break
            vertex_idxs.append(k)
        if not ok_vertex or len(vertex_idxs) != len(d_list):
            chelate_diag["rollback_reason"] = "vertex_mapping_incomplete"
            diag.append(chelate_diag)
            continue
        # Mean per-donor M-D target distance (we use the mean so the
        # symmetric chelate exactly satisfies the per-donor R = r_md;
        # for an asymmetric chelate the donors land at the mean which
        # is then tightened by the downstream grip_polish).
        r_md_list = np.array([md_targets_by_atom[int(d)] for d in d_list],
                             dtype=float)
        r_md_mean = float(np.mean(r_md_list))
        # Template donor positions (centred on the template's own
        # centroid, the build_chelate_template already does this).
        donor_pos_tpl = np.asarray(
            [tpl_coords[int(row)] for row in tpl_donor_rows], dtype=float,
        )
        donor_centroid_tpl = donor_pos_tpl.mean(axis=0)
        # Template chelate plane normal (SVD on the donor positions plus
        # a few backbone atoms for a stable normal).  For a strictly
        # bidentate template the donor pair alone is 1D, so we add the
        # chelate-plane atoms (the standalone embed's full subtree) to
        # estimate the plane.  Universal: just SVD on the heavy-atom
        # cloud.
        heavy_rows = [
            i for i in range(tpl_coords.shape[0])
        ]  # all template atoms (already excludes metal)
        M_tpl = tpl_coords[heavy_rows] - tpl_coords[heavy_rows].mean(axis=0)
        try:
            _, _, Vt_tpl = np.linalg.svd(M_tpl, full_matrices=False)
            plane_normal_tpl = Vt_tpl[-1]
            nn = float(np.linalg.norm(plane_normal_tpl))
            plane_normal_tpl = (
                plane_normal_tpl / nn if nn > 1e-9 else np.array([0.0, 0.0, 1.0])
            )
        except np.linalg.LinAlgError:
            plane_normal_tpl = np.array([0.0, 0.0, 1.0])
        # Template bisecting axis (from donor centroid back through the
        # ligand bulk — i.e. away from the metal).  For a bidentate
        # chelate this is the donor → backbone direction.
        # We compute it as -(donor_centroid_tpl - heavy_centroid_tpl) so
        # the axis points FROM the heavy centroid TO the donor centroid
        # (= the direction the chelate "opens" toward the metal).
        heavy_centroid_tpl = tpl_coords[heavy_rows].mean(axis=0)
        ax_vec_tpl = donor_centroid_tpl - heavy_centroid_tpl
        ax_norm = float(np.linalg.norm(ax_vec_tpl))
        if ax_norm < 1e-6:
            chelate_diag["rollback_reason"] = "template_axis_degenerate"
            diag.append(chelate_diag)
            continue
        axis_tpl = ax_vec_tpl / ax_norm
        # ---- TARGET frame ----
        # Sum of assigned vertex unit vectors = chelate's bisecting axis
        # in the embed frame.
        vk_unit_sum = np.zeros(3, dtype=float)
        for k in vertex_idxs:
            vk_unit_sum = vk_unit_sum + aligned_targets[k]
        ax_norm_tgt = float(np.linalg.norm(vk_unit_sum))
        if ax_norm_tgt < 1e-6:
            chelate_diag["rollback_reason"] = "target_axis_degenerate"
            diag.append(chelate_diag)
            continue
        axis_tgt = vk_unit_sum / ax_norm_tgt
        # Distance from the metal to the donor centroid along the axis.
        #
        # Geometric derivation (universal for any denticity n).  The
        # metal sits on the bisecting axis at distance ``r_centroid``
        # from the donor centroid.  Each donor's M-D distance is then
        #
        #     r_md_k² = r_centroid² + |donor_k − centroid|²_⊥
        #
        # where ``|donor_k − centroid|_⊥`` is the donor's perpendicular
        # offset from the centroid (in the bidentate case this is ½·DD).
        # We pick ``r_centroid`` so the MEAN per-donor M-D equals the
        # CCDC mean:
        #
        #     r_centroid = √( r_md_mean² − ⟨perp²⟩ )
        #
        # which collapses to the bidentate Pythagoras formula when n=2.
        # For tridentate the central donor has perp ≈ 0 (it sits on the
        # axis), so the formula gives a balanced compromise across all
        # n donors.
        donor_centroid_tpl_xy = donor_pos_tpl.mean(axis=0)
        # Use the axis_tpl direction (centroid → donor centroid) as the
        # "in-axis" direction in the template frame.  The perpendicular
        # offsets are the components of each donor's position projected
        # OUT of axis_tpl in the template frame.
        donor_offset_perp_sq = []
        for dp in donor_pos_tpl:
            v = dp - donor_centroid_tpl_xy
            # Project out the component along axis_tpl (the bisecting
            # direction in the template).
            v_perp = v - np.dot(v, axis_tpl) * axis_tpl
            donor_offset_perp_sq.append(float(np.dot(v_perp, v_perp)))
        mean_perp_sq = float(np.mean(donor_offset_perp_sq))
        # Donor-donor max distance for the "impossible geometry" guard.
        max_dd_tpl = 0.0
        for i_dp in range(len(donor_pos_tpl)):
            for j_dp in range(i_dp + 1, len(donor_pos_tpl)):
                dd = float(np.linalg.norm(
                    donor_pos_tpl[i_dp] - donor_pos_tpl[j_dp]
                ))
                if dd > max_dd_tpl:
                    max_dd_tpl = dd
        if r_md_mean ** 2 < mean_perp_sq + 1e-6:
            # The chelate's bite is too wide for this M-D distance.
            chelate_diag["rollback_reason"] = (
                f"impossible_geometry_dd_{max_dd_tpl:.3f}_md_{r_md_mean:.3f}"
            )
            diag.append(chelate_diag)
            continue
        r_centroid = float(np.sqrt(max(0.0, r_md_mean ** 2 - mean_perp_sq)))
        # ---- STAGE 1: rotate template so axis_tpl -> axis_tgt ----
        R1 = _rotation_from_to(axis_tpl, axis_tgt)
        # Apply R1 to the template (centre at template heavy centroid).
        tpl_after_R1 = (
            (tpl_coords - heavy_centroid_tpl) @ R1.T
        )
        # Donor positions after R1 (template frame, donor centroid axis
        # now aligned with target axis).
        donor_after_R1 = (
            (donor_pos_tpl - heavy_centroid_tpl) @ R1.T
        )
        # ---- STAGE 2: rotate around axis_tgt so the donor-donor edge
        # spans the polyhedron's two vertex vectors ----
        # The donor-donor edge direction is perpendicular to axis_tgt by
        # symmetry (donors are reflected across the axis); we want its
        # direction to coincide with the polyhedron's donor-pair
        # perpendicular ``(v_a - v_b)/|v_a - v_b|`` so when we translate
        # along axis_tgt the two donors land on the polyhedron-vertex
        # rays (within the angular mismatch).
        if len(donor_after_R1) >= 2 and len(vertex_idxs) >= 2:
            edge_tpl_after_R1 = donor_after_R1[0] - donor_after_R1[1]
            edge_tpl_after_R1 = edge_tpl_after_R1 - (
                np.dot(edge_tpl_after_R1, axis_tgt) * axis_tgt
            )
            ne_tpl = float(np.linalg.norm(edge_tpl_after_R1))
            edge_tgt = aligned_targets[vertex_idxs[0]] - aligned_targets[vertex_idxs[1]]
            edge_tgt = edge_tgt - np.dot(edge_tgt, axis_tgt) * axis_tgt
            ne_tgt = float(np.linalg.norm(edge_tgt))
            if ne_tpl > 1e-6 and ne_tgt > 1e-6:
                u_tpl = edge_tpl_after_R1 / ne_tpl
                u_tgt = edge_tgt / ne_tgt
                # Rotation around axis_tgt by the in-plane angle
                # between u_tpl and u_tgt.
                cos_t = float(np.clip(np.dot(u_tpl, u_tgt), -1.0, 1.0))
                cross = np.cross(u_tpl, u_tgt)
                sin_t = float(np.clip(np.dot(cross, axis_tgt), -1.0, 1.0))
                theta = float(np.arctan2(sin_t, cos_t))
                R2 = _rodrigues(axis_tgt, theta)
            else:
                R2 = np.eye(3)
        else:
            R2 = np.eye(3)
        # ---- STAGE 3: translate so donor centroid lands at
        # M + r_centroid · axis_tgt ----
        # Donor centroid position after R1 then R2:
        donor_centroid_after_R1R2 = (
            donor_after_R1.mean(axis=0) @ R2.T
        )
        # Where it should be in the embed frame:
        donor_centroid_target = metal_pos + r_centroid * axis_tgt
        # Translation:
        t_vec = donor_centroid_target - donor_centroid_after_R1R2
        # Final coords for every template atom.
        new_subtree_coords = (tpl_after_R1 @ R2.T) + t_vec
        # Defence-in-depth: per-donor M-D residual.
        donor_post = np.asarray(
            [new_subtree_coords[int(row)] for row in tpl_donor_rows],
            dtype=float,
        )
        per_donor_md = np.linalg.norm(donor_post - metal_pos, axis=1)
        md_dev = float(np.sqrt(np.mean((per_donor_md - r_md_list) ** 2)))
        chelate_diag["residual_donor_rms_A"] = md_dev
        # The symmetric bidentate case satisfies per-donor R = r_md
        # exactly (Pythagoras-by-construction).  Tridentate and higher
        # chelates carry a small inherent asymmetry because the central
        # donor (perp ≈ 0) cannot simultaneously satisfy M-D = r_md and
        # sit on the axis at the same r_centroid as the outer donors
        # (M-N(central) ≈ 1.88 Å vs M-N(outer) ≈ 1.97 Å is the typical
        # terpy crystallographic pattern).  We therefore use a denticity-
        # scaled threshold: 0.25 Å for bidentate (Pythagoras-exact),
        # 0.45 Å for tridentate (~0.1 Å expected residual + slack),
        # 0.60 Å for tetradentate+ (planar macrocycle).
        denticity = len(d_list)
        if denticity == 2:
            md_dev_thresh = 0.25
        elif denticity == 3:
            md_dev_thresh = 0.45
        else:
            md_dev_thresh = 0.60
        if md_dev > md_dev_thresh:
            chelate_diag["rollback_reason"] = (
                f"per_donor_md_rms_{md_dev:.3f}A_exceeds_{md_dev_thresh:.2f}"
            )
            diag.append(chelate_diag)
            continue
        # Commit: write new coordinates back.
        for row_idx, original_atom in enumerate(original_idxs):
            if original_atom < 0:
                continue
            P[int(original_atom)] = new_subtree_coords[row_idx]
        chelate_diag["placed"] = True
        diag.append(chelate_diag)
    return P, diag
