"""delfin.fffree.converter_backend — adapter wiring the metal-FF-free foundation
into the converter's smiles_to_xyz_isomers contract.

_fffree_isomers(smiles) -> [(xyz_string, label), ...]  or  None (-> legacy fallback).

v1 handles mononuclear Werner complexes with explicit metal-donor bonds, all
monodentate donors, CN 4/5/6 (decompose.py).  Anything else returns None so the
caller falls through to the existing pipeline.  Deterministic.
"""
from __future__ import annotations
import os
from collections import Counter
from typing import List, Optional, Tuple
import numpy as np
from rdkit import Chem

from delfin.fffree import decompose as DEC
from delfin.fffree import polya_isomer_count as PIC
from delfin.fffree import assemble_complex as AC
from delfin.fffree import polyhedra as PLY
from delfin.fffree import ligand_relax as LR
from delfin.fffree import backbone_reembed as _BR
import delfin._bond_decollapse as _bd


def _multibond_enabled() -> bool:
    """#279/#281: bond-order-aware collapse self-gate.  When enabled, genuine SHORT
    multiple/aromatic bonds (C≡O carbonyl, C≡N nitrile, N=N/N≡N azo, aromatic/imine
    C=N) are LENGTH-GATED-exempted from the `collapsed_bond` check so an EXCELLENT
    FF-free build is no longer false-flagged as collapsed and dropped to the distorted
    legacy-UFF fallback (root cause of the 9% native rate on the hard pool).  Default
    OFF => no exempt_pairs passed => exact historic gate (byte-identical)."""
    return os.environ.get("DELFIN_FFFREE_MULTIBOND_EXEMPT", "0") == "1"


def _local_multibond_ideals(mol):
    """LOCAL heavy-heavy multiple/aromatic bonds of a ligand ``mol`` -> mapping
    {(a1, a2): ideal_multibond_length} in the mol's OWN atom-index space (heavy-atom
    indices are preserved under AddHs, so they map 1:1 onto the assembled block).
    The ideal length is the bond-order-appropriate Pyykkö covalent-radius sum (#305
    PYYKKO_DOUBLE / PYYKKO_TRIPLE), so the self-gate can verify a flagged-short bond
    is genuinely a multiple bond (within tolerance) rather than a true collapse.
    Graph-only, deterministic; never raises."""
    from rdkit.Chem import BondType
    out = {}
    try:
        for b in mol.GetBonds():
            bt = b.GetBondType()
            if b.GetIsAromatic() or bt == BondType.AROMATIC:
                order = 1.5
            elif bt == BondType.DOUBLE:
                order = 2.0
            elif bt == BondType.TRIPLE:
                order = 3.0
            else:
                continue                              # single bond -> never exempt
            aa, ab = b.GetBeginAtom(), b.GetEndAtom()
            if aa.GetAtomicNum() <= 1 or ab.GetAtomicNum() <= 1:
                continue                              # heavy-heavy only
            sa, sb = aa.GetSymbol(), ab.GetSymbol()
            ideal = PLY._pyykko_radius(sa, order) + PLY._pyykko_radius(sb, order)
            i, j = aa.GetIdx(), ab.GetIdx()
            out[(min(i, j), max(i, j))] = float(ideal)
    except Exception:
        return {}
    return out


def _exempt_from_blocks(block_mols_offsets):
    """Assemble the GLOBAL length-gated exempt-pair mapping for a built complex from
    a list of (ligand_mol, block_offset) — the block_offset is the global index where
    the ligand's AddHs block begins (metal at 0).  Local heavy-atom indices map to
    global as offset + local (AddHs preserves heavy-atom order).  Returns
    {(gi, gj): ideal_multibond_length}; empty when the flag is off."""
    if not _multibond_enabled():
        return {}
    ex = {}
    for mol, off in block_mols_offsets:
        for (i, j), ideal in _local_multibond_ideals(mol).items():
            gi, gj = off + int(i), off + int(j)
            ex[(min(gi, gj), max(gi, gj))] = ideal
    return ex


def _heteroleptic_block_offsets(vertex_specs):
    """Deterministic per-ligand block offsets for the monodentate heteroleptic build
    order (assemble_heteroleptic_from_mols / _ensemble / _enumerate_geometry): metal at
    0, then one AddHs block per vertex_spec in order.  Returns
    [(frag_mol, offset), ...] matching the placement layout exactly."""
    out = []
    pos = 1
    for frag, _di in vertex_specs:
        out.append((frag, pos))
        pos += Chem.AddHs(frag).GetNumAtoms()
    return out


def _config_block_offsets(config, ligands):
    """Deterministic per-ligand block offsets for the chelate-config build order
    (assemble_from_config places ligands in first-appearance order of the config
    dict, one AddHs block per ligand instance).  Returns [(lig_mol, offset), ...]
    matching the placement layout exactly (same dict-iteration order as the builder)."""
    out = []
    pos = 1
    seen = set()
    for _v, (li, _arm) in config.items():
        if li in seen:
            continue
        seen.add(li)
        mol = ligands[li]["mol"]
        out.append((mol, pos))
        pos += Chem.AddHs(mol).GetNumAtoms()
    return out


def _maybe_relax(syms, P):
    """#38: env-gated COD-loss torsional/rigid-body ligand relaxer (default OFF).
    Relieves van-der-Waals clashes by rotating distal sub-trees about rotatable bonds —
    coordination (metal + donors within the coord_geom detection sphere) frozen, rigid
    fragments preserved, multi-axis never-worse firewall.  Validated on smoke: net +11,
    0 severe (hanom -22%, inter-ligand -29%, h-clash -27%, coord_geom unchanged).
    Enable via DELFIN_FFFREE_LIGAND_RELAX=1."""
    if os.environ.get("DELFIN_FFFREE_LIGAND_RELAX", "0") != "1":
        return syms, P
    try:
        P2 = np.asarray(P, dtype=float)
        mi = [i for i, s in enumerate(syms) if _bd._is_metal(s)]
        fixed = set(mi)
        for m in mi:
            for j in range(len(syms)):
                if j != m and syms[j] != "H" and float(np.linalg.norm(P2[j] - P2[m])) \
                        < 1.45 * _bd._ideal_bond(syms[m], syms[j]):
                    fixed.add(j)
        return list(syms), LR.relax(list(syms), P2, fixed)
    except Exception:
        return syms, P

_GEOM_TO_POLYA = {
    "L-2 linear": "linear",                        # iter-32f (DELFIN_FFFREE_CN_EXTEND): CN2
    "SP-3 trigonal planar": "trigonal_planar",     # iter-32c (User 2026-05-28 ADUMOD): CN3
    "T-3 T-shape": "tshape",
    "OC-6 octahedron": "octahedron",
    "SP-4 square planar": "square_planar",
    "T-4 tetrahedron": "tetrahedron",
    "TBP-5 trigonal bipyramid": "trigonal_bipyramid",
    "SPY-5 square pyramid": "square_pyramid",
    "TPR-6 trigonal prism": "trigonal_prism",      # iter-31 (User 2026-05-28): CN6 dual
    "PB-7 pentagonal bipyramid": "pentagonal_bipyramid",
    "SQAP-8 square antiprism": "square_antiprism",
    "TTP-9 tricapped trigonal prism": "tricapped_trigonal_prism",
}


# antipodal vertex pairs per geometry (polya vertex ordering) — for universal
# cis/trans/fac/mer classification from the coloring.  Iter-32d (User 2026-05-28
# GUVZIH "fac koord fehlt"): extended to CN5 (TBP/SPY) + CN3 (T-shape).
# Octahedron: opposite pairs along x/y/z axes.
# Square-planar: opposite pairs across the square.
# TBP-5: axials 0↔1 are antipodes (trans); equatorials 2,3,4 are all-cis to each other.
# SPY-5: apical (0) is "trans" to no basal vertex (the opposite is empty); basal pairs
#        2↔4 (diagonal across square base) are trans, 1↔3 also; adjacent are cis.
# T-shape: the two trans vertices are 0↔1 (the "T arms" 180°); 2 is the cis stem (90°).
_ANTIPODE = {
    "octahedron": {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4},
    "square_planar": {0: 2, 1: 3, 2: 0, 3: 1},
    "trigonal_bipyramid": {0: 1, 1: 0},        # only axial pair has a trans partner
    "square_pyramid": {1: 3, 3: 1, 2: 4, 4: 2},  # basal diagonals (apical has no trans)
    "tshape": {0: 1, 1: 0},                    # T-arms (vertex 2 = stem has no trans)
}


def _classify_coloring(geom_key, vertex_elems) -> str:
    """Universal scientific isomer name from per-vertex DONOR ELEMENTS + the
    polyhedron antipode structure (element-based trans-pair analysis, matching the
    project's _classify_isomer_label scheme — no system-specific rules):
      MA4B2 -> cis/trans · MA3B3 -> fac/mer · MA2B2C2 -> all-cis/all-trans/El-trans
      (square-planar MA2B2 -> cis/trans).  Returns '' for single-isomer cases."""
    from collections import Counter
    anti = _ANTIPODE.get(geom_key)
    if anti is None:
        return ""
    cnt = Counter(vertex_elems)
    n = len(vertex_elems)

    def is_trans(el):
        v = [i for i, e in enumerate(vertex_elems) if e == el]
        return any(anti.get(v[a]) == v[b] for a in range(len(v)) for b in range(a + 1, len(v)))

    pairs2 = [el for el, c in cnt.items() if c == 2]
    threes = [el for el, c in cnt.items() if c == 3]
    if n == 6:
        if threes:                                    # MA3B3 / MA3B2C
            return "mer" if is_trans(threes[0]) else "fac"
        if len(pairs2) == 3:                           # MA2B2C2
            trans_els = sorted(el for el in pairs2 if is_trans(el))
            if not trans_els:
                return "all-cis"
            if len(trans_els) == 3:
                return "all-trans"
            return "-".join(f"{e}trans" for e in trans_els)
        if len(pairs2) == 1:                           # MA4B2
            return "trans" if is_trans(pairs2[0]) else "cis"
    elif n == 4:                                       # SP-4 / T-4 MA2B2
        if len(pairs2) == 1:
            return "trans" if is_trans(pairs2[0]) else "cis"
    elif n == 5:                                       # TBP-5 / SPY-5
        if threes:                                     # MA3B2 (or MA3B1C1)
            return "mer" if is_trans(threes[0]) else "fac"
        if len(pairs2) == 1:                           # MA4B (degenerate; MA3B2 if 2nd singleton)
            return "trans" if is_trans(pairs2[0]) else "cis"
    elif n == 3:                                       # SP-3 trigonal-planar / T-3
        if len(pairs2) == 1:                           # MA2B
            return "trans" if is_trans(pairs2[0]) else "cis"
    return ""


def _xyz(syms, P) -> str:
    # HEADER-LESS atom block in the CANONICAL converter format ("{sym:4s}
    # {x:12.6f}..."), byte-identical to every other pool so viewers (Avogadro,
    # etc.) treat fffree output exactly like all other archives.  The pool
    # evaluator prepends "{count}\n{comment}".
    return "\n".join(f"{s:4s} {float(x):12.6f} {float(y):12.6f} {float(z):12.6f}"
                     for s, (x, y, z) in zip(syms, P))


def _lig_groups_from_vertex_specs(vertex_specs):
    """Construction-order ligand layout for the MONODENTATE heteroleptic path.

    The assembled frame is [metal] + AddHs(frag) blocks in vertex_specs order, so
    ligand i occupies the contiguous global-index block starting at 1 + sum(prev
    block sizes); its donor sits at block_start + donor_local_idx.  Returns the
    lig_groups list backbone_reembed.reembed_complex consumes, or None on failure.
    Universal, deterministic, graph-only (no coordinates needed)."""
    from rdkit import Chem as _Chem
    groups = []
    pos = 1
    for (frag, di) in vertex_specs:
        try:
            n = _Chem.AddHs(frag).GetNumAtoms()
        except Exception:
            return None
        groups.append({"mol": frag, "global_idxs": list(range(pos, pos + n)),
                       "donor_local": [int(di)]})
        pos += n
    return groups


def _lig_groups_from_config(config, ligands):
    """Construction-order ligand layout for the CHELATE assemble_from_config path.

    assemble_from_config iterates `by_lig` (unique ligand-index order from config)
    and appends AddHs(lg.mol) blocks; the donors are the ligand's donor_local_idxs.
    Reproduce that ordering exactly so each ligand's global-index block + donor
    locals map onto the native frame.  Returns lig_groups or None."""
    from rdkit import Chem as _Chem
    by_lig = {}
    for v, (li, arm) in config.items():
        by_lig.setdefault(li, []).append((v, arm))
    groups = []
    pos = 1
    for li in by_lig:                       # dict preserves insertion order (py3.7+)
        lg = ligands[li]
        try:
            n = _Chem.AddHs(lg["mol"]).GetNumAtoms()
        except Exception:
            return None
        dons = [int(x) for x in lg["donor_local_idxs"]]
        groups.append({"mol": lg["mol"], "global_idxs": list(range(pos, pos + n)),
                       "donor_local": dons})
        pos += n
    return groups


# vdW-level inter-ligand clash floor (Å) for the NEW-frame never-worse gate below.
# A re-embedded / re-seated conformer must not introduce a non-bonded heavy-heavy
# contact below this (a backbone folding into a NEIGHBOUR ligand collapses well
# inside the vdW shell long before the 0.60·Σcov gross-overlap floor _build_is_clean
# uses fires — measured: ACEQUY Fe(N(Dipp)(SiMe3))3 reembed frames at 1.98-2.01 A
# inter-ligand C-C, base frame 2.07-2.38 A).
_INTERLIG_VDW_FLOOR = 2.0
_INTERLIG_VDW_TOL = 0.05


def _interlig_vdw_gate_enabled() -> bool:
    """vdW-level inter-ligand clash filter for the ADDITIONAL conformer frames
    (backbone re-embed / conformer re-seating).  Active by default whenever those
    frames exist; toggle off with DELFIN_FFFREE_INTERLIG_VDW_GATE=0.  Byte-identical
    to the candidate when the reembed/seating flags themselves are off, because no
    extra frame is produced for it to filter (no-op)."""
    return os.environ.get("DELFIN_FFFREE_INTERLIG_VDW_GATE", "1") == "1"


def _min_nonbonded_heavy(syms, P) -> float:
    """Minimum NON-BONDED heavy-heavy distance (Å) in a frame.

    A pair is BONDED (and skipped) when d < 1.30·(rcov_i + rcov_j) -- the same
    graph-free covalent criterion as ``_bd._geometric_bonds`` -- so a genuine bond
    is never mistaken for a clash.  Metal atoms (and any metal-donor pair) are
    excluded: the metal is never a heavy-heavy partner here, so M-D contacts count
    as bonded by construction.  Hydrogens are skipped (heavy-heavy only, matching the
    vdW detector).  Returns +inf when no non-bonded heavy pair exists.  Pure geometry,
    deterministic, never raises on finite input."""
    P = np.asarray(P, dtype=float)
    n = len(syms)
    best = float("inf")
    for i in range(n):
        if syms[i] == "H" or _bd._is_metal(syms[i]):
            continue
        pi = P[i]
        for j in range(i + 1, n):
            if syms[j] == "H" or _bd._is_metal(syms[j]):
                continue
            d = float(np.linalg.norm(pi - P[j]))
            if d >= 1.30 * _bd._ideal_bond(syms[i], syms[j]):   # non-bonded only
                if d < best:
                    best = d
    return best


def _interlig_clash_ok(syms, P, base_min) -> bool:
    """NEVER-WORSE inter-ligand vdW gate for a NEW conformer frame (re-embed/re-seat).

    Reject the frame if its minimum non-bonded heavy-heavy distance is below
    ``max(_INTERLIG_VDW_FLOOR, base_min·(1 - _INTERLIG_VDW_TOL))`` -- i.e. the new
    conformer may NOT introduce an inter-ligand contact worse than the base frame
    already has (within a 5 % tolerance), and never below the hard 2.0 A vdW floor.
    ``base_min`` is the base (accepted) frame's own min non-bonded heavy-heavy
    distance; when unavailable (None / non-finite) the hard floor alone applies.
    This catches the per-ligand backbone-reembed collapse where one ligand folds into
    a neighbour (the reembed step freezes the core + re-folds each ligand WITHOUT
    inter-ligand awareness).  ``_build_is_clean``'s 0.60·Σcov gross-overlap floor
    (~0.9 A for C-C) sits far below the vdW shell and never fires for these ~2 A
    contacts.  Deterministic, geometry-only."""
    new_min = _min_nonbonded_heavy(syms, P)
    floor = _INTERLIG_VDW_FLOOR
    if base_min is not None and np.isfinite(base_min):
        floor = max(floor, float(base_min) * (1.0 - _INTERLIG_VDW_TOL))
    return new_min >= floor


def _append_reembed(results, metal, lig_groups, base_syms, base_P, base_label,
                    cn=None, geom=None, donors=None):
    """Backbone re-embed source (Task 2026-06-18, env DELFIN_FFFREE_BACKBONE_REEMBED).

    Given an ACCEPTED native base frame (base_syms, base_P) and its construction-order
    lig_groups, generate core-preserving GLOBAL backbone re-embed frames (metal + ALL
    donors frozen at their native positions; only the ligand BACKBONE re-folded via a
    fresh ETKDG embed grafted by rigid donor-Kabsch) and APPEND the self-gate-clean
    ones to `results`.  Additive: a re-embedded frame is only added on top of the
    native frame, never replaces it, and each is self-gated so it is never worse than
    legacy.  Default OFF -> no-op (byte-identical).  Never raises."""
    if not _BR.enabled() or lig_groups is None:
        return
    try:
        frames = _BR.reembed_complex(metal, lig_groups, (list(base_syms), base_P))
    except Exception:
        return
    # NEVER-WORSE inter-ligand vdW gate: each re-embedded frame re-folds every ligand
    # independently with the core frozen, WITHOUT inter-ligand awareness, so a backbone
    # can collapse into a neighbour ligand (~2 A C-C) — far above _build_is_clean's
    # 0.60·Σcov gross-overlap floor (~0.9 A), so the self-gate misses it.  base_min is
    # the accepted base frame's own min non-bonded heavy-heavy distance; a re-embedded
    # frame must not introduce a worse (closer) inter-ligand contact.
    _gate = _interlig_vdw_gate_enabled()
    base_min = _min_nonbonded_heavy(base_syms, base_P) if _gate else None
    for fi, (syms, P) in enumerate(frames):
        try:
            syms, P = _maybe_relax(syms, P)
            if not _build_is_clean(syms, P, cn=cn, geom=geom, donors=donors):
                continue
            if _gate and not _interlig_clash_ok(syms, P, base_min):
                continue                    # new conformer collapses inter-ligand -> drop
            results.append((_xyz(syms, P), f"{base_label}-reembed{fi+1}"))
        except Exception:
            continue


def _seating_enabled() -> bool:
    """Conformer-aware seating for large ligands (env DELFIN_FFFREE_CONFORMER_SEATING,
    default OFF -> byte-identical).  When ON, the decompose heavy-cap is raised so
    large-ligand complexes reach this build, and large ligands whose rigid placement
    fails the self-gate are re-seated by sampling conformers (core frozen ±0.05 A) and
    keeping the first clean fold (else legacy fallback)."""
    return os.environ.get("DELFIN_FFFREE_CONFORMER_SEATING", "0") == "1"


def _has_large_ligand(lig_groups) -> bool:
    """True if any ligand carries more heavy atoms / donor arm than the DEFAULT cap (8)
    — i.e. it is a complex that only reached this build because seating raised the cap.
    Conformer re-seating is engaged ONLY for these (cheap ligands seat fine rigidly)."""
    if not lig_groups:
        return False
    for lg in lig_groups:
        try:
            nheavy = sum(1 for a in lg["mol"].GetAtoms() if a.GetAtomicNum() > 1)
            dent = max(len(lg.get("donor_local", [1])), 1)
        except Exception:
            continue
        if nheavy / dent > DEC._HEAVY_CAP_DEFAULT:
            return True
    return False


def _seat_via_conformers(metal, lig_groups, base_syms, base_P,
                         cn=None, geom=None, donors=None):
    """Conformer-aware seating fallback for a large-ligand build that FAILED the
    self-gate under rigid placement (Task 2026-06-18, the dominant reach lever).

    Reuses the core-preserving conformer machinery (backbone_reembed.reembed_complex):
    metal + ALL donors are FROZEN on their native (ideal-polyhedron) vertices and only
    the ligand backbone is re-folded via a fresh global ETKDG/DG embed grafted by rigid
    donor-Kabsch.  Returns the FIRST re-seated fold that passes the self-gate (clash-
    free, non-collapsed, in-shell), or None if no conformer seats cleanly -> the caller
    bails to legacy (never-worse).  FF-free (geometry sampling, no metal-core relax),
    deterministic (fixed ETKDG seeds), core frozen ±0.05 A (reembed_complex hard guard).
    Never raises."""
    if lig_groups is None:
        return None
    try:
        frames = _BR.reembed_complex(metal, lig_groups, (list(base_syms), base_P))
    except Exception:
        return None
    # The rigid build FAILED the self-gate, so there is no accepted base frame whose
    # inter-ligand contact a re-seated fold must merely match -> apply the hard vdW
    # floor alone (base_min unavailable).  Same per-ligand core-frozen re-fold as
    # reembed: a seated fold can still drop a backbone into a neighbour ligand.
    _gate = _interlig_vdw_gate_enabled()
    for syms, P in frames:
        try:
            syms, P = _maybe_relax(syms, P)
            if not _build_is_clean(syms, P, cn=cn, geom=geom, donors=donors):
                continue
            if _gate and not _interlig_clash_ok(syms, P, None):
                continue                    # seated fold collapses inter-ligand -> skip
            return syms, P
        except Exception:
            continue
    return None


def _build_is_clean(syms, P, cn=None, geom=None, donors=None, exempt_pairs=None) -> bool:
    """Self-gate: reject a build that is destroyed — non-finite coordinates,
    any collapsed heavy-heavy bond, gross steric overlap, or OVER-COORDINATION
    (a non-coordinating atom intruding into the metal's first shell) — so fffree
    NEVER emits a structure worse than the legacy fallback would.  A failing build
    makes the whole complex fall back to the legacy pipeline (return None),
    guaranteeing fffree is never worse than UFF on its addressable subset.
    Universal, geometry-only.  Disable via DELFIN_FFFREE_SELFGATE=0.

    ``donors`` (optional, global indices of the cn constructed donor atoms): when
    fffree KNOWS the coordinating atoms (it built them), the coordination-shape and
    over-coordination checks use the donors directly instead of the cn-closest-heavy
    heuristic.  This is essential for CHELATES, whose ring backbone legitimately sits
    ~2.4-2.9 A from the metal: the heuristic miscounts a ring carbon as a donor (wrong
    cshm) or as over-coordination, falsely rejecting correct chelate geometry.

    ``exempt_pairs`` (optional): heavy-heavy bonds whose SHORT length is chemically
    correct (genuine triple/multiple/aromatic bonds such as a C≡O carbonyl ~1.12 A,
    a C≡N nitrile ~1.16 A, an azo N=N ~1.08 A, or an aromatic/imine C=N ~1.14 A), so
    they are NOT counted as collapsed.  Accepts either a set/iterable of (i,j) pairs
    (UNCONDITIONAL exemption; the historic hapto path) or a mapping
    {(i,j): ideal_multibond_length} (LENGTH-GATED exemption: only honoured when the
    observed length is within DELFIN_FFFREE_MULTIBOND_TOL of the ideal, so a real
    sub-ideal collapse is still caught).  Default None = byte-identical."""
    if os.environ.get("DELFIN_FFFREE_SELFGATE", "1") == "0":
        return True
    P = np.asarray(P, dtype=float)
    if P.size == 0 or not np.all(np.isfinite(P)):
        return False
    syms = list(syms)
    bonds = _bd._geometric_bonds(syms, P)
    # X-H collapse calibration (#306/#281): X-ray C-H/N-H/O-H bonds are legitimately
    # SHORT (~0.65-0.95 A); the uniform 0.82*ideal floor (=0.88 A for C-H) FALSELY
    # flags them as "collapsed", which drops an otherwise-clean geometry isomer to the
    # legacy fallback (measured: AFOWOH's trans-Cl isomer rejected over 2 aryl C-H at
    # 0.81/0.87 A — 0.01 A under threshold).  When enabled, H-involving bonds use a
    # lower floor (XH_FRAC*ideal, default 0.55 ~= the CCDC X-H metric floor) so a real
    # H-on-atom collapse (< ~0.6 A) is still caught while X-ray-short H pass.  Heavy-
    # heavy bonds keep the 0.82 floor unchanged.  Env-gated, default OFF => byte-id
    # (floor is 0.82*ideal for every bond, exactly the historic _count_collapsed /
    # exempt-branch behaviour).
    _xh = os.environ.get("DELFIN_FFFREE_XH_COLLAPSE", "0") == "1"
    _xh_frac = float(os.environ.get("DELFIN_FFFREE_XH_COLLAPSE_FRAC", "0.55"))

    def _coll_floor(a, b):
        if _xh and (a == "H" or b == "H"):
            return _xh_frac * _bd._ideal_bond(a, b)
        return 0.82 * _bd._ideal_bond(a, b)

    # exempt_pairs (#279/#281, DELFIN_FFFREE_MULTIBOND_EXEMPT): heavy-heavy bonds whose
    # SHORT length is chemically correct (genuine double/triple/aromatic bonds such as a
    # metal-carbonyl C≡O ~1.12 A, nitrile C≡N ~1.16 A, azo N=N ~1.08 A, or aromatic/imine
    # C=N ~1.14 A).  Two accepted shapes, both byte-identical when unused:
    #   * plain iterable of (i, j)            -> UNCONDITIONAL exemption (the historic
    #                                            hapto-path behaviour; kept unchanged).
    #   * mapping {(i, j): ideal_multibond}   -> LENGTH-GATED exemption: the bond is
    #                                            exempted ONLY when its observed length is
    #                                            within MULTIBOND_TOL (default 0.15 A) of
    #                                            the bond-order-appropriate ideal.  A real
    #                                            embedding COLLAPSE (e.g. a C-O at 0.76 A,
    #                                            far below even the triple-bond ideal) is
    #                                            NOT exempted and is STILL caught.  SAFETY.
    _mb_tol = float(os.environ.get("DELFIN_FFFREE_MULTIBOND_TOL", "0.15"))
    if isinstance(exempt_pairs, dict):
        _ex_len = {(min(i, j), max(i, j)): float(t) for (i, j), t in exempt_pairs.items()}
        _ex = set()                                   # length-gated, evaluated per bond
    else:
        _ex_len = {}
        _ex = {(min(i, j), max(i, j)) for i, j in (exempt_pairs or ())}

    def _is_exempt(i, j, d):
        key = (min(i, j), max(i, j))
        if key in _ex:                                # unconditional (hapto-path set form)
            return True
        ideal_mb = _ex_len.get(key)                   # length-gated (multibond dict form)
        return ideal_mb is not None and d >= ideal_mb - _mb_tol

    n_coll = 0
    for i, j in bonds:
        if _bd._is_metal(syms[i]) or _bd._is_metal(syms[j]):
            continue
        d_ij = float(np.linalg.norm(P[i] - P[j]))
        if d_ij >= _coll_floor(syms[i], syms[j]):
            continue
        if _is_exempt(i, j, d_ij):                    # genuine short multiple bond -> pass
            continue
        n_coll += 1
    if n_coll > 0:
        return False
    bset = {(min(i, j), max(i, j)) for i, j in bonds}
    n = len(syms)
    for i in range(n):
        if syms[i] == "H" or _bd._is_metal(syms[i]):
            continue
        for j in range(i + 1, n):
            if syms[j] == "H" or _bd._is_metal(syms[j]):
                continue
            if (i, j) in bset:
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < 0.60 * _bd._ideal_bond(syms[i], syms[j]):   # gross overlap
                return False
    mi = next((i for i in range(n) if _bd._is_metal(syms[i])), None)
    donor_set = set(donors) if donors else None
    # UNDER-coordination / decoordination guard (#324b, env DELFIN_FFFREE_COORD_INTEGRITY,
    # default OFF -> byte-identical to the candidate; 2/1197 candidate base frames carry a
    # decoordinated donor, so it MUST be opt-in to keep the OFF build bit-identical).
    # A conformer/rotamer/backbone-reembed/seating frame that rotates a backbone torsion
    # CARRYING a coordinating donor swings that donor OFF the metal (measured: ALAXOA
    # scorpionate S 2.4->5.8 A in EVERY torsion frame; ALAHEB Ir-N 2.1->3.0-3.4 A in 11/23
    # frames; ~14 % of all expansion frames pool-wide).  The over-coordination test below
    # never catches it (a donor LEAVING the shell intrudes on nothing) and the cshm shape
    # test silently re-selects a backbone carbon for the departed donor, so the broken frame
    # passes and leaks into the manifold.  Reject when a KNOWN donor sits farther from the
    # metal than its donor-TYPE ideal M-D (polyhedra.md_distance) + slack.  Donor-type-aware
    # by construction: the same 3.0 A is decoordination for a 2.1 A Fe-N donor but legitimate
    # stretch for a 2.7 A W-S donor.  The slack is generous (legit conformer spread is
    # <= ideal+0.6 A at p95, decoordination is ideal+1.7 A at p98) so a valid frame is never
    # rejected (verified: keeps AKEMUY's crystal-matching frame #13 and all clean frames,
    # drops only the donor-off frames).  Geometry-only, deterministic, never raises.
    if (donor_set is not None and mi is not None
            and os.environ.get("DELFIN_FFFREE_COORD_INTEGRITY", "0") == "1"):
        _coord_slack = float(os.environ.get("DELFIN_FFFREE_COORD_INTEGRITY_SLACK", "0.85"))
        for _d in donor_set:
            try:
                _ideal_md = float(PLY.md_distance(syms[mi], syms[_d]))
            except Exception:
                _ideal_md = 2.2
            if float(np.linalg.norm(P[_d] - P[mi])) > _ideal_md + _coord_slack:
                return False
    # over-coordination / spurious intrusion into the metal's first shell.
    if cn and mi is not None:
        if donor_set is not None:
            # donor-aware: reject only a NON-donor heavy atom that sits IN FRONT of
            # the coordination shell (closer than the donors) = real collapse; a
            # chelate ring backbone atom at normal ring distance (>= ~donor shell)
            # is legitimate and must pass.
            md = [float(np.linalg.norm(P[d] - P[mi])) for d in donor_set]
            md_min = min(md) if md else 0.0
            for j in range(n):
                if j == mi or j in donor_set or syms[j] == "H":
                    continue
                if float(np.linalg.norm(P[j] - P[mi])) < 0.92 * md_min:
                    return False
        else:
            close = 0
            for j in range(n):
                if j == mi or syms[j] == "H":
                    continue
                cutoff = max(1.45 * _bd._ideal_bond(syms[mi], syms[j]), 2.7)
                if float(np.linalg.norm(P[j] - P[mi])) < cutoff:
                    close += 1
            if close > cn + 1:                          # +1 slack for borderline
                return False
    # #39: reject catastrophic coordination-SHAPE outliers (CShM >> typical sets the
    # worst-case poly_max/cshm_max above UFF; legacy is better for that tail).
    # Threshold sits deep in the valley (p75 0.14 <-> p90 10.7).  Env DELFIN_FFFREE_SHAPE_MAX.
    if cn and geom and mi is not None:
        _shmax = float(os.environ.get("DELFIN_FFFREE_SHAPE_MAX", "20.0"))
        # High-CN (CN>=7) placement is less reliable than CN4-6, so a build that
        # only passes the loose CN4-6 threshold can still be worse than the legacy
        # fallback there.  A tighter high-CN shape gate (default 5.0) makes the
        # high-CN subset cleanly net-better than legacy (measured: net +11 vs +9 at
        # 20, regressions 10->6).  Deterministic, CN-keyed; env-tunable.
        if cn >= 7:
            _shmax = min(_shmax, float(os.environ.get("DELFIN_FFFREE_SHAPE_MAX_HIGHCN", "5.0")))
        if donor_set is not None and len(donor_set) == cn:
            sel = list(donor_set)                       # the KNOWN constructed donors
        else:
            ds = sorted((float(np.linalg.norm(P[j] - P[mi])), j)
                        for j in range(n) if j != mi and syms[j] != "H")
            sel = [j for _, j in ds[:cn]]
        if len(sel) >= cn:
            obs = np.array([(P[j] - P[mi]) / (max(float(np.linalg.norm(P[j] - P[mi])), 1e-9))
                            for j in sel])
            try:
                if PLY.cshm(obs, geom) > _shmax:
                    return False
            except Exception:
                pass
    return True


def _fffree_chelate_isomers(d, geom_key, max_isomers):
    """Build all distinct isomers of a chelate-containing complex (mixed bi-/
    monodentate) via the universal chelate-config enumerator + per-config
    geometric assembly.  Returns [(xyz, label), ...] or None."""
    ligands = d["ligands"]
    if any(lg["denticity"] >= 4 for lg in ligands):
        return None        # kappa>=4 (porphyrin/salen/DTPA) not yet supported -> legacy
    # Aromatic donors on the NEWLY-enabled CN5 chelate geometries (TBP-5/SPY-5) would be
    # placed face-on (ring-normal perp to M-N): _vsepr_reconstruct skips ring donors and
    # there is no in-plane orientation yet (deferred to the aromatic-N-in-plane iter).  The
    # self-gate catches collapse/over-coord/shape but NOT face-on, so route aromatic CN5
    # chelates to legacy (never-worse).  Scoped to TBP-5/SPY-5 only -> existing OC-6/SP-4
    # chelate coverage (e.g. M(bipy)3) stays byte-identical.
    if geom_key in ("trigonal_bipyramid", "square_pyramid"):
        for lg in ligands:
            lmol = lg["mol"]
            if any(lmol.GetAtomWithIdx(i).GetIsAromatic()
                   for i in lg.get("donor_local_idxs", [])):
                return None
    specs = []
    for lg in ligands:
        specs.append({
            "type": Chem.MolToSmiles(lg["mol"]),
            "denticity": lg["denticity"],
            "asym": len(set(lg.get("donor_elems", []))) > 1,
            # geometry-aware meridional restriction (default OFF: rigid_planar is
            # always False unless DELFIN_FFFREE_PLANAR_MER=1, so byte-identical)
            "rigid_planar": bool(lg.get("rigid_planar")),
        })
    try:
        configs = PIC.enumerate_chelate_configs(geom_key, specs)
    except Exception:
        return None
    if not configs:
        return None
    geom_tag = d["geometry"].split()[0]
    results = []
    # SIGMA-ensemble (Task A.1): the chelate σ sub-path emits ONE frame per config,
    # but the legacy converter sprays 4-25 frames per refcode (chelate-ring PUCKER +
    # ligand conformers) and best-of-ensemble MIN crystal-recall rewards the larger
    # spray purely on shot count (SARKOM/QEPFAS/UWUJAY trail legacy with the FIRST
    # coordination shell already correct).  Generalise the SAME proven lever (CN2
    # +5.9pp, hapto +5.8pp): per config emit a deterministic, RMSD-deduped ensemble.
    # The metallacycle / ETKDG conformer pool already samples chelate-ring pucker
    # (Cremer-Pople) across conformers, so assemble_from_config(n_frames>1) RETAINS
    # the distinct low-clash conformer combinations instead of collapsing to one.
    # Env-gated DELFIN_FFFREE_SIGMA_ENSEMBLE, default OFF (=> byte-identical, single
    # frame per config) ; size capped near the legacy spray.
    _sig_ens = os.environ.get("DELFIN_FFFREE_SIGMA_ENSEMBLE", "0") == "1"
    _n_chel = int(os.environ.get("DELFIN_FFFREE_SIGMA_CHELATE_NFRAMES", "8"))
    for k, config in enumerate(configs[:max_isomers]):
        # per-config pruning (generate-gate-floor): a geometrically infeasible config
        # (e.g. a fac vertex-triple for a mer pincer) is SKIPPED, not bailed -- so one
        # bad isomer no longer drops the whole complex to legacy.  Clean isomers survive.
        # #279/#281: genuine short multiple/aromatic bonds (global, length-gated) for the
        # collapse self-gate (e.g. AWELOD's aromatic C=N at 1.14 A).  The chelate builder
        # places ligands in first-appearance config order; offsets mirror that exactly.
        # Empty when DELFIN_FFFREE_MULTIBOND_EXEMPT unset -> byte-identical.
        _ex = _exempt_from_blocks(_config_block_offsets(config, ligands))
        if _sig_ens:
            # NEVER-WORSE GUARD: emit the ensemble for a config ONLY IF its
            # single-frame (clash-minimal) build ALSO passes the self-gate.  This
            # keeps the ENSEMBLE strictly ADDITIVE to the single-frame path's
            # accepted-config set: a config that the single-frame path would have
            # SKIPPED (-> potentially the whole complex falls back to legacy, which
            # may hold the crystallised structure) must NOT be rescued by the
            # larger conformer pool, or a refcode that legacy was winning (e.g.
            # KADZUL: single-frame None -> legacy 0.74Å) would lose its fallback to
            # a geometrically inferior FF-free frame (2.24Å).  Same self-gate, same
            # decision; the ensemble only adds conformer DIVERSITY on top.
            try:
                single = AC.assemble_from_config(d["metal"], d["geometry"], config, ligands)
            except Exception:
                single = None
            if single is None:
                continue
            ssyms, sP, sdonors = single
            ssyms, sP = _maybe_relax(ssyms, sP)
            if not _build_is_clean(ssyms, sP, cn=d.get("cn"), geom=d.get("geometry"),
                                   donors=sdonors, exempt_pairs=_ex):
                continue                                 # single-frame would skip -> skip
            try:
                built = AC.assemble_from_config(d["metal"], d["geometry"], config,
                                                ligands, n_frames=_n_chel)
            except Exception:
                built = None
            if not built:                                # ensemble failed -> use the single frame
                results.append((_xyz(ssyms, sP), f"{geom_tag}-chelate-{k+1}"))
                _append_reembed(results, d["metal"],
                                _lig_groups_from_config(config, ligands),
                                ssyms, sP, f"{geom_tag}-chelate-{k+1}",
                                cn=d.get("cn"), geom=d.get("geometry"), donors=sdonors)
                continue
            # re-embed off the canonical single frame (clash-minimal) for this config
            _append_reembed(results, d["metal"],
                            _lig_groups_from_config(config, ligands),
                            ssyms, sP, f"{geom_tag}-chelate-{k+1}",
                            cn=d.get("cn"), geom=d.get("geometry"), donors=sdonors)
            for fi, (syms, P, donors) in enumerate(built):
                syms, P = _maybe_relax(syms, P)
                if not _build_is_clean(syms, P, cn=d.get("cn"), geom=d.get("geometry"),
                                       donors=donors, exempt_pairs=_ex):   # per-frame self-gate
                    continue                             # skip a bad frame, keep clean ones
                lab = f"{geom_tag}-chelate-{k+1}" + (f"-conf{fi+1}" if fi else "")
                results.append((_xyz(syms, P), lab))
            continue
        try:
            built = AC.assemble_from_config(d["metal"], d["geometry"], config, ligands)
        except Exception:
            continue
        if built is None:
            continue
        syms, P, donors = built
        syms, P = _maybe_relax(syms, P)
        _clg = _lig_groups_from_config(config, ligands)
        if not _build_is_clean(syms, P, cn=d.get("cn"), geom=d.get("geometry"),
                               donors=donors, exempt_pairs=_ex):   # donor-aware self-gate
            # Conformer-aware seating (DELFIN_FFFREE_CONFORMER_SEATING, default OFF):
            # a large-arm chelate config whose rigid build fails the self-gate is
            # re-seated (donors frozen on the native vertices; backbone re-folded).
            # On success the clean fold is used for THIS config; on failure the config
            # is SKIPPED as before (never-worse).  Byte-identical when off (skip).
            reseated = None
            if _seating_enabled() and _has_large_ligand(_clg):
                reseated = _seat_via_conformers(d["metal"], _clg, syms, P,
                                                cn=d.get("cn"), geom=d.get("geometry"),
                                                donors=donors)
            if reseated is None:
                continue                          # skip this config
            syms, P = reseated
        _lab = f"{geom_tag}-chelate-{k+1}"
        results.append((_xyz(syms, P), _lab))
        # Backbone re-embed (env DELFIN_FFFREE_BACKBONE_REEMBED, default OFF): add
        # core-preserving global-fold variants of this accepted chelate frame.
        _append_reembed(results, d["metal"], _clg,
                        syms, P, _lab, cn=d.get("cn"), geom=d.get("geometry"),
                        donors=donors)
    return results or None


def _hapto_subst_rotamers(xyz, n_per=2):
    """Add substituent / co-ligand rotamers to ONE assembled rigid-hapto frame using
    the project FF-free σ-rotamer machinery (_rotamer_diversity.apply): it rotates
    distal sub-trees about rotatable bonds with the coordination shell frozen + a
    never-worse / M-D-invariant firewall, so the η-ring + donors stay rigid while
    peripheral methyls / phenyls / co-ligand arms sample distinct clock positions.
    Returns [xyz, rot1, ...] (base always first); falls back to [xyz] on any error.
    Deterministic (the rotamer grid + seeds are fixed inside the module)."""
    if os.environ.get("DELFIN_FFFREE_HAPTO_NO_SUBROT", "0") == "1" or n_per < 1:
        return [xyz]
    # Size guard: the OB rotor-tree rotamer search is O(atoms) per DOF; on very large
    # complexes (huge oligo-aryl substituents) it dominates wall-clock for little gain
    # (the η-ring + donors are frozen, so distal rotamers barely move heavy atoms near
    # the crystal core).  Skip it past a generous cap -> bounded per-SMILES time.
    if xyz.count("\n") + 1 > 160:
        return [xyz]
    try:
        from delfin import _rotamer_diversity as RD
        outs = RD.apply(xyz, n_per_isomer=int(n_per), n_states=3, max_dofs=4)
        return outs if outs else [xyz]
    except Exception:
        return [xyz]


def _fffree_hapto_isomers(d, max_isomers):
    """Build a hapto complex (≥1 η-face) on the FF-free path: rigid-η-unit
    construction keeps each Cp/arene/diene/allyl ring at its crystallographic
    metal→centroid distance instead of collapsing the ring carbons onto the metal
    (the legacy hapto defect).  Emits a deterministic, RMSD-deduplicated rigid
    ENSEMBLE (η-ring rotamers via symmetry-reduced spin, valence-gated η/σ ring-slip
    isomers, Cremer-Pople pucker for non-aromatic faces, plus substituent/co-ligand
    rotamers) comparable in size to the ~30-frame legacy spray, every member fully
    rigid (no collapse).  Returns [(xyz, label), ...] or None (-> legacy fallback).
    Env-gated upstream (the dict only carries 'has_eta' when
    DELFIN_FFFREE_RIGID_HAPTO=1).  Determinism: PYTHONHASHSEED=0 + fixed seeds."""
    geom_tag = d["geometry"].split()[0]
    cn = d.get("cn")
    geom = d.get("geometry")

    def _accept(built, label):
        if built is None:
            return None
        syms, P, donors, exempt_pairs = built
        syms, P = _maybe_relax(syms, P)
        if not _build_is_clean(syms, P, cn=cn, geom=geom,
                               donors=donors, exempt_pairs=exempt_pairs):
            return None                                 # never worse than legacy
        return (syms, P, label)

    # ENSEMBLE path (default for the rigid-hapto flag): enumerate distinct rigid
    # builds, gate each, then expand each accepted build with substituent rotamers.
    cap = max(1, int(max_isomers))
    try:
        builds = AC.assemble_hapto_ensemble(d["metal"], geom, d, max_builds=cap)
    except Exception:
        builds = None

    results = []
    seen_keys = set()
    if builds:
        n_subrot = 1 if len(builds) >= 8 else 2        # keep total near legacy size
        for bi, b in enumerate(builds):
            acc = _accept(b, None)
            if acc is None:
                continue
            syms, P, _ = acc
            base_xyz = _xyz(syms, P)
            tag = "hapto" if bi == 0 else f"hapto-iso{bi}"
            for ri, rot in enumerate(_hapto_subst_rotamers(base_xyz, n_per=n_subrot)):
                key = rot
                if key in seen_keys:
                    continue
                seen_keys.add(key)
                lab = f"{geom_tag}-{tag}" + (f"-r{ri}" if ri else "")
                results.append((rot, lab))
                if len(results) >= cap:
                    break
            if len(results) >= cap:
                break

    # --- Hebel A (DELFIN_FFFREE_HAPTO_AXIS_ROT, default OFF) --------------------
    # STRICTLY ADDITIVE: the baseline `results` above are emitted UNCHANGED (byte-id
    # when the flag is off).  With the flag on, APPEND extra frames that rotate the
    # CO-tripod / co-ligands / 2nd ring about the primary η M→centroid axis (a single
    # discrete rotamer DOF the base ensemble does not vary).  Every appended frame is
    # gated by the same never-worse ruler; the axis runs through the metal so M-D +
    # M-centroid distances are invariant.  Appended AFTER the cap so the baseline
    # (and its best-of-ensemble MIN) can never regress.
    if results and os.environ.get("DELFIN_FFFREE_HAPTO_AXIS_ROT", "0") == "1":
        n_axis = int(os.environ.get("DELFIN_FFFREE_HAPTO_AXIS_NROT", "8"))
        axis_cap = int(os.environ.get("DELFIN_FFFREE_HAPTO_AXIS_CAP", str(2 * cap)))
        try:
            ax_builds = AC.assemble_hapto_axis_rotants(
                d["metal"], geom, d, n_axis=n_axis, max_builds=axis_cap)
        except Exception:
            ax_builds = []
        for ai, b in enumerate(ax_builds):
            acc = _accept(b, None)
            if acc is None:
                continue
            syms, P, _ = acc
            rot = _xyz(syms, P)
            if rot in seen_keys:
                continue
            seen_keys.add(rot)
            results.append((rot, f"{geom_tag}-hapto-axis{ai+1}"))
            if len(results) >= cap + axis_cap:
                break

    if results:
        return results

    # Fallback: the single canonical build (historic v1 behaviour) so the rigid path
    # still fires even if the ensemble enumerator yields nothing clean.
    try:
        built = AC.assemble_hapto(d["metal"], geom, d)
    except Exception:
        return None
    acc = _accept(built, None)
    if acc is None:
        return None
    syms, P, _ = acc
    return [(_xyz(syms, P), f"{geom_tag}-hapto-1")]


def _coord_filter(results):
    """UNIVERSAL final decoordination filter (#324b coverage closure).

    The per-build self-gate (_build_is_clean coord-integrity branch) only covers the
    paths that route through it with the constructed donor set.  The DENSE generation
    paths -- backbone re-embed, conformer-coverage, conformer-seating, chelate-backbone,
    Cremer-Pople pucker, hapto-axis rotants, sigma-ensemble -- emit frames that bypass
    it, so under the full dense stack ~43% of frames can carry a decoordinated donor
    (a backbone torsion that carries a coordinating donor swung it off the metal).
    This final pass re-checks EVERY emitted frame regardless of which path produced it:
    a frame is dropped if any donor (the consistently-coordinated set, taken from the
    BEST-coordinated frame as reference) sits farther than its donor-type ideal M-D
    (polyhedra.md_distance) + slack.  Crystals pack TIGHT but stay COORDINATED, so this
    NEVER drops a crystal-like frame (verified: keeps AKEMUY #13) -- it only removes the
    ligand-flew-off junk (AGIPIT/ALAXOA).  Gated by DELFIN_FFFREE_COORD_INTEGRITY
    (default OFF -> returns results unchanged -> byte-identical).  Geometry-only, never raises."""
    if not results or os.environ.get("DELFIN_FFFREE_COORD_INTEGRITY", "0") != "1":
        return results
    try:
        slack = float(os.environ.get("DELFIN_FFFREE_COORD_INTEGRITY_SLACK", "0.85"))

        def _parse(xyz):
            ls = xyz.splitlines(); n = int(ls[0]); sy = []; P = []
            for ln in ls[2:2 + n]:
                p = ln.split()
                sy.append(p[0]); P.append([float(p[1]), float(p[2]), float(p[3])])
            return sy, np.asarray(P, dtype=float)

        # reference = the frame with the MOST donors within 2.75 A of the metal (the
        # best-coordinated frame), so a decoordinated base frame cannot hide a donor.
        ref_mi = ref_don = ref_n = None; best = -1
        parsed = []
        for xyz, lab in results:
            try:
                sy, P = _parse(xyz); parsed.append((sy, P, xyz, lab))
            except Exception:
                parsed.append((None, None, xyz, lab)); continue
            mi = next((i for i in range(len(sy)) if _bd._is_metal(sy[i])), None)
            if mi is None:
                continue
            don = [k for k in range(len(sy)) if k != mi and sy[k] != "H"
                   and float(np.linalg.norm(P[mi] - P[k])) < 2.75]
            if len(don) > best:
                best = len(don); ref_mi = mi; ref_don = don; ref_n = len(sy)
        if ref_mi is None or not ref_don:
            return results
        # threshold per donor (donor-type ideal M-D + slack), from the reference frame
        rsy = next(s for s, _, _, _ in parsed if s is not None and len(s) == ref_n)
        thr = {}
        for k in ref_don:
            try: thr[k] = float(PLY.md_distance(rsy[ref_mi], rsy[k])) + slack
            except Exception: thr[k] = 3.05
        kept = []
        for sy, P, xyz, lab in parsed:
            if sy is None or len(sy) != ref_n:
                kept.append((xyz, lab)); continue          # can't check -> keep
            if any(float(np.linalg.norm(P[ref_mi] - P[k])) > thr[k] for k in ref_don):
                continue                                    # a donor decoordinated -> drop
            kept.append((xyz, lab))
        return kept or [results[0]]                         # never empty -> keep base
    except Exception:
        return results


def _fffree_isomers(smiles: str, max_isomers: int = 50
                    ) -> Optional[List[Tuple[str, str]]]:
    d = DEC.decompose(smiles)
    if d is None:
        return None
    geom_key = _GEOM_TO_POLYA.get(d["geometry"])
    if geom_key is None or geom_key not in PIC._GROUPS:
        return None
    if d.get("has_eta"):
        return _coord_filter(_fffree_hapto_isomers(d, max_isomers))
    if d.get("has_chelate"):
        return _coord_filter(_fffree_chelate_isomers(d, geom_key, max_isomers))
    # ligand identity = canonical SMILES of each fragment; group by it
    lig_label, lig_ref, lab_elem = [], {}, {}
    for lg in d["ligands"]:
        try:
            lab = Chem.MolToSmiles(lg["mol"])
        except Exception:
            return None
        lig_label.append(lab)
        lig_ref.setdefault(lab, (lg["mol"], lg["donor_local_idx"]))
        lab_elem[lab] = lg["donor_elem"]
    spec = dict(Counter(lig_label))
    try:
        colorings = PIC.enumerate_isomers(geom_key, spec)
    except Exception:
        return None
    if not colorings:
        return None
    results: List[Tuple[str, str]] = []
    # CN2-ensemble (iter-32g): the rigid FF-free CN2 path emits ONE frame per coloring,
    # but the legacy multi-frame path sprays ~3-7 conformers, and best-of-ensemble MIN
    # rewards the larger spray on flexible-ligand CN2 (JEVDUJ/NESLED/YUXVUJ) — the same
    # pattern an ensemble already fixed for rigid-hapto.  So for CN2 emit a small,
    # RMSD-deduped conformer/rotamer ENSEMBLE per coloring (ligand internal conformers
    # + co-ligand orientation vary; the 180° linear core stays rigid), reusing the same
    # FF-free Layer-3 machinery (_ligand_confs_from_mol ETKDG pool + clash scoring +
    # refine).  Gated behind CN_EXTEND (already required for CN2 to exist); opt-out via
    # DELFIN_FFFREE_CN2_ENSEMBLE=0.  Default OFF overall (CN_EXTEND off => byte-id).
    _cn2_ens = (d.get("cn") == 2
                and os.environ.get("DELFIN_FFFREE_CN2_ENSEMBLE", "1") == "1")
    # SIGMA-ensemble (Task A.2): the monodentate CN4/5/6 path emits ONE frame per
    # coloring via assemble_heteroleptic_from_mols, but the legacy multi-frame
    # converter sprays ~8-13 conformers (cis/trans donor perms + ligand internal
    # conformers + co-ligand orientation) and best-of-ensemble MIN crystal-recall
    # rewards the larger spray (FANYAW/LEYLAC/CPOCEM trail legacy by frame count
    # alone — the first coordination shell is already correct).  This is the SAME
    # lever already proven for CN2 (+5.9pp) and rigid-hapto (+5.8pp).  Reuse the
    # identical CN2 ensemble machinery (assemble_heteroleptic_ensemble: ETKDG
    # conformer pool + axial-spin rotamers + inter-ligand combos, RMSD-deduped,
    # clash-aware, every frame self-gated).  Env-gated DELFIN_FFFREE_SIGMA_ENSEMBLE,
    # default OFF (=> byte-identical when unset).  Scoped to MONODENTATE CN4/5/6
    # (chelate + hapto + CN2/3 reach this branch only via their own paths/gates).
    _sigma_ens = (d.get("cn") in (4, 5, 6)
                  and os.environ.get("DELFIN_FFFREE_SIGMA_ENSEMBLE", "0") == "1")
    _n_ens = int(os.environ.get("DELFIN_FFFREE_CN2_NFRAMES", "8"))
    _n_sigma = int(os.environ.get("DELFIN_FFFREE_SIGMA_NFRAMES", "10"))
    for k, coloring in enumerate(colorings[:max_isomers]):
        vertex_specs = [lig_ref[lab] for lab in coloring]
        vertex_elems = [lab_elem[lab] for lab in coloring]
        name = _classify_coloring(geom_key, vertex_elems)
        geom_tag = d["geometry"].split()[0]
        base_label = f"{name}-{geom_tag}-{k+1}" if name else f"{geom_tag}-{k+1}"
        # #279/#281: genuine short multiple/aromatic bonds (global, length-gated) for the
        # collapse self-gate.  Empty when DELFIN_FFFREE_MULTIBOND_EXEMPT unset -> byte-id.
        _ex = _exempt_from_blocks(_heteroleptic_block_offsets(vertex_specs))
        if _cn2_ens or _sigma_ens:
            _nf = _n_ens if _cn2_ens else _n_sigma
            try:
                ens = AC.assemble_heteroleptic_ensemble(
                    d["metal"], d["geometry"], vertex_specs, n_frames=_nf)
            except Exception:
                ens = None
            if not ens:
                return None
            kept = 0
            for fi, (syms, P) in enumerate(ens):
                syms, P = _maybe_relax(syms, P)
                if not _build_is_clean(syms, P, cn=d.get("cn"), geom=d.get("geometry"),
                                       exempt_pairs=_ex):
                    continue            # skip a bad frame; keep the clean ones
                lbl = base_label if fi == 0 else f"{base_label}-conf{fi+1}"
                results.append((_xyz(syms, P), lbl))
                if fi == 0:             # re-embed off the canonical (clash-minimal) frame
                    _append_reembed(results, d["metal"],
                                    _lig_groups_from_vertex_specs(vertex_specs),
                                    syms, P, lbl, cn=d.get("cn"), geom=d.get("geometry"))
                kept += 1
            if kept == 0:               # whole coloring unbuildable -> legacy (never-worse)
                return None
            continue
        try:
            built = AC.assemble_heteroleptic_from_mols(d["metal"], d["geometry"], vertex_specs)
        except Exception:
            return None
        if built is None:
            return None
        syms, P = built
        syms, P = _maybe_relax(syms, P)
        _lg = _lig_groups_from_vertex_specs(vertex_specs)
        if not _build_is_clean(syms, P, cn=d.get("cn"), geom=d.get("geometry"),
                               exempt_pairs=_ex):   # self-gate: destroyed/over-coord/shape-outlier
            # Conformer-aware seating (DELFIN_FFFREE_CONFORMER_SEATING, default OFF):
            # large ligands (raised heavy-cap) often FAIL the rigid placement self-gate
            # because their backbone folds into the coordination shell.  Re-seat them by
            # sampling conformers with the metal + donors FROZEN on the ideal vertices
            # (±0.05 A guard) and keep the first clean fold; only large-ligand complexes
            # are re-seated (cheap ligands seat fine rigidly).  No clean fold -> legacy
            # (never-worse).  Byte-identical when the flag is off (this branch returns).
            if not (_seating_enabled() and _has_large_ligand(_lg)):
                return None
            reseated = _seat_via_conformers(d["metal"], _lg, syms, P,
                                            cn=d.get("cn"), geom=d.get("geometry"))
            if reseated is None:
                return None                 # no conformer seats cleanly -> legacy
            syms, P = reseated
        label = base_label
        results.append((_xyz(syms, P), label))
        # Backbone re-embed (env DELFIN_FFFREE_BACKBONE_REEMBED, default OFF): add
        # core-preserving global-fold variants of THIS accepted native frame.
        _append_reembed(results, d["metal"], _lg, syms, P, label,
                        cn=d.get("cn"), geom=d.get("geometry"))
    # CN5 polytopal completeness (#coverage): decompose defaults CN5 -> TBP-5, but SPY-5
    # is the Berry-pseudorotation partner — real CN5 complexes split between the two.
    # Additively enumerate SPY-5 too (best-effort; never bails the TBP-5 result).
    if d.get("cn") == 5:
        results += _enumerate_geometry(d, "square_pyramid", "SPY-5 square pyramid",
                                       lig_ref, lab_elem, spec, max_isomers)
    # Iter-31 (User 2026-05-28): CN6 dual OC-6 / TPR-6.  decompose defaults CN6 -> OC-6
    # but early-TM Mo/W CN6 prefer TPR — coverage gap previously missed (no TPR-6 in
    # the FF-free Pólya enumerator).  Additive, env-gated default OFF (byte-identical
    # when unset).  Same pattern as CN5 SPY-5: best-effort, never bails OC-6 result.
    if d.get("cn") == 6 and os.environ.get("DELFIN_FFFREE_TPR6", "0") == "1" \
            and d["geometry"] != "TPR-6 trigonal prism":
        results += _enumerate_geometry(d, "trigonal_prism", "TPR-6 trigonal prism",
                                       lig_ref, lab_elem, spec, max_isomers)
    # Iter-31 (User 2026-05-28): CN4 dual SP-4 / T-4.  decompose picks ONE per metal via
    # _PREFERRED_CN4_GEOMETRY ('SQ' or 'TH') but real CN4 complexes can be either — esp.
    # Cu²⁺ where both T-4 (Cu(I)-like) and SP-4 (Cu(II) Jahn-Teller) exist.  Additive,
    # env-gated default OFF.  Adds the OPPOSITE of whichever the primary picked.
    if d.get("cn") == 4 and os.environ.get("DELFIN_FFFREE_DUAL_CN4", "0") == "1":
        if d["geometry"] == "SP-4 square planar":
            results += _enumerate_geometry(d, "tetrahedron", "T-4 tetrahedron",
                                           lig_ref, lab_elem, spec, max_isomers)
        elif d["geometry"] == "T-4 tetrahedron":
            results += _enumerate_geometry(d, "square_planar", "SP-4 square planar",
                                           lig_ref, lab_elem, spec, max_isomers)
    # Iter-32c: CN3 dual SP-3 trigonal-planar / T-3 T-shape (mirror of dual-CN4).
    # decompose picks ONE (d⁸ → T-shape, else SP-3); dual flag adds the other.
    if d.get("cn") == 3 and os.environ.get("DELFIN_FFFREE_DUAL_CN3", "0") == "1":
        if d["geometry"] == "SP-3 trigonal planar":
            results += _enumerate_geometry(d, "tshape", "T-3 T-shape",
                                           lig_ref, lab_elem, spec, max_isomers)
        elif d["geometry"] == "T-3 T-shape":
            results += _enumerate_geometry(d, "trigonal_planar", "SP-3 trigonal planar",
                                           lig_ref, lab_elem, spec, max_isomers)
    # generate-gate-floor: never return zero isomers if the decomposition succeeded
    return _coord_filter(results) or None


def _enumerate_geometry(d, geom_key, geom_name, lig_ref, lab_elem, spec, max_isomers):
    """Build all clean isomers of `d`'s ligand set on a SPECIFIC polyhedron (geom_name).
    Best-effort: skips isomers that fail to build / fail the self-gate, returns [] on any
    enumeration error.  Used to add SPY-5 alongside TBP-5 for CN5 (polytopal completeness)."""
    out: List[Tuple[str, str]] = []
    try:
        colorings = PIC.enumerate_isomers(geom_key, spec)
    except Exception:
        return out
    geom_tag = geom_name.split()[0]
    for k, coloring in enumerate(colorings[:max_isomers]):
        vertex_specs = [lig_ref[lab] for lab in coloring]
        _ex = _exempt_from_blocks(_heteroleptic_block_offsets(vertex_specs))   # #279/#281
        try:
            built = AC.assemble_heteroleptic_from_mols(d["metal"], geom_name, vertex_specs)
        except Exception:
            continue
        if built is None:
            continue
        syms, P = built
        syms, P = _maybe_relax(syms, P)
        if not _build_is_clean(syms, P, cn=d.get("cn"), geom=geom_name, exempt_pairs=_ex):
            continue
        vertex_elems = [lab_elem[lab] for lab in coloring]
        name = _classify_coloring(geom_key, vertex_elems)
        label = f"{name}-{geom_tag}-{k+1}" if name else f"{geom_tag}-{k+1}"
        out.append((_xyz(syms, P), label))
    return out


if __name__ == "__main__":
    for label, smi in [("cisplatin", "N[Pt](N)(Cl)Cl"),
                       ("[CoCl3(NH3)3]", "[NH3][Co]([NH3])([NH3])([Cl])([Cl])[Cl]"),
                       ("hexammineCo", "[NH3][Co]([NH3])([NH3])([NH3])([NH3])[NH3]")]:
        r = _fffree_isomers(smi)
        if r is None:
            print(f"{label:<16} -> None (legacy)")
        else:
            print(f"{label:<16} -> {len(r)} isomers: {[lab for _, lab in r]}")
            print("   first xyz head:", r[0][0].splitlines()[0], r[0][0].splitlines()[1])
