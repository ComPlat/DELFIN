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
from typing import Dict, List, Optional, Tuple
import numpy as np
from rdkit import Chem

from delfin.fffree import decompose as DEC
from delfin.fffree import polya_isomer_count as PIC
from delfin.fffree import assemble_complex as AC
from delfin.fffree import polyhedra as PLY
from delfin.fffree import ligand_relax as LR
import delfin._bond_decollapse as _bd


# F1 coverage forensik: post-decompose stage failure recorder. Same TSV sink
# as ``decompose._F1_REASON_LOG``; instruments the polya / assemble / self-gate
# failure paths.  Default OFF (byte-identical).
def _f1_backend_log(smi: str, reason: str, extra: str = "") -> None:
    _path = os.environ.get("DELFIN_FFFREE_DECOMPOSE_REASON_LOG", "").strip()
    if not _path:
        return
    try:
        with open(_path, "a") as _fh:
            _fh.write(f"{smi}\t{reason}\t{extra}\n")
    except Exception:
        pass


def _maybe_relax(syms, P):
    """#38: env-gated COD-loss torsional/rigid-body ligand relaxer (default OFF).
    Relieves van-der-Waals clashes by rotating distal sub-trees about rotatable bonds —
    coordination (metal + donors within the coord_geom detection sphere) frozen, rigid
    fragments preserved, multi-axis never-worse firewall.  Validated on smoke: net +11,
    0 severe (hanom -22%, inter-ligand -29%, h-clash -27%, coord_geom unchanged).
    Enable via DELFIN_FFFREE_LIGAND_RELAX=1.

    Phase G9 (User 2026-05-31 'fundamental universell'): under PURE_TRACK3, ALSO
    run the FF-free defect refiner (delfin.fffree.refine.refine) to fix
    collapsed heavy-heavy bonds + clashes. This is the universal solution to
    the 75% of selfgate-rejected cases that have collapsed bonds from rigid
    fitting. Metal + donors frozen, ligand periphery relaxed. Pure geometric.

    Phase G12 (User 2026-05-31, Mogul TOP-blindspot on fb1ae9a-PT3): aromatic
    rings (PPh3 / PR3-aryl etc.) buckle after rigid placement + refine, leading
    to C sp2 bonds stretched to 1.835 Å vs COD 1.397 Å (+31 %, total severity
    26 307 on 41/500 files).  After the defect refiner, snap each detected
    aromatic ring to its SVD best-fit plane (perpendicular projection,
    minimal-movement), with ring-H rigidly riding the parent and M-coordinated
    atoms frozen.  Universal geometric, FF-free, auto under PURE_TRACK3.
    """
    out_syms = list(syms)
    P_curr = np.asarray(P, dtype=float)
    # Phase G9: universal FF-free defect refiner (auto under PURE_TRACK3)
    _PT3_REFINE = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
    if _PT3_REFINE:
        try:
            from delfin.fffree import refine as RF
            mi = [i for i, s in enumerate(syms) if _bd._is_metal(s)]
            fixed = set(mi)
            for m in mi:
                for j in range(len(syms)):
                    if j != m and syms[j] != "H" and float(np.linalg.norm(P_curr[j] - P_curr[m])) \
                            < 1.45 * _bd._ideal_bond(syms[m], syms[j]):
                        fixed.add(j)
            P_curr = RF.refine(out_syms, P_curr, fixed_idx=fixed)
        except Exception:
            pass
    # Phase G12 REVERTED (User 2026-05-31, smoke 461 vs fb1ae9a equal-n):
    # Aromatic-ring planarity snap was wired here but REGRESSED 37 severe axes
    # (net -59, pi_planar +115 %, hapto_geom +214 %, stereo +181 %) -- the SVD
    # plane projection over-corrected on chelate-attached rings and the
    # metal-coordinated-atom freeze did NOT propagate to the rest of the
    # rotational subtree, so the projection drift broke the linkage between
    # ring atoms and their non-ring substituents. Module + tests retained as
    # standby; an env-gated re-entry will require a per-ring "ring + bonded
    # substituent" rigid-body rotation (not a free atom-by-atom projection)
    # and additional validation before wire-on. See task #84 ("KEEP OFF").
    # Default OFF -- explicit no-op even under PURE_TRACK3.
    _AROMSNAP_ON = os.environ.get("DELFIN_FFFREE_AROMSNAP_FORCE", "0") == "1"
    if _AROMSNAP_ON:
        try:
            from delfin.fffree.aromatic_snap import snap_aromatic_rings
            out_syms, P_curr = snap_aromatic_rings(out_syms, P_curr)
        except Exception:
            pass
    if os.environ.get("DELFIN_FFFREE_LIGAND_RELAX", "0") != "1":
        return out_syms, P_curr
    try:
        mi = [i for i, s in enumerate(syms) if _bd._is_metal(s)]
        fixed = set(mi)
        for m in mi:
            for j in range(len(syms)):
                if j != m and syms[j] != "H" and float(np.linalg.norm(P_curr[j] - P_curr[m])) \
                        < 1.45 * _bd._ideal_bond(syms[m], syms[j]):
                    fixed.add(j)
        return list(syms), LR.relax(list(syms), P_curr, fixed)
    except Exception:
        return out_syms, P_curr

_GEOM_TO_POLYA = {
    "L-2 linear": "linear",                         # Phase G (2026-05-31): CN2 Cu(I)/Ag(I)/Au(I)
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
    # Task #44 / Mission A7 (2026-06-05): sandwich + piano-stool + half-sandwich.
    # Mission B1 (2026-06-05): use the EFFECTIVE per-site Pólya keys (cn=2/4/4
    # not 10/8/9), because the assembler treats a hapto ring as one
    # coordination site.  See ``sandwich_piano_polyhedra.effective_ref_vectors_sandwich``.
    "SANDWICH-10 bis-eta5-Cp": "sandwich_10_eff",
    "PIANO-STOOL-8 eta5-Cp+L3": "piano_stool_8_eff",
    "HALF-SANDWICH-9 eta6+L3": "half_sandwich_9_eff",
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


def _parse_xyz_block(xyz_block: str):
    """Inverse of :func:`_xyz` -- parse the headerless atom block back into
    ``(syms, P)``.  Used by the conformer-enumeration + RMSD-dedup
    post-processing path to round-trip XYZ strings without touching the
    upstream emit shape.  Returns ``(None, None)`` on any failure.
    """
    try:
        syms: List[str] = []
        P_rows: List[List[float]] = []
        for line in xyz_block.splitlines():
            parts = line.split()
            if len(parts) < 4:
                continue
            syms.append(parts[0])
            P_rows.append([float(parts[1]), float(parts[2]), float(parts[3])])
        if not syms:
            return None, None
        return syms, np.asarray(P_rows, dtype=float)
    except Exception:
        return None, None


def _conformer_enum_post(results, d):
    """Optional post-processing: enumerate single-bond rotamers and full
    Cremer-Pople ring puckers for each existing fffree result, then
    deduplicate the union by Kabsch-RMSD.

    Three independently env-gated knobs (default OFF -> byte-identical
    to HEAD):

      * ``DELFIN_FFFREE_ENUMERATE_ROTAMERS=1``  -> single-bond rotamer enum
      * ``DELFIN_FFFREE_RING_PUCKER_ALL=1``     -> Cremer-Pople enum on every
                                                   non-aromatic 5/6/7-ring
      * ``DELFIN_FFFREE_RMSD_DEDUP=1``          -> Kabsch-RMSD dedup of
                                                   the final union

    The two enumerators are ADDITIVE -- they extend the result list with
    new variants while keeping the originals.  When dedup is on, the
    extended list is collapsed by Butina clustering down to the essential
    conformer set (severity proxy = position in the input list, so the
    original isomers always win their cluster).
    """
    rot_on = os.environ.get("DELFIN_FFFREE_ENUMERATE_ROTAMERS", "0") == "1"
    pucker_on = os.environ.get("DELFIN_FFFREE_RING_PUCKER_ALL", "0") == "1"
    dedup_on = os.environ.get("DELFIN_FFFREE_RMSD_DEDUP", "0") == "1"
    if not (rot_on or pucker_on or dedup_on):
        return results
    if not results:
        return results

    # Hard caps to keep post-processing bounded per SMILES.  Without these,
    # a high-isomer-count input (chelate Pólya enumeration up to 64) crossed
    # with a flexible ligand (5 rotors -> 243 each) crossed with multi-ring
    # systems (Cremer-Pople -> 32 per isomer) can blow up to 100k+ frames
    # whose pairwise RMSD matrix becomes O(N^2 * M) = minutes per SMILES.
    # The caps are env-tunable.
    try:
        _ROT_PER = int(os.environ.get(
            "DELFIN_FFFREE_ENUMERATE_ROTAMERS_PER_ISOMER", "8"
        ))
    except (TypeError, ValueError):
        _ROT_PER = 8
    try:
        _PUCKER_PER = int(os.environ.get(
            "DELFIN_FFFREE_RING_PUCKER_PER_ISOMER", "8"
        ))
    except (TypeError, ValueError):
        _PUCKER_PER = 8
    try:
        _MAX_TOTAL = int(os.environ.get(
            "DELFIN_FFFREE_CONFORMER_POST_MAX_TOTAL", "200"
        ))
    except (TypeError, ValueError):
        _MAX_TOTAL = 200

    # ------------------------------------------------------------------
    # 1)  Parse each emitted XYZ block back into (syms, P).
    # ------------------------------------------------------------------
    parsed: List[Tuple[List[str], np.ndarray, str]] = []
    for xyz_block, label in results:
        syms_p, P_p = _parse_xyz_block(xyz_block)
        if syms_p is None:
            continue
        parsed.append((syms_p, P_p, label))
    if not parsed:
        return results

    # Determine universal CN / geometry context for the self-gate when
    # we have a decompose dict.  These are used to reject new variants
    # whose post-rotation/post-pucker geometry no longer passes
    # _build_is_clean (the same gate the original emit went through).
    _ctx_cn = None
    _ctx_geom = None
    _ctx_hapto = False
    if isinstance(d, dict):
        _ctx_cn = d.get("cn")
        _ctx_geom = d.get("geometry")
        _ctx_hapto = any(
            lg.get("is_hapto", False) for lg in d.get("ligands", []) or []
        )

    # ------------------------------------------------------------------
    # 2)  Build a tiny RDKit mol that we can use as graph topology for
    #     the rotamer enumerator (single bonds, ring info).  We re-use
    #     :func:`delfin._bond_decollapse._geometric_bonds` which gives a
    #     deterministic heavy-atom bond graph from coordinates -- the same
    #     graph the self-gate uses.
    # ------------------------------------------------------------------
    extended: List[Tuple[List[str], np.ndarray, str]] = list(parsed)

    if rot_on or pucker_on:
        try:
            from delfin.fffree.single_bond_rotamers import (
                enumerate_single_bond_rotamers,
            )
        except Exception:
            enumerate_single_bond_rotamers = None
        try:
            from delfin.fffree.ring_pucker_integration import (
                enumerate_mol_ring_conformers,
            )
        except Exception:
            enumerate_mol_ring_conformers = None

        for syms_p, P_p, label in parsed:
            mol_g = None
            try:
                from rdkit import Chem as _Chem
                rw = _Chem.RWMol()
                for s in syms_p:
                    rw.AddAtom(_Chem.Atom(s))
                # Add geometric bonds (heavy-heavy by distance) so the
                # rotamer + ring-pucker detection has a real bond graph.
                bonds = _bd._geometric_bonds(syms_p, P_p)
                added = set()
                for ii, jj in bonds:
                    key = (min(ii, jj), max(ii, jj))
                    if key in added:
                        continue
                    try:
                        rw.AddBond(int(ii), int(jj), _Chem.BondType.SINGLE)
                        added.add(key)
                    except Exception:
                        pass
                # Also add metal->donor bonds explicitly (geometric_bonds
                # skips them).  This lets the rotor / ring detector know
                # which atoms are coordinated -- metal-incident bonds are
                # excluded inside the rotamer enumerator, ring-pucker
                # detector keeps chelate rings.
                nm = len(syms_p)
                for im in range(nm):
                    if not _bd._is_metal(syms_p[im]):
                        continue
                    for jm in range(nm):
                        if jm == im or syms_p[jm] == "H":
                            continue
                        if _bd._is_metal(syms_p[jm]):
                            continue
                        try:
                            ideal = _bd._ideal_bond(syms_p[im], syms_p[jm])
                            dij = float(np.linalg.norm(P_p[jm] - P_p[im]))
                            if dij < 1.30 * ideal:
                                key = (min(im, jm), max(im, jm))
                                if key not in added:
                                    rw.AddBond(int(im), int(jm),
                                                _Chem.BondType.SINGLE)
                                    added.add(key)
                        except Exception:
                            pass
                mol_g = rw.GetMol()
                # Sanitize best-effort (ring detection needs SSSR).
                try:
                    _Chem.SanitizeMol(mol_g, catchErrors=True)
                except Exception:
                    pass
                # Explicit ring perception so GetRingInfo().AtomRings()
                # works on the geometric-bond mol even if SanitizeMol
                # bailed out partway.
                try:
                    _Chem.GetSSSR(mol_g)
                except Exception:
                    pass
                try:
                    _Chem.FastFindRings(mol_g)
                except Exception:
                    pass
            except Exception:
                mol_g = None
            if mol_g is None:
                continue

            if len(extended) >= _MAX_TOTAL:
                break

            # Subagent #129 follow-up: pre-polish inter-ligand clash gate.
            # When DELFIN_FFFREE_PRE_POLISH_CLASH_GATE (or auto-on under
            # ENUMERATE_ROTAMERS) is active, we count quick inter-ligand
            # clashes on the un-polished rotamer candidate BEFORE invoking
            # _maybe_relax / _build_is_clean.  Cheap up-front reject saves
            # the GRIP-polish compute for tractable cases.
            try:
                from delfin.fffree.inter_ligand_clash_gate import (
                    gate_enabled as _ilg_on,
                    count_inter_ligand_clashes_quick,
                    env_clash_threshold as _ilg_thr,
                    env_clash_vdw_fraction as _ilg_vdw,
                )
            except Exception:
                _ilg_on = None
                count_inter_ligand_clashes_quick = None
            _CLASH_GATE_ON = bool(
                _ilg_on(enumerate_rotamers_on=rot_on)
            ) if _ilg_on is not None else False

            # Identify ligand subgraphs ONCE per (parsed) isomer so the
            # clash gate doesn't re-walk the bond graph per rotamer.
            _ligand_subgraphs = None
            _metal_idx_local = None
            if _CLASH_GATE_ON:
                try:
                    from delfin.fffree.grip_ensemble import (
                        identify_ligand_subgraphs,
                    )
                    # Pick metal index from syms_p.
                    _metal_idx_local = next(
                        (i for i, s in enumerate(syms_p) if _bd._is_metal(s)),
                        None,
                    )
                    # Donors heuristic: all heavy near-metal atoms within
                    # 1.45 * ideal_bond -- same as _maybe_relax fixed-set.
                    _donor_idx: List[int] = []
                    if _metal_idx_local is not None:
                        for j in range(len(syms_p)):
                            if j == _metal_idx_local or syms_p[j] == "H":
                                continue
                            try:
                                dij = float(np.linalg.norm(
                                    P_p[j] - P_p[_metal_idx_local]
                                ))
                                if dij < 1.45 * _bd._ideal_bond(
                                    syms_p[_metal_idx_local], syms_p[j]
                                ):
                                    _donor_idx.append(j)
                            except Exception:
                                continue
                    _ligand_subgraphs = identify_ligand_subgraphs(
                        mol_g, _metal_idx_local or 0, _donor_idx,
                    )
                except Exception:
                    _ligand_subgraphs = None

            def _accept_variant(syms_v, P_v, new_lab):
                """Push the variant through _maybe_relax + _build_is_clean
                + the soft-polyhedron polish so a new rotamer/pucker frame
                receives the SAME post-processing the original isomer did.
                Reject the variant if the self-gate fails -- this prevents
                the dedup step from emitting unpolished broken geometries
                that inflate per-frame defect rates.

                Pre-polish clash gate (subagent #129 follow-up): when the
                env-flag is on, reject the variant BEFORE _maybe_relax when
                the quick inter-ligand clash count exceeds
                ``DELFIN_FFFREE_PRE_POLISH_CLASH_MAX`` (default 5).

                Returns the polished (syms, P, label) tuple on success or
                ``None`` if the variant is rejected.
                """
                # Pre-polish clash gate (cheap, runs BEFORE polish).
                if (_CLASH_GATE_ON
                    and _ligand_subgraphs is not None
                    and count_inter_ligand_clashes_quick is not None):
                    try:
                        _n_clash = count_inter_ligand_clashes_quick(
                            np.asarray(P_v, dtype=float),
                            list(syms_v),
                            _ligand_subgraphs,
                            metal_idx=_metal_idx_local,
                            threshold=_ilg_vdw(),
                        )
                    except Exception:
                        _n_clash = 0
                    if _n_clash > _ilg_thr():
                        return None
                try:
                    s_l, P_l = _maybe_relax(list(syms_v), np.asarray(P_v,
                                                                      dtype=float))
                except Exception:
                    return None
                if not np.all(np.isfinite(P_l)):
                    return None
                try:
                    if not _build_is_clean(
                        s_l, P_l, cn=_ctx_cn, geom=_ctx_geom,
                        has_hapto=_ctx_hapto,
                    ):
                        return None
                except Exception:
                    pass
                # G16 soft polyhedron polish (auto under PURE_TRACK3).
                if (
                    os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
                    or os.environ.get("DELFIN_FFFREE_SOFT_POLY", "0") == "1"
                ):
                    try:
                        s_l, P_l = _g16_soft_polyhedron_polish(
                            s_l, P_l, cn=_ctx_cn, geom=_ctx_geom,
                            has_hapto=_ctx_hapto,
                        )
                    except Exception:
                        pass
                return list(s_l), np.asarray(P_l, dtype=float), str(new_lab)

            # Rotamer enumeration -- skip the very first variant which is
            # the identity = original.  Capped at _ROT_PER per isomer.
            if rot_on and enumerate_single_bond_rotamers is not None:
                try:
                    rot_iter = enumerate_single_bond_rotamers(
                        mol_g, P_p, max_configs=_ROT_PER + 1
                    )
                    n_added = 0
                    for vk, (Pv, rot_lab) in enumerate(rot_iter):
                        if vk == 0:
                            continue
                        if n_added >= _ROT_PER:
                            break
                        if len(extended) >= _MAX_TOTAL:
                            break
                        if not np.all(np.isfinite(Pv)):
                            continue
                        polished = _accept_variant(
                            syms_p, Pv, f"{label}__{rot_lab}"
                        )
                        if polished is None:
                            continue
                        extended.append(polished)
                        n_added += 1
                except Exception:
                    pass

            # Cremer-Pople ring-pucker enumeration over EVERY non-aromatic
            # 5-7 ring.  The integrator already deduplicates internally,
            # caps the per-ring/total variants, and yields the original
            # geometry as its first frame (which we drop).
            if pucker_on and enumerate_mol_ring_conformers is not None:
                try:
                    pk_iter = enumerate_mol_ring_conformers(
                        mol_g, P_p,
                        max_per_ring=None,
                        skip_aromatic=True,
                        chelate_only=False,
                        rmsd_dedup_tol=0.15,
                        max_total_variants=_PUCKER_PER + 1,
                    )
                    n_added = 0
                    for vk, Pv in enumerate(pk_iter):
                        if vk == 0:
                            continue
                        if n_added >= _PUCKER_PER:
                            break
                        if len(extended) >= _MAX_TOTAL:
                            break
                        if not np.all(np.isfinite(Pv)):
                            continue
                        polished = _accept_variant(
                            syms_p, Pv, f"{label}__pucker{vk}"
                        )
                        if polished is None:
                            continue
                        extended.append(polished)
                        n_added += 1
                except Exception:
                    pass

    # ------------------------------------------------------------------
    # 3)  RMSD dedup the union.  Severity proxy: position in the extended
    #     list so the originals (added first) always win their cluster.
    # ------------------------------------------------------------------
    if dedup_on:
        try:
            from delfin.fffree.conformer_dedup import (
                dedup_by_rmsd,
                dedup_by_rmsd_preserve_originals,
                _env_pre_cluster_emit,
            )
            # Frames as (label, P, severity); shared syms (all atoms align
            # via the parse step above -- one isomer per input result).
            #
            # We dedup PER (atom-count + atom-symbol fingerprint) bucket so
            # that frames from different isomers with different atom
            # orderings don't get compared (Kabsch demands matching shape).
            #
            # Sub-bucket bookkeeping: ``n_orig_per_bucket`` keeps the number
            # of ORIGINAL parsed frames (i.e. those present in ``parsed``)
            # per bucket so the pre-cluster-emit path can preserve them.
            n_orig_total = len(parsed)
            buckets: Dict[Tuple[str, int], List[Tuple[str, np.ndarray, float]]] = {}
            n_orig_per_bucket: Dict[Tuple[str, int], int] = {}
            for k, (syms_p_e, P_p_e, label) in enumerate(extended):
                key = ("".join(syms_p_e), len(syms_p_e))
                buckets.setdefault(key, []).append((label, P_p_e, float(k)))
                if k < n_orig_total:
                    n_orig_per_bucket[key] = n_orig_per_bucket.get(key, 0) + 1
            collapsed: List[Tuple[List[str], np.ndarray, str]] = []
            _PRE_CLUSTER_EMIT = _env_pre_cluster_emit()
            # Iterate buckets in insertion order to keep output deterministic.
            for key, group in buckets.items():
                # syms reconstructed from the bucket key
                syms_b = list(key[0])  # may differ from list of element symbols
                # Use the actual syms from the first frame instead.
                first_idx = extended.index(
                    next(e for e in extended if "".join(e[0]) == key[0])
                )
                syms_b = list(extended[first_idx][0])
                if _PRE_CLUSTER_EMIT:
                    n_orig_in_bucket = n_orig_per_bucket.get(key, 0)
                    kept_frames = dedup_by_rmsd_preserve_originals(
                        group, n_originals=n_orig_in_bucket, syms=syms_b,
                    )
                else:
                    kept_frames = dedup_by_rmsd(group, syms=syms_b)
                for lab, P_k, _sev in kept_frames:
                    collapsed.append((list(syms_b), P_k, lab))
            extended = collapsed
        except Exception:
            pass

    # ------------------------------------------------------------------
    # 4)  Re-emit as the (xyz_string, label) shape the caller expects.
    # ------------------------------------------------------------------
    out: List[Tuple[str, str]] = []
    for syms_p, P_p, label in extended:
        try:
            out.append((_xyz(syms_p, P_p), str(label)))
        except Exception:
            continue
    return out or results


def _build_is_clean(syms, P, cn=None, geom=None, donors=None, has_hapto=False) -> bool:
    """Self-gate: reject a build that is destroyed — non-finite coordinates,
    any collapsed heavy-heavy bond, gross steric overlap, or OVER-COORDINATION
    (a non-coordinating atom intruding into the metal's first shell) — so fffree
    NEVER emits a structure worse than the legacy fallback would.  A failing build
    makes the whole complex fall back to the legacy pipeline (return None),
    guaranteeing fffree is never worse than UFF on its addressable subset.
    Universal, geometry-only.  Disable via DELFIN_FFFREE_SELFGATE=0.

    Phase G7 (2026-05-31): has_hapto=True relaxes CShM threshold for hapto
    structures (η-Cp/arene/allyl/diene placements have donors at ring-centroid
    distance, not σ-bond distance, which inflates CShM). Universal across η3-η8.

    ``donors`` (optional, global indices of the cn constructed donor atoms): when
    fffree KNOWS the coordinating atoms (it built them), the coordination-shape and
    over-coordination checks use the donors directly instead of the cn-closest-heavy
    heuristic.  This is essential for CHELATES, whose ring backbone legitimately sits
    ~2.4-2.9 A from the metal: the heuristic miscounts a ring carbon as a donor (wrong
    cshm) or as over-coordination, falsely rejecting correct chelate geometry."""
    if os.environ.get("DELFIN_FFFREE_SELFGATE", "1") == "0":
        return True
    P = np.asarray(P, dtype=float)
    if P.size == 0 or not np.all(np.isfinite(P)):
        return False
    syms = list(syms)
    bonds = _bd._geometric_bonds(syms, P)
    if _bd._count_collapsed(syms, P, bonds) > 0:        # any collapsed heavy bond
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
        # MISSION F1 (2026-06-05) coverage heal #3: under
        # DELFIN_FFFREE_F1_COVERAGE the self-gate CShM threshold is relaxed
        # significantly so high-cshm rough builds emit as native (not as
        # silent UFF fallback). Internals are downstream's job (GRIP,
        # UFF post-relax, post_grip_corrector run after).
        _F1_COV = (
            os.environ.get("DELFIN_FFFREE_F1_COVERAGE", "0") == "1"
        )
        _shmax = float(os.environ.get("DELFIN_FFFREE_SHAPE_MAX", "20.0"))
        if _F1_COV:
            _shmax = float(os.environ.get("DELFIN_FFFREE_SHAPE_MAX_F1", "200.0"))
        # Phase G7: hapto structures have donors at ring-centroid distance,
        # which inflates CShM measured at the σ-donor-shell. Relax to 60 Å.
        if has_hapto:
            _shmax = float(os.environ.get("DELFIN_FFFREE_SHAPE_MAX_HAPTO", "60.0"))
            if _F1_COV:
                _shmax = max(_shmax,
                             float(os.environ.get("DELFIN_FFFREE_SHAPE_MAX_HAPTO_F1", "400.0")))
        # High-CN (CN>=7) placement is less reliable than CN4-6, so a build that
        # only passes the loose CN4-6 threshold can still be worse than the legacy
        # fallback there.  A tighter high-CN shape gate (default 5.0) makes the
        # high-CN subset cleanly net-better than legacy (measured: net +11 vs +9 at
        # 20, regressions 10->6).  Deterministic, CN-keyed; env-tunable.
        if cn >= 7 and not has_hapto:
            _shmax = min(_shmax, float(os.environ.get("DELFIN_FFFREE_SHAPE_MAX_HIGHCN", "5.0")))
            if _F1_COV:
                # F1 mission: relax high-CN too -- internals downstream
                _shmax = max(_shmax,
                             float(os.environ.get("DELFIN_FFFREE_SHAPE_MAX_HIGHCN_F1", "50.0")))
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


def _g16_soft_polyhedron_polish(syms, P, cn=None, geom=None, donors=None,
                                  has_hapto=False):
    """Phase G16/G16b polish: soft polyhedron M-D relaxation AFTER the build
    has passed the clean check. Monodentate donors: per-atom subtree drag.
    Chelating donors: rigid-body translation of the chelate ligand body.
    Fallback to pre-snap on any failure or post-snap _build_is_clean fail.

    Targets the 28 iter-gate HEAL-FIRST axes (F3_bond +163 %, vdw +137 %,
    funcgrp +266 %) remaining after G15b eliminated the catastrophic-overlap
    class -- the construction-vs-relaxation trade-off documented in the
    Mogul -94.7 % paper headline.
    """
    try:
        from delfin.fffree.soft_polyhedron_relax import (
            relax_md_stretches, is_enabled as _soft_on
        )
        if not _soft_on():
            return list(syms), np.asarray(P, dtype=float)
        syms2, P2 = relax_md_stretches(list(syms), np.asarray(P, dtype=float))
        if not _build_is_clean(syms2, P2, cn=cn, geom=geom,
                                donors=donors, has_hapto=has_hapto):
            return list(syms), np.asarray(P, dtype=float)
        return syms2, P2
    except Exception:
        return list(syms), np.asarray(P, dtype=float)


def _g13_ring_scale_polish(syms, P, cn=None, geom=None, donors=None,
                            has_hapto=False):
    """Phase G13 polish: rigid uniform aromatic-ring scaling AFTER the build
    has passed the clean check.  Returns the polished (syms, P) if the polish
    also passes _build_is_clean, otherwise the input unchanged.  Falls back
    silently on any exception (the polish must never break a clean build).

    Distinct from the failed G12 plane-snap [[feedback_g12_aromsnap_failed_isomer_selection]]:
    G12 ran BEFORE the clean gate and perturbed isomer selection -> regressed
    pi_planar +115 % on smoke 461 via 17 A fallback isomers. G13 runs AFTER
    the gate; the worst case is a no-op revert to the already-clean input.
    """
    try:
        from delfin.fffree.aromatic_ring_scale import (
            scale_aromatic_rings, is_enabled as _ars_on
        )
        if not _ars_on():
            return list(syms), np.asarray(P, dtype=float)
        syms2, P2 = scale_aromatic_rings(list(syms), np.asarray(P, dtype=float))
        if not _build_is_clean(syms2, P2, cn=cn, geom=geom,
                                donors=donors, has_hapto=has_hapto):
            # post-polish broke the clean check -> revert
            return list(syms), np.asarray(P, dtype=float)
        return syms2, P2
    except Exception:
        return list(syms), np.asarray(P, dtype=float)


def _fffree_chelate_isomers(d, geom_key, max_isomers, smiles=""):
    """Build all distinct isomers of a chelate-containing complex (mixed bi-/
    monodentate) via the universal chelate-config enumerator + per-config
    geometric assembly.  Returns [(xyz, label), ...] or None."""
    ligands = d["ligands"]
    # Phase G5 (2026-05-31): allow ≥4-denticity ligands when PURE_TRACK3 enabled.
    # Previously blocked porphyrin (κ4)/salen (κ4)/DTPA (κ6). New path: enumerate
    # without rotation since high-denticity ligands have constrained orientation.
    _PT3_HIGH_DENT = (
        os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
        or os.environ.get("DELFIN_FFFREE_HIGH_DENTICITY", "0") == "1"
    )
    if any(lg["denticity"] >= 4 for lg in ligands) and not _PT3_HIGH_DENT:
        max_d = max(lg["denticity"] for lg in ligands)
        _f1_backend_log(smiles, "chelate_denticity_ge4_disabled",
                        f"max_denticity={max_d}")
        return None        # kappa>=4 (porphyrin/salen/DTPA) not yet supported -> legacy
    # Phase G8 (User 2026-05-31): fundamental + universal → remove the
    # class-specific aromatic-CN5 chelate pre-rejection. Let the universal
    # self-gate (_build_is_clean) decide based on geometry, not on chemical-
    # class hard-coding.  When PURE_TRACK3 enabled, skip this pre-filter and
    # rely on the universal mechanisms downstream.
    # Original gate kept for env-OFF byte-identity (existing OC-6/SP-4 untouched).
    _PT3_NO_AROMATIC_GATE = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
    if (not _PT3_NO_AROMATIC_GATE
        and geom_key in ("trigonal_bipyramid", "square_pyramid")):
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
        })
    try:
        configs = PIC.enumerate_chelate_configs(geom_key, specs)
    except Exception as _e:
        _f1_backend_log(smiles, "chelate_polya_exception", str(_e)[:80])
        return None
    if not configs:
        _f1_backend_log(smiles, "chelate_polya_no_configs", geom_key)
        return None
    geom_tag = d["geometry"].split()[0]
    results = []
    # MISSION F1 (2026-06-05) coverage heal #4a: hold one best-effort
    # build that PASSED assembly but FAILED the self-gate so we can emit
    # it as last-resort native when no clean isomer survives.
    _F1_COV = os.environ.get("DELFIN_FFFREE_F1_COVERAGE", "0") == "1"
    _f1_holdover = None
    for k, config in enumerate(configs[:max_isomers]):
        # per-config pruning (generate-gate-floor): a geometrically infeasible config
        # (e.g. a fac vertex-triple for a mer pincer) is SKIPPED, not bailed -- so one
        # bad isomer no longer drops the whole complex to legacy.  Clean isomers survive.
        try:
            built = AC.assemble_from_config(d["metal"], d["geometry"], config, ligands)
        except Exception:
            continue
        if built is None:
            continue
        syms, P, donors = built
        syms, P = _maybe_relax(syms, P)
        # Phase G7: propagate hapto status to self-gate (relaxes CShM threshold).
        _has_hapto = any(lg.get("is_hapto") for lg in ligands)
        if not _build_is_clean(syms, P, cn=d.get("cn"), geom=d.get("geometry"),
                               donors=donors, has_hapto=_has_hapto):   # donor-aware self-gate -> skip config
            if _F1_COV and _f1_holdover is None:
                _f1_holdover = (syms, P, k, geom_tag)
            continue
        # Phase G13b: aromatic ring scaling polish with HAPTO-PI EXPLICIT FREEZE
        # + SUBSTITUENT SUBTREE DRAG. Iter G13 (first variant) regressed
        # hapto_geom +162 % and frame_pct_element_list_change +98 %; the fixes
        # are (a) freeze any metal near-neighbour cluster of >=3 same-element
        # heavy atoms within 3.0 A (covers eta-3 / eta-5 / eta-6 hapto-pi
        # regardless of cov-sum mismatch), (b) translate every substituent
        # subtree atom by the same vector as its parent ring atom (preserves
        # ring_C / substituent_root attachment + all subtree-internal bonds).
        # Default OFF; PT3 auto-enable via env flag.
        # G13b REVERTED (HEAL-FIRST x2): substituent-subtree-drag broke
        # element-list (+99 %), decoord (+126 %), vdW (+57 %). Standby behind
        # explicit FORCE flag only. Need build-time fix, not post-hoc polish.
        _polish_on = (os.environ.get("DELFIN_FFFREE_RING_SCALE_FORCE", "0") == "1")
        if _polish_on:
            syms, P = _g13_ring_scale_polish(
                syms, P, cn=d.get("cn"), geom=d.get("geometry"),
                donors=donors, has_hapto=_has_hapto,
            )
        # G16/G16b: soft polyhedron M-D relaxation (chelate-aware rigid-body +
        # monodentate subtree drag). Auto-enabled under PT3.
        _soft_poly_on = (
            os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
            or os.environ.get("DELFIN_FFFREE_SOFT_POLY", "0") == "1"
        )
        if _soft_poly_on:
            syms, P = _g16_soft_polyhedron_polish(
                syms, P, cn=d.get("cn"), geom=d.get("geometry"),
                donors=donors, has_hapto=_has_hapto,
            )
        results.append((_xyz(syms, P), f"{geom_tag}-chelate-{k+1}"))
    # MISSION F1 (2026-06-05) coverage heal #4b: emit last-resort build
    # when nothing passed the self-gate (env-gated).
    if _F1_COV and not results and _f1_holdover is not None:
        _syms, _P, _k, _gt = _f1_holdover
        results.append((_xyz(_syms, _P), f"{_gt}-chelate-{_k+1}-F1lastresort"))
        _f1_backend_log(smiles, "f1_chelate_last_resort_emit", geom_key)
    # Hebel #101 (2026-06-02): optional rotamer + ring-pucker enumeration
    # post-processed by Kabsch-RMSD dedup.  All three env-flags default OFF
    # -> this call is a no-op byte-identical to HEAD when none are set.
    results = _conformer_enum_post(results, d)
    if not results:
        _f1_backend_log(smiles, "chelate_all_configs_failed", geom_key)
    return results or None


def _fffree_isomers(smiles: str, max_isomers: int = 50
                    ) -> Optional[List[Tuple[str, str]]]:
    # MOGUL-PRIMARY hook (2026-06-07, feat-mogul-primary-2026-06-07):
    # when ``DELFIN_FFFREE_MOGUL_PRIMARY=1`` is set we delegate to the
    # CCDC-bounds-matrix-driven construction described in the project
    # draft manuscript (sections "Closed-form coverage of the
    # configuration space" and "Empirical bounds from a 2.45M-key
    # crystallographic manifold").
    #
    # This path REPLACES the idealised polyhedron + Pyykkö + VSEPR
    # cascade with ONE constrained DG embed using the CCDC empirical
    # distributions for every 1-2 / 1-3 / 1-4 / M-D / D-D / ring pair
    # AND the seven dedicated TM categories (NHC carbene, hapto-η²/⁵/⁶,
    # μ-bridge, agostic, oxidative-addition) detected universally from
    # the molecular graph.
    #
    # On ANY failure inside the primary path we fall through to the
    # legacy V16 cascade so production stays safe.  Default OFF -> the
    # entire block is skipped and HEAD bytes are unchanged.
    if os.environ.get("DELFIN_FFFREE_MOGUL_PRIMARY", "0") == "1":
        try:
            from delfin.fffree import assemble_via_mogul as _MP
            _res = _MP.assemble_complex_mogul_primary(smiles)
            if _res is not None:
                _syms, _P = _res
                # The Pólya enumeration over isomers can be layered on
                # top later; for the architectural pivot we emit ONE
                # primary structure (the CCDC-bounds-matrix-driven
                # embed) which matches the draft's "single constrained
                # DG embed" contract.  Label tag matches the legacy
                # naming convention so the pool tooling treats this
                # entry as a normal fffree result.
                _label = "MOGUL-PRIMARY-1"
                return [(_xyz(_syms, _P), _label)]
        except Exception:
            # silent fall-through to legacy on primary failure
            pass

    # GRIP-Ensemble hook (Hebel #101, 2026-06-02).  When the env flag
    # ``DELFIN_FFFREE_GRIP_ENSEMBLE=1`` is set we delegate to the
    # deterministic-completeness enumerator (Pólya × Cremer-Pople ×
    # GRIP-polish × ranked).  Default OFF -> byte-identical to HEAD
    # (this entire block is skipped when the flag is unset).  On
    # ANY failure inside the ensemble we fall through to the legacy
    # path so production stays safe.
    if os.environ.get("DELFIN_FFFREE_GRIP_ENSEMBLE", "0") == "1":
        try:
            _eres = _ensemble_isomers(smiles, max_isomers)
            if _eres is not None and len(_eres) > 0:
                # Hebel #101 conformer post-processing also applies to
                # the GRIP-ensemble path -- decompose first so the helper
                # has access to the metal / cn metadata; fall back to a
                # minimal stub on failure.
                try:
                    _d_ens = DEC.decompose(smiles)
                except Exception:
                    _d_ens = None
                _eres = _conformer_enum_post(_eres, _d_ens or {})
                return _eres
        except Exception:
            pass
        # silent fall-through to legacy on ensemble failure

    d = DEC.decompose(smiles)
    if d is None:
        return None
    geom_key = _GEOM_TO_POLYA.get(d["geometry"])
    if geom_key is None or geom_key not in PIC._GROUPS:
        _f1_backend_log(smiles, "geom_not_in_polya_groups", d.get("geometry", "?"))
        return None
    if d.get("has_chelate"):
        return _fffree_chelate_isomers(d, geom_key, max_isomers, smiles=smiles)
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
    except Exception as _e:
        _f1_backend_log(smiles, "polya_enumerate_exception", str(_e)[:80])
        return None
    if not colorings:
        _f1_backend_log(smiles, "polya_no_colorings", geom_key)
        return None
    results: List[Tuple[str, str]] = []
    # MISSION F1 (2026-06-05) coverage heal #4c: same last-resort mechanism
    # for non-chelate. Under F1_COVERAGE, exceptions / None / self-gate
    # failures CONTINUE to next coloring instead of bailing the whole SMILES.
    _F1_COV_NC = os.environ.get("DELFIN_FFFREE_F1_COVERAGE", "0") == "1"
    _f1_holdover_nc = None
    for k, coloring in enumerate(colorings[:max_isomers]):
        vertex_specs = [lig_ref[lab] for lab in coloring]
        try:
            built = AC.assemble_heteroleptic_from_mols(d["metal"], d["geometry"], vertex_specs)
        except Exception as _e:
            _f1_backend_log(smiles, "assemble_exception", str(_e)[:80])
            if _F1_COV_NC:
                continue
            return None
        if built is None:
            _f1_backend_log(smiles, "assemble_returned_none", d.get("geometry", "?"))
            if _F1_COV_NC:
                continue
            return None
        syms, P = built
        syms, P = _maybe_relax(syms, P)
        # Phase G7: propagate hapto status to self-gate
        _has_hapto = any(lg.get("is_hapto") for lg in d.get("ligands", []))
        if not _build_is_clean(syms, P, cn=d.get("cn"), geom=d.get("geometry"),
                                has_hapto=_has_hapto):   # self-gate: destroyed/over-coord/shape-outlier -> legacy
            _f1_backend_log(smiles, "self_gate_failed",
                            f"cn={d.get('cn')} geom={d.get('geometry')}")
            if _F1_COV_NC:
                if _f1_holdover_nc is None:
                    _f1_holdover_nc = (syms, P, k, coloring)
                continue
            return None
        # Phase G13b polish (see chelate path above for full rationale).
        # G13b REVERTED (HEAL-FIRST x2): substituent-subtree-drag broke
        # element-list (+99 %), decoord (+126 %), vdW (+57 %). Standby behind
        # explicit FORCE flag only. Need build-time fix, not post-hoc polish.
        _polish_on = (os.environ.get("DELFIN_FFFREE_RING_SCALE_FORCE", "0") == "1")
        if _polish_on:
            syms, P = _g13_ring_scale_polish(
                syms, P, cn=d.get("cn"), geom=d.get("geometry"),
                has_hapto=_has_hapto,
            )
        # G16/G16b soft polyhedron polish (non-chelate path)
        _soft_poly_on = (
            os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
            or os.environ.get("DELFIN_FFFREE_SOFT_POLY", "0") == "1"
        )
        if _soft_poly_on:
            syms, P = _g16_soft_polyhedron_polish(
                syms, P, cn=d.get("cn"), geom=d.get("geometry"),
                has_hapto=_has_hapto,
            )
        vertex_elems = [lab_elem[lab] for lab in coloring]
        name = _classify_coloring(geom_key, vertex_elems)
        geom_tag = d["geometry"].split()[0]
        label = f"{name}-{geom_tag}-{k+1}" if name else f"{geom_tag}-{k+1}"
        results.append((_xyz(syms, P), label))
    # MISSION F1 (2026-06-05) coverage heal #4d: last-resort emit for
    # non-chelate when no clean coloring made it through.
    if _F1_COV_NC and not results and _f1_holdover_nc is not None:
        _syms, _P, _k, _coloring = _f1_holdover_nc
        vertex_elems = [lab_elem[lab] for lab in _coloring]
        name = _classify_coloring(geom_key, vertex_elems)
        geom_tag = d["geometry"].split()[0]
        label = (f"{name}-{geom_tag}-{_k+1}-F1lastresort" if name
                 else f"{geom_tag}-{_k+1}-F1lastresort")
        results.append((_xyz(_syms, _P), label))
        _f1_backend_log(smiles, "f1_nonchelate_last_resort_emit", geom_key)
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
    # Hebel #101 (2026-06-02): optional rotamer + ring-pucker enumeration
    # post-processed by Kabsch-RMSD dedup.  All three env-flags default OFF
    # -> this call is a no-op byte-identical to HEAD when none are set.
    results = _conformer_enum_post(results, d)
    # generate-gate-floor: never return zero isomers if the decomposition succeeded
    return results or None


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
        try:
            built = AC.assemble_heteroleptic_from_mols(d["metal"], geom_name, vertex_specs)
        except Exception:
            continue
        if built is None:
            continue
        syms, P = built
        syms, P = _maybe_relax(syms, P)
        if not _build_is_clean(syms, P, cn=d.get("cn"), geom=geom_name):
            continue
        vertex_elems = [lab_elem[lab] for lab in coloring]
        name = _classify_coloring(geom_key, vertex_elems)
        label = f"{name}-{geom_tag}-{k+1}" if name else f"{geom_tag}-{k+1}"
        out.append((_xyz(syms, P), label))
    return out


# ---------------------------------------------------------------------------
# GRIP-Ensemble adapter (Hebel #101, 2026-06-02)
# ---------------------------------------------------------------------------
def _ensemble_isomers(smiles: str, max_isomers: int = 50
                       ) -> Optional[List[Tuple[str, str]]]:
    """Bridge :func:`grip_ensemble_enumerate` -> ``[(xyz, label), ...]``.

    Called only when ``DELFIN_FFFREE_GRIP_ENSEMBLE=1``.  Returns the
    top-K (or full) ranked ensemble as the same ``(xyz_string, label)``
    list shape :func:`_fffree_isomers` produces, so downstream
    smiles_converter handling is unchanged.

    Determinism
    -----------
    The ensemble enumerator iterates Pólya configs + Cremer-Pople states
    in fixed sorted order, with a deterministic score-based ranking and
    lex tiebreak by ``(isomer_id, conformer_id, label)``.  Two calls with
    the same SMILES + same env produce the same XYZ list byte-for-byte.

    Failure mode
    ------------
    Returns ``None`` on ANY enumerator failure (no candidates, decompose
    skip, exception); the caller (``_fffree_isomers``) then falls through
    to the legacy non-ensemble path so production stays safe.
    """
    try:
        from .grip_ensemble import (
            grip_ensemble_enumerate,
            ensemble_emit_full,
            ensemble_topk,
            DEFAULT_MAX_CONFORMERS_PER_ISOMER,
            DEFAULT_MAX_TOTAL,
        )
    except Exception:
        return None
    # Cremer-Pople knob: opt-in via env.  When the GRIP-Ensemble is active
    # (the gate that gets us into this function in the first place) we
    # default to 4 ring-pucker conformers per isomer so the enumeration
    # actually covers chair/boat/twist + envelope variants on common
    # 5/6-rings — matches the spec in
    # ``project_grip_ensemble_complete_conformer_coverage_2026_06_02``.
    # Set ``DELFIN_FFFREE_GRIP_ENSEMBLE_MAX_CONF=1`` to fall back to a
    # single conformer per isomer for forensic A/B.
    try:
        max_conf = int(os.environ.get(
            "DELFIN_FFFREE_GRIP_ENSEMBLE_MAX_CONF", "4"
        ).strip())
    except (TypeError, ValueError):
        max_conf = 4
    if max_conf < 1:
        max_conf = 1
    try:
        max_tot = int(os.environ.get(
            "DELFIN_FFFREE_GRIP_ENSEMBLE_MAX_TOTAL",
            str(DEFAULT_MAX_TOTAL),
        ).strip())
    except (TypeError, ValueError):
        max_tot = DEFAULT_MAX_TOTAL
    try:
        res = grip_ensemble_enumerate(
            smiles,
            max_isomers=int(max_isomers),
            max_conformers_per_isomer=int(max_conf),
            max_total=int(max_tot),
        )
    except Exception:
        return None
    if res is None or res.skip_reason or not res.candidates:
        return None
    # Pick top-K (or full ensemble when env DELFIN_FFFREE_GRIP_ENSEMBLE_FULL=1).
    if ensemble_emit_full():
        winners = res.candidates
    else:
        winners = res.top_k or res.candidates[: max(1, ensemble_topk())]
    out: List[Tuple[str, str]] = []
    for cand in winners:
        # Reuse the canonical _xyz formatter so emitted blocks match every
        # other archive byte-for-byte.
        try:
            xyz = _xyz(list(cand.syms), cand.P)
        except Exception:
            continue
        out.append((xyz, cand.label))
    return out or None


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
