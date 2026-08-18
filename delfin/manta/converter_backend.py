"""delfin.manta.converter_backend — adapter wiring the metal-FF-free foundation
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

from delfin.manta import decompose as DEC
from delfin.manta import polya_isomer_count as PIC
from delfin.manta import assemble_complex as AC
from delfin.manta import polyhedra as PLY
from delfin.manta import ligand_relax as LR
from delfin.manta import backbone_reembed as _BR
import delfin.manta._bond_decollapse as _bd


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


def _local_heavy_bonds(mol):
    """Local heavy-heavy bond index pairs of one ligand fragment."""
    out = []
    for b in mol.GetBonds():
        a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if mol.GetAtomWithIdx(a1).GetAtomicNum() > 1 and mol.GetAtomWithIdx(a2).GetAtomicNum() > 1:
            out.append((min(a1, a2), max(a1, a2)))
    return out


def _graph_bonds_from_blocks(block_mols_offsets):
    """GLOBAL heavy-heavy bonds the SMILES graph REQUIRES, from the same
    (ligand_mol, block_offset) layout _exempt_from_blocks uses.

    This is the missing half of the self-gate.  `_build_is_clean` perceives its bonds
    GEOMETRICALLY (`bonds = _bd._geometric_bonds(syms, P)`), so it can only judge bonds
    that ARE there -- it asks "is this contact too short?" and never "is a bond the
    molecule requires MISSING?".  A torn ligand is therefore invisible BY CONSTRUCTION:
    once two bonded atoms drift apart, perception simply stops reporting the bond and
    there is nothing left to fail.  Measured consequence: `smiles_topology` fires on
    816 of 941 champion systems (87 %) and `core_torn` on 699, while the construction
    self-gate reports the build clean.  Anchoring against the graph closes the only
    direction the gate was blind in."""
    out = set()
    for mol, off in block_mols_offsets:
        for (i, j) in _local_heavy_bonds(mol):
            gi, gj = off + int(i), off + int(j)
            out.add((min(gi, gj), max(gi, gj)))
    return out


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


# --- CN4 dual-geometry (Td <-> SP-4) completeness ----------------------------
# For CN4 the crystal can be tetrahedral (T-4, ~109.5 deg donor-M-donor) OR square
# planar (SP-4, 90 deg x4 / 180 deg x2).  decompose() picks exactly ONE per metal
# (_default_geometry: SP-4 for the d8 set, else T-4), so the manifold previously held
# only one of the two CN4 shapes -- the crystal's actual geometry was MISSING whenever
# it disagreed with the metal's default (eye: QAKTOO is Td but Au->SP-4; INICIR misses
# both Td and a clean SP-4).  The legacy DELFIN_FFFREE_DUAL_CN4 flag added the opposite
# geometry on the MONODENTATE path only; CN4 CHELATE complexes (the bulk of real CN4
# phosphine/phosphite cages, incl. QAKTOO/INICIR) had NO dual-geometry pass at all and
# always emitted just the single decompose-chosen shape.
#
# DELFIN_FFFREE_CN4_BOTH=1 (default OFF -> byte-identical) makes CN4 enumeration ALWAYS
# emit BOTH the T-4 tetrahedron AND the SP-4 square-planar isomers whenever CN==4 and the
# donor set permits both -- covering the monodentate AND the chelate paths -- so the
# manifold contains whichever geometry the crystal picks.  Strictly ADDITIVE (never
# removes the primary geometry); the opposite-geometry pass is best-effort and never
# bails the primary result.  Universal / graph-only, deterministic, never raises.
_CN4_GEOMS = ("T-4 tetrahedron", "SP-4 square planar")


def _cn4_both_enabled() -> bool:
    """CN4 dual-geometry completeness: emit BOTH Td and SP-4 for plausible CN4
    (default OFF -> byte-identical when unset)."""
    return os.environ.get("DELFIN_FFFREE_CN4_BOTH", "0") == "1"


def _cn4_opposite_geometry(geom_name: str) -> Optional[str]:
    """The OTHER CN4 polyhedron name (T-4 <-> SP-4), or None if `geom_name` is not a
    CN4 geometry.  Used to additively enumerate the partner shape."""
    if geom_name == "T-4 tetrahedron":
        return "SP-4 square planar"
    if geom_name == "SP-4 square planar":
        return "T-4 tetrahedron"
    return None


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


def _config_template_mol(metal, lig_groups, syms):
    """An RDKit mol in the FRAME'S OWN atom order: metal at 0, then AddHs(ligand)
    blocks in construction order -- the same layout _lig_groups_from_config derives
    and backbone_reembed already relies on.

    WHY THIS EXISTS.  Every xyz->mol bridge in the tree (_xyz_to_rdkit_conformer)
    demands the element symbol to match at EVERY index, and the mol parsed from the
    SMILES carries RDKit's own order -- so it never matches a frame of ours and every
    consumer silently gets nothing.  That is exactly how the first attempt at wiring
    ring pucker into the FF-free path measured 185 of 187 systems byte-identical.

    Built here rather than reused from _finish_config_frame's ``cm`` on purpose: that
    one is deliberately lossy (atoms recreated from atomic number alone, so charge,
    aromaticity and chirality are dropped, and M-D is a SINGLE bond), because its only
    consumer is a bond-graph-only sp2 flattener.  CombineMols keeps all of it, and the
    M-D bond is DATIVE -- a SINGLE bond would push a neutral 3-coordinate donor N to
    valence 4, sanitisation would stop at SANITIZE_PROPERTIES, RingInfo would never be
    initialised, and every ring consumer would read zero rings.

    Returns None unless the result matches ``syms`` symbol-for-symbol -- an aggregate
    count check would pass a compensating pair of errors."""
    from rdkit import Chem as _Chem
    try:
        m = _Chem.RWMol()
        m.AddAtom(_Chem.Atom(str(metal)))
        for g in lig_groups:
            gh = _Chem.AddHs(g["mol"])
            base = m.GetNumAtoms()
            m.InsertMol(gh)
            for d in g["donor_local"]:
                m.AddBond(0, base + int(d), _Chem.BondType.DATIVE)
        mm = m.GetMol()
        try:
            _Chem.SanitizeMol(mm, catchErrors=True)
        except Exception:
            pass
        try:
            _Chem.FastFindRings(mm)          # belt: RingInfo even if sanitize stopped early
        except Exception:
            pass
        # SAY WHY IT FAILED.  A bare `return None` here is a SILENT null lever, and that is
        # the exact failure mode that cost two wrong conclusions on 2026-08-02: the first
        # ring-pucker wiring returned 0 for every system and nothing recorded that it had.
        # Trace-gated, so the default path stays silent and byte-identical.
        if mm.GetNumAtoms() != len(syms):
            return _scope_no("PUCKER_TEMPLATE_NATOMS",
                             "tmpl=%d frame=%d" % (mm.GetNumAtoms(), len(syms)))
        for i in range(mm.GetNumAtoms()):
            if mm.GetAtomWithIdx(i).GetSymbol() != syms[i]:
                return _scope_no("PUCKER_TEMPLATE_ORDER",
                                 "i=%d tmpl=%s frame=%s" % (
                                     i, mm.GetAtomWithIdx(i).GetSymbol(), syms[i]))
        return mm
    except Exception as _e:
        return _scope_no("PUCKER_TEMPLATE_RAISED", type(_e).__name__)




def _append_ffree_ring_puckers(results, metal, lig_groups, base_syms, base_P, base_label,
                               cn=None, geom=None, donors=None, exempt_pairs=None,
                               graph_bonds=None, max_isomers=0, budget=48):
    """Append the Cremer-Pople ring-pucker SIBLINGS of one accepted FF-free frame.

    Named apart from smiles_converter._append_ring_puckers on purpose: that one drives
    the METAL-FREE organic pool off a sanitised ETKDG mol.  This one drives a metal
    complex off OUR frame, in OUR atom order, with the coordination sphere frozen.

    THE conformer lever, and the biggest by reach: over 1000 systems, 655 are torsionally
    RIGID (one torsional state), so for two thirds of the space the ring pucker IS the
    entire conformer manifold -- and the FF-free path has none of it (2.23 frames per
    system, ZERO with a conformer suffix, against 30.6 elsewhere).  ETKDG cannot supply
    it either: UFF does not cross the ~10 kcal/mol chair-boat barrier, so every ring stays
    in whatever basin the random seed happened to hit.

    THE ONE PLACE DELFIN_FFFREE_RING_PUCKER IS READ (default OFF -> byte-identical).
    An earlier attempt read it in smiles_converter at the FF-free return and measured 185
    of 187 systems byte-identical: the legacy emitter bridges XYZ to a mol through a
    symbol-per-index match, and the SMILES-parsed mol carries RDKit's order while our
    frames carry metal-at-0 plus AddHs(ligand) blocks.  They never coincide.

    Additive by construction: the base frame is never touched, every sibling must pass the
    SAME self-gate the base passed, and a failing pucker is dropped rather than allowed to
    replace anything.  ``donors`` are the GLOBAL indices the assembler itself returned, so
    the frozen set is exact -- the legacy emitter has to guess it back from metal neighbours.
    """
    if os.environ.get("DELFIN_FFFREE_RING_PUCKER", "0") != "1" or not lig_groups:
        return
    try:
        from delfin.manta import _ring_pucker as _rpuck
        from rdkit import Chem as _Chem
        import numpy as _np
    except Exception:
        return
    m = _config_template_mol(metal, lig_groups, base_syms)
    if m is None:
        return
    try:
        conf = _Chem.Conformer(m.GetNumAtoms())
        for i in range(m.GetNumAtoms()):
            conf.SetAtomPosition(i, [float(base_P[i][0]), float(base_P[i][1]),
                                     float(base_P[i][2])])
        m.RemoveAllConformers()
        m.AddConformer(conf, assignId=True)
        frozen = {0} | {int(x) for x in (donors or [])}
        # the PRIMARY's own scores -- the bar every sibling has to clear (see below)
        _dloc = sorted(int(x) for x in (donors or []))
        try:
            from delfin.manta import assemble_complex as _AC
            _base_bad = (bool(_AC._collapsed_heavy_bonds_strict(list(base_syms), base_P)),
                         float(_AC._beta_score(list(base_syms), base_P, _dloc)))
        except Exception:
            _AC, _base_bad = None, None
        try:
            _base_min = _min_nonbonded_heavy(base_syms, base_P)
        except Exception:
            _base_min = None
        # angle_skip = the METAL alone.  Its angles come from the polyhedron, not from
        # hybridisation: the VSEPR gate sees nh == 4 on a CN4 centre and demands 109.5 deg,
        # so a square-planar d8's two 180 deg trans pairs read as a 70.5 deg error that no
        # pucker caused and none can fix.  Without this, EVERY combination of EVERY SP-4
        # and T-3 complex is rejected -- which reads as "no reach" for entirely the wrong
        # reason.  Default off in _ring_pucker, so no existing caller changes.
        _out = _rpuck.generate(m, frozen=frozen, budget=int(budget), angle_skip={0})
    except Exception:
        return
    # ===== ZAEHLUNG DER VERWERFUNGSGRUENDE (17.08.2026) ==============================
    # GEMESSEN am 16.08. (`folds`, 965 Systeme): **1055 von 1766 Ringidentitaeten =
    # 59,7 % tragen ueber den GANZEN Manifold nur EINE Faltung** -- und das MIT diesem
    # laufenden Emitter (RING_PUCKER ist Champion-Flag #24, und der FF-freie Pfad ist
    # 98,8 % der Faelle, Feuerzensus 14.08.).  Der Bauer existiert also, er laeuft, und
    # die Luecke bleibt.
    #
    # Die Frage ist damit nicht "wie bauen wir die zweite Faltung", sondern **an welchem
    # der sechs Tore sie stirbt**.  Von aussen sieht "erzeugt und verworfen" genauso aus
    # wie "nie gebaut" -- derselbe Fehlschluss wie beim Feuerzensus am 14.08.
    #
    # Reine INSTRUMENTIERUNG: zaehlt und meldet einmal, aendert am Verhalten nichts.
    # Die Meldung haengt an DELFIN_FFFREE_PUCKER_TRACE (Vorgabe aus), damit sie in
    # Produktionslaeufen nicht mitlaeuft.
    _pk = {"erzeugt": len(_out or []), "cap": 0, "anzahl": 0, "reihenfolge": 0,
           "unclean": 0, "kollaps": 0, "beta": 0, "clash": 0, "akzeptiert": 0}
    for _px, _plab in (_out or []):
        if max_isomers and len(results) >= max_isomers:
            _pk["cap"] += 1
            break
        try:
            _lines = [ln.split() for ln in _px.splitlines() if ln.strip()]
            if len(_lines) != len(base_syms):
                _pk["anzahl"] += 1
                continue
            _ps = [t[0] for t in _lines]
            if _ps != list(base_syms):
                _pk["reihenfolge"] += 1
                continue                      # order must survive; never guess a mapping
            _pP = _np.array([[float(t[1]), float(t[2]), float(t[3])] for t in _lines])
            _ps, _pP = _maybe_relax(_ps, _pP)
            if not _build_is_clean(_ps, _pP, cn=cn, geom=geom, donors=donors,
                                   exempt_pairs=exempt_pairs, graph_bonds=graph_bonds):
                _pk["unclean"] += 1
                continue
            # NEVER-WORSE PER SIBLING.  _build_is_clean asks "is this buildable", which is a
            # LOWER bar than "is this good".  Measured 2026-08-02 (ffpuck2): the pass was
            # strictly additive -- the primary frame stayed BYTE-IDENTICAL on every system --
            # and it still failed the gate, on exactly two: CAZJEW pyramid_frame_regressed,
            # YAGQIG quality_agg.  Those are arithmetic: 5 frames instead of 1, so a per-frame
            # count or mean moves even though nothing that existed got worse.
            #
            # The metric is not what is wrong.  The goal is "EVERY frame without an anomaly",
            # so a sibling that carries a defect the primary does not have is a real cost, and
            # arguing with the aggregate would be exactly the kind of metric-gaming this
            # project forbids.  So the sibling has to clear the SAME bar as the frame it hangs
            # off: no new collapsed bond, and no worse out-of-plane.  Both predicates already
            # exist and are the ones the gate itself judges by.
            if _AC is not None and _base_bad is not None:
                try:
                    if (_AC._collapsed_heavy_bonds_strict(_ps, _pP)
                            and not _base_bad[0]):
                        _pk["kollaps"] += 1; continue   # introduces a collapse the primary lacks
                    if _AC._beta_score(_ps, _pP, _dloc) > _base_bad[1] + 1e-9:
                        _pk["beta"] += 1; continue      # flatter donors were the point; worse is not
                except Exception:
                    pass
            # ... and the same never-worse test _append_reembed already applies to ITS extra
            # frames: a sibling must not bring the closest non-bonded heavy contact in tighter
            # than the primary has it.  Measured (ffpuck3): the beta+collapse filter removed
            # CAZJEW from the blockers exactly as intended and left ONE system, YAGQIG, failing
            # on tier2 and quality_agg -- neither beta nor collapse, so a third quantity.  This
            # is the cheapest one the eye also judges by, and reusing the pattern beats
            # inventing a fourth.
            if _base_min is not None:
                try:
                    if not _interlig_clash_ok(_ps, _pP, _base_min):
                        _pk["clash"] += 1; continue
                except Exception:
                    pass
            _pk["akzeptiert"] += 1
            results.append((_xyz(_ps, _pP), f"{base_label}-{_plab}"))
        except Exception:
            continue
    # EINMAL melden, welches Tor die Faltungen kostet.  Ohne diese Zeile ist "erzeugt und
    # verworfen" von "nie gebaut" nicht zu unterscheiden -- genau der Fehlschluss, der am
    # 14.08. den Feuerzensus wertlos gemacht hat.  Vorgabe AUS.
    if os.environ.get("DELFIN_FFFREE_PUCKER_TRACE", "0") == "1" and _pk["erzeugt"]:
        try:
            import logging as _lg
            _lg.getLogger(__name__).warning(
                "[pucker-trace] %s: erzeugt=%d akzeptiert=%d | verworfen: anzahl=%d "
                "reihenfolge=%d unclean=%d kollaps=%d beta=%d clash=%d cap=%d",
                base_label, _pk["erzeugt"], _pk["akzeptiert"], _pk["anzahl"],
                _pk["reihenfolge"], _pk["unclean"], _pk["kollaps"], _pk["beta"],
                _pk["clash"], _pk["cap"])
        except Exception:
            pass


# ---------------------------------------------------------------------------
# THE COUNTABLE CONFORMER AXIS
# ---------------------------------------------------------------------------
# For ISOMERS this project has Polya: a theory set, hence a denominator, hence
# completeness as a CHECKABLE quantity.  For CONFORMERS it has nothing of the kind -- the
# builder embeds, keeps one, and nobody can say how many it should have had.  Measured
# consequence: the FF-free chelate path emits 2.23 frames per system with ZERO conformer
# suffixes, against 30.6 on everything else, and ccdc_pucker_realized is false for a
# quarter of the census.
#
# THE LAW THAT MAKES IT COUNTABLE.  A molecule's geometry is fixed by 1-2 (bonds,
# chemistry), 1-3 (angles, hybridisation) and 1-4 (torsions) -- and the torsions are the
# ONLY free part, and they are DISCRETE.  Their minima follow from the LOCAL environment
# of the bond, not from a uniform raster:
#
#     sp3-sp3              3    gauche+, anti, gauche-      (syn ~0 is a MAXIMUM)
#     sp3-sp2              2    the two eclipsing minima at the sp2 centre
#     sp2-sp2 conjugated   2    s-cis, s-trans
#     ring bond            0    NO degree of freedom -- ring closure, not torsion
#     terminal/symmetric   1    a methyl turns but changes nothing distinguishable
#
# Choosing a well therefore COMPLETES the specification: after it, nothing is left
# underdetermined.  A conformer stops being something one SAMPLES and becomes something one
# INDEXES -- and the index set is countable.  That is the whole point.
#
# WHY IT IS TRACTABLE, measured over 1000 systems (conformer_theory.py): 655 are torsionally
# RIGID (n_theory == 1) -- for two thirds of the space the ring pucker IS the entire
# conformer manifold, which is why that lever came first.  Of the rest the median is 1 and
# p90 is 36; only 20 systems explode, and for those the standing rule is to rank by ENERGY,
# never by RMSD.
#
# WHY IT IS BUILT AS SIBLINGS.  Every flag that ever landed in this project ADDS
# (D8_SQ_ADD "a PURELY ADDITIVE sibling ... the PRIMARY frame is untouched", CN6_OH_ADD,
# STEREOCENTER_ENUM, CN4_BOTH); everything measured on 2026-08-02 that CHOSE instead of
# added, died -- 17 A/Bs.  The reason turned out to be measurable: cap_LOST can only be paid
# by a system we ALREADY build perfectly (a system with broken_frac 1.0 has no capability to
# lose), so the recurring losers are the CLEANEST systems we have -- manifold_clean 44 % vs
# 19 % over 953 systems.  The only safe change to a perfect system is one that does not
# touch it.

_WELL_SP3_SP3 = (60.0, 180.0, 300.0)     # gauche+, anti, gauche-
_WELL_CONJ = (0.0, 180.0)                # s-cis / s-trans, and the sp2 eclipsing pair
_WELL_MAX_SIBLINGS = 12                  # bounded: the product is exponential by nature


def _sp2_like(atom) -> bool:
    """sp2 for the purpose of torsional wells: aromatic, or carrying a multiple bond."""
    try:
        if atom.GetIsAromatic():
            return True
        for b in atom.GetBonds():
            if b.GetBondTypeAsDouble() > 1.4:
                return True
    except Exception:
        pass
    return False


def _torsion_wells(mol, i, j):
    """The absolute dihedral minima of the i-j bond, from its LOCAL environment."""
    try:
        a, b = mol.GetAtomWithIdx(int(i)), mol.GetAtomWithIdx(int(j))
    except Exception:
        return ()
    s1, s2 = _sp2_like(a), _sp2_like(b)
    if s1 and s2:
        return _WELL_CONJ                 # conjugated: s-cis / s-trans
    if s1 or s2:
        return _WELL_CONJ                 # sp3-sp2: the two eclipsing minima at the sp2
    return _WELL_SP3_SP3                  # sp3-sp3


def _distal_side(mol, i, j):
    """Atom indices on j's side of the i-j bond (the half a rotation moves)."""
    seen, stack = {int(j)}, [int(j)]
    while stack:
        cur = stack.pop()
        for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
            k = nb.GetIdx()
            if k == int(i) or k in seen:
                continue
            seen.add(k)
            stack.append(k)
    return seen


def _append_ffree_torsion_wells(results, metal, lig_groups, base_syms, base_P, base_label,
                                cn=None, geom=None, donors=None, exempt_pairs=None,
                                graph_bonds=None, max_isomers=0):
    """Append one SIBLING per distinct torsional well vector of the accepted frame.

    THE ONE PLACE DELFIN_FFFREE_TORSION_WELLS IS READ (default OFF -> byte-identical).

    Only bonds whose MOVING half contains no frozen atom are turned -- the metal and every
    donor stay exactly put, so the coordination sphere is carried through untouched and only
    pendant backbone turns.  That is deliberately the same restriction torsion_relax applies,
    and it is why a rigid chelate contributes nothing here while a large flexible ligand --
    the 464-system LIGAND_TOO_LARGE class -- contributes most.

    Every sibling must clear the SAME bar as the frame it hangs off: the self-gate, no new
    collapsed bond, no worse beta, and no tighter closest contact.  Nothing that already
    exists is altered or dropped."""
    if os.environ.get("DELFIN_FFFREE_TORSION_WELLS", "0") != "1" or not lig_groups:
        return
    try:
        from delfin.manta import assemble_complex as _AC
        from rdkit import Chem as _Chem
        from rdkit.Chem import rdMolTransforms as _RT
        import numpy as _np
        import itertools as _it
    except Exception:
        return
    m = _config_template_mol(metal, lig_groups, base_syms)
    if m is None:
        return
    frozen = {0} | {int(x) for x in (donors or [])}
    _dloc = sorted(int(x) for x in (donors or []))
    try:
        base_bad = (bool(_AC._collapsed_heavy_bonds_strict(list(base_syms), base_P)),
                    float(_AC._beta_score(list(base_syms), base_P, _dloc)))
        base_min = _min_nonbonded_heavy(base_syms, base_P)
    except Exception:
        return
    # rotatable bonds whose moving half is free of metal AND donors
    dofs = []
    try:
        for b in m.GetBonds():
            if b.GetBondType() != _Chem.BondType.SINGLE or b.IsInRing():
                continue
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if m.GetAtomWithIdx(i).GetAtomicNum() == 1 or m.GetAtomWithIdx(j).GetAtomicNum() == 1:
                continue
            h1 = sum(1 for n in m.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() > 1)
            h2 = sum(1 for n in m.GetAtomWithIdx(j).GetNeighbors() if n.GetAtomicNum() > 1)
            if h1 < 2 or h2 < 2:
                continue                      # terminal spin: nothing distinguishable turns
            for (a, c) in ((i, j), (j, i)):
                if _distal_side(m, a, c) & frozen:
                    continue                  # would move the metal or a donor
                r1 = next((n.GetIdx() for n in m.GetAtomWithIdx(a).GetNeighbors()
                           if n.GetIdx() != c and n.GetAtomicNum() > 1), None)
                r2 = next((n.GetIdx() for n in m.GetAtomWithIdx(c).GetNeighbors()
                           if n.GetIdx() != a and n.GetAtomicNum() > 1), None)
                if r1 is None or r2 is None:
                    continue
                w = _torsion_wells(m, a, c)
                if len(w) > 1:
                    dofs.append((r1, a, c, r2, w))
                break
    except Exception:
        return
    if not dofs:
        return _scope_no("WELLS_NO_DOF", "nrot=0")
    dofs = dofs[:4]                           # bound the product; 3^4 = 81 before dedup
    try:
        conf0 = _Chem.Conformer(m.GetNumAtoms())
        for i in range(m.GetNumAtoms()):
            conf0.SetAtomPosition(i, [float(base_P[i][0]), float(base_P[i][1]),
                                      float(base_P[i][2])])
        m.RemoveAllConformers()
        m.AddConformer(conf0, assignId=True)
    except Exception:
        return
    n_added = 0
    for combo in _it.product(*[d[4] for d in dofs]):
        if n_added >= _WELL_MAX_SIBLINGS:
            break
        if max_isomers and len(results) >= max_isomers:
            break
        try:
            w = _Chem.Mol(m)
            c = w.GetConformer()
            for (r1, a, b2, r2, _wl), ang in zip(dofs, combo):
                _RT.SetDihedralDeg(c, int(r1), int(a), int(b2), int(r2), float(ang))
            _pP = _np.array(c.GetPositions(), float)
            _ps = list(base_syms)
            if _np.allclose(_pP, _np.asarray(base_P, float), atol=1e-6):
                continue                      # the base frame already sits in this well
            _ps, _pP = _maybe_relax(_ps, _pP)
            if not _build_is_clean(_ps, _pP, cn=cn, geom=geom, donors=donors,
                                   exempt_pairs=exempt_pairs, graph_bonds=graph_bonds):
                continue
            if _AC._collapsed_heavy_bonds_strict(_ps, _pP) and not base_bad[0]:
                continue
            if _AC._beta_score(_ps, _pP, _dloc) > base_bad[1] + 1e-9:
                continue
            if base_min is not None and not _interlig_clash_ok(_ps, _pP, base_min):
                continue
            results.append((_xyz(_ps, _pP), "%s-well%d" % (base_label, n_added + 1)))
            n_added += 1
        except Exception:
            continue
# vdW-level inter-ligand clash floor (Å) for the NEW-frame never-worse gate below.
# A re-embedded / re-seated conformer must not introduce a non-bonded heavy-heavy
# contact below this (a backbone folding into a NEIGHBOUR ligand collapses well
# inside the vdW shell long before the 0.60·Σcov gross-overlap floor _build_is_clean
# uses fires — measured: ACEQUY Fe(N(Dipp)(SiMe3))3 reembed frames at 1.98-2.01 Å
# inter-ligand C-C, base frame 2.07-2.38 Å).
#
# ⚠ WIEDERHERGESTELLT 17.08.2026 — DIESE ZWEI ZEILEN WAREN FUENFZEHN TAGE WEG.
# Commit a9b82598 (02.08., "Ring-Pucker an die Stelle, die die Atomreihenfolge des
# Frames kennt") loeschte sie und liess die zwei Benutzungen in `_interlig_clash_ok`
# stehen.  Seither warf die Funktion bei JEDEM Aufruf einen NameError, den das
# `except Exception: continue` in `_append_reembed` verschluckte -- und damit fiel
# JEDER Frame, der das Selbstgate bestanden hatte, lautlos heraus.
#
# GEMESSEN am 17.08. mit zwei Feuerspuren ueber `pool_ffonly` (187 Systeme):
#   :938 reembed gerufen ........ 367 Treffer / 183 Systeme
#   :951 Frame in der Schleife .. 2446 / 182
#   :952 Selbstgate gerufen ..... 2446 / 182
#   :953 Selbstgate verwirft .... 1156 / 99
#   :954 Selbstgate BESTANDEN ... 1290 / 146     <- 2446-1156 = 1290, exakt
#   :955 vdW verwirft ...........    0 / 0
#   :956 angehaengt .............    0 / 0       <- kein einziger, je
# 1290 Frames verschwanden zwischen 954 und 956, also IN dieser Funktion.
#
# WAS ES GEKOSTET HAT: `BACKBONE_REEMBED` hatte am 01.08. Reichweite 135/187 = 72 %
# (Byte-Vergleich der actsweep-Arme) und am 17.08. 0/24.  Der Mechanismus zielt auf
# `ccdc_backbone`, die ZWEITGROESSTE Fehlmasse (411 von 840).  Er war nie kaputt --
# er erzeugt weiterhin 2446 Frames auf 182 von 187 Systemen.  Es wurde nur nichts
# davon angenommen.
#
# ⚠ DIE LEHRE STEHT SCHON IM REGISTER, und sie hat hier fuenfzehn Tage gekostet:
# ein Fehler unter einem breiten `except` ist unsichtbar, bis jemand ZAEHLT.  Ein
# Zensus ueber die Reichweite haette es am 02.08. gefunden; es gab keinen Verlauf,
# gegen den ein Abfall haette auffallen koennen.  Genau dafuer gibt es seit heute
# `harness/reach_watch.py` und `results/reach_history.jsonl`.
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


_SP2_PLANAR_BAND = 3.0        # deg of angle-sum deficit a sibling may add; same width as the
                              # gate's own pyramid band (harness/loop.py:1255)


def _sp2_planarity_worst(syms, P) -> float:
    """Worst deviation from planarity (deg) over every THREE-COORDINATE heavy atom.

    THE QUANTITY THE GATE ACTUALLY READS.  pyramid_frame_regressed compares
    ``pyramid_over_worst`` -- the WORST frame's sp2 Walsh excess -- and it has NO allowance for
    added frames, deliberately: a manifold is only complete up to frames that are themselves
    realistic, so a conformer that pyramidalises an sp2 centre is a defect and not coverage.

    The sibling bars first used ``_beta_score`` for this and it did not work: measured on 187
    systems, adding that test to both ensemble paths left the blocked systems EXACTLY unchanged
    (LIBCEH, TIBMEX, TIQFAB, XUYXOE before and after).  The two are different quantities.
    ``_beta_score`` asks whether the METAL lies in a DONOR's plane, scores only donors that
    have a plane, and counts only the excess over the crystal band; the gate asks whether ANY
    sp2 centre in the frame -- backbone included -- has been bent out of its own plane.  A
    conformer can be flawless at both donors and still fold a backbone sp2.

    Geometry only, no mol and no hybridisation label: an atom with exactly three bonded
    neighbours is planar when its three bond angles sum to 360 deg, and the deficit from 360
    is the deviation (0 for planar, about 31.5 for a perfect tetrahedral centre).  Atom-
    specific by construction -- no functional group is named, and a centre nobody has
    classified reads the same as any other.  Metals are skipped, so the metal-donor bond is
    never one of the three.
    """
    P = np.asarray(P, float)
    n = len(syms)
    heavy = [i for i in range(n) if not _bd._is_metal(syms[i])]
    nbrs = {i: [] for i in heavy}
    for ai in range(len(heavy)):
        i = heavy[ai]
        for bj in range(ai + 1, len(heavy)):
            j = heavy[bj]
            if float(np.linalg.norm(P[i] - P[j])) <= 1.30 * _bd._ideal_bond(syms[i], syms[j]):
                nbrs[i].append(j)
                nbrs[j].append(i)
    worst = 0.0
    for i in heavy:
        v = nbrs[i]
        if len(v) != 3 or syms[i] == "H":
            continue
        u = []
        for j in v:
            d = P[j] - P[i]
            nd = float(np.linalg.norm(d))
            if nd < 1e-9:
                u = []
                break
            u.append(d / nd)
        if len(u) != 3:
            continue
        s = 0.0
        for a, b in ((0, 1), (0, 2), (1, 2)):
            s += float(np.degrees(np.arccos(max(-1.0, min(1.0, float(np.dot(u[a], u[b])))))))
        dev = 360.0 - s
        if dev > worst:
            worst = dev
    return worst


def _org_bond_worst(syms, P) -> float:
    """Worst RELATIVE deviation of a bonded heavy-heavy pair from its covalent ideal.

    WHY THIS AND NOT THE COLLAPSE TEST.  The collapse tests fire below 0.82 x ideal, which is a
    destroyed bond.  The eye's org_bond axis is graded and fires far earlier, and that is what
    actually blocked the beta sibling: on AVUNUC02 -- the ONE system standing between that
    lever and the gate -- the added frame carried n_root_defects 0 -> 2 with no collapse and no
    contact, and the rows say why:

        org_bond_worst_sev   2.22 -> 2.71        (worse)
        coord_angle_best_maxdev  28.5 -> 22.8    (better)
        arrangement_worst_dev    87.4 -> 85.1    (better)

    The beta-optimal pick is a DIFFERENT conformer out of the pool, so it brings its own
    internal bond lengths; it was chosen for being flat at the donor and paid for it with a
    strained backbone bond.  Flatness does not buy a stretched bond, so a sibling that is worse
    here is not added.  Symmetric and relative, so it reads the same for any element pair.
    """
    P = np.asarray(P, float)
    n = len(syms)
    worst = 0.0
    for i in range(n):
        if syms[i] == "H" or _bd._is_metal(syms[i]):
            continue
        for j in range(i + 1, n):
            if syms[j] == "H" or _bd._is_metal(syms[j]):
                continue
            ideal = _bd._ideal_bond(syms[i], syms[j])
            d = float(np.linalg.norm(P[i] - P[j]))
            if d > 1.30 * ideal or ideal <= 0:
                continue                       # not a bond
            dev = abs(d - ideal) / ideal
            if dev > worst:
                worst = dev
    return worst


_ORG_BOND_BAND = 0.02         # 2 % of the covalent ideal: deterministic wobble, not a defect


def _org_bond_ok(syms, P, base_worst) -> bool:
    """A sibling may not strain a ligand bond further than the frame it hangs off already does."""
    try:
        return _org_bond_worst(syms, P) <= float(base_worst) + _ORG_BOND_BAND
    except Exception:
        return False


def _sp2_planarity_ok(syms, P, base_worst) -> bool:
    """A sibling may not bend an sp2 centre further than the frame it hangs off already does."""
    try:
        return _sp2_planarity_worst(syms, P) <= float(base_worst) + _SP2_PLANAR_BAND
    except Exception:
        return False


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


def _decollapse_enabled() -> bool:
    """Post-seating heavy-bond decollapse (env DELFIN_FFFREE_CHELATE_DECOLLAPSE,
    default OFF -> byte-identical).  The chelate seating (per-donor radial rescale
    in _orient_chelate_to_vertices) can collapse a backbone bond NEAR the metal
    (diagnosed: 87% of severe collapses are chelate-dispatch, 77% within 3 A of the
    metal).  When ON, a config that would FAIL the collapse self-gate is first run
    through the safe firewall mover before being skipped."""
    return os.environ.get("DELFIN_FFFREE_CHELATE_DECOLLAPSE", "0") == "1"


def _maybe_decollapse(syms, P):
    """Pull collapsed heavy-heavy bonds apart via the proven firewall mover
    (delfin.manta._bond_decollapse.correct_xyz: metals + donors FROZEN, collapse must
    STRICTLY drop, and no measured axis — vdW clashes, F20 H-planarity, F3 bond
    distortion, M-D invariant — may worsen, else the relaxed frame is rejected and
    the input returned bit-exact).  The collapsed backbone atoms (C/H, neither metal
    nor donor) are free to move, so a ring/backbone bond squeezed into the metal's
    shell during seating is relaxed back to its ideal length while the coordination
    geometry is held fixed.  Returns the SAME (syms, P) objects when the flag is off
    OR when the mover makes no change (no collapse, or firewall rejected) -> the
    caller's identity check keeps the OFF path byte-identical and skips a pointless
    re-gate when nothing moved."""
    if not _decollapse_enabled():
        return syms, P
    try:
        from delfin.manta import _bond_decollapse as _BD
        xin = _xyz(syms, P)
        xout = _BD.correct_xyz(None, xin)
        if xout == xin:
            return syms, P                       # bit-exact: nothing moved
        s2, P2, _ = _BD._parse(xout)
        if len(s2) == len(syms):
            return s2, np.asarray(P2, dtype=float)
    except Exception:
        pass
    return syms, P


def _trilat_rescue(build_fn, *, cn=None, geom=None, exempt_pairs=None, graph_bonds=None):
    """LAST RUNG of the seating ladder: re-assemble with trilaterated donor targets.

    Reached only after the rigid build failed the self-gate AND the conformer ladder is
    exhausted, so whatever this returns replaces a DISCARD -- it can never displace a clean
    frame.  That makes it additive BY CONSTRUCTION: never-worse holds structurally instead of
    having to be re-measured on every pool.

    Why a rung and not the primary path (trilatAB2, 995 systems, 2026-08-02): as the primary
    path trilateration lost 12 systems that were ALL topo_correct before and gained 11 that
    ALL had no valid frame before -- not one borderline case in either direction.  It repairs
    what is broken and damages what is whole; the ladder is the shape that fits that.

    Returns (syms, P) or None.  None whenever anything is off, missing, throws, or still
    fails the gate -- the caller then proceeds exactly as it did before.
    """
    if not AC.trilat_rescue_enabled():
        return None
    try:
        with AC.trilaterate_rescue():
            _b = build_fn()
    except Exception:
        return None
    if not _b:
        return None
    # The chelate builder returns (syms, P, donors); the heteroleptic one (syms, P).  Gate the
    # rescue against the donors IT produced, not the ones the rejected build had.
    if len(_b) == 3:
        _s, _P, _don = _b
    else:
        _s, _P = _b
        _don = None
    try:
        _s, _P = _maybe_relax(_s, _P)
        if not _build_is_clean(_s, _P, cn=cn, geom=geom, donors=_don,
                               exempt_pairs=exempt_pairs, graph_bonds=graph_bonds):
            return None
    except Exception:
        return None
    return _s, _P


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


# Ideal heavy-heavy multiple/aromatic bond lengths (A) per element pair, LONGEST-BOND-ORDER FIRST:
# (triple, double, aromatic) where known.  Single source of truth -- the collapse self-gate reads it
# for its length-based exemption, and the hapto path (assemble_complex._collect_exempt) reads it to
# emit a LENGTH-GATED exemption instead of an unconditional one.
MULTIBOND_IDEALS = {("C", "C"): (1.20, 1.34, 1.39), ("C", "N"): (1.16, 1.28, 1.34),
                    ("C", "O"): (1.13, 1.21, 1.28), ("N", "N"): (1.10, 1.25),
                    ("N", "O"): (1.21, 1.24), ("C", "S"): (1.55, 1.60),
                    ("O", "O"): (1.21,), ("N", "S"): (1.54,), ("C", "P"): (1.66, 1.55)}


def multibond_ideal(e1, e2, order):
    """Ideal length (A) for a bond of `order` between elements e1/e2, or None if unknown.

    order >= 2.5 -> triple, >= 1.75 -> double, else aromatic (falling back to the double
    entry when the pair has no separate aromatic ideal).  Returning None lets the caller
    keep the historic UNCONDITIONAL exemption for pairs we have no ideal for -- never-worse
    by construction: an unknown pair behaves exactly as it does today."""
    tup = MULTIBOND_IDEALS.get((min(e1, e2), max(e1, e2)))
    if not tup:
        return None
    idx = 0 if order >= 2.5 else (1 if order >= 1.75 else 2)
    if idx >= len(tup):
        idx = len(tup) - 1
    return float(tup[idx])


def _build_is_clean(syms, P, cn=None, geom=None, donors=None, exempt_pairs=None,
                    graph_bonds=None, block_bounds=None) -> bool:
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
        # 14.08.2026: war ein NACKTES `return False`.  Sechs Verwerfungsgruende tragen einen
        # Namen und erscheinen im Verwurfszensus, dieser nicht -- er fiel damit stillschweigend
        # in die Restmenge und war von "kein Isomer enumeriert" nicht zu unterscheiden.
        # Das ist der teuerste Grund von allen, den man NICHT sehen will: eine nicht-endliche
        # Koordinate ist ein BAUFEHLER, kein Chemieurteil, und gehoert getrennt gezaehlt.
        return _gate_no("NONFINITE_COORDS")
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
        return _bd.COLLAPSE_FLOOR * _bd._ideal_bond(a, b)   # war Literal 0.82; siehe _bond_decollapse

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

    # LENGTH-BASED multibond exemption (DELFIN_FFFREE_MULTIBOND_LENGTH_EXEMPT, default OFF
    # -> byte-identical).  ~20% of the FF-free "collapse" self-gate rejections are FALSE
    # FLAGS: a genuine multiple/aromatic bond (a metal-carbonyl C≡O ~1.13, a nitrile C≡N
    # ~1.16, an imine C=N ~1.28, an aromatic C~C ~1.39) whose bond-ORDER was lost in the
    # constructed mol, so the order-based _is_exempt never fires and the 0.82×single-bond
    # floor flags it as collapsed (measured: ~4 pp of build-coverage lost this way).  This
    # check exempts a flagged bond whose LENGTH itself matches a known multiple/aromatic
    # ideal for its element pair (no bond-order needed) -- a real collapse (e.g. C-C at
    # 0.38 A) sits far below every multibond ideal and is still caught.  Tolerance reuses
    # DELFIN_FFFREE_MULTIBOND_TOL.
    _mb_len_exempt = os.environ.get("DELFIN_FFFREE_MULTIBOND_LENGTH_EXEMPT", "0") == "1"
    _MB_LEN = MULTIBOND_IDEALS          # hoisted to module level so the hapto path can share it

    def _len_is_multibond(a, b, d):
        key = (min(a, b), max(a, b))
        for ideal in _MB_LEN.get(key, ()):  # legit short multiple/aromatic bond length?
            if abs(d - ideal) <= _mb_tol:
                return True
        return False

    n_coll = 0
    for i, j in bonds:
        if _bd._is_metal(syms[i]) or _bd._is_metal(syms[j]):
            continue
        d_ij = float(np.linalg.norm(P[i] - P[j]))
        if d_ij >= _coll_floor(syms[i], syms[j]):
            continue
        if _is_exempt(i, j, d_ij):                    # genuine short multiple bond -> pass
            continue
        if _mb_len_exempt and _len_is_multibond(syms[i], syms[j], d_ij):
            continue                                  # length-matches a multibond ideal -> pass
        n_coll += 1
    if n_coll > 0:
        return _gate_no("COLLAPSED_BOND")
    # GRAPH ANCHOR (DELFIN_FFFREE_TORN_GATE, default OFF -> byte-identical).  Everything above
    # judges bonds that geometric perception FOUND -- "is this contact too short?".  Nothing asks
    # the opposite question: "is a bond the MOLECULE REQUIRES missing?".  A torn ligand is therefore
    # invisible by construction: once two bonded atoms drift apart, _geometric_bonds stops reporting
    # them and there is nothing left to fail.  Measured: smiles_topology fires on 816 of 941 champion
    # systems (87 %) and core_torn on 699, while this gate calls the same builds clean -- it is the
    # single largest blind direction we have.  The bound mirrors the eye's own torn criterion
    # (_isolated_reseat._frame_topology_valid: broken above 1.4x the covalent sum); metal bonds are
    # excluded because M-D distances legitimately vary far more than covalent radii predict.
    if graph_bonds and os.environ.get("DELFIN_FFFREE_TORN_GATE", "0") == "1":
        _torn_f = float(os.environ.get("DELFIN_FFFREE_TORN_FACTOR", "1.4"))
        for i, j in graph_bonds:
            if i >= len(syms) or j >= len(syms):
                continue
            if _bd._is_metal(syms[i]) or _bd._is_metal(syms[j]):
                continue
            if float(np.linalg.norm(P[i] - P[j])) > _torn_f * _bd._ideal_bond(syms[i], syms[j]):
                return _gate_no("TORN_BOND")                              # the graph requires this bond; it is torn
    # ⚠ 2026-08-07 VERENGT.  Die erste Fassung schlug bei JEDER wahrgenommenen Bindung an, die
    # der Graph nicht fuehrt -- und war damit zu grob: tpr6spur verlor 7 Isomere und riss
    # ccdc_isomer_lost, ccdc_arrangement_lost und ccdc_backbone_lost, also SCHLECHTER als ohne.
    # Grund: die geometrische Wahrnehmung sieht auch LIGANDINTERN Bindungen, die der
    # Blockgraph nicht auffuehrt (Ringschluesse, Perzeptionsrand), und die sind harmlos.
    # Die physikalische Ursache im Prisma ist enger: die Dreiecksflaechen stehen auf Deckung,
    # also kollidieren VERSCHIEDENE LIGANDEN miteinander.  Genau darauf wird jetzt geprueft --
    # eine wahrgenommene Bindung zwischen zwei verschiedenen Ligandbloecken kann es chemisch
    # nicht geben, denn jeder Block ist ein eigenes Molekuel.
    if (block_bounds and graph_bonds
            and os.environ.get("DELFIN_FFFREE_SPURIOUS_BOND", "0") == "1"):
        _req = {(min(i, j), max(i, j)) for i, j in graph_bonds}

        def _blk(a):
            for _b, (_s, _e) in enumerate(block_bounds):
                if _s <= a < _e:
                    return _b
            return -1
        for i, j in bonds:
            if i >= len(syms) or j >= len(syms):
                continue
            if _bd._is_metal(syms[i]) or _bd._is_metal(syms[j]):
                continue
            if syms[i] == "H" or syms[j] == "H":
                continue
            if (min(i, j), max(i, j)) in _req:
                continue
            _bi, _bj = _blk(i), _blk(j)
            if _bi >= 0 and _bj >= 0 and _bi != _bj:
                return _gate_no("SPURIOUS_INTERLIG_BOND")   # zwei getrennte Molekuele, verschmolzen
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
                return _gate_no("GROSS_OVERLAP")
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
                return _gate_no("DONOR_DECOORD")
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
                    return _gate_no("SHELL_INTRUDER")
        else:
            close = 0
            for j in range(n):
                if j == mi or syms[j] == "H":
                    continue
                cutoff = max(1.45 * _bd._ideal_bond(syms[mi], syms[j]), 2.7)
                if float(np.linalg.norm(P[j] - P[mi])) < cutoff:
                    close += 1
            if close > cn + 1:                          # +1 slack for borderline
                return _gate_no("OVERCOORD")
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
                    # 14.08.2026: war ein NACKTES `return False` -- der EINZIGE Verwurf des
                    # Selbstgates, der auf die Polyedergestalt zielt, und der einzige ohne
                    # Namen.  Genau der fehlt in jeder Verwurfsrangliste, an der entschieden
                    # wird, welcher Defekt als naechstes drankommt.
                    return _gate_no("SHAPE_CSHM")
            except Exception:
                pass
    return True


def _gate_no(reason):
    """Say WHICH self-gate criterion rejected a build, then reject it.

    Behaviour is unchanged: it returns False exactly as the bare `return False` did, and says
    nothing unless the trace is on -- it delegates that decision to _iso_trace so the flag is
    read in ONE place (the alternative was a second copy of the env lookup, which is how this
    codebase ended up with nine different metal predicates).

    WHY THE NAME MATTERS.  The self-gate discards one enumerated isomer in three (measured
    2026-08-01: 66 enumerated, 24 dropped by the gate, 0 build failures), and that gap IS the
    distance between the paper's "complete by construction" and the 73 % the FF-free builder
    realises.  But "the gate said no" is not a root: COLLAPSED_BOND, TORN_BOND, GROSS_OVERLAP,
    DONOR_DECOORD, SHELL_INTRUDER and OVERCOORD are six different defects with six different
    fixes, and picking one without the counts would be guessing.
    """
    _iso_trace("GATE_" + reason, -1, "-")
    return False


def _iso_trace(reason, k, geom_tag):
    """Say WHY an enumerated coordination isomer never became a frame.

    DELFIN_FFFREE_ISO_TRACE=1, default silent -> byte-identical.

    THE MEASUREMENT THIS EXISTS FOR.  Against Burnside-Polya theory the FF-free builder
    realises 73 % of the predicted isomers where it fires, while the legacy path realises
    97 % (measured 2026-08-01 over 187 systems: 522 vs 692 of 706).  The paper claims
    "complete by construction", so that 27 % is the gap between the claim and the code.
    The enumerator is NOT the hole -- enumerate_chelate_configs produces the configs and
    they are then dropped one by one, either because the build returns None or because the
    self-gate rejects the geometry.  Which of the two decides where the work goes:
    a failing BUILD is a constructor problem, a failing SELF-GATE is a SEATING problem --
    and the seating roots are already measured (metal out of the donor plane at 15.9 deg
    for tetradentates, bent sp centres, hydrogens left behind by the rescale).
    Counting is the cheapest way to tell them apart, and nothing counted before.
    """
    _ff_trace_write("[ISO_DROP] %s config=%d geom=%s" % (reason, k, geom_tag))


def _ff_trace_on():
    """The ONE place the FF-free trace flag is read (DELFIN_FFFREE_ISO_TRACE)."""
    _v = os.environ.get("DELFIN_FFFREE_ISO_TRACE", "0")
    return bool(_v) and _v != "0"


def _ff_trace_write(line):
    """Write one trace line -- to a FILE when the flag names a path, else stderr.

    ⚠ 2026-08-09, DER GRUND.  Dieser Trace wurde am 01.08. gebaut, um genau EINE Frage zu
    beantworten: warum faellt die Mehrheit der Systeme aus dem FF-freien Bauer heraus.  Der
    Docstring von _scope_no sagt es woertlich ("808 anonymous fall-throughs into a ranked
    work list").  Er hat nie eine Zeile geliefert -- weil er nach STDERR schreibt, und
    loop.py verwirft den stderr der Bau-Worker (er landet nur in results/debug_<rid>.log,
    und nur wenn --debug auf genau dieses System zeigt).

    Fuer den CLEANGATE-Trace ist derselbe Fehler am 04.08. schon einmal behoben worden --
    dort steht der Kommentar "A FILE, not stderr ... measured: 139 of 142 systems built, ZERO
    trace lines".  Die Korrektur wurde nur nicht auf den Nachbartrace uebertragen, und so hat
    das wichtigste Diagnoseinstrument des Projekts eine Woche lang ins Leere geschrieben.

    DELFIN_FFFREE_ISO_TRACE=1        wie bisher: stderr
    DELFIN_FFFREE_ISO_TRACE=<pfad>   O_APPEND, eine kurze Zeile pro Aufruf; das ist ueber
                                     parallele Worker hinweg atomar genug.
    DELFIN_FFFREE_TRACE_RID          optional, wird vorangestellt, damit aus der Verteilung
                                     ein POOL werden kann und nicht nur eine Rangliste.
    """
    _v = os.environ.get("DELFIN_FFFREE_ISO_TRACE", "0")
    if not _v or _v == "0":
        return
    _rid = os.environ.get("DELFIN_FFFREE_TRACE_RID", "")
    _out = ("%s %s" % (_rid, line)) if _rid else line
    try:
        if _v == "1":
            os.write(2, (_out + "\n").encode())
        else:
            with open(_v, "a") as _fh:
                _fh.write(_out + "\n")
    except Exception:
        pass


def _scope_no(reason, detail=""):
    """Say WHY the FF-free builder declined a whole system, then decline it.

    Returns None exactly as the bare `return None` did, and is silent unless the trace is on.

    THE MEASUREMENT THIS EXISTS FOR (User 2026-08-01: "ausrollen ... bis alle fffree laufen
    und besser als legacy sind").  On the shipped champion the FF-free builder fires for 187
    of 995 systems -- the other 808 fall through to legacy, and NOTHING recorded why.  That
    matters more than the quality gap does: measured on 934 identical systems with an
    identical denominator and an isomer counted only when a CLEAN frame realises it, FF-free
    reaches 50.5 % against legacy's 54.0 % and holds EIGHT TIMES as many defect-free manifolds
    (57.2 % vs 6.9 %).  So the distance to dropping legacy is not quality, it is SCOPE -- and
    a scope you cannot see cannot be rolled out class by class, largest first.

    Every `return None` in _fffree_isomers is a silent decline today; naming them turns 808
    anonymous fall-throughs into a ranked work list.
    """
    _ff_trace_write("[FFREE_SCOPE] %s%s" % (reason, (" " + detail) if detail else ""))
    return None


def _rescue_first_config_enabled() -> bool:
    """Darf das ERSTE Config eines Systems die letzte Rettungssprosse benutzen?

    (DELFIN_FFFREE_RESCUE_FIRST_CONFIG, Vorgabe AUS -> byte-identisch.)  Wirkt nur, wenn der
    Aufrufer zugleich `union=True` durchreicht -- die Begruendung steht an der Sperre selbst.
    """
    return os.environ.get("DELFIN_FFFREE_RESCUE_FIRST_CONFIG", "0") == "1"


def _fffree_chelate_isomers(d, geom_key, max_isomers, union: bool = False):
    """Build all distinct isomers of a chelate-containing complex (mixed bi-/
    monodentate) via the universal chelate-config enumerator + per-config
    geometric assembly.  Returns [(xyz, label), ...] or None."""
    ligands = d["ligands"]
    if (os.environ.get("DELFIN_FFFREE_KAPPA4", "0") != "1"
            and any(lg["denticity"] >= 4 for lg in ligands)):
        return None        # kappa>=4 (porphyrin/salen/DTPA): default legacy; KAPPA4=1 enables
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
    # the denominator: how many isomers the enumerator OFFERED, before any is dropped
    _iso_trace("ENUMERATED", len(configs[:max_isomers]), geom_tag)
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
        _gb = _graph_bonds_from_blocks(_config_block_offsets(config, ligands))
        # THE ENSEMBLE USED TO SHORT-CIRCUIT THE WHOLE PATH BELOW, AND THAT IS WHAT BROKE IT.
        #
        # It built its own canonical frame with assemble_from_config and emitted built[0] under
        # the PLAIN label -- while the single-frame path below reaches its frame through
        # _build_config_never_worse plus decollapse plus conformer seating.  Two different
        # builders under one label: switching the ensemble on REPLACED the primary frame rather
        # than adding to it, and it could even shift which config came first (CEJPIQ:
        # OC-6-chelate-2 -> OC-6-chelate-3-conf2).  Its own comment said the ensemble "only adds
        # conformer DIVERSITY on top"; the code did not do that.
        #
        # Measured on the 187-system pool: frame0 changed on 49 systems.  Fixing the sigma path
        # alone (further down this file) brought it to 22 -- and all 22 remaining were chelates,
        # i.e. exactly this branch.  ccdc_backbone_lost, which takes a MAX over frames and so
        # CANNOT fall for a truly additive lever, went 2 -> 0 with that first half of the fix.
        #
        # So the branch is gone.  The ordinary path below runs untouched, and the ensemble
        # conformers are appended as SIBLINGS next to the ring-pucker, torsion-well and beta
        # siblings, which is where every other conformer axis in this file already lives.  The
        # primary frame is then byte-identical to the flag-off build BY CONSTRUCTION, not by a
        # check that has to be trusted.
        built = _build_config_never_worse(d, config, ligands, geom_key)
        if built is None:
            _iso_trace("BUILD_NONE", k, geom_tag)
            continue
        syms, P, donors = built
        syms, P = _maybe_relax(syms, P)
        _clg = _lig_groups_from_config(config, ligands)
        if not _build_is_clean(syms, P, cn=d.get("cn"), geom=d.get("geometry"),
                               donors=donors, exempt_pairs=_ex, graph_bonds=_gb):   # donor-aware self-gate
            # Decollapse (DELFIN_FFFREE_CHELATE_DECOLLAPSE, default OFF, byte-id):
            # the seating can squeeze a backbone bond into the metal's shell; the
            # safe firewall mover pulls it back to ideal length (metals+donors frozen,
            # collapse must strictly drop, no axis worsens).  Re-gate after — accept
            # ONLY if the moved frame now passes the FULL self-gate (shape/overcoord
            # re-checked too), else fall through to seating/skip as before.
            _dc_s, _dc_P = _maybe_decollapse(syms, P)
            if (_dc_P is not P) and _build_is_clean(
                    _dc_s, _dc_P, cn=d.get("cn"), geom=d.get("geometry"),
                    donors=donors, exempt_pairs=_ex, graph_bonds=_gb):
                syms, P = _dc_s, _dc_P
            else:
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
                # Last rung -- but ONLY once this system already has an accepted config.
                # `results` non-empty means FF-free is definitely keeping the system, so a
                # rescued config can only ADD an isomer; it cannot flip the system away from
                # legacy.  That is the additivity the heteroleptic branch could not have,
                # and the difference is exactly what trilatresc measured.
                #
                # ===== DIE PRAEMISSE DIESER SPERRE FAELLT IM UNION-MODUS WEG (2026-08-09) =====
                # Der Satz oben nennt den Grund selbst: sie schuetzt davor, dass ein gerettetes
                # ERSTES Config das System von legacy WEGZIEHT.  Laeuft der Aufrufer im
                # Vereinigungsmodus, gibt es dieses Wegziehen nicht -- legacy baut zu Ende und
                # seine Frames bleiben stehen, unsere kommen daneben.  Die Sperre verteidigt
                # dann gegen eine Gefahr, die nicht mehr existiert.
                #
                # WAS SIE KOSTET, gemessen am 09.08. aus der Bauspur ueber 40 Systeme mit
                # Chelat-Enumeration: 309 Isomere enumeriert, 236 vom Bauer selbst verworfen
                # = 76 Prozent.  Und die Sprossen darueber stehen im Champion alle auf AUS
                # (DECOLLAPSE aus, CONFORMER_SEATING aus und global negativ mit cap_lost 78).
                # Fuer das erste Config eines Systems gibt es damit heute UEBERHAUPT KEINE
                # Rettung; scheitern alle Configs daran, verliert das System den FF-freien
                # Bauer ganz -- CHELATE_EMPTY, gemessen 35 von 136 Systemen.
                #
                # ⛔ Die Bedingung ist eine UND-Verknuepfung, nicht bloss Dokumentation: ohne
                # Union ist die Praemisse echt und die Sperre richtig.  Dieselbe Lehre wie am
                # 07.08. -- MULTIBOND_LENGTH_EXEMPT allein cap_lost 40, mit Union cap_lost 0
                # und gained 26.  Nicht der Hebel war der Schaden, das Entweder-Oder war es.
                if reseated is None and (results
                                         or (union and _rescue_first_config_enabled())):
                    reseated = _trilat_rescue(
                        lambda: _build_config_never_worse(d, config, ligands, geom_key),
                        cn=d.get("cn"), geom=d.get("geometry"),
                        exempt_pairs=_ex, graph_bonds=_gb)
                if reseated is None:
                    _iso_trace("SELFGATE", k, geom_tag)
                    continue                          # skip this config
                syms, P = reseated
        _lab = f"{geom_tag}-chelate-{k+1}"
        results.append((_xyz(syms, P), _lab))
        # Saat fuer das Kreuzprodukt Konformer x Faltung -- siehe die lange Begruendung
        # am Ende dieser Schleife.  Pro akzeptiertem Frame frisch, nie ueber Frames
        # hinweg gesammelt: der Faltungspass braucht die Atomreihenfolge SEINES Frames.
        _prod_on = os.environ.get("DELFIN_FFFREE_MANIFOLD_PRODUCT", "0") == "1"
        _prod_seeds = []
        # BETA AS A SIBLING, NOT AS A REPLACEMENT (DELFIN_FFFREE_BETA_SIBLING, default OFF).
        #
        # Measured 2026-08-02: choosing the flat-beta conformer INSTEAD of the clash-minimal
        # one lowers pyramidal_sp2 from 19.05 % to 12.70 % (and 37.07 -> 12.77 on the worst
        # systems) -- but costs 3 to 4 capabilities, because the frame it replaces was better
        # somewhere else.  Every flag that ever LANDED here adds instead of replacing, and
        # every one measured that night which replaced, died.  So: build the beta-optimal
        # frame as well and append it, leaving the primary byte-identical.  The eye's crystal
        # floors read the BEST frame over the manifold, so a second, flatter frame can only
        # help them -- and it cannot take anything away, because nothing was removed.
        #
        # A SIBLING MUST CLEAR THE SAME BAR AS THE FRAME IT HANGS OFF.  The self-gate alone is
        # not that bar: it asks "is this buildable", not "is this as good as what we already
        # have".  Measured 2026-08-03 on the 25-system A/B -- with the eye correction and the
        # round-trip axis both in place, EVERY term was zero except one, and that one was a
        # single system:
        #     AVUNUC02   frames 1 -> 2   n_root_defects 0 -> 2   broken_frac 0.0 -> 0.5
        # The added beta frame was itself BROKEN, carrying two independent root causes, and the
        # gate was right to refuse it (the term already forgives a rise up to the number of
        # frames ADDED; this one exceeded it).  The ring-pucker siblings pass the same gate
        # cleanly because they already carry these three checks -- so the beta sibling gets
        # them too: no NEW collapsed bond, no TIGHTER closest contact than the primary already
        # has.  ("No worse beta" is trivially true here: this frame IS the beta-optimal one.)
        if os.environ.get("DELFIN_FFFREE_BETA_SIBLING", "0") == "1":
            try:
                _bb = AC.assemble_from_config(d["metal"], d["geometry"], config, ligands,
                                              prefer_beta=True)
            except Exception:
                _bb = None
            if _bb is not None:
                _bs, _bP, _bd = _bb
                _bs, _bP = _maybe_relax(_bs, _bP)
                _bx = _xyz(_bs, _bP)
                if (_bx != _xyz(syms, P)                      # a DIFFERENT frame, or nothing
                        and (not max_isomers or len(results) < max_isomers)
                        and _build_is_clean(_bs, _bP, cn=d.get("cn"), geom=d.get("geometry"),
                                            donors=_bd, exempt_pairs=_ex, graph_bonds=_gb)):
                    # NO CONTACT TEST HERE, AND THAT IS MEASURED, NOT ASSUMED.  The bar first
                    # also demanded "no closer contact than the primary" -- but the primary is
                    # SELECTED as the clash-minimal frame, so that is a bar almost nothing can
                    # clear.  With it, betastrict measured affected=0 on 187 systems and the
                    # ensemble siblings measured 0 of 3 on a probe where the same code without
                    # it measured 3 of 3.  A sibling is judged on ITS OWN properties -- no new
                    # collapse -- and the clash question stays where it belongs, in the
                    # self-gate that every frame passes anyway.
                    _ok = True
                    if os.environ.get("DELFIN_FFFREE_BETA_SIBLING_STRICT", "0") == "1":
                        try:
                            if (AC._collapsed_heavy_bonds_strict(_bs, _bP)
                                    and not AC._collapsed_heavy_bonds_strict(syms, P)):
                                _ok = False                   # a collapse the primary does not have
                            # ... and the graded bond axis, which is what AVUNUC02 actually
                            # failed on: no collapse, no contact, a STRAINED ligand bond
                            # (org_bond_worst_sev 2.22 -> 2.71).  See _org_bond_worst.
                            if _ok and not _org_bond_ok(_bs, _bP, _org_bond_worst(syms, P)):
                                _ok = False
                            if _ok and not _sp2_planarity_ok(_bs, _bP,
                                                             _sp2_planarity_worst(syms, P)):
                                _ok = False                   # bends an sp2 the primary keeps flat
                        except Exception:
                            _ok = False                       # cannot prove it is as good -> do not add
                    if _ok:
                        results.append((_bx, f"{_lab}-beta"))
        # LONE-PAIR ORIENTATION AS A SIBLING (DELFIN_FFFREE_LP_SIBLING, default OFF).
        #
        # Seating two donors onto two vertices leaves exactly ONE rotational freedom -- about
        # the donor-donor axis, which both donors lie ON, so it moves neither of them.  That
        # one angle decides whether the metal sits in a conjugated donor's pi plane, i.e. beta,
        # our largest realism gap (18-19 % of frames against 0.98 % of crystals).  Today it is
        # decided by whatever the Kabsch SVD returns.
        #
        # Setting it in the SEATING was measured and is negative: reach 49 of 187, but
        # cap_LOST 2, valid 28 -> 27, mean +0.181 -- the turned backbone collides, the config
        # is dropped and the system hands over to legacy.  Guarding that with a collapse test
        # made it byte-identical instead (affected 0), which the source had already recorded
        # once as "why the first version changed almost nothing".
        #
        # So it is built the way things land here: an EXTRA frame.  The primary keeps the
        # orientation it has, a second frame carries the lone-pair-aligned one, and the eye's
        # crystal floors -- which read the BEST frame over the manifold -- can find it.  A
        # collision in the sibling costs nothing, because nothing was taken away.
        if os.environ.get("DELFIN_FFFREE_LP_SIBLING", "0") == "1":
            try:
                _lb = AC.assemble_from_config(d["metal"], d["geometry"], config, ligands,
                                              lp_orient=True)
            except Exception:
                _lb = None
            if _lb is not None:
                _ls, _lP, _ld = _lb
                _ls, _lP = _maybe_relax(_ls, _lP)
                _lx = _xyz(_ls, _lP)
                if (_lx != _xyz(syms, P)
                        and (not max_isomers or len(results) < max_isomers)
                        and _build_is_clean(_ls, _lP, cn=d.get("cn"), geom=d.get("geometry"),
                                            donors=_ld, exempt_pairs=_ex, graph_bonds=_gb)):
                    try:
                        _dl = sorted(int(x) for x in (_ld or []))
                        # it only earns its place if it is FLATTER; an equal one adds nothing
                        _better = (AC._beta_score(_ls, _lP, _dl)
                                   < AC._beta_score(list(syms), P, _dl) - 1e-9)
                        _nocoll = not (AC._collapsed_heavy_bonds_strict(_ls, _lP)
                                       and not AC._collapsed_heavy_bonds_strict(syms, P))
                    except Exception:
                        _better = _nocoll = False
                    if _better and _nocoll:
                        results.append((_lx, f"{_lab}-lp"))
        # SIGMA-ENSEMBLE CONFORMERS, now as siblings of the accepted frame rather than in
        # place of it (see the long note where the old short-circuit branch used to be).
        # Every one clears the same per-frame self-gate as before; the one that reproduces
        # the primary exactly is dropped instead of being emitted twice.
        if _sig_ens:
            try:
                _sens = AC.assemble_from_config(d["metal"], d["geometry"], config,
                                                ligands, n_frames=_n_chel)
            except Exception:
                _sens = None
            _pxyz = _xyz(syms, P)
            _pworst = None                            # primary's worst sp2 bend, on first need
            for _sfi, _sfr in enumerate(_sens or []):
                if max_isomers and len(results) >= max_isomers:
                    break
                _ss, _sP, _sd = _sfr
                _ss, _sP = _maybe_relax(_ss, _sP)
                _sxyz = _xyz(_ss, _sP)
                if _sxyz == _pxyz:
                    continue                          # this IS the primary frame
                if not _build_is_clean(_ss, _sP, cn=d.get("cn"), geom=d.get("geometry"),
                                       donors=_sd, exempt_pairs=_ex, graph_bonds=_gb):
                    continue                          # skip a bad frame, keep the clean ones
                # SAME BAR AS THE PRIMARY.  pyramid_frame_regressed reads the WORST frame's
                # sp2 excess with no allowance for added frames, and deliberately so -- a
                # manifold is only complete up to frames that are themselves realistic.  So a
                # conformer that pyramidalises a donor's sp2 centre, collapses a bond or sits
                # closer than the primary is not completeness, it is a defect with a label.
                # These are the checks the ring-pucker siblings already carry.
                try:
                    if (AC._collapsed_heavy_bonds_strict(_ss, _sP)
                            and not AC._collapsed_heavy_bonds_strict(syms, P)):
                        continue
                    _dloc = sorted(int(x) for x in (_sd or []))
                    if (AC._beta_score(_ss, _sP, _dloc)
                            > AC._beta_score(list(syms), P, _dloc) + 1e-9):
                        continue
                    if _pworst is None:
                        _pworst = _sp2_planarity_worst(syms, P)
                    if not _sp2_planarity_ok(_ss, _sP, _pworst):
                        continue                      # bends an sp2 the primary keeps flat
                except Exception:
                    continue                          # cannot prove equivalence -> do not add
                results.append((_sxyz, f"{_lab}-conf{_sfi+1}"))
                if _prod_on:
                    _prod_seeds.append((_ss, _sP, f"{_lab}-conf{_sfi+1}", _sd))
        # Ring-pucker siblings of this accepted frame (default OFF -> byte-identical).
        # Runs HERE, next to the frame it belongs to, because this is where the frame's
        # own atom order is known -- see the function for why the same call from the
        # smiles_converter return point measured 185 of 187 systems byte-identical.
        _append_ffree_ring_puckers(results, d["metal"], _clg, syms, P, _lab,
                                   cn=d.get("cn"), geom=d.get("geometry"),
                                   donors=donors, exempt_pairs=_ex, graph_bonds=_gb,
                                   max_isomers=max_isomers)
        # ... and the TORSIONAL half of the same axis.  The ring pucker covers the 655 of
        # 1000 systems that are torsionally rigid; this covers the rest.  Together they are
        # the countable conformer space: pucker basins x torsional wells, both enumerated,
        # both emitted as siblings, the primary untouched.
        _append_ffree_torsion_wells(results, d["metal"], _clg, syms, P, _lab,
                                    cn=d.get("cn"), geom=d.get("geometry"),
                                    donors=donors, exempt_pairs=_ex, graph_bonds=_gb,
                                    max_isomers=max_isomers)
        # Backbone re-embed (env DELFIN_FFFREE_BACKBONE_REEMBED, default OFF): add
        # core-preserving global-fold variants of this accepted chelate frame.
        _append_reembed(results, d["metal"], _clg,
                        syms, P, _lab, cn=d.get("cn"), geom=d.get("geometry"),
                        donors=donors)
        # ===== DER MANIFOLD WAR EINE SUMME, KEIN PRODUKT ==========================
        # Gemessen 2026-08-18 an 143904 Frame-Etiketten aus fuenf Archiven: die
        # Kombination "Konformer UND Ringfaltung" existiert NULL mal, obwohl beide
        # Achsen einzeln reichlich vertreten sind (6518 Faltungs- gegen 94191
        # Konformer-Etiketten).  107 Systeme bedienen beide Achsen -- keines baut ein
        # einziges Kreuzprodukt:
        #     CECWEM  1 primaer + 1 conf + 9 pucker = 11 gebaut,  Produkt waere 20
        #     QILQUX  4 + 17 + 10                   = 39 gebaut,  Produkt waere 231
        #     XIZTOS  2 + 15 + 6                    = 23 gebaut,  Produkt waere 119
        # Die Ursache steht drei Aufrufe weiter oben: die Geschwister-Erzeuger bekommen
        # alle `syms, P`, also das UNVERAENDERTE Primaerframe.  Sie schreiben in
        # `results`, lesen es aber nie.  Die Nachpaesse in smiles_converter sind dagegen
        # eine echte Kette (_ff = f(_ff)) -- deshalb komponieren nur die.
        #
        # DELFIN_FFFREE_MANIFOLD_PRODUCT (Vorgabe 0 -> byte-identisch) laesst die
        # Faltungs- und Torsionspaesse zusaetzlich ueber die AKZEPTIERTEN
        # Konformer-Geschwister laufen.  Aus 1 + 1 + 9 wird 1 + 1 + 9 + 9.
        #
        # ⚠ WARUM NUR DIE KONFORMER-SAAT UND NICHT AUCH -beta/-lp/-reembed:
        # der gemessene Befund ist "pucker x conf = 0"; -beta und -lp sind klein
        # (Einzelframes), und _append_reembed ist gemessen SCHMUTZIG (bbrefix: +11,9 pp
        # harte Frames).  Ein Produkt ueber eine schmutzige Achse vervielfacht den
        # Schmutz.  Die Erweiterung auf weitere Saaten gehoert hinter das
        # Zusicherungsprotokoll, nicht hierhin.
        #
        # ⚠ KEINE STILLE KAPPUNG.  Der Deckel meldet sich, wenn er bindet -- eine stille
        # Kappung liest sich hinterher als "mehr gab es nicht", und genau dieser Fehler
        # steckt schon in dofs[:4] und _WELL_MAX_SIBLINGS.
        if _prod_on and _prod_seeds:
            # Der Deckel ist ein BACKSTOP, kein Auswahlmittel: 24 liegt ueber der
            # groessten beobachteten Saatzahl (QILQUX, 17).  Er meldet sich ueber den
            # Trace-Kanal dieses Moduls -- also sichtbar, sobald
            # DELFIN_FFFREE_ISO_TRACE auf eine Datei zeigt, und sonst nicht.  Das ist
            # bewusst hier notiert: eine Kappung, von der man nur unter einem zweiten
            # Schalter erfaehrt, ist halb still, und wer die Zahlen liest, muss das
            # wissen.  Die eigentliche Schranke bleibt max_isomers.
            _pmax = max(1, int(os.environ.get(
                "DELFIN_FFFREE_MANIFOLD_PRODUCT_MAX", "24")))
            if len(_prod_seeds) > _pmax:
                _ff_trace_write(
                    "[PRODUCT_CAP] seeds=%d cap=%d dropped=%d lab=%s"
                    % (len(_prod_seeds), _pmax, len(_prod_seeds) - _pmax, _lab))
            for _ps, _pP, _plab, _pd in _prod_seeds[:_pmax]:
                if max_isomers and len(results) >= max_isomers:
                    break
                _append_ffree_ring_puckers(results, d["metal"], _clg, _ps, _pP, _plab,
                                           cn=d.get("cn"), geom=d.get("geometry"),
                                           donors=(_pd if _pd is not None else donors),
                                           exempt_pairs=_ex, graph_bonds=_gb,
                                           max_isomers=max_isomers)
                _append_ffree_torsion_wells(results, d["metal"], _clg, _ps, _pP, _plab,
                                            cn=d.get("cn"), geom=d.get("geometry"),
                                            donors=(_pd if _pd is not None else donors),
                                            exempt_pairs=_ex, graph_bonds=_gb,
                                            max_isomers=max_isomers)
    return results or None


def _shell_cshm(d, built):
    """Full coordination-shell CShM of a built frame vs the target geometry (lower =
    closer to the ideal polyhedron).  Returns +inf on failure / unclean build, so a
    failed alternative never wins the never-worse comparison."""
    try:
        syms, P, donors = built
        M = P[0]
        vecs = [P[i] - M for i in donors]
        return float(PLY.cshm(vecs, d["geometry"]))
    except Exception:
        return float("inf")


def _planar_polydentate_place_enabled() -> bool:
    """Gate for the in-plane (metal-COPLANAR) PLACEMENT of a rigid planar polydentate
    (DELFIN_FFFREE_PLANAR_POLYDENTATE_PLACE, default OFF => byte-identical)."""
    return os.environ.get("DELFIN_FFFREE_PLANAR_POLYDENTATE_PLACE", "0") == "1"


def _coplanar_not_worse(d, base, cand):
    """NEVER-WORSE accept test for a coplanar-placed candidate frame ``cand`` vs the
    historic ``base`` frame.  Accept only if the candidate (a) builds, (b) does NOT
    increase full coordination-shell CShM beyond a hair (1e-3 tolerance), and (c) does
    NOT worsen the minimum non-bonded heavy-heavy (inter-ligand) distance below the
    base's own (within the standard 5 % vdW tolerance / 2.0 A floor).  Geometry-only,
    deterministic; any failure => reject (keep base)."""
    if cand is None:
        return False
    if base is None:
        # no historic frame at all (the folded build itself failed the gate later) —
        # accept the coplanar candidate only on its own merit (finite CShM).
        return np.isfinite(_shell_cshm(d, cand))
    try:
        if _shell_cshm(d, cand) > _shell_cshm(d, base) + 1e-3:
            return False
        bsyms, bP, _ = base
        csyms, cP, _ = cand
        base_min = _min_nonbonded_heavy(bsyms, bP)
        return _interlig_clash_ok(csyms, cP, base_min)
    except Exception:
        return False


def _build_config_never_worse(d, config, ligands, geom_key):
    """Build one chelate-isomer config.  For a RIGID PLANAR tridentate on CN5 (TBP-5 /
    SPY-5) with DELFIN_FFFREE_PLANAR_MER_CN5=1, build the frame BOTH ways — with the
    meridional bite constraint (planar_bite=True) and without (planar_bite=False, the
    historic folded seating) — gate both, and KEEP the lower full-shell CShM.  The
    meridional build wins where it opens the bite toward the ideal (ANUCOE 12.7->2.1);
    the folded build survives where the forced bite would distort the shape WORSE
    (strict never-worse).

    RIGID PLANAR polydentate IN-PLANE PLACEMENT (DELFIN_FFFREE_PLANAR_POLYDENTATE_PLACE,
    default OFF): the folded metallacycle embed lifts the metal ~1.0 A OUT of the
    donors' own plane (a flat tridentate physically requires the metal IN that plane).
    When enabled and the config carries a rigid_planar dent>=3 ligand, ALSO build a
    metal-COPLANAR frame (planar_coplanar=True) and keep it over the historic frame iff
    it does not worsen CShM or the minimum inter-ligand distance (``_coplanar_not_worse``
    — strict never-worse).  Universal (any CN/geometry the chelate path reaches), graph
    /geometry-only.  Otherwise (flags off / not rigid-planar) this is the plain single
    build -> byte-identical."""
    metal, geom = d["metal"], d["geometry"]
    _has_rp = any(lg.get("rigid_planar") and lg.get("denticity") >= 3 for lg in ligands)
    _cn5_rp = (
        os.environ.get("DELFIN_FFFREE_PLANAR_MER_CN5", "0") == "1"
        and geom in ("TBP-5 trigonal bipyramid", "SPY-5 square pyramid")
        and any(lg.get("rigid_planar") and lg.get("denticity") == 3 for lg in ligands))
    _cop = _planar_polydentate_place_enabled() and _has_rp
    if not _cn5_rp and not _cop:
        try:
            return AC.assemble_from_config(metal, geom, config, ligands)
        except Exception:
            return None
    # base build: either the CN5 meridional-vs-folded never-worse pick, or the plain
    # historic single build (when only the coplanar flag is on).
    # The base build must be the HISTORIC folded build regardless of the env flag, so
    # explicitly pass planar_coplanar=False (else assemble_from_config would honour the
    # env flag and the "base" would itself be coplanar -> no never-worse reference).
    if _cn5_rp:
        best = None
        best_cshm = float("inf")
        for pb in (True, False):                  # meridional (bite-forced) vs folded
            try:
                b = AC.assemble_from_config(metal, geom, config, ligands,
                                            planar_bite=pb, planar_coplanar=False)
            except Exception:
                b = None
            if b is None:
                continue
            c = _shell_cshm(d, b)
            if c < best_cshm:
                best_cshm, best = c, b
    else:
        try:
            best = AC.assemble_from_config(metal, geom, config, ligands,
                                           planar_coplanar=False)
        except Exception:
            best = None
    if not _cop:
        return best
    # coplanar in-plane placement candidate (metal solved INTO the rigid donor plane)
    try:
        cand = AC.assemble_from_config(metal, geom, config, ligands, planar_coplanar=True)
    except Exception:
        cand = None
    if _coplanar_not_worse(d, best, cand):
        return cand
    return best


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
        from delfin.manta import _rotamer_diversity as RD
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
        # True transition / f-block metal Z (the genuine coordination centre). Heavy
        # MAIN-GROUP "metals" (Sb/Sn/Bi/Pb/Ga/In/Tl/Ge/As/Te) are DONORS when bound
        # to such a centre, NOT a second centre -- so they must (a) not be picked as
        # the centre and (b) be recognised as donors despite their longer M-E bond.
        _TF_Z = (set(range(21, 31)) | set(range(39, 49)) | set(range(57, 81))
                 | set(range(89, 104)))
        _Z = {"Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27,
              "Ni": 28, "Cu": 29, "Zn": 30, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42,
              "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "La": 57,
              "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78,
              "Au": 79, "Hg": 80}

        def _is_true_metal(s):
            return _Z.get(s, 0) in _TF_Z

        def _parse(xyz):
            # robust to BOTH standard XYZ (count + comment header) AND the
            # atom-lines-only format the public emitter uses (no header). Take any
            # line whose last 3 tokens are floats as an atom; skip count/comment.
            sy = []; P = []
            for ln in str(xyz).splitlines():
                p = ln.split()
                if len(p) < 4:
                    continue
                try:
                    xyz3 = [float(p[1]), float(p[2]), float(p[3])]
                except ValueError:
                    continue
                sy.append(p[0]); P.append(xyz3)
            if not sy:
                raise ValueError("no atoms")
            return sy, np.asarray(P, dtype=float)

        # donor-detection cutoff: 2.75 A for ordinary donors, but heavy main-group
        # donors sit at a longer M-E bond (our build ~3.3-3.6 A), so use a generous
        # element-aware cutoff = 1.45*(vdW_m + vdW_d)/2 floored at 2.75 -- recognises
        # a bonded heavy donor without admitting a flown-off one.
        _HEAVY_MG = {"Sb", "Sn", "Bi", "Pb", "Ga", "In", "Tl", "Ge", "As", "Te"}

        def _det_cut(ms, ds):
            # heavy main-group donors are built with a longer M-E bond (~3.3-3.6 A);
            # recognise them out to 4.0 A so a seated-but-long donor counts as a donor
            # (the reference-relative threshold then drops only the flown-off frames).
            if ds in _HEAVY_MG:
                return 4.0
            return max(2.75, 0.725 * (_bd._VDW.get(ms, 2.1) + _bd._VDW.get(ds, 1.7)))

        # reference = frame with the MOST donors within the (element-aware) cutoff,
        # using the TRUE metal centre.
        ref_mi = ref_don = ref_n = None; best = -1
        parsed = []
        for xyz, lab in results:
            try:
                sy, P = _parse(xyz); parsed.append((sy, P, xyz, lab))
            except Exception:
                parsed.append((None, None, xyz, lab)); continue
            mi = next((i for i in range(len(sy)) if _is_true_metal(sy[i])), None)
            if mi is None:
                mi = next((i for i in range(len(sy)) if _bd._is_metal(sy[i])), None)
            if mi is None:
                continue
            don = [k for k in range(len(sy)) if k != mi and sy[k] != "H"
                   and float(np.linalg.norm(P[mi] - P[k])) < _det_cut(sy[mi], sy[k])]
            if len(don) > best:
                best = len(don); ref_mi = mi; ref_don = don; ref_n = len(sy)
        if ref_mi is None or not ref_don:
            return results
        # per-donor threshold = REFERENCE (best-coordinated) M-D distance + slack
        # (reference-relative, so robust to a mis-calibrated absolute ideal: a heavy
        # donor seated at 3.6 A is the reference, decoordination is 3.6 + slack).
        # Bounded below by the donor-type ideal + slack so an already-stretched
        # reference frame cannot raise the bar arbitrarily.
        rsy = rP = None
        for s, p, _, _ in parsed:
            if s is not None and len(s) == ref_n:
                rsy, rP = s, p; break
        if rsy is None:
            return results
        thr = {}
        for k in ref_don:
            ref_d = float(np.linalg.norm(rP[ref_mi] - rP[k]))
            try:
                ideal = float(PLY.md_distance(rsy[ref_mi], rsy[k]))
            except Exception:
                ideal = ref_d
            thr[k] = max(ref_d, ideal) + slack
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


def _fffree_isomers(smiles: str, max_isomers: int = 50, union: bool = False
                    ) -> Optional[List[Tuple[str, str]]]:
    # `union`: der Aufrufer sagt, ob er unsere Frames NEBEN die von legacy stellt statt
    # statt ihrer.  Bewusst ein Argument und keine zweite Lesung des Vereinigungs-Schalters
    # -- dessen Lesestelle in smiles_converter traegt die Zusage, die EINZIGE zu sein, und
    # ein Repo mit neun verschiedenen Metall-Praedikaten hat sich diese Regel verdient.
    # Wirkung an genau EINER Stelle: der letzten Rettungssprosse im Chelat-Pfad.
    d = DEC.decompose(smiles)
    if d is None:
        return _scope_no("DECOMPOSE_NONE")
    # The coordination number is one number per complex and constant for the whole build,
    # and it is what actually moves a metal-donor distance: measured, Cd-N runs 2.283 /
    # 2.342 / 2.357 at CN 4 / 5 / 6.  Record it once here rather than threading it through
    # all twelve md_distance() call sites; read only behind DELFIN_FFREE_MD_MEASURED, and
    # a stale or missing value merely falls back to the coarser element-pair band.
    try:
        from delfin.manta import polyhedra as _PLY
        _PLY.set_current_cn(d.get("cn"))
    except Exception:
        pass
    geom_key = _GEOM_TO_POLYA.get(d["geometry"])
    if geom_key is None or geom_key not in PIC._GROUPS:
        return _scope_no("GEOM_NOT_IN_POLYA",
                         "geom=%s cn=%s key=%s" % (d.get("geometry"), d.get("cn"), geom_key))
    if d.get("has_eta"):
        return (_coord_filter(_fffree_hapto_isomers(d, max_isomers))
                or _scope_no("HAPTO_EMPTY", "cn=%s geom=%s" % (d.get("cn"), d.get("geometry"))))
    if d.get("has_chelate"):
        chel = _fffree_chelate_isomers(d, geom_key, max_isomers, union=union) or []
        # CN4 dual-geometry completeness (DELFIN_FFFREE_CN4_BOTH, default OFF ->
        # byte-identical): the chelate path builds only on the single decompose-chosen
        # CN4 shape, so the partner geometry (the one the crystal may actually have --
        # e.g. QAKTOO is Td but Au defaults to SP-4) was never emitted.  Additively
        # enumerate the OTHER CN4 polyhedron's chelate isomers on top.  Best-effort:
        # the opposite-geometry pass never bails the primary result.
        # NEVER-WORSE GUARD (2026-07-01): only ADD the opposite-CN4 chelate isomers when
        # the DEFAULT-geometry chelate path itself produced isomers (`chel` non-empty).
        # If the default chelate returned nothing (e.g. DIPZOV: T-4 chelate = 0 isomers),
        # this branch must return None so the caller falls through to the binding-mode /
        # non-chelate fallback that builds the CORRECT isomers (DIPZOV cis/see-saw = 0.14).
        # Without this guard the SP-4 `extra` alone masks the empty default, SUPPRESSES the
        # good fallback, and REGRESSES the system (DIPZOV 0.14 -> 0.56).  Purely additive:
        # where the default chelate succeeds, the opposite geometry is still added (YEGGUO
        # wins via the separate non-chelate path, unaffected).
        if (_cn4_both_enabled() and d.get("cn") == 4 and chel
                and _cn4_opposite_geometry(d["geometry"]) is not None):
            opp = _cn4_opposite_geometry(d["geometry"])
            opp_key = _GEOM_TO_POLYA.get(opp)
            if opp_key is not None and opp_key in PIC._GROUPS:
                d_opp = dict(d)
                d_opp["geometry"] = opp
                try:
                    extra = _fffree_chelate_isomers(d_opp, opp_key, max_isomers) or []
                except Exception:
                    extra = []
                chel = chel + extra
                if os.environ.get("DELFIN_CN4_DEBUG", "0") == "1":
                    os.write(2, ("[CN4_BOTH] default_chel=%d extra=%d total=%d geom=%s opp=%s\n"
                                 % (len(chel)-len(extra), len(extra), len(chel),
                                    d["geometry"], opp)).encode())
        _r = _coord_filter(chel or None)
        if os.environ.get("DELFIN_CN4_DEBUG", "0") == "1":
            os.write(2, ("[CN4_BOTH] after _coord_filter: %d (was %d)\n"
                         % (len(_r or []), len(chel))).encode())
        # A chelate complex that produced nothing falls through to legacy here.  That is the
        # single most consequential silent decline in the file -- chelates are the bulk of the
        # corpus -- so it is named like every other one.
        return _r or _scope_no("CHELATE_EMPTY", "cn=%s geom=%s nchel=%d"
                               % (d.get("cn"), d.get("geometry"), len(chel)))
    # ligand identity = canonical SMILES of each fragment; group by it
    lig_label, lig_ref, lab_elem = [], {}, {}
    for lg in d["ligands"]:
        try:
            lab = Chem.MolToSmiles(lg["mol"])
        except Exception:
            return _scope_no("LIGAND_SMILES_FAIL", "cn=%s" % d.get("cn"))
        lig_label.append(lab)
        lig_ref.setdefault(lab, (lg["mol"], lg["donor_local_idx"]))
        lab_elem[lab] = lg["donor_elem"]
    spec = dict(Counter(lig_label))
    try:
        colorings = PIC.enumerate_isomers(geom_key, spec)
    except Exception:
        return _scope_no("ENUM_ERROR", "geom=%s cn=%s" % (geom_key, d.get("cn")))
    if not colorings:
        return _scope_no("NO_COLORINGS", "geom=%s cn=%s nlig=%d"
                         % (geom_key, d.get("cn"), len(lig_label)))
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
        # THE ENSEMBLE SHORT-CIRCUIT USED TO SIT HERE, AND IT WAS NEVER ADDITIVE.
        #
        # It emitted ens[0] under the PLAIN base_label -- the label the single-frame path below
        # gives the frame it builds with assemble_heteroleptic_from_mols.  Two different
        # builders under one label, so switching the ensemble on REPLACED the primary frame
        # instead of adding to it; and where ens[0] failed the self-gate, the first SURVIVING
        # frame took the primary slot outright (QAYZUL: OC-6-1 -> OC-6-1-conf8, at an unchanged
        # frame count of one).  Measured on the 187-system pool: frame0 changed on 49 systems.
        #
        # Of the 16 systems blocking that A/B, 9 had a changed primary frame -- a real
        # regression on an existing frame, nothing an eye correction may excuse -- and only 7
        # had an untouched primary.  ccdc_backbone_lost proves it independently: that fraction
        # takes a MAX over frames (weddell/detectors/find_conformer_coverage.py:504-519), so a
        # genuinely additive lever CANNOT lower it, yet it dropped on 2 systems.
        #
        # A first fix kept the branch and rebuilt a canonical frame inside it.  That was the
        # wrong shape: it left 22 systems still changing, because a canonical frame assembled
        # here is not the frame the path below reaches through its self-gate and conformer
        # seating.  So the branch is gone entirely.  The ordinary path runs untouched and the
        # ensemble is appended as SIBLINGS after it -- the primary is byte-identical to the
        # flag-off build BY CONSTRUCTION rather than by a check that has to be trusted.
        try:
            built = AC.assemble_heteroleptic_from_mols(d["metal"], d["geometry"], vertex_specs)
        except Exception:
            return _scope_no("ASSEMBLE_EXC", "cn=%s geom=%s k=%d"
                             % (d.get("cn"), geom_tag, k))
        if built is None:
            return _scope_no("ASSEMBLE_NONE", "cn=%s geom=%s k=%d"
                             % (d.get("cn"), geom_tag, k))
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
            # NO TRILATERATION RESCUE HERE, and the reason is measured (trilatresc, 187
            # systems, 2026-08-02): this branch RETURNS, i.e. one bad coloring sends the
            # WHOLE system to legacy.  A rescue here does not fill a discard -- it CANCELS
            # THE HANDOVER, and legacy was building those systems better (FUBJAP 4 isomers
            # -> 1, TAJFAP 6 -> 2, REYFEI 3 -> 2).  My "additive by construction" claim was
            # false in exactly one nameable way: FF-free's failure is not a discard, it is a
            # handover.  Widening reach here trades coverage for quality, which is the whole
            # roll-out question and must never be decided by a seating fallback.
            _tried_seating = _seating_enabled() and _has_large_ligand(_lg)
            reseated = _seat_via_conformers(d["metal"], _lg, syms, P,
                                            cn=d.get("cn"), geom=d.get("geometry")) \
                if _tried_seating else None
            if reseated is None:
                return _scope_no("RESEAT_FAILED" if _tried_seating else "GATE_NO_RESEAT",
                                 "cn=%s geom=%s k=%d seating=%d"
                                 % (d.get("cn"), geom_tag, k, int(_seating_enabled())))
            syms, P = reseated
        label = base_label
        results.append((_xyz(syms, P), label))
        # Ensemble conformers as SIBLINGS of the accepted frame (see the note where the old
        # short-circuit branch used to be).  Same per-frame self-gate as before; the frame that
        # reproduces the primary exactly is dropped rather than emitted twice.
        if _cn2_ens or _sigma_ens:
            try:
                _ens = AC.assemble_heteroleptic_ensemble(
                    d["metal"], d["geometry"], vertex_specs,
                    n_frames=(_n_ens if _cn2_ens else _n_sigma))
            except Exception:
                _ens = None
            _pxyz = _xyz(syms, P)
            _hdons = None                   # donor indices, derived once on first need
            _hworst = None                  # primary's worst sp2 bend, on first need
            for _efi, _efr in enumerate(_ens or []):
                if max_isomers and len(results) >= max_isomers:
                    break
                _es, _eP = _maybe_relax(_efr[0], _efr[1])
                _exyz = _xyz(_es, _eP)
                if _exyz == _pxyz:
                    continue                # this IS the primary frame
                if not _build_is_clean(_es, _eP, cn=d.get("cn"), geom=d.get("geometry"),
                                       exempt_pairs=_ex):
                    continue                # skip a bad frame; keep the clean ones
                # Same bar as the primary, for the same reason as in the chelate path above.
                # The ensemble returns symbols and coordinates only, but the donors ARE known
                # here: the frame is [metal] + one AddHs block per vertex_spec in order, so
                # donor i sits at its block start plus its own local index -- the very layout
                # _heteroleptic_block_offsets already computes.  Measured why it matters: of
                # the four systems still blocking the additive ensemble on pyramid_frame_
                # regressed, THREE (LIBCEH, TIQFAB, XUYXOE) build no chelate frame at all and
                # come through here, where the beta test was missing.
                try:
                    if (AC._collapsed_heavy_bonds_strict(_es, _eP)
                            and not AC._collapsed_heavy_bonds_strict(syms, P)):
                        continue
                    if _hdons is None:
                        _hdons = []
                        _hp = 1
                        for _hfrag, _hdi in vertex_specs:
                            _hdons.append(_hp + int(_hdi))
                            _hp += Chem.AddHs(_hfrag).GetNumAtoms()
                        _hdons = sorted(_hdons)
                    if (AC._beta_score(_es, _eP, _hdons)
                            > AC._beta_score(list(syms), P, _hdons) + 1e-9):
                        continue            # a conformer that pyramidalises a donor is a defect
                    if _hworst is None:
                        _hworst = _sp2_planarity_worst(syms, P)
                    if not _sp2_planarity_ok(_es, _eP, _hworst):
                        continue            # ... and one that bends a BACKBONE sp2 likewise
                except Exception:
                    continue                # cannot prove equivalence -> do not add
                results.append((_exyz, f"{base_label}-conf{_efi+1}"))
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
    # GELTUNGSBEREICH (DELFIN_FFFREE_TPR6_EARLY_TM, default OFF -> byte-identisch).
    # Der Hebel oben feuert auf JEDEM CN6-System, dessen Primaergeometrie nicht schon TPR ist.
    # Seine eigene Begruendung ist viel enger: "early-TM Mo/W CN6 prefer TPR".  Beim spaeten
    # Uebergangsmetall ist das trigonale Prisma chemisch unrealistisch -- und genau dort reisst
    # der Bau.  GEMESSEN an tpr6cn6 (69 CN6-Systeme, FF-frei gebaut): 40 von 40 Systemen besser,
    # keines schlechter, mean -8,45, capability_lost 0 -- blockiert allein daran, dass die
    # hinzugefuegten Prismen-Frames zu 28 % hart sind gegen einen Boden von 4,7 %, mit
    # smiles_topology 17 gegen 1 und core_torn 3 gegen 0 als Defekttypen.
    #
    # Beide vorhandenen Gates scheiden aus: TORN_GATE sieht nur die FEHLENDE Bindung (zweimal
    # REACH 0/24), der Konsens-Gate TOPOLOGY_GATE verwirft die Prismen KOMPLETT (isomers_lost 9,
    # jedes betroffene System 3 -> 2 Isomere) -- er kann "anderes Isomer" nicht von "Artefakt"
    # unterscheiden.  Bleibt: das Prisma dort NICHT bauen, wo es chemisch nicht vorkommt.
    #
    # Die Menge ist elementbasiert und universell -- kein SMILES, kein Refcode, kein System.
    # d0-d2-Uebergangsmetalle der Gruppen 3-7, fuer die trigonal-prismatisches CN6 dokumentiert
    # ist (klassisch Mo/W-Dithiolene, dazu Nb/Ta/V/Zr/Hf/Re).  Eine Ladung steht am
    # Zerlegungs-Dict nicht zur Verfuegung, deshalb Element statt d-Zahl.
    _TPR6_EARLY_TM = frozenset((
        "Sc", "Y", "La",       # Gruppe 3
        "Ti", "Zr", "Hf",      # Gruppe 4
        "V", "Nb", "Ta",       # Gruppe 5
        "Cr", "Mo", "W",       # Gruppe 6
        "Mn", "Tc", "Re",      # Gruppe 7
    ))
    _tpr6_on = os.environ.get("DELFIN_FFFREE_TPR6", "0") == "1"
    if _tpr6_on and os.environ.get("DELFIN_FFFREE_TPR6_EARLY_TM", "0") == "1":
        _tpr6_on = str(d.get("metal") or "") in _TPR6_EARLY_TM
    if d.get("cn") == 6 and _tpr6_on \
            and d["geometry"] != "TPR-6 trigonal prism":
        results += _enumerate_geometry(d, "trigonal_prism", "TPR-6 trigonal prism",
                                       lig_ref, lab_elem, spec, max_isomers)
    # Iter-31 (User 2026-05-28): CN4 dual SP-4 / T-4.  decompose picks ONE per metal via
    # _PREFERRED_CN4_GEOMETRY ('SQ' or 'TH') but real CN4 complexes can be either — esp.
    # Cu²⁺ where both T-4 (Cu(I)-like) and SP-4 (Cu(II) Jahn-Teller) exist.  Additive,
    # env-gated default OFF.  Adds the OPPOSITE of whichever the primary picked.
    # DELFIN_FFFREE_CN4_BOTH (default OFF) generalises the same additive dual-geometry to
    # ALWAYS emit both CN4 shapes (and is also wired on the chelate path above) so the
    # crystal's geometry is never absent from the manifold; either flag triggers the add.
    if d.get("cn") == 4 and (os.environ.get("DELFIN_FFFREE_DUAL_CN4", "0") == "1"
                             or _cn4_both_enabled()):
        opp = _cn4_opposite_geometry(d["geometry"])
        if opp is not None:
            opp_key = _GEOM_TO_POLYA.get(opp)
            results += _enumerate_geometry(d, opp_key, opp,
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
    # Iter-32g (User 2026-06-19, eye-flagged ATENET): CN3 trigonal-PYRAMIDAL isomer.
    # decompose only ever emits SP-3 (planar, 120deg) or T-3 (T-shape, 90/180deg) for
    # CN3 — the third real CN3 geometry, the trigonal PYRAMID (3 donors on one
    # hemisphere, metal at the apex above the donor plane, ~107deg donor-M-donor; the
    # "vacant tetrahedron" / NH3 lone-pair shape), was never enumerated.  Additively
    # build it on the TPY-3 polyhedron alongside whatever the primary picked.  Env-gated
    # default OFF -> byte-identical when unset (the new geometry never appears).  Same
    # additive, best-effort, never-bails pattern as CN5 SPY-5 / CN6 TPR-6 / dual-CN4.
    if d.get("cn") == 3 and os.environ.get("DELFIN_FFFREE_CN3_PYRAMIDAL", "0") == "1" \
            and d["geometry"] != "TPY-3 trigonal pyramidal":
        results += _enumerate_geometry(d, "trigonal_pyramidal",
                                       "TPY-3 trigonal pyramidal",
                                       lig_ref, lab_elem, spec, max_isomers)
    # generate-gate-floor: never return zero isomers if the decomposition succeeded
    return (_coord_filter(results)
            or _scope_no("COORD_FILTER_EMPTY", "cn=%s geom=%s nres=%d"
                         % (d.get("cn"), d.get("geometry"), len(results))))


def _enumerate_geometry(d, geom_key, geom_name, lig_ref, lab_elem, spec, max_isomers):
    _topo_env_on = os.environ.get("DELFIN_FFFREE_TOPO_ENV", "0") == "1"
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
        # GRAPH-ANKER FUER DEN ADDITIVEN PFAD (2026-08-06).  Der Selbst-Gate bekam hier
        # `exempt_pairs`, aber NIE `graph_bonds` -- und sein Riss-Test haengt an genau dieser
        # Vorbedingung (`if graph_bonds and DELFIN_FFFREE_TORN_GATE`).  Auf allen additiv
        # enumerierten Frames war der Gate damit STRUKTURELL unerreichbar, obwohl sein eigener
        # Kommentar ihn "the single largest blind direction we have" nennt.
        # GEMESSEN: tpr6cn6 (TPR6 auf 69 CN6-Systemen) verbesserte 40 von 40 Systemen,
        # mean -8,45, cap_lost 0 -- und scheiterte allein daran, dass die hinzugefuegten
        # Prismen-Frames zu 28 % hart sind gegen einen Boden von 4,7 %.  Die Defekttypen
        # sagen warum: smiles_topology 17 gegen 1, core_torn 3 gegen 0.  Die Liganden REISSEN.
        # tpr6torn (TPR6+TORN_GATE gegen TPR6) meldete dann 0/24 Reichweite -- der Gate kam
        # gar nicht an.  Diese Zeile ist der Grund.
        # `_graph_bonds_from_blocks` nimmt laut eigenem Docstring genau dasselbe
        # (ligand_mol, block_offset)-Layout wie `_exempt_from_blocks` -- gleiche Quelle,
        # gleiche Zeile, keine neue Annahme.  Ohne TORN_GATE bleibt alles byte-identisch,
        # denn der Gate liest den Parameter nur unter seinem eigenen Env-Schalter.
        _gb = _graph_bonds_from_blocks(_heteroleptic_block_offsets(vertex_specs))
        # Blockgrenzen fuer die Interligand-Pruefung: jeder Ligand belegt einen
        # zusammenhaengenden Atombereich, Metall auf 0.  Gleiche Quelle wie _gb.
        _bb = []
        for _frag, _off in _heteroleptic_block_offsets(vertex_specs):
            _bb.append((_off, _off + Chem.AddHs(_frag).GetNumAtoms()))
        try:
            built = AC.assemble_heteroleptic_from_mols(d["metal"], geom_name, vertex_specs)
        except Exception:
            continue
        if built is None:
            continue
        syms, P = built
        syms, P = _maybe_relax(syms, P)
        if not _build_is_clean(syms, P, cn=d.get("cn"), geom=geom_name, exempt_pairs=_ex,
                               graph_bonds=_gb, block_bounds=_bb):
            continue
        # ===== TOPOLOGIE-FILTER AUF DEM ERGAENZTEN FRAME (2026-08-08) =====
        # ⛔ GEMESSEN UND WIDERLEGT -- NICHT WIEDER EINSCHALTEN OHNE NEUES KRITERIUM.
        #
        #   topoenvsolo (TOPO_ENV ALLEIN gegen den Champion, 40er-Sonde, CN6-Pool):
        #     ccdc_arrangement_lost 3   LIBNAO, LIBNES, URUTEH
        #     isomers_lost 4            HOQVAN, LIBNAO, LIBNES, URUTEH
        #     1 besser / 4 schlechter
        #   tpr6topoenv (mit TPR6 zusammen): exakt DIESELBEN Systeme, exakt dieselben Terme.
        #
        # Der ganze Schaden kommt vom Filter ALLEIN.  Meine erste Erklaerung -- er treffe ueber
        # die gemeinsame Funktion auch die CN4_BOTH-Ergaenzungen des Champions -- war FALSCH:
        # die Isolation zeigt, dass es nicht der Ort ist, sondern das KRITERIUM.  Der Vergleich
        # der Nachbar-Element-MENGE feuert auf Frames, die reale Arrangements und Isomere
        # tragen; die geometrische Perzeption sieht dort Nachbarschaften, die der Blockgraph
        # nicht auffuehrt, ohne dass etwas kaputt waere.
        #
        # Zum Vergleich, ohne diesen Filter:  tpr6final = BLOCKER 0, 17 von 17 Systemen besser,
        # historischer Boden bestanden.  Der Filter macht aus null Blockern zwei.
        #
        # Der Code bleibt stehen (Projektregel: nie loeschen), aber die Tarnung ist weg: wer
        # ihn einschaltet, weiss jetzt, dass er gemessen und schlechter ist.
        # DELFIN_FFFREE_TOPO_ENV, default OFF -> byte-identisch.
        #
        # WARUM HIER UND NICHT IM SELBST-GATE.  Ein Versuch, dasselbe in _build_is_clean zu
        # pruefen (SPURIOUS_BOND), traf AUCH die primaeren Frames und verlor Isomere:
        # tpr6spur2 meldete isomers_lost 4 und ccdc_arrangement_lost 3, also SCHLECHTER als
        # ohne.  Hier laeuft der Test ausschliesslich auf dem ERGAENZTEN Frame -- er kann
        # damit nur eine Ergaenzung verwerfen, nie ein bestehendes Isomer.  Das ist dieselbe
        # Bauform wie jede Landung dieses Projekts: ADD, never replace.
        #
        # WAS ER PRUEFT.  Der gemessene Rest von tpr6final war ausschliesslich das: 187
        # ergaenzte Frames, davon 53 hart, Defekttyp smiles_topology 17 und core_torn 3.
        # Das Auge definiert smiles_topology so: aus dem Molekuel steht fest, welche
        # SCHWER-NACHBARSCHAFT ein Atom haben MUSS; ein Frame-Atom, dessen wahrgenommene
        # Nachbarschaft davon abweicht, traegt einen Topologiebruch.
        #
        # graph_bonds liegt hier bereits in FRAME-Indizes vor (aus denselben Bloecken wie
        # exempt_pairs), also entfaellt das Zuordnungsproblem, an dem der Mechanismus in
        # smiles_converter:31880 einmal gestorben ist ("they never coincide, so it returned 0
        # every time").  Verglichen wird pro Atom die MENGE der Nachbar-Elemente; eine
        # Perzeptionsrandbedingung, die nur die ANZAHL gleicher Elemente aendert, bleibt
        # damit unauffaellig, waehrend ein abgeloester Substituent oder ein verschmolzener
        # Kontakt die Menge veraendert.
        if _topo_env_on and _gb:
            try:
                from collections import Counter as _te_C
                _req_nb = {}
                for _i, _j in _gb:
                    _req_nb.setdefault(_i, []).append(syms[_j])
                    _req_nb.setdefault(_j, []).append(syms[_i])
                _got_nb = {}
                for _i, _j in _bd._geometric_bonds(syms, P):
                    if syms[_i] == "H" or syms[_j] == "H":
                        continue
                    if _bd._is_metal(syms[_i]) or _bd._is_metal(syms[_j]):
                        continue
                    _got_nb.setdefault(_i, []).append(syms[_j])
                    _got_nb.setdefault(_j, []).append(syms[_i])
                _broken = False
                for _i, _s in enumerate(syms):
                    if _s == "H" or _bd._is_metal(_s):
                        continue
                    if _i not in _req_nb:
                        continue          # kein Sollwert -> nicht beurteilbar
                    if set(_te_C(_got_nb.get(_i, []))) != set(_te_C(_req_nb[_i])):
                        _broken = True
                        break
                if _broken:
                    continue              # ergaenzter Frame mit gebrochener Topologie
            except Exception:
                pass                      # nicht beurteilbar -> alte Sicherungen gelten
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
