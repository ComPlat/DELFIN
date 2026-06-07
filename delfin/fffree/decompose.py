"""delfin.fffree.decompose — SMILES → coordination decomposition for the
metal-FF-free backend.

Scope (falls back to None otherwise): mononuclear complex with explicit
metal-donor bonds in the graph, CN 4/5/6.  Returns the metal, CN, default
geometry, and per-ligand fragments (mol + local donor indices + donor elements +
denticity) — monodentate and chelating (bi-/tridentate) ligands are both
detected; ``has_chelate`` flags whether any ligand is polydentate.  CN outside
4-6, >1 metal, >tridentate ligands, or unparseable input -> return None.
"""
from __future__ import annotations
import os
from typing import Optional, Dict, List
from rdkit import Chem
import delfin._bond_decollapse as bd

# d8 square-planar-preferring metals (else CN4 -> tetrahedral)
_D8 = {"Pt", "Pd", "Ni", "Au", "Rh", "Ir"}


# ===== F1 coverage forensik: fine-grained failure reason recorder =====
# Default OFF (byte-identical). When DELFIN_FFFREE_DECOMPOSE_REASON_LOG points
# at a writable path, every decompose() return-None path records a one-line
# TSV reason for the SMILES.  Format: ``<smi>\t<reason>\t<extra>\n``.
# Used by the F1 coverage forensik harness to bucket the ~80 % UFF-fallback
# population by root cause without changing dispatch behaviour.


def _f1_log(smi: str, reason: str, extra: str = "") -> None:
    _path = os.environ.get("DELFIN_FFFREE_DECOMPOSE_REASON_LOG", "").strip()
    if not _path:
        return
    try:
        with open(_path, "a") as _fh:
            _fh.write(f"{smi}\t{reason}\t{extra}\n")
    except Exception:
        pass


def _default_geometry(metal: str, cn: int) -> Optional[str]:
    # Phase C f-block CN8-12 dispatch (Task #64, 2026-06-04).
    # Env-gated default OFF: only when DELFIN_FFFREE_FBLOCK_CN8_12=1 (or
    # PURE_TRACK3=1) AND metal is Ln/An AND CN in 8-12 do we route to the
    # dedicated f-block polyhedron (SAP-8, TTP-9, BICAP-10, CAP-11, IH-12).
    # Otherwise the legacy dispatch below runs unchanged (byte-identical).
    if (os.environ.get("DELFIN_FFFREE_FBLOCK_CN8_12", "0") == "1"
            or os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"):
        try:
            from delfin.fffree import f_block_polyhedra as _FBP
            if _FBP.is_f_block(metal):
                g = _FBP.default_geometry_fblock(metal, cn)
                if g is not None:
                    return g
        except ImportError:
            pass
    # Main-group LP-aware dispatch (2026-06-07, hmaximilian).  When the
    # metal carries a stereo-active lone pair (Sn²⁺, Pb²⁺, Sb³⁺, Bi³⁺,
    # In⁺, Tl⁺, Ge²⁺, As³⁺, Po⁴⁺) AND the observed CN is in the LP-aware
    # table (CN 1-5), the LP-aware polyhedron replaces the legacy choice.
    # Default ON under DELFIN_FFFREE_MAIN_GROUP_LP=1 (auto-on under
    # MOGUL_PRIMARY=1); byte-identical to legacy when unset.
    try:
        from delfin.fffree import main_group_polyhedron as _MGP
        if _MGP.main_group_lp_enabled():
            _mg = _MGP.main_group_polyhedron_for(metal, cn, None)
            if _mg is not None:
                return _mg
    except ImportError:
        pass
    if cn == 2:
        # Phase G: linear coordination (Cu(I), Ag(I), Au(I), Hg(II)).
        # Simple D∞h: 2 donors 180° apart on metal axis.
        return "L-2 linear"
    if cn == 6:
        return "OC-6 octahedron"
    if cn == 5:
        return "TBP-5 trigonal bipyramid"
    if cn == 4:
        return "SP-4 square planar" if metal in _D8 else "T-4 tetrahedron"
    if cn == 3:
        # iter-32c (User 2026-05-28 ADUMOD: Pd CN3 was built as Td/linear).
        # d⁸ Pt/Pd/Ni/Au/Rh/Ir CN3 prefer T-shape (square-planar with one vacancy);
        # all other metals get trigonal-planar SP-3.  Both isomers can be additively
        # enumerated via DELFIN_FFFREE_DUAL_CN3=1 (mirror of dual-CN4 / dual-CN6 wiring).
        return "T-3 T-shape" if metal in _D8 else "SP-3 trigonal planar"
    if cn == 7:
        return "PB-7 pentagonal bipyramid"
    if cn == 8:
        return "SQAP-8 square antiprism"
    if cn == 9:
        return "TTP-9 tricapped trigonal prism"
    # Mission A2 (2026-06-05): non-f-block CN10 polyhedra (BICAP-10 default).
    # Env-gated DELFIN_FFFREE_CN10_POLYHEDRA=1 (or PURE_TRACK3=1).  When unset
    # this falls through to ``return None`` so HEAD behaviour is byte-identical.
    if cn == 10:
        if (os.environ.get("DELFIN_FFFREE_CN10_POLYHEDRA", "0") == "1"
                or os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"):
            return "BICAP-10 bicapped square antiprism"
    return None


def decompose(smiles: str) -> Optional[Dict]:
    # Reuse the converter's full organometallic mol-preparation (stk / dative-bond
    # conversion / charge+H perception) so cleaved ligands have correct chemistry
    # (NH3 stays NH3, Cl⁻ stays Cl⁻ — not HCl).  Lazy import avoids the import
    # cycle (smiles_converter -> converter_backend -> decompose).
    mol = None
    try:
        from delfin.smiles_converter import _prepare_mol_for_embedding
        mol = _prepare_mol_for_embedding(smiles)
    except Exception:
        mol = None
    if mol is None:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        _f1_log(smiles, "smiles_unparseable")
        return None
    metals = [a.GetIdx() for a in mol.GetAtoms() if bd._is_metal(a.GetSymbol())]
    if len(metals) == 0:
        _f1_log(smiles, "no_metal")
        return None
    # Phase G4 (2026-05-31): multi-metal extension.
    # Previously: mononuclear only. Now: under PURE_TRACK3 or DELFIN_FFFREE_MULTI_METAL,
    # we pick the PRIMARY metal (max neighbours) and treat the other metal + its
    # local ligands as a single "metal-containing ligand" fragment.
    # Use case: HUCPIH (Sn-Ru bond), Mn2(CO)10, Re2(Cl)10 = pick the more
    # densely-coordinated metal as primary.
    # Expected CCDC build-rate gain: +10%.
    _PT3_MULTI = (os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
                  or os.environ.get("DELFIN_FFFREE_MULTI_METAL", "0") == "1")
    if len(metals) > 1:
        if not _PT3_MULTI:
            _f1_log(smiles, "multi_metal_disabled", f"n_metals={len(metals)}")
            return None                               # mononuclear only (legacy)
        # Pick primary metal = highest CN (most neighbors)
        m = max(metals, key=lambda mi: mol.GetAtomWithIdx(mi).GetDegree())
    else:
        m = metals[0]
    matom = mol.GetAtomWithIdx(m)
    donor_idx = [n.GetIdx() for n in matom.GetNeighbors()]
    cn = len(donor_idx)
    # High-CN (7-9) support is env-gated: default coverage stays CN 4-6 (byte-
    # identical), CN 7-9 polyhedra (PB-7/SQAP-8/TTP-9) only when DELFIN_FFFREE_HIGHCN=1.
    # Low-CN (3) support is env-gated: DELFIN_FFFREE_CN3=1 enables SP-3/T-3 (iter-32c).
    _allowed = set()
    _allowed.update({4, 5, 6})
    # Phase G3 (2026-05-31): under PURE_TRACK3, auto-enable CN2 (linear),
    # CN3 (SP-3/T-3), and CN7-9 (PB-7/SQAP-8/TTP-9). The hapto-π detection
    # below can collapse CN=9 (η6-arene + sigma) to effective CN=4-5, so
    # high CN must reach the fragment-analysis stage to be reduced.
    _PT3_AUTO = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
    if _PT3_AUTO or os.environ.get("DELFIN_FFFREE_HIGHCN", "0") == "1":
        _allowed.update({7, 8, 9})
    if _PT3_AUTO or os.environ.get("DELFIN_FFFREE_CN3", "0") == "1":
        _allowed.add(3)
    if _PT3_AUTO or os.environ.get("DELFIN_FFFREE_CN2", "0") == "1":
        _allowed.add(2)
    # Phase C f-block CN10-12 (Task #64): only enabled when explicit env
    # flag is set AND the metal is f-block (checked AFTER metal is known).
    # Setting only CN10-12 here is safe — _default_geometry returns None
    # for non-f-block at CN10-12, which bails to legacy (None return).
    if _PT3_AUTO or os.environ.get("DELFIN_FFFREE_FBLOCK_CN8_12", "0") == "1":
        _allowed.update({10, 11, 12})
    # Mission A2 (2026-06-05): non-f-block CN10 via BICAP-10/CSAP-10/SAP-10.
    # Independent env flag so the CN10 path can be exercised without enabling
    # the full f-block CN8-12 dispatch.
    if _PT3_AUTO or os.environ.get("DELFIN_FFFREE_CN10_POLYHEDRA", "0") == "1":
        _allowed.add(10)
    if cn not in _allowed:
        _f1_log(smiles, "cn_out_of_range", f"cn={cn}")
        return None
    metal = matom.GetSymbol()
    geometry = _default_geometry(metal, cn)
    if geometry is None:
        _f1_log(smiles, "no_geometry_for_cn", f"metal={metal} cn={cn}")
        return None

    donor_elem = {d: mol.GetAtomWithIdx(d).GetSymbol() for d in donor_idx}
    # break metal-donor bonds (keep atom indices stable), then split into
    # sanitized fragment mols (sanitize fixes valences -> N regains its H etc.).
    em = Chem.RWMol(mol)
    for d in donor_idx:
        if em.GetBondBetweenAtoms(m, d) is not None:
            em.RemoveBond(m, d)
    mapping: List = []
    # Phase G14-PIVOT (User 2026-05-31 night, 70.5 % legacy fallback discovery):
    # sanitized GetMolFrags fails on charged-aromatic NHC, carbonyl carbene,
    # iminophosphorane and other systems whose fragment valence becomes
    # invalid after the M-X bond is cleaved (e.g. [N+] in imidazolium losing
    # its C-M coordination -> N now has invalid coordination). Try sanitised
    # first (preserves existing behaviour); on failure under PT3, retry with
    # sanitizeFrags=False so the fragment graph survives. Expected native
    # rate gain: ~21/200 = 10.5 pp.
    # G14-PIVOT GetMolFrags-sanitize-False fallback REVERTED -- it rescued
    # 21/200 SMILES into fffree-native, but downstream assembly produced
    # rougher internals than the UFF-relaxed legacy fallback those SMILES
    # would have used, regressing the iter_gate basket by 21 severe axes.
    # Returns to the original sanitize=True only path.
    try:
        frags = Chem.GetMolFrags(em, asMols=True, sanitizeFrags=True,
                                 fragsMolAtomMapping=mapping)
    except Exception as _e:
        # MISSION F1 (2026-06-05) coverage heal #1: retry with
        # sanitizeFrags=False under DELFIN_FFFREE_F1_COVERAGE.  Charged-
        # aromatic NHC/carbene/iminophosphorane fragments have invalid
        # valence after M-X cleavage; sanitize-True rejects entirely.
        # G14-PIVOT pattern was reverted for INTERNALS; F1 is COVERAGE
        # mission -- internals are downstream's responsibility (GRIP,
        # UFF post-relax, post_grip_corrector).  Env-gated default OFF.
        if os.environ.get("DELFIN_FFFREE_F1_COVERAGE", "0") == "1":
            mapping = []
            try:
                frags = Chem.GetMolFrags(em, asMols=True, sanitizeFrags=False,
                                         fragsMolAtomMapping=mapping)
            except Exception as _e2:
                _f1_log(smiles, "frag_sanitize_failed_twice", str(_e2)[:80])
                return None
        else:
            _f1_log(smiles, "frag_sanitize_failed", str(_e)[:80])
            return None
    donor_set = set(donor_idx)
    ligands: List[Dict] = []
    n_chelate_bonds = 0
    for fmol, orig_idxs in zip(frags, mapping):
        orig = list(orig_idxs)
        if m in orig:
            continue                                  # the metal's own fragment
        fdonors = [o for o in orig if o in donor_set]
        if len(fdonors) == 0:
            # Phase G: bridging / spectator → previously hard fail. With
            # DELFIN_FFFREE_DECOMPOSE_EXTENDED=1, skip these fragments instead
            # (treat as auxiliary). Expected CCDC build-rate gain: +20%.
            if (os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
                or os.environ.get("DELFIN_FFFREE_DECOMPOSE_EXTENDED", "0") == "1"):
                continue
            _f1_log(smiles, "bridging_spectator_fragment")
            return None                               # bridging / spectator -> legacy
        # Phase G3 (User 2026-05-31): Hapto-π detection.
        # If a fragment has ≥3 carbon donors that are all part of the same
        # aromatic ring or conjugated π-system, classify as hapto (η3/η4/η5/η6/η7/η8).
        # Map to a single hapto-donor at the ring centroid. This admits
        # ferrocene/Cp/arene/COT complexes that previously hit the >tridentate gate.
        # Expected CCDC build-rate gain: +15%.
        _PT3_AUTO = (os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
                     or os.environ.get("DELFIN_FFFREE_HAPTO_DETECT", "0") == "1")
        is_hapto = False
        hapto_eta = 0
        if _PT3_AUTO and len(fdonors) >= 3:
            # Check if fdonors are all C (Cp/arene/allyl/diene) OR all N (porphyrin core)
            all_carbon = all(donor_elem[d] == "C" for d in fdonors)
            # Phase G3: ring-based hapto detection (works on SMILES that don't mark
            # aromaticity, e.g. [C+] cation Cp notation).
            try:
                local_donor_idxs = [orig.index(o) for o in fdonors]
                ring_info = fmol.GetRingInfo()
                # Find smallest ring containing ALL local-donor indices.
                for r in ring_info.AtomRings():
                    if all(li in r for li in local_donor_idxs):
                        # All donors lie in this ring → hapto-π
                        if all_carbon and len(fdonors) in (3, 4, 5, 6, 7, 8):
                            is_hapto = True
                            hapto_eta = len(fdonors)
                            break
                if not is_hapto:
                    # Try open-chain hapto (allyl η3, diene η4) via bond adjacency:
                    # donors form an unbroken chain in the fragment graph
                    if all_carbon and len(fdonors) in (3, 4):
                        # Check chain: each donor (except endpoints) bonded to 2 others
                        adj = {li: set() for li in local_donor_idxs}
                        for b in fmol.GetBonds():
                            a, bb = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
                            if a in adj and bb in adj:
                                adj[a].add(bb); adj[bb].add(a)
                        # Count degrees: chain = 2 endpoints (deg 1) + rest deg 2
                        degs = sorted(len(adj[li]) for li in local_donor_idxs)
                        if degs == [1, 1] + [2] * (len(fdonors) - 2):
                            is_hapto = True
                            hapto_eta = len(fdonors)
            except Exception:
                pass

        if len(fdonors) > 6 and not is_hapto:
            # Phase G: relax from > 3 to > 6 — allows tetra/penta/hexadentate
            # ligands (porphyrin κ4, EDTA κ6, salen κ4).
            # If hapto (η7 cycloheptatrienyl, η8 COT): allow.
            _f1_log(smiles, "hyperdentate_nonhapto", f"denticity={len(fdonors)}")
            return None                               # >hexadentate (very rare) -> legacy
        local_donors = [orig.index(o) for o in fdonors]
        if is_hapto:
            ligands.append({
                "mol": fmol,
                "donor_local_idx": local_donors[0],       # primary (back-compat; centroid effective)
                "donor_local_idxs": local_donors,         # all hapto-C atoms
                "denticity": 1,                            # Phase G3: hapto-π = 1 coord site
                "donor_elem": "C_hapto",                   # marker for downstream
                "donor_elems": ["C_hapto"],
                "is_hapto": True,
                "hapto_eta": hapto_eta,
            })
            n_chelate_bonds += 1                            # hapto = 1 coord site
        else:
            ligands.append({
                "mol": fmol,
                "donor_local_idx": local_donors[0],       # primary (back-compat)
                "donor_local_idxs": local_donors,         # all donors (chelate)
                "denticity": len(fdonors),
                "donor_elem": donor_elem[fdonors[0]],
                "donor_elems": [donor_elem[o] for o in fdonors],
                "is_hapto": False,
                "hapto_eta": 0,
            })
            n_chelate_bonds += len(fdonors)
    if n_chelate_bonds != cn:
        # Phase G3: hapto detection may reduce effective CN. Compute new effective CN.
        _PT3_AUTO = (os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
                     or os.environ.get("DELFIN_FFFREE_HAPTO_DETECT", "0") == "1")
        if _PT3_AUTO and any(lg.get("is_hapto") for lg in ligands):
            # CN is now sum of effective denticities (hapto=1, sigma=fdonors)
            new_cn = sum(lg["denticity"] for lg in ligands)
            # Phase C: extend upper cap to 12 when f-block env-gate is set.
            _CN_MAX = 12 if (
                os.environ.get("DELFIN_FFFREE_FBLOCK_CN8_12", "0") == "1"
                or os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
            ) else 9
            if 2 <= new_cn <= _CN_MAX:
                cn = new_cn
                # Re-compute geometry for new CN
                geometry = _default_geometry(metal, cn)
                if geometry is None:
                    _f1_log(smiles, "no_geometry_after_hapto_collapse",
                            f"metal={metal} cn={cn}")
                    return None
            else:
                _f1_log(smiles, "hapto_collapse_cn_out_of_range",
                        f"new_cn={new_cn}")
                return None
        else:
            _f1_log(smiles, "denticity_mismatch_cn",
                    f"cn={cn} sum_denticity={n_chelate_bonds}")
            return None
    has_chelate = any(lg["denticity"] >= 2 for lg in ligands)
    # Mission B1 (2026-06-05): sandwich / piano-stool / half-sandwich dispatch.
    # The legacy hapto path collapses the η-ring to a single point-donor on the
    # generic polyhedron (TBP-5, OC-6, TTP-9 …) which is geometrically wrong:
    # ferrocene becomes a CN2 dumbbell with M-Cp_centroid ≈ 1.0 Å instead of
    # 1.65 Å, and (η⁶-arene)Ru(L)₃ becomes a SQAP-8 with M-arene ≈ 1.31 Å
    # instead of ≈ 1.62 Å (Mission A7 forensik on pool 2792332: 414 half-sand
    # structures, median M-arene 1.313 Å vs ideal 1.62-1.69 Å).
    #
    # Under DELFIN_FFFREE_SANDWICH_DISPATCH=1 (auto-on under PURE_TRACK3=1) we
    # match three hapto-templates and override the generic geometry with one
    # of the dedicated sandwich polyhedra:
    #
    #   2× η⁵-Cp, nothing else                      → SANDWICH-10 (CN=10)
    #   1× η⁵-Cp + 3 σ-donors                       → PIANO-STOOL-8 (CN=8)
    #   1× η⁶-arene + 3 σ-donors                    → HALF-SANDWICH-9 (CN=9)
    #
    # Vertex set + Pólya groups for these names are already registered in
    # ``delfin.fffree.sandwich_piano_polyhedra`` and ``polya_isomer_count``;
    # ``polyhedra.ref_vectors`` / ``geometries_for_cn`` know how to dispatch
    # to them under the same env-gate.  Default OFF byte-identical with HEAD.
    _SANDWICH_DISPATCH = (
        os.environ.get("DELFIN_FFFREE_SANDWICH_DISPATCH", "0") == "1"
        or os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
    )
    if _SANDWICH_DISPATCH:
        hapto_etas = sorted(lg.get("hapto_eta", 0) for lg in ligands
                            if lg.get("is_hapto"))
        n_sigma = sum(1 for lg in ligands if not lg.get("is_hapto"))
        # The assembler counts ONE coordination site per hapto ring (see
        # ``hapto_modes.m_centroid_distance`` and the piano-stool placement at
        # ``assemble_complex.assemble_from_config`` lines 838-865), so the
        # effective ``cn`` exposed to the assembler is the EFFECTIVE site
        # count, NOT the atom-vertex count of the full polyhedron.  The
        # corresponding ``ref_vectors`` lookup in ``polyhedra`` returns the
        # per-site effective vertex array (1 ring axis + N σ vertices); see
        # ``sandwich_piano_polyhedra.effective_ref_vectors_sandwich``.
        #
        # 2× η⁵-Cp + nothing else → SANDWICH-10 (effective cn = 2)
        if hapto_etas == [5, 5] and n_sigma == 0:
            cn = 2
            geometry = "SANDWICH-10 bis-eta5-Cp"
        # 1× η⁵-Cp + 3 σ-donors → PIANO-STOOL-8 (effective cn = 4)
        elif hapto_etas == [5] and n_sigma == 3:
            cn = 4
            geometry = "PIANO-STOOL-8 eta5-Cp+L3"
        # 1× η⁶-arene + 3 σ-donors → HALF-SANDWICH-9 (effective cn = 4)
        elif hapto_etas == [6] and n_sigma == 3:
            cn = 4
            geometry = "HALF-SANDWICH-9 eta6+L3"
    # Ligand-complexity gate, measured PER DONOR ARM (heavy atoms / denticity):
    # small/simple ligands place cleanly; very large conjugated ligands need
    # conformer work the rigid placement doesn't do, so bail to None for those.
    # Per-arm (not per-ligand) so a polydentate chelate (e.g. citrate, 13 heavy
    # over 6 donors ~ 2/arm) is not unfairly rejected for its total size.
    #
    # Phase G (User 2026-05-31): with macrocycle.py + ring_pucker integration,
    # large polydentate ligands (porphyrin, salen) become handlable. When
    # DELFIN_FFFREE_DECOMPOSE_EXTENDED=1 (auto-on under PURE_TRACK3), raise the
    # per-donor heavy-atom cap from 8 → 20 to admit macrocyclic ligands.
    # Expected CCDC build-rate gain: +10% (porphyrin/salen/calix-class).
    _EXTENDED = (
        os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
        or os.environ.get("DELFIN_FFFREE_DECOMPOSE_EXTENDED", "0") == "1"
    )
    # G14-PIVOT REVERTED (smoke 520 vs eqn 517 = net -27, 21 severe regressions
    # including frame_pct_element_list_change +102 %, F3_bond +56 %,
    # lig_pct_realistic -45 %). Lifting decompose coverage admits more SMILES
    # into fffree-native quality territory, replacing the UFF-quality legacy
    # output for those SMILES in the comparison archive. The iter_gate then
    # measures fffree-native (rougher internals) vs UFF-legacy (relaxed)
    # on those rescued SMILES and finds fffree worse. This confirms the
    # diagnosis that the 63-axis gap is structural to fffree's construction
    # rougher-than-UFF internals; pure coverage shift cannot close it.
    # Keeping previous cap 20.
    # MISSION F1 (2026-06-05) coverage heal #2: raise cap further under
    # DELFIN_FFFREE_F1_COVERAGE. PT3 cap=20 already admits porphyrin/salen;
    # F1 cap=30 admits very large bidentate aryl/phosphine ligands (BINAP,
    # biphephos, fluorinated arylphosphines) common in catalysis.  Forensik
    # on smoke500/PT3=1 found 53 SMILES rejected at the cap=20 boundary --
    # cap=30 rescues most. Downstream relaxers handle the internals.
    if os.environ.get("DELFIN_FFFREE_F1_COVERAGE", "0") == "1":
        MAX_HEAVY_PER_DONOR = 30
    elif _EXTENDED:
        MAX_HEAVY_PER_DONOR = 20
    else:
        MAX_HEAVY_PER_DONOR = 8
    for lg in ligands:
        nheavy = sum(1 for a in lg["mol"].GetAtoms() if a.GetAtomicNum() > 1)
        if nheavy / max(lg["denticity"], 1) > MAX_HEAVY_PER_DONOR:
            _f1_log(smiles, "ligand_too_large",
                    f"nheavy={nheavy} denticity={lg['denticity']} cap={MAX_HEAVY_PER_DONOR}")
            return None
    return {"metal": metal, "cn": cn, "geometry": geometry,
            "has_chelate": has_chelate,
            "donor_elems": [lg["donor_elem"] for lg in ligands],
            "ligands": ligands}


if __name__ == "__main__":
    for label, smi in [("cisplatin", "N[Pt](N)(Cl)Cl"),
                       ("hexammineCo", "[NH3][Co]([NH3])([NH3])([NH3])([NH3])[NH3]"),
                       ("[CoCl3(NH3)3]", "[NH3][Co]([NH3])([NH3])([Cl])([Cl])[Cl]")]:
        d = decompose(smi)
        if d:
            print(f"{label:<16} metal={d['metal']} CN={d['cn']} geom={d['geometry']} "
                  f"donors={d['donor_elems']}")
        else:
            print(f"{label:<16} -> None (legacy fallback)")
