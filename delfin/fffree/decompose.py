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


def _default_geometry(metal: str, cn: int) -> Optional[str]:
    if cn == 2:
        # iter-32f (DELFIN_FFFREE_CN_EXTEND): linear two-coordinate — the canonical
        # d10 geometry (Cu(I)/Ag(I)/Au(I)/Hg(II)).  No metal-specific branch: CN2 is
        # essentially always linear.
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
    return None


def _rigid_hapto_enabled() -> bool:
    """Rigid-η construction on the FF-free path (default OFF -> byte-identical)."""
    return os.environ.get("DELFIN_FFFREE_RIGID_HAPTO", "0") == "1"


# Heavy-atom per-donor-arm complexity cap (decompose ligand-complexity gate).  At the
# default cap (8) a ligand with >8 heavy atoms / donor is REJECTED -> the whole complex
# bails to legacy (0 enumerated isomers).  The isomer-completeness audit
# (ISOMER_COMPLETENESS_AUDIT_2026_06_18 §4/§6) measured this single gate as the DOMINANT
# reach gap (~51.5 % of the pool: large/conjugated ligands — phosphines with big
# substituents, polypyridyl, large arene systems).  Conformer-aware SEATING
# (DELFIN_FFFREE_CONFORMER_SEATING, default OFF) raises this cap so those complexes
# reach the FF-free build, where the conformer/backbone-re-embed machinery seats the
# large ligand and the self-gate keeps the never-worse guarantee.
_HEAVY_CAP_DEFAULT = 8


def _seating_enabled() -> bool:
    """Conformer-aware seating for large ligands (default OFF -> byte-identical):
    when ON, the per-arm heavy-atom cap is raised so large-ligand complexes enter the
    FF-free path instead of bailing to legacy at decompose()."""
    return os.environ.get("DELFIN_FFFREE_CONFORMER_SEATING", "0") == "1"


def _heavy_cap() -> int:
    """Effective per-donor-arm heavy-atom cap.  Default 8 (byte-identical).  Raised to
    DELFIN_FFFREE_SEATING_HEAVY_CAP (default 24) when conformer seating is enabled, so
    realistic TMC ligands (polypyridyl, substituted phosphines, large arenes) reach the
    FF-free build; the conformer seating + self-gate enforce never-worse downstream."""
    if not _seating_enabled():
        return _HEAVY_CAP_DEFAULT
    try:
        return int(os.environ.get("DELFIN_FFFREE_SEATING_HEAVY_CAP", "24"))
    except Exception:
        return 24


def _eta_groups(mol, m: int) -> List[List[int]]:
    """Detect η (hapto) groups among the metal's carbon neighbours: maximal sets
    of ≥2 metal-bound carbons that are mutually contiguous through C-C bonds
    (a Cp / arene / diene / allyl π-face).  Graph-only, deterministic; returns a
    list of sorted carbon-index lists.  A lone metal-bound carbon (σ M-C, e.g. a
    carbonyl C or a methyl) is NOT an η group and is left out.

    Mirrors smiles_converter._find_hapto_groups but scoped to ONE metal and used
    only by the rigid-hapto FF-free path (env-gated)."""
    matom = mol.GetAtomWithIdx(m)
    c_nbrs = [n.GetIdx() for n in matom.GetNeighbors() if n.GetAtomicNum() == 6]
    cset = set(c_nbrs)
    if len(cset) < 2:
        return []
    seen: set = set()
    groups: List[List[int]] = []
    for start in c_nbrs:
        if start in seen:
            continue
        comp: List[int] = []
        stack = [start]
        seen.add(start)
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
                ni = nb.GetIdx()
                if ni in cset and ni not in seen:
                    seen.add(ni)
                    stack.append(ni)
        if len(comp) >= 2:                            # contiguous π-face = η group
            groups.append(sorted(comp))
    return groups


def _restore_eta_ring_hydrogens(fmol, local_eta: List[int]):
    """Restore the hydrogens dropped from η-ring carbons after the M-C cleave.

    The dative-bond SMILES encoding writes π-face carbons as ``[C+]`` cations
    (degree-3 with the metal bond).  When ``_decompose_hapto`` cleaves the M-C
    bonds the bare ring carbons keep ``noImplicit=True`` + a radical electron, so
    RDKit's ``AddHs`` adds ZERO hydrogen to them — the η-ring CH hydrogens vanish
    from every emitted frame (eye-caught: CEBQEI 11 H vs 17, BOCMIR 6 vs 10).

    Real Cp / arene / diene / allyl ring carbons are sp2 (valence 3): a ring
    carbon with no substituent (2 heavy neighbours) bears exactly ONE H; a
    substituted ring carbon (3 heavy neighbours) bears none.  We therefore drop
    the encoding-artifact formal charge + radical on each η-carbon and set its
    explicit-H count to ``max(0, 3 - heavy_neighbours)``.  The H land exo / in the
    ring plane via RDKit's geometric embed (handled downstream in
    ``_ligand_confs_from_mol``).  Graph-only, deterministic, idempotent; returns a
    sanitized RWMol-derived mol or the original fmol on any failure (never raises,
    never produces a non-finite molecule)."""
    try:
        rw = Chem.RWMol(fmol)
        for i in local_eta:
            a = rw.GetAtomWithIdx(int(i))
            if a.GetAtomicNum() != 6:
                continue                                  # η-faces are carbon
            n_heavy = sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() > 1)
            nh = max(0, 3 - n_heavy)                       # sp2 ring carbon, valence 3
            a.SetFormalCharge(0)                           # drop the [C+] artifact
            a.SetNumRadicalElectrons(0)
            a.SetNoImplicit(True)
            a.SetNumExplicitHs(nh)
        m2 = rw.GetMol()
        Chem.SanitizeMol(m2)
        return m2
    except Exception:
        return fmol


def _nhc_carbene_enabled() -> bool:
    """Free-carbene repair of metal-bound carbene-C after the M-C cleave (default
    OFF -> byte-identical: the repair branch is never entered)."""
    return os.environ.get("DELFIN_FFFREE_NHC_CARBENE", "0") == "1"


def _carbene_carbon_idxs(mol, m: int, matom) -> List[int]:
    """Graph-only detection of carbene-carbon DONORS on the metal.

    A carbene C (NHC = imidazol-2-yliden and its saturated/benzannulated kin, but
    also any divalent C-donor) is encoded in the dative-bond SMILES as ``[C+]``: a
    ring carbon that bonds the metal AND sits between heteroatoms (the canonical
    N-C-N of an NHC; generalised to ≥2 ring heteroatom neighbours, N/O/S/P).  Once
    the M-C bond is cleaved this carbon drops to degree 2 (two σ bonds to the ring
    heteroatoms) while still carrying the ``[C+]``/aromatic encoding artifact, which
    leaves RDKit unable to kekulise the ring (KekulizeException) -> the whole
    fragment split fails -> the system falls back to legacy.

    Returns the sorted local indices of such carbons (metal-bound, in a ring, with
    ≥2 ring heteroatom neighbours).  A lone σ-C (methyl, carbonyl C, aryl-C) has
    fewer than two heteroatom neighbours and is never flagged.  Deterministic."""
    out: List[int] = []
    for n in matom.GetNeighbors():
        if n.GetAtomicNum() != 6:
            continue
        if not n.IsInRing():
            continue
        het = 0
        for x in n.GetNeighbors():
            if x.GetIdx() == m:
                continue
            if x.GetAtomicNum() in (7, 8, 15, 16) and x.IsInRing():
                het += 1
        if het >= 2:
            out.append(n.GetIdx())
    return sorted(out)


def _repair_carbene_carbons(em, carbene_idxs: List[int]):
    """Repair cleaved carbene carbons IN-PLACE on the bond-removed RWMol ``em`` to
    the RDKit-valid representation of a FREE singlet carbene.

    After the M-C cleave the carbene carbon is a degree-2 carbon still tagged
    ``[C+]`` (charge +1, valence 3, the dative-bond encoding artifact) — RDKit
    cannot kekulise the surrounding ring, so ``SanitizeMol``/``sanitizeFrags``
    raises.  The chemically correct free species (after heterolytic M-C cleavage the
    σ lone pair stays on carbon) is the neutral divalent singlet carbene ``:C:`` —
    an sp² carbon with a σ lone pair and an empty p orbital, i.e. RDKit valence 2
    with 2 (non-bonding) radical electrons, charge 0, no implicit/explicit H.  This
    is exactly the representation RDKit accepts and round-trips, and embeds to the
    real ~101-106° N-C-N carbene angle.  Idempotent; only touches the flagged
    carbons.  Mutates ``em`` in place (the caller owns the RWMol)."""
    cset = set(int(i) for i in carbene_idxs)
    for i in cset:
        a = em.GetAtomWithIdx(i)
        if a.GetAtomicNum() != 6:
            continue
        a.SetFormalCharge(0)            # drop the [C+] dative-encoding artifact
        a.SetNumExplicitHs(0)
        a.SetNumRadicalElectrons(2)     # divalent singlet carbene :C: (val 2)
        a.SetNoImplicit(True)
        a.SetIsAromatic(False)          # carbene C is no longer ring-aromatic


def _decompose_hapto(smiles: str, mol, m: int, matom) -> Optional[Dict]:
    """Hapto-aware decomposition (DELFIN_FFFREE_RIGID_HAPTO=1 only).

    Treats each contiguous metal-bound π-face (Cp / arene / diene / allyl) as ONE
    coordination site occupying ONE polyhedron vertex (the ring centroid), so the
    effective CN is (#σ-donors) + (#η-groups) rather than the inflated η-carbon
    count.  Returns a decompose dict with per-ligand 'eta' metadata that
    assemble_complex builds as a rigid ring on the centroid vertex, or None to
    fall back to the legacy hapto path.  Universal, graph-only, deterministic."""
    egroups = _eta_groups(mol, m)
    if not egroups:
        return None                                   # no η-face -> not our case
    eta_atoms = set()
    for g in egroups:
        eta_atoms.update(g)
    nbr_idx = [n.GetIdx() for n in matom.GetNeighbors()]
    sigma_idx = [d for d in nbr_idx if d not in eta_atoms]
    cn = len(sigma_idx) + len(egroups)                # effective coordination number
    _allowed = {4, 5, 6}
    if os.environ.get("DELFIN_FFFREE_HIGHCN", "0") == "1":
        _allowed.update({7, 8, 9})
    if os.environ.get("DELFIN_FFFREE_CN3", "0") == "1":
        _allowed.add(3)
    if cn not in _allowed:
        return None
    metal = matom.GetSymbol()
    geometry = _default_geometry(metal, cn)
    if geometry is None:
        return None
    donor_elem = {d: mol.GetAtomWithIdx(d).GetSymbol() for d in sigma_idx}
    # Cleave: break every M-σ-donor bond AND every M-(η-carbon) bond, then split
    # into sanitized fragment mols (indices stay stable inside each fragment).
    em = Chem.RWMol(mol)
    for d in nbr_idx:
        if em.GetBondBetweenAtoms(m, d) is not None:
            em.RemoveBond(m, d)
    mapping: List = []
    try:
        frags = Chem.GetMolFrags(em, asMols=True, sanitizeFrags=True,
                                 fragsMolAtomMapping=mapping)
    except Exception:
        return None
    sigma_set = set(sigma_idx)
    # map global η-atom -> its group id (for per-fragment grouping)
    eta_gid = {}
    for gid, g in enumerate(egroups):
        for a in g:
            eta_gid[a] = gid
    ligands: List[Dict] = []
    n_sites = 0
    for fmol, orig_idxs in zip(frags, mapping):
        orig = list(orig_idxs)
        if m in orig:
            continue                                  # the metal's own fragment
        f_sigma = [o for o in orig if o in sigma_set]
        f_eta_gids = sorted({eta_gid[o] for o in orig if o in eta_gid})
        if not f_sigma and not f_eta_gids:
            return None                               # bridging / spectator -> legacy
        # A fragment may carry both σ-donors and η-faces (rare); but each must map
        # cleanly onto its own vertex.  Keep v1 simple+safe: a fragment is EITHER a
        # set of σ-donors (handled like the Werner path) OR exactly ONE η-face.
        if f_eta_gids:
            if len(f_eta_gids) != 1 or f_sigma:
                return None                           # mixed/multi η in one frag -> legacy
            gid = f_eta_gids[0]
            g = egroups[gid]
            local_eta = [orig.index(o) for o in g]
            # Restore the η-ring CH hydrogens the M-C cleave stripped (the [C+]
            # artifact otherwise has AddHs add zero H to the π-face carbons).
            fmol = _restore_eta_ring_hydrogens(fmol, local_eta)
            ligands.append({
                "mol": fmol,
                "is_eta": True,
                "eta_local_idxs": local_eta,          # ring carbons (local indices)
                "eta_n": len(g),                      # η hapticity (5=Cp, 6=arene, ...)
                "denticity": 1,                       # occupies ONE vertex (centroid)
                "donor_elem": "C",
            })
            n_sites += 1
        else:
            if len(f_sigma) > 3:
                return None                           # >tridentate (rare) -> legacy
            local_donors = [orig.index(o) for o in f_sigma]
            ligands.append({
                "mol": fmol,
                "is_eta": False,
                "donor_local_idx": local_donors[0],
                "donor_local_idxs": local_donors,
                "denticity": len(f_sigma),
                "donor_elem": donor_elem[f_sigma[0]],
                "donor_elems": [donor_elem[o] for o in f_sigma],
            })
            n_sites += len(f_sigma)
    if n_sites != cn:
        return None
    # Ligand-complexity gate (mirror of the Werner path), but EXEMPT η-faces: an
    # η-ring is placed rigidly as one unit, so its heavy-atom count is irrelevant.
    MAX_HEAVY_PER_DONOR = _heavy_cap()    # raised under conformer-seating (default 8)
    for lg in ligands:
        if lg.get("is_eta"):
            continue
        nheavy = sum(1 for a in lg["mol"].GetAtoms() if a.GetAtomicNum() > 1)
        if nheavy / max(lg["denticity"], 1) > MAX_HEAVY_PER_DONOR:
            return None
    has_chelate = any((not lg.get("is_eta")) and lg["denticity"] >= 2
                      for lg in ligands)
    has_eta = any(lg.get("is_eta") for lg in ligands)
    return {"metal": metal, "cn": cn, "geometry": geometry,
            "has_chelate": has_chelate, "has_eta": has_eta,
            "donor_elems": [lg["donor_elem"] for lg in ligands],
            "ligands": ligands}


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
        return None
    metals = [a.GetIdx() for a in mol.GetAtoms() if bd._is_metal(a.GetSymbol())]
    if len(metals) != 1:
        return None                                   # mononuclear only (v1)
    m = metals[0]
    matom = mol.GetAtomWithIdx(m)
    # Rigid-hapto path (env-gated, default OFF): if the metal carries a contiguous
    # π-face (Cp / arene / diene / allyl), collapse each face to ONE coordination
    # site so the effective CN passes the 4-6 gate and reaches the FF-free build.
    # Byte-identical when the flag is off (the branch is never entered).
    if _rigid_hapto_enabled():
        d_hap = _decompose_hapto(smiles, mol, m, matom)
        if d_hap is not None:
            return d_hap
        # no η-face (or unhandled hapto topology) -> fall through to the Werner path,
        # which is byte-identical to the flag-off behaviour for non-hapto inputs.
    donor_idx = [n.GetIdx() for n in matom.GetNeighbors()]
    cn = len(donor_idx)
    # High-CN (7-9) support is env-gated: default coverage stays CN 4-6 (byte-
    # identical), CN 7-9 polyhedra (PB-7/SQAP-8/TTP-9) only when DELFIN_FFFREE_HIGHCN=1.
    # Low-CN (3) support is env-gated: DELFIN_FFFREE_CN3=1 enables SP-3/T-3 (iter-32c).
    # CN-coverage extension (iter-32f): DELFIN_FFFREE_CN_EXTEND=1 accepts the most
    # impactful out-of-range buckets — CN2 (linear, the big one: d10 Cu/Ag/Au/Hg) plus
    # CN7 (PB-7) and CN8 (SQAP-8) — routing them to the polyhedra that already exist.
    # Default OFF -> byte-identical (the branch is never entered).
    _allowed = set()
    _allowed.update({4, 5, 6})
    if os.environ.get("DELFIN_FFFREE_HIGHCN", "0") == "1":
        _allowed.update({7, 8, 9})
    if os.environ.get("DELFIN_FFFREE_CN3", "0") == "1":
        _allowed.add(3)
    if os.environ.get("DELFIN_FFFREE_CN_EXTEND", "0") == "1":
        _allowed.update({2, 7, 8})
    if cn not in _allowed:
        return None
    metal = matom.GetSymbol()
    geometry = _default_geometry(metal, cn)
    if geometry is None:
        return None

    donor_elem = {d: mol.GetAtomWithIdx(d).GetSymbol() for d in donor_idx}
    # break metal-donor bonds (keep atom indices stable), then split into
    # sanitized fragment mols (sanitize fixes valences -> N regains its H etc.).
    em = Chem.RWMol(mol)
    for d in donor_idx:
        if em.GetBondBetweenAtoms(m, d) is not None:
            em.RemoveBond(m, d)
    mapping: List = []
    try:
        frags = Chem.GetMolFrags(em, asMols=True, sanitizeFrags=True,
                                 fragsMolAtomMapping=mapping)
    except Exception:
        # NHC / carbene-C donor repair (DELFIN_FFFREE_NHC_CARBENE=1 only; default OFF
        # -> this branch is never entered and the failure stays a legacy fallback,
        # byte-identical).  Defect-driven: the fragment split only fails (kekulize/
        # valence) when a cleaved carbene-C kept its [C+] encoding; repair those
        # carbons to the free singlet carbene :C: and retry the split ONCE.  A
        # genuinely unsplittable input (no carbene, real valence error) still has no
        # carbene candidates -> no repair -> return None unchanged.
        frags = None
        if _nhc_carbene_enabled():
            carbene_idxs = _carbene_carbon_idxs(mol, m, matom)
            if carbene_idxs:
                _repair_carbene_carbons(em, carbene_idxs)
                mapping = []
                try:
                    frags = Chem.GetMolFrags(em, asMols=True, sanitizeFrags=True,
                                             fragsMolAtomMapping=mapping)
                except Exception:
                    frags = None
        if frags is None:
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
            return None                               # bridging / spectator -> legacy
        if len(fdonors) > 3:
            return None                               # >tridentate (rare) -> legacy
        local_donors = [orig.index(o) for o in fdonors]
        ligands.append({
            "mol": fmol,
            "donor_local_idx": local_donors[0],       # primary (back-compat)
            "donor_local_idxs": local_donors,         # all donors (chelate)
            "denticity": len(fdonors),
            "donor_elem": donor_elem[fdonors[0]],
            "donor_elems": [donor_elem[o] for o in fdonors],
        })
        n_chelate_bonds += len(fdonors)
    if n_chelate_bonds != cn:
        return None
    has_chelate = any(lg["denticity"] >= 2 for lg in ligands)
    # Ligand-complexity gate, measured PER DONOR ARM (heavy atoms / denticity):
    # small/simple ligands place cleanly; very large conjugated ligands need
    # conformer work the rigid placement doesn't do, so the DEFAULT bails to None
    # for those.  Per-arm (not per-ligand) so a polydentate chelate (e.g. citrate,
    # 13 heavy over 6 donors ~ 2/arm) is not unfairly rejected for its total size.
    # Conformer seating (DELFIN_FFFREE_CONFORMER_SEATING) raises the cap so the large
    # ligand reaches the build, where the conformer-aware seating (+ backbone re-embed)
    # places it and the self-gate enforces never-worse.  Default OFF -> cap 8 (byte-id).
    MAX_HEAVY_PER_DONOR = _heavy_cap()
    #
    # JOINT-DECLASH gate-lift (DELFIN_FFFREE_JOINT_DECLASH, default OFF -> byte-id):
    # large MONODENTATE ligands (per-arm 9-12 heavy: PMePh2, mesityloxide, ArO,
    # N(SiMe3)2, ...) are the cleanest "class-B" sub-population — their coordination
    # core builds IDEAL and only the ligand BODIES clash, which the new joint
    # inter-ligand declash pass now relieves before the self-gate.  So raise the
    # per-arm ceiling for MONODENTATE arms only (8 -> 12).  CHELATE arms keep the
    # historic cap (their backbone needs ring/metallacycle work the declash does not
    # do); denticity>3 / kappa4 already bailed above (class-C, separate).
    _jd = os.environ.get("DELFIN_FFFREE_JOINT_DECLASH", "0") == "1"
    _mono_cap = int(os.environ.get("DELFIN_FFFREE_MONO_HEAVY_CAP", "12")) if _jd else 8
    for lg in ligands:
        nheavy = sum(1 for a in lg["mol"].GetAtoms() if a.GetAtomicNum() > 1)
        # union of both flags: monodentate arms may be raised by EITHER the
        # seating heavy-cap OR the joint-declash mono-cap (whichever is larger);
        # chelate arms keep the seating-aware base cap.  Both flags OFF => 8 (byte-id).
        cap = max(MAX_HEAVY_PER_DONOR, _mono_cap) if lg["denticity"] == 1 else MAX_HEAVY_PER_DONOR
        if nheavy / max(lg["denticity"], 1) > cap:
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
