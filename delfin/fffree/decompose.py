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
        return None
    metals = [a.GetIdx() for a in mol.GetAtoms() if bd._is_metal(a.GetSymbol())]
    if len(metals) == 0:
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
            if 2 <= new_cn <= 9:
                cn = new_cn
                # Re-compute geometry for new CN
                geometry = _default_geometry(metal, cn)
                if geometry is None:
                    return None
            else:
                return None
        else:
            return None
    has_chelate = any(lg["denticity"] >= 2 for lg in ligands)
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
    MAX_HEAVY_PER_DONOR = 20 if _EXTENDED else 8
    for lg in ligands:
        nheavy = sum(1 for a in lg["mol"].GetAtoms() if a.GetAtomicNum() > 1)
        if nheavy / max(lg["denticity"], 1) > MAX_HEAVY_PER_DONOR:
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
