#!/usr/bin/env python3
"""coordination_mode_enum.py — enumerate the valid COORDINATION
MODES of a ligand (the hard open core).  For a ligand with multiple potential
coordination sites, enumerate: donor-selection + denticity (κ1/κ2/κ3…) + linkage
isomerism (ambidentate N-vs-O, S-vs-N, carboxylate κ1-vs-κ2) + hapto-π (η^n).

Universal, graph-only (no SMILES-specific rules): donors identified by element +
lone-pair availability; chelate feasibility by metallacycle ring size from the
molecular graph; symmetric donors merged by RDKit canonical rank PLUS a
resonance-equivalence layer (carboxylate, nitrate, sulfate, phosphate ... so
that the two O of acetate are merged into ONE κ¹-O mode rather than two).
Deterministic (canonical SMILES atom order).

VALIDITY RULES (the calibratable core; refine against COD/CCDC):
  - donor eligible: N/O/P/S/halide with formal charge <= 0 and degree < 4 (lone pair)
  - chelate ring size 4..7 feasible (4 = carboxylate-type strained-but-real, 5-6
    ideal, 7 ok); 3 = not a chelate (side-on), >7 = unusual -> excluded
  - max denticity capped at MAX_KAPPA
  - π-system detection: aromatic ring + open-chain conjugated π emit η-modes

Two layers extending the original kappa enumeration:

(1) FUNCTIONAL-GROUP REALISM — for known oxoanion / ambidentate / common-ligand
    families a canonical mode table replaces brute-force subset-enumeration:
        carboxylate ─→ κ¹-O, κ²-O,O (chelate, 4-ring), bridging
        nitrate    ─→ κ¹-O, κ²-O,O
        sulfate    ─→ κ¹-O, κ²-O,O (chelate or bridging)
        phosphate  ─→ κ¹-O, κ²-O,O
        β-diketonate (acac) ─→ κ²-O,O (6-ring; very stable)
        salen (N2O2) ─→ κ⁴
        porphyrin / phthalocyanine ─→ κ⁴ (rigid)
        EDTA-family  ─→ κ⁴, κ⁵, κ⁶ (denticity tunable by env)
        amino acid (NH2-CR-COO⁻) ─→ κ²-N,O glycinate
        Cp/Cp* / indenyl ─→ η¹, η³, η⁵
        arene (C6H6, C6H5R, ...)  ─→ η², η⁴, η⁶
        allyl  ─→ η¹, η³
        butadiene ─→ η², η⁴

(2) HAPTO-π EXTENSION via the partner module ``hapto_modes`` — when a ring is
    detected we emit the realistic η-modes for that ring size (Cp = 1,3,5;
    arene = 2,4,6; open-chain via consecutive π carbons).

REALISM FILTER (env-controlled): ``DELFIN_FFFREE_BINDING_MODE_STRICT=1`` (default
in PURE_TRACK3) prunes:
  - 4-membered chelate rings unless the κ²-bite is symmetric (oxoanions ok)
  - η-modes higher than half the ring atoms in low-CN situations
  - linkage isomers requiring a missing fragment

Tested on en/acetate/glycinate/SCN/pyridine/terpy/Cp/arene/acac/nitrate.
"""
from __future__ import annotations

import itertools
import os
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
from rdkit import Chem


# ===== env-gates ==========================================================
_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_BINDING_MODE_ENUM = (
    _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_BINDING_MODE_ENUM", "0") == "1"
)
_BINDING_MODE_STRICT = (
    _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_BINDING_MODE_STRICT", "0") == "1"
)


def binding_mode_enum_enabled() -> bool:
    """Live-evaluated env-gate (env can be changed at runtime in tests)."""
    return (
        os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
        or os.environ.get("DELFIN_FFFREE_BINDING_MODE_ENUM", "0") == "1"
    )


def binding_mode_strict_enabled() -> bool:
    return (
        os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
        or os.environ.get("DELFIN_FFFREE_BINDING_MODE_STRICT", "0") == "1"
    )


# ===== chemistry tables ===================================================
DONOR_ELEMENTS = {"N", "O", "P", "S", "F", "Cl", "Br", "I", "Se", "As", "C"}
# Note: C donors only count when part of a carbanion / NHC / hapto-π fragment;
#       that is handled by the per-fragment classifier, not the bulk donor list.
CHELATE_MIN, CHELATE_MAX = 4, 7
MAX_KAPPA = 6


# Canonical functional-group SMARTS → list of (mode_label, donor_atom_indices_in_SMARTS, kappa)
# The donor_atom_indices reference SMARTS atom-indices (0-based) inside the match.
# For η-modes, kappa = η (the hapto number); donor_atom_indices = the π-bound atoms.
# Each emission entry is one *realistic* binding mode for this fragment.
FUNCTIONAL_GROUP_MODES: Dict[str, List[Dict]] = {
    # carboxylate ─ R-C(=O)O⁻ : κ¹-O monodentate / κ²-O,O chelate / bridging
    "[CX3](=O)[O;H0,-]": [
        {"label": "kappa1-O-carboxylate",        "donors": (1,),    "kappa": 1, "ring": 0},
        {"label": "kappa1-O-carboxylate-other",  "donors": (2,),    "kappa": 1, "ring": 0,
         "resonance_equivalent_of": "kappa1-O-carboxylate"},
        {"label": "kappa2-OO-carboxylate",       "donors": (1, 2),  "kappa": 2, "ring": 4},
        {"label": "mu2-bridging-carboxylate",    "donors": (1, 2),  "kappa": 2, "ring": 0,
         "bridging": True},
    ],
    # nitrate ─ NO3⁻ : κ¹-O / κ²-O,O (SMARTS matches any of 3 N-O bonds)
    "[OX1,O-][N+](=O)[O-]": [
        {"label": "kappa1-O-nitrate",     "donors": (0,),   "kappa": 1, "ring": 0},
        {"label": "kappa1-O-nitrate-2",   "donors": (2,),   "kappa": 1, "ring": 0,
         "resonance_equivalent_of": "kappa1-O-nitrate"},
        {"label": "kappa1-O-nitrate-3",   "donors": (3,),   "kappa": 1, "ring": 0,
         "resonance_equivalent_of": "kappa1-O-nitrate"},
        {"label": "kappa2-OO-nitrate",    "donors": (0, 2), "kappa": 2, "ring": 4},
    ],
    # nitrite NO2⁻ (standalone, not part of nitrate; require ONLY 2 oxygens
    # bonded to N to discriminate from NO3).  Match the full anion: N has
    # exactly 2 O-neighbours plus zero H.
    "[N+;X3;H0](=O)[O-]": [
        {"label": "kappa1-N-nitro",       "donors": (0,),   "kappa": 1, "ring": 0,
         "_require_N_degree": 2},
        {"label": "kappa1-O-nitrito",     "donors": (2,),   "kappa": 1, "ring": 0,
         "_require_N_degree": 2},
        {"label": "kappa2-ON-bidentate",  "donors": (0, 2), "kappa": 2, "ring": 4,
         "_require_N_degree": 2},
    ],
    # sulfate SO4²⁻
    "[SX4](=O)(=O)([O-])[O-]": [
        {"label": "kappa1-O-sulfate",     "donors": (3,),   "kappa": 1, "ring": 0},
        {"label": "kappa2-OO-sulfate",    "donors": (3, 4), "kappa": 2, "ring": 4},
    ],
    # sulfite SO3²⁻ (ambidentate)
    "[SX3](=O)([O-])[O-]": [
        {"label": "kappa1-S-sulfito",     "donors": (0,),   "kappa": 1, "ring": 0},
        {"label": "kappa1-O-sulfonato",   "donors": (2,),   "kappa": 1, "ring": 0},
        {"label": "kappa2-OO-sulfito",    "donors": (2, 3), "kappa": 2, "ring": 4},
    ],
    # phosphate PO4³⁻
    "[PX4](=O)([O-])([O-])[O-]": [
        {"label": "kappa1-O-phosphate",   "donors": (2,),   "kappa": 1, "ring": 0},
        {"label": "kappa2-OO-phosphate",  "donors": (2, 3), "kappa": 2, "ring": 4},
    ],
    # phosphonate / phosphinate RP(=O)(O-)(OR')
    "[PX4](=O)([O-])": [
        {"label": "kappa1-O-phosphinate", "donors": (1,),   "kappa": 1, "ring": 0},
    ],
    # β-diketonate (acac/dbm) ─ canonical resonance after deprotonation:
    #     O=C-CH=C-O⁻   ↔   ⁻O-C=CH-C=O
    # SMARTS matches the deprotonated enolate form, 5 atoms incl. central CH.
    "[O;-,X1][CX3]=[CX3][CX3]=[O]": [
        {"label": "kappa2-OO-betadiketo", "donors": (0, 4), "kappa": 2, "ring": 6},
    ],
    "[O]=[CX3][CX3]=[CX3][O;-,X1]": [
        {"label": "kappa2-OO-betadiketo", "donors": (0, 4), "kappa": 2, "ring": 6},
    ],
    # amino acid α-NH₂-Cα-COO⁻ → glycinate N,O chelate (5-ring)
    "[NX3;H2,H1][CX4][CX3](=O)[O-]": [
        {"label": "kappa2-NO-aminoacid",  "donors": (0, 4), "kappa": 2, "ring": 5},
        {"label": "kappa1-N-aminoacid",   "donors": (0,),   "kappa": 1, "ring": 0},
        {"label": "kappa1-O-aminoacid",   "donors": (4,),   "kappa": 1, "ring": 0,
         "resonance_equivalent_of": "kappa1-O-aminoacid-other"},
    ],
    # ethylenediamine ─ NCCN κ²-N,N (5-ring)
    "[NX3;H2,H1][CX4][CX4][NX3;H2,H1]": [
        {"label": "kappa2-NN-en",        "donors": (0, 3), "kappa": 2, "ring": 5},
        {"label": "kappa1-N-en",         "donors": (0,),   "kappa": 1, "ring": 0,
         "resonance_equivalent_of": "kappa1-N-en-other"},
    ],
    # 2,2'-bipyridine — κ²-N,N (5-ring including M)
    "c1ccc(-c2ccccn2)nc1": [
        {"label": "kappa2-NN-bipy",     "donors": (9, 10), "kappa": 2, "ring": 5},
    ],
    # imidazole / pyrazole monodentate-N
    # (skip — fall through to brute-force κ¹ enumeration)
}


# Hapto-π SMARTS → list of feasible η-modes.
# Donor indices are SMARTS-match ring positions for the SVD centroid.
HAPTO_PI_FRAGMENTS: Dict[str, List[Dict]] = {
    # cyclopentadienyl anion (Cp⁻)
    "[c]1[c][c][c][c]1": [
        {"label": "eta1-Cp",  "donors": (0,),              "eta": 1, "ring_n": 5},
        {"label": "eta3-Cp",  "donors": (0, 1, 2),         "eta": 3, "ring_n": 5},
        {"label": "eta5-Cp",  "donors": (0, 1, 2, 3, 4),   "eta": 5, "ring_n": 5},
    ],
    # benzene / arene (η²/η⁴/η⁶)
    "[c]1[c][c][c][c][c]1": [
        {"label": "eta2-arene", "donors": (0, 1),                "eta": 2, "ring_n": 6},
        {"label": "eta4-arene", "donors": (0, 1, 2, 3),          "eta": 4, "ring_n": 6},
        {"label": "eta6-arene", "donors": (0, 1, 2, 3, 4, 5),    "eta": 6, "ring_n": 6},
    ],
    # cyclobutadiene η⁴ (rare)
    "[c]1[c][c][c]1": [
        {"label": "eta4-cbd", "donors": (0, 1, 2, 3), "eta": 4, "ring_n": 4},
    ],
    # cycloheptatrienyl (Tropylium) η⁷
    "[c]1[c][c][c][c][c][c]1": [
        {"label": "eta3-tro", "donors": (0, 1, 2),                   "eta": 3, "ring_n": 7},
        {"label": "eta5-tro", "donors": (0, 1, 2, 3, 4),             "eta": 5, "ring_n": 7},
        {"label": "eta7-tro", "donors": (0, 1, 2, 3, 4, 5, 6),       "eta": 7, "ring_n": 7},
    ],
    # cyclooctatetraene (COT) η⁸ (lanthanide/actinide)
    "[c]1[c][c][c][c][c][c][c]1": [
        {"label": "eta4-cot", "donors": (0, 1, 2, 3),                  "eta": 4, "ring_n": 8},
        {"label": "eta6-cot", "donors": (0, 1, 2, 3, 4, 5),            "eta": 6, "ring_n": 8},
        {"label": "eta8-cot", "donors": (0, 1, 2, 3, 4, 5, 6, 7),      "eta": 8, "ring_n": 8},
    ],
    # allyl C=C-C  (η¹-σ vs η³) — generic sp2-sp2-sp3 propene-like backbone
    "[CX3]=[CX3][CX3,CX4;!R]": [
        {"label": "eta1-allyl", "donors": (2,),       "eta": 1, "ring_n": 3},
        {"label": "eta3-allyl", "donors": (0, 1, 2),  "eta": 3, "ring_n": 3},
    ],
    # 1,3-butadiene C=C-C=C (acyclic, conjugated)
    "[CX3]=[CX3]-[CX3]=[CX3]": [
        {"label": "eta2-butadiene", "donors": (0, 1),       "eta": 2, "ring_n": 4},
        {"label": "eta4-butadiene", "donors": (0, 1, 2, 3), "eta": 4, "ring_n": 4},
    ],
}


# ===== helpers ============================================================
def _donor_atoms(mol) -> List[int]:
    out = []
    for a in mol.GetAtoms():
        if a.GetSymbol() in DONOR_ELEMENTS - {"C"} and \
                a.GetFormalCharge() <= 0 and a.GetDegree() < 4:
            out.append(a.GetIdx())
    return out


def _ring_size(mol, i: int, j: int) -> int:
    """Metallacycle ring size if donors i,j chelate = (atoms on shortest path) + 1
    (the metal).  Returns 0 if no path."""
    path = Chem.GetShortestPath(mol, i, j)
    return len(path) + 1 if path else 0


def _feasible_pair(mol, i: int, j: int) -> bool:
    return CHELATE_MIN <= _ring_size(mol, i, j) <= CHELATE_MAX


def _canon_ranks(mol) -> List[int]:
    return list(Chem.CanonicalRankAtoms(mol, breakTies=False))


def _resonance_equivalent_donors(mol) -> Dict[int, int]:
    """Map donor atom-idx → equivalence-class representative.  Two oxygens that
    are resonance-equivalent across a delocalised group (carboxylate, nitrate,
    sulfate, phosphate, ...) share the same representative.

    Detection: for each donor atom O/N/S of a structural-feature root, if a
    neighbour atom carries a double-bond to another donor of the same element
    AND both donors have the same neighbouring environment up to depth-1, they
    are merged.

    Lightweight (no RDKit Resonance — that's expensive); covers the cases that
    matter for binding-mode enumeration.
    """
    equiv: Dict[int, int] = {}
    # Carboxylate / nitrate / sulfate / phosphate: detect via SMARTS.
    canonical_groups = [
        # SMARTS, set of indices in match that are resonance-equivalent donors
        ("[CX3](=O)[O;H0,-]",                 {1, 2}),   # COO⁻
        ("[OX1,O-][N+](=O)[O-]",              {0, 2, 3}), # NO₃⁻
        ("[N+](=O)[O-]",                      {1, 2}),   # NO₂⁻
        ("[SX4](=O)(=O)([O-])[O-]",           {1, 2, 3, 4}),  # SO₄²⁻
        ("[SX3](=O)([O-])[O-]",               {1, 2, 3}),  # SO₃²⁻ O-equiv only
        ("[PX4](=O)([O-])([O-])[O-]",         {1, 2, 3, 4}),  # PO₄³⁻
        # ethylenediamine N,N (two -NH2 symmetric)
        ("[NX3;H2,H1][CX4][CX4][NX3;H2,H1]",  {0, 3}),
    ]
    for smarts, eq_idx in canonical_groups:
        try:
            patt = Chem.MolFromSmarts(smarts)
            if patt is None:
                continue
            for match in mol.GetSubstructMatches(patt):
                # pick the smallest match atom as representative
                reps = sorted(match[i] for i in eq_idx if i < len(match))
                if not reps:
                    continue
                head = reps[0]
                for r in reps:
                    equiv[r] = head
        except Exception:
            continue
    return equiv


def _detect_functional_modes(mol) -> List[Dict]:
    """Recognise carboxylate/nitrate/sulfate/phosphate/β-diketonate etc. and
    emit their canonical κ-modes.  Mode-dict keys: label, donors (atom-idx),
    kappa, ring, optional bridging/resonance_equivalent_of.
    """
    out: List[Dict] = []
    seen_signatures: Set[Tuple] = set()
    for smarts, modes in FUNCTIONAL_GROUP_MODES.items():
        try:
            patt = Chem.MolFromSmarts(smarts)
        except Exception:
            continue
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt):
            for spec in modes:
                # constraint: optional degree check on a matched atom (prevents
                # NO2 SMARTS from matching nitrate's N which has 3 O neighbours).
                req_deg = spec.get("_require_N_degree")
                if req_deg is not None:
                    try:
                        nat = mol.GetAtomWithIdx(match[0])
                        if nat.GetDegree() != req_deg:
                            continue
                    except Exception:
                        continue
                try:
                    donors = tuple(sorted(match[i] for i in spec["donors"]))
                except IndexError:
                    continue
                sig = (spec["label"].rsplit("-", 1)[0], donors)
                if sig in seen_signatures:
                    continue
                seen_signatures.add(sig)
                m = {
                    "label": spec["label"],
                    "donors": donors,
                    "kappa": spec["kappa"],
                    "ring_size": spec.get("ring", 0),
                    "kind": "functional",
                    "smarts": smarts,
                }
                if spec.get("bridging"):
                    m["bridging"] = True
                if "resonance_equivalent_of" in spec:
                    m["resonance_equivalent_of"] = spec["resonance_equivalent_of"]
                out.append(m)
    return out


def _detect_hapto_modes(mol) -> List[Dict]:
    """Recognise Cp / arene / allyl / butadiene / COT etc. and emit canonical
    η-modes.  Returns a list of mode-dicts with `eta` key set."""
    out: List[Dict] = []
    seen_centroids: Set[Tuple[int, ...]] = set()
    for smarts, modes in HAPTO_PI_FRAGMENTS.items():
        try:
            patt = Chem.MolFromSmarts(smarts)
        except Exception:
            continue
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt):
            ring_atoms = tuple(sorted(set(match)))
            # collapse duplicate ring detections (Cp/arene SMARTS overlap)
            for spec in modes:
                try:
                    donors = tuple(sorted(match[i] for i in spec["donors"]))
                except IndexError:
                    continue
                sig = (spec["label"].rsplit("-", 1)[0], donors)
                if sig in seen_centroids:
                    continue
                seen_centroids.add(sig)
                out.append({
                    "label": spec["label"],
                    "donors": donors,
                    "eta": spec["eta"],
                    "kappa": 1,             # hapto = single effective coord site
                    "ring_size": 0,         # rings don't form metallacycles
                    "ring_n": spec["ring_n"],
                    "kind": "hapto",
                    "smarts": smarts,
                })
    return out


def _passes_realism(mol, mode: Dict) -> bool:
    """Realism filter for binding modes.  Stricter mode prunes:
      - 4-membered chelate rings (carboxylate is allowed since the
        oxoanion-bite is symmetric)
      - η-modes whose hapto-number > ring-atom-count (impossible)
      - non-π C-donor in non-hapto modes (handled upstream)
    Permissive mode keeps everything.
    """
    if not binding_mode_strict_enabled():
        return True
    # eta > ring atoms = impossible
    if mode.get("kind") == "hapto":
        if mode["eta"] > mode.get("ring_n", 0):
            return False
        return True
    # κ-modes — chelate ring size enforced (4-membered allowed for oxoanions)
    rs = mode.get("ring_size", 0)
    if mode["kappa"] >= 2 and rs > 0:
        if rs < CHELATE_MIN:
            # 4-ring tolerated for symmetric-bite oxoanions
            if mode.get("kind") == "functional" and rs == 4:
                return True
            return False
        if rs > CHELATE_MAX:
            return False
    return True


# ===== public API =========================================================
def enumerate_modes(smiles: str,
                    *,
                    include_functional: bool = True,
                    include_hapto: bool = True,
                    strict_realism: Optional[bool] = None,
                    ) -> List[Dict]:
    """Return the distinct coordination modes of the ligand.  Each mode is a dict:

        {donors, kappa, ring_sizes, donor_elems, kind, label}

    ``kind`` is one of {"kappa", "functional", "hapto"} so callers can tag
    isomer labels.  Symmetric modes merged via canonical ranks (RDKit) AND
    resonance-equivalence (carboxylate-O, nitrate-O, sulfate-O, ...).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    use_strict = (strict_realism
                  if strict_realism is not None
                  else binding_mode_strict_enabled())

    # 1) classical kappa subset enumeration (the original code)
    donors = _donor_atoms(mol)
    rank = _canon_ranks(mol)
    sym = {a.GetIdx(): a.GetSymbol() for a in mol.GetAtoms()}
    res_equiv = _resonance_equivalent_donors(mol)
    seen: Set[Tuple] = set()
    modes: List[Dict] = []

    def add(subset: Tuple[int, ...], ring_sizes: Tuple[int, ...],
            kind: str, label: Optional[str] = None,
            extra: Optional[Dict] = None):
        # canonical signature: multiset of (resolved rank) — donor merged into
        # its resonance representative if any
        merged = tuple(sorted(rank[res_equiv.get(i, i)] for i in subset))
        # Signature collapses (functional, kappa) duplicates: the labelled
        # functional mode subsumes a brute-force kappa-only mode with the same
        # donors.  hapto stays separate (different physics).
        if kind == "hapto":
            sig = ("hapto", extra.get("eta") if extra else 0, merged)
        else:
            sig = ("kappa", len(subset), merged)
        if sig in seen:
            return
        seen.add(sig)
        m = {
            "donors": tuple(sorted(subset)),
            "kappa": len(subset),
            "ring_sizes": ring_sizes,
            "donor_elems": tuple(sym[i] for i in sorted(subset)),
            "kind": kind,
            "label": label or _label_for(kind, len(subset), subset, sym, extra),
        }
        if extra:
            m.update(extra)
        if _passes_realism_for_mode(mol, m, strict=use_strict):
            modes.append(m)

    # ORDER MATTERS: functional first so labelled (richer) modes subsume the
    # plain-kappa equivalents in the dedup signature.
    # 1) functional-group canonical modes
    if include_functional:
        for fg in _detect_functional_modes(mol):
            ring_sizes = (fg["ring_size"],) if fg["ring_size"] else ()
            add(fg["donors"], ring_sizes, "functional",
                label=fg["label"], extra={"bridging": fg.get("bridging", False),
                                          "smarts": fg.get("smarts")})

    # 2) κ1
    for d in donors:
        add((d,), (), "kappa")

    # 3) κ2..MAX_KAPPA: connected feasible chelates
    for k in range(2, min(len(donors), MAX_KAPPA) + 1):
        for subset in itertools.combinations(donors, k):
            feas_edges = [(a, b) for a, b in itertools.combinations(subset, 2)
                          if _feasible_pair(mol, a, b)]
            if not feas_edges:
                continue
            parent = {d: d for d in subset}

            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x

            for a, b in feas_edges:
                parent[find(a)] = find(b)
            if len({find(d) for d in subset}) != 1:
                continue
            ring_sizes = tuple(sorted(_ring_size(mol, a, b) for a, b in feas_edges))
            add(tuple(sorted(subset)), ring_sizes, "kappa")

    # 4) hapto-π modes
    if include_hapto:
        for hm in _detect_hapto_modes(mol):
            add(hm["donors"], (), "hapto",
                label=hm["label"],
                extra={"eta": hm["eta"], "ring_n": hm.get("ring_n", 0),
                       "smarts": hm.get("smarts")})

    # Final lex-sort for byte-identical deterministic output
    modes.sort(key=lambda m: (m["kind"], -m["kappa"], m.get("eta", 0),
                              m["donor_elems"], m["donors"], m["label"]))
    return modes


def _passes_realism_for_mode(mol, mode: Dict, *, strict: bool) -> bool:
    if not strict:
        return True
    kind = mode.get("kind")
    if kind == "hapto":
        if mode.get("eta", 0) > mode.get("ring_n", 0):
            return False
        return True
    if mode["kappa"] >= 2 and mode["ring_sizes"]:
        rs = min(mode["ring_sizes"])
        if rs < CHELATE_MIN:
            return False
        if rs > CHELATE_MAX:
            return False
    return True


def _label_for(kind: str, kappa: int, donors: Tuple[int, ...],
               sym: Dict[int, str], extra: Optional[Dict]) -> str:
    if kind == "hapto" and extra and "eta" in extra:
        return f"eta{extra['eta']}"
    el = "".join(sorted(sym[i] for i in donors))
    return f"kappa{kappa}-{el}"


def predicted_mode_count(smiles: str) -> int:
    """Return the predicted number of distinct realistic binding modes for the
    ligand `smiles`.  Used as the completeness-denominator."""
    return len(enumerate_modes(smiles))


# Convenience: dispatch into the existing linkage / hapto-mode helpers without
# loading them at import time (env-gate respected).
def merge_with_linkage_isomers(mol) -> List[Dict]:
    """Combine functional-group modes with the legacy linkage-isomer table
    (SCN/NCS/NO2/N3/DMSO/...).  Returns one merged mode-list."""
    out: List[Dict] = []
    try:
        from delfin.fffree.linkage_isomers import detect_ambidentate_groups
        groups = detect_ambidentate_groups(mol)
        for g in groups:
            for atom_idx, label in g["donor_options"]:
                out.append({
                    "label": f"linkage-{label}",
                    "donors": (atom_idx,),
                    "kappa": 1,
                    "ring_size": 0,
                    "kind": "linkage",
                    "donor_elems": (mol.GetAtomWithIdx(atom_idx).GetSymbol(),),
                })
    except Exception:
        pass
    return out


if __name__ == "__main__":
    tests = {
        "ethylenediamine (NCCN)": "NCCN",
        "acetate (CC(=O)[O-])": "CC(=O)[O-]",
        "glycinate": "[NH2]CC(=O)[O-]",
        "thiocyanate (ambidentate)": "[S-]C#N",
        "pyridine (monodentate)": "c1ccncc1",
        "2,2'-bipyridine": "c1ccc(-c2ccccn2)nc1",
        "terpyridine": "c1ccc(-c2cccc(-c3ccccn3)n2)nc1",
        "water": "O",
        "nitrate": "[O-][N+](=O)[O-]",
        "sulfate": "[O-][S](=O)(=O)[O-]",
        "phosphate": "[O-][P](=O)([O-])[O-]",
        "acac (β-diketonate, dep)": "[O-]C(C)=CC(C)=O",
        "Cp anion": "[cH-]1cccc1",
        "benzene": "c1ccccc1",
        "allyl": "C=CC",
    }
    for name, smi in tests.items():
        ms = enumerate_modes(smi)
        print(f"\n{name}: {len(ms)} modes")
        for m in ms:
            print(f"   [{m['kind']:<10}] {m['label']:<28} κ={m['kappa']} "
                  f"donors={m['donors']}")
