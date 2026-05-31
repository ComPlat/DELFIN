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
import delfin._bond_decollapse as _bd


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


def _build_is_clean(syms, P, cn=None, geom=None, donors=None) -> bool:
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
        })
    try:
        configs = PIC.enumerate_chelate_configs(geom_key, specs)
    except Exception:
        return None
    if not configs:
        return None
    geom_tag = d["geometry"].split()[0]
    results = []
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
        if not _build_is_clean(syms, P, cn=d.get("cn"), geom=d.get("geometry"),
                               donors=donors):   # donor-aware self-gate -> skip config
            continue
        results.append((_xyz(syms, P), f"{geom_tag}-chelate-{k+1}"))
    return results or None


def _fffree_isomers(smiles: str, max_isomers: int = 50
                    ) -> Optional[List[Tuple[str, str]]]:
    d = DEC.decompose(smiles)
    if d is None:
        return None
    geom_key = _GEOM_TO_POLYA.get(d["geometry"])
    if geom_key is None or geom_key not in PIC._GROUPS:
        return None
    if d.get("has_chelate"):
        return _fffree_chelate_isomers(d, geom_key, max_isomers)
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
    for k, coloring in enumerate(colorings[:max_isomers]):
        vertex_specs = [lig_ref[lab] for lab in coloring]
        try:
            built = AC.assemble_heteroleptic_from_mols(d["metal"], d["geometry"], vertex_specs)
        except Exception:
            return None
        if built is None:
            return None
        syms, P = built
        syms, P = _maybe_relax(syms, P)
        if not _build_is_clean(syms, P, cn=d.get("cn"), geom=d.get("geometry")):   # self-gate: destroyed/over-coord/shape-outlier -> legacy
            return None
        vertex_elems = [lab_elem[lab] for lab in coloring]
        name = _classify_coloring(geom_key, vertex_elems)
        geom_tag = d["geometry"].split()[0]
        label = f"{name}-{geom_tag}-{k+1}" if name else f"{geom_tag}-{k+1}"
        results.append((_xyz(syms, P), label))
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
