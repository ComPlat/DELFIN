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
    for fi, (syms, P) in enumerate(frames):
        try:
            syms, P = _maybe_relax(syms, P)
            if not _build_is_clean(syms, P, cn=cn, geom=geom, donors=donors):
                continue
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
    for syms, P in frames:
        try:
            syms, P = _maybe_relax(syms, P)
            if _build_is_clean(syms, P, cn=cn, geom=geom, donors=donors):
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

    ``exempt_pairs`` (optional set of (min,max) global index pairs): heavy-heavy
    bonds whose SHORT length is chemically correct (genuine triple/multiple bonds
    such as a C≡O carbonyl ~1.12 A or a C≡N nitrile), so they are NOT counted as
    collapsed.  Default None = byte-identical to the historic gate."""
    if os.environ.get("DELFIN_FFFREE_SELFGATE", "1") == "0":
        return True
    P = np.asarray(P, dtype=float)
    if P.size == 0 or not np.all(np.isfinite(P)):
        return False
    syms = list(syms)
    bonds = _bd._geometric_bonds(syms, P)
    if exempt_pairs:
        _ex = {(min(i, j), max(i, j)) for i, j in exempt_pairs}
        n_coll = sum(1 for i, j in bonds
                     if not (_bd._is_metal(syms[i]) or _bd._is_metal(syms[j]))
                     and (min(i, j), max(i, j)) not in _ex
                     and float(np.linalg.norm(P[i] - P[j])) < 0.82 * _bd._ideal_bond(syms[i], syms[j]))
        if n_coll > 0:
            return False
    elif _bd._count_collapsed(syms, P, bonds) > 0:      # any collapsed heavy bond
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
                                   donors=sdonors):
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
                                       donors=donors):   # per-frame self-gate
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
                               donors=donors):   # donor-aware self-gate
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


def _fffree_isomers(smiles: str, max_isomers: int = 50
                    ) -> Optional[List[Tuple[str, str]]]:
    d = DEC.decompose(smiles)
    if d is None:
        return None
    geom_key = _GEOM_TO_POLYA.get(d["geometry"])
    if geom_key is None or geom_key not in PIC._GROUPS:
        return None
    if d.get("has_eta"):
        return _fffree_hapto_isomers(d, max_isomers)
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
                if not _build_is_clean(syms, P, cn=d.get("cn"), geom=d.get("geometry")):
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
        if not _build_is_clean(syms, P, cn=d.get("cn"), geom=d.get("geometry")):   # self-gate: destroyed/over-coord/shape-outlier
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
