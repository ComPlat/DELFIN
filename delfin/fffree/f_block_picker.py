"""delfin.fffree.f_block_picker — f-block-aware polyhedron picker for CN 8-12.

Lanthanide and actinide chemistry differs from d-block in three load-bearing
ways that the existing first-candidate / chelate-aware picker cannot capture:

  1. **Ionic bonding.**  Ln³⁺ / An^n+ M–X interactions are ~70-90 % ionic
     (Shannon ionic-radius model fits within ±0.05 Å).  There is far less
     directional preference than in d-block covalent bonding, so the SAME
     CN can adopt several near-degenerate polyhedra (e.g. CN8 SAP-8 vs DD-8,
     CN9 TTP-9 vs CSAP-9) depending on donor identity and packing.

  2. **High CN dominant.**  CN 8-12 are typical (CN8 ≈ 35 %, CN9 ≈ 40 %,
     CN10-12 ≈ 20 % of structurally-characterised Ln complexes).  The naive
     ``geometries_for_cn(CN)[0]`` rule picks ONE polyhedron per CN and
     silently discards the alternative — wrong for the majority of Ln cases.

  3. **Multi-modal M–D distributions.**  Lanthanide contraction (La→Lu r
     shrinks by ~0.15 Å) PLUS multi-CN coexistence in the same crystal
     (mixed CN8/CN9 sites in Ln nitrates) yields characteristic D-D vertex
     distributions that are diagnostic of the polyhedron — capture this via
     a Mahalanobis match against the CCDC pair-bond library.

The existing module ``f_block_polyhedra`` provides the *vertex math* for
SAP-8, DD-8, TTP-9, CSAP-9, BICAP-10, CAP-11, IH-12, the metal-classification
``is_f_block``, and a Shannon-radii M–D fallback distance.  It does NOT pick
among multiple candidates — it returns the FIRST entry of ``FBLOCK_GEOM_BY_CN``
deterministically.  This module adds the picking layer.

Universal contract
------------------

* NO SMILES patterns.
* NO per-metal override tables — the only metal-aware quantity is the
  Pyykkö covalent radius (universal table for all Ln/An).
* The score is purely geometric: D-D pair distance distribution match
  against the CCDC pair-bond library, plus a symmetry-preference term
  for ionic bonding (high-symmetry polyhedra preferred at parity).
* All f-block detection is by atomic-number range (Z 57-71 Ln, Z 89-103 An).
* The picker is also defined for non-f-block high-CN metals (Y, Sc, Cd, Mo
  at CN 7-12) — these get the SAP-8 / DD-8 / TTP-9 / CSAP-9 / BICAP-10
  candidates added to their list since the chemistry is the same family
  (large early-/late-d cations with weakly directional bonding).

Env-flag::

    DELFIN_FFFREE_F_BLOCK_PICK = "1"   # default ON under MOGUL_PRIMARY=1
    DELFIN_FFFREE_F_BLOCK_PICK = "0"   # legacy first-candidate rule

Default-OFF byte-identical with HEAD when MOGUL_PRIMARY is off.

Author: hmaximilian <hmaximilian496@gmail.com>
Branch: feat-mogul-primary-2026-06-07
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np


__all__ = [
    "PYYKKO_COV_RADIUS",
    "HIGH_CN_NON_FBLOCK_CANDIDATES",
    "f_block_picker_enabled",
    "is_high_cn_non_fblock_target",
    "f_block_candidates_for_cn",
    "score_polyhedron_dd_distribution",
    "f_block_polyhedron_picker",
    "pyykko_delta",
]


# ---------------------------------------------------------------------------
# Pyykkö (2009) covalent radii — universal table for Δ(Pyykkö) computation
# ---------------------------------------------------------------------------


# Source: Pyykkö & Atsumi, Chem. Eur. J. 15 (2009) 186–197.  Single-bond
# covalent radii (Å) for all elements.  Used for the "Pyykkö-Δ" metric
# computed in the smoke 100 report (mean |r_obs - (r_M + r_X)|).
PYYKKO_COV_RADIUS: Dict[str, float] = {
    # Period 1-3 donors / non-metals
    "H": 0.32, "B": 0.85, "C": 0.75, "N": 0.71, "O": 0.63, "F": 0.64,
    "Si": 1.16, "P": 1.11, "S": 1.03, "Cl": 0.99,
    "As": 1.21, "Se": 1.16, "Br": 1.14,
    "Te": 1.35, "I": 1.33,
    # First-row TM (target subset)
    "Sc": 1.48, "Ti": 1.36, "V": 1.34, "Cr": 1.22, "Mn": 1.19, "Fe": 1.16,
    "Co": 1.11, "Ni": 1.10, "Cu": 1.12, "Zn": 1.18,
    # Second-row TM (target subset)
    "Y": 1.63, "Zr": 1.54, "Nb": 1.47, "Mo": 1.38, "Tc": 1.28, "Ru": 1.25,
    "Rh": 1.25, "Pd": 1.20, "Ag": 1.28, "Cd": 1.36,
    # Third-row TM
    "Hf": 1.52, "Ta": 1.46, "W": 1.37, "Re": 1.31, "Os": 1.29, "Ir": 1.22,
    "Pt": 1.23, "Au": 1.24, "Hg": 1.33,
    # Lanthanides (Z = 57-71)
    "La": 1.80, "Ce": 1.63, "Pr": 1.76, "Nd": 1.74, "Pm": 1.73, "Sm": 1.72,
    "Eu": 1.68, "Gd": 1.69, "Tb": 1.68, "Dy": 1.67, "Ho": 1.66, "Er": 1.65,
    "Tm": 1.64, "Yb": 1.70, "Lu": 1.62,
    # Actinides (Z = 89-103)
    "Ac": 1.86, "Th": 1.75, "Pa": 1.69, "U": 1.70, "Np": 1.71, "Pu": 1.72,
    "Am": 1.66, "Cm": 1.66, "Bk": 1.68, "Cf": 1.68, "Es": 1.65, "Fm": 1.67,
    "Md": 1.73, "No": 1.76, "Lr": 1.61,
}


# ---------------------------------------------------------------------------
# High-CN non-f-block target metals — get the same candidate set as f-block
# ---------------------------------------------------------------------------


# Universal: group-3 (Y, Sc) and late-d / heavy non-d-block at CN >= 8 share
# the f-block chemistry of large, weakly-directional cations.  This set is
# defined by atomic number ranges, NOT a per-metal hardcode:
#   * Z = 21 (Sc),  Z = 39 (Y),  Z = 57-71 (Ln, true f-block)
#   * Z = 89-103 (An, true f-block)
#   * Z = 40-44 (4d early: Zr, Nb, Mo, Tc, Ru) at CN >= 8 (rare but real,
#     e.g. [Mo(CN)8]4-, K4[Mo(CN)8])
#   * Z = 48 (Cd) at CN >= 7 (Cd-uracil mixed oxide-amine sites)
# The atomic-number table is universal — every metal whose CN >= 7 and Z in
# {21, 39-48 d-block early/heavy, 57-71 Ln, 72-79 5d, 89-103 An}.
_HIGH_CN_TARGET_Z: frozenset = frozenset(
    {21, 39}                                          # Sc, Y
    | set(range(40, 49))                              # Zr-Cd (4d-block + Cd)
    | set(range(57, 72))                              # Ln (Z 57-71)
    | set(range(72, 80))                              # Hf-Au (5d-block)
    | set(range(89, 104))                             # An (Z 89-103)
)

# Symbol → Z table for the metals in the picker scope (subset used above).
_SYM_TO_Z: Dict[str, int] = {
    "Sc": 21, "Y": 39,
    "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44,
    "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48,
    "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63,
    "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70,
    "Lu": 71,
    "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78,
    "Au": 79,
    "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95,
    "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100, "Md": 101, "No": 102,
    "Lr": 103,
}


# Canonical high-CN candidate set (CN -> ordered list of geometries),
# derived purely from ``f_block_polyhedra.FBLOCK_GEOM_BY_CN`` plus the legacy
# names already present in ``polyhedra.GEOM_BY_CN`` (SQAP-8 alias of SAP-8,
# PB-7).  The picker considers this entire set for any metal in scope.
HIGH_CN_NON_FBLOCK_CANDIDATES: Dict[int, List[str]] = {
    7: ["PB-7 pentagonal bipyramid"],
    8: [
        "SQAP-8 square antiprism",       # legacy SAP-8 alias
        "SAP-8 square antiprism",
        "DD-8 dodecahedron",
    ],
    9: [
        "TTP-9 tricapped trigonal prism",
        "CSAP-9 capped square antiprism",
    ],
    10: ["BICAP-10 bicapped square antiprism"],
    11: ["CAP-11 monocapped pentagonal antiprism"],
    12: ["IH-12 icosahedron"],
}


# ---------------------------------------------------------------------------
# Env-flag
# ---------------------------------------------------------------------------


def f_block_picker_enabled() -> bool:
    """True iff the f-block-aware picker is wired on.

    Default-ON when the Mogul-PRIMARY path is active so f-block / high-CN
    SMILES are routed through this picker by default.  The user can switch
    back to the naive first-candidate rule by setting
    ``DELFIN_FFFREE_F_BLOCK_PICK=0`` explicitly.

    Default-OFF byte-identical with HEAD when MOGUL_PRIMARY itself is off,
    since this picker is only consulted from the MOGUL_PRIMARY code paths
    (and via explicit force=True callers).
    """
    raw = os.environ.get("DELFIN_FFFREE_F_BLOCK_PICK", "").strip()
    if raw == "":
        # Default ON only when MOGUL_PRIMARY is active.
        return os.environ.get("DELFIN_FFFREE_MOGUL_PRIMARY", "0") == "1"
    return raw == "1"


# ---------------------------------------------------------------------------
# Metal classification helpers
# ---------------------------------------------------------------------------


def is_high_cn_non_fblock_target(metal_sym: str, cn: int) -> bool:
    """True if ``metal_sym`` at coordination number ``cn`` is in the picker's
    high-CN non-f-block target set.

    The picker covers:
      * any f-block metal (Z 57-71 Ln, Z 89-103 An) at CN 8-12
      * any group-3 / 4d-early / 4d-late / 5d-block / Cd metal at CN >= 7
        — these chemically resemble Ln in their large ionic radius and
        weak directional preference, so they benefit from the same
        multi-candidate Mahalanobis ranking.

    Universal: pure atomic-number range test, no SMILES patterns, no
    per-metal hardcodes.
    """
    z = _SYM_TO_Z.get(str(metal_sym))
    if z is None:
        return False
    if int(cn) < 7:
        return False
    return z in _HIGH_CN_TARGET_Z


def f_block_candidates_for_cn(metal_sym: str, cn: int) -> List[str]:
    """Return the candidate polyhedron list for ``(metal_sym, cn)``.

    Merges the legacy ``GEOM_BY_CN[cn]`` (so the canonical first-rule
    candidate is preserved at index 0 even when this picker is OFF) with
    the high-CN f-block / Wells canonical set
    (``HIGH_CN_NON_FBLOCK_CANDIDATES[cn]``).  Returns an empty list if
    ``cn`` is outside the supported range (7-12) or the metal is outside
    the picker scope.

    The returned list is deterministic and lex-stable (no platform-
    dependent dict ordering).
    """
    try:
        from delfin.fffree import polyhedra as _polyhedra
    except ImportError:
        return []
    legacy = list(_polyhedra.GEOM_BY_CN.get(int(cn), []))
    extra = list(HIGH_CN_NON_FBLOCK_CANDIDATES.get(int(cn), []))
    # Merge preserving order; legacy first so byte-identical when picker OFF.
    out: List[str] = list(legacy)
    for g in extra:
        if g not in out:
            out.append(g)
    return out


# ---------------------------------------------------------------------------
# Pyykkö-Δ utility (for the validation report and as a candidate score term)
# ---------------------------------------------------------------------------


def pyykko_delta(
    metal_sym: str,
    donor_syms: Sequence[str],
    md_distances: Sequence[float],
) -> Optional[float]:
    """Mean |r_obs - (r_M + r_X)| (Å) for an observed metal-donor shell.

    Returns ``None`` when ``metal_sym`` or any donor element is not in the
    Pyykkö table (so the caller can decide between skip / covalent-radii
    fallback / ionic-radii fallback).

    Universal — pure table lookup, no per-class branches.
    """
    r_m = PYYKKO_COV_RADIUS.get(str(metal_sym))
    if r_m is None:
        return None
    diffs: List[float] = []
    for sym, d in zip(donor_syms, md_distances):
        r_d = PYYKKO_COV_RADIUS.get(str(sym))
        if r_d is None:
            return None
        diffs.append(abs(float(d) - (r_m + r_d)))
    if not diffs:
        return None
    return float(np.mean(diffs))


# ---------------------------------------------------------------------------
# CCDC pair-bond library helper (lazy-loaded; cached)
# ---------------------------------------------------------------------------


_PAIRBOND_CACHE: Optional[Dict[Tuple[str, str], Tuple[float, float, int]]] = None


def _load_pairbond_library() -> Dict[Tuple[str, str], Tuple[float, float, int]]:
    """Return cached ``{(elem_a, elem_b): (mu, sigma, n)}`` from the GRIP lib.

    The GRIP TM library (``grip_lib_v5_TM.npz``) stores element-pair bond
    distance statistics keyed as JSON strings
    ``["A", hyb_a, "B", hyb_b]`` — we collapse over hybridisation by taking
    the n-weighted mean / pooled-σ across all matching (A, B) entries.

    Returns an empty dict if the library is missing or unparseable, so the
    score term degrades gracefully to "no distribution match" (uniform).
    """
    global _PAIRBOND_CACHE
    if _PAIRBOND_CACHE is not None:
        return _PAIRBOND_CACHE
    out: Dict[Tuple[str, str], List[Tuple[float, float, int]]] = {}
    lib_path = os.environ.get(
        "DELFIN_GRIP_LIB_PATH",
        "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5_TM.npz",
    )
    try:
        import json
        lib = np.load(lib_path, allow_pickle=True)
        pkeys = lib["pair_bond_keys"]
        mu = lib["pair_bond_mu"]
        sg = lib["pair_bond_sigma"]
        nn = lib["pair_bond_n"]
        for i, k in enumerate(pkeys):
            try:
                parsed = json.loads(str(k))
                if not isinstance(parsed, list) or len(parsed) < 3:
                    continue
                a = str(parsed[0])
                b = str(parsed[2])
            except Exception:
                continue
            key = tuple(sorted([a, b]))
            out.setdefault(key, []).append((float(mu[i]), float(sg[i]), int(nn[i])))
    except Exception:
        _PAIRBOND_CACHE = {}
        return _PAIRBOND_CACHE
    collapsed: Dict[Tuple[str, str], Tuple[float, float, int]] = {}
    for key, entries in out.items():
        total_n = sum(e[2] for e in entries)
        if total_n <= 0:
            continue
        mean_mu = sum(e[0] * e[2] for e in entries) / total_n
        # Pooled σ² = (Σ n_i (σ_i² + μ_i²) ) / Σ n_i  - mean_mu²
        var_pool = (
            sum(e[2] * (e[1] * e[1] + e[0] * e[0]) for e in entries) / total_n
            - mean_mu * mean_mu
        )
        sigma_pool = math.sqrt(max(1e-6, var_pool))
        collapsed[key] = (float(mean_mu), float(sigma_pool), int(total_n))
    _PAIRBOND_CACHE = collapsed
    return _PAIRBOND_CACHE


def _expected_dd_distance(donor_a: str, donor_b: str) -> Optional[Tuple[float, float]]:
    """Return ``(mu, sigma)`` for a donor-donor pair from the CCDC library,
    or ``None`` if not tabulated."""
    if not donor_a or not donor_b:
        return None
    lib = _load_pairbond_library()
    key = tuple(sorted([donor_a, donor_b]))
    e = lib.get(key)
    if e is None:
        return None
    return (float(e[0]), float(e[1]))


# ---------------------------------------------------------------------------
# Polyhedron scoring — D-D distribution match + symmetry preference
# ---------------------------------------------------------------------------


# Approximate "high-symmetry preference" weights for ionic-bonded cations:
# higher = more preferred at parity.  Pure heuristic — these numbers ONLY
# break score ties (the dominant term is the empirical D-D match).
_SYMMETRY_WEIGHT: Dict[str, float] = {
    # CN 7
    "PB-7 pentagonal bipyramid": 1.00,          # D5h
    # CN 8
    "SQAP-8 square antiprism": 0.95,            # legacy SAP-8 alias
    "SAP-8 square antiprism": 1.00,             # D4d
    "DD-8 dodecahedron": 0.90,                  # D2d
    # CN 9
    "TTP-9 tricapped trigonal prism": 1.00,     # D3h
    "CSAP-9 capped square antiprism": 0.85,     # C4v
    # CN 10
    "BICAP-10 bicapped square antiprism": 1.00,  # D4d
    "CSAP-10 capped square antiprism": 0.85,
    "SAP-10 pentagonal antiprism": 0.90,
    # CN 11
    "CAP-11 monocapped pentagonal antiprism": 1.00,  # C5v
    # CN 12
    "IH-12 icosahedron": 1.00,                  # Ih
}


def _polyhedron_dd_distances(
    geometry: str,
    cn: int,
    md_targets: Sequence[float],
) -> Optional[np.ndarray]:
    """Compute the D-D pair distances (Å) of ``geometry`` placed at the
    metal with each donor at its ideal M-D distance.

    Returns an (n*(n-1)/2,) array ordered by (i, j) with i < j, or
    ``None`` if the geometry cannot be loaded.

    Universal: pure geometric — vertex unit vectors * per-donor M-D yield
    the donor positions; pairwise distances are Euclidean.
    """
    try:
        from delfin.fffree import polyhedra as _polyhedra
    except ImportError:
        return None
    try:
        V = _polyhedra.ref_vectors(str(geometry))
    except (KeyError, Exception):
        return None
    if V is None:
        return None
    V = np.asarray(V, dtype=float)
    if V.ndim != 2 or V.shape[1] != 3 or V.shape[0] < int(cn):
        return None
    V = V[: int(cn)]
    norms = np.linalg.norm(V, axis=1, keepdims=True)
    norms = np.where(norms < 1e-9, 1.0, norms)
    Vu = V / norms
    md = np.asarray(md_targets, dtype=float).reshape(-1, 1)
    if len(md) != int(cn):
        return None
    P = Vu * md  # donor positions in metal frame
    out: List[float] = []
    for i in range(int(cn)):
        for j in range(i + 1, int(cn)):
            out.append(float(np.linalg.norm(P[i] - P[j])))
    return np.asarray(out, dtype=float)


def score_polyhedron_dd_distribution(
    geometry: str,
    metal_sym: str,
    donor_syms: Sequence[str],
    md_targets: Optional[Sequence[float]] = None,
) -> Optional[float]:
    """Score (lower = better) of ``geometry`` for the given donor set.

    Score = mean squared (z-distance) over all donor-donor pairs, where
    z-distance = |d_pred(geom) - mu_CCDC| / sigma_CCDC.  ``mu_CCDC`` and
    ``sigma_CCDC`` are the empirical pair-bond statistics from the GRIP
    library (collapsed over hybridisation).  Pairs not tabulated in the
    CCDC library are skipped.

    The symmetry-preference weight is folded in as a small multiplicative
    bonus (score / weight), so at parity (equal Mahalanobis sum) the
    higher-symmetry polyhedron wins — matches the ionic-bonding chemistry
    of Ln / heavy non-d cations.

    Parameters
    ----------
    geometry
        Polyhedron canonical name (e.g. ``"SAP-8 square antiprism"``).
    metal_sym
        Metal symbol (used to compute M-D targets if not supplied).
    donor_syms
        Per-vertex donor element symbols (``len = cn``).
    md_targets
        Optional per-donor target M-D distances (Å).  When ``None``,
        we compute them via :func:`polyhedra.md_distance`.

    Returns
    -------
    float or None
        Lower-is-better score.  ``None`` on degenerate input (geometry
        not loadable, donor list empty, or NO pair is tabulated in CCDC
        — in that case the caller should fall back to symmetry-only).
    """
    cn = int(len(donor_syms))
    if cn < 2:
        return None
    if md_targets is None:
        try:
            from delfin.fffree import polyhedra as _polyhedra
        except ImportError:
            return None
        md_targets = [
            float(_polyhedra.md_distance(str(metal_sym), str(s)))
            for s in donor_syms
        ]
    pred = _polyhedron_dd_distances(geometry, cn, md_targets)
    if pred is None:
        return None

    # Pair index → (donor_a, donor_b) lookup, same (i, j) order as pred.
    idx_pairs: List[Tuple[str, str]] = []
    for i in range(cn):
        for j in range(i + 1, cn):
            idx_pairs.append((str(donor_syms[i]), str(donor_syms[j])))

    z2_sum = 0.0
    n_tab = 0
    for k, (a, b) in enumerate(idx_pairs):
        ms = _expected_dd_distance(a, b)
        if ms is None:
            continue
        mu, sigma = ms
        sig = max(0.05, float(sigma))  # floor σ to avoid blow-up
        z = (float(pred[k]) - float(mu)) / sig
        z2_sum += z * z
        n_tab += 1

    if n_tab == 0:
        return None
    base = z2_sum / float(n_tab)
    weight = max(0.5, float(_SYMMETRY_WEIGHT.get(str(geometry), 0.80)))
    return float(base) / float(weight)


# ---------------------------------------------------------------------------
# Main picker entry point
# ---------------------------------------------------------------------------


def f_block_polyhedron_picker(
    metal_sym: str,
    cn: int,
    donors: Sequence[str],
    md_targets: Optional[Sequence[float]] = None,
    candidates: Optional[Sequence[str]] = None,
    force: bool = False,
) -> Optional[str]:
    """f-block-aware polyhedron selector for CN 7-12.

    Enumerates all candidate polyhedra for ``(metal_sym, cn)``, scores each
    by the donor-donor distance distribution match against the CCDC
    pair-bond library (with a high-symmetry preference weight for ionic
    bonding), and returns the lowest-score name.

    Universal contract
    ------------------
    * NO SMILES patterns.
    * No per-metal hardcoded preference — only the universal Pyykkö table
      (for the validation report) and the per-element CCDC pair-bond
      statistics drive the score.
    * Falls back gracefully to the first candidate when:
       - the picker env-flag is OFF and ``force`` is False,
       - the CN is outside 7-12,
       - the metal is outside the picker scope,
       - no candidate scores successfully (e.g. all pair-bond keys absent).

    Parameters
    ----------
    metal_sym : str
        Metal element symbol (e.g. ``"La"``, ``"Y"``, ``"U"``).
    cn : int
        Coordination number.
    donors : sequence of str
        Per-vertex donor element symbols (length = ``cn``).
    md_targets : sequence of float, optional
        Per-donor M-D target distances (Å).  When ``None`` we compute them
        via :func:`polyhedra.md_distance`, which already routes f-block
        metals through the Shannon-radii table when the gate is on.
    candidates : sequence of str, optional
        Explicit candidate list.  When ``None`` we use
        :func:`f_block_candidates_for_cn(metal_sym, cn)`.
    force : bool, default False
        Bypass the env-flag gate.  When True, the picker runs even if
        ``DELFIN_FFFREE_F_BLOCK_PICK`` is unset — useful for tests and for
        the chelate-aware picker which delegates to us conditionally.

    Returns
    -------
    str or None
        The chosen polyhedron canonical name.  ``None`` if no candidate
        could be loaded at all (very rare; only for unsupported CN).
    """
    cn = int(cn)
    cand_list: List[str]
    if candidates is not None:
        cand_list = list(candidates)
    else:
        cand_list = f_block_candidates_for_cn(metal_sym, cn)
    if not cand_list:
        return None

    # Default-OFF byte-identical: when the gate is OFF AND not forced,
    # return the first candidate (this is the legacy first-rule).
    if not force and not f_block_picker_enabled():
        return cand_list[0]

    # Out-of-scope metals fall back to the first candidate (legacy).
    if not is_high_cn_non_fblock_target(metal_sym, cn):
        # is_high_cn_non_fblock_target returns False for e.g. Cu/Fe at CN<=6
        # AND for metals not in our atomic-number ranges.  Check if metal
        # is f-block specifically — if it is, we still proceed (the chemistry
        # match is exact); if it isn't AND CN<7, return legacy.
        try:
            from delfin.fffree import f_block_polyhedra as _FBP
            if not _FBP.is_f_block(metal_sym):
                return cand_list[0]
        except ImportError:
            return cand_list[0]

    # Score every candidate.  Ties broken by candidate input order (i.e.
    # the legacy first-candidate rule is the deterministic tie-break).
    best_name = cand_list[0]
    best_score = float("inf")
    any_scored = False
    for name in cand_list:
        score = score_polyhedron_dd_distribution(
            geometry=str(name),
            metal_sym=str(metal_sym),
            donor_syms=donors,
            md_targets=md_targets,
        )
        if score is None:
            # No CCDC pair-bond coverage for ANY donor-donor pair under
            # this geometry — fall through to symmetry-only ranking
            # (apply a symmetry-only penalty so a higher-sym candidate
            # wins at parity).
            weight = float(_SYMMETRY_WEIGHT.get(str(name), 0.80))
            score = 100.0 / max(0.5, weight)
        any_scored = True
        # Strict less-than for deterministic tie-break (input order wins).
        if score < best_score - 1e-9:
            best_score = score
            best_name = name
    if not any_scored:
        return cand_list[0]
    return best_name
