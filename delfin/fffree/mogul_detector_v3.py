"""Mogul-style per-feature anomaly detector v3 (CCDC-calibrated).

This module is the TOP-1 calibration fix produced by the CCDC Mogul cross-
validation pass (``scripts/ccdc_mogul_validation.py``,
``iters/CCDC_MOGUL_DETECTOR_CALIBRATION_2026_06_04.md``).

Motivation
----------
Our legacy aggregator
``agent_workspace/quality_framework/scripts/metric_cod_fragment_lookup.py``
emits anomaly records **per central atom**: each record describes one heavy
atom whose fragment-context matches a CSD-derived L2/L1/L0 fragment bucket
and whose bond/angle/improper channels deviate above a configurable MAD
threshold. The strict thresholds (``sev_mad >= 5`` and a ``0.12 Å`` cap on
bonds, ``20 deg`` cap on angles) give very high precision but poor recall
relative to CCDC's authoritative Mogul (``ccdc.conformer.GeometryAnalyser``,
``|z| >= 2``):

* bond recall ~13 %
* angle recall ~4 %
* torsion recall 0 %  (channel not emitted at all)

The detector v3 fixes the **emission semantics** without changing any
default behaviour:

1. **Per-feature records** — each anomaly carries the *full participating
   atom set* (both atoms of a bond, the three atoms of an angle, the four
   atoms of a torsion) so downstream calibration / comparison tools can
   attribute hits to the correct CCDC feature rather than guessing from the
   central atom only.
2. **CCDC-calibrated soft threshold** — when ``mode='ccdc'`` we drop the
   geometric cap and lower the MAD threshold to ``2.5`` (≈ CCDC's
   ``|z| >= 2``) so the per-feature recall matches the published Mogul
   grade.
3. **Torsion channel** — a new ``per_torsion`` mode walks each
   ``i-j-k-l`` quadruplet and computes the deviation against the **same
   COD index** that already feeds bond and angle channels. The pooled-MAD
   statistics from neighbour records satisfy the torsion comparison; we
   cap MAD with ``5 deg`` to guard against degenerate distributions.

Default state
-------------
The module is **default-OFF**: the legacy aggregator is not modified, and
``DELFIN_MOGUL_V3_DETECTOR`` (and the function-level ``mode`` parameter)
must be set to ``"ccdc"`` to use the new emission semantics. Byte-identity
of the legacy aggregator output is preserved.

Env flag
--------
``DELFIN_MOGUL_V3_DETECTOR=1`` selects the CCDC-calibrated per-feature
emission when callers use :func:`detect_anomalies_v3` with default
``mode=None``.

Library
-------
The detector reads the same fragment index as
``metric_cod_fragment_lookup`` (``COD_FRAGMENT_INDEX`` env, default
``cod_fragment_index_v3.json``). This module deliberately re-implements the
fragment walk so it is a stand-alone, testable, dependency-free unit;
identical statistical primitives (``_mad``, ``_classify_bond_order``,
``_hybridization``) are shared with the legacy script via direct import,
falling back to local copies when the QF scripts directory is not on the
path.
"""
from __future__ import annotations

import json
import math
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "detect_anomalies_v3",
    "DEFAULT_INDEX_PATH",
    "MOGUL_V3_ENV",
    "CCDC_MAD_THRESHOLD",
    "LEGACY_MAD_THRESHOLD",
]

# ----------------------------------------------------------------- #
#  Shared with metric_cod_fragment_lookup
# ----------------------------------------------------------------- #
_QF_SCRIPTS = "/home/qmchem_max/agent_workspace/quality_framework/scripts"
if _QF_SCRIPTS not in sys.path:
    sys.path.insert(0, _QF_SCRIPTS)

# These imports are pure utilities used by the legacy aggregator. They live
# in QF/scripts; if QF is not on disk we provide local fallbacks so the
# detector still imports (tests use these fallbacks).
try:
    import metric_cod_fragment_lookup as _MCFL  # type: ignore

    _legacy_load_index = _MCFL._load_index
    _legacy_classify_bo = _MCFL._classify_bond_order
    _legacy_hybridization = _MCFL._hybridization
    _legacy_signed_dihedral = _MCFL._signed_dihedral
    _legacy_atom_hyb = _MCFL._atom_hyb
    _legacy_l1_pool = _MCFL._l1_pool
    _legacy_l0_pool = _MCFL._l0_pool
    _legacy_mad = _MCFL._mad
    _MIN_MATCHES = _MCFL.MIN_MATCHES
    _MAX_NBRS = _MCFL.MAX_NBRS_FOR_KEY
    _MAX_2ND = _MCFL.MAX_SECOND_SPHERE
    _LEGACY_BOND_DELTA_CAP = _MCFL.BOND_DELTA_CAP_ANG
    _LEGACY_ANGLE_DELTA_CAP = _MCFL.ANGLE_DELTA_CAP_DEG
except Exception:  # pragma: no cover - fallback path for fresh CI
    _MCFL = None  # type: ignore
    _legacy_load_index = lambda: {}  # type: ignore
    _legacy_classify_bo = lambda d, a, b: 1  # type: ignore
    _legacy_hybridization = lambda *args, **kw: "sp3"  # type: ignore

    def _legacy_signed_dihedral(p0, p1, p2, p3):  # type: ignore
        b1 = np.asarray(p1) - np.asarray(p0)
        b2 = np.asarray(p2) - np.asarray(p1)
        b3 = np.asarray(p3) - np.asarray(p2)
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)
        m = np.cross(n1, b2 / max(np.linalg.norm(b2), 1e-12))
        x = float(np.dot(n1, n2))
        y = float(np.dot(m, n2))
        return math.degrees(math.atan2(y, x))

    def _legacy_atom_hyb(syms, P, adj, a):  # type: ignore
        return 0, 109.5, "sp3"

    def _legacy_l1_pool(idx, prefix):  # type: ignore
        return []

    def _legacy_l0_pool(idx, elem, hyb):  # type: ignore
        return []

    def _legacy_mad(arr):  # type: ignore
        arr = np.asarray(arr, dtype=float)
        if not len(arr):
            return 0.0, 0.0
        med = float(np.median(arr))
        return float(np.median(np.abs(arr - med))), med

    _MIN_MATCHES = 15
    _MAX_NBRS = 5
    _MAX_2ND = 12
    _LEGACY_BOND_DELTA_CAP = 0.12
    _LEGACY_ANGLE_DELTA_CAP = 20.0

# bond_decollapse helpers (parse XYZ, metal test, covalent radii)
sys.path.insert(0, "/home/qmchem_max/ComPlat/DELFIN")
import delfin._bond_decollapse as _bd  # noqa: E402

# Public constants ------------------------------------------------- #
DEFAULT_INDEX_PATH = os.environ.get(
    "COD_FRAGMENT_INDEX",
    "/home/qmchem_max/agent_workspace/quality_framework/reports/cod_fragment_index_v3.json",
)
MOGUL_V3_ENV = "DELFIN_MOGUL_V3_DETECTOR"

# CCDC Mogul flags features with |z| >= 2; equivalent MAD severity is ~2.5
# (a normal distribution has 1 sigma = 1.4826 * MAD).
CCDC_MAD_THRESHOLD = 2.5
# Legacy threshold (high precision, low recall) — keeps the existing
# behaviour byte-identical when ``mode='legacy'``.
LEGACY_MAD_THRESHOLD = 5.0
# Torsion MAD-floor to guard against degenerate distributions (rings give
# zero spread).
TORSION_MAD_FLOOR_DEG = 5.0


# ----------------------------------------------------------------- #
#  Adjacency + per-atom local context (mirrors legacy)
# ----------------------------------------------------------------- #
def _build_adjacency(syms: Sequence[str], P: np.ndarray):
    """Return ``(adj, h_count, hyb_cache)`` as the legacy walker would.

    ``adj``:    heavy-atom-only adjacency, dict {a: [(b, d), ...]}
    ``h_count``: number of bonded hydrogens per heavy atom
    ``hyb_cache``: per heavy atom hybridisation label
    """
    from delfin.fffree.polyhedra import COV

    def _r(s):
        return COV.get(s, _bd._COV.get(s, 0.9))

    n = len(syms)
    adj: Dict[int, List[Tuple[int, float]]] = {
        i: [] for i in range(n) if syms[i] != "H" and not _bd._is_metal(syms[i])
    }
    for a in adj:
        for b in adj:
            if a < b:
                d = float(np.linalg.norm(P[a] - P[b]))
                if d < 1.30 * (_r(syms[a]) + _r(syms[b])):
                    adj[a].append((b, d))
                    adj[b].append((a, d))
    h_count: Dict[int, int] = {}
    for h in range(n):
        if syms[h] != "H":
            continue
        best = (1e9, -1)
        for k in range(n):
            if k == h or syms[k] == "H" or _bd._is_metal(syms[k]):
                continue
            d = float(np.linalg.norm(P[h] - P[k]))
            if d < best[0]:
                best = (d, k)
        if best[1] >= 0 and best[0] < 1.30 * (_r("H") + _r(syms[best[1]])):
            h_count[best[1]] = h_count.get(best[1], 0) + 1
    hyb_cache: Dict[int, str] = {}
    for a in adj:
        _, _, hy = _legacy_atom_hyb(syms, P, adj, a)
        hyb_cache[a] = hy
    return adj, h_count, hyb_cache


# ----------------------------------------------------------------- #
#  Local fragment lookup that mirrors the legacy cascade
# ----------------------------------------------------------------- #
def _fragment_lookup(idx, syms, a, nbr_full, hyb, hyb_cache, adj):
    """Mirrors ``metric_cod_fragment_lookup.cod_fragment_lookup_violations``
    cascade. Returns ``(cod_records, level)`` or ``([], None)``."""
    nbr_desc = [[t[0], t[1], t[2], t[3]] for t in nbr_full]
    second = []
    seen = {a} | {t[5] for t in nbr_full}
    for t in nbr_full:
        b = t[5]
        for c, d_nc in adj.get(b, []):
            if c in seen:
                continue
            bo = _legacy_classify_bo(d_nc, syms[b], syms[c])
            second.append((syms[c], bo, hyb_cache.get(c, "sp3")))
    second.sort()
    second_capped = second[:_MAX_2ND]
    l2_key = json.dumps(
        [syms[a], hyb, [list(d) for d in nbr_desc], [list(s2) for s2 in second_capped]],
        separators=(",", ":"),
    )
    cod = idx.get(l2_key, [])
    level = "L2"
    if len(cod) < _MIN_MATCHES:
        l1_prefix = [syms[a], hyb, [list(d) for d in nbr_desc]]
        cod = _legacy_l1_pool(idx, l1_prefix)
        level = "L1"
    if len(cod) < _MIN_MATCHES:
        cod = _legacy_l0_pool(idx, syms[a], hyb)
        level = "L0"
    if len(cod) < _MIN_MATCHES:
        return [], None
    return cod, level


# ----------------------------------------------------------------- #
#  Threshold helpers
# ----------------------------------------------------------------- #
def _is_unusual_bond(obs, med, mad, threshold, *, apply_cap):
    """Return ``(bool, severity)`` for the bond test."""
    mad = max(float(mad), 0.01)
    delta = abs(obs - med)
    sev = delta / mad
    if apply_cap:
        return (sev > threshold and delta > _LEGACY_BOND_DELTA_CAP), sev
    return (sev > threshold), sev


def _is_unusual_angle(obs, med, mad, threshold, *, apply_cap):
    mad = max(float(mad), 1.0)
    delta = abs(obs - med)
    sev = delta / mad
    if apply_cap:
        return (sev > threshold and delta > _LEGACY_ANGLE_DELTA_CAP), sev
    return (sev > threshold), sev


# ----------------------------------------------------------------- #
#  Public API
# ----------------------------------------------------------------- #
def detect_anomalies_v3(
    syms: Sequence[str],
    P,
    *,
    mode: Optional[str] = None,
    index: Optional[Dict] = None,
) -> List[Dict]:
    """Return per-feature anomaly records for ``(syms, P)``.

    Parameters
    ----------
    syms, P
        Atom symbols and Nx3 coordinates (Å).
    mode
        ``"legacy"``: byte-identical reproduction of the legacy
        ``cod_fragment_lookup_violations`` thresholds (sev >= 5, geometric
        cap). Default when ``DELFIN_MOGUL_V3_DETECTOR`` is unset.
        ``"ccdc"``: CCDC-calibrated thresholds (sev >= 2.5, no geometric
        cap), torsion channel emitted. Selected when
        ``DELFIN_MOGUL_V3_DETECTOR=1`` or passed explicitly.
    index
        Pre-loaded fragment index. If ``None`` the module-level cache
        from ``metric_cod_fragment_lookup`` is reused.

    Returns
    -------
    list of dict
        Each dict has keys:
        * ``axis``       — ``"bond"``, ``"angle"``, or ``"torsion"``
        * ``atoms``      — list of participating atom indices (sorted-set)
        * ``atom_syms``  — element symbols at those indices
        * ``center``     — central atom for the fragment lookup
        * ``obs``        — observed value (Å or deg)
        * ``median``     — fragment median value
        * ``mad``        — fragment MAD
        * ``sev_mad``    — severity (delta / MAD)
        * ``level``      — fragment lookup level ("L0", "L1", "L2")
        * ``n_cod``      — number of matching CSD records
        * ``threshold``  — MAD threshold applied
        * ``mode``       — selected mode

    Records are deterministic: ordered by ``(axis, sorted atoms)``.
    """
    if mode is None:
        mode = "ccdc" if os.environ.get(MOGUL_V3_ENV, "") == "1" else "legacy"
    if mode not in ("legacy", "ccdc"):
        raise ValueError(f"unknown mode {mode!r}")
    threshold = LEGACY_MAD_THRESHOLD if mode == "legacy" else CCDC_MAD_THRESHOLD
    apply_cap = mode == "legacy"
    emit_torsion = mode == "ccdc"

    P = np.asarray(P, dtype=float)
    n = len(syms)
    if n < 3:
        return []

    idx = index if index is not None else _legacy_load_index()
    if not idx:
        return []

    adj, h_count, hyb_cache = _build_adjacency(syms, P)

    out: List[Dict] = []

    # Track which (axis, atom-set) tuples we have already emitted so a
    # bond/angle anomaly is reported exactly once even though both end
    # atoms host the same fragment-lookup attempt.
    emitted: set = set()

    for a in adj:
        nbrs = adj[a][: _MAX_NBRS * 2]
        if not nbrs:
            continue
        nbr_full = [
            (
                syms[b],
                _legacy_classify_bo(d, syms[a], syms[b]),
                int(h_count.get(b, 0)),
                hyb_cache.get(b, "sp3"),
                float(d),
                int(b),
            )
            for b, d in nbrs
        ]
        nbr_full.sort(key=lambda t: (t[0], t[1], t[2], t[3]))
        nbr_full = nbr_full[: _MAX_NBRS]
        my_bonds = [t[4] for t in nbr_full]

        # 3-ring detection (legacy semantics).
        in_3ring_pair = set()
        for i in range(len(nbr_full)):
            bi = nbr_full[i][5]
            bi_nbrs = {bn for (bn, _) in adj.get(bi, [])}
            for j in range(i + 1, len(nbr_full)):
                bj = nbr_full[j][5]
                if bj in bi_nbrs:
                    in_3ring_pair.add((i, j))

        # Angle list with per-pair atom mapping.
        angle_records: List[Tuple[int, int, float]] = []
        angles = []
        for i in range(len(nbr_full)):
            for j in range(i + 1, len(nbr_full)):
                if (i, j) in in_3ring_pair:
                    angles.append(None)
                    continue
                v1 = P[nbr_full[i][5]] - P[a]
                v2 = P[nbr_full[j][5]] - P[a]
                n1 = float(np.linalg.norm(v1))
                n2 = float(np.linalg.norm(v2))
                if n1 < 1e-6 or n2 < 1e-6:
                    angles.append(None)
                    continue
                cosv = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
                ang = float(math.degrees(math.acos(cosv)))
                angles.append(ang)
                angle_records.append((nbr_full[i][5], nbr_full[j][5], ang))
        avg_ang = (
            sum(x for x in angles if x is not None)
            / max(1, sum(1 for x in angles if x is not None))
            if angles
            else 109.5
        )
        hyb = _legacy_hybridization(syms[a], len(nbr_full), avg_ang)

        # Fragment lookup
        cod, level = _fragment_lookup(idx, syms, a, nbr_full, hyb, hyb_cache, adj)
        if not cod:
            continue

        n_pos = min(len(my_bonds), _MAX_NBRS)
        # --- bonds ------------------------------------------------ #
        if level == "L2":
            # position-aware
            for pi in range(n_pos):
                cod_pos = [
                    r.get("bonds", [None] * n_pos)[pi]
                    for r in cod
                    if r.get("bonds") and len(r["bonds"]) > pi
                ]
                cod_pos = [v for v in cod_pos if v is not None]
                if len(cod_pos) < _MIN_MATCHES:
                    continue
                mad_b, med_b = _legacy_mad(cod_pos)
                unusual, sev = _is_unusual_bond(
                    my_bonds[pi], med_b, mad_b, threshold, apply_cap=apply_cap
                )
                if not unusual:
                    continue
                partner = nbr_full[pi][5]
                atom_set = tuple(sorted((a, partner)))
                key = ("bond", atom_set)
                if key in emitted:
                    continue
                emitted.add(key)
                out.append(
                    {
                        "axis": "bond",
                        "atoms": list(atom_set),
                        "atom_syms": [syms[atom_set[0]], syms[atom_set[1]]],
                        "center": a,
                        "obs": round(float(my_bonds[pi]), 3),
                        "median": round(float(med_b), 3),
                        "mad": round(float(max(mad_b, 0.01)), 4),
                        "sev_mad": round(float(sev), 2),
                        "level": level,
                        "n_cod": len(cod_pos),
                        "threshold": threshold,
                        "mode": mode,
                    }
                )
        else:
            cod_bonds_pool = []
            for r in cod:
                cod_bonds_pool.extend(r.get("bonds", []))
            if len(cod_bonds_pool) >= _MIN_MATCHES:
                mad_b, med_b = _legacy_mad(cod_bonds_pool)
                for pi, d in enumerate(my_bonds):
                    unusual, sev = _is_unusual_bond(
                        d, med_b, mad_b, threshold, apply_cap=apply_cap
                    )
                    if not unusual:
                        continue
                    partner = nbr_full[pi][5]
                    atom_set = tuple(sorted((a, partner)))
                    key = ("bond", atom_set)
                    if key in emitted:
                        continue
                    emitted.add(key)
                    out.append(
                        {
                            "axis": "bond",
                            "atoms": list(atom_set),
                            "atom_syms": [syms[atom_set[0]], syms[atom_set[1]]],
                            "center": a,
                            "obs": round(float(d), 3),
                            "median": round(float(med_b), 3),
                            "mad": round(float(max(mad_b, 0.01)), 4),
                            "sev_mad": round(float(sev), 2),
                            "level": level,
                            "n_cod": len(cod_bonds_pool),
                            "threshold": threshold,
                            "mode": mode,
                        }
                    )

        # --- angles ----------------------------------------------- #
        cod_angles_all = []
        for r in cod:
            cod_angles_all.extend(r.get("angles", []))
        if len(cod_angles_all) >= _MIN_MATCHES and angle_records:
            mad_a, med_a = _legacy_mad(cod_angles_all)
            for (b, c, ang) in angle_records:
                unusual, sev = _is_unusual_angle(
                    ang, med_a, mad_a, threshold, apply_cap=apply_cap
                )
                if not unusual:
                    continue
                atom_set = tuple(sorted((a, b, c)))
                key = ("angle", atom_set)
                if key in emitted:
                    continue
                emitted.add(key)
                out.append(
                    {
                        "axis": "angle",
                        "atoms": list(atom_set),
                        "atom_syms": [syms[i] for i in atom_set],
                        "center": a,
                        "obs": round(float(ang), 2),
                        "median": round(float(med_a), 2),
                        "mad": round(float(max(mad_a, 1.0)), 2),
                        "sev_mad": round(float(sev), 2),
                        "level": level,
                        "n_cod": len(cod_angles_all),
                        "threshold": threshold,
                        "mode": mode,
                    }
                )

        # --- torsions (CCDC mode only) ---------------------------- #
        if emit_torsion:
            # i-a-j is an angle; a torsion includes one more outer atom k
            # adjacent to one of the angle wings.
            cod_tor_all = []
            for r in cod:
                cod_tor_all.extend([abs(v) for v in r.get("torsions", []) if v is not None])
            if len(cod_tor_all) < _MIN_MATCHES:
                # Fall back to abs(improper) distribution which is also
                # collected per-record in the v3 index.
                for r in cod:
                    v = r.get("improper")
                    if v is not None:
                        cod_tor_all.append(abs(v))
            if len(cod_tor_all) >= _MIN_MATCHES:
                mad_t, med_t = _legacy_mad(cod_tor_all)
                mad_t = max(mad_t, TORSION_MAD_FLOOR_DEG)
                for (b, c, _ang) in angle_records:
                    # Walk one extra atom from b (skipping a and c).
                    for (d_atom, _dd) in adj.get(b, []):
                        if d_atom in (a, c, b):
                            continue
                        try:
                            tor = abs(
                                _legacy_signed_dihedral(P[d_atom], P[b], P[a], P[c])
                            )
                        except Exception:
                            continue
                        delta = abs(tor - med_t)
                        sev = delta / mad_t
                        if sev <= threshold:
                            continue
                        atom_set = tuple(sorted((d_atom, b, a, c)))
                        if len(set(atom_set)) < 4:
                            continue
                        key = ("torsion", atom_set)
                        if key in emitted:
                            continue
                        emitted.add(key)
                        out.append(
                            {
                                "axis": "torsion",
                                "atoms": list(atom_set),
                                "atom_syms": [syms[i] for i in atom_set],
                                "center": a,
                                "obs": round(float(tor), 2),
                                "median": round(float(med_t), 2),
                                "mad": round(float(mad_t), 2),
                                "sev_mad": round(float(sev), 2),
                                "level": level,
                                "n_cod": len(cod_tor_all),
                                "threshold": threshold,
                                "mode": mode,
                            }
                        )

    # Deterministic order
    out.sort(key=lambda r: (r["axis"], tuple(r["atoms"])))
    return out


# ----------------------------------------------------------------- #
#  CLI for quick manual checks
# ----------------------------------------------------------------- #
def _cli_main() -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("xyz")
    ap.add_argument("--mode", choices=("legacy", "ccdc"), default=None)
    ap.add_argument("--frame", type=int, default=0)
    args = ap.parse_args()
    text = open(args.xyz).read().splitlines()
    na = int(text[0].strip())
    start = args.frame * (na + 2)
    frame_lines = text[start : start + na + 2]
    if len(frame_lines) < na + 2:
        print("frame out of range", file=sys.stderr)
        return 2
    s, P, _ = _bd._parse("\n".join(frame_lines))
    recs = detect_anomalies_v3(s, P, mode=args.mode)
    print(f"# {args.xyz} frame={args.frame} mode={args.mode or 'env'} "
          f"anomalies={len(recs)}")
    for r in recs:
        print(json.dumps(r))
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli_main())
