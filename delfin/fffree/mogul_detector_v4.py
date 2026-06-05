"""Mogul-style per-feature anomaly detector v4 (CCDC + v5-GMM-aware).

v4 extends ``mogul_detector_v3`` with three additions, all driven by the
v5 mogul library (``grip_lib_v5.npz``):

1. **X-H bond awareness** — v3 excludes H from heavy-atom adjacency, so
   X-H bonds (M-H, N-H, O-H, C-H, etc.) are never emitted. v4 emits
   per-X-H records using the v5 ``pair_bond`` table for
   element-and-hybridisation-specific (mu, sigma) priors.
2. **M-D-aware angle channel** — v3 also excludes metals from the
   adjacency walk so D-M-D and M-D-X angles are never tested. v4 emits
   per-metal-angle records using v5 ``tm_triple_angle_block`` (preferred,
   metal-block aware) with fallback to ``triple_angle``.
3. **Torsion multimodal acceptance** — v3 uses a single MAD over a
   pool of torsion samples that may be tri-modal (gauche/anti/eclipsed
   for sp3, or syn/anti for hindered sp2). v4 looks the torsion key up
   in v5's GMM (1-/2-/3-component) and uses the *minimum* z-score across
   components weighted by their mixture weight. A "match" against any
   well-populated component (pi >= 0.05) drops the anomaly.

Default state
-------------
**Default-OFF**, byte-identical to v3 unless ``DELFIN_USE_MOGUL_V4=1``.
When the env-flag is set (or ``use_v4=True`` is passed explicitly), v4
extensions are applied. Otherwise v3 is invoked unchanged.

Library path
------------
``DELFIN_MOGUL_V5_LIB`` env (default
``/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5.npz``).

Public API
----------
``detect_anomalies_v4(syms, P, ...)``
    Same return type as ``detect_anomalies_v3`` plus additional records
    for X-H bonds and M-D angles. Includes a ``source`` key
    (``"v3"``/``"v4_xh"``/``"v4_md_angle"``) and, where the GMM was
    consulted, the chosen component.
"""
from __future__ import annotations

import json
import math
import os
import sys
import threading
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

sys.path.insert(0, "/home/qmchem_max/ComPlat/DELFIN")
import delfin._bond_decollapse as _bd  # noqa: E402

from . import mogul_detector_v3 as _v3  # noqa: E402

__all__ = [
    "detect_anomalies_v4",
    "MOGUL_V4_ENV",
    "MOGUL_V5_LIB_ENV",
    "DEFAULT_V5_LIB_PATH",
    "V5Lookup",
    "load_v5_lookup",
    "METALS",
]

MOGUL_V4_ENV = "DELFIN_USE_MOGUL_V4"
MOGUL_V5_LIB_ENV = "DELFIN_MOGUL_V5_LIB"
DEFAULT_V5_LIB_PATH = "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5.npz"

# CCDC equivalent: |z| >= 2 (Gaussian sigma is direct; no MAD inflation needed)
V4_GAUSSIAN_THRESHOLD = 2.0
# Minimum component weight to count as a real GMM mode
V4_GMM_PI_FLOOR = 0.05
# Min sample count to trust a v5 entry
V4_MIN_N = 15

METALS = set(
    [
        # alkali / alkaline earth
        "Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba",
        # d-block 3d
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        # d-block 4d
        "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        # d-block 5d
        "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        # p-block metals occasionally treated as TMC centres
        "Al", "Ga", "In", "Sn", "Tl", "Pb", "Bi",
        # f-block
        "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
        "Er", "Tm", "Yb", "Lu",
    ]
)


def _metal_block(sym: str) -> str:
    if sym in ("Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"):
        return "3d"
    if sym in ("Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd"):
        return "4d"
    if sym in ("La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"):
        return "5d"
    if sym in ("Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
               "Er", "Tm", "Yb", "Lu"):
        return "4f"
    return "*"


# --------------------------------------------------------------- #
#  v5 lookup wrapper (lazy + cached)
# --------------------------------------------------------------- #
class V5Lookup:
    """Cached lookup over the v5 lib's TM-aware tables."""

    def __init__(self, path: Optional[str] = None):
        self.path = path or os.environ.get(MOGUL_V5_LIB_ENV, DEFAULT_V5_LIB_PATH)
        self._loaded = False
        self._pair_bond: Dict[str, Tuple[float, float, int]] = {}
        self._triple_angle: Dict[str, Tuple[float, float, int]] = {}
        self._tm_pair: Dict[str, Tuple[float, float, int]] = {}
        self._tm_triple: Dict[str, Tuple[float, float, int]] = {}
        self._torsion_simple: Dict[str, int] = {}
        self._torsion_n_components = None
        self._torsion_pi = None
        self._torsion_mu = None
        self._torsion_sigma = None

    def _ensure_loaded(self) -> bool:
        if self._loaded:
            return True
        if not os.path.exists(self.path):
            self._loaded = True  # mark as loaded-empty to avoid re-attempts
            return False
        lib = np.load(self.path, allow_pickle=True)

        def _b(keys_k, mu_k, sigma_k, n_k, out):
            keys = lib[keys_k]
            mu = lib[mu_k]
            sigma = lib[sigma_k]
            n = lib[n_k]
            for i in range(len(keys)):
                out[str(keys[i])] = (float(mu[i]), float(sigma[i]), int(n[i]))

        _b("pair_bond_keys", "pair_bond_mu", "pair_bond_sigma", "pair_bond_n",
           self._pair_bond)
        _b("triple_angle_keys", "triple_angle_mu", "triple_angle_sigma",
           "triple_angle_n", self._triple_angle)
        try:
            _b("tm_pair_bond_block_keys", "tm_pair_bond_block_mu",
               "tm_pair_bond_block_sigma", "tm_pair_bond_block_n",
               self._tm_pair)
        except Exception:
            pass
        try:
            _b("tm_triple_angle_block_keys", "tm_triple_angle_block_mu",
               "tm_triple_angle_block_sigma", "tm_triple_angle_block_n",
               self._tm_triple)
        except Exception:
            pass
        # torsion GMM (sorted-element fallback)
        try:
            tk = lib["torsion_keys"]
            self._torsion_n_components = np.asarray(lib["torsion_n_components"])
            self._torsion_pi = np.asarray(lib["torsion_pi"])
            self._torsion_mu = np.asarray(lib["torsion_mu"])
            self._torsion_sigma = np.asarray(lib["torsion_sigma"])
            for i, k in enumerate(tk):
                try:
                    arr = json.loads(k)
                    elems = sorted([arr[0], arr[2], arr[4], arr[5]])
                    kk = "-".join(elems)
                    if kk not in self._torsion_simple:
                        self._torsion_simple[kk] = i
                except Exception:
                    continue
        except Exception:
            pass
        self._loaded = True
        return True

    # ---------- public lookups ---------- #
    def pair_bond(self, a_sym: str, a_hyb: str, b_sym: str, b_hyb: str):
        if not self._ensure_loaded():
            return None
        a, b = sorted([(a_sym, a_hyb), (b_sym, b_hyb)])
        for hyb_a in (a[1], "*"):
            for hyb_b in (b[1], "*"):
                k = json.dumps([a[0], hyb_a, b[0], hyb_b], separators=(",", ":"))
                if k in self._pair_bond:
                    return self._pair_bond[k]
        return None

    def triple_angle(self, center: str, center_hyb: str, w1: str, w2: str):
        if not self._ensure_loaded():
            return None
        wings = sorted([w1, w2])
        for hyb in (center_hyb, "*"):
            k = json.dumps([center, hyb, wings[0], wings[1]],
                           separators=(",", ":"))
            if k in self._triple_angle:
                return self._triple_angle[k]
        return None

    def tm_pair_bond(self, m: str, x_sym: str, x_hyb: str):
        if not self._ensure_loaded():
            return None
        block = _metal_block(m)
        a, b = sorted([(m, block), (x_sym, x_hyb)])
        for hyb_a in (a[1], "*"):
            for hyb_b in (b[1], "*"):
                k = json.dumps([a[0], hyb_a, b[0], hyb_b],
                               separators=(",", ":"))
                if k in self._tm_pair:
                    return self._tm_pair[k]
        return None

    def tm_triple_angle(self, m: str, w1: str, w2: str):
        if not self._ensure_loaded():
            return None
        block = _metal_block(m)
        wings = sorted([w1, w2])
        k = json.dumps([m, block, wings[0], wings[1]],
                       separators=(",", ":"))
        return self._tm_triple.get(k)

    def torsion_gmm(self, e0: str, e1: str, e2: str, e3: str):
        if not self._ensure_loaded() or self._torsion_n_components is None:
            return None
        kk = "-".join(sorted([e0, e1, e2, e3]))
        idx = self._torsion_simple.get(kk)
        if idx is None:
            return None
        return (
            int(self._torsion_n_components[idx]),
            self._torsion_pi[idx],
            self._torsion_mu[idx],
            self._torsion_sigma[idx],
        )


# Module-level shared lookup (lazy + cached + thread-safe init).
_GLOBAL_LOOKUP: Optional[V5Lookup] = None
_GLOBAL_LOCK = threading.Lock()


def load_v5_lookup(path: Optional[str] = None, force: bool = False) -> V5Lookup:
    """Return the module-level v5 lookup, building it if needed."""
    global _GLOBAL_LOOKUP
    with _GLOBAL_LOCK:
        if _GLOBAL_LOOKUP is None or force:
            _GLOBAL_LOOKUP = V5Lookup(path=path)
            # eager-load so first detection call doesn't pay the cost
            _GLOBAL_LOOKUP._ensure_loaded()
        return _GLOBAL_LOOKUP


# --------------------------------------------------------------- #
#  X-H + M-D detection helpers (used by detect_anomalies_v4)
# --------------------------------------------------------------- #
def _build_full_adjacency(syms, P) -> Dict[int, List[Tuple[int, float]]]:
    """Heavy-atom + H + metal adjacency. Symmetric, distance-based."""
    from delfin.fffree.polyhedra import COV

    def _r(s: str) -> float:
        return COV.get(s, _bd._COV.get(s, 0.9))

    n = len(syms)
    adj: Dict[int, List[Tuple[int, float]]] = {i: [] for i in range(n)}
    cutoff_factor = 1.30
    for a in range(n):
        for b in range(a + 1, n):
            d = float(np.linalg.norm(P[a] - P[b]))
            ra, rb = _r(syms[a]), _r(syms[b])
            # for metal-X use a slightly more generous cutoff
            f = cutoff_factor
            if _bd._is_metal(syms[a]) or _bd._is_metal(syms[b]):
                f = 1.35
            if d < f * (ra + rb):
                adj[a].append((b, d))
                adj[b].append((a, d))
    return adj


def _heavy_hyb_for_atom(syms, P, full_adj, a) -> str:
    """Rough sp/sp2/sp3 by ALL non-H neighbours (incl. metals).

    Mirrors v5 lib's hybridisation logic: total non-H connections matter
    for sp/sp2/sp3 classification (e.g. carbonyl C bound to Fe + O is sp,
    not sp3 — the metal-coordination DOES contribute one connection).
    """
    nh_nbrs = [b for (b, _) in full_adj.get(a, []) if syms[b] != "H"]
    n = len(nh_nbrs)
    if n == 0:
        return "*"
    if n >= 4:
        return "sp3"
    if n == 1:
        # Terminal heavy bond — use element-level heuristic
        s_a = syms[a]
        if s_a in ("C", "N"):
            return "sp"  # likely C#X / N#X / C=X terminus
        return "sp3"
    # n == 2 or 3 — use mean angle
    angles = []
    for i in range(len(nh_nbrs)):
        for j in range(i + 1, len(nh_nbrs)):
            v1 = P[nh_nbrs[i]] - P[a]
            v2 = P[nh_nbrs[j]] - P[a]
            n1 = float(np.linalg.norm(v1))
            n2 = float(np.linalg.norm(v2))
            if n1 < 1e-6 or n2 < 1e-6:
                continue
            cv = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
            angles.append(math.degrees(math.acos(cv)))
    if not angles:
        return "sp3"
    avg = sum(angles) / len(angles)
    if n == 2 and avg > 150:
        return "sp"
    if n == 3 and avg > 115:
        return "sp2"
    if n == 2 and avg > 115:
        return "sp2"
    return "sp3"


def _emit_xh_records(syms, P, full_adj, v5: V5Lookup, threshold: float
                     ) -> List[Dict]:
    """Test every X-H bond against v5 pair_bond statistics."""
    out: List[Dict] = []
    seen: set = set()
    for h in range(len(syms)):
        if syms[h] != "H":
            continue
        # find the unique heavy neighbour (closest)
        best = (1e9, -1)
        for (nbr, d) in full_adj.get(h, []):
            if syms[nbr] == "H":
                continue
            if d < best[0]:
                best = (d, nbr)
        if best[1] < 0:
            continue
        x = best[1]
        d_xh = best[0]
        x_sym = syms[x]
        x_hyb = _heavy_hyb_for_atom(syms, P, full_adj, x)
        # try TM-specific table first for M-H bonds
        rec = None
        if x_sym in METALS:
            rec = v5.tm_pair_bond(x_sym, "H", "*")
        if rec is None:
            rec = v5.pair_bond("H", "*", x_sym, x_hyb)
        if rec is None:
            continue
        mu, sigma, n = rec
        if n < V4_MIN_N or sigma < 1e-6:
            continue
        sev = abs(d_xh - mu) / sigma
        if sev < threshold:
            continue
        atoms = tuple(sorted((h, x)))
        if ("bond", atoms) in seen:
            continue
        seen.add(("bond", atoms))
        out.append({
            "axis": "bond",
            "atoms": list(atoms),
            "atom_syms": [syms[atoms[0]], syms[atoms[1]]],
            "center": x,
            "obs": round(float(d_xh), 3),
            "median": round(float(mu), 3),
            "mad": round(float(sigma), 4),  # report sigma in mad slot
            "sev_mad": round(float(sev), 2),
            "level": "v5_pair_bond",
            "n_cod": int(n),
            "threshold": float(threshold),
            "mode": "v4",
            "source": "v4_xh",
        })
    out.sort(key=lambda r: (r["axis"], tuple(r["atoms"])))
    return out


def _emit_md_angle_records(syms, P, full_adj, v5: V5Lookup, threshold: float
                           ) -> List[Dict]:
    """Test D-M-D and M-D-X angles against v5 tm_triple_angle / triple_angle.

    D-M-D: center is metal, both wings are heavy-non-metal donors.
    M-D-X: center is donor (non-metal), one wing is metal, other wing
           is a heavy non-metal neighbour.
    """
    out: List[Dict] = []
    seen: set = set()
    n = len(syms)
    for c in range(n):
        c_sym = syms[c]
        if c_sym == "H":
            continue
        heavy_nbrs = [b for (b, _) in full_adj.get(c, []) if syms[b] != "H"]
        if len(heavy_nbrs) < 2:
            continue
        c_is_metal = c_sym in METALS
        for i in range(len(heavy_nbrs)):
            for j in range(i + 1, len(heavy_nbrs)):
                w1, w2 = heavy_nbrs[i], heavy_nbrs[j]
                w1_sym, w2_sym = syms[w1], syms[w2]
                w1_metal = w1_sym in METALS
                w2_metal = w2_sym in METALS

                # classify
                if c_is_metal:
                    if w1_metal or w2_metal:
                        # skip M-M-M for now (not in v5 priors)
                        continue
                    source = "v4_dmd"
                elif (w1_metal and not w2_metal) or (w2_metal and not w1_metal):
                    source = "v4_mdx"
                else:
                    continue

                v_a = P[w1] - P[c]
                v_b = P[w2] - P[c]
                na = float(np.linalg.norm(v_a))
                nb = float(np.linalg.norm(v_b))
                if na < 1e-6 or nb < 1e-6:
                    continue
                cv = float(np.clip(np.dot(v_a, v_b) / (na * nb), -1.0, 1.0))
                ang = math.degrees(math.acos(cv))

                rec = None
                if source == "v4_dmd":
                    rec = v5.tm_triple_angle(c_sym, w1_sym, w2_sym)
                if rec is None:
                    c_hyb = _heavy_hyb_for_atom(syms, P, full_adj, c) if not c_is_metal else "*"
                    rec = v5.triple_angle(c_sym, c_hyb, w1_sym, w2_sym)
                if rec is None:
                    continue
                mu, sigma, n_obs = rec
                if n_obs < V4_MIN_N or sigma < 1e-6:
                    continue
                sev = abs(ang - mu) / sigma
                if sev < threshold:
                    continue
                atoms = tuple(sorted((c, w1, w2)))
                # D-M-D and M-D-X have the same atom set but different
                # vertex, so the dedup key must include the center.
                key = ("angle", atoms, c, source)
                if key in seen:
                    continue
                seen.add(key)
                out.append({
                    "axis": "angle",
                    "atoms": list(atoms),
                    "atom_syms": [syms[i] for i in atoms],
                    "center": c,
                    "obs": round(float(ang), 2),
                    "median": round(float(mu), 2),
                    "mad": round(float(sigma), 3),
                    "sev_mad": round(float(sev), 2),
                    "level": "v5_tm_angle" if source == "v4_dmd" else "v5_triple_angle",
                    "n_cod": int(n_obs),
                    "threshold": float(threshold),
                    "mode": "v4",
                    "source": source,
                })
    out.sort(key=lambda r: (r["axis"], tuple(r["atoms"])))
    return out


def _reclassify_torsions_with_gmm(records: List[Dict], syms, v5: V5Lookup,
                                  threshold: float) -> Tuple[List[Dict], int]:
    """Pass-through v3 torsion records but DROP those that match a v5 GMM
    component within ``threshold``.

    Returns (kept_records, n_dropped).
    """
    kept: List[Dict] = []
    dropped = 0
    for rec in records:
        if rec.get("axis") != "torsion":
            kept.append(rec)
            continue
        atom_syms = rec.get("atom_syms") or []
        if len(atom_syms) != 4:
            kept.append(rec)
            continue
        gmm = v5.torsion_gmm(*atom_syms)
        if gmm is None:
            kept.append(rec)
            continue
        n_comp, pi, mu, sigma = gmm
        obs = float(rec.get("obs", 0.0))
        # v3 emits abs(tor) — try both signs since true sign is unknown
        match = False
        best_sev = 1e9
        best_comp = -1
        for sign in (1.0, -1.0):
            ob = sign * obs
            for c in range(n_comp):
                if not np.isfinite(mu[c]) or not np.isfinite(sigma[c]):
                    continue
                if sigma[c] < 1e-9:
                    continue
                if pi[c] < V4_GMM_PI_FLOOR:
                    continue
                s = abs(ob - mu[c]) / sigma[c]
                if s < best_sev:
                    best_sev = s
                    best_comp = c
                if s < threshold:
                    match = True
                    break
            if match:
                break
        if match:
            dropped += 1
            continue
        # update rec with best component context, but still keep it
        rec2 = dict(rec)
        rec2["v5_gmm_best_sev"] = round(float(best_sev), 2)
        rec2["v5_gmm_best_comp"] = int(best_comp)
        rec2["v5_gmm_n_components"] = int(n_comp)
        rec2["source"] = "v3"
        kept.append(rec2)
    return kept, dropped


# --------------------------------------------------------------- #
#  Public detector
# --------------------------------------------------------------- #
def detect_anomalies_v4(
    syms: Sequence[str],
    P,
    *,
    use_v4: Optional[bool] = None,
    index: Optional[Dict] = None,
    v5_lib_path: Optional[str] = None,
) -> List[Dict]:
    """Mogul v4 detector — v3 + (X-H bonds, M-D angles, GMM-torsions).

    Parameters
    ----------
    syms, P
        Atom symbols and Nx3 coordinates (Å).
    use_v4
        If ``True``, apply v4 extensions. If ``False``, return raw v3
        records. If ``None`` (default), consult ``DELFIN_USE_MOGUL_V4``
        env: '1' enables v4; anything else falls back to v3.
    index
        Pre-loaded COD fragment index for the v3 pass.
    v5_lib_path
        Override v5 library path (test only).

    Returns
    -------
    list of dict
        Combined v3 + v4 records, sorted by ``(axis, atoms)``.
    """
    if use_v4 is None:
        use_v4 = os.environ.get(MOGUL_V4_ENV, "") == "1"
    if not use_v4:
        # byte-identical pass-through: respect existing v3 env semantics
        return _v3.detect_anomalies_v3(syms, P, index=index)

    # v4 mode — start from v3 CCDC anomalies
    v3_recs = _v3.detect_anomalies_v3(syms, P, mode="ccdc", index=index)
    # tag v3 records with source
    for r in v3_recs:
        r.setdefault("source", "v3")

    P = np.asarray(P, dtype=float)
    if len(syms) < 3:
        return sorted(v3_recs, key=lambda r: (r["axis"], tuple(r["atoms"])))

    v5 = load_v5_lookup(path=v5_lib_path)
    if not v5._ensure_loaded():
        # v5 not on disk — return v3 only, but still tagged
        return sorted(v3_recs, key=lambda r: (r["axis"], tuple(r["atoms"])))

    threshold = V4_GAUSSIAN_THRESHOLD

    full_adj = _build_full_adjacency(syms, P)

    # 1. X-H bonds (new channel)
    xh_recs = _emit_xh_records(syms, P, full_adj, v5, threshold)

    # 2. M-D angles (new channel)
    md_recs = _emit_md_angle_records(syms, P, full_adj, v5, threshold)

    # 3. Torsions: reclassify v3 torsions using v5 GMM
    reclassified, _dropped = _reclassify_torsions_with_gmm(
        v3_recs, syms, v5, threshold
    )

    combined = reclassified + xh_recs + md_recs
    # de-duplicate by (axis, sorted-atom-tuple, center) — different
    # centres represent different physical angles (D-M-D vs M-D-X).
    # Within a single (axis, atoms, center) we prefer v4 records.
    by_key: Dict[Tuple, Dict] = {}
    for r in combined:
        key = (r["axis"], tuple(sorted(r["atoms"])), r.get("center"))
        existing = by_key.get(key)
        if existing is None:
            by_key[key] = r
            continue
        # tie-break: v4 source wins over v3
        if existing.get("source", "") == "v3" and r.get("source", "v3") != "v3":
            by_key[key] = r
    out = list(by_key.values())
    out.sort(key=lambda r: (r["axis"], tuple(r["atoms"]),
                             r.get("center") or -1))
    return out
