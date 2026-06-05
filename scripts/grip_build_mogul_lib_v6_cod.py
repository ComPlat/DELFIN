#!/usr/bin/env python3
"""GRIP v6: COD fragment library (open-access twin of v5 CCDC).

User-Direktive 2026-06-05: "beides nutzten" — both COD and CCDC.

v6 = open-access analogue of v5, built from the Crystallography Open Database
(COD) instead of CCDC.  Preserves the v5 schema byte-for-byte so the merged
lookup-chain in :mod:`delfin.fffree.grip_mogul_lookup` can splice COD and
CCDC libraries together without branch divergence.

WHY COD (in addition to CCDC)
-----------------------------
1. **Open-access provenance** — paper claim "we use COD" is reproducible by
   reviewers without a CCDC licence.
2. **CCDC legal-separation doctrine** (feedback_ccdc_legal_separation_doctrine_2026_06_04.md):
   CCDC for us/build, NOT for adopters/runtime.  v6 is the open answer.
3. **Cross-validation** — different curators, different methodology.
   Fragments where COD and CCDC disagree = data-quality flags.
4. **Recent open entries** that CCDC has not yet indexed.

DATA SOURCE
-----------
``/home/qmchem_max/agent_workspace/quality_framework/COD_all/clean/`` =
51,196 cleaned COD entries already exported as XYZ (Cartesian) with a
header comment carrying COD-id / metal / formula / hapticity.  The dataset
is pre-curated:
  * largest metal-bearing component already extracted
  * disorder already collapsed (no occupancy<1 atoms)
  * counter-ions / solvents already dropped
This matches v5's Lever 1 cleaning intent natively, so we do NOT need to
re-implement the CCDC mol.components / secondary-disorder logic.

v6 LEVERS (mirror v5 perfectly where applicable)
-------------------------------------------------
Lever 1 — cleaning : N/A (clean COD is already curated)
Lever 2 — TM categories : SAME as v5
    (carbene, hapto_eta2/5/6, mu_bridge, agostic, ox_addition)
Lever 3 — Mogul native validation : N/A (no Mogul on COD; mogul_validated=False)
Lever 4 — failed-entry recovery : SAME structural intent
    (skip bad-coord atoms, keep good ones)
Lever 5 — TM disaggregation 3d/4d/5d/f : SAME as v5

DETERMINISM
-----------
* PYTHONHASHSEED=0 set at import.
* COD files processed in lex-sorted filename order.
* Same robust aggregation (median + 1.4826·MAD).
* Same GMM torsion fits (random_state=42, n_init=3, 100k subsample cap).

USAGE
-----
::

    PYTHONHASHSEED=0 \\
    /home/qmchem_max/micromamba/envs/delfin/bin/python \\
        scripts/grip_build_mogul_lib_v6_cod.py \\
        --workers 32 \\
        --out reports/grip_lib_v6_cod.npz

OUTPUT SCHEMA
-------------
Schema-identical to v5 (so the v5 GripLibrary class loads it unchanged),
plus a few v6-specific provenance fields::

    version                     int32 = 6
    source                      str   = "COD"
    cod_n_entries_scanned       int32
    cod_n_extracted             int32
    cleaning_applied            int32 = 1   (clean COD is pre-cleaned)
    (... all v5 tables identical ...)

The lookup module (:mod:`grip_mogul_lookup`) gains an env flag
``DELFIN_FFFREE_GRIP_LIB_COD_PATH``.  When set, the lib is loaded ALONGSIDE
the CCDC v5 lib and the lookup chain becomes:
    merged → ccdc-only → cod-only → fallback
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import multiprocessing as mp
import os
import random
import sys
import time
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# Determinism BEFORE other imports
os.environ.setdefault("PYTHONHASHSEED", "0")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
np.random.seed(42)
random.seed(42)
warnings.filterwarnings("ignore")

ROOT = Path("/home/qmchem_max/agent_workspace/quality_framework")
COD_CLEAN_DIR = ROOT / "COD_all" / "clean"
DEFAULT_OUT = ROOT / "reports" / "grip_lib_v6_cod.npz"
DEFAULT_LOG = Path("/tmp/grip_lib_v6_cod_build.log")


def _setup_logging(log_path: Path) -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[
            logging.FileHandler(log_path, mode="w"),
            logging.StreamHandler(sys.stdout),
        ],
        force=True,
    )
    return logging.getLogger("grip_build_v6")


log = logging.getLogger("grip_build_v6")


# ============================================================================
# Conventions — verbatim from v5 to keep keys byte-identical
# ============================================================================

_METALS_D = frozenset({
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
})
_METALS_F = frozenset({
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
    "Ac", "Th", "Pa", "U", "Np", "Pu",
})
_METALS = _METALS_D | _METALS_F

_METAL_BLOCK = {}
for s in ("Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"):
    _METAL_BLOCK[s] = "3d"
for s in ("Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd"):
    _METAL_BLOCK[s] = "4d"
for s in ("Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"):
    _METAL_BLOCK[s] = "5d"
for s in _METALS_F:
    _METAL_BLOCK[s] = "f"


def _block_of(sym: str) -> str:
    return _METAL_BLOCK.get(sym, "*")


_D8_SP_CANDIDATES = frozenset({"Ni", "Pd", "Pt", "Rh", "Ir"})


def _infer_hyb(symbol: str, heavy_degree: int, total_degree: int) -> str:
    """Identical to v5 ``_infer_hyb`` (metals -> '*')."""
    if symbol in _METALS:
        return "*"
    deg = int(total_degree)
    if symbol == "C":
        if deg <= 2:
            return "sp"
        if deg == 3:
            return "sp2"
        return "sp3"
    if symbol == "N":
        if deg <= 1:
            return "sp"
        if deg == 2:
            return "sp2"
        if deg == 3:
            return "sp3"
        return "sp3"
    if symbol == "O":
        if deg <= 1:
            return "sp2"
        if deg == 2:
            return "sp3"
        return "sp3"
    if symbol == "S":
        if deg <= 1:
            return "sp2"
        if deg == 2:
            return "sp3"
        if deg <= 4:
            return "sp3"
        return "sp3"
    if symbol == "P":
        if deg <= 3:
            return "sp3"
        if deg == 4:
            return "sp3"
        if deg == 5:
            return "sp3d"
        return "sp3d2"
    if symbol in ("F", "Cl", "Br", "I"):
        return "sp3"
    if symbol == "H":
        return "*"
    if symbol == "B":
        if deg == 3:
            return "sp2"
        return "sp3"
    if symbol == "Si":
        return "sp3"
    return "sp3"


def _to_key_str(parsed) -> str:
    return json.dumps(parsed, separators=(",", ":"), ensure_ascii=False)


def _fallback_levels(parsed) -> list[str]:
    """Mirror grip_mogul_lookup._fallback_levels."""
    if not isinstance(parsed, list) or len(parsed) < 4:
        return [_to_key_str(parsed)]
    Zc, hyb_c, neighbors, second_shell = parsed[0], parsed[1], parsed[2], parsed[3]
    levels = [
        _to_key_str([Zc, hyb_c, neighbors, second_shell]),
        _to_key_str([Zc, hyb_c, neighbors, []]),
    ]
    neigh_no_ring = [[n[0], n[1], -1, n[3]] if len(n) >= 4 else n for n in neighbors]
    levels.append(_to_key_str([Zc, hyb_c, neigh_no_ring, []]))
    neigh_no_hyb = [[n[0], n[1], -1, "*"] if len(n) >= 4 else n for n in neighbors]
    levels.append(_to_key_str([Zc, hyb_c, neigh_no_hyb, []]))
    neigh_no_Z = [["*", n[1], -1, "*"] if len(n) >= 4 else n for n in neighbors]
    levels.append(_to_key_str([Zc, hyb_c, neigh_no_Z, []]))
    levels.append(_to_key_str([Zc, hyb_c, [["*", -1, -1, "*"]] * len(neighbors), []]))
    seen: set[str] = set()
    out: list[str] = []
    for s in levels:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


def _torsion_canonicalise(Z_a, Z_b, hyb_b, Z_c, hyb_c, Z_d,
                          ring_bc: int, arom_b: bool, arom_c: bool) -> list:
    bc1 = (str(Z_b), str(hyb_b))
    bc2 = (str(Z_c), str(hyb_c))
    if bc1 <= bc2:
        return [str(Z_b), str(hyb_b), str(Z_c), str(hyb_c),
                str(Z_a), str(Z_d), int(ring_bc),
                bool(arom_b), bool(arom_c)]
    return [str(Z_c), str(hyb_c), str(Z_b), str(hyb_b),
            str(Z_d), str(Z_a), int(ring_bc),
            bool(arom_c), bool(arom_b)]


def _to_torsion_key_str(parsed: list) -> str:
    return json.dumps(parsed, separators=(",", ":"), ensure_ascii=False)


def _torsion_fallback_levels(parsed: list) -> list[str]:
    if not isinstance(parsed, list) or len(parsed) != 9:
        return [_to_torsion_key_str(parsed)]
    Zb, hyb_b, Zc, hyb_c, Za, Zd, ring_bc, arom_b, arom_c = parsed
    levels: list[str] = []
    levels.append(_to_torsion_key_str([Zb, hyb_b, Zc, hyb_c, Za, Zd, ring_bc, arom_b, arom_c]))
    levels.append(_to_torsion_key_str([Zb, hyb_b, Zc, hyb_c, Za, Zd, -1, arom_b, arom_c]))
    levels.append(_to_torsion_key_str([Zb, hyb_b, Zc, hyb_c, Za, Zd, -1, False, False]))
    levels.append(_to_torsion_key_str([Zb, hyb_b, Zc, hyb_c, "*", "*", -1, False, False]))
    levels.append(_to_torsion_key_str([Zb, "*", Zc, "*", "*", "*", -1, False, False]))
    seen: set[str] = set()
    out: list[str] = []
    for s in levels:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


def _pair_bond_key(z1: str, hyb1: str, z2: str, hyb2: str) -> str:
    a = (str(z1), str(hyb1))
    b = (str(z2), str(hyb2))
    lo, hi = (a, b) if a <= b else (b, a)
    return json.dumps([lo[0], lo[1], hi[0], hi[1]],
                      separators=(",", ":"), ensure_ascii=False)


def _triple_angle_key(zc: str, hyb_c: str, z_l: str, z_r: str) -> str:
    lo, hi = (z_l, z_r) if z_l <= z_r else (z_r, z_l)
    return json.dumps([str(zc), str(hyb_c), str(lo), str(hi)],
                      separators=(",", ":"), ensure_ascii=False)


def _improper_pair_key(zc: str, hyb_c: str, nbr_zs: List[str]) -> str:
    nzs = sorted(str(z) for z in nbr_zs)
    return json.dumps([str(zc), str(hyb_c), nzs],
                      separators=(",", ":"), ensure_ascii=False)


# ============================================================================
# Lever 2 — TM fragment category detector tolerances (identical to v5)
# ============================================================================
_AGOSTIC_MH_LO, _AGOSTIC_MH_HI = 1.7, 2.6
_AGOSTIC_MC_LO, _AGOSTIC_MC_HI = 2.0, 3.1
_HAPTO_MC_LO, _HAPTO_MC_HI = 1.9, 2.8


# ============================================================================
# Worker — read clean-COD XYZ, recover bonds with RDKit, emit v5-schema keys
# ============================================================================

def _worker_init():
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    np.random.seed(42)
    random.seed(42)


def _read_cod_xyz(path: str) -> Optional[Tuple[List[str], np.ndarray]]:
    """Parse a clean-COD XYZ file.

    Header layout (line 1 = n, line 2 = ``COD:<id> metal:<m> formula:<f> hapto:<h>``,
    lines 3..n+2 = ``<sym> x y z``).  Returns ``(symbols, positions)`` or
    ``None`` on parse failure.
    """
    try:
        with open(path) as fh:
            n_line = fh.readline().strip()
            try:
                n = int(n_line)
            except ValueError:
                return None
            _ = fh.readline()  # metadata
            syms: List[str] = []
            coords: List[Tuple[float, float, float]] = []
            for _ in range(n):
                line = fh.readline()
                if not line:
                    break
                parts = line.split()
                if len(parts) < 4:
                    continue
                sym = parts[0]
                try:
                    x = float(parts[1]); y = float(parts[2]); z = float(parts[3])
                except ValueError:
                    continue
                if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)):
                    continue
                syms.append(sym)
                coords.append((x, y, z))
        if len(syms) < 4 or len(syms) > 5000:
            return None
        return syms, np.asarray(coords, dtype=np.float64)
    except Exception:
        return None


# Covalent radii table (Angstrom) for fast distance-based bond detection.
# Source: revised values used by ASE / OpenBabel.  We bond when
# d <= (r_i + r_j) * 1.20 and d >= 0.4 (i.e. avoid duplicate sites).
_COV_R = {
    "H": 0.31, "He": 0.28,
    "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66,
    "F": 0.57, "Ne": 0.58,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05,
    "Cl": 1.02, "Ar": 1.06,
    "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39,
    "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
    "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20, "Kr": 1.16,
    "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54,
    "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
    "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.40,
    "Cs": 2.44, "Ba": 2.15,
    "La": 2.07, "Ce": 2.04, "Pr": 2.03, "Nd": 2.01, "Pm": 1.99, "Sm": 1.98,
    "Eu": 1.98, "Gd": 1.96, "Tb": 1.94, "Dy": 1.92, "Ho": 1.92, "Er": 1.89,
    "Tm": 1.90, "Yb": 1.87, "Lu": 1.87,
    "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41,
    "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "Tl": 1.45, "Pb": 1.46, "Bi": 1.48,
    "Ac": 2.15, "Th": 2.06, "Pa": 2.00, "U": 1.96, "Np": 1.90, "Pu": 1.87,
}
_BOND_SCALE = 1.20


def _detect_bonds(symbols: List[str], pos: np.ndarray) -> List[List[Tuple[int, np.ndarray, float]]]:
    """Distance-based bond detection using scaled covalent radii.

    Returns ``nbrs`` such that ``nbrs[i]`` is a list of
    ``(j, displacement_vector_i_to_j, distance)``.  Filters d>=0.4 and
    d<=r_sum*BOND_SCALE (typical 1.6-3.2 A for organics, up to 3.5 A for M-X).
    Also hard-caps d<=4.5 A as in v5 to avoid spurious neighbours.
    """
    n = len(symbols)
    nbrs: List[List[Tuple[int, np.ndarray, float]]] = [[] for _ in range(n)]
    # Build covrad array (default 0.77 if unknown)
    rcov = np.array([_COV_R.get(s, 0.77) for s in symbols], dtype=np.float64)
    # Use scipy KD-tree for O(n log n)
    try:
        from scipy.spatial import cKDTree
        tree = cKDTree(pos)
        # max distance any element pair ~ 2*max(rcov)*BOND_SCALE  ~ 6 A safely
        max_d = float(2 * rcov.max() * _BOND_SCALE)
        pairs = tree.query_pairs(r=max_d, output_type="ndarray")
    except Exception:
        # Fallback: brute force (only for small molecules)
        pairs = []
        for i in range(n):
            for j in range(i + 1, n):
                d = float(np.linalg.norm(pos[i] - pos[j]))
                if d <= 4.5:
                    pairs.append((i, j))
        pairs = np.asarray(pairs, dtype=np.int64) if pairs else np.zeros((0, 2), dtype=np.int64)
    for ij in pairs:
        i, j = int(ij[0]), int(ij[1])
        disp = pos[j] - pos[i]
        d = float(np.linalg.norm(disp))
        if d < 0.4 or d > 4.5:
            continue
        cutoff = (rcov[i] + rcov[j]) * _BOND_SCALE
        if d > cutoff:
            continue
        nbrs[i].append((j, disp, d))
        nbrs[j].append((i, -disp, d))
    return nbrs


def _extract_one_entry(path: str) -> Optional[dict]:
    """COD analogue of v5 ``_extract_one_entry``.

    Differences from CCDC version:
        * no entry-level disorder / polymeric / no_3d flags (pre-cleaned)
        * no mol.components (already main component only)
        * neighbours derived from covalent-radius distance test, not ccdc.bonds
    All emitted keys / fragment categories are byte-identical to v5.
    """
    audit = {
        "had_disorder": False,
        "had_partial_coord": False,
        "n_components_total": 1,
        "n_components_metal": 1,
        "was_cleaned": True,  # source already cleaned
        "fail_reason": "",
    }
    parsed_xyz = _read_cod_xyz(path)
    if parsed_xyz is None:
        audit["fail_reason"] = "read_fail"
        return {"audit": audit}
    symbols, positions_np = parsed_xyz
    n = len(symbols)

    nbrs = _detect_bonds(symbols, positions_np)

    heavy_deg = [
        sum(1 for (j, _, _) in nbrs[i] if symbols[j] != "H")
        for i in range(n)
    ]
    total_deg = [len(nbrs[i]) for i in range(n)]
    hyb = [_infer_hyb(symbols[i], heavy_deg[i], total_deg[i]) for i in range(n)]

    nonmetal_idx = [i for i in range(n)
                    if symbols[i] not in _METALS and symbols[i] != "H"]
    if len(nonmetal_idx) < 3:
        audit["fail_reason"] = "too_few_nonmetals"
        return {"audit": audit}
    nm_ok = sum(1 for i in nonmetal_idx if total_deg[i] <= 6)
    if nm_ok < int(0.8 * len(nonmetal_idx)):
        audit["fail_reason"] = "high_degree_filter"
        return {"audit": audit}

    # Heavy adjacency for ring detection (verbatim v5)
    adj_heavy: List[List[int]] = [[] for _ in range(n)]
    for i in range(n):
        if symbols[i] == "H":
            continue
        for (j, _, _) in nbrs[i]:
            if symbols[j] == "H":
                continue
            adj_heavy[i].append(j)
        adj_heavy[i] = sorted(set(adj_heavy[i]))

    MAX_RING = 8

    def _shortest_path_avoiding_edge(src: int, dst: int,
                                     forbid_uv: Tuple[int, int]) -> int:
        if src == dst:
            return 0
        dist = {src: 0}
        frontier = [src]
        for _depth in range(MAX_RING):
            if not frontier:
                return -1
            next_frontier = []
            for u in frontier:
                for v in adj_heavy[u]:
                    if (u == forbid_uv[0] and v == forbid_uv[1]) or \
                       (u == forbid_uv[1] and v == forbid_uv[0]):
                        continue
                    if v in dist:
                        continue
                    dist[v] = dist[u] + 1
                    if v == dst:
                        return dist[v]
                    if dist[v] < MAX_RING:
                        next_frontier.append(v)
            frontier = next_frontier
        return -1

    ring_min_atom = [-1] * n
    bond_ring_min: Dict[Tuple[int, int], int] = {}
    for i in range(n):
        if symbols[i] == "H":
            continue
        for j in adj_heavy[i]:
            if j <= i:
                continue
            sp = _shortest_path_avoiding_edge(i, j, (i, j))
            if sp > 0:
                ring_len = sp + 1
                if 3 <= ring_len <= MAX_RING:
                    bond_ring_min[(i, j)] = ring_len
    for (a, b), rl in bond_ring_min.items():
        if ring_min_atom[a] == -1 or rl < ring_min_atom[a]:
            ring_min_atom[a] = rl
        if ring_min_atom[b] == -1 or rl < ring_min_atom[b]:
            ring_min_atom[b] = rl

    def _bond_rmin(a: int, b: int) -> int:
        key = (min(a, b), max(a, b))
        return bond_ring_min.get(key, -1)

    rings_of_size: Dict[int, List[List[int]]] = {5: [], 6: []}
    if n <= 2000:
        seen_rings: set = set()
        for start in range(n):
            if symbols[start] == "H":
                continue
            stack = [(start, [start])]
            while stack:
                u, path = stack.pop()
                if len(path) > 6:
                    continue
                for v in adj_heavy[u]:
                    if symbols[v] == "H":
                        continue
                    if len(path) >= 3 and v == path[0] and len(path) in (5, 6):
                        canon = tuple(sorted(path))
                        if canon not in seen_rings:
                            seen_rings.add(canon)
                            rings_of_size[len(path)].append(list(path))
                    elif v not in path:
                        stack.append((v, path + [v]))

    is_aromatic = [False] * n
    for i in range(n):
        if symbols[i] == "H":
            continue
        if ring_min_atom[i] not in (5, 6):
            continue
        if hyb[i] != "sp2":
            continue
        nb_sp2 = sum(1 for j in adj_heavy[i] if hyb[j] == "sp2")
        if nb_sp2 >= 2:
            is_aromatic[i] = True

    # ---- v3 + v4 + v5 fragment emission ----
    bonds_out: Dict[str, List[float]] = defaultdict(list)
    angles_out: Dict[str, List[float]] = defaultdict(list)
    impropers_out: Dict[str, List[float]] = defaultdict(list)
    torsions_out: Dict[str, List[float]] = defaultdict(list)
    pair_bonds_out: Dict[str, List[float]] = defaultdict(list)
    triple_angles_out: Dict[str, List[float]] = defaultdict(list)
    improper_pairs_out: Dict[str, List[float]] = defaultdict(list)
    tm_carbene_out: Dict[str, List[float]] = defaultdict(list)
    tm_eta2_out: Dict[str, List[float]] = defaultdict(list)
    tm_eta5_out: Dict[str, List[float]] = defaultdict(list)
    tm_eta6_out: Dict[str, List[float]] = defaultdict(list)
    tm_mu_bridge_out: Dict[str, List[float]] = defaultdict(list)
    tm_agostic_out: Dict[str, List[float]] = defaultdict(list)
    tm_ox_addition_out: Dict[str, List[float]] = defaultdict(list)
    tm_pair_bond_block_out: Dict[str, List[float]] = defaultdict(list)
    tm_triple_angle_block_out: Dict[str, List[float]] = defaultdict(list)

    def _is_metal(sym: str) -> bool:
        return sym in _METALS

    centre_keys: List[Optional[str]] = [None] * n
    for i in range(n):
        sym_i = symbols[i]
        if sym_i == "H" or _is_metal(sym_i):
            continue
        hyb_i = hyb[i]
        if hyb_i == "*":
            continue
        first_shell = []
        for (j, _, _) in nbrs[i]:
            sym_j = symbols[j]
            hyb_j = hyb[j] if sym_j != "H" else "*"
            ring_j = ring_min_atom[j]
            first_shell.append([sym_j, 1, int(ring_j), str(hyb_j)])
        first_shell.sort(key=lambda v: (str(v[0]), str(v[3]), int(v[1]), int(v[2])))
        second_shell = []
        for (j, _, _) in nbrs[i]:
            if symbols[j] == "H":
                continue
            for (k, _, _) in nbrs[j]:
                if k == i:
                    continue
                sym_k = symbols[k]
                if sym_k == "H" or _is_metal(sym_k):
                    continue
                second_shell.append([sym_k, total_deg[k], hyb[k]])
        second_shell.sort(key=lambda v: (str(v[0]), int(v[1]), str(v[2])))
        centre_keys[i] = _to_key_str([sym_i, hyb_i, first_shell, second_shell])

    # ---- v3 bonds ----
    for i in range(n):
        key = centre_keys[i]
        if key is None:
            continue
        for (j, _, d) in nbrs[i]:
            if 0.5 < d < 3.5:
                bonds_out[key].append(float(d))

    # ---- v4 pair_bonds + v5 block-disagg ----
    for i in range(n):
        sym_i = symbols[i]
        if sym_i == "H":
            continue
        hyb_i = hyb[i]
        for (j, _, d) in nbrs[i]:
            if j <= i:
                continue
            if not (0.5 < d < 3.5):
                continue
            sym_j = symbols[j]
            hyb_j = hyb[j]
            k = _pair_bond_key(sym_i, hyb_i, sym_j, hyb_j)
            pair_bonds_out[k].append(float(d))
            if _is_metal(sym_i) or _is_metal(sym_j):
                hyb_i_b = _block_of(sym_i) if _is_metal(sym_i) else hyb_i
                hyb_j_b = _block_of(sym_j) if _is_metal(sym_j) else hyb_j
                kb = _pair_bond_key(sym_i, hyb_i_b, sym_j, hyb_j_b)
                tm_pair_bond_block_out[kb].append(float(d))

    # ---- v3 angles + v4 triple_angles + v5 block-disagg ----
    for b in range(n):
        key = centre_keys[b]
        nbr_local = nbrs[b]
        if len(nbr_local) < 2:
            continue
        for ia in range(len(nbr_local)):
            for ic in range(ia + 1, len(nbr_local)):
                ja, vec_a, _da = nbr_local[ia]
                jc, vec_c, _dc = nbr_local[ic]
                na = float(np.linalg.norm(vec_a))
                nc = float(np.linalg.norm(vec_c))
                if na < 1e-6 or nc < 1e-6:
                    continue
                cos_a = float(np.dot(vec_a, vec_c) / (na * nc))
                cos_a = max(-1.0, min(1.0, cos_a))
                angle_deg = math.degrees(math.acos(cos_a))
                if not (15.0 <= angle_deg <= 175.0):
                    continue
                if key is not None:
                    angles_out[key].append(angle_deg)
                sym_b = symbols[b]
                if sym_b == "H":
                    continue
                k_t = _triple_angle_key(sym_b, hyb[b], symbols[ja], symbols[jc])
                triple_angles_out[k_t].append(angle_deg)
                if _is_metal(sym_b):
                    k_tb = _triple_angle_key(
                        sym_b, _block_of(sym_b),
                        symbols[ja], symbols[jc]
                    )
                    tm_triple_angle_block_out[k_tb].append(angle_deg)

    # ---- v3 impropers + v4 improper_pairs ----
    for c_idx in range(n):
        if symbols[c_idx] == "H":
            continue
        hyb_c = hyb[c_idx]
        is_metal_c = _is_metal(symbols[c_idx])
        if not is_metal_c and hyb_c not in ("sp2", "sp3"):
            continue
        heavy_nbrs = [(j, v, d) for (j, v, d) in nbrs[c_idx] if symbols[j] != "H"]
        if len(heavy_nbrs) != 3:
            continue
        nb_data = [(symbols[j], hyb[j] if symbols[j] != "H" else "*",
                    j, v, d) for (j, v, d) in heavy_nbrs]
        nb_data.sort(key=lambda t: (t[0], t[1]))
        vec_n = [t[3] for t in nb_data]
        r1, r2, r3 = vec_n[0], vec_n[1], vec_n[2]
        b1 = -r1
        b2 = r2
        b3 = r3 - r2
        n1v = np.cross(b1, b2)
        n2v = np.cross(b2, b3)
        b2_norm = float(np.linalg.norm(b2))
        if b2_norm < 1e-6:
            continue
        nn1 = float(np.linalg.norm(n1v))
        nn2 = float(np.linalg.norm(n2v))
        if nn1 < 1e-6 or nn2 < 1e-6:
            continue
        cosphi = float(np.dot(n1v, n2v) / (nn1 * nn2))
        sinphi = float(np.dot(np.cross(n1v, n2v), b2 / b2_norm) / (nn1 * nn2))
        phi_deg = math.degrees(math.atan2(sinphi, cosphi))
        if not is_metal_c:
            key = centre_keys[c_idx]
            if key is not None:
                impropers_out[key].append(float(phi_deg))
        nbr_zs = [t[0] for t in nb_data]
        k_ip = _improper_pair_key(symbols[c_idx], hyb_c, nbr_zs)
        improper_pairs_out[k_ip].append(float(phi_deg))

    # ---- v3 torsions ----
    for b_idx in range(n):
        if symbols[b_idx] == "H" or _is_metal(symbols[b_idx]):
            continue
        nb_b = nbrs[b_idx]
        for (c_idx, _vc, _dc) in nb_b:
            if c_idx <= b_idx:
                continue
            if symbols[c_idx] == "H" or _is_metal(symbols[c_idx]):
                continue
            vec_bc = None
            for (jx, vx, _dx) in nb_b:
                if jx == c_idx:
                    vec_bc = vx
                    break
            if vec_bc is None:
                continue
            for (a_idx, vec_ba_raw, _da) in nb_b:
                if a_idx == c_idx:
                    continue
                if symbols[a_idx] == "H" or _is_metal(symbols[a_idx]):
                    continue
                vec_ba = -vec_ba_raw
                for (d_idx, vec_cd, _dd) in nbrs[c_idx]:
                    if d_idx == b_idx or d_idx == a_idx:
                        continue
                    if symbols[d_idx] == "H" or _is_metal(symbols[d_idx]):
                        continue
                    b1 = vec_ba
                    b2 = vec_bc
                    b3 = vec_cd
                    nn1v = np.cross(b1, b2)
                    nn2v = np.cross(b2, b3)
                    b2_norm = float(np.linalg.norm(b2))
                    if b2_norm < 1e-6:
                        continue
                    nn1m = float(np.linalg.norm(nn1v))
                    nn2m = float(np.linalg.norm(nn2v))
                    if nn1m < 1e-6 or nn2m < 1e-6:
                        continue
                    cosphi = float(np.dot(nn1v, nn2v) / (nn1m * nn2m))
                    sinphi = float(np.dot(np.cross(nn1v, nn2v),
                                           b2 / b2_norm) / (nn1m * nn2m))
                    phi_deg = math.degrees(math.atan2(sinphi, cosphi))
                    ring_bc = _bond_rmin(b_idx, c_idx)
                    arom_b = bool(is_aromatic[b_idx])
                    arom_c = bool(is_aromatic[c_idx])
                    parsed = _torsion_canonicalise(
                        symbols[a_idx], symbols[b_idx], hyb[b_idx],
                        symbols[c_idx], hyb[c_idx], symbols[d_idx],
                        int(ring_bc), arom_b, arom_c,
                    )
                    torsions_out[_to_torsion_key_str(parsed)].append(float(phi_deg))

    # ---- Lever 2 — TM categories (verbatim v5 logic) ----
    metal_indices = [i for i in range(n) if _is_metal(symbols[i])]

    for mi in metal_indices:
        m_sym = symbols[mi]
        m_block = _block_of(m_sym)
        for (j, _, d) in nbrs[mi]:
            if symbols[j] != "C":
                continue
            if not (1.7 <= d <= 2.15):
                continue
            if is_aromatic[j]:
                continue
            if ring_min_atom[j] in (5, 6):
                continue
            c_heavy_nonM = sum(
                1 for (k, _, _) in nbrs[j]
                if symbols[k] != "H" and not _is_metal(symbols[k])
            )
            if heavy_deg[j] == 3 and c_heavy_nonM == 2:
                k = json.dumps([m_sym, m_block, "C", "carbene"],
                               separators=(",", ":"), ensure_ascii=False)
                tm_carbene_out[k].append(float(d))

    for mi in metal_indices:
        m_sym = symbols[mi]
        m_block = _block_of(m_sym)
        c_nbrs = [(j, d) for (j, _, d) in nbrs[mi]
                  if symbols[j] == "C" and _HAPTO_MC_LO <= d <= _HAPTO_MC_HI]
        for ia in range(len(c_nbrs)):
            for ib in range(ia + 1, len(c_nbrs)):
                ja, da = c_nbrs[ia]
                jb, db = c_nbrs[ib]
                bonded = any(k == jb for (k, _, _) in nbrs[ja])
                if not bonded:
                    continue
                if hyb[ja] != "sp2" or hyb[jb] != "sp2":
                    continue
                ring_ja = ring_min_atom[ja]
                ring_jb = ring_min_atom[jb]
                if ring_ja in (5, 6) and ring_jb == ring_ja:
                    continue
                k = json.dumps([m_sym, m_block, "C", "eta2"],
                               separators=(",", ":"), ensure_ascii=False)
                tm_eta2_out[k].append(float((da + db) / 2.0))

    for ring_len in (5, 6):
        for ring in rings_of_size[ring_len]:
            ring_set = set(ring)
            if not all(symbols[r] == "C" for r in ring):
                continue
            if not all(hyb[r] == "sp2" for r in ring):
                continue
            for mi in metal_indices:
                m_sym = symbols[mi]
                m_block = _block_of(m_sym)
                d_list = []
                ok = True
                for r in ring:
                    d_mr = None
                    for (j, _, d) in nbrs[mi]:
                        if j == r:
                            d_mr = d
                            break
                    if d_mr is None or not (_HAPTO_MC_LO <= d_mr <= _HAPTO_MC_HI):
                        ok = False
                        break
                    d_list.append(d_mr)
                if ok and d_list:
                    cat = "eta5" if ring_len == 5 else "eta6"
                    k = json.dumps([m_sym, m_block, "C", cat],
                                   separators=(",", ":"), ensure_ascii=False)
                    if cat == "eta5":
                        tm_eta5_out[k].extend(float(x) for x in d_list)
                    else:
                        tm_eta6_out[k].extend(float(x) for x in d_list)

    for xi in range(n):
        s = symbols[xi]
        if s not in ("F", "Cl", "Br", "I", "O", "S", "N"):
            continue
        metal_nbrs = [(j, d) for (j, _, d) in nbrs[xi] if _is_metal(symbols[j])]
        if len(metal_nbrs) < 2:
            continue
        for (mj, d_mx) in metal_nbrs:
            m_sym = symbols[mj]
            m_block = _block_of(m_sym)
            k = json.dumps([m_sym, m_block, s, "mu"],
                           separators=(",", ":"), ensure_ascii=False)
            tm_mu_bridge_out[k].append(float(d_mx))

    for mi in metal_indices:
        m_sym = symbols[mi]
        m_block = _block_of(m_sym)
        for (hj, _, d_mh) in nbrs[mi]:
            if symbols[hj] != "H":
                continue
            if not (_AGOSTIC_MH_LO <= d_mh <= _AGOSTIC_MH_HI):
                continue
            for (cj, _, d_hc) in nbrs[hj]:
                if cj == mi:
                    continue
                if symbols[cj] != "C":
                    continue
                d_mc = None
                for (k, _, d_test) in nbrs[mi]:
                    if k == cj:
                        d_mc = d_test
                        break
                if d_mc is None:
                    diff = positions_np[cj] - positions_np[mi]
                    d_mc = float(np.linalg.norm(diff))
                if not (_AGOSTIC_MC_LO <= d_mc <= _AGOSTIC_MC_HI):
                    continue
                k = json.dumps([m_sym, m_block, "H", "agostic"],
                               separators=(",", ":"), ensure_ascii=False)
                tm_agostic_out[k].append(float(d_mh))

    for mi in metal_indices:
        m_sym = symbols[mi]
        if m_sym not in _D8_SP_CANDIDATES:
            continue
        heavy_nbrs = [(j, v, d) for (j, v, d) in nbrs[mi] if symbols[j] != "H"]
        if len(heavy_nbrs) < 4:
            continue
        m_block = _block_of(m_sym)
        heavy_nbrs.sort(key=lambda t: t[2])
        plane_4 = heavy_nbrs[:4]
        D = np.array([v for (_, v, _) in plane_4])
        try:
            _, _, Vt = np.linalg.svd(D, full_matrices=False)
        except Exception:
            continue
        normal = Vt[-1]
        proj = D @ normal
        if float(np.mean(proj * proj)) > 0.5:
            continue
        for (j, v, d) in heavy_nbrs[4:]:
            v_norm = float(np.linalg.norm(v))
            if v_norm < 1e-6:
                continue
            cos_axial = float(np.dot(v, normal) / v_norm)
            if abs(cos_axial) < 0.85:
                continue
            k = json.dumps([m_sym, m_block, symbols[j], "ax"],
                           separators=(",", ":"), ensure_ascii=False)
            tm_ox_addition_out[k].append(float(d))

    audit["fail_reason"] = "ok"
    return {
        "bonds": dict(bonds_out),
        "angles": dict(angles_out),
        "impropers": dict(impropers_out),
        "torsions": dict(torsions_out),
        "pair_bonds": dict(pair_bonds_out),
        "triple_angles": dict(triple_angles_out),
        "improper_pairs": dict(improper_pairs_out),
        "tm_carbene": dict(tm_carbene_out),
        "tm_eta2": dict(tm_eta2_out),
        "tm_eta5": dict(tm_eta5_out),
        "tm_eta6": dict(tm_eta6_out),
        "tm_mu_bridge": dict(tm_mu_bridge_out),
        "tm_agostic": dict(tm_agostic_out),
        "tm_ox_addition": dict(tm_ox_addition_out),
        "tm_pair_bond_block": dict(tm_pair_bond_block_out),
        "tm_triple_angle_block": dict(tm_triple_angle_block_out),
        "audit": audit,
        "stats": {"n_atoms": n,
                  "n_heavy_bonds": sum(len(a) for a in adj_heavy) // 2,
                  "n_metals": len(metal_indices)},
    }


def _process_chunk(paths_chunk: List[str]):
    bonds: Dict[str, List[float]] = defaultdict(list)
    angles: Dict[str, List[float]] = defaultdict(list)
    impropers: Dict[str, List[float]] = defaultdict(list)
    torsions: Dict[str, List[float]] = defaultdict(list)
    pair_bonds: Dict[str, List[float]] = defaultdict(list)
    triple_angles: Dict[str, List[float]] = defaultdict(list)
    improper_pairs: Dict[str, List[float]] = defaultdict(list)
    tm_carbene: Dict[str, List[float]] = defaultdict(list)
    tm_eta2: Dict[str, List[float]] = defaultdict(list)
    tm_eta5: Dict[str, List[float]] = defaultdict(list)
    tm_eta6: Dict[str, List[float]] = defaultdict(list)
    tm_mu_bridge: Dict[str, List[float]] = defaultdict(list)
    tm_agostic: Dict[str, List[float]] = defaultdict(list)
    tm_ox_addition: Dict[str, List[float]] = defaultdict(list)
    tm_pair_bond_block: Dict[str, List[float]] = defaultdict(list)
    tm_triple_angle_block: Dict[str, List[float]] = defaultdict(list)
    n_done = 0
    audit_counts: Dict[str, int] = defaultdict(int)
    for p in paths_chunk:
        out = _extract_one_entry(p)
        n_done += 1
        if out is None:
            audit_counts["read_fail"] += 1
            continue
        a = out.get("audit", {})
        reason = a.get("fail_reason", "")
        if reason and reason != "ok":
            audit_counts[reason] += 1
            continue
        audit_counts["ok"] += 1
        for k, vs in out["bonds"].items(): bonds[k].extend(vs)
        for k, vs in out["angles"].items(): angles[k].extend(vs)
        for k, vs in out["impropers"].items(): impropers[k].extend(vs)
        for k, vs in out["torsions"].items(): torsions[k].extend(vs)
        for k, vs in out["pair_bonds"].items(): pair_bonds[k].extend(vs)
        for k, vs in out["triple_angles"].items(): triple_angles[k].extend(vs)
        for k, vs in out["improper_pairs"].items(): improper_pairs[k].extend(vs)
        for k, vs in out["tm_carbene"].items(): tm_carbene[k].extend(vs)
        for k, vs in out["tm_eta2"].items(): tm_eta2[k].extend(vs)
        for k, vs in out["tm_eta5"].items(): tm_eta5[k].extend(vs)
        for k, vs in out["tm_eta6"].items(): tm_eta6[k].extend(vs)
        for k, vs in out["tm_mu_bridge"].items(): tm_mu_bridge[k].extend(vs)
        for k, vs in out["tm_agostic"].items(): tm_agostic[k].extend(vs)
        for k, vs in out["tm_ox_addition"].items(): tm_ox_addition[k].extend(vs)
        for k, vs in out["tm_pair_bond_block"].items(): tm_pair_bond_block[k].extend(vs)
        for k, vs in out["tm_triple_angle_block"].items(): tm_triple_angle_block[k].extend(vs)
    return (
        dict(bonds), dict(angles), dict(impropers), dict(torsions),
        dict(pair_bonds), dict(triple_angles), dict(improper_pairs),
        dict(tm_carbene), dict(tm_eta2), dict(tm_eta5), dict(tm_eta6),
        dict(tm_mu_bridge), dict(tm_agostic), dict(tm_ox_addition),
        dict(tm_pair_bond_block), dict(tm_triple_angle_block),
        n_done, dict(audit_counts),
    )


# ============================================================================
# Aggregation (verbatim from v5)
# ============================================================================

def _robust_aggregate(values: List[float]):
    if not values:
        return (np.nan, np.nan, 0, np.nan, np.nan, np.nan)
    arr = np.asarray(values, dtype=np.float64)
    n = arr.size
    p5, p50, p95 = np.percentile(arr, [5.0, 50.0, 95.0])
    mu = float(p50)
    if n >= 2:
        mad = float(np.median(np.abs(arr - mu)))
        sigma = max(1.4826 * mad, 1e-4)
        if mad == 0.0:
            sigma = max(float(np.std(arr)), 1e-4)
    else:
        sigma = 1e-3
    return (mu, float(sigma), int(n), float(p5), float(p50), float(p95))


def _fit_torsion_gmm(values: np.ndarray, max_components: int = 3) -> dict:
    from sklearn.mixture import GaussianMixture
    GMM_SUBSAMPLE_MAX = 100_000
    arr = np.asarray(values, dtype=np.float64).reshape(-1, 1)
    n = arr.shape[0]
    if n < 2:
        return {
            "n_components": 1,
            "pi": [1.0, np.nan, np.nan][:max_components],
            "mu": [float(arr.flatten()[0]), np.nan, np.nan][:max_components],
            "sigma": [1.0, np.nan, np.nan][:max_components],
        }
    if n > GMM_SUBSAMPLE_MAX:
        rng = np.random.default_rng(42)
        idx = rng.choice(n, size=GMM_SUBSAMPLE_MAX, replace=False)
        idx.sort()
        fit_arr = arr[idx]
    else:
        fit_arr = arr
    n_fit = fit_arr.shape[0]
    best_bic = np.inf
    best_gm = None
    for k in range(1, max_components + 1):
        if n_fit < k * 2:
            break
        try:
            gm = GaussianMixture(
                n_components=k,
                covariance_type="full",
                random_state=42,
                n_init=3,
                max_iter=200,
                reg_covar=1e-4,
            )
            gm.fit(fit_arr)
            bic = gm.bic(fit_arr)
        except Exception:
            continue
        if bic < best_bic:
            best_bic = bic
            best_gm = gm
    if best_gm is None:
        return {
            "n_components": 1,
            "pi": [1.0] + [np.nan] * (max_components - 1),
            "mu": [float(np.median(arr))] + [np.nan] * (max_components - 1),
            "sigma": [max(float(np.std(arr)), 1.0)] + [np.nan] * (max_components - 1),
        }
    mus = best_gm.means_.flatten()
    pis = best_gm.weights_.flatten()
    cov = best_gm.covariances_
    sigmas = np.sqrt(cov.reshape(-1)).flatten()
    order = np.argsort(mus)
    mus = mus[order]
    pis = pis[order]
    sigmas = sigmas[order]
    n_k = len(mus)
    pad_pi = list(pis) + [np.nan] * (max_components - n_k)
    pad_mu = list(mus) + [np.nan] * (max_components - n_k)
    pad_sig = list(sigmas) + [np.nan] * (max_components - n_k)
    return {
        "n_components": int(n_k),
        "pi": pad_pi[:max_components],
        "mu": pad_mu[:max_components],
        "sigma": pad_sig[:max_components],
    }


def _agg_table(vals_dict: Dict[str, List[float]]):
    keys = sorted(vals_dict.keys())
    n_k = len(keys)
    mu = np.full(n_k, np.nan, dtype=np.float64)
    sigma = np.full(n_k, np.nan, dtype=np.float64)
    n_arr = np.zeros(n_k, dtype=np.int32)
    for i, k in enumerate(keys):
        m, s, nn, _, _, _ = _robust_aggregate(vals_dict[k])
        mu[i] = m; sigma[i] = s; n_arr[i] = nn
    return (np.array(keys, dtype=np.object_), mu, sigma, n_arr)


# ============================================================================
# Driver
# ============================================================================

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--cod-dir", type=Path, default=COD_CLEAN_DIR,
                        help="Directory of clean COD XYZ files")
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--limit", type=int, default=0,
                        help="Use only the first N entries (smoke build).")
    parser.add_argument("--chunksize", type=int, default=128)
    parser.add_argument("--log", type=Path, default=DEFAULT_LOG)
    args = parser.parse_args()

    global log
    log = _setup_logging(args.log)

    log.info("=" * 72)
    log.info("GRIP-LIB v6 BUILD (COD, open-access twin of v5)")
    log.info("=" * 72)

    log.info(f"Listing COD entries from {args.cod_dir}")
    t0 = time.time()
    cod_paths = sorted(str(p) for p in args.cod_dir.iterdir() if p.suffix == ".xyz")
    n_total = len(cod_paths)
    log.info(f"  found {n_total} COD .xyz files in {time.time()-t0:.1f}s")
    if args.limit > 0:
        cod_paths = cod_paths[:args.limit]
        log.info(f"  --limit {args.limit} -> using {len(cod_paths)} entries")
    n_to_use = len(cod_paths)

    log.info(f"GRIP v6 build start workers={args.workers} "
             f"chunksize={args.chunksize} out={args.out}")

    bond_vals: Dict[str, List[float]] = defaultdict(list)
    angle_vals: Dict[str, List[float]] = defaultdict(list)
    improper_vals: Dict[str, List[float]] = defaultdict(list)
    torsion_vals: Dict[str, List[float]] = defaultdict(list)
    pair_bond_vals: Dict[str, List[float]] = defaultdict(list)
    triple_angle_vals: Dict[str, List[float]] = defaultdict(list)
    improper_pair_vals: Dict[str, List[float]] = defaultdict(list)
    tm_carbene_vals: Dict[str, List[float]] = defaultdict(list)
    tm_eta2_vals: Dict[str, List[float]] = defaultdict(list)
    tm_eta5_vals: Dict[str, List[float]] = defaultdict(list)
    tm_eta6_vals: Dict[str, List[float]] = defaultdict(list)
    tm_mu_bridge_vals: Dict[str, List[float]] = defaultdict(list)
    tm_agostic_vals: Dict[str, List[float]] = defaultdict(list)
    tm_ox_addition_vals: Dict[str, List[float]] = defaultdict(list)
    tm_pair_bond_block_vals: Dict[str, List[float]] = defaultdict(list)
    tm_triple_angle_block_vals: Dict[str, List[float]] = defaultdict(list)

    audit_counts_global: Dict[str, int] = defaultdict(int)

    chunks = [
        cod_paths[i:i + args.chunksize]
        for i in range(0, n_to_use, args.chunksize)
    ]
    n_chunks = len(chunks)
    log.info(f"  total chunks: {n_chunks} (~{args.chunksize} entries each)")

    t0 = time.time()
    last_log = t0
    n_done = 0

    ctx = mp.get_context("forkserver")
    with ctx.Pool(processes=args.workers, initializer=_worker_init) as pool:
        try:
            iterator = pool.imap_unordered(_process_chunk, chunks, chunksize=1)
            for chunk_i, result in enumerate(iterator, start=1):
                (bonds, angles, impropers, torsions,
                 pair_bonds, triple_angles, improper_pairs,
                 tm_carb, tm_e2, tm_e5, tm_e6,
                 tm_mub, tm_ago, tm_ox,
                 tm_pbb, tm_tab,
                 c_done, c_audit) = result
                for k, vs in bonds.items(): bond_vals[k].extend(vs)
                for k, vs in angles.items(): angle_vals[k].extend(vs)
                for k, vs in impropers.items(): improper_vals[k].extend(vs)
                for k, vs in torsions.items(): torsion_vals[k].extend(vs)
                for k, vs in pair_bonds.items(): pair_bond_vals[k].extend(vs)
                for k, vs in triple_angles.items(): triple_angle_vals[k].extend(vs)
                for k, vs in improper_pairs.items(): improper_pair_vals[k].extend(vs)
                for k, vs in tm_carb.items(): tm_carbene_vals[k].extend(vs)
                for k, vs in tm_e2.items(): tm_eta2_vals[k].extend(vs)
                for k, vs in tm_e5.items(): tm_eta5_vals[k].extend(vs)
                for k, vs in tm_e6.items(): tm_eta6_vals[k].extend(vs)
                for k, vs in tm_mub.items(): tm_mu_bridge_vals[k].extend(vs)
                for k, vs in tm_ago.items(): tm_agostic_vals[k].extend(vs)
                for k, vs in tm_ox.items(): tm_ox_addition_vals[k].extend(vs)
                for k, vs in tm_pbb.items(): tm_pair_bond_block_vals[k].extend(vs)
                for k, vs in tm_tab.items(): tm_triple_angle_block_vals[k].extend(vs)
                n_done += c_done
                for k, v in c_audit.items():
                    audit_counts_global[k] += v

                if (chunk_i % 50 == 0) or (time.time() - last_log) > 60.0:
                    rate = n_done / max(time.time() - t0, 1e-6)
                    eta_s = (n_to_use - n_done) / max(rate, 1e-6)
                    log.info(
                        f"  processed {n_done}/{n_to_use} "
                        f"({100*n_done/max(n_to_use,1):.1f}%) "
                        f"chunks={chunk_i}/{n_chunks} "
                        f"rate={rate:.0f} entries/s "
                        f"ETA={eta_s/60.0:.1f}min "
                        f"ok={audit_counts_global.get('ok',0)} "
                        f"pb_keys={len(pair_bond_vals)} "
                        f"carbene={len(tm_carbene_vals)} "
                        f"eta5={len(tm_eta5_vals)} "
                        f"mubr={len(tm_mu_bridge_vals)} "
                        f"agost={len(tm_agostic_vals)}"
                    )
                    last_log = time.time()
        except KeyboardInterrupt:
            log.warning("Interrupted; saving partial state")
            pool.terminate()

    log.info(f"Extraction done in {time.time()-t0:.1f}s: {n_done} entries scanned")
    log.info("Audit breakdown:")
    for r, c in sorted(audit_counts_global.items(), key=lambda kv: -kv[1]):
        log.info(f"  {r:30s} {c:8d}  ({100*c/max(n_done,1):.2f}%)")

    n_ok = audit_counts_global.get("ok", 0)

    # ---- v3 fragment fallback expansion ----
    log.info("Expanding v3 fragment values into fallback bins")
    t0 = time.time()
    def _expand(values_dict: Dict[str, List[float]]) -> Dict[str, List[float]]:
        out = defaultdict(list)
        for k, vals in values_dict.items():
            try:
                parsed = json.loads(k)
            except Exception:
                out[k].extend(vals)
                continue
            levels = _fallback_levels(parsed)
            for lvl in levels:
                out[lvl].extend(vals)
        return out

    bond_bins = _expand(bond_vals)
    angle_bins = _expand(angle_vals)
    improper_bins = _expand(improper_vals)
    log.info(f"  v3 expansion done in {time.time()-t0:.1f}s")

    log.info("Expanding torsion values into fallback bins")
    t0 = time.time()
    torsion_bins: Dict[str, List[float]] = defaultdict(list)
    for k, vals in torsion_vals.items():
        try:
            parsed = json.loads(k)
        except Exception:
            torsion_bins[k].extend(vals)
            continue
        for lvl in _torsion_fallback_levels(parsed):
            torsion_bins[lvl].extend(vals)
    log.info(f"  torsion expansion done in {time.time()-t0:.1f}s "
             f"({len(torsion_bins)} bins)")

    # ---- v3 master aggregation ----
    log.info("Aggregating master v3 bins")
    t0 = time.time()
    master_keys = sorted(
        set(bond_bins.keys()) | set(angle_bins.keys()) | set(improper_bins.keys())
    )
    n_master = len(master_keys)
    key_to_idx = {k: i for i, k in enumerate(master_keys)}

    bond_mu = np.full(n_master, np.nan, dtype=np.float64)
    bond_sigma = np.full(n_master, np.nan, dtype=np.float64)
    bond_n = np.zeros(n_master, dtype=np.int32)
    bond_p5 = np.full(n_master, np.nan, dtype=np.float32)
    bond_p50 = np.full(n_master, np.nan, dtype=np.float32)
    bond_p95 = np.full(n_master, np.nan, dtype=np.float32)
    angle_mu = np.full(n_master, np.nan, dtype=np.float64)
    angle_sigma = np.full(n_master, np.nan, dtype=np.float64)
    angle_n = np.zeros(n_master, dtype=np.int32)
    angle_p5 = np.full(n_master, np.nan, dtype=np.float32)
    angle_p50 = np.full(n_master, np.nan, dtype=np.float32)
    angle_p95 = np.full(n_master, np.nan, dtype=np.float32)
    improper_mu = np.full(n_master, np.nan, dtype=np.float64)
    improper_sigma = np.full(n_master, np.nan, dtype=np.float64)
    improper_n = np.zeros(n_master, dtype=np.int32)
    improper_p5 = np.full(n_master, np.nan, dtype=np.float32)
    improper_p50 = np.full(n_master, np.nan, dtype=np.float32)
    improper_p95 = np.full(n_master, np.nan, dtype=np.float32)

    for i, k in enumerate(master_keys):
        if k in bond_bins:
            mu, sigma, n, p5, p50, p95 = _robust_aggregate(bond_bins[k])
            bond_mu[i] = mu; bond_sigma[i] = sigma; bond_n[i] = n
            bond_p5[i] = p5; bond_p50[i] = p50; bond_p95[i] = p95
        if k in angle_bins:
            mu, sigma, n, p5, p50, p95 = _robust_aggregate(angle_bins[k])
            angle_mu[i] = mu; angle_sigma[i] = sigma; angle_n[i] = n
            angle_p5[i] = p5; angle_p50[i] = p50; angle_p95[i] = p95
        if k in improper_bins:
            mu, sigma, n, p5, p50, p95 = _robust_aggregate(improper_bins[k])
            improper_mu[i] = mu; improper_sigma[i] = sigma; improper_n[i] = n
            improper_p5[i] = p5; improper_p50[i] = p50; improper_p95[i] = p95
    log.info(f"  v3 master aggregation done in {time.time()-t0:.1f}s "
             f"({n_master} keys)")

    # ---- v4 pair tables ----
    log.info("Aggregating v4 pair tables")
    t0 = time.time()
    (pair_bond_keys, pair_bond_mu, pair_bond_sigma, pair_bond_n) = _agg_table(pair_bond_vals)
    (triple_angle_keys, triple_angle_mu, triple_angle_sigma, triple_angle_n) = _agg_table(triple_angle_vals)
    (improper_pair_keys, improper_pair_mu, improper_pair_sigma, improper_pair_n) = _agg_table(improper_pair_vals)
    log.info(f"  v4 pair aggregation done in {time.time()-t0:.1f}s: "
             f"pb={len(pair_bond_keys)} ta={len(triple_angle_keys)} "
             f"ip={len(improper_pair_keys)}")

    # ---- v5 TM-category and block-disagg tables ----
    log.info("Aggregating v5 TM-category tables")
    t0 = time.time()
    (tm_carb_keys, tm_carb_mu, tm_carb_sigma, tm_carb_n) = _agg_table(tm_carbene_vals)
    (tm_e2_keys, tm_e2_mu, tm_e2_sigma, tm_e2_n) = _agg_table(tm_eta2_vals)
    (tm_e5_keys, tm_e5_mu, tm_e5_sigma, tm_e5_n) = _agg_table(tm_eta5_vals)
    (tm_e6_keys, tm_e6_mu, tm_e6_sigma, tm_e6_n) = _agg_table(tm_eta6_vals)
    (tm_mub_keys, tm_mub_mu, tm_mub_sigma, tm_mub_n) = _agg_table(tm_mu_bridge_vals)
    (tm_ago_keys, tm_ago_mu, tm_ago_sigma, tm_ago_n) = _agg_table(tm_agostic_vals)
    (tm_ox_keys, tm_ox_mu, tm_ox_sigma, tm_ox_n) = _agg_table(tm_ox_addition_vals)
    (tm_pbb_keys, tm_pbb_mu, tm_pbb_sigma, tm_pbb_n) = _agg_table(tm_pair_bond_block_vals)
    (tm_tab_keys, tm_tab_mu, tm_tab_sigma, tm_tab_n) = _agg_table(tm_triple_angle_block_vals)
    log.info(f"  v5 TM-cat aggregation done in {time.time()-t0:.1f}s: "
             f"carb={len(tm_carb_keys)} e2={len(tm_e2_keys)} "
             f"e5={len(tm_e5_keys)} e6={len(tm_e6_keys)} "
             f"mu_br={len(tm_mub_keys)} agost={len(tm_ago_keys)} "
             f"ox={len(tm_ox_keys)} "
             f"pb_block={len(tm_pbb_keys)} ta_block={len(tm_tab_keys)}")

    # ---- v3 torsion aggregation (with GMM) ----
    log.info("Aggregating torsion bins (with per-bin GMM fits)")
    t0 = time.time()
    torsion_keys = sorted(torsion_bins.keys())
    n_torsion = len(torsion_keys)
    log.info(f"  torsion bin count: {n_torsion}")
    MAX_COMPONENTS = 3
    torsion_n = np.zeros(n_torsion, dtype=np.int32)
    torsion_n_components = np.ones(n_torsion, dtype=np.int32)
    torsion_pi = np.full((n_torsion, MAX_COMPONENTS), np.nan, dtype=np.float32)
    torsion_mu = np.full((n_torsion, MAX_COMPONENTS), np.nan, dtype=np.float32)
    torsion_sigma = np.full((n_torsion, MAX_COMPONENTS), np.nan, dtype=np.float32)

    fits_n5 = 0
    progress_every = max(1, n_torsion // 20)
    for i, k in enumerate(torsion_keys):
        vals = np.asarray(torsion_bins[k], dtype=np.float64)
        n_val = int(vals.size)
        torsion_n[i] = n_val
        if n_val >= 5:
            gmm = _fit_torsion_gmm(vals, MAX_COMPONENTS)
            fits_n5 += 1
        else:
            mu_med = float(np.median(vals)) if n_val > 0 else 0.0
            sig = max(float(np.std(vals)), 5.0) if n_val > 1 else 30.0
            gmm = {"n_components": 1, "pi": [1.0, np.nan, np.nan],
                   "mu": [mu_med, np.nan, np.nan], "sigma": [sig, np.nan, np.nan]}
        torsion_n_components[i] = int(gmm["n_components"])
        for j in range(MAX_COMPONENTS):
            torsion_pi[i, j] = float(gmm["pi"][j])
            torsion_mu[i, j] = float(gmm["mu"][j])
            torsion_sigma[i, j] = float(gmm["sigma"][j])
        if (i + 1) % progress_every == 0:
            elapsed = time.time() - t0
            log.info(f"  ...GMM fit {i+1}/{n_torsion} "
                     f"({100*(i+1)/n_torsion:.1f}%) elapsed={elapsed:.1f}s")
    log.info(f"  torsion aggregation done in {time.time()-t0:.1f}s "
             f"({fits_n5}/{n_torsion} n>=5 GMM fits)")

    # ---- Fragment fallback CSR encoding (v3 schema) ----
    log.info("Encoding fragment fallback chains")
    t0 = time.time()
    orig_frag_keys = sorted(set(bond_vals.keys()) | set(angle_vals.keys()) |
                            set(improper_vals.keys()))
    n_orig = len(orig_frag_keys)
    MAX_CHAIN_LEN = 6
    fb_flat = np.full(n_orig * MAX_CHAIN_LEN, -1, dtype=np.int32)
    fb_ptr = np.zeros(n_orig + 1, dtype=np.int32)
    write_pos = 0
    for i, k in enumerate(orig_frag_keys):
        fb_ptr[i] = write_pos
        try:
            parsed = json.loads(k)
        except Exception:
            chain = [k]
        else:
            chain = _fallback_levels(parsed)
        for lvl in chain:
            j = key_to_idx.get(lvl, -1)
            if j >= 0:
                fb_flat[write_pos] = j
                write_pos += 1
    fb_ptr[n_orig] = write_pos
    fb_flat = fb_flat[:write_pos]
    log.info(f"  frag fallback chains: n_orig={n_orig} entries={write_pos} "
             f"in {time.time()-t0:.1f}s")

    log.info("Encoding torsion fallback chains")
    t0 = time.time()
    orig_torsion_keys = sorted(torsion_vals.keys())
    n_orig_t = len(orig_torsion_keys)
    torsion_key_to_idx = {k: i for i, k in enumerate(torsion_keys)}
    MAX_T_CHAIN = 6
    torsion_fb_flat = np.full(n_orig_t * MAX_T_CHAIN, -1, dtype=np.int32)
    torsion_fb_ptr = np.zeros(n_orig_t + 1, dtype=np.int32)
    write_pos_t = 0
    for i, k in enumerate(orig_torsion_keys):
        torsion_fb_ptr[i] = write_pos_t
        try:
            parsed = json.loads(k)
        except Exception:
            chain = [k]
        else:
            chain = _torsion_fallback_levels(parsed)
        for lvl in chain:
            j = torsion_key_to_idx.get(lvl, -1)
            if j >= 0:
                torsion_fb_flat[write_pos_t] = j
                write_pos_t += 1
    torsion_fb_ptr[n_orig_t] = write_pos_t
    torsion_fb_flat = torsion_fb_flat[:write_pos_t]
    log.info(f"  torsion fallback chains in {time.time()-t0:.1f}s")

    mogul_validated = np.zeros(n_master, dtype=bool)

    # ---- Coverage report ----
    bonds_present = int((bond_n > 0).sum())
    bonds_ge5 = int((bond_n >= 5).sum())
    angles_present = int((angle_n > 0).sum())
    angles_ge5 = int((angle_n >= 5).sum())
    imp_present = int((improper_n > 0).sum())
    imp_ge5 = int((improper_n >= 5).sum())
    tor_present = int((torsion_n > 0).sum())
    tor_ge5 = int((torsion_n >= 5).sum())
    pb_ge5 = int((pair_bond_n >= 5).sum())
    ta_ge5 = int((triple_angle_n >= 5).sum())
    ip_ge5 = int((improper_pair_n >= 5).sum())
    carb_ge5 = int((tm_carb_n >= 5).sum())
    e2_ge5 = int((tm_e2_n >= 5).sum())
    e5_ge5 = int((tm_e5_n >= 5).sum())
    e6_ge5 = int((tm_e6_n >= 5).sum())
    mub_ge5 = int((tm_mub_n >= 5).sum())
    ago_ge5 = int((tm_ago_n >= 5).sum())
    ox_ge5 = int((tm_ox_n >= 5).sum())
    pbb_ge5 = int((tm_pbb_n >= 5).sum())
    tab_ge5 = int((tm_tab_n >= 5).sum())

    log.info("=" * 72)
    log.info("GRIP-LIB v6 (COD) COVERAGE REPORT")
    log.info("=" * 72)
    log.info(f"COD scanned: {n_done}  ok: {n_ok}")
    log.info(f"v3 master keys: {n_master}  v3 orig frag keys: {n_orig}")
    log.info(f"  bonds   obs={bonds_present}/{n_master} n>=5={bonds_ge5}/{n_master}")
    log.info(f"  angles  obs={angles_present}/{n_master} n>=5={angles_ge5}/{n_master}")
    log.info(f"  imps    obs={imp_present}/{n_master} n>=5={imp_ge5}/{n_master}")
    log.info(f"v3 torsion bins: {n_torsion} obs={tor_present} n>=5={tor_ge5}")
    log.info(f"v4 pair tables: pb={len(pair_bond_keys)} n>=5={pb_ge5}, "
             f"ta={len(triple_angle_keys)} n>=5={ta_ge5}, "
             f"ip={len(improper_pair_keys)} n>=5={ip_ge5}")
    log.info(f"v5 TM categories (n>=5):")
    log.info(f"  carbene     {len(tm_carb_keys):5d} keys ({carb_ge5} n>=5)")
    log.info(f"  hapto_eta2  {len(tm_e2_keys):5d} keys ({e2_ge5} n>=5)")
    log.info(f"  hapto_eta5  {len(tm_e5_keys):5d} keys ({e5_ge5} n>=5)")
    log.info(f"  hapto_eta6  {len(tm_e6_keys):5d} keys ({e6_ge5} n>=5)")
    log.info(f"  mu_bridge   {len(tm_mub_keys):5d} keys ({mub_ge5} n>=5)")
    log.info(f"  agostic     {len(tm_ago_keys):5d} keys ({ago_ge5} n>=5)")
    log.info(f"  ox_addition {len(tm_ox_keys):5d} keys ({ox_ge5} n>=5)")
    log.info(f"v5 block-disagg: pb_block={len(tm_pbb_keys)} n>=5={pbb_ge5}, "
             f"ta_block={len(tm_tab_keys)} n>=5={tab_ge5}")
    log.info("=" * 72)

    log.info(f"Writing {args.out}")
    t0 = time.time()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        version=np.int32(6),
        source=np.array("COD", dtype=np.object_),
        cleaning_applied=np.int32(1),
        cod_n_entries_scanned=np.int32(n_done),
        cod_n_extracted=np.int32(n_ok),
        # v5-compatible audit fields (zero for v6 since COD is pre-cleaned)
        meta_n_total_scanned=np.int32(n_done),
        meta_n_extracted=np.int32(n_ok),
        meta_n_polymeric=np.int32(0),
        meta_n_no_3d=np.int32(0),
        meta_n_too_small=np.int32(
            audit_counts_global.get("too_few_nonmetals", 0)
            + audit_counts_global.get("too_few_atoms", 0)
        ),
        meta_n_recovered_disorder=np.int32(0),
        meta_n_recovered_partial_coord=np.int32(0),
        # v3 master
        n_master=np.int32(n_master),
        n_orig=np.int32(n_orig),
        n_torsion=np.int32(n_torsion),
        n_orig_torsion=np.int32(n_orig_t),
        keys=np.array(master_keys, dtype=np.object_),
        orig_keys=np.array(orig_frag_keys, dtype=np.object_),
        orig_idx=np.zeros(n_orig, dtype=np.int32),
        bond_mu=bond_mu, bond_sigma=bond_sigma, bond_n=bond_n,
        bond_p5=bond_p5, bond_p50=bond_p50, bond_p95=bond_p95,
        angle_mu=angle_mu, angle_sigma=angle_sigma, angle_n=angle_n,
        angle_p5=angle_p5, angle_p50=angle_p50, angle_p95=angle_p95,
        improper_mu=improper_mu, improper_sigma=improper_sigma, improper_n=improper_n,
        improper_p5=improper_p5, improper_p50=improper_p50, improper_p95=improper_p95,
        fb_flat=fb_flat, fb_ptr=fb_ptr,
        torsion_keys=np.array(torsion_keys, dtype=np.object_),
        torsion_orig_keys=np.array(orig_torsion_keys, dtype=np.object_),
        torsion_n=torsion_n, torsion_n_components=torsion_n_components,
        torsion_pi=torsion_pi, torsion_mu=torsion_mu, torsion_sigma=torsion_sigma,
        torsion_fb_flat=torsion_fb_flat, torsion_fb_ptr=torsion_fb_ptr,
        # v4 pair tables
        n_pair_bond=np.int32(len(pair_bond_keys)),
        n_triple_angle=np.int32(len(triple_angle_keys)),
        n_improper_pair=np.int32(len(improper_pair_keys)),
        pair_bond_keys=pair_bond_keys,
        pair_bond_mu=pair_bond_mu, pair_bond_sigma=pair_bond_sigma, pair_bond_n=pair_bond_n,
        triple_angle_keys=triple_angle_keys,
        triple_angle_mu=triple_angle_mu, triple_angle_sigma=triple_angle_sigma, triple_angle_n=triple_angle_n,
        improper_pair_keys=improper_pair_keys,
        improper_pair_mu=improper_pair_mu, improper_pair_sigma=improper_pair_sigma, improper_pair_n=improper_pair_n,
        # v5 TM-category tables
        tm_carbene_keys=tm_carb_keys,
        tm_carbene_mu=tm_carb_mu, tm_carbene_sigma=tm_carb_sigma, tm_carbene_n=tm_carb_n,
        tm_hapto_eta2_keys=tm_e2_keys,
        tm_hapto_eta2_mu=tm_e2_mu, tm_hapto_eta2_sigma=tm_e2_sigma, tm_hapto_eta2_n=tm_e2_n,
        tm_hapto_eta5_keys=tm_e5_keys,
        tm_hapto_eta5_mu=tm_e5_mu, tm_hapto_eta5_sigma=tm_e5_sigma, tm_hapto_eta5_n=tm_e5_n,
        tm_hapto_eta6_keys=tm_e6_keys,
        tm_hapto_eta6_mu=tm_e6_mu, tm_hapto_eta6_sigma=tm_e6_sigma, tm_hapto_eta6_n=tm_e6_n,
        tm_mu_bridge_keys=tm_mub_keys,
        tm_mu_bridge_mu=tm_mub_mu, tm_mu_bridge_sigma=tm_mub_sigma, tm_mu_bridge_n=tm_mub_n,
        tm_agostic_keys=tm_ago_keys,
        tm_agostic_mu=tm_ago_mu, tm_agostic_sigma=tm_ago_sigma, tm_agostic_n=tm_ago_n,
        tm_ox_addition_keys=tm_ox_keys,
        tm_ox_addition_mu=tm_ox_mu, tm_ox_addition_sigma=tm_ox_sigma, tm_ox_addition_n=tm_ox_n,
        # v5 block-disaggregated tables
        tm_pair_bond_block_keys=tm_pbb_keys,
        tm_pair_bond_block_mu=tm_pbb_mu, tm_pair_bond_block_sigma=tm_pbb_sigma, tm_pair_bond_block_n=tm_pbb_n,
        tm_triple_angle_block_keys=tm_tab_keys,
        tm_triple_angle_block_mu=tm_tab_mu, tm_triple_angle_block_sigma=tm_tab_sigma, tm_triple_angle_block_n=tm_tab_n,
        mogul_validated=mogul_validated,
    )
    log.info(f"Wrote {args.out} ({args.out.stat().st_size/1e6:.1f} MB) "
             f"in {time.time()-t0:.1f}s")

    log.info("v6 BUILD DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
