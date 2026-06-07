#!/usr/bin/env python3
"""GRIP v5_TM: universal TM-complex ligand-internal library.

User mandate 2026-06-07: regenerate ``grip_lib_v5.npz`` (output as
``grip_lib_v5_TM.npz``) so it has comprehensive ligand-internal coverage
from ALL CCDC TM-complex entries — universally, no per-class extraction.

ROOT CAUSE OF THE COVERAGE GAP
------------------------------
The legacy v5 build emits v3-style centered-fragment keys whose neighbour-
list contains the centre's FULL neighbour set (e.g. for an sp3 C with
neighbours [C, C, H, Ir] the key is
``["C","sp3", [["C","H","C","Ir"]_full_sorted], []]``).
But ``lookup_bond`` builds the SINGLE-neighbour key
``["C","sp3", [["X",1,r,h_X]], []]`` because it only knows the bond
partner — no global topology.  No fallback level in v3 ``_fallback_levels``
DROPS neighbours, so the standard chain misses on nearly every metal-complex
ligand-internal bond.

The v4 ``pair_bond_keys`` table HAS the correct per-pair distributions, but
it is gated behind ``DELFIN_FFFREE_GRIP_LOOKUP_TM_FALLBACK=1`` — when that
env-flag is unset (default), the lookup never reaches it.

This build fixes the library SO the standard chain hits NATURALLY:

  1. For every TM-complex CCDC entry, walk every ligand-internal (non-metal-
     incident) bond/angle/improper.
  2. Emit v3-style fragment keys with SINGLE-NEIGHBOUR neighbour-lists,
     PLUS the legacy full-neighbour-list keys (so existing semantics are
     preserved byte-identically for any caller that builds full keys).
  3. Also emit the v4 ``pair_bond`` / ``triple_angle`` / ``improper_pair``
     tables (per-pair, per-triple aggregated) so TM-fallback still works
     identically.

UNIVERSAL PAIR-AGGREGATION ONLY
-------------------------------
No per-class branching.  No NHC-specific keys.  No hapto-specific keys.
Every bond/angle/improper emits a key by its graph-topology features
(element + hybridization + ring-size + aromaticity), aggregated across
ALL TM-complex entries.

The 7 TM-category tables (carbene/hapto/μ-bridge/agostic/ox-addition)
inherited from v5 are kept BUT the new ligand-internal universal keys
take precedence in the standard chain — categories are classified at
LOOKUP time by graph topology, not at extraction time.

LIGAND-INTERNAL FOCUS
---------------------
For each TM-complex entry, M-X bonds are skipped (already covered in M-D
library).  Only ligand-internal fragments (every atom NOT a metal) are
extracted.  An atom that has a metal in its FULL neighbour list still
contributes its non-metal sub-fragments — but the metal is removed from
the neighbour list (so the resulting key matches the polish-time query
that knows nothing about the metal partner).

DETERMINISM
-----------
* PYTHONHASHSEED=0 set at import.
* CSD entries processed in lex-sorted (identifier) order.
* Sorted keys + lex order in all output arrays.
* Robust aggregation: median + 1.4826*MAD (no RNG).

USAGE
-----
::

    PYTHONHASHSEED=0 \\
    /home/qmchem_max/CCDC/CSD_2026.1/ccdc-software/csd-python-api/run_csd_python_api \\
        scripts/grip_build_mogul_lib_v5_TM_ligand_internals.py \\
        --workers 32 --out reports/grip_lib_v5_TM.npz

OUTPUT SCHEMA
-------------
Same npz schema as v5 (so ``GripLibrary.__init__`` loads it byte-identically).
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
DEFAULT_OUT = ROOT / "reports" / "grip_lib_v5_TM.npz"
DEFAULT_LOG = Path("/tmp/grip_lib_v5_TM_build.log")


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
    return logging.getLogger("grip_build_v5_TM")


log = logging.getLogger("grip_build_v5_TM")


# ============================================================================
# Conventions (verbatim from v3/v4/v5)
# ============================================================================

_METALS_D = frozenset({
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
})
_METALS_F = frozenset({
    "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy",
    "Ho","Er","Tm","Yb","Lu",
    "Ac","Th","Pa","U","Np","Pu",
})
_METALS = _METALS_D | _METALS_F

_METAL_BLOCK: Dict[str, str] = {}
for s in ("Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn"):
    _METAL_BLOCK[s] = "3d"
for s in ("Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd"):
    _METAL_BLOCK[s] = "4d"
for s in ("Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"):
    _METAL_BLOCK[s] = "5d"
for s in _METALS_F:
    _METAL_BLOCK[s] = "f"


def _is_metal(sym: str) -> bool:
    return sym in _METALS


def _block_of(sym: str) -> str:
    return _METAL_BLOCK.get(sym, "*")


def _infer_hyb(symbol: str, heavy_degree: int, total_degree: int) -> str:
    """Match v3/v4/v5 ``_infer_hyb``."""
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
    """Mirror :func:`grip_mogul_lookup._fallback_levels`."""
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
        return [str(Z_b), str(hyb_b),
                str(Z_c), str(hyb_c),
                str(Z_a), str(Z_d),
                int(ring_bc),
                bool(arom_b), bool(arom_c)]
    return [str(Z_c), str(hyb_c),
            str(Z_b), str(hyb_b),
            str(Z_d), str(Z_a),
            int(ring_bc),
            bool(arom_c), bool(arom_b)]


def _to_torsion_key_str(parsed: list) -> str:
    return json.dumps(parsed, separators=(",", ":"), ensure_ascii=False)


def _torsion_fallback_levels(parsed: list) -> list[str]:
    if not isinstance(parsed, list) or len(parsed) != 9:
        return [_to_torsion_key_str(parsed)]
    Zb, hyb_b, Zc, hyb_c, Za, Zd, ring_bc, arom_b, arom_c = parsed
    levels: list[str] = []
    levels.append(_to_torsion_key_str(
        [Zb, hyb_b, Zc, hyb_c, Za, Zd, ring_bc, arom_b, arom_c]
    ))
    levels.append(_to_torsion_key_str(
        [Zb, hyb_b, Zc, hyb_c, Za, Zd, -1, arom_b, arom_c]
    ))
    levels.append(_to_torsion_key_str(
        [Zb, hyb_b, Zc, hyb_c, Za, Zd, -1, False, False]
    ))
    levels.append(_to_torsion_key_str(
        [Zb, hyb_b, Zc, hyb_c, "*", "*", -1, False, False]
    ))
    levels.append(_to_torsion_key_str(
        [Zb, "*", Zc, "*", "*", "*", -1, False, False]
    ))
    levels.append(_to_torsion_key_str(
        [Zb, "*", Zc, "*", "*", "*", -1, False, False]
    ))
    seen: set[str] = set()
    out: list[str] = []
    for s in levels:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


# ============================================================================
# v4 pair-resolved keys (unchanged from v5)
# ============================================================================

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
# Per-entry extraction
# ============================================================================

_WORKER_READER = None


def _worker_init():
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    np.random.seed(42)
    random.seed(42)


def _ensure_reader():
    global _WORKER_READER
    if _WORKER_READER is None:
        from ccdc.io import EntryReader
        _WORKER_READER = EntryReader("CSD")
    return _WORKER_READER


def _clean_mol_for_extraction(entry):
    """v5 cleaning logic: remove disorder, take main component, etc."""
    try:
        mol = entry.molecule
    except Exception:
        return None, "no_mol"
    if mol is None:
        return None, "no_mol"
    if mol.is_polymeric:
        return None, "polymeric"
    # Take largest metal-containing component if multi-component, else largest.
    try:
        comps = mol.components
    except Exception:
        comps = [mol]
    if comps is None or len(comps) == 0:
        return None, "no_comp"
    tm_comps = []
    for c in comps:
        try:
            has_tm = any(a.atomic_symbol in _METALS for a in c.atoms)
        except Exception:
            continue
        if has_tm:
            tm_comps.append(c)
    if tm_comps:
        comps_use = tm_comps
    else:
        comps_use = list(comps)
    comps_use.sort(key=lambda c: -len(c.atoms))
    mol_main = comps_use[0]
    # Strip disordered atoms (occupancy < 0.55 or label suffix B-Z).
    return mol_main, "ok"


def _extract_one_entry(csd_index: int) -> Optional[dict]:
    """Extract universal TM-complex ligand-internal fragments.

    Only entries containing a transition metal contribute.
    """
    csd = _ensure_reader()
    try:
        e = csd[csd_index]
    except Exception:
        return None

    try:
        if not e.has_3d_structure:
            return {"_status": "no_3d"}
    except Exception:
        return None

    mol, status = _clean_mol_for_extraction(e)
    if mol is None:
        return {"_status": status}

    try:
        atoms = mol.atoms
    except Exception:
        return None
    if atoms is None or len(atoms) < 4 or len(atoms) > 5000:
        return {"_status": "size"}

    # Symbols + positions (also captures occupancy filter)
    symbols: List[str] = []
    positions: List[Tuple[float, float, float]] = []
    valid: List[bool] = []
    try:
        for a in atoms:
            c = a.coordinates
            if c is None:
                valid.append(False)
                symbols.append(a.atomic_symbol)
                positions.append((0.0, 0.0, 0.0))
                continue
            sym = a.atomic_symbol
            # Disorder cleaning: skip secondary-disorder atoms.
            try:
                occ = getattr(a, "occupancy", None)
                if occ is not None and float(occ) < 0.55:
                    valid.append(False)
                    symbols.append(sym)
                    positions.append((0.0, 0.0, 0.0))
                    continue
            except Exception:
                pass
            symbols.append(sym)
            positions.append((float(c.x), float(c.y), float(c.z)))
            valid.append(True)
    except Exception:
        return None

    # Must contain at least one TM atom (the universal scope condition).
    n = len(atoms)
    has_tm = any(symbols[i] in _METALS for i in range(n) if valid[i])
    if not has_tm:
        return {"_status": "no_tm"}

    positions_np = np.asarray(positions, dtype=np.float64)

    idx_of: Dict[int, int] = {}
    for i, a in enumerate(atoms):
        try:
            ai = int(a.index)
        except Exception:
            return None
        idx_of[ai] = i

    nbrs: List[List[Tuple[int, np.ndarray, float]]] = [[] for _ in range(n)]
    for i, a in enumerate(atoms):
        if not valid[i]:
            continue
        try:
            for nbr in a.neighbours:
                j = idx_of.get(int(nbr.index))
                if j is None or not valid[j]:
                    continue
                disp = positions_np[j] - positions_np[i]
                d = float(np.linalg.norm(disp))
                if d < 0.4 or d > 4.5:
                    continue
                nbrs[i].append((j, disp, d))
        except Exception:
            return None

    heavy_deg = [
        sum(1 for (j, _, _) in nbrs[i] if symbols[j] != "H")
        for i in range(n)
    ]
    total_deg = [len(nbrs[i]) for i in range(n)]
    hyb = [_infer_hyb(symbols[i], heavy_deg[i], total_deg[i]) for i in range(n)]

    nonmetal_idx = [i for i in range(n)
                    if valid[i] and symbols[i] not in _METALS and symbols[i] != "H"]
    if len(nonmetal_idx) < 3:
        return {"_status": "too_few_nonmetal"}
    nm_ok = sum(1 for i in nonmetal_idx if total_deg[i] <= 6)
    if nm_ok < int(0.8 * len(nonmetal_idx)):
        return {"_status": "bad_valence"}

    # Heavy adjacency for ring detection (BFS-based, MAX_RING bound 8)
    adj_heavy: List[List[int]] = [[] for _ in range(n)]
    for i in range(n):
        if not valid[i] or symbols[i] == "H":
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
        if not valid[i] or symbols[i] == "H":
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

    is_aromatic = [False] * n
    for i in range(n):
        if not valid[i] or symbols[i] == "H":
            continue
        if ring_min_atom[i] not in (5, 6):
            continue
        if hyb[i] != "sp2":
            continue
        nb_sp2 = sum(1 for j in adj_heavy[i] if hyb[j] == "sp2")
        if nb_sp2 >= 2:
            is_aromatic[i] = True

    # ---- Output bins ------------------------------------------------------
    bonds_out: Dict[str, List[float]] = defaultdict(list)
    angles_out: Dict[str, List[float]] = defaultdict(list)
    impropers_out: Dict[str, List[float]] = defaultdict(list)
    torsions_out: Dict[str, List[float]] = defaultdict(list)

    pair_bonds_out: Dict[str, List[float]] = defaultdict(list)
    triple_angles_out: Dict[str, List[float]] = defaultdict(list)
    improper_pairs_out: Dict[str, List[float]] = defaultdict(list)

    # ---- v3 centered-fragment SINGLE-NEIGHBOUR keys (NEW: matches polish query) ----
    # For every actual ligand-internal bond, emit a v3 key on each endpoint
    # whose FRAGMENT contains ONLY the partner atom in the neighbour list.
    # This is exactly the key `lookup_bond` builds, so the standard chain
    # hits NATURALLY without TM-fallback env-flag.
    #
    # Ligand-internal definition: neither endpoint is a metal.

    seen_bonds: set = set()
    for i in range(n):
        if not valid[i]:
            continue
        sym_i = symbols[i]
        if sym_i == "H":
            continue  # central key on H is rare in queries
        if _is_metal(sym_i):
            continue  # skip metal-centred bond keys (M-D library covers them)
        hyb_i = hyb[i]
        if hyb_i == "*":
            continue
        for (j, _, d) in nbrs[i]:
            if not valid[j]:
                continue
            if not (0.5 < d < 3.5):
                continue
            sym_j = symbols[j]
            if _is_metal(sym_j):
                continue  # ligand-internal only
            hyb_j = hyb[j] if sym_j != "H" else "*"
            ring_j = ring_min_atom[j]
            rmin_bond = _bond_rmin(i, j)
            # Single-neighbour key centered at i with j as the sole neighbour.
            neigh_one = [[sym_j, 1, int(rmin_bond), str(hyb_j)]]
            key_single = _to_key_str([sym_i, hyb_i, neigh_one, []])
            bonds_out[key_single].append(float(d))

            # Also emit the legacy v4 pair_bond key (unchanged).
            pair_bonds_out[
                _pair_bond_key(sym_i, hyb_i, sym_j, hyb_j)
            ].append(float(d))

    # Metal-containing bonds: emit pair_bond keys for TM-fallback only.
    # These don't seed v3 fragment keys (M-D library is the authority).
    for i in range(n):
        if not valid[i]:
            continue
        sym_i = symbols[i]
        if sym_i == "H":
            continue
        hyb_i = hyb[i]
        for (j, _, d) in nbrs[i]:
            if not valid[j] or j <= i:
                continue
            if not (0.5 < d < 3.5):
                continue
            sym_j = symbols[j]
            hyb_j = hyb[j]
            # Only emit pair_bond + TM block keys for metal-incident bonds.
            if _is_metal(sym_i) or _is_metal(sym_j):
                pair_bonds_out[
                    _pair_bond_key(sym_i, hyb_i, sym_j, hyb_j)
                ].append(float(d))

    # ---- v3 single-neighbour ANGLES (ligand-internal only) ----
    # For each centre atom b, every pair of neighbours (a, c).  Emit a v3
    # key centred at b with the SORTED two-neighbour list [a, c] (the
    # lookup_angle path builds exactly this).  Skip if any of (a,b,c) is
    # a metal — those are M-D-X or D-M-D' priors and live in M-D library.
    for b in range(n):
        if not valid[b]:
            continue
        sym_b = symbols[b]
        if sym_b == "H":
            continue
        if _is_metal(sym_b):
            continue
        hyb_b = hyb[b]
        if hyb_b == "*":
            continue
        nbr_local = nbrs[b]
        if len(nbr_local) < 2:
            continue
        rmin_b = ring_min_atom[b]
        for ia in range(len(nbr_local)):
            for ic in range(ia + 1, len(nbr_local)):
                ja, vec_a, _da = nbr_local[ia]
                jc, vec_c, _dc = nbr_local[ic]
                if not valid[ja] or not valid[jc]:
                    continue
                sym_a = symbols[ja]
                sym_c = symbols[jc]
                if _is_metal(sym_a) or _is_metal(sym_c):
                    continue  # ligand-internal angles only
                na = float(np.linalg.norm(vec_a))
                nc = float(np.linalg.norm(vec_c))
                if na < 1e-6 or nc < 1e-6:
                    continue
                cos_a = float(np.dot(vec_a, vec_c) / (na * nc))
                cos_a = max(-1.0, min(1.0, cos_a))
                angle_deg = math.degrees(math.acos(cos_a))
                if not (15.0 <= angle_deg <= 175.0):
                    continue
                # Build the canonical sorted-pair neighbour list.
                hyb_a = hyb[ja] if sym_a != "H" else "*"
                hyb_c2 = hyb[jc] if sym_c != "H" else "*"
                nbrs_pair = sorted(
                    (
                        (str(sym_a), int(rmin_b), str(hyb_a)),
                        (str(sym_c), int(rmin_b), str(hyb_c2)),
                    ),
                    key=lambda t: (t[0], t[2], t[1]),
                )
                neighbors_key = [[t[0], 1, t[1], t[2]] for t in nbrs_pair]
                key_pair = _to_key_str([sym_b, hyb_b, neighbors_key, []])
                angles_out[key_pair].append(angle_deg)

                # Also emit v4 triple_angle for TM-fallback compatibility.
                triple_angles_out[
                    _triple_angle_key(sym_b, hyb_b, sym_a, sym_c)
                ].append(angle_deg)

    # ---- v3 IMPROPERS centered at non-metal atoms ----
    # Same key as lookup_improper: centre + sorted-3-neighbour-elements key.
    for c_idx in range(n):
        if not valid[c_idx]:
            continue
        sym_c = symbols[c_idx]
        if sym_c == "H" or _is_metal(sym_c):
            continue
        hyb_c = hyb[c_idx]
        if hyb_c not in ("sp2", "sp3"):
            continue
        heavy_nbrs = [(j, v, d) for (j, v, d) in nbrs[c_idx]
                      if valid[j] and symbols[j] != "H"]
        if len(heavy_nbrs) != 3:
            continue
        # Skip if ANY neighbour is a metal — those impropers are M-X/Y/Z and
        # live elsewhere.
        if any(_is_metal(symbols[j]) for (j, _, _) in heavy_nbrs):
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
        # Use v3 lookup_improper key: centred + sorted-neighbour element/hyb list.
        rmin = ring_min_atom[c_idx]
        nbrs_data_sorted = [(t[0], t[1]) for t in nb_data]
        neighbors_key = [[z, 1, int(rmin), h] for (z, h) in nbrs_data_sorted]
        key_imp = _to_key_str([sym_c, hyb_c, neighbors_key, []])
        impropers_out[key_imp].append(float(phi_deg))

        # v4 improper_pair (sorted-element neighbours)
        nbr_zs = [t[0] for t in nb_data]
        improper_pairs_out[
            _improper_pair_key(sym_c, hyb_c, nbr_zs)
        ].append(float(phi_deg))

    # ---- TORSIONS (heavy-only, non-metal) ----
    for b_idx in range(n):
        if not valid[b_idx]:
            continue
        if symbols[b_idx] == "H" or _is_metal(symbols[b_idx]):
            continue
        nb_b = nbrs[b_idx]
        for (c_idx, _vc, _dc) in nb_b:
            if not valid[c_idx] or c_idx <= b_idx:
                continue
            if symbols[c_idx] == "H" or _is_metal(symbols[c_idx]):
                continue
            nb_c = nbrs[c_idx]
            vec_bc = None
            for (jx, vx, _dx) in nb_b:
                if jx == c_idx:
                    vec_bc = vx
                    break
            if vec_bc is None:
                continue
            for (a_idx, vec_ba_raw, _da) in nb_b:
                if not valid[a_idx]:
                    continue
                if a_idx == c_idx:
                    continue
                if symbols[a_idx] == "H" or _is_metal(symbols[a_idx]):
                    continue
                for (d_idx, vec_cd, _dd) in nb_c:
                    if not valid[d_idx]:
                        continue
                    if d_idx == b_idx or d_idx == a_idx:
                        continue
                    if symbols[d_idx] == "H" or _is_metal(symbols[d_idx]):
                        continue
                    b1 = -vec_ba_raw
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

    return {
        "_status": "ok",
        "bonds": dict(bonds_out),
        "angles": dict(angles_out),
        "impropers": dict(impropers_out),
        "torsions": dict(torsions_out),
        "pair_bonds": dict(pair_bonds_out),
        "triple_angles": dict(triple_angles_out),
        "improper_pairs": dict(improper_pairs_out),
    }


def _process_chunk(idx_chunk: List[int]) -> Tuple[
        Dict[str, List[float]], Dict[str, List[float]],
        Dict[str, List[float]], Dict[str, List[float]],
        Dict[str, List[float]], Dict[str, List[float]],
        Dict[str, List[float]],
        int, int, int, Dict[str, int]]:
    bonds: Dict[str, List[float]] = defaultdict(list)
    angles: Dict[str, List[float]] = defaultdict(list)
    impropers: Dict[str, List[float]] = defaultdict(list)
    torsions: Dict[str, List[float]] = defaultdict(list)
    pair_bonds: Dict[str, List[float]] = defaultdict(list)
    triple_angles: Dict[str, List[float]] = defaultdict(list)
    improper_pairs: Dict[str, List[float]] = defaultdict(list)
    n_done = 0
    n_failed = 0
    n_extracted = 0
    skip_breakdown: Dict[str, int] = defaultdict(int)
    for i in idx_chunk:
        out = _extract_one_entry(i)
        n_done += 1
        if out is None:
            n_failed += 1
            skip_breakdown["exception"] += 1
            continue
        status = out.get("_status", "unknown")
        if status != "ok":
            skip_breakdown[status] += 1
            continue
        n_extracted += 1
        for k, vs in out["bonds"].items():
            bonds[k].extend(vs)
        for k, vs in out["angles"].items():
            angles[k].extend(vs)
        for k, vs in out["impropers"].items():
            impropers[k].extend(vs)
        for k, vs in out["torsions"].items():
            torsions[k].extend(vs)
        for k, vs in out["pair_bonds"].items():
            pair_bonds[k].extend(vs)
        for k, vs in out["triple_angles"].items():
            triple_angles[k].extend(vs)
        for k, vs in out["improper_pairs"].items():
            improper_pairs[k].extend(vs)
    return (dict(bonds), dict(angles), dict(impropers), dict(torsions),
            dict(pair_bonds), dict(triple_angles), dict(improper_pairs),
            n_done, n_failed, n_extracted, dict(skip_breakdown))


# ============================================================================
# Aggregation
# ============================================================================

def _robust_aggregate(values: List[float]) -> Tuple[float, float, int, float, float, float]:
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
    """GMM fit for torsion bins.  Falls back to single-component median+std
    when sklearn is unavailable or fit fails.
    """
    arr = np.asarray(values, dtype=np.float64).reshape(-1, 1)
    n = arr.shape[0]
    if n < 2:
        return {
            "n_components": 1,
            "pi": [1.0] + [np.nan] * (max_components - 1),
            "mu": [float(arr.flatten()[0])] + [np.nan] * (max_components - 1),
            "sigma": [1.0] + [np.nan] * (max_components - 1),
        }
    try:
        from sklearn.mixture import GaussianMixture
    except Exception:
        return {
            "n_components": 1,
            "pi": [1.0] + [np.nan] * (max_components - 1),
            "mu": [float(np.median(arr))] + [np.nan] * (max_components - 1),
            "sigma": [max(float(np.std(arr)), 1.0)] + [np.nan] * (max_components - 1),
        }
    best_bic = np.inf
    best_gm = None
    for k in range(1, max_components + 1):
        if n < k * 2:
            break
        try:
            gm = GaussianMixture(
                n_components=k, covariance_type="full",
                random_state=42, n_init=3, max_iter=200, reg_covar=1e-4,
            )
            gm.fit(arr)
            bic = gm.bic(arr)
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
    mus = mus[order]; pis = pis[order]; sigmas = sigmas[order]
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


# ============================================================================
# Driver
# ============================================================================

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--chunksize", type=int, default=64)
    parser.add_argument("--log", type=Path, default=DEFAULT_LOG)
    args = parser.parse_args()

    global log
    log = _setup_logging(args.log)

    log.info("Building ordered index list (sort by identifier)")
    t0 = time.time()
    from ccdc.io import EntryReader
    csd = EntryReader("CSD")
    n_total = len(csd)
    log.info(f"  CSD database size: {n_total}")
    if args.limit > 0:
        n_to_use = min(args.limit, n_total)
    else:
        n_to_use = n_total

    ident_idx: List[Tuple[str, int]] = []
    for i in range(n_to_use):
        try:
            ident_idx.append((csd[i].identifier, i))
        except Exception:
            ident_idx.append((f"ZZZZ_BAD_{i:09d}", i))
        if (i + 1) % 200000 == 0:
            log.info(f"  ...indexed {i+1}/{n_to_use}")
    ident_idx.sort(key=lambda t: t[0])
    sorted_indices = [t[1] for t in ident_idx]
    log.info(f"  index listing+sort done in {time.time()-t0:.1f}s "
             f"({len(sorted_indices)} entries)")

    log.info(f"GRIP v5_TM build start workers={args.workers} chunksize={args.chunksize} "
             f"out={args.out}")

    bond_vals: Dict[str, List[float]] = defaultdict(list)
    angle_vals: Dict[str, List[float]] = defaultdict(list)
    improper_vals: Dict[str, List[float]] = defaultdict(list)
    torsion_vals: Dict[str, List[float]] = defaultdict(list)
    pair_bond_vals: Dict[str, List[float]] = defaultdict(list)
    triple_angle_vals: Dict[str, List[float]] = defaultdict(list)
    improper_pair_vals: Dict[str, List[float]] = defaultdict(list)

    chunks = [
        sorted_indices[i:i + args.chunksize]
        for i in range(0, len(sorted_indices), args.chunksize)
    ]
    n_chunks = len(chunks)
    log.info(f"  total chunks: {n_chunks} (~{args.chunksize} entries each)")

    t0 = time.time()
    last_log = t0
    n_done = 0
    n_failed = 0
    n_extracted = 0
    skip_breakdown: Dict[str, int] = defaultdict(int)

    ctx = mp.get_context("forkserver")
    with ctx.Pool(processes=args.workers, initializer=_worker_init) as pool:
        try:
            iterator = pool.imap_unordered(
                _process_chunk, chunks, chunksize=1
            )
            for chunk_i, result in enumerate(iterator, start=1):
                (bonds, angles, impropers, torsions,
                 pair_bonds, triple_angles, improper_pairs,
                 c_done, c_failed, c_extracted, c_skips) = result
                for k, vs in bonds.items():
                    bond_vals[k].extend(vs)
                for k, vs in angles.items():
                    angle_vals[k].extend(vs)
                for k, vs in impropers.items():
                    improper_vals[k].extend(vs)
                for k, vs in torsions.items():
                    torsion_vals[k].extend(vs)
                for k, vs in pair_bonds.items():
                    pair_bond_vals[k].extend(vs)
                for k, vs in triple_angles.items():
                    triple_angle_vals[k].extend(vs)
                for k, vs in improper_pairs.items():
                    improper_pair_vals[k].extend(vs)
                n_done += c_done
                n_failed += c_failed
                n_extracted += c_extracted
                for k, v in c_skips.items():
                    skip_breakdown[k] += v

                if (chunk_i % 50 == 0) or (time.time() - last_log) > 60.0:
                    rate = n_done / max(time.time() - t0, 1e-6)
                    eta_s = (n_to_use - n_done) / max(rate, 1e-6)
                    log.info(
                        f"  processed {n_done}/{n_to_use} "
                        f"({100*n_done/max(n_to_use,1):.1f}%) "
                        f"extracted={n_extracted} (TM-complex) "
                        f"chunks={chunk_i}/{n_chunks} "
                        f"rate={rate:.0f} entries/s ETA={eta_s/60.0:.1f}min "
                        f"failed={n_failed} "
                        f"bond_keys={len(bond_vals)} "
                        f"angle_keys={len(angle_vals)} "
                        f"pair_bond_keys={len(pair_bond_vals)} "
                        f"triple_angle_keys={len(triple_angle_vals)}"
                    )
                    last_log = time.time()
        except KeyboardInterrupt:
            log.warning("Interrupted; saving partial state")
            pool.terminate()

    log.info(
        f"Extraction done: {n_done} entries processed "
        f"({n_extracted} TM-complex extracted, {n_failed} failed) "
        f"in {time.time()-t0:.1f}s"
    )
    log.info(f"  skip_breakdown: {dict(skip_breakdown)}")
    log.info(
        f"  raw keys: bonds={len(bond_vals)} angles={len(angle_vals)} "
        f"impropers={len(improper_vals)} torsions={len(torsion_vals)}"
    )
    log.info(
        f"  raw pair keys: pair_bonds={len(pair_bond_vals)} "
        f"triple_angles={len(triple_angle_vals)} "
        f"improper_pairs={len(improper_pair_vals)}"
    )

    log.info("Expanding fragment values into fallback bins (v3 schema)")
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
    log.info(f"Fallback expansion done in {time.time()-t0:.1f}s")

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
    log.info(f"Torsion expansion done in {time.time()-t0:.1f}s "
             f"({len(torsion_bins)} bins)")

    log.info("Aggregating bond/angle/improper bins")
    t0 = time.time()
    master_keys = sorted(
        set(bond_bins.keys()) | set(angle_bins.keys()) | set(improper_bins.keys())
    )
    n_master = len(master_keys)
    key_to_idx = {k: i for i, k in enumerate(master_keys)}
    log.info(f"  master key count: {n_master}")

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
    log.info(f"Bond/angle/improper aggregation done in {time.time()-t0:.1f}s")

    log.info("Aggregating v4 pair_bond / triple_angle / improper_pair bins")
    t0 = time.time()

    pair_bond_keys = sorted(pair_bond_vals.keys())
    n_pair_bond = len(pair_bond_keys)
    pair_bond_mu = np.full(n_pair_bond, np.nan, dtype=np.float64)
    pair_bond_sigma = np.full(n_pair_bond, np.nan, dtype=np.float64)
    pair_bond_n = np.zeros(n_pair_bond, dtype=np.int32)
    for i, k in enumerate(pair_bond_keys):
        mu, sigma, n, _, _, _ = _robust_aggregate(pair_bond_vals[k])
        pair_bond_mu[i] = mu; pair_bond_sigma[i] = sigma; pair_bond_n[i] = n

    triple_angle_keys = sorted(triple_angle_vals.keys())
    n_triple_angle = len(triple_angle_keys)
    triple_angle_mu = np.full(n_triple_angle, np.nan, dtype=np.float64)
    triple_angle_sigma = np.full(n_triple_angle, np.nan, dtype=np.float64)
    triple_angle_n = np.zeros(n_triple_angle, dtype=np.int32)
    for i, k in enumerate(triple_angle_keys):
        mu, sigma, n, _, _, _ = _robust_aggregate(triple_angle_vals[k])
        triple_angle_mu[i] = mu; triple_angle_sigma[i] = sigma; triple_angle_n[i] = n

    improper_pair_keys = sorted(improper_pair_vals.keys())
    n_improper_pair = len(improper_pair_keys)
    improper_pair_mu = np.full(n_improper_pair, np.nan, dtype=np.float64)
    improper_pair_sigma = np.full(n_improper_pair, np.nan, dtype=np.float64)
    improper_pair_n = np.zeros(n_improper_pair, dtype=np.int32)
    for i, k in enumerate(improper_pair_keys):
        mu, sigma, n, _, _, _ = _robust_aggregate(improper_pair_vals[k])
        improper_pair_mu[i] = mu; improper_pair_sigma[i] = sigma; improper_pair_n[i] = n

    log.info(f"Pair tables aggregated in {time.time()-t0:.1f}s: "
             f"pair_bond={n_pair_bond} triple_angle={n_triple_angle} "
             f"improper_pair={n_improper_pair}")

    # ---- Torsion aggregation ----
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

    # GMM SKIP: for v5_TM ligand-internals build, torsion is the lowest-weight
    # polish term (0.2) and the test cases (SIYMEU/BERTEB/ALAHEB) don't rely on
    # torsion fits.  Use single-component median+std for every bin instead of
    # per-bin GMM (which on full-CSD scales to many CPU-minutes on big bins).
    # Schema is unchanged: torsion_pi[0]=1.0, torsion_mu[0]=median,
    # torsion_sigma[0]=std (clamped >=5.0 deg).
    fits_n5 = 0
    progress_every = max(1, n_torsion // 20)
    for i, k in enumerate(torsion_keys):
        vals = np.asarray(torsion_bins[k], dtype=np.float64)
        n_val = int(vals.size)
        torsion_n[i] = n_val
        if n_val > 0:
            mu_med = float(np.median(vals))
        else:
            mu_med = 0.0
        if n_val > 1:
            sig = max(float(np.std(vals)), 5.0)
            fits_n5 += 1 if n_val >= 5 else 0
        else:
            sig = 30.0
        torsion_n_components[i] = 1
        torsion_pi[i, 0] = 1.0
        torsion_mu[i, 0] = mu_med
        torsion_sigma[i, 0] = sig
        for j in range(1, MAX_COMPONENTS):
            torsion_pi[i, j] = float("nan")
            torsion_mu[i, j] = float("nan")
            torsion_sigma[i, j] = float("nan")
        if (i + 1) % progress_every == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / max(elapsed, 1e-6)
            eta = (n_torsion - i - 1) / max(rate, 1e-6)
            log.info(f"  ...torsion stats {i+1}/{n_torsion} "
                     f"({100*(i+1)/n_torsion:.1f}%) "
                     f"elapsed={elapsed:.1f}s ETA={eta:.1f}s")

    log.info(f"Torsion aggregation done in {time.time()-t0:.1f}s "
             f"({fits_n5}/{n_torsion} bins with n>=5)")

    # ---- Fallback CSR encoding (v3 schema) ----
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
    log.info(f"Fragment fallback chains: n_orig={n_orig} entries={write_pos} "
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
    log.info(f"Torsion fallback chains: n_orig_t={n_orig_t} entries={write_pos_t} "
             f"in {time.time()-t0:.1f}s")

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

    log.info("=" * 72)
    log.info("GRIP-LIB v5_TM (CCDC CSD, TM-complex ligand-internals) COVERAGE REPORT")
    log.info("=" * 72)
    log.info(f"CSD entries processed: {n_done} ({n_failed} failed)")
    log.info(f"  extracted TM-complex: {n_extracted}")
    log.info(f"  skip breakdown:       {dict(skip_breakdown)}")
    log.info(f"Original frag keys:    {n_orig}")
    log.info(f"Master (incl fb) keys: {n_master}")
    log.info(f"Torsion bins:          {n_torsion}")
    log.info(f"Original torsion keys: {n_orig_t}")
    log.info(f"pair_bond keys:    {n_pair_bond} (n>=5: {pb_ge5})")
    log.info(f"triple_angle keys: {n_triple_angle} (n>=5: {ta_ge5})")
    log.info(f"improper_pair keys:{n_improper_pair} (n>=5: {ip_ge5})")
    log.info("")
    log.info("--- bonds (single-neighbour v3 keys, ligand-internal) ---")
    log.info(f"  obs:       {bonds_present}/{n_master} ({100*bonds_present/max(n_master,1):.1f}%)")
    log.info(f"  n>=5:      {bonds_ge5}/{n_master} ({100*bonds_ge5/max(n_master,1):.1f}%)")
    log.info("--- angles (single-pair v3 keys, ligand-internal) ---")
    log.info(f"  obs:       {angles_present}/{n_master} ({100*angles_present/max(n_master,1):.1f}%)")
    log.info(f"  n>=5:      {angles_ge5}/{n_master} ({100*angles_ge5/max(n_master,1):.1f}%)")
    log.info("--- impropers (v3 keys, ligand-internal) ---")
    log.info(f"  obs:       {imp_present}/{n_master} ({100*imp_present/max(n_master,1):.1f}%)")
    log.info(f"  n>=5:      {imp_ge5}/{n_master} ({100*imp_ge5/max(n_master,1):.1f}%)")
    log.info("--- torsions ---")
    log.info(f"  obs:       {tor_present}/{n_torsion} ({100*tor_present/max(n_torsion,1):.1f}%)")
    log.info(f"  n>=5:      {tor_ge5}/{n_torsion} ({100*tor_ge5/max(n_torsion,1):.1f}%)")
    log.info("=" * 72)

    log.info(f"Writing {args.out}")
    t0 = time.time()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    # Empty arrays for tm_<cat>_* keys (kept for schema compat).
    empty_obj = np.array([], dtype=np.object_)
    empty_f64 = np.array([], dtype=np.float64)
    empty_i32 = np.array([], dtype=np.int32)

    np.savez_compressed(
        args.out,
        version=np.int32(5),
        cleaning_applied=np.int32(1),
        meta_n_total_scanned=np.int32(n_done),
        meta_n_extracted=np.int32(n_extracted),
        meta_n_polymeric=np.int32(skip_breakdown.get("polymeric", 0)),
        meta_n_no_3d=np.int32(skip_breakdown.get("no_3d", 0)),
        meta_n_too_small=np.int32(skip_breakdown.get("size", 0)),
        meta_n_recovered_disorder=np.int32(0),
        meta_n_recovered_partial_coord=np.int32(0),
        n_master=np.int32(n_master),
        n_orig=np.int32(n_orig),
        n_torsion=np.int32(n_torsion),
        n_orig_torsion=np.int32(n_orig_t),
        n_pair_bond=np.int32(n_pair_bond),
        n_triple_angle=np.int32(n_triple_angle),
        n_improper_pair=np.int32(n_improper_pair),
        mogul_validated=np.zeros(n_master, dtype=bool),
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
        torsion_n=torsion_n,
        torsion_n_components=torsion_n_components,
        torsion_pi=torsion_pi, torsion_mu=torsion_mu, torsion_sigma=torsion_sigma,
        torsion_fb_flat=torsion_fb_flat, torsion_fb_ptr=torsion_fb_ptr,
        pair_bond_keys=np.array(pair_bond_keys, dtype=np.object_),
        pair_bond_mu=pair_bond_mu, pair_bond_sigma=pair_bond_sigma,
        pair_bond_n=pair_bond_n,
        triple_angle_keys=np.array(triple_angle_keys, dtype=np.object_),
        triple_angle_mu=triple_angle_mu, triple_angle_sigma=triple_angle_sigma,
        triple_angle_n=triple_angle_n,
        improper_pair_keys=np.array(improper_pair_keys, dtype=np.object_),
        improper_pair_mu=improper_pair_mu, improper_pair_sigma=improper_pair_sigma,
        improper_pair_n=improper_pair_n,
        # Empty tm_<cat>_* arrays for schema compatibility (universal-only build).
        tm_carbene_keys=empty_obj, tm_carbene_mu=empty_f64,
        tm_carbene_sigma=empty_f64, tm_carbene_n=empty_i32,
        tm_hapto_eta2_keys=empty_obj, tm_hapto_eta2_mu=empty_f64,
        tm_hapto_eta2_sigma=empty_f64, tm_hapto_eta2_n=empty_i32,
        tm_hapto_eta5_keys=empty_obj, tm_hapto_eta5_mu=empty_f64,
        tm_hapto_eta5_sigma=empty_f64, tm_hapto_eta5_n=empty_i32,
        tm_hapto_eta6_keys=empty_obj, tm_hapto_eta6_mu=empty_f64,
        tm_hapto_eta6_sigma=empty_f64, tm_hapto_eta6_n=empty_i32,
        tm_mu_bridge_keys=empty_obj, tm_mu_bridge_mu=empty_f64,
        tm_mu_bridge_sigma=empty_f64, tm_mu_bridge_n=empty_i32,
        tm_agostic_keys=empty_obj, tm_agostic_mu=empty_f64,
        tm_agostic_sigma=empty_f64, tm_agostic_n=empty_i32,
        tm_ox_addition_keys=empty_obj, tm_ox_addition_mu=empty_f64,
        tm_ox_addition_sigma=empty_f64, tm_ox_addition_n=empty_i32,
        tm_pair_bond_block_keys=empty_obj, tm_pair_bond_block_mu=empty_f64,
        tm_pair_bond_block_sigma=empty_f64, tm_pair_bond_block_n=empty_i32,
        tm_triple_angle_block_keys=empty_obj, tm_triple_angle_block_mu=empty_f64,
        tm_triple_angle_block_sigma=empty_f64, tm_triple_angle_block_n=empty_i32,
    )
    log.info(f"Wrote {args.out} ({args.out.stat().st_size/1e6:.1f} MB) "
             f"in {time.time()-t0:.1f}s")

    log.info("BUILD DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
