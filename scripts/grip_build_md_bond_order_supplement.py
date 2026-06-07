#!/usr/bin/env python3
"""GRIP M-D bond-order supplement library.

Audit verdict (2026-06-07): the v5_TM library has *excellent* W/Mo/Tc/Re donor
coverage (n>>1000 per donor element) and returns realistic empirical means for
single-bond cases (W-N(sp3)=2.224, W-O(sp3)=1.926, W-Cl=2.411).

THE problem: the ``pair_bond_keys`` table key is ``(Z1, hyb1, Z2, hyb2)``
which carries NO bond-order context.  This collapses chemically distinct
populations into a single mean:

  * W-N(sp2):  M=N imido/nitride (~1.75 Å, BO≈2-3) pooled with
               κ¹-pyridine N (~2.20 Å, BO=1, aromatic ring N).
               Library returns μ=1.779 (n=2962) — wrong for pyridine.
  * Re/Tc/Mo analogues: same problem.
  * W-C(sp):   M=C carbenoid + W-CO carbonyl + W-C alkynyl pooled.
  * W-O(sp2):  M=O oxo (~1.7 Å) pooled with κ¹-carboxylate-O.

Universal fix: extend the M-D pair_bond key to include bond-order context.
NEW key format (M-D only):
    ``[Z1, hyb1, "@bo<N>", Z2, hyb2]``
where ``<N>`` is ``1``, ``2``, ``3``, ``4`` (CCDC native bond_type → integer)
or ``a`` (aromatic).  All non-M-D bonds keep the legacy 4-tuple key — this
is an additive supplement.

This script ONLY emits M-D pair_bond keys (Z1 or Z2 is a metal); the rest of
the library is consumed unchanged.  Output: a supplement npz with a single
new array set ``md_bo_pair_*`` that the loader can merge into the existing
library at runtime.

Scope: every CSD entry that contains at least one of {W, Mo, Tc, Re} -- this
is the user's mission focus for 5d-mid donors.  About 75-100k entries (vs
1.4M total), so the extraction completes in 10-30 min vs 6-12 h for a full
re-build.

USAGE
-----
::

    PYTHONHASHSEED=0 \\
    /home/qmchem_max/CCDC/CSD_2026.1/ccdc-software/csd-python-api/run_csd_python_api \\
        scripts/grip_build_md_bond_order_supplement.py \\
        --workers 16 \\
        --metals W,Mo,Tc,Re \\
        --out /home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_md_bo_supplement_5d_mid.npz

OUTPUT SCHEMA
-------------
::

    version : int = 1
    metals  : array[str] -- which metals were included
    md_bo_pair_keys : object array of JSON strings
                      "[Z1,hyb1,'@bo<X>',Z2,hyb2]" canonical-sorted on (Z, hyb)
                      with the metal endpoint always second (Z2 is the metal).
    md_bo_pair_mu, md_bo_pair_sigma, md_bo_pair_p5, md_bo_pair_p50,
    md_bo_pair_p95, md_bo_pair_n  : per-key aggregates (robust median/MAD).

Author: hmaximilian <hmaximilian496@gmail.com>
"""
from __future__ import annotations

import argparse
import json
import logging
import multiprocessing as mp
import os
import random
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# Determinism BEFORE any other import
os.environ.setdefault("PYTHONHASHSEED", "0")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
np.random.seed(42)
random.seed(42)


_METALS_D = frozenset({
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
})


def _infer_hyb(symbol: str, total_degree: int) -> str:
    """Mirror grip_build_mogul_lib_v5_TM_ligand_internals.py _infer_hyb."""
    if symbol in _METALS_D:
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
        return "sp3"
    if symbol == "O":
        if deg <= 1:
            return "sp2"
        return "sp3"
    if symbol == "S":
        if deg <= 1:
            return "sp2"
        return "sp3"
    if symbol == "P":
        if deg <= 4:
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
    return "sp3"


def _bo_tag(bond_type, is_aromatic: bool) -> str:
    """Convert CCDC bond_type to a stable tag character.

    Returns one of: ``'a'`` (aromatic), ``'1'``, ``'2'``, ``'3'``, ``'4'``,
    ``'?'`` (unknown).
    """
    if is_aromatic:
        return "a"
    bt = str(bond_type).lower()
    if "aromatic" in bt:
        return "a"
    if bt == "single" or bt == "1":
        return "1"
    if bt == "double" or bt == "2":
        return "2"
    if bt == "triple" or bt == "3":
        return "3"
    if bt == "quadruple" or bt == "4":
        return "4"
    if "delocalized" in bt or "pi" in bt:
        return "a"
    return "?"


def _md_bo_key(z_metal: str, z_donor: str, hyb_donor: str, bo_tag: str) -> str:
    """Canonical M-D key with bond-order context.

    Format: ``[Z_donor, hyb_donor, '@bo<X>', Z_metal, '*']``
    Metal always second; metal hyb always '*' (transition metals have no
    Hill-style hyb).  Donor element + hybridization keep their natural form.
    """
    return json.dumps(
        [str(z_donor), str(hyb_donor), f"@bo{bo_tag}", str(z_metal), "*"],
        separators=(",", ":"),
        ensure_ascii=False,
    )


# ---------------------------------------------------------------------------
# CCDC worker
# ---------------------------------------------------------------------------
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


def _extract_md_bo_for_entry(args):
    """Return list of (key_str, distance) for ALL metal-donor bonds in entry.

    Quality gates (mirror v5_TM build):
        * Entry must have 3D coordinates.
        * Skip polymeric/disordered entries (use main asymmetric component).
        * Bond distance must be in (0.5, 3.5) Å -- harvested CCDC's own bonds,
          so contact-vs-bond decision is the CCDC bond model's.

    The returned list is partial -- aggregation happens in the parent process.
    """
    identifier, target_metals = args
    out: List[Tuple[str, float]] = []
    try:
        r = _ensure_reader()
        e = r.entry(identifier)
        if e is None:
            return out
        if not e.has_3d_structure:
            return out
        mol = e.molecule
        if mol is None or len(mol.atoms) < 2:
            return out
        # Drop disorder: take heaviest component when polymeric
        if mol.is_polymeric:
            return out
        # Build atom degree counts from CCDC bonds (more reliable than RDKit
        # for organometallics; matches CSD's native graph).  Use the CCDC
        # ``atom.index`` attribute (stable across atom-proxy fetches) -- using
        # ``id(a)`` is unsafe because CCDC returns a new proxy per access.
        atoms = mol.atoms
        n_atoms = len(atoms)
        degrees = [0] * n_atoms
        for b in mol.bonds:
            try:
                a1, a2 = b.atoms
                i1 = int(a1.index)
                i2 = int(a2.index)
            except Exception:
                continue
            if 0 <= i1 < n_atoms and 0 <= i2 < n_atoms:
                degrees[i1] += 1
                degrees[i2] += 1
        # Iterate every bond; emit entry only if one endpoint is a target metal
        # and the other is a non-metal donor.
        for b in mol.bonds:
            try:
                a1, a2 = b.atoms
            except Exception:
                continue
            s1 = a1.atomic_symbol
            s2 = a2.atomic_symbol
            m_a1 = s1 in target_metals
            m_a2 = s2 in target_metals
            if m_a1 == m_a2:
                continue  # both metal or both non-metal -- not M-D
            # Filter out other-metal endpoint that isn't a target metal
            if m_a1:
                metal_sym = s1
                metal_atom = a1
                donor_sym = s2
                donor_atom = a2
            else:
                metal_sym = s2
                metal_atom = a2
                donor_sym = s1
                donor_atom = a1
            # Skip M-M bonds (other-metal donor) -- those are not "donor" bonds
            if donor_sym in _METALS_D:
                continue
            # Coordinates → distance
            try:
                c_m = metal_atom.coordinates
                c_d = donor_atom.coordinates
                if c_m is None or c_d is None:
                    continue
                dist = float(np.linalg.norm([
                    c_m.x - c_d.x, c_m.y - c_d.y, c_m.z - c_d.z
                ]))
            except Exception:
                continue
            if not (0.5 < dist < 3.5):
                continue
            # Donor degree from CCDC graph (drops M-edge influence: include
            # all incident bonds -- the lookup-time donor hyb is from the
            # full graph too, so this is consistent).
            try:
                idx_d = int(donor_atom.index)
                deg_d = degrees[idx_d] if 0 <= idx_d < n_atoms else 1
            except Exception:
                deg_d = 1
            hyb_d = _infer_hyb(donor_sym, deg_d)
            # Bond-order tag from CCDC native model.
            try:
                is_arom = bool(getattr(b, "is_conjugated", False)) or \
                          bool(getattr(b, "is_in_aromatic_ring", False))
            except Exception:
                is_arom = False
            bo = _bo_tag(b.bond_type, is_arom)
            key = _md_bo_key(metal_sym, donor_sym, hyb_d, bo)
            out.append((key, dist))
    except Exception:
        return out
    return out


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------
def _robust_aggregate(values: List[float]) -> Tuple[float, float, float, float, float, int]:
    """Return (mu, sigma, p5, p50, p95, n) via median + 1.4826*MAD."""
    arr = np.asarray(values, dtype=np.float64)
    n = int(len(arr))
    if n == 0:
        return (float("nan"), float("nan"), float("nan"), float("nan"), float("nan"), 0)
    med = float(np.median(arr))
    mad = float(np.median(np.abs(arr - med)))
    sigma = 1.4826 * mad
    p5 = float(np.percentile(arr, 5))
    p95 = float(np.percentile(arr, 95))
    return (med, sigma, p5, med, p95, n)


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=str, required=True,
                        help="Output npz path")
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--metals", type=str, default="W,Mo,Tc,Re",
                        help="Comma-separated metal symbols to include")
    parser.add_argument("--limit", type=int, default=None,
                        help="Process only the first N candidate identifiers (testing)")
    parser.add_argument("--log", type=str, default="/tmp/grip_md_bo_supplement_build.log")
    args = parser.parse_args()

    target_metals = frozenset(s.strip() for s in args.metals.split(",") if s.strip())
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[
            logging.FileHandler(args.log, mode="w"),
            logging.StreamHandler(sys.stdout),
        ],
        force=True,
    )
    log = logging.getLogger("grip_md_bo_supplement")
    log.info(f"Target metals: {sorted(target_metals)}")
    log.info(f"Output: {out_path}")

    # First pass: collect identifiers containing target metals.  Filter via
    # the CSD substructure query for speed when possible; else stream-iterate.
    log.info("Pass 1: collecting candidate identifiers (substructure pre-filter).")
    t0 = time.time()
    candidates: List[str] = []
    try:
        from ccdc.search import SubstructureSearch, SMARTSSubstructure
        # SMARTS-OR for each target metal
        for m in sorted(target_metals):
            sm = SubstructureSearch()
            sm.add_substructure(SMARTSSubstructure(f"[{m}]"))
            log.info(f"  Searching SMARTS [{m}] ...")
            hits = sm.search()
            for h in hits:
                candidates.append(h.entry.identifier)
            log.info(f"    [{m}]: {len(hits)} hits (running total {len(candidates)})")
    except Exception as ex:
        log.warning(f"Substructure pre-filter failed ({ex!r}); stream-iterating full CSD.")
        from ccdc.io import EntryReader
        r = EntryReader("CSD")
        for i, e in enumerate(r):
            if i % 50000 == 0:
                log.info(f"  stream {i}/{len(r)} (candidates so far: {len(candidates)})")
            try:
                mol = e.molecule
                if mol is None:
                    continue
                syms = set(a.atomic_symbol for a in mol.atoms)
                if syms & target_metals:
                    candidates.append(e.identifier)
            except Exception:
                continue

    # Lex-sort + dedup for determinism
    candidates = sorted(set(candidates))
    log.info(f"Pass 1 done in {time.time()-t0:.1f}s: {len(candidates)} candidate identifiers.")
    if args.limit:
        candidates = candidates[: args.limit]
        log.info(f"--limit applied: processing first {len(candidates)} only.")

    # Pass 2: extract M-D bonds in parallel
    log.info(f"Pass 2: extracting M-D bonds with {args.workers} workers.")
    t0 = time.time()
    bonds_by_key: Dict[str, List[float]] = defaultdict(list)
    args_iter = ((ident, target_metals) for ident in candidates)
    if args.workers > 1:
        with mp.Pool(args.workers, initializer=_worker_init) as pool:
            chunksize = max(1, min(64, len(candidates) // (args.workers * 8)))
            for n_done, hits in enumerate(
                pool.imap_unordered(_extract_md_bo_for_entry, args_iter, chunksize=chunksize),
                start=1,
            ):
                for key, dist in hits:
                    bonds_by_key[key].append(dist)
                if n_done % 2000 == 0:
                    log.info(f"  {n_done}/{len(candidates)} entries, "
                             f"keys-so-far={len(bonds_by_key)}, "
                             f"elapsed={time.time()-t0:.0f}s")
    else:
        _worker_init()
        for n_done, a in enumerate(args_iter, start=1):
            for key, dist in _extract_md_bo_for_entry(a):
                bonds_by_key[key].append(dist)
            if n_done % 500 == 0:
                log.info(f"  serial {n_done}/{len(candidates)}, keys-so-far={len(bonds_by_key)}")
    log.info(f"Pass 2 done in {time.time()-t0:.1f}s: {len(bonds_by_key)} unique keys.")

    # Aggregate
    log.info("Aggregating per-key distributions.")
    keys_sorted = sorted(bonds_by_key.keys())
    n_keys = len(keys_sorted)
    mu = np.zeros(n_keys, dtype=np.float64)
    sigma = np.zeros(n_keys, dtype=np.float64)
    p5 = np.zeros(n_keys, dtype=np.float64)
    p50 = np.zeros(n_keys, dtype=np.float64)
    p95 = np.zeros(n_keys, dtype=np.float64)
    n_arr = np.zeros(n_keys, dtype=np.int32)
    for i, k in enumerate(keys_sorted):
        m, s, q5, q50, q95, n = _robust_aggregate(bonds_by_key[k])
        mu[i] = m
        sigma[i] = s
        p5[i] = q5
        p50[i] = q50
        p95[i] = q95
        n_arr[i] = n

    # Some stats
    nh = int((n_arr >= 5).sum())
    log.info(f"Keys with n>=5: {nh}/{n_keys}")

    # Write
    log.info(f"Writing npz: {out_path}")
    np.savez(
        out_path,
        version=np.int32(1),
        metals=np.array(sorted(target_metals), dtype=np.object_),
        md_bo_pair_keys=np.array(keys_sorted, dtype=np.object_),
        md_bo_pair_mu=mu,
        md_bo_pair_sigma=sigma,
        md_bo_pair_p5=p5,
        md_bo_pair_p50=p50,
        md_bo_pair_p95=p95,
        md_bo_pair_n=n_arr,
    )
    log.info(f"Done.  Output: {out_path} ({out_path.stat().st_size/1e6:.2f} MB)")


if __name__ == "__main__":
    main()
