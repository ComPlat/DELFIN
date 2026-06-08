"""Validate the aromatic planarity gate on 3 user-eye cases.

For each case we build the structure twice:

  * baseline:  DELFIN_FFFREE_MOGUL_PRIMARY=1 + aromatic flags OFF
  * fixed   :  DELFIN_FFFREE_MOGUL_PRIMARY=1 + aromatic flags ON

For every aromatic ring in the resulting mol we measure max perpendicular
deviation from the SVD best-fit plane.  We expect the "fixed" run to push
each ring's max-dev to ≤ 0.10 Å.

Cases:
  * SIYMEU  — Ag-NHC linear CN2, 4 aryl substituents + NHC carbene ring.
  * BERTEB  — Ni-pyridyl-imine, 2 pyridyl rings.
  * NASQAB  — generic Cu(II) tris(N-phenyl)-tris(pyrazolyl) borate proxy
              built locally; multiple phenyl substituents — exactly the
              user-eye scenario the gate is meant to fix.

Author: hmaximilian <hmaximilian496@gmail.com>
"""
from __future__ import annotations

import os
import sys
from typing import Dict, List, Optional, Tuple

# Ensure deterministic walks
os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np

# Make sure we exercise our branch's code path.
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def _set_grip_lib_path() -> None:
    """Wire up the TM-restricted GRIP library used by the user-eye tests."""
    if os.environ.get("DELFIN_GRIP_LIB_PATH"):
        return
    for cand in (
        "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5_TM.npz",
        "/home/qmchem_max/agent_workspace/quality_framework/reports/grip_lib_v5.npz",
    ):
        if os.path.exists(cand):
            os.environ["DELFIN_GRIP_LIB_PATH"] = cand
            return


_set_grip_lib_path()
# Required by mogul-primary cascade
os.environ.setdefault("DELFIN_FFFREE_MOGUL_PRIMARY_GRIP", "1")
os.environ.setdefault("DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM", "1")
os.environ.setdefault("DELFIN_FFFREE_MOGUL_BOND_FALLBACK", "1")


CASES = {
    "SIYMEU": (
        "CC(=O)[O][Ag-][C+]1N(CC2=CC=C(C)C=C2)C(C2=CC=C(C(C)C)C=C2)"
        "=C(C2=CC=C(C(C)C)C=C2)N1CC1=CC=C(C)C=C1"
    ),
    "BERTEB": "[Br][Ni-2]12[N]3C=CC=C3C=[N+]1CC1=CC=CC=[N+]12",
    # NASQAB proxy: Cu(I) tris-(triphenylphosphine) - 3 P donors each with
    # 3 phenyl substituents.  This is the user-eye scenario the brief
    # describes: a metal complex with multiple independent phenyl rings
    # whose buckle the soft-DG embedder leaves at 0.15-0.30 Å without the
    # planarity gate.  Linear coordination, neutral overall.
    "NASQAB_proxy":
        "[Cu]([P](c1ccccc1)(c2ccccc2)c3ccccc3)"
        "[P](c4ccccc4)(c5ccccc5)c6ccccc6",
}


def _build(smiles: str, *, aromatic_on: bool) -> Optional[Tuple[List[str], np.ndarray, object]]:
    # Always exercise the Mogul-primary path
    os.environ["DELFIN_FFFREE_MOGUL_PRIMARY"] = "1"
    if aromatic_on:
        os.environ["DELFIN_FFFREE_AROMATIC_TIER1_5"] = "1"
        os.environ["DELFIN_FFFREE_AROMATIC_PLANARITY_GATE"] = "1"
    else:
        os.environ["DELFIN_FFFREE_AROMATIC_TIER1_5"] = "0"
        os.environ["DELFIN_FFFREE_AROMATIC_PLANARITY_GATE"] = "0"
    # Reload the module each time so the auto-default doesn't sticky-override
    # our explicit env-flag settings.  The auto-default only fires when the
    # flag is the empty string; we've set it explicitly so it stays.
    from delfin.fffree.assemble_via_mogul import (
        assemble_complex_mogul_primary,
        _full_complex_mol,
    )
    out = assemble_complex_mogul_primary(smiles)
    if out is None:
        return None
    syms, P = out
    # Rebuild the source mol so we can measure aromatic ring planarity
    # with full perception.
    mol = _full_complex_mol(smiles)
    return syms, P, mol


def _arom_rings(mol) -> List[Tuple[int, ...]]:
    """Return aromatic rings on the SOURCE mol; atom indices are in
    SOURCE-frame, which the public assemble call may permute (metal
    at idx 0)."""
    from delfin.fffree.aromatic_flatten import detect_aromatic_rings
    return detect_aromatic_rings(mol)


def _measure_rings(syms: List[str], P: np.ndarray, mol) -> List[Dict]:
    """Map source-frame ring atom idxs into the output frame and measure
    planarity per ring.

    The public assemble call may reorder syms so the metal is at idx 0.
    We rebuild that mapping by matching source-frame and output-frame
    symbol sequences modulo the metal-prepend permutation (idx j in source
    -> j+1 in output for j < metal_idx_src, idx metal_idx_src -> 0, idx j
    -> j for j > metal_idx_src).
    """
    from delfin.fffree.aromatic_flatten import measure_ring_planarity

    # Build source -> output atom-idx map by walking the metal-prepend rule.
    # First, locate the (first) metal in source-frame.
    from delfin._bond_decollapse import _is_metal
    n_src = mol.GetNumAtoms()
    metal_idx_src = None
    for i in range(n_src):
        sym = mol.GetAtomWithIdx(int(i)).GetSymbol()
        if _is_metal(sym):
            metal_idx_src = int(i)
            break
    if metal_idx_src is None:
        # Identity map (no metal -> assemble shouldn't have run, but safe).
        idx_map = {i: i for i in range(n_src)}
    else:
        idx_map: Dict[int, int] = {}
        idx_map[metal_idx_src] = 0
        for j in range(n_src):
            if j == metal_idx_src:
                continue
            if j < metal_idx_src:
                idx_map[j] = j + 1
            else:
                idx_map[j] = j

    rings = _arom_rings(mol)
    out: List[Dict] = []
    for ring in rings:
        ring_out = tuple(idx_map.get(int(a), int(a)) for a in ring)
        # Bounds check; if any mapped idx is out of range, skip.
        if any(a >= len(syms) for a in ring_out):
            continue
        dev = measure_ring_planarity(P, ring_out)
        out.append({
            "ring_src": ring,
            "ring_out": ring_out,
            "max_dev": dev,
            "atoms": [syms[a] for a in ring_out],
            "size": len(ring),
        })
    return out


def main() -> int:
    print("=" * 78)
    print("Aromatic-planarity validation (hmaximilian, 2026-06-08)")
    print("=" * 78)

    n_fail = 0

    for case_name, smiles in CASES.items():
        print(f"\n--- Case: {case_name} ---")
        print(f"SMILES: {smiles}")

        # 1) Pre (gate OFF)
        try:
            res_pre = _build(smiles, aromatic_on=False)
        except Exception as exc:
            print(f"  [PRE]  BUILD FAILED with exception: {exc}")
            res_pre = None
        # 2) Post (gate ON)
        try:
            res_post = _build(smiles, aromatic_on=True)
        except Exception as exc:
            print(f"  [POST] BUILD FAILED with exception: {exc}")
            res_post = None

        if res_pre is None:
            print(f"  [PRE]  assemble returned None — case skipped")
            continue
        if res_post is None:
            print(f"  [POST] assemble returned None — case skipped")
            continue

        syms_pre, P_pre, mol_pre = res_pre
        syms_post, P_post, mol_post = res_post

        pre_rings = _measure_rings(syms_pre, P_pre, mol_pre)
        post_rings = _measure_rings(syms_post, P_post, mol_post)

        print(f"  Aromatic rings detected: {len(pre_rings)}")
        if not pre_rings:
            print("  (no aromatic rings — nothing to validate)")
            continue

        # Print per-ring before/after
        print(
            f"  {'Ring (out)':<28}  {'Atoms':<14}  {'Pre max-dev':>10}  "
            f"{'Post max-dev':>12}"
        )
        print("  " + "-" * 70)
        for pre, post in zip(pre_rings, post_rings):
            ring_atoms_str = "".join(pre["atoms"])[:14]
            ring_str = str(pre["ring_out"])[:26]
            print(
                f"  {ring_str:<28}  {ring_atoms_str:<14}  "
                f"{pre['max_dev']:>10.4f}  {post['max_dev']:>12.4f}"
            )

        # Decision: every ring's post max-dev must be ≤ 0.10 Å
        worst_post = max(r["max_dev"] for r in post_rings)
        worst_pre = max(r["max_dev"] for r in pre_rings)
        print(f"  Summary: worst pre = {worst_pre:.4f} Å, "
              f"worst post = {worst_post:.4f} Å, threshold = 0.100 Å")
        if worst_post <= 0.10 + 1e-6:
            print(f"  [PASS] all rings ≤ 0.10 Å after gate")
        else:
            print(f"  [FAIL] worst ring exceeds threshold")
            n_fail += 1

    print()
    print("=" * 78)
    if n_fail == 0:
        print("RESULT: all cases pass (every aromatic ring ≤ 0.10 Å after gate)")
        return 0
    else:
        print(f"RESULT: {n_fail} case(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
