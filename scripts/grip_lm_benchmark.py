"""GRIP L-BFGS vs LM (TRF) — 50-SMILES smoke benchmark.

Usage::

    PYTHONHASHSEED=0 micromamba run -n delfin python scripts/grip_lm_benchmark.py

Reports per-method:
* Convergence count (returned non-trivially polished)
* Wall time (sum / mean)
* Final mogul severity
* Final position RMSD between methods
* Function evaluation count (LM nfev; L-BFGS nit)

Deterministic when ``PYTHONHASHSEED=0`` is set (the benchmark itself uses
RDKit's ``randomSeed=42`` for embedding so the input geometries are pinned).
"""
from __future__ import annotations

import os

os.environ.setdefault("PYTHONHASHSEED", "0")

import statistics
import time
from typing import List, Tuple

import numpy as np


# A representative 50-SMILES sample covering organics, simple chelates,
# aromatics, and small molecules.  Goal is a sanity check, not a coverage
# benchmark; the production benchmark is a 500-SMILES pool smoke against
# the full assemble_complex pipeline.
SMILES_50: List[str] = [
    "CCO", "CCN", "CCC", "c1ccccc1", "Cc1ccccc1", "Oc1ccccc1", "Nc1ccccc1",
    "CC(=O)O", "CCOC", "CCOCC", "CC(=O)C", "CC#N", "CN(C)C", "CCNCC",
    "CC(N)C", "OCC(O)CO", "OC1CCCCC1", "c1ccc2ccccc2c1", "c1ccncc1",
    "c1ccsc1", "c1ccoc1", "C1CCNCC1", "C1CCOCC1", "C1CC1", "C1CCCCC1",
    "C1CCC1", "C1=CC=CC=C1", "CC(C)C", "CC(C)O", "CC(C)N",
    "CCCCO", "CCCCN", "CCCCCC", "c1ccc(O)cc1", "c1ccc(N)cc1",
    "C(=O)O", "C(=O)N", "CCC(=O)CC", "CCCC=O", "CC=CC",
    "OCC=CCO", "NCCN", "OCCO", "NCC(=O)O", "SCC=O",
    "Cc1ccc(C)cc1", "Cc1ccncc1", "CC1=CC=CC=C1O",
    "CCOC(=O)C", "C[N+](C)(C)C",
]


def _embed(smiles: str):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass
    return mol


def _coords(mol) -> np.ndarray:
    conf = mol.GetConformer()
    return np.array(
        [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
        dtype=np.float64,
    )


def _safe_donor(mol) -> int:
    """Pick a deterministic 'donor' atom index for the M-D constraint."""
    if mol.GetNumAtoms() < 2:
        return 0
    return 1


def _run_one_smiles(smiles: str, method: str) -> Tuple[bool, float, float, int, float]:
    """Returns (success, wall_time_s, final_sev, n_iter, rmsd_to_input)."""
    from delfin.fffree.grip_polish import grip_polish, _GRIP_METHOD_ENV
    if method == "lm":
        os.environ[_GRIP_METHOD_ENV] = "lm"
    else:
        os.environ.pop(_GRIP_METHOD_ENV, None)

    mol = _embed(smiles)
    if mol is None:
        return False, 0.0, float("nan"), 0, 0.0
    P0 = _coords(mol)
    # Slight perturbation so the polish has work to do.
    rng = np.random.default_rng(7)
    P_in = P0 + rng.standard_normal(P0.shape) * 0.03

    t0 = time.perf_counter()
    res = grip_polish(
        P_in.copy(), mol, metal=0, donors=[_safe_donor(mol)], geom="",
        clash_weight=5.0, return_diagnostics=True,
    )
    t1 = time.perf_counter()
    rmsd = float(np.sqrt(np.mean(np.sum((res.P - P_in) ** 2, axis=1))))
    return True, (t1 - t0), float(res.severity_after), int(res.n_iter), rmsd


def main():
    print("=" * 78)
    print("GRIP — L-BFGS vs LM (TRF) Benchmark")
    print(f"Sample: {len(SMILES_50)} SMILES (deterministic embeds @ seed=42)")
    print("=" * 78)

    # Per-method tallies
    results: dict = {}
    for method in ("lbfgs", "lm"):
        times = []
        sevs = []
        nits = []
        rmsds = []
        polished_positions = []
        n_ok = 0
        n_total = 0
        for s in SMILES_50:
            ok, t, sev, nit, rmsd = _run_one_smiles(s, method)
            if not ok:
                continue
            n_total += 1
            if np.isfinite(sev):
                n_ok += 1
                times.append(t)
                sevs.append(sev)
                nits.append(nit)
                rmsds.append(rmsd)
        results[method] = {
            "n_ok": n_ok,
            "n_total": n_total,
            "mean_time_ms": 1000.0 * statistics.mean(times) if times else 0.0,
            "total_time_s": sum(times),
            "mean_sev": statistics.mean(sevs) if sevs else float("nan"),
            "mean_nit": statistics.mean(nits) if nits else 0.0,
            "mean_rmsd_from_input": statistics.mean(rmsds) if rmsds else 0.0,
        }

    # Print
    print()
    print(f"{'Metric':36}  {'L-BFGS':>14}  {'LM (TRF)':>14}")
    print("-" * 78)
    keys = [
        ("Successful polishes (n)", "n_ok"),
        ("Total polishes attempted", "n_total"),
        ("Mean wall time (ms)", "mean_time_ms"),
        ("Total wall time (s)", "total_time_s"),
        ("Mean final severity", "mean_sev"),
        ("Mean iterations/nfev", "mean_nit"),
        ("Mean RMSD from input (Å)", "mean_rmsd_from_input"),
    ]
    for label, k in keys:
        v_l = results["lbfgs"][k]
        v_lm = results["lm"][k]
        if isinstance(v_l, float):
            print(f"{label:36}  {v_l:>14.4f}  {v_lm:>14.4f}")
        else:
            print(f"{label:36}  {v_l:>14}  {v_lm:>14}")

    # Compare positions across methods (per-SMILES RMSD)
    print()
    print("Cross-method position comparison (LM vs L-BFGS):")
    cross_rmsds = []
    for s in SMILES_50:
        from delfin.fffree.grip_polish import grip_polish, _GRIP_METHOD_ENV
        mol = _embed(s)
        if mol is None:
            continue
        P0 = _coords(mol)
        rng = np.random.default_rng(7)
        P_in = P0 + rng.standard_normal(P0.shape) * 0.03
        os.environ.pop(_GRIP_METHOD_ENV, None)
        out_lbfgs = grip_polish(P_in.copy(), mol, metal=0, donors=[_safe_donor(mol)],
                                geom="", clash_weight=5.0)
        os.environ[_GRIP_METHOD_ENV] = "lm"
        out_lm = grip_polish(P_in.copy(), mol, metal=0, donors=[_safe_donor(mol)],
                             geom="", clash_weight=5.0)
        os.environ.pop(_GRIP_METHOD_ENV, None)
        rmsd = float(np.sqrt(np.mean(np.sum((out_lbfgs - out_lm) ** 2, axis=1))))
        cross_rmsds.append((s, rmsd))
    n_close = sum(1 for _, r in cross_rmsds if r < 0.01)
    print(f"  pairs with RMSD < 0.01 Å : {n_close}/{len(cross_rmsds)}")
    print(f"  max cross-method RMSD    : {max((r for _, r in cross_rmsds), default=0.0):.4f}")
    print(f"  mean cross-method RMSD   : "
          f"{statistics.mean(r for _, r in cross_rmsds) if cross_rmsds else 0.0:.4f}")

    print()
    print("=" * 78)
    print("HEADLINE:")
    sev_ratio = results["lm"]["mean_sev"] / max(results["lbfgs"]["mean_sev"], 1e-12)
    time_ratio = results["lm"]["mean_time_ms"] / max(results["lbfgs"]["mean_time_ms"], 1e-12)
    print(f"  LM/L-BFGS mean severity ratio: {sev_ratio:.3f}  (<1 = LM lower loss)")
    print(f"  LM/L-BFGS mean time ratio    : {time_ratio:.3f}  (<1 = LM faster)")
    print("=" * 78)


if __name__ == "__main__":
    main()
