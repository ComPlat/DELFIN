"""Smoke validation: GRIP-Healing-Mode on real broken structures.

Reads a detector jsonl (default: VOLLPOOL bondlen.jsonl), picks 20 files
with the highest ``n_outliers``, loads the XYZ frames, and applies the
healer.  Reports the fraction of broken bonds that drop below the
σ-threshold after healing.

Usage::

    PYTHONHASHSEED=0 micromamba run -n delfin python scripts/smoke_grip_healing.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

os.environ.setdefault("PYTHONHASHSEED", "0")

# Ensure the WORKTREE's delfin package is imported, not the parent-repo's
# installed delfin (which lacks grip_healing).
_HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_HERE))

import numpy as np

# Bond-length table — covalent radii (sum) for the most common pairs.  This
# is the same fallback used by the heal when the GRIP library has no entry.
_COV = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
    "B": 0.84, "Si": 1.11, "Se": 1.20,
}
_METALS = set(
    "Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os "
    "Ir Pt Au Hg La Ce Lu Sn Pb Ge Sb Bi".split()
)


def _ideal_bond(a: str, b: str) -> float:
    return _COV.get(a, 0.9) + _COV.get(b, 0.9)


def _is_metal(s: str) -> bool:
    return s in _METALS


def parse_xyz(text: str) -> Tuple[List[str], np.ndarray]:
    lines = text.splitlines()
    syms: List[str] = []
    pts: List[List[float]] = []
    for ln in lines:
        p = ln.split()
        if len(p) == 4:
            try:
                xyz = [float(p[1]), float(p[2]), float(p[3])]
            except ValueError:
                continue
            syms.append(p[0])
            pts.append(xyz)
    return syms, np.array(pts, dtype=float)


def geometric_bonds(
    syms: List[str], P: np.ndarray,
) -> List[Tuple[int, int]]:
    """Heavy-heavy + X-H pairs within 1.3 * Σcov are treated as bonds."""
    n = len(syms)
    bonds: List[Tuple[int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(P[i] - P[j]))
            if d < 1.30 * _ideal_bond(syms[i], syms[j]):
                bonds.append((i, j))
    return bonds


def count_broken_bonds(
    syms: List[str], P: np.ndarray, bonds: List[Tuple[int, int]],
) -> int:
    """Count bonds whose length deviates > 25% from the element-pair ideal."""
    cnt = 0
    for (i, j) in bonds:
        d = float(np.linalg.norm(P[i] - P[j]))
        ideal = _ideal_bond(syms[i], syms[j])
        if ideal <= 0:
            continue
        rel = abs(d - ideal) / ideal
        if rel > 0.25:
            cnt += 1
    return cnt


def find_metal_donor_set(
    syms: List[str], P: np.ndarray,
) -> Tuple[int, List[int]]:
    """Return ``(metal_idx, donor_indices)`` -- first metal + atoms within
    1.3 * (cov_M + cov_X) of it.  Returns ``(-1, [])`` if no metal found."""
    metal_idx = -1
    for i, s in enumerate(syms):
        if _is_metal(s):
            metal_idx = i
            break
    if metal_idx < 0:
        return -1, []
    donors: List[int] = []
    for j in range(len(syms)):
        if j == metal_idx:
            continue
        d = float(np.linalg.norm(P[metal_idx] - P[j]))
        if d < 1.30 * _ideal_bond(syms[metal_idx], syms[j]):
            donors.append(j)
    return metal_idx, sorted(donors)


def heal_structure(
    syms: List[str], P: np.ndarray,
) -> Tuple[np.ndarray, "object"]:
    """Run the healer on one XYZ.  Returns ``(P_healed, diagnostics)``.

    Uses the GRIP library if available; otherwise relies on covalent-radius
    fallback constants.  Freezes metal + first-coord donors during heal.
    """
    from delfin.fffree.grip_healing import iterative_topology_repositioning

    bonds = geometric_bonds(syms, P)
    metal_idx, donors = find_metal_donor_set(syms, P)
    frozen = [metal_idx, *donors] if metal_idx >= 0 else []

    # Build ideal table from covalent radii.  We do not require the
    # GRIP library here -- the σ-threshold uses a wide fallback (1.0)
    # so the heal still triggers for severely broken bonds (> 30% off).
    ideal = {}
    for (i, j) in bonds:
        if _is_metal(syms[i]) or _is_metal(syms[j]):
            continue  # don't heal M-D bonds (frozen anyway)
        mu = _ideal_bond(syms[i], syms[j])
        # Use σ = 0.10 (typical bond stiffness band, ~ 7% of bond length).
        ideal[(min(i, j), max(i, j))] = (float(mu), 0.10)

    # Filter bonds to non-metal-incident ones.
    heal_bonds = [
        (a, b) for (a, b) in bonds
        if not _is_metal(syms[a]) and not _is_metal(syms[b])
    ]

    P_healed, diag = iterative_topology_repositioning(
        P.copy(), heal_bonds,
        ideal_lengths=ideal,
        symbols=syms,
        frozen_atoms=frozen,
        sigma_threshold=3.0,
        max_iter=50,
        tol=0.01,
        return_diagnostics=True,
    )
    return P_healed, diag


def main():
    pool_dir = Path(
        "/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/"
        "b00f9a0-full7-VOLLPOOL"
    )
    bondlen_jsonl = Path(
        "/home/qmchem_max/agent_workspace/quality_framework/xyz_archive/"
        "b00f9a0-full7-VOLLPOOL_bondlen.jsonl"
    )
    if not pool_dir.exists():
        print(f"Pool dir missing: {pool_dir}")
        return 1
    if not bondlen_jsonl.exists():
        print(f"bondlen jsonl missing: {bondlen_jsonl}")
        return 1

    # Aggregate worst-n_outliers per file (over all frames).
    worst_outliers: Dict[str, int] = {}
    with bondlen_jsonl.open() as fh:
        for line in fh:
            try:
                rec = json.loads(line)
            except Exception:
                continue
            f = str(rec.get("file", ""))
            n = int(rec.get("n_outliers", 0))
            if n > worst_outliers.get(f, 0):
                worst_outliers[f] = n
    if not worst_outliers:
        print("No bondlen records found")
        return 1

    # Pick top 20 files with the most outliers.
    top = sorted(worst_outliers.items(), key=lambda t: -t[1])[:20]
    print(f"# Smoke validation: GRIP-Healing-Mode on {len(top)} broken structures")
    print(f"# Pool: {pool_dir.name}")
    print(f"# Detector source: {bondlen_jsonl.name}")
    print()
    print(f"{'file':<40} {'bonds':>6} {'broken0':>8} {'broken1':>8} {'iters':>6} {'conv':>5} {'reduced%':>9}")
    print("-" * 90)

    n_files = 0
    n_files_healed = 0
    total_broken_before = 0
    total_broken_after = 0

    for fname, _n_out in top:
        xyz_path = pool_dir / fname
        if not xyz_path.exists():
            continue
        try:
            text = xyz_path.read_text()
        except Exception:
            continue
        syms, P = parse_xyz(text)
        if len(syms) < 3:
            continue

        bonds = geometric_bonds(syms, P)
        broken_before = count_broken_bonds(syms, P, bonds)
        if broken_before == 0:
            continue

        try:
            P_healed, diag = heal_structure(syms, P)
        except Exception as exc:
            print(f"{fname:<40} HEAL FAILED: {exc!r}")
            continue
        broken_after = count_broken_bonds(syms, P_healed, bonds)
        n_files += 1
        if broken_after < broken_before:
            n_files_healed += 1
        total_broken_before += broken_before
        total_broken_after += broken_after
        reduced_pct = 0.0
        if broken_before > 0:
            reduced_pct = 100.0 * (broken_before - broken_after) / broken_before
        conv = "Y" if diag.converged else "N"
        print(
            f"{fname:<40} {len(bonds):>6d} {broken_before:>8d} "
            f"{broken_after:>8d} {diag.n_iter:>6d} {conv:>5s} "
            f"{reduced_pct:>9.1f}"
        )

    print()
    print(f"# Files processed: {n_files}")
    print(f"# Files with > 0 broken bonds: {n_files}")
    print(f"# Files where broken count decreased: {n_files_healed}")
    if n_files > 0:
        pct_healed = 100.0 * n_files_healed / n_files
        print(f"# % structures successfully (partially) healed: {pct_healed:.1f}%")
    if total_broken_before > 0:
        reduction = 100.0 * (total_broken_before - total_broken_after) / total_broken_before
        print(f"# Total broken bonds before: {total_broken_before}")
        print(f"# Total broken bonds after:  {total_broken_after}")
        print(f"# Aggregate bond-break reduction: {reduction:.1f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
