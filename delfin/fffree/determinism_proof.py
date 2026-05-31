"""delfin.fffree.determinism_proof — 100% byte-identical builds verification.

Tests determinism of the fffree pipeline by running the same SMILES set TWICE
sequentially and checking byte-identity of the output XYZ files.

Algorithm:
  1. Run pipeline on SMILES set #1 → archive A
  2. Run pipeline on SMILES set #1 (same) → archive B
  3. SHA-256 each XYZ file in both archives
  4. Report: byte-identical fraction

For 100% determinism: every XYZ file must be byte-identical between runs.

Universal validation script. FF-free.
"""
from __future__ import annotations

import hashlib
import os
import subprocess
import sys
from pathlib import Path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def compare_archives(archive_a: Path, archive_b: Path) -> dict:
    """Compare two archives for byte-identity.

    Returns: {total_a, total_b, common_keys, byte_identical, mismatches}
    """
    files_a = {f.name: f for f in archive_a.glob("*.xyz")}
    files_b = {f.name: f for f in archive_b.glob("*.xyz")}
    common = set(files_a.keys()) & set(files_b.keys())
    identical = 0
    mismatches = []
    for name in sorted(common):
        h_a = sha256_file(files_a[name])
        h_b = sha256_file(files_b[name])
        if h_a == h_b:
            identical += 1
        else:
            mismatches.append((name, h_a[:8], h_b[:8]))
    return {
        "total_a": len(files_a),
        "total_b": len(files_b),
        "common": len(common),
        "byte_identical": identical,
        "identical_pct": 100.0 * identical / max(len(common), 1),
        "n_mismatches": len(mismatches),
        "first_mismatches": mismatches[:10],
    }


def run_determinism_check(
    smiles_file: str,
    work_dir: str,
    label_a: str,
    label_b: str,
    env_vars: dict,
    parallel: int = 20,
    timeout: int = 180,
) -> dict:
    """Run pipeline twice + compare byte-identity.

    Returns full comparison report.
    """
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    for label in (label_a, label_b):
        cmd = [
            "python3",
            "scripts/pool_evaluator.py",
            smiles_file,
            "--xyz-archive", str(work / "archives"),
            "--commit-label", label,
            "--parallel", str(parallel),
            "--timeout", str(timeout),
        ]
        env = os.environ.copy()
        env.update(env_vars)
        log = work / f"{label}.log"
        print(f"Running {label} ...")
        with open(log, "w") as f:
            subprocess.run(cmd, env=env, stdout=f, stderr=subprocess.STDOUT)
    archive_a = work / "archives" / label_a
    archive_b = work / "archives" / label_b
    return compare_archives(archive_a, archive_b)


if __name__ == "__main__":
    # Example usage; actual invocation done by user with master_v3_plus
    if len(sys.argv) < 3:
        print("Usage: determinism_proof.py <archive_a> <archive_b>")
        sys.exit(1)
    a = Path(sys.argv[1])
    b = Path(sys.argv[2])
    result = compare_archives(a, b)
    print(f"\n=== Determinism Check ===")
    print(f"Archive A: {result['total_a']} XYZ files")
    print(f"Archive B: {result['total_b']} XYZ files")
    print(f"Common files: {result['common']}")
    print(f"Byte-identical: {result['byte_identical']} ({result['identical_pct']:.2f}%)")
    print(f"Mismatches: {result['n_mismatches']}")
    if result["first_mismatches"]:
        print("First mismatches:")
        for name, a, b in result["first_mismatches"][:5]:
            print(f"  {name}: {a}... vs {b}...")
