"""Cluster helper: SMILES -> isomer XYZ generation with parallel UFF in tray/.

This module is intended to be executed inside a submitted batch job.
It reads one SMILES from input.txt, generates coordination isomers, applies
Open Babel UFF refinement in parallel, and writes all outputs into tray/.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import re
from pathlib import Path
from typing import List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.smiles_converter import smiles_to_xyz_isomers, _optimize_xyz_openbabel

logger = get_logger(__name__)


def _read_first_smiles(path: Path) -> Optional[str]:
    if not path.exists():
        return None
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        return line
    return None


def _safe_label(label: str, idx: int) -> str:
    base = (label or f"isomer_{idx}").strip()
    base = re.sub(r"[^A-Za-z0-9._-]+", "_", base)
    return base[:80] if base else f"isomer_{idx}"


def _xyz_block(xyz_body: str, comment: str) -> str:
    lines = [ln.rstrip() for ln in xyz_body.splitlines() if ln.strip()]
    return f"{len(lines)}\n{comment}\n" + "\n".join(lines) + "\n"


def _optimize_entry(entry: Tuple[int, str, str]) -> Tuple[int, str, str]:
    idx, xyz, label = entry
    try:
        return idx, _optimize_xyz_openbabel(xyz), label
    except Exception as exc:
        logger.debug("UFF optimization failed for isomer %s (%s), keeping original.", idx, exc)
        return idx, xyz, label


def run_cluster_smiles_convert(
    input_file: Path,
    out_dir: Path,
    max_isomers: int,
    parallel_jobs: int,
    apply_uff: bool,
) -> int:
    smiles = _read_first_smiles(input_file)
    if not smiles:
        logger.error("No SMILES found in %s", input_file)
        return 1

    logger.info("Starting SMILES conversion: %s", smiles)
    # Generate full isomer pool first without UFF so we can refine in parallel.
    isomers, error = smiles_to_xyz_isomers(
        smiles,
        max_isomers=max_isomers,
        apply_uff=False,
    )
    if error or not isomers:
        logger.error("SMILES conversion failed: %s", error or "No isomers generated")
        return 2

    if apply_uff:
        jobs = max(1, int(parallel_jobs))
        workers = min(jobs, len(isomers))
        entries = [(i, xyz, label) for i, (xyz, label) in enumerate(isomers, 1)]
        if workers > 1:
            logger.info("Running parallel UFF refinement with %d workers", workers)
            with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
                refined = list(pool.map(_optimize_entry, entries))
            refined.sort(key=lambda x: x[0])
            isomers = [(xyz, label) for _i, xyz, label in refined]
        else:
            isomers = [(_optimize_entry((i, xyz, label))[1], label) for i, (xyz, label) in enumerate(isomers, 1)]

    out_dir.mkdir(parents=True, exist_ok=True)
    all_blocks: List[str] = []
    labels: List[str] = []

    for idx, (xyz, label) in enumerate(isomers, 1):
        label_txt = (label or f"isomer_{idx}").strip()
        labels.append(f"{idx}\t{label_txt}")
        safe = _safe_label(label_txt, idx)
        block = _xyz_block(xyz, f"SMILES cluster convert | isomer {idx} | {label_txt}")
        all_blocks.append(block)
        (out_dir / f"isomer_{idx:03d}_{safe}.xyz").write_text(block, encoding="utf-8")

    if all_blocks:
        (out_dir / "SMILES_convert_best.xyz").write_text(all_blocks[0], encoding="utf-8")
    (out_dir / "SMILES_convert_all.xyz").write_text("".join(all_blocks), encoding="utf-8")
    (out_dir / "SMILES_convert_labels.txt").write_text("\n".join(labels) + "\n", encoding="utf-8")
    (out_dir / "SMILES_convert_input.txt").write_text(smiles + "\n", encoding="utf-8")

    logger.info("Wrote %d isomers to %s", len(isomers), out_dir)
    logger.info("Best structure: %s", out_dir / "SMILES_convert_best.xyz")
    logger.info("All structures: %s", out_dir / "SMILES_convert_all.xyz")
    return 0


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="delfin-smiles-cluster-convert",
        description="Generate SMILES isomers with optional parallel UFF refinement.",
    )
    parser.add_argument("input_file", nargs="?", default="input.txt")
    parser.add_argument("--out-dir", default="tray")
    parser.add_argument("--max-isomers", type=int, default=400)
    parser.add_argument("--parallel-jobs", type=int, default=40)
    parser.add_argument("--no-uff", action="store_true", help="Disable UFF refinement")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    return run_cluster_smiles_convert(
        input_file=Path(args.input_file),
        out_dir=Path(args.out_dir),
        max_isomers=max(1, int(args.max_isomers)),
        parallel_jobs=max(1, int(args.parallel_jobs)),
        apply_uff=not args.no_uff,
    )


if __name__ == "__main__":
    raise SystemExit(main())
