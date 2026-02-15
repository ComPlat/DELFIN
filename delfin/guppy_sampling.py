"""Generate and rank multiple XTB-optimized structures from one SMILES input.

Workflow:
1. Read SMILES from input file (default: input.txt)
2. Repeat N times (default: 20):
   - convert SMILES to a 3D start structure
   - run ORCA XTB2 OPT in run_XX/XTB2/
   - extract energy from output_XTB.out
   - collect optimized geometry from XTB.xyz
3. Sort successful runs by energy
4. Write trajectory to GUPPY_try.xyz (XYZ multi-frame)

Comment line format in trajectory:
    run_XX <energy_in_Eh>
so the energy is the second column as requested.
"""

from __future__ import annotations

import argparse
import os
import re
import threading
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from delfin.common.logging import get_logger
from delfin.dynamic_pool import JobPriority, PoolJob
from delfin.global_manager import get_global_manager
from delfin.orca import run_orca
from delfin.smiles_converter import RDKIT_AVAILABLE, smiles_to_xyz

if RDKIT_AVAILABLE:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from delfin.smiles_converter import _mol_to_xyz, _prepare_mol_for_embedding

logger = get_logger(__name__)

_TOTAL_ENERGY_RE = re.compile(
    r"total\s+energy\s+([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s+Eh\b",
    re.IGNORECASE,
)
_FINAL_SP_ENERGY_RE = re.compile(
    r"final\s+single\s+point\s+energy\s+([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\b",
    re.IGNORECASE,
)
_BRACKET_TOKEN_RE = re.compile(r"\[([^\]]+)\]")


def _derive_charge_from_smiles(smiles: str) -> int:
    """Derive total charge from full SMILES (metal + ligands).

    Examples:
    - [Fe+2] -> +2
    - [CH-]  -> -1
    - [N+]   -> +1
    - [Cu++] -> +2
    """
    if RDKIT_AVAILABLE:
        try:
            parser_params = Chem.SmilesParserParams()
            parser_params.sanitize = False
            parser_params.removeHs = False
            parser_params.strictParsing = False
            mol = Chem.MolFromSmiles(smiles, parser_params)
            if mol is not None:
                return int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms()))
        except Exception as exc:  # noqa: BLE001
            logger.debug("RDKit formal-charge extraction failed, falling back to token parser: %s", exc)

    # Fallback parser: explicit +/- annotations in bracket atoms.
    total = 0
    for token in _BRACKET_TOKEN_RE.findall(smiles):
        i = 0
        n = len(token)
        while i < n:
            ch = token[i]
            if ch not in "+-":
                i += 1
                continue

            sign = 1 if ch == "+" else -1

            # Handle repeated signs like "++" / "--"
            j = i
            while j < n and token[j] == ch:
                j += 1
            repeated = j - i

            # Handle optional magnitude digits after sign(s): +2 / -3
            k = j
            while k < n and token[k].isdigit():
                k += 1

            if k > j:
                magnitude = int(token[j:k])
                total += sign * magnitude
                i = k
            else:
                total += sign * repeated
                i = j
    return total


def _read_first_smiles_line(input_path: Path) -> str:
    """Return first non-empty, non-comment line from input file."""
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    for line in input_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("*"):
            continue
        return stripped
    raise ValueError(f"No SMILES line found in {input_path}")


def _convert_smiles_with_seed(smiles: str, seed: int) -> Tuple[Optional[str], Optional[str]]:
    """Convert SMILES to XYZ coordinates (DELFIN coordinate format, no header).

    Uses deterministic per-run seeds for diverse but reproducible starts.
    Falls back to regular smiles_to_xyz if seeded embedding is unavailable.
    """
    if RDKIT_AVAILABLE:
        try:
            mol = _prepare_mol_for_embedding(smiles)
        except Exception as exc:  # noqa: BLE001
            logger.debug("Seeded embedding prep failed: %s", exc)
            mol = None

        if mol is not None:
            try:
                mol.RemoveAllConformers()
                params = AllChem.ETKDGv3()
                params.useRandomCoords = True
                params.randomSeed = int(seed)
                params.enforceChirality = False
                result = AllChem.EmbedMolecule(mol, params)
                if result == 0:
                    return _mol_to_xyz(mol), None
            except Exception as exc:  # noqa: BLE001
                logger.debug("Seeded embedding failed: %s", exc)

    return smiles_to_xyz(smiles)


def _write_xtb_input(
    inp_path: Path,
    coords_lines: List[str],
    *,
    charge: int,
    multiplicity: int,
    pal: int,
    maxcore: int,
    method: str,
) -> None:
    """Write ORCA XTB optimization input file."""
    blocks = [
        f"!{method} OPT",
        f"%maxcore {maxcore}",
        f"%pal nprocs {pal} end",
        f"*xyz {charge} {multiplicity}",
    ]
    blocks.extend(coords_lines)
    blocks.append("*")
    inp_path.write_text("\n".join(blocks) + "\n", encoding="utf-8")


def _extract_total_energy_eh(output_path: Path) -> Optional[float]:
    """Extract last energy value from ORCA output.

    Supports both:
    - '... total energy ... Eh ...' (XTB printout)
    - 'FINAL SINGLE POINT ENERGY ...' (ORCA summary)
    """
    if not output_path.exists():
        return None
    energy = None
    for line in output_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        for pattern in (_TOTAL_ENERGY_RE, _FINAL_SP_ENERGY_RE):
            match = pattern.search(line)
            if not match:
                continue
            try:
                energy = float(match.group(1))
            except ValueError:
                continue
    return energy


def _read_xyz_coordinates(xyz_path: Path) -> Tuple[int, List[str]]:
    """Read XYZ file and return (natoms, coordinate_lines_without_header)."""
    if not xyz_path.exists():
        raise FileNotFoundError(f"Missing XYZ file: {xyz_path}")

    lines = [ln.rstrip() for ln in xyz_path.read_text(encoding="utf-8", errors="ignore").splitlines()]
    non_empty = [ln for ln in lines if ln.strip()]
    if not non_empty:
        raise ValueError(f"Empty XYZ file: {xyz_path}")

    natoms = None
    coord_start = 0
    try:
        natoms = int(non_empty[0].split()[0])
        coord_start = 2
    except (ValueError, IndexError):
        natoms = None
        coord_start = 0

    coords = non_empty[coord_start:]
    if natoms is None:
        natoms = len(coords)
    else:
        coords = coords[:natoms]

    if natoms <= 0 or len(coords) < natoms:
        raise ValueError(f"Invalid XYZ content in {xyz_path}")

    return natoms, coords


def _write_ranked_trajectory(
    output_path: Path,
    ranked_results: List[Tuple[float, int, List[str], int]],
) -> None:
    """Write sorted structures as multi-frame XYZ trajectory."""
    with output_path.open("w", encoding="utf-8") as handle:
        for energy, natoms, coords, run_idx in ranked_results:
            handle.write(f"{natoms}\n")
            handle.write(f"run_{run_idx:02d} {energy:.12f}\n")
            for line in coords[:natoms]:
                handle.write(f"{line}\n")


def _execute_single_sampling_run(
    *,
    run_idx: int,
    run_seed: int,
    smiles: str,
    resolved_charge: int,
    multiplicity: int,
    pal: int,
    maxcore: int,
    method: str,
    workdir: Path,
) -> Tuple[bool, Optional[Tuple[float, int, List[str], int]], Optional[str]]:
    """Execute one SMILES->XTB2 run and return (ok, result, error)."""
    run_dir = workdir / f"run_{run_idx:02d}"
    xtb_dir = run_dir / "XTB2"
    xtb_dir.mkdir(parents=True, exist_ok=True)

    xyz_text, error = _convert_smiles_with_seed(smiles, run_seed)
    if error or not xyz_text:
        return False, None, f"SMILES conversion failed: {error}"

    coords_lines = [ln.rstrip() for ln in xyz_text.splitlines() if ln.strip()]
    if not coords_lines:
        return False, None, "Converted XYZ is empty"

    start_xyz = run_dir / "start_converted.xyz"
    start_xyz.write_text(
        f"{len(coords_lines)}\nrun_{run_idx:02d} start\n" + "\n".join(coords_lines) + "\n",
        encoding="utf-8",
    )

    inp_path = xtb_dir / "XTB.inp"
    out_path = xtb_dir / "output_XTB.out"
    xyz_path = xtb_dir / "XTB.xyz"
    _write_xtb_input(
        inp_path,
        coords_lines,
        charge=resolved_charge,
        multiplicity=multiplicity,
        pal=pal,
        maxcore=maxcore,
        method=method,
    )

    ok = run_orca(
        str(inp_path),
        str(out_path),
        working_dir=xtb_dir,
        isolate=True,
    )
    if not ok:
        return False, None, "ORCA XTB run failed"

    energy = _extract_total_energy_eh(out_path)
    if energy is None:
        return False, None, f"Could not extract total energy from {out_path}"

    try:
        natoms, opt_coords = _read_xyz_coordinates(xyz_path)
    except Exception as exc:  # noqa: BLE001
        return False, None, f"Could not read optimized XYZ: {exc}"

    return True, (energy, natoms, opt_coords, run_idx), None


def run_sampling(
    *,
    input_file: Path,
    runs: int,
    charge: Optional[int],
    pal: int,
    maxcore: int,
    parallel_jobs: int,
    method: str,
    output_file: Path,
    workdir: Path,
    seed: int,
    allow_partial: bool,
) -> int:
    """Execute repeated SMILES->XTB2 workflow and write ranked trajectory."""
    smiles = _read_first_smiles_line(input_file)
    resolved_charge = charge if charge is not None else _derive_charge_from_smiles(smiles)
    resolved_parallel_jobs = max(1, min(parallel_jobs, runs, pal))
    per_job_pal = max(1, pal // resolved_parallel_jobs)

    logger.info("Using SMILES from %s", input_file)
    logger.info("Sampling runs: %d", runs)
    logger.info("Using charge: %d (net charge from whole SMILES: metal + ligands)", resolved_charge)
    logger.info("Total PAL budget: %d", pal)
    logger.info("Maxcore per core: %d MB", maxcore)
    logger.info("Parallel jobs: %d", resolved_parallel_jobs)
    logger.info("Target PAL per run: %d", per_job_pal)
    if charge is None:
        logger.info("Charge was auto-derived from SMILES explicit +/- annotations.")

    workdir.mkdir(parents=True, exist_ok=True)
    results: List[Tuple[float, int, List[str], int]] = []
    failed_runs: List[int] = []
    results_lock = threading.Lock()

    # Closed-shell only as requested.
    multiplicity = 1
    logger.info("Using multiplicity: %d (fixed closed-shell)", multiplicity)

    def _record_run(run_idx: int, run_seed: int, assigned_pal: int) -> None:
        try:
            ok, result, error = _execute_single_sampling_run(
                run_idx=run_idx,
                run_seed=run_seed,
                smiles=smiles,
                resolved_charge=resolved_charge,
                multiplicity=multiplicity,
                pal=max(1, assigned_pal),
                maxcore=maxcore,
                method=method,
                workdir=workdir,
            )
            with results_lock:
                if ok and result is not None:
                    energy = result[0]
                    logger.info("[run %02d] energy = %.12f Eh (PAL=%d)", run_idx, energy, max(1, assigned_pal))
                    results.append(result)
                else:
                    logger.error("[run %02d] %s", run_idx, error or "Unknown run failure")
                    failed_runs.append(run_idx)
        except Exception as exc:  # noqa: BLE001
            with results_lock:
                logger.error("[run %02d] Unexpected error: %s", run_idx, exc)
                failed_runs.append(run_idx)

    use_pool = False
    pool = None
    try:
        manager = get_global_manager()
        manager.ensure_initialized(
            {
                "PAL": pal,
                "maxcore": maxcore,
                "pal_jobs": resolved_parallel_jobs,
                "parallel_workflows": "enable",
            }
        )
        pool = manager.get_pool()
        use_pool = True
        logger.info("Using global manager pool scheduling for GUPPY runs.")
    except Exception as exc:  # noqa: BLE001
        logger.warning("Global manager unavailable for GUPPY (%s). Falling back to sequential runs.", exc)

    if use_pool and pool is not None:
        estimated_runtime = float(max(300, int(os.environ.get("GUPPY_EST_RUNTIME_S", "1800"))))
        for run_idx in range(1, runs + 1):
            run_seed = seed + (run_idx - 1) * 1009
            run_dir = workdir / f"run_{run_idx:02d}"

            def runner(*_args, cur_idx=run_idx, cur_seed=run_seed, **kwargs) -> None:
                allocated = kwargs.get("cores", per_job_pal)
                try:
                    assigned = max(1, int(allocated))
                except (TypeError, ValueError):
                    assigned = per_job_pal
                _record_run(cur_idx, cur_seed, assigned)

            pool_job = PoolJob(
                job_id=f"GUPPY_RUN_{run_idx:02d}",
                cores_min=1,
                cores_optimal=per_job_pal,
                cores_max=per_job_pal,
                memory_mb=max(256, per_job_pal * maxcore),
                priority=JobPriority.NORMAL,
                execute_func=runner,
                args=(),
                kwargs={},
                estimated_duration=estimated_runtime,
                working_dir=run_dir,
            )
            pool_job.suppress_pool_logs = True
            pool.submit_job(pool_job)

        pool.wait_for_completion()
    else:
        for run_idx in range(1, runs + 1):
            run_seed = seed + (run_idx - 1) * 1009
            _record_run(run_idx, run_seed, per_job_pal)

    if not results:
        logger.error("No successful XTB runs. Nothing to write.")
        return 1

    results.sort(key=lambda item: (item[0], item[3]))
    _write_ranked_trajectory(output_file, results)

    logger.info("Wrote ranked trajectory: %s", output_file)
    logger.info("Successful runs: %d / %d", len(results), runs)
    if failed_runs:
        logger.warning("Failed runs: %s", ", ".join(str(i) for i in failed_runs))

    if failed_runs and not allow_partial:
        logger.error("Partial result detected and --allow-partial not set.")
        return 1
    return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="delfin-guppy",
        description="Run repeated SMILES->XTB2 OPT sampling and rank structures by energy.",
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        default="input.txt",
        help="Input file containing SMILES in first non-empty line (default: input.txt)",
    )
    parser.add_argument("--runs", type=int, default=int(os.environ.get("GUPPY_RUNS", "20")))
    parser.add_argument(
        "--charge",
        type=int,
        default=(int(os.environ["GUPPY_CHARGE"]) if "GUPPY_CHARGE" in os.environ else None),
        help="Total charge override. If omitted, net charge is derived from full SMILES (metal + ligands).",
    )
    parser.add_argument(
        "--pal",
        type=int,
        default=int(os.environ.get("GUPPY_PAL", os.environ.get("SLURM_CPUS_PER_TASK", "40"))),
        help="Total core budget for all GUPPY runs (default: $GUPPY_PAL or $SLURM_CPUS_PER_TASK or 40)",
    )
    parser.add_argument(
        "--maxcore",
        type=int,
        default=int(os.environ.get("GUPPY_MAXCORE", os.environ.get("DELFIN_MAXCORE", "6000"))),
        help="Memory budget per core in MB for pool scheduling (default: $GUPPY_MAXCORE or $DELFIN_MAXCORE or 6000)",
    )
    parser.add_argument(
        "--parallel-jobs",
        type=int,
        default=int(os.environ.get("GUPPY_PARALLEL_JOBS", "4")),
        help="Maximum number of parallel GUPPY runs sharing PAL/maxcore (default: 4)",
    )
    parser.add_argument("--method", default=os.environ.get("GUPPY_XTB_METHOD", "XTB2"))
    parser.add_argument("--output", default="GUPPY_try.xyz")
    parser.add_argument("--workdir", default="GUPPY")
    parser.add_argument("--seed", type=int, default=int(os.environ.get("GUPPY_SEED", "31")))
    parser.add_argument(
        "--allow-partial",
        action="store_true",
        help="Return success even if some of the runs fail.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.runs <= 0:
        parser.error("--runs must be > 0")
    if args.pal <= 0:
        parser.error("--pal must be > 0")
    if args.maxcore <= 0:
        parser.error("--maxcore must be > 0")
    if args.parallel_jobs <= 0:
        parser.error("--parallel-jobs must be > 0")

    return run_sampling(
        input_file=Path(args.input_file),
        runs=args.runs,
        charge=args.charge,
        pal=args.pal,
        maxcore=args.maxcore,
        parallel_jobs=args.parallel_jobs,
        method=str(args.method).strip() or "XTB2",
        output_file=Path(args.output),
        workdir=Path(args.workdir),
        seed=args.seed,
        allow_partial=args.allow_partial,
    )


if __name__ == "__main__":
    raise SystemExit(main())
