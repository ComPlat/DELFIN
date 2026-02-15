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
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from delfin.common.logging import get_logger
from delfin.orca import run_orca
from delfin.smiles_converter import RDKIT_AVAILABLE, smiles_to_xyz

if RDKIT_AVAILABLE:
    from rdkit.Chem import AllChem
    from delfin.smiles_converter import _mol_to_xyz, _prepare_mol_for_embedding

logger = get_logger(__name__)

_TOTAL_ENERGY_RE = re.compile(
    r"^\s*total\s+energy\s+([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s+Eh\b",
    re.IGNORECASE,
)
_BRACKET_TOKEN_RE = re.compile(r"\[([^\]]+)\]")


def _derive_charge_from_smiles(smiles: str) -> int:
    """Derive total charge from explicit +/- markers inside bracket atoms.

    Examples:
    - [Fe+2] -> +2
    - [CH-]  -> -1
    - [N+]   -> +1
    - [Cu++] -> +2
    """
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
    method: str,
) -> None:
    """Write ORCA XTB optimization input file."""
    blocks = [
        f"!{method} OPT",
        f"%pal nprocs {pal} end",
        f"*xyz {charge} {multiplicity}",
    ]
    blocks.extend(coords_lines)
    blocks.append("*")
    inp_path.write_text("\n".join(blocks) + "\n", encoding="utf-8")


def _extract_total_energy_eh(output_path: Path) -> Optional[float]:
    """Extract last 'total energy ... Eh' value from ORCA output."""
    if not output_path.exists():
        return None
    energy = None
    for line in output_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        match = _TOTAL_ENERGY_RE.search(line)
        if match:
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


def run_sampling(
    *,
    input_file: Path,
    runs: int,
    charge: Optional[int],
    pal: int,
    method: str,
    output_file: Path,
    workdir: Path,
    seed: int,
    allow_partial: bool,
) -> int:
    """Execute repeated SMILES->XTB2 workflow and write ranked trajectory."""
    smiles = _read_first_smiles_line(input_file)
    resolved_charge = charge if charge is not None else _derive_charge_from_smiles(smiles)
    logger.info("Using SMILES from %s", input_file)
    logger.info("Sampling runs: %d", runs)
    logger.info("Using charge: %d", resolved_charge)
    if charge is None:
        logger.info("Charge was auto-derived from SMILES explicit +/- annotations.")

    workdir.mkdir(parents=True, exist_ok=True)
    results: List[Tuple[float, int, List[str], int]] = []
    failed_runs: List[int] = []

    # Closed-shell only as requested.
    multiplicity = 1
    logger.info("Using multiplicity: %d (fixed closed-shell)", multiplicity)

    for run_idx in range(1, runs + 1):
        run_seed = seed + (run_idx - 1) * 1009
        run_dir = workdir / f"run_{run_idx:02d}"
        xtb_dir = run_dir / "XTB2"
        xtb_dir.mkdir(parents=True, exist_ok=True)

        xyz_text, error = _convert_smiles_with_seed(smiles, run_seed)
        if error or not xyz_text:
            logger.error("[run %02d] SMILES conversion failed: %s", run_idx, error)
            failed_runs.append(run_idx)
            continue

        coords_lines = [ln.rstrip() for ln in xyz_text.splitlines() if ln.strip()]
        if not coords_lines:
            logger.error("[run %02d] Converted XYZ is empty", run_idx)
            failed_runs.append(run_idx)
            continue

        # Keep the generated start geometry for traceability.
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
            method=method,
        )

        ok = run_orca(str(inp_path), str(out_path))
        if not ok:
            logger.error("[run %02d] ORCA XTB run failed", run_idx)
            failed_runs.append(run_idx)
            continue

        energy = _extract_total_energy_eh(out_path)
        if energy is None:
            logger.error("[run %02d] Could not extract total energy from %s", run_idx, out_path)
            failed_runs.append(run_idx)
            continue

        try:
            natoms, opt_coords = _read_xyz_coordinates(xyz_path)
        except Exception as exc:  # noqa: BLE001
            logger.error("[run %02d] Could not read optimized XYZ: %s", run_idx, exc)
            failed_runs.append(run_idx)
            continue

        logger.info("[run %02d] energy = %.12f Eh", run_idx, energy)
        results.append((energy, natoms, opt_coords, run_idx))

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
        help="Total charge override. If omitted, charge is derived from SMILES +/- markers.",
    )
    parser.add_argument(
        "--pal",
        type=int,
        default=int(os.environ.get("GUPPY_PAL", os.environ.get("SLURM_CPUS_PER_TASK", "40"))),
        help="ORCA PAL for XTB runs (default: $GUPPY_PAL or $SLURM_CPUS_PER_TASK or 40)",
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

    return run_sampling(
        input_file=Path(args.input_file),
        runs=args.runs,
        charge=args.charge,
        pal=args.pal,
        method=str(args.method).strip() or "XTB2",
        output_file=Path(args.output),
        workdir=Path(args.workdir),
        seed=args.seed,
        allow_partial=args.allow_partial,
    )


if __name__ == "__main__":
    raise SystemExit(main())
