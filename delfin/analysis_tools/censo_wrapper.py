"""Python wrapper for CENSO conformer ensemble sorting.

CENSO (Commandline ENergetic SOrting) refines conformer ensembles from CREST
through successive DFT levels, optionally computing thermochemistry (mRRHO),
NMR, UV-Vis, or optical rotation.

Install::

    conda install -c conda-forge censo
    # or: pip install git+https://github.com/grimme-lab/CENSO.git

Requires xTB and/or ORCA for the QM calculations.

Usage::

    from delfin.analysis_tools.censo_wrapper import run_censo, quick_sort

    # Full CENSO run on a CREST ensemble
    result = run_censo("crest_conformers.xyz", charge=0, solvent="chcl3")

    # Quick energy-only sort (prescreening only)
    ranked = quick_sort("crest_conformers.xyz", charge=0)
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Optional

from . import censo_available


def run_censo(
    ensemble_file: str | Path,
    *,
    charge: int = 0,
    uhf: int = 0,
    solvent: str = "",
    functional: str = "r2scan-3c",
    basis: str = "automatic",
    prescreening: bool = True,
    screening: bool = True,
    optimization: bool = False,
    refinement: bool = False,
    nprocs: int = 4,
    omp: int = 2,
    threshold: float = 4.0,
    bhess: bool = True,
    extra_args: Optional[list[str]] = None,
    timeout: int = 7200,
    working_dir: Optional[str | Path] = None,
) -> dict:
    """Run CENSO on a conformer ensemble.

    Parameters
    ----------
    ensemble_file : path to XYZ ensemble (e.g. crest_conformers.xyz)
    charge : molecular charge
    uhf : number of unpaired electrons
    solvent : solvent name (e.g. "chcl3", "water", "dmso") or "" for gas phase
    functional : DFT functional for screening/optimization
    basis : basis set (or "automatic")
    prescreening : enable prescreening step (fast, xTB-level)
    screening : enable screening step (DFT single-points)
    optimization : enable geometry optimization step
    refinement : enable final refinement step
    nprocs : number of parallel CENSO processes (-P)
    omp : OpenMP threads per process (-O)
    threshold : energy threshold in kcal/mol for ensemble pruning
    bhess : use SPH + structure constraints for mRRHO
    extra_args : additional CLI arguments
    timeout : max seconds
    working_dir : directory to run CENSO in (default: ensemble file's directory)

    Returns
    -------
    dict with keys: returncode, stdout, stderr, output_dir
    """
    if not censo_available():
        raise RuntimeError(
            "CENSO is not installed. Install with:\n"
            "  conda install -c conda-forge censo\n"
            "  or: pip install git+https://github.com/grimme-lab/CENSO.git"
        )

    ensemble_path = Path(ensemble_file).resolve()
    if not ensemble_path.is_file():
        raise FileNotFoundError(f"Ensemble file not found: {ensemble_path}")

    work_dir = Path(working_dir) if working_dir else ensemble_path.parent

    cmd = ["censo", "-inp", str(ensemble_path)]
    cmd.extend(["-chrg", str(charge)])
    if uhf:
        cmd.extend(["-uhf", str(uhf)])
    if solvent:
        cmd.extend(["-solvent", solvent])
    cmd.extend(["-func", functional, "-basis", basis])
    cmd.extend(["-P", str(nprocs), "-O", str(omp)])
    cmd.extend(["-thrpart1", str(threshold)])

    if prescreening:
        cmd.extend(["-part1", "on"])
    if screening:
        cmd.extend(["-part2", "on"])
    if optimization:
        cmd.extend(["-part3", "on"])
    if refinement:
        cmd.extend(["-part4", "on"])
    if bhess:
        cmd.append("--bhess")

    if extra_args:
        cmd.extend(extra_args)

    result = subprocess.run(
        cmd,
        cwd=str(work_dir),
        capture_output=True,
        text=True,
        timeout=timeout,
    )

    return {
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_dir": str(work_dir),
    }


def quick_sort(
    ensemble_file: str | Path,
    *,
    charge: int = 0,
    uhf: int = 0,
    solvent: str = "",
    nprocs: int = 4,
    threshold: float = 6.0,
    timeout: int = 3600,
) -> dict:
    """Quick prescreening-only CENSO run for fast ensemble sorting.

    Uses only xTB-level prescreening (part1) — fast but less accurate.
    Good for initial filtering before expensive DFT calculations.
    """
    return run_censo(
        ensemble_file,
        charge=charge,
        uhf=uhf,
        solvent=solvent,
        prescreening=True,
        screening=False,
        optimization=False,
        refinement=False,
        nprocs=nprocs,
        threshold=threshold,
        timeout=timeout,
    )


def parse_censo_energies(censo_output_dir: str | Path) -> list[dict]:
    """Parse CENSO output to extract ranked conformer energies.

    Parameters
    ----------
    censo_output_dir : directory containing CENSO output files

    Returns
    -------
    list of dicts with keys: conformer, energy_au, population, dG_kcal
    sorted by energy ascending
    """
    output_dir = Path(censo_output_dir)
    results = []

    # Try to parse the final enso.json or results files
    enso_json = output_dir / "enso.json"
    if enso_json.is_file():
        import json

        data = json.loads(enso_json.read_text())
        for conf_name, conf_data in data.items():
            if not isinstance(conf_data, dict):
                continue
            energy = conf_data.get("energy")
            if energy is not None:
                results.append({
                    "conformer": conf_name,
                    "energy_au": float(energy),
                    "population": conf_data.get("population"),
                    "dG_kcal": conf_data.get("rel_free_energy"),
                })

    # Fallback: parse censo.out text
    if not results:
        censo_out = output_dir / "censo.out"
        if censo_out.is_file():
            import re

            text = censo_out.read_text()
            # Look for energy table lines
            pattern = re.compile(
                r"(CONF\d+)\s+(-?\d+\.\d+)\s+"
            )
            for match in pattern.finditer(text):
                results.append({
                    "conformer": match.group(1),
                    "energy_au": float(match.group(2)),
                    "population": None,
                    "dG_kcal": None,
                })

    results.sort(key=lambda r: r["energy_au"])
    return results
