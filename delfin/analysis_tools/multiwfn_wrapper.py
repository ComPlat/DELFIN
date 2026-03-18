"""Python wrapper for Multiwfn wavefunction analysis.

Multiwfn is a standalone binary driven via interactive stdin menus.
This module automates common analyses by piping menu selections via subprocess.

Requires:
  - Multiwfn binary in PATH (download from http://sobereva.com/multiwfn/)
  - Input files: .molden, .wfn, .wfx, or .fch (generate from ORCA via orca_2mkl)

ORCA workflow::

    # Generate .molden from ORCA .gbw
    orca_2mkl basename -molden

    # Then analyse
    from delfin.analysis_tools.multiwfn_wrapper import bond_order_analysis
    result = bond_order_analysis("basename.molden.input")
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path
from typing import Optional

from . import multiwfn_available


def _run_multiwfn(
    input_file: str | Path,
    commands: str,
    *,
    timeout: int = 300,
) -> str:
    """Run Multiwfn with the given menu commands piped to stdin.

    Parameters
    ----------
    input_file : path to .molden / .wfn / .wfx / .fch file
    commands : multi-line string of menu selections (one per line)
    timeout : max seconds to wait

    Returns
    -------
    Combined stdout+stderr output from Multiwfn.
    """
    if not multiwfn_available():
        raise RuntimeError(
            "Multiwfn binary not found in PATH. "
            "Download from http://sobereva.com/multiwfn/ and add to PATH."
        )

    input_path = Path(input_file)
    if not input_path.is_file():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    result = subprocess.run(
        ["Multiwfn", str(input_path)],
        input=commands,
        capture_output=True,
        text=True,
        timeout=timeout,
        cwd=str(input_path.parent),
    )
    return result.stdout + result.stderr


def convert_orca_gbw_to_molden(
    gbw_file: str | Path,
    *,
    orca_2mkl: str = "orca_2mkl",
) -> Path:
    """Convert ORCA .gbw to .molden.input using orca_2mkl.

    Returns the path to the generated .molden.input file.
    """
    gbw_path = Path(gbw_file)
    if not gbw_path.is_file():
        raise FileNotFoundError(f"GBW file not found: {gbw_path}")

    basename = gbw_path.stem
    subprocess.run(
        [orca_2mkl, basename, "-molden"],
        cwd=str(gbw_path.parent),
        capture_output=True,
        text=True,
        check=True,
        timeout=120,
    )
    molden_file = gbw_path.parent / f"{basename}.molden.input"
    if not molden_file.is_file():
        raise FileNotFoundError(
            f"orca_2mkl did not produce expected file: {molden_file}"
        )
    return molden_file


def bond_order_analysis(
    input_file: str | Path,
    *,
    method: str = "mayer",
) -> dict:
    """Run bond order analysis.

    Parameters
    ----------
    input_file : .molden / .wfn / .wfx file
    method : "mayer" (default), "wiberg", or "fuzzy"

    Returns
    -------
    dict with key "output" (raw text) and "bond_orders" (list of dicts)
    """
    method_map = {"mayer": "1", "wiberg": "2", "fuzzy": "3"}
    method_key = method_map.get(method.lower(), "1")

    # Menu: 9 (Bond order analysis) -> method -> q
    commands = f"9\n{method_key}\nq\n"
    output = _run_multiwfn(input_file, commands)

    bond_orders = []
    # Parse lines like: "   1(C )  --   2(C )    1.2345"
    pattern = re.compile(
        r"\s*(\d+)\(\s*(\w+)\s*\)\s*--\s*(\d+)\(\s*(\w+)\s*\)\s+([\d.]+)"
    )
    for match in pattern.finditer(output):
        bond_orders.append({
            "atom1_idx": int(match.group(1)),
            "atom1_elem": match.group(2),
            "atom2_idx": int(match.group(3)),
            "atom2_elem": match.group(4),
            "bond_order": float(match.group(5)),
        })

    return {"output": output, "bond_orders": bond_orders, "method": method}


def population_analysis(
    input_file: str | Path,
    *,
    method: str = "hirshfeld",
) -> dict:
    """Run population analysis (atomic charges).

    Parameters
    ----------
    method : "hirshfeld", "mulliken", "lowdin", or "becke"

    Returns
    -------
    dict with "output" and "charges" (list of dicts with atom_idx, elem, charge)
    """
    method_map = {
        "hirshfeld": "1",
        "mulliken": "5",
        "lowdin": "6",
        "becke": "10",
    }
    method_key = method_map.get(method.lower(), "1")

    # Menu: 7 (Population analysis) -> method -> q
    commands = f"7\n{method_key}\nq\n"
    output = _run_multiwfn(input_file, commands)

    charges = []
    # Parse charge output
    pattern = re.compile(
        r"\s*(\d+)\(\s*(\w+)\s*\)\s+.*?(-?[\d.]+)\s*$", re.MULTILINE
    )
    for match in pattern.finditer(output):
        charges.append({
            "atom_idx": int(match.group(1)),
            "elem": match.group(2),
            "charge": float(match.group(3)),
        })

    return {"output": output, "charges": charges, "method": method}


def orbital_composition(
    input_file: str | Path,
    *,
    orbital_indices: Optional[list[int]] = None,
) -> dict:
    """Run orbital composition analysis (HOMO, LUMO, etc.).

    Parameters
    ----------
    orbital_indices : list of 1-based orbital indices to analyse.
        If None, analyses HOMO and LUMO.

    Returns
    -------
    dict with "output" (raw text)
    """
    # Menu: 8 (Orbital composition) -> 1 (Mulliken) -> q
    commands = "8\n1\nq\n"
    output = _run_multiwfn(input_file, commands)
    return {"output": output}


def electrostatic_potential(
    input_file: str | Path,
    *,
    output_cube: Optional[str] = None,
) -> dict:
    """Calculate molecular electrostatic potential (ESP) on van der Waals surface.

    Returns
    -------
    dict with "output" and optionally path to generated cube file
    """
    # Menu: 12 (Quantitative ESP) -> 1 (on vdW surface) -> q
    commands = "12\n1\nq\n"
    output = _run_multiwfn(input_file, commands)
    return {"output": output}
