"""Python-level wrapper for Genarris integration in DELFIN workflows.

Provides convenience functions for crystal structure generation
without requiring users to manage Genarris internals directly.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from delfin.csp_tools import require_genarris
from delfin.common.logging import get_logger

logger = get_logger(__name__)


def run_genarris_cli(
    config_path: str | Path,
    *,
    np: int = 1,
    cwd: str | Path | None = None,
    mpirun: str = "mpirun",
    capture_output: bool = False,
    timeout: int | None = None,
) -> subprocess.CompletedProcess:
    """Run Genarris via its CLI (mpirun -np N gnrs --config FILE).

    Parameters
    ----------
    config_path : path to Genarris INI config file
    np : number of MPI processes
    cwd : working directory
    mpirun : mpirun executable (default: "mpirun")
    capture_output : capture stdout/stderr
    timeout : timeout in seconds

    Returns
    -------
    subprocess.CompletedProcess
    """
    from shutil import which

    gnrs_bin = which("gnrs")
    if gnrs_bin is None:
        from delfin.csp_tools import get_csp_tools_root

        candidate = get_csp_tools_root() / "bin" / "gnrs"
        if candidate.is_file():
            gnrs_bin = str(candidate)
        else:
            raise FileNotFoundError(
                "gnrs not found. Install Genarris: bash delfin/csp_tools/install_csp_tools.sh"
            )

    cmd = [mpirun, "-np", str(np), gnrs_bin, "--config", str(config_path)]
    logger.info("Running Genarris: %s", " ".join(cmd))

    return subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        capture_output=capture_output,
        text=True,
        check=False,
        timeout=timeout,
    )


def generate_random_structures(
    molecule_path: str | Path,
    *,
    num_structures: int = 100,
    space_groups: list[int] | None = None,
    Z: int = 4,
    output_dir: str | Path | None = None,
    np: int = 1,
    extra_config: dict[str, str] | None = None,
) -> Path:
    """Generate random molecular crystal structures using Genarris.

    This is a convenience wrapper that creates a Genarris config file
    and runs the generation.

    Parameters
    ----------
    molecule_path : path to molecular geometry (xyz, json, etc.)
    num_structures : number of structures to generate
    space_groups : list of space group numbers (default: common organic groups)
    Z : molecules per unit cell
    output_dir : where to write output (default: tempdir)
    np : MPI processes
    extra_config : additional Genarris config key-value pairs

    Returns
    -------
    Path to the output directory containing generated structures.
    """
    require_genarris("Random crystal structure generation")

    molecule_path = Path(molecule_path).resolve()
    if not molecule_path.is_file():
        raise FileNotFoundError(f"Molecule file not found: {molecule_path}")

    if output_dir is None:
        output_dir = Path(tempfile.mkdtemp(prefix="genarris_"))
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    if space_groups is None:
        # Common space groups for organic crystals
        space_groups = [2, 4, 14, 15, 19, 33, 61]

    sg_str = " ".join(str(sg) for sg in space_groups)

    config_lines = [
        "[Genarris_master]",
        "procedures = random_crystal_generation",
        "",
        "[random_crystal_generation]",
        f"molecule_path = {molecule_path}",
        f"num_structures = {num_structures}",
        f"space_groups = {sg_str}",
        f"Z = {Z}",
        f"output_dir = {output_dir}",
    ]

    if extra_config:
        for key, value in extra_config.items():
            config_lines.append(f"{key} = {value}")

    config_file = output_dir / "genarris_config.ini"
    config_file.write_text("\n".join(config_lines) + "\n", encoding="utf-8")

    logger.info(
        "Generating %d random crystal structures in %s", num_structures, output_dir
    )
    result = run_genarris_cli(config_file, np=np, cwd=output_dir, capture_output=True)

    if result.returncode != 0:
        logger.error("Genarris failed (rc=%d): %s", result.returncode, result.stderr)
        raise RuntimeError(
            f"Genarris crystal generation failed (rc={result.returncode}). "
            f"stderr: {result.stderr[:500]}"
        )

    logger.info("Crystal structure generation complete: %s", output_dir)
    return output_dir
