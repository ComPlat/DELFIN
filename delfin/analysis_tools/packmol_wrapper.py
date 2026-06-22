"""Wrapper around Packmol for building initial configurations of molecular systems.

Packmol packs molecules into defined spatial regions, useful for:
  - Solvation boxes for MD simulations
  - Mixtures and interfaces
  - Random initial configurations for sampling
"""

from __future__ import annotations

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Optional


def _find_packmol() -> Optional[str]:
    """Return the Packmol executable path or None."""
    import shutil
    return shutil.which("packmol")


def solvate(
    solute_xyz: str | Path,
    solvent_xyz: str | Path,
    n_solvent: int = 50,
    box_side: float = 20.0,
    output_xyz: Optional[str | Path] = None,
    tolerance: float = 2.0,
    working_dir: Optional[str | Path] = None,
) -> Optional[Path]:
    """Place a solute in a box of solvent molecules.

    Parameters
    ----------
    solute_xyz : path to solute XYZ file
    solvent_xyz : path to solvent molecule XYZ file
    n_solvent : number of solvent molecules
    box_side : box side length in Angstrom
    output_xyz : output file path (default: working_dir/packed.xyz)
    tolerance : minimum distance between molecules (Angstrom)
    working_dir : directory for Packmol I/O (default: temp dir)

    Returns
    -------
    Path to the output XYZ file, or None on failure.
    """
    packmol = _find_packmol()
    if not packmol:
        logging.error("Packmol executable not found in PATH.")
        return None

    solute_xyz = Path(solute_xyz).resolve()
    solvent_xyz = Path(solvent_xyz).resolve()

    if not solute_xyz.exists():
        logging.error("Solute file not found: %s", solute_xyz)
        return None
    if not solvent_xyz.exists():
        logging.error("Solvent file not found: %s", solvent_xyz)
        return None

    if working_dir:
        wd = Path(working_dir)
        wd.mkdir(parents=True, exist_ok=True)
    else:
        wd = Path(tempfile.mkdtemp(prefix="packmol_"))

    if output_xyz:
        out = Path(output_xyz).resolve()
    else:
        out = wd / "packed.xyz"

    half = box_side / 2.0
    inp_text = f"""\
tolerance {tolerance}
filetype xyz
output {out}

structure {solute_xyz}
  number 1
  fixed 0.0 0.0 0.0 0.0 0.0 0.0
  center
end structure

structure {solvent_xyz}
  number {n_solvent}
  inside box -{half} -{half} -{half} {half} {half} {half}
end structure
"""
    inp_file = wd / "packmol.inp"
    inp_file.write_text(inp_text, encoding="utf-8")
    log_file = wd / "packmol.log"

    try:
        with log_file.open("w", encoding="utf-8") as lf:
            subprocess.run(
                [packmol],
                input=inp_text,
                stdout=lf,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=300,
                check=True,
            )
    except subprocess.CalledProcessError as exc:
        logging.error("Packmol failed (exit %d). See %s", exc.returncode, log_file)
        return None
    except subprocess.TimeoutExpired:
        logging.error("Packmol timed out after 300s.")
        return None

    if out.exists():
        logging.info("Packmol output: %s", out)
        return out

    logging.error("Packmol did not produce output file: %s", out)
    return None


def pack_molecules(
    structures: list[dict],
    output_xyz: str | Path,
    tolerance: float = 2.0,
    working_dir: Optional[str | Path] = None,
) -> Optional[Path]:
    """General-purpose packing of multiple molecule types.

    Parameters
    ----------
    structures : list of dicts, each with keys:
        - xyz: path to XYZ file
        - number: how many copies
        - constraint: spatial constraint, e.g.
          "inside box -10 -10 -10 10 10 10"
          "fixed 0 0 0 0 0 0" (with optional "center")
    output_xyz : output file path
    tolerance : minimum distance (Angstrom)
    working_dir : directory for Packmol I/O

    Returns
    -------
    Path to output file, or None on failure.
    """
    packmol = _find_packmol()
    if not packmol:
        logging.error("Packmol executable not found in PATH.")
        return None

    if working_dir:
        wd = Path(working_dir)
        wd.mkdir(parents=True, exist_ok=True)
    else:
        wd = Path(tempfile.mkdtemp(prefix="packmol_"))

    out = Path(output_xyz).resolve()
    lines = [f"tolerance {tolerance}", "filetype xyz", f"output {out}", ""]

    for s in structures:
        xyz = Path(s["xyz"]).resolve()
        if not xyz.exists():
            logging.error("Structure file not found: %s", xyz)
            return None
        lines.append(f"structure {xyz}")
        lines.append(f"  number {s['number']}")
        lines.append(f"  {s['constraint']}")
        lines.append("end structure")
        lines.append("")

    inp_text = "\n".join(lines)
    log_file = wd / "packmol.log"

    try:
        with log_file.open("w", encoding="utf-8") as lf:
            subprocess.run(
                [packmol],
                input=inp_text,
                stdout=lf,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=300,
                check=True,
            )
    except subprocess.CalledProcessError as exc:
        logging.error("Packmol failed (exit %d). See %s", exc.returncode, log_file)
        return None
    except subprocess.TimeoutExpired:
        logging.error("Packmol timed out after 300s.")
        return None

    if out.exists():
        logging.info("Packmol output: %s", out)
        return out

    logging.error("Packmol did not produce output file: %s", out)
    return None
