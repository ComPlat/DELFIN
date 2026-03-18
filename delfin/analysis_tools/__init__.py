"""Analysis tools for DELFIN.

Provides unified access to:
  - Multiwfn   — wavefunction analysis (orbital composition, bond orders, ESP, ELF, population analysis)
  - CENSO      — conformer ensemble sorting and refinement with thermochemistry
  - morfeus    — steric descriptors (buried volume, cone angle, bite angle, Sterimol)
  - cclib      — computational chemistry log file parser (Gaussian, ORCA, GAMESS, …)
  - nglview    — interactive 3D molecular viewer for Jupyter notebooks
  - Packmol    — initial configurations for molecular dynamics simulations
"""

from __future__ import annotations

import importlib
import os
import shutil
from pathlib import Path
from typing import Optional


def get_analysis_tools_root() -> Path:
    """Return the root directory of the analysis_tools package."""
    env_root = os.environ.get("DELFIN_ANALYSIS_TOOLS_ROOT")
    if env_root:
        return Path(env_root).expanduser().resolve()
    return Path(__file__).resolve().parent


# ── Multiwfn ──────────────────────────────────────────────────────────

def multiwfn_available() -> bool:
    """Check whether the Multiwfn binary is in PATH."""
    return shutil.which("Multiwfn") is not None


def get_multiwfn_path() -> Optional[str]:
    """Return the full path to the Multiwfn binary, or None."""
    return shutil.which("Multiwfn")


def get_multiwfn_version() -> Optional[str]:
    """Return the Multiwfn version string, or None if not available."""
    if not multiwfn_available():
        return None
    try:
        import subprocess

        result = subprocess.run(
            ["Multiwfn", "--version"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        # Multiwfn prints version info to stdout on startup
        output = result.stdout + result.stderr
        for line in output.splitlines():
            if "Version" in line or "Multiwfn" in line:
                return line.strip()
        return "installed"
    except Exception:
        return "installed"


# ── CENSO ─────────────────────────────────────────────────────────────

def censo_available() -> bool:
    """Check whether CENSO is importable or available as CLI."""
    try:
        if importlib.util.find_spec("censo") is not None:
            return True
    except (ModuleNotFoundError, ValueError):
        pass
    return shutil.which("censo") is not None


def get_censo_version() -> Optional[str]:
    """Return the CENSO version string, or None."""
    if not censo_available():
        return None
    try:
        from importlib.metadata import version

        return version("censo")
    except Exception:
        pass
    try:
        import subprocess

        result = subprocess.run(
            ["censo", "--version"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        output = (result.stdout + result.stderr).strip()
        if output:
            return output.splitlines()[0]
    except Exception:
        pass
    return "installed"


# ── morfeus ───────────────────────────────────────────────────────────

def morfeus_available() -> bool:
    """Check whether morfeus-ml is importable."""
    try:
        return importlib.util.find_spec("morfeus") is not None
    except (ModuleNotFoundError, ValueError):
        return False


def get_morfeus_version() -> Optional[str]:
    """Return the morfeus version string, or None."""
    if not morfeus_available():
        return None
    try:
        from importlib.metadata import version

        return version("morfeus-ml")
    except Exception:
        pass
    try:
        import morfeus

        return getattr(morfeus, "__version__", "installed")
    except Exception:
        return "installed"


# ── cclib ─────────────────────────────────────────────────────────────

def cclib_available() -> bool:
    """Check whether cclib is importable."""
    try:
        return importlib.util.find_spec("cclib") is not None
    except (ModuleNotFoundError, ValueError):
        return False


def get_cclib_version() -> Optional[str]:
    """Return the cclib version string, or None."""
    if not cclib_available():
        return None
    try:
        from importlib.metadata import version
        return version("cclib")
    except Exception:
        return "installed"


# ── nglview ──────────────────────────────────────────────────────────

def nglview_available() -> bool:
    """Check whether nglview is importable."""
    try:
        return importlib.util.find_spec("nglview") is not None
    except (ModuleNotFoundError, ValueError):
        return False


def get_nglview_version() -> Optional[str]:
    """Return the nglview version string, or None."""
    if not nglview_available():
        return None
    try:
        from importlib.metadata import version
        return version("nglview")
    except Exception:
        return "installed"


# ── Packmol ──────────────────────────────────────────────────────────

def packmol_available() -> bool:
    """Check whether the Packmol binary is in PATH."""
    return shutil.which("packmol") is not None


def get_packmol_version() -> Optional[str]:
    """Return the Packmol version string, or None."""
    if not packmol_available():
        return None
    try:
        import re
        import subprocess
        result = subprocess.run(
            ["packmol"],
            input="",
            capture_output=True,
            text=True,
            timeout=5,
        )
        output = result.stdout + result.stderr
        for line in output.splitlines():
            if "Version" in line:
                # Extract version number like "20.15.3"
                m = re.search(r"(\d+\.\d+(?:\.\d+)?)", line)
                if m:
                    return m.group(1)
                return line.strip()
        return "installed"
    except Exception:
        return "installed"


# ── summary ───────────────────────────────────────────────────────────

def available_tools() -> list[str]:
    """Return a list of available analysis tool names."""
    tools = []
    if multiwfn_available():
        tools.append("multiwfn")
    if censo_available():
        tools.append("censo")
    if morfeus_available():
        tools.append("morfeus")
    if cclib_available():
        tools.append("cclib")
    if nglview_available():
        tools.append("nglview")
    if packmol_available():
        tools.append("packmol")
    return tools


def collect_analysis_summary() -> dict:
    """Return a dict summarising analysis tool status for dashboard display."""
    tools_info = []
    for label, avail_fn, ver_fn, description in [
        ("Multiwfn", multiwfn_available, get_multiwfn_version,
         "Wavefunction analysis: orbitals, bond orders, ESP, ELF, population analysis"),
        ("CENSO", censo_available, get_censo_version,
         "Conformer ensemble sorting with DFT refinement and thermochemistry"),
        ("morfeus", morfeus_available, get_morfeus_version,
         "Steric descriptors: buried volume, cone angle, bite angle, Sterimol"),
        ("cclib", cclib_available, get_cclib_version,
         "Log file parser for ORCA, Gaussian, GAMESS, Turbomole, NWChem and more"),
        ("nglview", nglview_available, get_nglview_version,
         "Interactive 3D molecular viewer for Jupyter/Voilà notebooks"),
        ("Packmol", packmol_available, get_packmol_version,
         "Packing molecules for MD simulation initial configurations"),
    ]:
        ok = avail_fn()
        tools_info.append({
            "name": label,
            "installed": ok,
            "version": ver_fn() or "" if ok else "",
            "description": description,
        })
    return {"tools": tools_info}
