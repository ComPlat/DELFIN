"""Wrapper around cclib for parsing computational chemistry output files.

cclib supports ORCA, Gaussian, GAMESS, Turbomole, NWChem, Q-Chem, and many more.
This wrapper provides convenient high-level functions for the most common use cases
in DELFIN workflows.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional


def parse_output(filepath: str | Path) -> Optional[dict]:
    """Parse a computational chemistry output file and return key results.

    Returns a dict with available fields:
      - scfenergies: list of SCF energies (eV)
      - moenergies: MO energies per spin channel
      - homos: HOMO indices per spin channel
      - charge: formal charge
      - mult: spin multiplicity
      - natom: number of atoms
      - atomcoords: final geometry (Nx3 array)
      - atomnos: atomic numbers
      - vibfreqs: vibrational frequencies (cm⁻¹)
      - enthalpy: enthalpy (hartree)
      - freeenergy: free energy (hartree)
    """
    try:
        import cclib
    except ImportError:
        logging.error("cclib is not installed. Run: pip install cclib")
        return None

    filepath = Path(filepath)
    if not filepath.exists():
        logging.error("File not found: %s", filepath)
        return None

    try:
        data = cclib.io.ccread(str(filepath))
    except Exception as exc:
        logging.error("cclib failed to parse %s: %s", filepath.name, exc)
        return None

    if data is None:
        logging.warning("cclib returned no data for %s", filepath.name)
        return None

    result = {}
    for attr in (
        "scfenergies", "moenergies", "homos", "charge", "mult",
        "natom", "atomcoords", "atomnos", "vibfreqs",
        "enthalpy", "freeenergy",
    ):
        val = getattr(data, attr, None)
        if val is not None:
            result[attr] = val

    return result


def extract_energies(filepath: str | Path) -> Optional[dict]:
    """Extract final SCF energy, enthalpy, and free energy from an output file.

    Returns dict with keys: scf_ev, scf_hartree, enthalpy, freeenergy (all floats).
    """
    parsed = parse_output(filepath)
    if not parsed:
        return None

    result = {}
    if "scfenergies" in parsed:
        scf_ev = float(parsed["scfenergies"][-1])
        result["scf_ev"] = scf_ev
        result["scf_hartree"] = scf_ev / 27.211386245988
    if "enthalpy" in parsed:
        result["enthalpy"] = float(parsed["enthalpy"])
    if "freeenergy" in parsed:
        result["freeenergy"] = float(parsed["freeenergy"])
    return result if result else None


def extract_vibrations(filepath: str | Path) -> Optional[dict]:
    """Extract vibrational frequencies from an output file.

    Returns dict with: freqs (list[float]), n_imaginary (int).
    """
    parsed = parse_output(filepath)
    if not parsed or "vibfreqs" not in parsed:
        return None

    freqs = [float(f) for f in parsed["vibfreqs"]]
    return {
        "freqs": freqs,
        "n_imaginary": sum(1 for f in freqs if f < 0),
    }
