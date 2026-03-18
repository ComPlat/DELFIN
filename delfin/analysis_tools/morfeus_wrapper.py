"""Python wrapper for morfeus steric descriptor calculations.

morfeus computes steric descriptors commonly used for characterising
transition-metal complexes and ligands:

  - Buried volume (%Vbur)
  - Cone angle
  - Bite angle (bidentate ligands)
  - Sterimol parameters (L, B1, B5)
  - Solid angle

Install::

    pip install morfeus-ml

Usage::

    from delfin.analysis_tools.morfeus_wrapper import (
        buried_volume, cone_angle, bite_angle, sterimol, full_steric_report,
    )

    # From XYZ file — specify metal center atom index (0-based)
    bv = buried_volume("complex.xyz", metal_index=0)
    print(f"Buried volume: {bv['percent_buried_volume']:.1f}%")

    ca = cone_angle("complex.xyz", metal_index=0)
    print(f"Cone angle: {ca['cone_angle']:.1f}°")

    # Full report for a metal complex
    report = full_steric_report("complex.xyz", metal_index=0)
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from . import morfeus_available


def _require_morfeus():
    if not morfeus_available():
        raise ImportError(
            "morfeus-ml is not installed. Install with:\n"
            "  pip install morfeus-ml"
        )


def _read_structure(xyz_file: str | Path) -> tuple:
    """Read XYZ file and return (elements, coordinates)."""
    _require_morfeus()
    from morfeus import read_xyz

    return read_xyz(str(xyz_file))


def buried_volume(
    xyz_file: str | Path,
    metal_index: int,
    *,
    radius: float = 3.5,
    include_hydrogens: bool = True,
    excluded_atoms: Optional[list[int]] = None,
) -> dict:
    """Calculate percent buried volume (%Vbur) around a metal center.

    Parameters
    ----------
    xyz_file : path to XYZ file
    metal_index : 0-based index of the metal center atom
    radius : sphere radius in Angstrom (default 3.5)
    include_hydrogens : whether to include H atoms
    excluded_atoms : list of 0-based atom indices to exclude

    Returns
    -------
    dict with keys: percent_buried_volume, free_volume, buried_volume,
                    total_volume, metal_index, radius
    """
    _require_morfeus()
    from morfeus import BuriedVolume

    elements, coordinates = _read_structure(xyz_file)
    # morfeus uses 1-based indexing
    metal_idx_1 = metal_index + 1
    excluded_1 = [i + 1 for i in (excluded_atoms or [])]

    bv = BuriedVolume(
        elements,
        coordinates,
        metal_idx_1,
        radius=radius,
        include_hs=include_hydrogens,
        excluded_atoms=excluded_1 or None,
    )

    return {
        "percent_buried_volume": bv.percent_buried_volume,
        "free_volume": bv.free_volume,
        "buried_volume": bv.buried_volume,
        "total_volume": 4 / 3 * 3.14159265 * radius**3,
        "metal_index": metal_index,
        "radius": radius,
    }


def cone_angle(
    xyz_file: str | Path,
    metal_index: int,
) -> dict:
    """Calculate exact cone angle of a ligand around a metal center.

    Parameters
    ----------
    xyz_file : path to XYZ file
    metal_index : 0-based index of the metal center atom

    Returns
    -------
    dict with keys: cone_angle, metal_index
    """
    _require_morfeus()
    from morfeus import ConeAngle

    elements, coordinates = _read_structure(xyz_file)
    metal_idx_1 = metal_index + 1

    ca = ConeAngle(elements, coordinates, metal_idx_1)

    return {
        "cone_angle": ca.cone_angle,
        "metal_index": metal_index,
    }


def bite_angle(
    xyz_file: str | Path,
    metal_index: int,
    donor_indices: list[int],
) -> dict:
    """Calculate bite angle for a bidentate ligand.

    Parameters
    ----------
    xyz_file : path to XYZ file
    metal_index : 0-based index of the metal center
    donor_indices : 0-based indices of the two donor atoms

    Returns
    -------
    dict with keys: bite_angle, metal_index, donor_indices
    """
    if len(donor_indices) != 2:
        raise ValueError("bite_angle requires exactly 2 donor atom indices")

    _require_morfeus()

    import numpy as np

    elements, coordinates = _read_structure(xyz_file)

    # Simple geometric bite angle: angle donor1-metal-donor2
    m = coordinates[metal_index]
    d1 = coordinates[donor_indices[0]]
    d2 = coordinates[donor_indices[1]]

    v1 = d1 - m
    v2 = d2 - m
    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle_deg = float(np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0))))

    return {
        "bite_angle": angle_deg,
        "metal_index": metal_index,
        "donor_indices": donor_indices,
    }


def sterimol(
    xyz_file: str | Path,
    atom1_index: int,
    atom2_index: int,
) -> dict:
    """Calculate Sterimol parameters (L, B1, B5) along an axis.

    Parameters
    ----------
    xyz_file : path to XYZ file
    atom1_index : 0-based index of the first atom defining the axis
    atom2_index : 0-based index of the second atom defining the axis

    Returns
    -------
    dict with keys: L, B1, B5 (in Angstrom)
    """
    _require_morfeus()
    from morfeus import Sterimol

    elements, coordinates = _read_structure(xyz_file)
    # morfeus uses 1-based
    s = Sterimol(elements, coordinates, atom1_index + 1, atom2_index + 1)

    return {
        "L": s.L_value,
        "B1": s.B_1_value,
        "B5": s.B_5_value,
        "atom1_index": atom1_index,
        "atom2_index": atom2_index,
    }


def full_steric_report(
    xyz_file: str | Path,
    metal_index: int,
    *,
    donor_indices: Optional[list[int]] = None,
    radius: float = 3.5,
) -> dict:
    """Generate a full steric report for a metal complex.

    Computes buried volume, cone angle, and optionally bite angle
    (if donor_indices are provided).

    Parameters
    ----------
    xyz_file : path to XYZ file
    metal_index : 0-based index of the metal center
    donor_indices : optional 0-based indices of donor atoms for bite angle
    radius : sphere radius for buried volume

    Returns
    -------
    dict with keys: buried_volume, cone_angle, bite_angle (if applicable)
    """
    report = {
        "xyz_file": str(xyz_file),
        "metal_index": metal_index,
    }

    try:
        report["buried_volume"] = buried_volume(
            xyz_file, metal_index, radius=radius
        )
    except Exception as exc:
        report["buried_volume"] = {"error": str(exc)}

    try:
        report["cone_angle"] = cone_angle(xyz_file, metal_index)
    except Exception as exc:
        report["cone_angle"] = {"error": str(exc)}

    if donor_indices and len(donor_indices) == 2:
        try:
            report["bite_angle"] = bite_angle(
                xyz_file, metal_index, donor_indices
            )
        except Exception as exc:
            report["bite_angle"] = {"error": str(exc)}

    return report
