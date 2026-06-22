"""
NMR Spectrum Parser for ORCA output files.

Parses 1H NMR chemical shieldings, coordinates, and J-coupling constants
from ORCA NMR calculations (EPRNMR block).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from delfin.common.logging import get_logger

logger = get_logger(__name__)

# TMS reference shielding for 1H (TPSS/pcSseg-1 level, typical value).
# Users can override via CLI flag.
TMS_SHIELDING_1H = 31.4


@dataclass
class NMRShielding:
    """One nucleus from the CHEMICAL SHIELDING SUMMARY."""
    atom_index: int        # 0-based index in the coordinate block
    orca_index: int        # ORCA nucleus index (as printed in the output)
    element: str
    isotropic_ppm: float
    anisotropy_ppm: float

    @property
    def chemical_shift(self) -> float:
        """delta = sigma_ref - sigma_calc  (using module-level TMS ref)."""
        return TMS_SHIELDING_1H - self.isotropic_ppm


@dataclass
class Atom:
    """Cartesian atom from the coordinate block."""
    index: int   # 0-based
    element: str
    x: float
    y: float
    z: float


@dataclass
class NMRResult:
    """Complete parsed NMR result from an ORCA output file."""
    atoms: List[Atom] = field(default_factory=list)
    shieldings: List[NMRShielding] = field(default_factory=list)
    j_couplings: Dict[Tuple[int, int], float] = field(default_factory=dict)
    source_file: Optional[str] = None

    @property
    def h_shieldings(self) -> List[NMRShielding]:
        """Return only hydrogen shieldings."""
        return [s for s in self.shieldings if s.element == "H"]

    def get_atom(self, orca_index: int) -> Optional[Atom]:
        """Get atom by its ORCA nucleus index."""
        if 0 <= orca_index < len(self.atoms):
            return self.atoms[orca_index]
        return None

    def h_bonded_to(self, h_orca_index: int) -> Optional[Atom]:
        """Find the heavy atom that a hydrogen is bonded to (nearest non-H)."""
        h_atom = self.get_atom(h_orca_index)
        if h_atom is None:
            return None
        best_dist = 999.0
        best_atom = None
        for atom in self.atoms:
            if atom.index == h_atom.index or atom.element == "H":
                continue
            dx = atom.x - h_atom.x
            dy = atom.y - h_atom.y
            dz = atom.z - h_atom.z
            dist = (dx*dx + dy*dy + dz*dz) ** 0.5
            if dist < best_dist:
                best_dist = dist
                best_atom = atom
        return best_atom


def parse_nmr_orca(output_file: Path, tms_ref: float = TMS_SHIELDING_1H) -> NMRResult:
    """
    Parse NMR data from an ORCA output file.

    Extracts:
      - Cartesian coordinates (ANGSTROEM block)
      - Chemical shielding summary (isotropic + anisotropy)
      - Isotropic J-coupling matrix

    Parameters
    ----------
    output_file : Path
        Path to ORCA .out file with NMR calculation.
    tms_ref : float
        TMS reference shielding in ppm for delta calculation.

    Returns
    -------
    NMRResult
    """
    global TMS_SHIELDING_1H
    TMS_SHIELDING_1H = tms_ref

    text = Path(output_file).read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    result = NMRResult(source_file=str(output_file))
    result.atoms = _parse_coordinates(lines)
    result.shieldings = _parse_shielding_summary(lines, len(result.atoms))
    result.j_couplings = _parse_j_coupling_summary(lines)

    logger.info(
        "Parsed %d atoms, %d shieldings (%d H), %d J-couplings from %s",
        len(result.atoms), len(result.shieldings), len(result.h_shieldings),
        len(result.j_couplings), output_file.name,
    )
    return result


# ---------------------------------------------------------------------------
# Internal parsers
# ---------------------------------------------------------------------------

def _parse_coordinates(lines: List[str]) -> List[Atom]:
    """Parse the first CARTESIAN COORDINATES (ANGSTROEM) block."""
    atoms: List[Atom] = []
    in_block = False
    idx = 0
    for line in lines:
        if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
            in_block = True
            atoms.clear()
            idx = 0
            continue
        if in_block:
            if line.startswith("---"):
                continue
            stripped = line.strip()
            if not stripped:
                if atoms:
                    break  # end of block
                continue
            parts = stripped.split()
            if len(parts) == 4:
                try:
                    atoms.append(Atom(
                        index=idx,
                        element=parts[0],
                        x=float(parts[1]),
                        y=float(parts[2]),
                        z=float(parts[3]),
                    ))
                    idx += 1
                except ValueError:
                    pass
    return atoms


_SHIELDING_RE = re.compile(
    r'^\s*(\d+)\s+([A-Za-z]+)\s+([-\d.]+)\s+([-\d.]+)\s*$'
)


def _parse_shielding_summary(lines: List[str], num_atoms: int) -> List[NMRShielding]:
    """Parse the CHEMICAL SHIELDING SUMMARY table."""
    shieldings: List[NMRShielding] = []
    in_summary = False
    for line in lines:
        if "CHEMICAL SHIELDING SUMMARY" in line:
            in_summary = True
            shieldings.clear()
            continue
        if in_summary:
            m = _SHIELDING_RE.match(line)
            if m:
                orca_idx = int(m.group(1))
                shieldings.append(NMRShielding(
                    atom_index=orca_idx,
                    orca_index=orca_idx,
                    element=m.group(2),
                    isotropic_ppm=float(m.group(3)),
                    anisotropy_ppm=float(m.group(4)),
                ))
            elif shieldings and line.strip() == "":
                break  # end of table
    return shieldings


def _parse_j_coupling_summary(lines: List[str]) -> Dict[Tuple[int, int], float]:
    """Parse the SUMMARY OF ISOTROPIC COUPLING CONSTANTS J (Hz) matrix."""
    couplings: Dict[Tuple[int, int], float] = {}

    in_summary = False
    col_indices: List[int] = []

    for line in lines:
        if "SUMMARY OF ISOTROPIC COUPLING CONSTANTS" in line:
            in_summary = True
            couplings.clear()
            col_indices = []
            continue

        if not in_summary:
            continue

        stripped = line.strip()
        if not stripped or stripped.startswith("---"):
            continue

        # Header row: "  46 H    47 H    48 H ..."
        # or data row: "  46 H   0.000   7.636 ..."
        parts = stripped.split()

        # Detect header row: all tokens alternate between int and element
        if _is_header_row(parts):
            col_indices = [int(parts[i]) for i in range(0, len(parts), 2)]
            continue

        # Data row: first two tokens are "idx element", rest are floats
        if len(parts) >= 3 and parts[1].isalpha():
            try:
                row_idx = int(parts[0])
            except ValueError:
                continue
            values = parts[2:]
            for ci, val_str in enumerate(values):
                if ci >= len(col_indices):
                    break
                col_idx = col_indices[ci]
                try:
                    j_val = float(val_str)
                except ValueError:
                    continue
                if abs(j_val) > 0.001 and row_idx != col_idx:
                    key = (min(row_idx, col_idx), max(row_idx, col_idx))
                    couplings[key] = j_val

    return couplings


def _is_header_row(parts: List[str]) -> bool:
    """Check if a split line is a J-coupling header row like '46 H  47 H  48 H'."""
    if len(parts) < 2 or len(parts) % 2 != 0:
        return False
    for i in range(0, len(parts), 2):
        if not parts[i].isdigit():
            return False
        if not parts[i + 1].isalpha():
            return False
    return True
