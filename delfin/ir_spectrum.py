"""
IR Spectrum Parser for ORCA output files.

Parses infrared vibrational spectra from ORCA frequency calculations.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import List

from delfin.common.logging import get_logger

logger = get_logger(__name__)


@dataclass
class IRMode:
    """Represents a single vibrational mode in an IR spectrum."""
    mode_number: int
    frequency_cm1: float
    epsilon: float  # L/(mol*cm)
    intensity_km_mol: float
    t_squared: float
    tx: float
    ty: float
    tz: float

    @property
    def transmittance_percent(self) -> float:
        """Calculate transmittance percentage (for plotting)."""
        # Simple Beer-Lambert approximation: higher intensity = lower transmittance
        # For typical plotting, we invert the intensity
        # This is a simplified model; adjust as needed
        if self.intensity_km_mol <= 0:
            return 100.0
        # Scale intensity to reasonable transmittance range
        # Typical strong peaks have intensity 100-500 km/mol -> ~10-50% transmittance
        # Weak peaks < 10 km/mol -> ~90-99% transmittance
        max_absorption = min(self.intensity_km_mol / 500.0, 1.0)
        return 100.0 * (1.0 - max_absorption)


def parse_ir_spectrum(output_file: Path) -> List[IRMode]:
    """
    Parse IR spectrum from ORCA frequency output file.

    Args:
        output_file: Path to ORCA .out file containing frequency calculation

    Returns:
        List of IRMode objects representing the vibrational modes

    Example output format:
        -----------
        IR SPECTRUM
        -----------

         Mode   freq       eps      Int      T**2         TX        TY        TZ
               cm**-1   L/(mol*cm) km/mol    a.u.
        ----------------------------------------------------------------------------
          6:     89.55   0.000058    0.30  0.000203  (-0.002326 -0.009616 -0.010276)
          7:    113.29   0.000190    0.96  0.000522  ( 0.012132 -0.002599 -0.019194)
    """
    if not output_file.exists():
        logger.error(f"Output file not found: {output_file}")
        return []

    try:
        with open(output_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        # Find IR SPECTRUM section
        ir_section_match = re.search(
            r'-+\s*\n\s*IR SPECTRUM\s*\n\s*-+\s*\n.*?'
            r'Mode\s+freq\s+eps\s+Int\s+T\*\*2\s+TX\s+TY\s+TZ.*?\n\s*-+\s*\n'
            r'(.*?)(?:\n\s*\n|\Z)',
            content,
            re.DOTALL | re.IGNORECASE
        )

        if not ir_section_match:
            logger.warning(f"No IR SPECTRUM section found in {output_file}")
            return []

        ir_data_block = ir_section_match.group(1)
        modes = []

        # Parse each mode line
        # Format: Mode: freq eps Int T**2 (TX TY TZ)
        mode_pattern = re.compile(
            r'^\s*(\d+):\s+'  # Mode number
            r'([-\d.]+)\s+'   # Frequency
            r'([-\d.eE+]+)\s+'  # Epsilon
            r'([-\d.eE+]+)\s+'  # Intensity
            r'([-\d.eE+]+)\s+'  # T**2
            r'\(\s*([-\d.eE+]+)\s+'  # TX
            r'([-\d.eE+]+)\s+'       # TY
            r'([-\d.eE+]+)\s*\)',    # TZ
            re.MULTILINE
        )

        for match in mode_pattern.finditer(ir_data_block):
            mode = IRMode(
                mode_number=int(match.group(1)),
                frequency_cm1=float(match.group(2)),
                epsilon=float(match.group(3)),
                intensity_km_mol=float(match.group(4)),
                t_squared=float(match.group(5)),
                tx=float(match.group(6)),
                ty=float(match.group(7)),
                tz=float(match.group(8))
            )
            modes.append(mode)

        logger.info(f"Parsed {len(modes)} IR modes from {output_file.name}")
        return modes

    except Exception as e:
        logger.error(f"Error parsing IR spectrum from {output_file}: {e}", exc_info=True)
        return []


def calculate_spectrum_range(modes: List[IRMode]) -> tuple[float, float]:
    """
    Calculate appropriate x-axis range for IR spectrum plot.

    Args:
        modes: List of IR modes

    Returns:
        Tuple of (min_wavenumber, max_wavenumber) with padding
    """
    if not modes:
        return 400.0, 4000.0  # Default IR range

    frequencies = [m.frequency_cm1 for m in modes]
    min_freq = min(frequencies)
    max_freq = max(frequencies)

    # Add 10% padding
    padding = (max_freq - min_freq) * 0.1

    # Round to nearest 100
    min_range = max(400, int((min_freq - padding) / 100) * 100)
    max_range = min(4000, int((max_freq + padding + 99) / 100) * 100)

    return min_range, max_range
