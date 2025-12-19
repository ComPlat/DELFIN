"""Parse and process UV-Vis absorption spectra from ORCA output files."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import numpy as np

from delfin.common.logging import get_logger

logger = get_logger(__name__)


@dataclass
class Transition:
    """Single electronic transition with associated properties."""

    from_state: str  # e.g., "S0", "0-1A"
    to_state: str    # e.g., "T1", "1-3A"
    energy_ev: float
    energy_cm1: float
    wavelength_nm: float
    fosc: float  # Oscillator strength
    d2: float    # Dipole moment squared
    dx: float
    dy: float
    dz: float

    @property
    def readable_transition(self) -> str:
        """Convert ORCA notation to readable format (S0 -> T1, etc.)."""
        return f"{self._translate_state(self.from_state)} → {self._translate_state(self.to_state)}"

    @staticmethod
    def _translate_state(orca_state: str) -> str:
        """Translate ORCA state notation to physical notation.

        Examples:
            0-1A → S0
            1-1A → S1
            1-3A → T1
            2-3A → T2
        """
        match = re.match(r'(\d+)-([13])A', orca_state)
        if not match:
            return orca_state

        root_number = int(match.group(1))
        multiplicity = match.group(2)

        if multiplicity == '1':
            # Singlet
            return f"S{root_number}"
        elif multiplicity == '3':
            # Triplet: 1-3A is T1, 2-3A is T2, etc.
            return f"T{root_number}"
        else:
            return orca_state


def parse_absorption_spectrum(output_file: Path) -> List[Transition]:
    """Parse ABSORPTION SPECTRUM section from ORCA output file.

    Args:
        output_file: Path to ORCA .out file

    Returns:
        List of Transition objects
    """
    transitions = []

    try:
        with open(output_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except FileNotFoundError:
        logger.warning(f"Output file not found: {output_file}")
        return transitions

    # Find ALL absorption spectrum sections and take the LAST one (after final geometry)
    pattern = r'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS.*?\n-+\n.*?\n-+\n(.*?)(?:\n\n|\Z)'
    matches = re.findall(pattern, content, re.DOTALL)

    if not matches:
        logger.info(f"No absorption spectrum found in {output_file}")
        return transitions

    # Use the LAST spectrum (after geometry optimization converged)
    spectrum_text = matches[-1]
    if len(matches) > 1:
        logger.info(f"Found {len(matches)} absorption spectra in {output_file}, using the last one")

    # Parse each transition line
    # Format: 0-1A  ->  1-3A    2.600831   20977.1   476.7   0.000000000   0.00000   0.00000   0.00000   0.00000
    line_pattern = r'(\S+)\s+->\s+(\S+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)'

    for line in spectrum_text.strip().split('\n'):
        match = re.search(line_pattern, line)
        if match:
            from_state = match.group(1)
            to_state = match.group(2)
            energy_ev = float(match.group(3))
            energy_cm1 = float(match.group(4))
            wavelength_nm = float(match.group(5))
            fosc = float(match.group(6))
            d2 = float(match.group(7))
            dx = float(match.group(8))
            dy = float(match.group(9))
            dz = float(match.group(10))

            transitions.append(Transition(
                from_state=from_state,
                to_state=to_state,
                energy_ev=energy_ev,
                energy_cm1=energy_cm1,
                wavelength_nm=wavelength_nm,
                fosc=fosc,
                d2=d2,
                dx=dx,
                dy=dy,
                dz=dz
            ))

    logger.info(f"Parsed {len(transitions)} transitions from {output_file.name}")
    return transitions


def gaussian_broadening(
    transitions: List[Transition],
    wavelength_range: Tuple[float, float] = (200, 800),
    num_points: int = 1000,
    fwhm: float = 20.0
) -> Tuple[np.ndarray, np.ndarray]:
    """Apply Gaussian broadening to create a continuous spectrum.

    Args:
        transitions: List of transitions to broaden
        wavelength_range: (min_nm, max_nm) range for spectrum
        num_points: Number of points in output spectrum
        fwhm: Full width at half maximum (nm) for Gaussian broadening

    Returns:
        (wavelengths, intensities) arrays
    """
    wavelengths = np.linspace(wavelength_range[0], wavelength_range[1], num_points)
    intensities = np.zeros(num_points)

    # Convert FWHM to standard deviation
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))

    # Add Gaussian contribution from each transition
    for trans in transitions:
        if trans.fosc > 0:  # Only include transitions with non-zero oscillator strength
            gaussian = trans.fosc * np.exp(-0.5 * ((wavelengths - trans.wavelength_nm) / sigma) ** 2)
            intensities += gaussian

    # No normalization - keep actual fosc values for scientific accuracy
    return wavelengths, intensities


def filter_significant_transitions(
    transitions: List[Transition],
    min_fosc: float = 0.001
) -> List[Transition]:
    """Filter transitions by minimum oscillator strength.

    Args:
        transitions: List of all transitions
        min_fosc: Minimum oscillator strength threshold

    Returns:
        Filtered list of significant transitions
    """
    return [t for t in transitions if t.fosc >= min_fosc]
