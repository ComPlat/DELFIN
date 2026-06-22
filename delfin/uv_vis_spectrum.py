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
    excitations: List[dict] = None  # List of orbital excitations with weights
    homo_number: int = None  # HOMO orbital number for this transition

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

    # Convert ORCA notation to S/T notation
    # ORCA numbers states by energy (all symmetries together), but we need to count separately:
    # - First X-1A (singlet, X>0) = S1, second = S2, etc.
    # - First X-3A (triplet) = T1, second = T2, etc.
    singlet_counter = 0
    triplet_counter = 0
    state_mapping = {"0-1A": "S0"}  # Ground state

    for trans in transitions:
        # Map to_state if not already mapped
        if trans.to_state not in state_mapping:
            if "-1A" in trans.to_state and trans.to_state != "0-1A":
                singlet_counter += 1
                state_mapping[trans.to_state] = f"S{singlet_counter}"
            elif "-3A" in trans.to_state:
                triplet_counter += 1
                state_mapping[trans.to_state] = f"T{triplet_counter}"
            else:
                # Other symmetries (keep ORCA notation)
                state_mapping[trans.to_state] = trans.to_state

        # Map from_state if not already mapped
        if trans.from_state not in state_mapping:
            state_mapping[trans.from_state] = trans.from_state

    # Apply mapping to all transitions
    for trans in transitions:
        trans.from_state = state_mapping.get(trans.from_state, trans.from_state)
        trans.to_state = state_mapping.get(trans.to_state, trans.to_state)

    # Parse orbital excitations from STATE sections
    _add_orbital_excitations(content, transitions)

    logger.info(f"Parsed {len(transitions)} transitions from {output_file.name}")
    return transitions


def _add_orbital_excitations(content: str, transitions: List[Transition]) -> None:
    """Parse STATE sections and add orbital excitation information to transitions.

    Args:
        content: Full ORCA output file content
        transitions: List of Transition objects to augment with excitation data
    """
    # Find HOMO orbital number from the last ORBITAL ENERGIES section
    homo_number = _find_homo_number(content)
    if homo_number is None:
        return  # Cannot parse excitations without HOMO number

    # Parse all STATE sections
    # Format: STATE  1:  E=   0.091053 au      2.478 eV    19983.7 cm**-1 <S**2> =   2.000000 Mult 3
    #             39a ->  42a  :     0.849876
    state_pattern = r'STATE\s+(\d+):\s+E=\s+[\d.]+\s+au\s+([\d.]+)\s+eV.*?\n((?:\s+\d+[ab]\s+->\s+\d+[ab]\s+:\s+[\d.]+\s*\n)+)'

    state_matches = re.findall(state_pattern, content, re.MULTILINE)

    # Build map from energy to excitations
    energy_to_excitations = {}
    for state_num, energy_ev_str, excitations_text in state_matches:
        energy_ev = float(energy_ev_str)

        # Parse individual excitations
        excitations = []
        exc_pattern = r'(\d+)([ab])\s+->\s+(\d+)([ab])\s+:\s+([\d.]+)'
        for exc_match in re.finditer(exc_pattern, excitations_text):
            from_orb = int(exc_match.group(1))
            from_spin = exc_match.group(2)
            to_orb = int(exc_match.group(3))
            to_spin = exc_match.group(4)
            weight = float(exc_match.group(5))

            excitations.append({
                "from_orbital": from_orb,
                "from_spin": from_spin,
                "to_orbital": to_orb,
                "to_spin": to_spin,
                "weight": weight
            })

        if excitations:
            energy_to_excitations[energy_ev] = excitations

    # Match transitions to STATE sections by energy (with tolerance)
    for trans in transitions:
        # Find closest matching energy in STATE sections
        closest_energy = None
        min_diff = float('inf')
        for state_energy in energy_to_excitations.keys():
            diff = abs(trans.energy_ev - state_energy)
            if diff < min_diff and diff < 0.01:  # 0.01 eV tolerance
                min_diff = diff
                closest_energy = state_energy

        if closest_energy is not None:
            trans.excitations = energy_to_excitations[closest_energy]
            trans.homo_number = homo_number


def _find_homo_number(content: str) -> int:
    """Find HOMO orbital number from ORBITAL ENERGIES section.

    Args:
        content: Full ORCA output file content

    Returns:
        HOMO orbital number or None if not found
    """
    # Find last ORBITAL ENERGIES section (after final optimization)
    orbital_sections = list(re.finditer(r'ORBITAL ENERGIES\s*\n-+', content))
    if not orbital_sections:
        return None

    last_section = orbital_sections[-1]
    orbital_text = content[last_section.end():][:500000]  # Take next 500k chars

    # Find HOMO (last orbital with OCC = 2.0000 or 1.0000)
    homo_num = None
    for line in orbital_text.split('\n'):
        orbital_match = re.match(r'\s*(\d+)\s+(2\.0000|1\.0000)\s+', line)
        if orbital_match:
            homo_num = int(orbital_match.group(1))

    return homo_num


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
