# energies.py
import re
from typing import Optional, Sequence, Dict, Any, Tuple

from delfin.common.logging import get_logger

logger = get_logger(__name__)

FLOAT_RE = r'([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)'

def _read_text(path: str) -> Optional[str]:
    """Read text content from file with error handling.

    Args:
        path: Path to the file to read

    Returns:
        File content as string, or None if file not found
    """
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            return f.read()
    except FileNotFoundError:
        logger.error(f"File {path} not found.")
        return None

def _search_last_float(path: str, patterns: Sequence[str]) -> Optional[float]:
    """Search for the last occurrence of a float matching any pattern in file.

    Args:
        path: Path to the file to search
        patterns: List of regex patterns to search for

    Returns:
        Last matching float value, or None if no match found
    """
    text = _read_text(path)
    if text is None:
        return None
    for pat in patterns:
        matches = re.findall(pat, text, flags=re.IGNORECASE | re.DOTALL)
        if matches:
            try:
                return float(matches[-1])
            except ValueError:
                logger.error(f"Could not convert '{matches[-1]}' to float from {path}.")
                return None
    return None

def find_gibbs_energy(filename: str) -> Optional[float]:
    """Extract Gibbs free energy from ORCA output file.

    Searches for 'Final Gibbs free energy' patterns in the output.

    Args:
        filename: Path to ORCA output file

    Returns:
        Gibbs free energy in Hartree, or None if not found
    """

    patterns = [
        rf"Final\s+Gibbs\s+free\s+energy.*?{FLOAT_RE}\s*(?:E[hH]|a\.u\.)?",
        rf"Gibbs\s+free\s+energy.*?{FLOAT_RE}\s*(?:E[hH]|a\.u\.)?",
    ]
    return _search_last_float(filename, patterns)

def find_ZPE(filename: str) -> Optional[float]:
    """Extract zero-point energy from ORCA output file.

    Searches for 'Zero point energy' or 'Zero point correction' patterns.

    Args:
        filename: Path to ORCA output file

    Returns:
        Zero-point energy in Hartree, or None if not found
    """
    patterns = [
        rf"Zero[\s-]?point\s+energy.*?{FLOAT_RE}\s*(?:E[hH]|a\.u\.)?",
        rf"Zero[\s-]?point\s+correction.*?{FLOAT_RE}\s*(?:E[hH]|a\.u\.)?",
    ]
    return _search_last_float(filename, patterns)

def find_electronic_energy(filename: str) -> Optional[float]:
    """Extract final single point electronic energy from ORCA output file.

    Searches for 'FINAL SINGLE POINT ENERGY' or equivalent patterns.

    Args:
        filename: Path to ORCA output file

    Returns:
        Electronic energy in Hartree, or None if not found
    """
    patterns = [
        rf"FINAL\s+SINGLE\s+POINT\s+ENERGY\s+{FLOAT_RE}",
        rf"Total\s+Energy\s*:\s*{FLOAT_RE}",
        rf"Electronic\s+energy.*?{FLOAT_RE}\s*(?:E[hH]|a\.u\.)?",
    ]
    return _search_last_float(filename, patterns)


def find_state1_ohne_SOC(filename3: str) -> Optional[float]:
    """Extract S1 excited state energy without spin-orbit coupling.

    Parses TD-DFT results from 'CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' section.

    Args:
        filename3: Path to ORCA output file with TD-DFT results

    Returns:
        S1 state energy in eV, or None if not found
    """
    search_text = "CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
    last_value = None
    try:
        with open(filename3, 'r') as file:
            lines = file.readlines()
            for i, line in enumerate(lines):
                if search_text in line:
                    target_line_index = i + 5  # Fifth line after the header
                    if target_line_index < len(lines):
                        target_line = lines[target_line_index]
                        parts = target_line[20:].split()
                        try:
                            last_value = float(parts[0])  # First numeric element
                        except ValueError:
                            logger.error(f"Could not convert '{parts[0]}' to a float.")
                            continue
    except FileNotFoundError:
        logger.error(f"File {filename3} not found.")
        return None
    return last_value

def find_state3_ohne_SOC(filename3: str) -> Optional[float]:
    """Extract S3 excited state energy without spin-orbit coupling.

    Parses TD-DFT results from 'CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' section.

    Args:
        filename3: Path to ORCA output file with TD-DFT results

    Returns:
        S3 state energy in eV, or None if not found
    """
    search_text = "CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
    last_value = None
    try:
        with open(filename3, 'r') as file:
            lines = file.readlines()
            for i, line in enumerate(lines):
                if search_text in line:
                    target_line_index = i + 7  # Seventh line after the header
                    if target_line_index < len(lines):
                        target_line = lines[target_line_index]
                        parts = target_line[20:].split()
                        try:
                            last_value = float(parts[0])  # First numeric element
                        except ValueError:
                            logger.error(f"Could not convert '{parts[0]}' to a float.")
                            continue
    except FileNotFoundError:
        logger.error(f"File {filename3} not found.")
        return None
    return last_value


def find_state1_mit_SOC(filename3: str) -> Optional[float]:
    """Extract S1 excited state energy with spin-orbit coupling corrections.

    Parses SOC-corrected TD-DFT results from ORCA output.

    Args:
        filename3: Path to ORCA output file with SOC-TD-DFT results

    Returns:
        SOC-corrected S1 state energy in eV, or None if not found
    """
    search_text = "SOC CORRECTED CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
    try:
        with open(filename3, 'r') as file:
            lines = file.readlines()
            for i, line in enumerate(lines):
                if search_text in line:
                    target_line_index = i + 5  # Fifth line after the header
                    if target_line_index < len(lines):
                        target_line = lines[target_line_index]
                        parts = target_line[20:].split()
                        try:
                            return float(parts[0])  # First numeric value
                        except ValueError:
                            logger.error(f"Could not convert '{parts[0]}' to a float.")
                            return None
    except FileNotFoundError:
        logger.error(f"File {filename3} not found.")
    return None

def find_state3_mit_SOC(filename3: str) -> Optional[float]:
    """Extract S3 excited state energy with spin-orbit coupling corrections.

    Parses SOC-corrected TD-DFT results from ORCA output.

    Args:
        filename3: Path to ORCA output file with SOC-TD-DFT results

    Returns:
        SOC-corrected S3 state energy in eV, or None if not found
    """
    search_text = "SOC CORRECTED CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
    try:
        with open(filename3, 'r') as file:
            lines = file.readlines()
            for i, line in enumerate(lines):
                if search_text in line:
                    target_line_index = i + 7  # Seventh line after the header
                    if target_line_index < len(lines):
                        target_line = lines[target_line_index]
                        parts = target_line[20:].split()
                        try:
                            return float(parts[0])  # First numeric value
                        except ValueError:
                            logger.error(f"Could not convert '{parts[0]}' to a float.")
                            return None
    except FileNotFoundError:
        logger.error(f"File {filename3} not found.")
    return None

def check_and_execute_SOC(filename3: str, config: Dict[str, Any]) -> Tuple[Optional[float], Optional[float]]:
    """Extract excited state energies with or without SOC based on configuration.

    Chooses between SOC-corrected or regular TD-DFT results based on config['DOSOC'].

    Args:
        filename3: Path to ORCA output file with TD-DFT results
        config: Configuration dictionary containing 'DOSOC' setting

    Returns:
        Tuple of (S1_energy, S3_energy) in eV, or (None, None) if not found
    """
    if config['DOSOC'] == "TRUE":
        state1_mit_SOC = find_state1_mit_SOC(filename3)  # Direct value
        state3_mit_SOC = find_state3_mit_SOC(filename3)
        return state1_mit_SOC, state3_mit_SOC
    if config['DOSOC'] == "FALSE":
        state1_ohne_SOC = find_state1_ohne_SOC(filename3)  # Direct value
        state3_ohne_SOC = find_state3_ohne_SOC(filename3)
        return state1_ohne_SOC, state3_ohne_SOC
