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
        logger.info(f"File {path} not found; skipping energy extraction.")
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
    For S0.out files with multiple jobs, extracts the last energy from Job 1.

    Args:
        filename: Path to ORCA output file

    Returns:
        Electronic energy in Hartree, or None if not found
    """
    text = _read_text(filename)
    if text is None:
        return None

    # Check if this is S0.out with multiple jobs
    if "S0.out" in filename and "JOB NUMBER  2" in text:
        # Find position of Job 2 marker
        job2_match = re.search(r'\$+\s+JOB\s+NUMBER\s+2\s+\$+', text, flags=re.IGNORECASE)
        if job2_match:
            # Search only in text before Job 2
            text_job1 = text[:job2_match.start()]
            # Find last FINAL SINGLE POINT ENERGY in Job 1
            pattern = rf"FINAL\s+SINGLE\s+POINT\s+ENERGY\s+{FLOAT_RE}"
            matches = re.findall(pattern, text_job1, flags=re.IGNORECASE)
            if matches:
                try:
                    return float(matches[-1])
                except ValueError:
                    logger.error(f"Could not convert '{matches[-1]}' to float from {filename}.")
                    return None

    # Default behavior for other files or S0.out without Job 2
    patterns = [
        rf"FINAL\s+SINGLE\s+POINT\s+ENERGY\s+{FLOAT_RE}",
        rf"Total\s+Energy\s*:\s*{FLOAT_RE}",
        rf"Electronic\s+energy.*?{FLOAT_RE}\s*(?:E[hH]|a\.u\.)?",
    ]
    return _search_last_float(filename, patterns)


def _extract_excited_state_energy(filename: str, state_number: int = 1, use_soc: bool = False) -> Optional[float]:
    """Extract excited state energy from ORCA TD-DFT output.

    Consolidated function to extract S1 or S3 energies with or without SOC corrections.

    Args:
        filename: Path to ORCA output file with TD-DFT results
        state_number: State to extract (1 for S1, 3 for S3)
        use_soc: Whether to use SOC-corrected results

    Returns:
        Excited state energy in eV, or None if not found
    """
    if use_soc:
        search_text = "SOC CORRECTED CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
    else:
        search_text = "CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"

    # S1 is at offset +5, S3 is at offset +7
    offset = 5 if state_number == 1 else 7

    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
            for i, line in enumerate(lines):
                if search_text in line:
                    target_line_index = i + offset
                    if target_line_index < len(lines):
                        target_line = lines[target_line_index]
                        parts = target_line[20:].split()
                        try:
                            return float(parts[0])
                        except ValueError:
                            logger.error(f"Could not convert '{parts[0]}' to a float.")
                            return None
    except FileNotFoundError:
        logger.info(f"File {filename} not found; skipping excited-state extraction.")
    return None


def find_state1_ohne_SOC(filename3: str) -> Optional[float]:
    """Extract S1 excited state energy without spin-orbit coupling."""
    return _extract_excited_state_energy(filename3, state_number=1, use_soc=False)


def find_state3_ohne_SOC(filename3: str) -> Optional[float]:
    """Extract S3 excited state energy without spin-orbit coupling."""
    return _extract_excited_state_energy(filename3, state_number=3, use_soc=False)


def find_state1_mit_SOC(filename3: str) -> Optional[float]:
    """Extract S1 excited state energy with spin-orbit coupling corrections."""
    return _extract_excited_state_energy(filename3, state_number=1, use_soc=True)


def find_state3_mit_SOC(filename3: str) -> Optional[float]:
    """Extract S3 excited state energy with spin-orbit coupling corrections."""
    return _extract_excited_state_energy(filename3, state_number=3, use_soc=True)

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
