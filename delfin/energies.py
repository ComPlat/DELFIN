# energies.py
import re
import logging
from typing import Optional, Sequence, Dict, Any, Tuple

FLOAT_RE = r'([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)'

def _read_text(path: str) -> Optional[str]:
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            return f.read()
    except FileNotFoundError:
        logging.error(f"File {path} not found.")
        return None

def _search_last_float(path: str, patterns: Sequence[str]) -> Optional[float]:
    text = _read_text(path)
    if text is None:
        return None
    for pat in patterns:
        matches = re.findall(pat, text, flags=re.IGNORECASE | re.DOTALL)
        if matches:
            try:
                return float(matches[-1])
            except ValueError:
                logging.error(f"Could not convert '{matches[-1]}' to float from {path}.")
                return None
    return None

def find_gibbs_energy(filename: str) -> Optional[float]:

    patterns = [
        rf"Final\s+Gibbs\s+free\s+energy.*?{FLOAT_RE}\s*(?:E[hH]|a\.u\.)?",
        rf"Gibbs\s+free\s+energy.*?{FLOAT_RE}\s*(?:E[hH]|a\.u\.)?",
    ]
    return _search_last_float(filename, patterns)

def find_ZPE(filename: str) -> Optional[float]:
    patterns = [
        rf"Zero[\s-]?point\s+energy.*?{FLOAT_RE}\s*(?:E[hH]|a\.u\.)?",
        rf"Zero[\s-]?point\s+correction.*?{FLOAT_RE}\s*(?:E[hH]|a\.u\.)?",
    ]
    return _search_last_float(filename, patterns)

def find_electronic_energy(filename: str) -> Optional[float]:
    patterns = [
        rf"FINAL\s+SINGLE\s+POINT\s+ENERGY\s+{FLOAT_RE}",
        rf"Total\s+Energy\s*:\s*{FLOAT_RE}",
        rf"Electronic\s+energy.*?{FLOAT_RE}\s*(?:E[hH]|a\.u\.)?",
    ]
    return _search_last_float(filename, patterns)


def find_state1_ohne_SOC(filename3: str) -> Optional[float]:
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
                            logging.error(f"Could not convert '{parts[0]}' to a float.")
                            continue
    except FileNotFoundError:
        logging.error(f"File {filename3} not found.")
        return None
    return last_value

def find_state3_ohne_SOC(filename3: str) -> Optional[float]:
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
                            logging.error(f"Could not convert '{parts[0]}' to a float.")
                            continue
    except FileNotFoundError:
        logging.error(f"File {filename3} not found.")
        return None
    return last_value


def find_state1_mit_SOC(filename3: str) -> Optional[float]:
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
                            logging.error(f"Could not convert '{parts[0]}' to a float.")
                            return None
    except FileNotFoundError:
        logging.error(f"File {filename3} not found.")
    return None

def find_state3_mit_SOC(filename3: str) -> Optional[float]:
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
                            logging.error(f"Could not convert '{parts[0]}' to a float.")
                            return None
    except FileNotFoundError:
        logging.error(f"File {filename3} not found.")
    return None

def check_and_execute_SOC(filename3: str, config: Dict[str, Any]) -> Tuple[Optional[float], Optional[float]]:
    if config['DOSOC'] == "TRUE":
        state1_mit_SOC = find_state1_mit_SOC(filename3)  # Direct value
        state3_mit_SOC = find_state3_mit_SOC(filename3)
        return state1_mit_SOC, state3_mit_SOC
    if config['DOSOC'] == "FALSE":
        state1_ohne_SOC = find_state1_ohne_SOC(filename3)  # Direct value
        state3_ohne_SOC = find_state3_ohne_SOC(filename3)
        return state1_ohne_SOC, state3_ohne_SOC
