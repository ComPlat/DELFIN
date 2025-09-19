import re
from pathlib import Path
from typing import Optional, Union, TextIO, Dict


def extract_last_uhf_deviation(
    file: Union[str, Path, TextIO],
    *,
    encoding: str = "utf-8",
    raise_on_missing: bool = False
) -> Optional[float]:
    """
    Parse the ORCA output for the last occurrence of the
    'UHF SPIN CONTAMINATION' block and return the final 'Deviation' value.

    Parameters
    ----------
    file : str | Path | TextIO
        Path to the ORCA output file or an open file handle.
    encoding : str
        File encoding when opening by path.
    raise_on_missing : bool
        If True, raise ValueError when no deviation is found.

    Returns
    -------
    float | None
        The last 'Deviation' value if found; otherwise None (or raises).
    """
    should_close = False
    if hasattr(file, "read"):
        fh = file
    else:
        fh = open(file, "r", encoding=encoding, errors="replace")
        should_close = True

    last_deviation: Optional[float] = None
    in_block = False

    try:
        for line in fh:
            if "UHF SPIN CONTAMINATION" in line:
                in_block = True
                continue
            if in_block:
                if "Deviation" in line:
                    m = re.search(r'([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)', line)
                    if m:
                        try:
                            last_deviation = float(m.group(1))
                        except ValueError:
                            pass
                    in_block = False
        if last_deviation is None and raise_on_missing:
            raise ValueError("No 'Deviation' found in 'UHF SPIN CONTAMINATION' section.")
        return last_deviation
    finally:
        if should_close:
            fh.close()


# --------------------------------------------------------------------
# Spin-Hamiltonian (Heisenberg–Dirac–van Vleck) J-Parser
# Always takes the LAST 'Spin-Hamiltonian Analysis' block in the file.
# We look for lines like:
#   J(3) =   ....... cm**-1  (from -(E[HS]-E[BS])/(<S**2>HS-<S**2>BS))
# --------------------------------------------------------------------

_J_LINE_RE = re.compile(
    r'J\(\s*([123])\s*\)\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)\s*cm\*\*-1',
    re.IGNORECASE
)

def extract_last_J_block(
    file: Union[str, Path, TextIO],
    *,
    encoding: str = "utf-8",
    raise_on_missing: bool = False
) -> Optional[Dict[str, float]]:
    """
    Parse the last 'Spin-Hamiltonian Analysis based on H(HDvV)= -2J*SA*SB'
    block and return a dict of J values: {'J1': ..., 'J2': ..., 'J3': ...}.

    Parameters
    ----------
    file : str | Path | TextIO
        Path to the ORCA output file or an open file handle.
    encoding : str
        File encoding when opening by path.
    raise_on_missing : bool
        If True, raise ValueError when no block is found.

    Returns
    -------
    dict | None
        Dict with any of J1/J2/J3 found in the last block, or None if none found.
    """
    should_close = False
    if hasattr(file, "read"):
        fh = file
    else:
        fh = open(file, "r", encoding=encoding, errors="replace")
        should_close = True

    last_block: Optional[Dict[str, float]] = None
    in_block = False
    have_any_j = False
    current: Dict[str, float] = {}

    try:
        for raw in fh:
            line = raw.rstrip("\n")

            # Detect the beginning of the block
            if ("Spin-Hamiltonian Analysis" in line) and ("H(HDvV)= -2J*SA*SB" in line):
                in_block = True
                have_any_j = False
                current = {}
                continue

            if in_block:
                # Capture J lines
                m = _J_LINE_RE.search(line)
                if m:
                    which = m.group(1)
                    val = float(m.group(2))
                    current[f"J{which}"] = val
                    have_any_j = True
                    continue

                # Heuristic: after capturing at least one J line, a long rule line
                # (or an empty line) usually indicates the end of the framed block.
                if have_any_j and (line.strip().startswith("---") or line.strip() == ""):
                    if current:
                        last_block = dict(current)
                    in_block = False
                    have_any_j = False
                    current = {}
                    continue

        # If file ended while still in a block, finalize it
        if in_block and have_any_j and current:
            last_block = dict(current)

        if last_block is None and raise_on_missing:
            raise ValueError("No Spin-Hamiltonian J block found.")
        return last_block
    finally:
        if should_close:
            fh.close()


def extract_last_J3(
    file: Union[str, Path, TextIO],
    *,
    encoding: str = "utf-8",
    raise_on_missing: bool = False
) -> Optional[float]:
    """
    Convenience helper: return only the last J(3) value found in the file.

    Parameters
    ----------
    file : str | Path | TextIO
        Path to the ORCA output file or an open file handle.
    encoding : str
        File encoding when opening by path.
    raise_on_missing : bool
        If True, raise ValueError when J(3) is not found.

    Returns
    -------
    float | None
        J(3) in cm^-1 if found; otherwise None (or raises).
    """
    block = extract_last_J_block(file, encoding=encoding, raise_on_missing=False)
    if block is None:
        if raise_on_missing:
            raise ValueError("No Spin-Hamiltonian J block found.")
        return None
    j3 = block.get("J3")
    if j3 is None and raise_on_missing:
        raise ValueError("J(3) not found in the last Spin-Hamiltonian block.")
    return j3
