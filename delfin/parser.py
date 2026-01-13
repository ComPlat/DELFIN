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


def parse_hyperpolarizability(
    file: Union[str, Path, TextIO],
    *,
    encoding: str = "utf-8",
) -> Optional[Dict[str, float]]:
    """
    Parse the static hyperpolarizability tensor from ORCA output.

    Searches for 'STATIC HYPERPOLARIZABILITY TENSOR' and extracts all 27
    Cartesian components β_ijk (i,j,k ∈ {x,y,z}).

    Parameters
    ----------
    file : str | Path | TextIO
        Path to the ORCA output file or an open file handle.
    encoding : str
        File encoding when opening by path.

    Returns
    -------
    dict | None
        Dictionary with keys like 'xxx', 'xxy', 'xyz', etc. containing
        the β_ijk values in atomic units, or None if not found.
    """
    should_close = False
    if hasattr(file, "read"):
        fh = file
    else:
        fh = open(file, "r", encoding=encoding, errors="replace")
        should_close = True

    tensor = {}
    in_tensor = False

    try:
        for line in fh:
            # Look for the tensor header
            if "STATIC HYPERPOLARIZABILITY TENSOR" in line:
                in_tensor = True
                continue

            if in_tensor:
                # Look for lines like "     ( x x x ):         -931.858682"
                match = re.match(r'\s*\(\s*([xyz])\s+([xyz])\s+([xyz])\s*\)\s*:\s*([+-]?\d+\.\d+)', line)
                if match:
                    i, j, k, value = match.groups()
                    key = f"{i}{j}{k}"
                    tensor[key] = float(value)

                # End of tensor section (empty line or next section)
                if line.strip() == "" and tensor:
                    break

        return tensor if tensor else None
    finally:
        if should_close:
            fh.close()


def calculate_beta_zzz_aligned(
    beta_tensor: Dict[str, float],
    dipole_x: float,
    dipole_y: float,
    dipole_z: float,
) -> float:
    """
    Calculate β'_zzz in the dipole-aligned coordinate system.

    This rotates the hyperpolarizability tensor so that the dipole moment
    points along the z-axis, then extracts the β'_zzz component.

    Parameters
    ----------
    beta_tensor : dict
        Dictionary with β_ijk components in atomic units
    dipole_x, dipole_y, dipole_z : float
        Dipole moment components in atomic units

    Returns
    -------
    float
        β'_zzz in the dipole-aligned coordinate system (a.u.)
    """
    import math

    # Calculate dipole magnitude
    dipole_mag = math.sqrt(dipole_x**2 + dipole_y**2 + dipole_z**2)

    # Handle zero dipole case
    if dipole_mag < 1e-10:
        # No preferred direction - return original beta_zzz
        return beta_tensor.get('zzz', 0.0)

    # Normalized dipole vector components (third row of rotation matrix)
    R_31 = dipole_x / dipole_mag
    R_32 = dipole_y / dipole_mag
    R_33 = dipole_z / dipole_mag

    # Helper function to get tensor components with permutation symmetry
    def get_beta(i: str, j: str, k: str) -> float:
        """Get β_ijk, trying all permutations due to tensor symmetry."""
        # Try direct lookup first
        key = i + j + k
        if key in beta_tensor:
            return beta_tensor[key]
        # Try all 6 permutations (tensor is symmetric under index exchange)
        for perm in [f"{i}{k}{j}", f"{j}{i}{k}", f"{j}{k}{i}", f"{k}{i}{j}", f"{k}{j}{i}"]:
            if perm in beta_tensor:
                return beta_tensor[perm]
        return 0.0

    # Calculate β'_zzz using tensor transformation
    # β'_zzz = Σ_ijk R_3i R_3j R_3k β_ijk
    # With symmetry: use reduced form with appropriate multiplicities

    beta_zzz_new = (
        R_31**3 * get_beta('x', 'x', 'x') +
        3 * R_31**2 * R_32 * get_beta('x', 'x', 'y') +
        3 * R_31**2 * R_33 * get_beta('x', 'x', 'z') +
        3 * R_31 * R_32**2 * get_beta('x', 'y', 'y') +
        6 * R_31 * R_32 * R_33 * get_beta('x', 'y', 'z') +
        3 * R_31 * R_33**2 * get_beta('x', 'z', 'z') +
        R_32**3 * get_beta('y', 'y', 'y') +
        3 * R_32**2 * R_33 * get_beta('y', 'y', 'z') +
        3 * R_32 * R_33**2 * get_beta('y', 'z', 'z') +
        R_33**3 * get_beta('z', 'z', 'z')
    )

    return beta_zzz_new


def calculate_beta_properties(
    beta_tensor: Dict[str, float],
    dipole_x: float,
    dipole_y: float,
    dipole_z: float,
) -> Dict[str, float]:
    """
    Calculate derived hyperpolarizability properties from the tensor.

    Parameters
    ----------
    beta_tensor : dict
        Dictionary with β_ijk components (e.g., {'xxx': -931.86, 'xxy': -210.84, ...})
    dipole_x, dipole_y, dipole_z : float
        Dipole moment components in atomic units

    Returns
    -------
    dict
        Dictionary containing:
        - 'beta_x_au', 'beta_y_au', 'beta_z_au': Contracted vector components (a.u.)
        - 'beta_tot_au': Total hyperpolarizability magnitude (a.u.)
        - 'beta_mu_au': Projection onto dipole vector (a.u.)
        - 'beta_zzz_au': β_zzz tensor component (a.u.)
        - 'beta_zzz_aligned_au': β'_zzz rotated to dipole-aligned frame (a.u.)
        - 'beta_x_esu', 'beta_y_esu', 'beta_z_esu': Converted to esu
        - 'beta_tot_esu': Total in esu
        - 'beta_mu_esu': Projection in esu
        - 'beta_zzz_esu': β_zzz tensor component in esu
        - 'beta_zzz_aligned_esu': β'_zzz in dipole-aligned frame (esu)
    """
    # Conversion factor: 1 a.u. = 8.6393 × 10⁻³³ esu
    AU_TO_ESU = 8.6393e-33

    # Calculate contracted vector components
    beta_x = (beta_tensor.get('xxx', 0) +
              beta_tensor.get('xyy', 0) +
              beta_tensor.get('xzz', 0))

    beta_y = (beta_tensor.get('yxx', 0) +
              beta_tensor.get('yyy', 0) +
              beta_tensor.get('yzz', 0))

    beta_z = (beta_tensor.get('zxx', 0) +
              beta_tensor.get('zyy', 0) +
              beta_tensor.get('zzz', 0))

    # Total hyperpolarizability
    import math
    beta_tot = math.sqrt(beta_x**2 + beta_y**2 + beta_z**2)

    # Projection onto dipole vector
    dipole_mag = math.sqrt(dipole_x**2 + dipole_y**2 + dipole_z**2)
    if dipole_mag > 0:
        beta_mu = (beta_x * dipole_x + beta_y * dipole_y + beta_z * dipole_z) / dipole_mag
    else:
        beta_mu = 0.0

    # Extract β_zzz tensor component
    beta_zzz = beta_tensor.get('zzz', 0.0)

    # Calculate β'_zzz in dipole-aligned coordinate system
    beta_zzz_aligned = calculate_beta_zzz_aligned(beta_tensor, dipole_x, dipole_y, dipole_z)

    return {
        'beta_x_au': beta_x,
        'beta_y_au': beta_y,
        'beta_z_au': beta_z,
        'beta_tot_au': beta_tot,
        'beta_mu_au': beta_mu,
        'beta_zzz_au': beta_zzz,
        'beta_zzz_aligned_au': beta_zzz_aligned,
        'beta_x_esu': beta_x * AU_TO_ESU,
        'beta_y_esu': beta_y * AU_TO_ESU,
        'beta_z_esu': beta_z * AU_TO_ESU,
        'beta_tot_esu': beta_tot * AU_TO_ESU,
        'beta_mu_esu': beta_mu * AU_TO_ESU,
        'beta_zzz_esu': beta_zzz * AU_TO_ESU,
        'beta_zzz_aligned_esu': beta_zzz_aligned * AU_TO_ESU,
    }


def parse_polarizability(
    file: Union[str, Path, TextIO],
    *,
    encoding: str = "utf-8",
) -> Optional[Dict[str, float]]:
    """
    Parse the static polarizability from ORCA output.

    Searches for 'STATIC POLARIZABILITY TENSOR' and extracts the
    isotropic polarizability value.

    Parameters
    ----------
    file : str | Path | TextIO
        Path to the ORCA output file or an open file handle.
    encoding : str
        File encoding when opening by path.

    Returns
    -------
    dict | None
        Dictionary with 'isotropic_au' and 'isotropic_angstrom3' keys,
        or None if not found.
    """
    should_close = False
    if hasattr(file, "read"):
        fh = file
    else:
        fh = open(file, "r", encoding=encoding, errors="replace")
        should_close = True

    isotropic = None

    try:
        for line in fh:
            # Look for "Isotropic polarizability :  12.86035"
            if "Isotropic polarizability" in line:
                match = re.search(r'Isotropic polarizability\s*:\s*([+-]?\d+\.\d+)', line)
                if match:
                    isotropic = float(match.group(1))
                    break

        if isotropic is not None:
            # Conversion: 1 a.u. = 0.1482 Å³
            AU_TO_ANGSTROM3 = 0.1482
            return {
                'isotropic_au': isotropic,
                'isotropic_angstrom3': isotropic * AU_TO_ANGSTROM3,
            }
        return None
    finally:
        if should_close:
            fh.close()
