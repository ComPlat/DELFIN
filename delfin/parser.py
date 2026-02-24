import math
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


def _compute_kleinman_score(beta_tensor: Dict[str, float]) -> float:
    """
    Measure deviation from Kleinman permutation symmetry.

    Groups all 27 β_ijk components by their sorted index triple (canonical
    form) and computes ε = RMS(β_ijk − group_mean) / RMS(β_ijk).

    Returns
    -------
    float
        Score ε ∈ [0, 1]:  0 = perfect Kleinman symmetry,
        near 1 = strongly broken symmetry.
    """
    groups: Dict[str, list] = {}
    for i in 'xyz':
        for j in 'xyz':
            for k in 'xyz':
                canon = ''.join(sorted([i, j, k]))
                groups.setdefault(canon, []).append(beta_tensor.get(i + j + k, 0.0))

    sum_sq_diff = 0.0
    sum_sq_val = 0.0
    for vals in groups.values():
        mean = sum(vals) / len(vals)
        for v in vals:
            sum_sq_diff += (v - mean) ** 2
            sum_sq_val += v ** 2

    if sum_sq_val < 1e-30:
        return 0.0
    return math.sqrt(sum_sq_diff / sum_sq_val)


def _symmetrize_tensor(beta_tensor: Dict[str, float]) -> Dict[str, float]:
    """
    Return a Kleinman-symmetrised copy of beta_tensor.

    Replaces each group of index-permutation-related components with their
    mean, enforcing β_ijk = β_ikj = β_jik = β_jki = β_kij = β_kji.
    The returned dict always contains all 27 keys.
    """
    totals: Dict[str, float] = {}
    counts: Dict[str, int] = {}
    for i in 'xyz':
        for j in 'xyz':
            for k in 'xyz':
                canon = ''.join(sorted([i, j, k]))
                val = beta_tensor.get(i + j + k, 0.0)
                totals[canon] = totals.get(canon, 0.0) + val
                counts[canon] = counts.get(canon, 0) + 1
    means = {c: totals[c] / counts[c] for c in totals}

    sym: Dict[str, float] = {}
    for i in 'xyz':
        for j in 'xyz':
            for k in 'xyz':
                canon = ''.join(sorted([i, j, k]))
                sym[i + j + k] = means[canon]
    return sym


def _resolve_kleinman(mode: str, score: float) -> bool:
    """
    Decide whether to apply Kleinman symmetrisation.

    Parameters
    ----------
    mode : str
        'on'   – always apply
        'off'  – never apply
        'auto' – apply only when score < 0.05
    score : float
        Kleinman asymmetry score from _compute_kleinman_score.

    Notes
    -----
    Kleinman symmetry (β_ijk = β_ikj = … for all permutations) holds in the
    static / off-resonant limit (Kleinman, Phys. Rev. 126, 1977, 1962).
    For ORCA static SCF/DFT calculations the score is typically < 1e-6.
    The threshold ε < 0.05 (5 % RMS asymmetry) is a conservative cut-off:
    at this level the symmetrisation changes individual components by at most
    a few percent, well within the error of the static approximation itself.
    For frequency-dependent β near resonances the score will exceed 0.05 and
    Kleinman will correctly be left unapplied.
    """
    if mode == 'on':
        return True
    if mode == 'off':
        return False
    # 'auto': apply when tensor is near-symmetric (static / off-resonant)
    return score < 0.05


def _get_beta_component(
    beta_tensor: Dict[str, float],
    i: str,
    j: str,
    k: str,
    *,
    kleinman: bool = True,
) -> float:
    """
    Return β_ijk.

    Parameters
    ----------
    kleinman : bool
        If True (default), fall back to permuted indices when the direct
        component is absent – valid under Kleinman / full permutation symmetry
        (Kleinman, Phys. Rev. 126, 1977, 1962).
        If False, return 0.0 for missing components without permutation lookup.

    Notes
    -----
    ORCA static calculations output all 27 Cartesian components explicitly,
    so the permutation fallback is rarely needed in practice.  It is kept as
    a safety net for partial tensors and for use in the β_HRS formula, where
    the formula structure assumes Kleinman-equivalent index groups.
    """
    key = i + j + k
    if key in beta_tensor:
        return beta_tensor[key]
    if kleinman:
        for perm in (f"{i}{k}{j}", f"{j}{i}{k}", f"{j}{k}{i}", f"{k}{i}{j}", f"{k}{j}{i}"):
            if perm in beta_tensor:
                return beta_tensor[perm]
    return 0.0


def calculate_beta_zzz_aligned(
    beta_tensor: Dict[str, float],
    dipole_x: float,
    dipole_y: float,
    dipole_z: float,
    *,
    kleinman_mode: str = 'auto',
) -> tuple:
    """
    Calculate β'_zzz in the dipole-aligned coordinate system.

    Rotates the hyperpolarizability tensor so that the dipole moment points
    along the z-axis, then extracts β'_zzz via the exhaustive 27-term sum:

        β'_zzz = Σ_{i,j,k ∈ {x,y,z}} R_i R_j R_k β_ijk

    where R = dipole / |dipole|.  No implicit multiplicity factors are used,
    so the result is correct regardless of whether permutation symmetry holds.

    Parameters
    ----------
    beta_tensor : dict
        β_ijk components in atomic units.
    dipole_x, dipole_y, dipole_z : float
        Dipole moment components in atomic units.
    kleinman_mode : str
        'auto' – symmetrize tensor before rotation when score ε < 0.05 (default)
        'on'   – always symmetrize before rotation
        'off'  – use raw tensor, no symmetrisation

    Returns
    -------
    tuple
        (beta_zzz_aligned: float  – β'_zzz in a.u.,
         kleinman_score: float    – symmetry deviation ε ∈ [0, 1],
         kleinman_applied: bool   – whether symmetrisation was applied)
    """
    score = _compute_kleinman_score(beta_tensor)
    use_kleinman = _resolve_kleinman(kleinman_mode, score)
    working = _symmetrize_tensor(beta_tensor) if use_kleinman else beta_tensor

    dipole_mag = math.sqrt(dipole_x ** 2 + dipole_y ** 2 + dipole_z ** 2)
    if dipole_mag < 1e-10:
        return working.get('zzz', 0.0), score, use_kleinman

    R = [dipole_x / dipole_mag, dipole_y / dipole_mag, dipole_z / dipole_mag]
    axes = ['x', 'y', 'z']

    # Full 27-term rotation – no implicit symmetry assumed
    result = 0.0
    for a, Ra in zip(axes, R):
        for b, Rb in zip(axes, R):
            for c, Rc in zip(axes, R):
                result += Ra * Rb * Rc * working.get(a + b + c, 0.0)

    return result, score, use_kleinman


def calculate_beta_properties(
    beta_tensor: Dict[str, float],
    dipole_x: float,
    dipole_y: float,
    dipole_z: float,
    *,
    kleinman_mode: str = 'auto',
) -> Dict[str, float]:
    """
    Calculate derived hyperpolarizability properties from the tensor.

    Parameters
    ----------
    beta_tensor : dict
        Dictionary with β_ijk components (e.g., {'xxx': -931.86, 'xxy': -210.84, ...})
    dipole_x, dipole_y, dipole_z : float
        Dipole moment components in atomic units.
    kleinman_mode : str
        Controls Kleinman permutation symmetry handling:
        'auto' – symmetrize tensor when score ε < 0.05 (default)
        'on'   – always symmetrize
        'off'  – never symmetrize (use raw ORCA tensor throughout)

    Returns
    -------
    dict
        Dictionary containing (all units noted):
        beta_x/y/z_au           – contracted vector components (a.u.)
        beta_tot_au             – total magnitude (a.u.)
        beta_mu_au              – dipole projection (a.u.)
        beta_zzz_au             – raw β_zzz component (a.u.)
        beta_zzz_aligned_au     – β'_zzz, dipole-aligned, active mode (a.u.)
        beta_zzz_aligned_raw_au – β'_zzz without Kleinman symmetrisation (a.u.)
        beta_zzz_aligned_sym_au – β'_zzz with Kleinman symmetrisation (a.u.)
        beta_HRS_au             – Hyper-Rayleigh scattering β_HRS (a.u.)
        kleinman_score          – asymmetry score ε ∈ [0, 1]
        kleinman_applied        – whether symmetrisation was used
        … (plus _esu and _esu_30 unit variants for all quantities above)
    """
    # Conversion factor: 1 a.u. = 8.6393 × 10⁻³³ esu
    AU_TO_ESU = 8.6393e-33

    # Compute Kleinman score and decide symmetrisation once
    score = _compute_kleinman_score(beta_tensor)
    use_kleinman = _resolve_kleinman(kleinman_mode, score)
    working = _symmetrize_tensor(beta_tensor) if use_kleinman else beta_tensor

    # ---------- contracted vector components ----------
    # The "vector part" of β is defined by contracting over two indices:
    #   βᵢ = Σⱼ β_ijj   (sum over j for fixed i)
    # This gives a polar vector (βₓ, βᵧ, β_z) whose magnitude β_tot and
    # projection onto the dipole β_μ are orientation-independent scalars
    # commonly reported for push-pull chromophores.
    # Reference: Champagne & Bishop, Adv. Chem. Phys. 126, 41 (2003), eq. 2.
    beta_x = (
        _get_beta_component(working, 'x', 'x', 'x')
        + _get_beta_component(working, 'x', 'y', 'y')
        + _get_beta_component(working, 'x', 'z', 'z')
    )
    beta_y = (
        _get_beta_component(working, 'y', 'x', 'x')
        + _get_beta_component(working, 'y', 'y', 'y')
        + _get_beta_component(working, 'y', 'z', 'z')
    )
    beta_z = (
        _get_beta_component(working, 'z', 'x', 'x')
        + _get_beta_component(working, 'z', 'y', 'y')
        + _get_beta_component(working, 'z', 'z', 'z')
    )

    beta_tot = math.sqrt(beta_x ** 2 + beta_y ** 2 + beta_z ** 2)

    # β_μ: projection of the β-vector onto the ground-state dipole moment.
    # This is what EFISH experiments measure (after orientational averaging).
    dipole_mag = math.sqrt(dipole_x ** 2 + dipole_y ** 2 + dipole_z ** 2)
    if dipole_mag > 0:
        beta_mu = (beta_x * dipole_x + beta_y * dipole_y + beta_z * dipole_z) / dipole_mag
    else:
        beta_mu = 0.0

    # β_zzz: the single component along the laboratory z-axis (not dipole-aligned).
    beta_zzz = _get_beta_component(working, 'z', 'z', 'z')

    # ---------- β'_zzz: three variants via full 27-term rotation ----------
    # Active mode (honors kleinman_mode)
    beta_zzz_aligned, _, _ = calculate_beta_zzz_aligned(
        beta_tensor, dipole_x, dipole_y, dipole_z, kleinman_mode=kleinman_mode
    )
    # Always raw (no symmetrisation)
    beta_zzz_aligned_raw, _, _ = calculate_beta_zzz_aligned(
        beta_tensor, dipole_x, dipole_y, dipole_z, kleinman_mode='off'
    )
    # Always Kleinman-symmetrised
    beta_zzz_aligned_sym, _, _ = calculate_beta_zzz_aligned(
        beta_tensor, dipole_x, dipole_y, dipole_z, kleinman_mode='on'
    )

    # ---------- β_HRS (isotropic rotational average) ----------
    #
    # Physical setup: 90° Hyper-Rayleigh Scattering (HRS) in isotropic solution.
    #   - Incident beam polarised along Z (lab frame).
    #   - Detection along X (perpendicular to beam), collecting both polarisations.
    #   - Molecular orientation fully randomised → orientational average over SO(3).
    #
    # The HRS intensity is proportional to ⟨|β_ZZZ|²⟩ + ⟨|β_XZZ|²⟩  (lab frame).
    # Carrying out the SO(3) average under Kleinman symmetry
    # (β_ijk = β_ikj = β_jik = …, valid for static / off-resonant β) gives:
    #
    #   β²_HRS = (6/35)   · Σ_i   β_iii²
    #           + (16/105) · Σ_{i≠j} β_iii · β_iij
    #           + (38/105) · Σ_{i≠j} β_iij²
    #           + (16/105) · Σ_cyclic β_iij · β_jkk      (cyclic: xyz, yzx, zxy)
    #           + (4/7)    · β_xyz²
    #
    # Note: 20/35 = 4/7.  The formula uses the active 'working' tensor
    # (Kleinman-symmetrised or raw, depending on kleinman_mode).
    # β_HRS is NOT divided by 2 – no SHG convention factor is applied here.
    axes = ("x", "y", "z")
    sum_iii_sq = sum(_get_beta_component(working, i, i, i) ** 2 for i in axes)

    sum_iii_iij = 0.0
    sum_iij_sq = 0.0
    for i in axes:
        beta_iii = _get_beta_component(working, i, i, i)
        for j in axes:
            if i == j:
                continue
            beta_iij = _get_beta_component(working, i, i, j)
            sum_iii_iij += beta_iii * beta_iij
            sum_iij_sq += beta_iij ** 2

    cyclic_triplets = (("x", "y", "z"), ("y", "z", "x"), ("z", "x", "y"))
    sum_cyclic_iij_jkk = sum(
        _get_beta_component(working, i, i, j) * _get_beta_component(working, j, k, k)
        for i, j, k in cyclic_triplets
    )

    beta_ijk = _get_beta_component(working, "x", "y", "z")
    beta_hrs_sq = (
        (6.0 / 35.0) * sum_iii_sq
        + (16.0 / 105.0) * sum_iii_iij
        + (38.0 / 105.0) * sum_iij_sq
        + (16.0 / 105.0) * sum_cyclic_iij_jkk
        + (20.0 / 35.0) * (beta_ijk ** 2)   # 20/35 = 4/7
    )
    beta_hrs = math.sqrt(max(beta_hrs_sq, 0.0))

    # ---------- unit conversions ----------
    beta_tot_esu = beta_tot * AU_TO_ESU
    beta_mu_esu = beta_mu * AU_TO_ESU
    beta_zzz_esu = beta_zzz * AU_TO_ESU
    beta_zzz_aligned_esu = beta_zzz_aligned * AU_TO_ESU
    beta_hrs_esu = beta_hrs * AU_TO_ESU

    beta_tot_esu_30 = beta_tot_esu * 1e30
    beta_mu_esu_30 = beta_mu_esu * 1e30
    beta_zzz_aligned_esu_30 = beta_zzz_aligned_esu * 1e30
    beta_zzz_aligned_raw_esu_30 = beta_zzz_aligned_raw * AU_TO_ESU * 1e30
    beta_zzz_aligned_sym_esu_30 = beta_zzz_aligned_sym * AU_TO_ESU * 1e30
    beta_hrs_esu_30 = beta_hrs_esu * 1e30

    # Convention variants (Willets et al., J. Chem. Phys. 97, 7590, 1992):
    #
    # T convention (Taylor, used by ORCA and most QC codes):
    #   μᵢ(F) = μᵢ⁰ + αᵢⱼFⱼ + ½ βᵢⱼₖ^T FⱼFₖ + ...
    #   → All _T_ values are the raw ORCA output.
    #
    # B convention (used by many experimental SHG/HRS papers):
    #   μᵢ(F) = μᵢ⁰ + αᵢⱼFⱼ + βᵢⱼₖ^B FⱼFₖ + ...   (no ½)
    #   → β^B = ½ β^T   →   All _B_ values = ORCA / 2.
    #
    # When comparing with experiment, check which convention the paper uses.
    beta_zzz_aligned_sym_T_esu_30 = beta_zzz_aligned_sym_esu_30          # ORCA / T
    beta_zzz_aligned_sym_B_esu_30 = beta_zzz_aligned_sym_esu_30 / 2.0   # exp  / B
    beta_hrs_T_esu_30 = beta_hrs_esu_30                                   # ORCA / T
    beta_hrs_B_esu_30 = beta_hrs_esu_30 / 2.0                            # exp  / B

    return {
        # Atomic units
        'beta_x_au': beta_x,
        'beta_y_au': beta_y,
        'beta_z_au': beta_z,
        'beta_tot_au': beta_tot,
        'beta_mu_au': beta_mu,
        'beta_zzz_au': beta_zzz,
        'beta_zzz_aligned_au': beta_zzz_aligned,
        'beta_zzz_aligned_raw_au': beta_zzz_aligned_raw,
        'beta_zzz_aligned_sym_au': beta_zzz_aligned_sym,
        # esu
        'beta_x_esu': beta_x * AU_TO_ESU,
        'beta_y_esu': beta_y * AU_TO_ESU,
        'beta_z_esu': beta_z * AU_TO_ESU,
        'beta_tot_esu': beta_tot_esu,
        'beta_mu_esu': beta_mu_esu,
        'beta_zzz_esu': beta_zzz_esu,
        'beta_zzz_aligned_esu': beta_zzz_aligned_esu,
        # 10^-30 esu
        'beta_tot_esu_30': beta_tot_esu_30,
        'beta_mu_esu_30': beta_mu_esu_30,
        'beta_zzz_aligned_esu_30': beta_zzz_aligned_esu_30,
        'beta_zzz_aligned_raw_esu_30': beta_zzz_aligned_raw_esu_30,
        'beta_zzz_aligned_sym_esu_30': beta_zzz_aligned_sym_esu_30,
        # backward-compat key
        'beta_zzz_aligned_esu_30_kleinman': beta_zzz_aligned_sym_esu_30,
        # T convention (ORCA / QC): μ = μ⁰ + αE + ½βE² + ...
        'beta_zzz_aligned_sym_T_esu_30': beta_zzz_aligned_sym_T_esu_30,
        'beta_HRS_T_esu_30': beta_hrs_T_esu_30,
        # B convention (experiment): μ = μ⁰ + αE + βE² + ...  → β^B = β^T / 2
        'beta_zzz_aligned_sym_B_esu_30': beta_zzz_aligned_sym_B_esu_30,
        'beta_HRS_B_esu_30': beta_hrs_B_esu_30,
        # β_HRS base values
        'beta_HRS': beta_hrs,
        'beta_HRS_au': beta_hrs,
        'beta_HRS_esu': beta_hrs_esu,
        'beta_HRS_esu_30': beta_hrs_esu_30,
        # Kleinman metadata
        'kleinman_score': score,
        'kleinman_applied': use_kleinman,
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
