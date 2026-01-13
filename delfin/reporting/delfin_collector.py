"""
DELFIN Data Collector - Collects all data from DELFIN calculations into JSON format.

This module parses:
- initial.inp (metadata, molecule info)
- initial.out (initial S0 optimization)
- ox_step_1.out, red_step_1.out (oxidized/reduced states)
- OCCUPIER/ folder (charge transfer analysis, if exists)
- ESD/ folder (all excited states, ISC, IC)
- All .xyz geometries
- All thermochemistry data (Gibbs, ZPE, etc.)
- All spectroscopic data (absorption, emission)
- All rate constants (ISC, IC, radiative)
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, Optional

from delfin.common.logging import get_logger
from delfin.uv_vis_spectrum import parse_absorption_spectrum
from delfin.utils import get_git_commit_info

logger = get_logger(__name__)

HARTREE_TO_EV = 27.211386245988
HARTREE_TO_KJ_MOL = 2625.499639


def parse_delfin_txt(project_dir: Path) -> Dict[str, Any]:
    """Parse DELFIN.txt summary file for redox potentials and other info."""
    delfin_txt = project_dir / "DELFIN.txt"
    if not delfin_txt.exists():
        return {}

    data = {
        "redox_potentials_vs_fc": {},
        "total_run_time": None,
    }

    try:
        with open(delfin_txt, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        # Parse redox potentials (V vs. Fc+/Fc)
        redox_patterns = [
            (r'E_red\s*=\s*([-\d.]+)', 'E_red'),
            (r'E_red_2\s*=\s*([-\d.]+)', 'E_red_2'),
            (r'E_red_3\s*=\s*([-\d.]+)', 'E_red_3'),
            (r'E_ox\s*=\s*([-\d.]+)', 'E_ox'),
            (r'E_ox_2\s*=\s*([-\d.]+)', 'E_ox_2'),
            (r'E_ox_3\s*=\s*([-\d.]+)', 'E_ox_3'),
            (r'E_ref\s*=\s*([-\d.]+)', 'E_ref'),
        ]

        for pattern, key in redox_patterns:
            match = re.search(pattern, content)
            if match:
                try:
                    data["redox_potentials_vs_fc"][key] = float(match.group(1))
                except ValueError:
                    pass

        # Parse total run time
        time_match = re.search(r'TOTAL RUN TIME:\s*(\d+)\s*hours?\s*(\d+)\s*minutes?\s*([\d.]+)\s*seconds?', content, re.IGNORECASE)
        if time_match:
            hours = int(time_match.group(1))
            minutes = int(time_match.group(2))
            seconds = float(time_match.group(3))
            data["total_run_time"] = {
                "hours": hours,
                "minutes": minutes,
                "seconds": seconds,
                "total_seconds": hours * 3600 + minutes * 60 + seconds
            }

    except Exception as e:
        logger.error(f"Error parsing DELFIN.txt: {e}")

    return data


def parse_control_flags(project_dir: Path) -> Dict[str, Any]:
    """Parse CONTROL.txt for simple flags (e.g., IMAG=yes)."""
    control = project_dir / "CONTROL.txt"
    flags: Dict[str, Any] = {}
    if not control.exists():
        return flags

    try:
        content = control.read_text(encoding="utf-8", errors="ignore")
        imag_match = re.search(r'^\s*IMAG\s*=\s*(\w+)', content, flags=re.IGNORECASE | re.MULTILINE)
        if imag_match:
            flags["imag"] = imag_match.group(1).strip().lower() == "yes"
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to parse CONTROL.txt for flags: %s", exc)

    return flags


def parse_initial_inp(project_dir: Path) -> Dict[str, Any]:
    """Parse initial.inp (or fallback S0.inp) file for metadata."""
    initial_inp = project_dir / "initial.inp"
    # Fallback: use S0.inp inside ESD if initial.inp is missing
    if not initial_inp.exists():
        s0_inp = project_dir / "ESD" / "S0.inp"
        if s0_inp.exists():
            initial_inp = s0_inp
            logger.info(f"initial.inp not found, using fallback {s0_inp}")
        else:
            logger.warning(f"initial.inp not found in {project_dir} and no fallback S0.inp")
            return {}

    metadata = {
        "functional": None,
        "basis_set": None,
        "auxiliary_basis": None,
        "ri_method": None,
        "dispersion_correction": None,
        "implicit_solvation": None,
        "solvent": None,
        "charge": 0,
        "multiplicity": 1
    }

    try:
        with open(initial_inp, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        # Extract from ! line
        keyword_match = re.search(r'!\s+(.+)', content)
        if keyword_match:
            keywords = keyword_match.group(1).split()
            for kw in keywords:
                kw_upper = kw.upper()
                if any(func in kw_upper for func in ['PBE0', 'PBE', 'B3LYP', 'WB97', 'CAM-B3LYP', 'M062X', 'TPSS']):
                    if not metadata['functional']:  # Don't overwrite
                        metadata['functional'] = kw
                elif any(basis in kw_upper for basis in ['DEF2-SVP', 'DEF2-TZVP', 'DEF2-TZVPP', '6-31G', 'CC-PVDZ']):
                    metadata['basis_set'] = kw
                elif '/J' in kw or '-J' in kw or 'JKFIT' in kw_upper or 'AUX' in kw_upper:
                    metadata['auxiliary_basis'] = kw
                elif any(token in kw_upper for token in ['RIJCOSX', 'RICOSX', 'RIJONX', 'RI']):
                    metadata['ri_method'] = kw
                elif 'D3' in kw_upper or 'D4' in kw_upper:
                    metadata['dispersion_correction'] = kw
                elif 'CPCM' in kw_upper:
                    metadata['implicit_solvation'] = 'CPCM'
                elif 'SMD' in kw_upper:
                    metadata['implicit_solvation'] = 'SMD'

        # Extract solvent
        solvent_match = re.search(r'(CPCM|SMD)\((\w+)\)', content, re.IGNORECASE)
        if solvent_match:
            metadata['solvent'] = solvent_match.group(2)

        # Extract charge and multiplicity from coordinate block
        coord_match = re.search(r'\*\s+xyz\s+(-?\d+)\s+(\d+)', content)
        if coord_match:
            metadata['charge'] = int(coord_match.group(1))
            metadata['multiplicity'] = int(coord_match.group(2))

    except Exception as e:
        logger.error(f"Error parsing initial.inp: {e}")

    return metadata


def parse_tddft_settings(inp_file: Path) -> Dict[str, Any]:
    """Extract key TDDFT keywords (e.g., NROOTS, MAXDIM, FOLLOWIROOT) from an ORCA input."""
    settings: Dict[str, Any] = {}
    if not inp_file.exists():
        return settings

    try:
        content = inp_file.read_text(encoding="utf-8", errors="ignore")
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to read %s for TDDFT settings: %s", inp_file, exc)
        return settings

    # Only scan the TDDFT block to avoid false positives elsewhere
    tddft_blocks = re.findall(r"%\s*TDDFT(.*?)(?:\n\s*end|\Z)", content, flags=re.IGNORECASE | re.DOTALL)
    block_text = "\n".join(tddft_blocks) if tddft_blocks else content

    def _extract_int(keyword: str) -> Optional[int]:
        match = re.search(rf"{keyword}\s*=?\s*([-\d]+)", block_text, flags=re.IGNORECASE)
        if match:
            try:
                return int(match.group(1))
            except Exception:
                return match.group(1)
        return None

    for key in ["NROOTS", "MAXDIM", "FOLLOWIROOT"]:
        val = _extract_int(key)
        if val is not None:
            settings[key.lower()] = val

    return settings


def parse_charge_multiplicity(inp_file: Path) -> Optional[Dict[str, int]]:
    """Parse charge and multiplicity from an ORCA input file."""
    if not inp_file.exists():
        return None

    try:
        with open(inp_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        coord_match = re.search(r'\*\s+xyz\s+(-?\d+)\s+(\d+)', content)
        if coord_match:
            return {
                "charge": int(coord_match.group(1)),
                "multiplicity": int(coord_match.group(2))
            }
    except Exception as e:
        logger.error(f"Error parsing charge/multiplicity from {inp_file}: {e}")

    return None


def parse_orbitals_from_occuper(out_file: Path) -> Optional[Dict[str, Any]]:
    """Parse orbitals from the matched OCCUPIER .out file if present (last block)."""
    return parse_orbitals(out_file) if out_file.exists() else None


def parse_scf_energy(output_file: Path) -> Optional[Dict[str, Any]]:
    """Parse SCF energy from ORCA output file.

    Uses the last FINAL SINGLE POINT ENERGY before 'OPTIMIZATION RUN DONE'.
    This excludes subsequent single-point calculations (e.g., TDDFT at optimized geometry).
    """
    if not output_file.exists():
        return None

    try:
        with open(output_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        # Find position of optimization completion marker
        opt_done_match = re.search(r'\*\*\* OPTIMIZATION RUN DONE \*\*\*', content)
        if opt_done_match:
            # Only consider energies before optimization completion
            search_content = content[:opt_done_match.start()]
        else:
            # No optimization marker - use entire file
            search_content = content

        # Find all FINAL SINGLE POINT ENERGY entries in search region
        energy_matches = re.findall(r'FINAL SINGLE POINT ENERGY\s+([-\d.]+)', search_content)

        if energy_matches:
            # Use last FSPE before OPTIMIZATION RUN DONE
            energy_hartree = float(energy_matches[-1])
            return {
                "hartree": energy_hartree,
                "eV": energy_hartree * HARTREE_TO_EV,
                "kJ_mol": energy_hartree * HARTREE_TO_KJ_MOL
            }
    except Exception as e:
        logger.error(f"Error parsing SCF energy from {output_file}: {e}")

    return None


def parse_thermochemistry(output_file: Path) -> Optional[Dict[str, Any]]:
    """Parse thermochemistry data from frequency calculation."""
    if not output_file.exists():
        return None

    try:
        with open(output_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        thermo = {}

        # Temperature and pressure
        temp_match = re.search(r'Temperature\s+\.\.\.\s+([\d.]+)\s+K', content)
        press_match = re.search(r'Pressure\s+\.\.\.\s+([\d.]+)\s+atm', content)

        if temp_match:
            thermo['temperature_K'] = float(temp_match.group(1))
        if press_match:
            thermo['pressure_atm'] = float(press_match.group(1))

        # Zero-point energy
        zpe_match = re.search(r'Zero point energy\s+\.\.\.\s+([\d.]+)\s+Eh', content)
        if zpe_match:
            zpe_hartree = float(zpe_match.group(1))
            thermo['zero_point_energy_hartree'] = zpe_hartree
            thermo['zero_point_energy_eV'] = zpe_hartree * HARTREE_TO_EV
            thermo['zero_point_energy_kJ_mol'] = zpe_hartree * HARTREE_TO_KJ_MOL

        # Thermal corrections / thermodynamic totals
        thermal_energy_match = re.search(r'Total thermal energy\s+([-\d.]+)\s+Eh', content)
        thermal_enthalpy_match = re.search(r'Total enthalpy\s+([-\d.]+)\s+Eh', content)
        thermal_gibbs_match = re.search(r'Final Gibbs free energy\s+\.\.\.\s+([-\d.]+)\s+Eh', content)

        if thermal_energy_match:
            thermo["total_thermal_energy_hartree"] = float(thermal_energy_match.group(1))
            thermo["total_thermal_energy_eV"] = thermo["total_thermal_energy_hartree"] * HARTREE_TO_EV

        if thermal_enthalpy_match:
            thermo["total_enthalpy_hartree"] = float(thermal_enthalpy_match.group(1))
            thermo["total_enthalpy_eV"] = thermo["total_enthalpy_hartree"] * HARTREE_TO_EV

        if thermal_gibbs_match:
            gibbs_hartree = float(thermal_gibbs_match.group(1))
            thermo['total_gibbs_free_energy_hartree'] = gibbs_hartree
            thermo['total_gibbs_free_energy_eV'] = gibbs_hartree * HARTREE_TO_EV

        # Entropy
        entropy_match = re.search(r'Total entropy correction\s+([-\d.]+)\s+Eh', content)
        if entropy_match:
            thermo['entropy_correction_hartree'] = float(entropy_match.group(1))

        return thermo if thermo else None

    except Exception as e:
        logger.error(f"Error parsing thermochemistry from {output_file}: {e}")

    return None


def parse_dipole_moment(output_file: Path) -> Optional[Dict[str, Any]]:
    """Parse permanent dipole moment from ORCA output."""
    if not output_file.exists():
        return None

    try:
        with open(output_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        # Find the LAST dipole moment (after optimization)
        # Pattern matches:
        # Total Dipole Moment    :     -0.286613982       1.115620361       0.186158294
        #                         -----------------------------------------
        # Magnitude (a.u.)       :      1.166795301
        # Magnitude (Debye)      :      2.965757963
        dipole_matches = re.findall(
            r'Total Dipole Moment\s+:\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+).*?'
            r'Magnitude \(a\.u\.\)\s+:\s+([-\d.]+).*?'
            r'Magnitude \(Debye\)\s+:\s+([-\d.]+)',
            content,
            re.DOTALL
        )

        if dipole_matches:
            # Take the last match
            x_au, y_au, z_au, mag_au, mag_debye = dipole_matches[-1]
            return {
                "x_au": float(x_au),
                "y_au": float(y_au),
                "z_au": float(z_au),
                "magnitude_au": float(mag_au),
                "magnitude_debye": float(mag_debye),
                "vector_debye": [float(x_au), float(y_au), float(z_au)]  # Actually in a.u. but kept for compatibility
            }
    except Exception as e:
        logger.error(f"Error parsing dipole moment from {output_file}: {e}")

    return None


def parse_orbitals(output_file: Path) -> Optional[Dict[str, Any]]:
    """Parse HOMO/LUMO energies from ORCA output (uses last orbital block)."""
    if not output_file.exists():
        return None

    try:
        with open(output_file, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()

        orbital_line_re = re.compile(
            r'^\s*(\d+)\s+([\d.]+)\s+([+\-]?\d*\.?\d+(?:[Ee][+\-]?\d+)?)\s+([+\-]?\d*\.?\d+(?:[Ee][+\-]?\d+)?)'
        )
        header_re = re.compile(r'^\s*NO\s+OCC\s+E\(Eh\)\s+E\(eV\)', re.IGNORECASE)

        blocks: list[list[Dict[str, Any]]] = []
        current_block: list[Dict[str, Any]] = []
        in_block = False
        current_spin: Optional[str] = None

        for line in lines:
            stripped = line.strip()

            if stripped.startswith("ORBITAL ENERGIES"):
                # Start new block
                in_block = True
                current_block = []
                current_spin = None
                continue

            if not in_block:
                continue

            # Track spin information if present (open-shell cases)
            upper_line = stripped.upper()
            if "SPIN UP" in upper_line:
                current_spin = "alpha"
                continue
            if "SPIN DOWN" in upper_line:
                current_spin = "beta"
                continue

            # Skip header/separator lines inside a block
            if header_re.match(line) or set(stripped) == {"-"} or stripped == "":
                continue

            match = orbital_line_re.match(line)
            if match:
                idx = int(match.group(1))
                occ = float(match.group(2))
                energy_eh = float(match.group(3))
                energy_ev = float(match.group(4))

                orbital_entry = {
                    "index": idx,
                    "occupancy": occ,
                    "energy_hartree": energy_eh,
                    "energy_eV": energy_ev
                }
                if current_spin:
                    orbital_entry["spin"] = current_spin

                current_block.append(orbital_entry)
                continue

            # End of current block detected
            if stripped.startswith("*Only") or stripped.startswith("MOLECULAR ORBITALS"):
                if current_block:
                    blocks.append(current_block)
                in_block = False
                current_block = []
                continue

            # Any other content ends the block once we have data
            if current_block:
                blocks.append(current_block)
            in_block = False
            current_block = []

        # Catch a trailing block if file ended right after data
        if in_block and current_block:
            blocks.append(current_block)

        if not blocks:
            return None

        orbital_list = blocks[-1]
        orbitals = []
        homo_energy = None
        lumo_energy = None

        for entry in orbital_list:
            orbitals.append(entry)
            occ = entry["occupancy"]
            energy_ev = entry["energy_eV"]

            if occ > 1e-3:
                homo_energy = energy_ev
            elif lumo_energy is None:
                lumo_energy = energy_ev

        if homo_energy is not None and lumo_energy is not None:
            return {
                "homo_eV": homo_energy,
                "lumo_eV": lumo_energy,
                "gap_eV": lumo_energy - homo_energy,
                "orbital_list": orbitals
            }

    except Exception as e:
        logger.error(f"Error parsing orbitals from {output_file}: {e}")

    return None


def parse_occupier_folder(folder: Path) -> Optional[Dict[str, Any]]:
    """Parse a single OCCUPIER folder (multiplicity scans etc.)."""
    txt = folder / "OCCUPIER.txt"
    if not txt.exists():
        return None

    try:
        with open(txt, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        data: Dict[str, Any] = {
            "source_folder": folder.name,
            "method": None,
            "charge": None,
            "preferred_index": None,
            "electron_number": None,
            "lowest_energy_hartree": None,
            "entries": []
        }

        method_match = next((re.search(r"^Method:\s*(.+)", ln) for ln in lines if ln.strip().startswith("Method:")), None)
        if method_match:
            data["method"] = method_match.group(1).strip()

        for ln in lines:
            if "Charge:" in ln:
                m = re.search(r"Charge:\s*([-]?\d+)", ln)
                if m:
                    data["charge"] = int(m.group(1))
            if "Electron number" in ln:
                m = re.search(r"Electron number:\s*([^\s)]+)", ln)
                if m:
                    data["electron_number"] = m.group(1).rstrip(')')
            if "Preferred Index" in ln:
                m = re.search(r"Preferred Index:\s*(\d+)", ln)
                if m:
                    data["preferred_index"] = int(m.group(1))
            if "LOWEST FINAL SINGLE POINT ENERGY" in ln:
                m = re.search(r'LOWEST FINAL SINGLE POINT ENERGY:\s*([-\d.]+)', ln)
                if m:
                    data["lowest_energy_hartree"] = float(m.group(1))

        block: Dict[str, Any] = {}
        for ln in lines:
            energy_match = re.search(r"FINAL SINGLE POINT ENERGY \((\d+)\)\s*=\s*([-\d.]+)\s*\(H\)(.*)", ln)
            if energy_match:
                if block:
                    data["entries"].append(block)
                    block = {}
                idx = int(energy_match.group(1))
                energy = float(energy_match.group(2))
                tail = energy_match.group(3) or ""
                block = {
                    "index": idx,
                    "energy_hartree": energy,
                    "preferred": "<--" in tail,
                }
                continue

            mult_match = re.search(r"multiplicity\s+(\d+)", ln, re.IGNORECASE)
            if mult_match and block:
                block["multiplicity"] = int(mult_match.group(1))
                continue

            spin_match = re.search(r"Spin Contamination.*?:\s*([-\d.]+|N/A)", ln, re.IGNORECASE)
            if spin_match and block:
                val = spin_match.group(1)
                block["spin_contamination"] = None if val.upper() == "N/A" else float(val)
                continue

            unpaired_match = re.search(r"Unpaired e.*?:\s*([-\d]+)\s*\|\s*([-\d]+)", ln)
            if unpaired_match and block:
                block["unpaired_alpha"] = int(unpaired_match.group(1))
                block["unpaired_beta"] = int(unpaired_match.group(2))
                continue

        if block:
            data["entries"].append(block)

        if data["entries"] and data["preferred_index"] is None:
            # Fallback: mark first preferred flag
            for entry in data["entries"]:
                if entry.get("preferred"):
                    data["preferred_index"] = entry["index"]
                    break

        return data
    except Exception as e:
        logger.error(f"Error parsing OCCUPIER folder {folder}: {e}")
        return None


def _calculate_total_isc_rate(ms_components: Dict[str, Dict[str, Any]]) -> float:
    """Calculate total ISC rate from Ms components with symmetry approximation.

    - 3 components (Ms=-1,0,+1): sum all three
    - 2 components (Ms=0,±1): rate(Ms=0) + 2×rate(Ms=±1)
    - 1 component (Ms=0 only): 3×rate(Ms=0)

    Args:
        ms_components: Dictionary of Ms component data (ms_0, ms_p1, ms_m1)

    Returns:
        Total ISC rate (s^-1)
    """
    # Get rates for each Ms component
    rate_ms0 = ms_components.get("ms_0", {}).get("rate_s1") or 0
    rate_msp1 = ms_components.get("ms_p1", {}).get("rate_s1") or 0
    rate_msm1 = ms_components.get("ms_m1", {}).get("rate_s1") or 0

    # Check which components are present
    has_ms0 = "ms_0" in ms_components
    has_msp1 = "ms_p1" in ms_components
    has_msm1 = "ms_m1" in ms_components
    num_components = sum([has_ms0, has_msp1, has_msm1])

    if num_components == 1:
        # Only 1 component: use 3× that value (approximates all 3 Ms components)
        if has_ms0:
            return 3 * rate_ms0
        elif has_msp1:
            return 3 * rate_msp1
        else:  # has_msm1
            return 3 * rate_msm1
    elif num_components == 2:
        # 2 components: use symmetry approximation
        if has_msp1 and not has_msm1:
            # Missing Ms=-1, use 2× Ms=+1
            return rate_ms0 + 2 * rate_msp1
        elif has_msm1 and not has_msp1:
            # Missing Ms=+1, use 2× Ms=-1
            return rate_ms0 + 2 * rate_msm1
        else:
            # Has both Ms=±1 but no Ms=0 (unusual), just sum
            return rate_ms0 + rate_msp1 + rate_msm1
    else:
        # All 3 present: simple sum
        return rate_ms0 + rate_msp1 + rate_msm1


def parse_isc_data(esd_dir: Path, state1: str, state2: str) -> Optional[Dict[str, Any]]:
    """Parse ISC data from ISC output files."""
    isc_data = {"ms_components": {}}

    # Parse all ms components
    for ms in ["ms0", "msp1", "msm1"]:
        isc_file = esd_dir / f"{state1}_{state2}_ISC_{ms}.out"
        if not isc_file.exists():
            continue

        try:
            with open(isc_file, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()

            # Extract SOC (complex), rate, and FC/HT percentages
            soc_match = re.search(
                r'Reference SOC \(Re and Im\):\s*([-\d.eE+]+),\s*([-\d.eE+]+)',
                content
            )
            rate_match = re.search(
                r'(?:k_ISC|ISC rate constant is)\s*=?\s*([\d.eE+-]+)\s*s-?1',
                content
            )
            fc_ht_match = re.search(
                r'with\s+([\d.]+)\s+from\s+FC\s+and\s+([\d.]+)\s+from\s+HT',
                content,
                re.IGNORECASE
            )
            delta_match = re.search(
                r'0-0 energy difference:\s*([-\d.]+)\s*cm-1',
                content
            )
            temp_match = re.search(r'Temperature used:\s*([-\d.]+)\s*K', content)

            ms_key = ms.replace("ms", "ms_")
            soc_re = float(soc_match.group(1)) if soc_match else None
            soc_im = float(soc_match.group(2)) if soc_match else None
            soc_abs = None
            if soc_re is not None and soc_im is not None:
                soc_abs = (soc_re ** 2 + soc_im ** 2) ** 0.5
            isc_data["ms_components"][ms_key] = {
                "soc_re_cm1": soc_re,
                "soc_im_cm1": soc_im,
                "soc_abs_cm1": soc_abs,
                "rate_s1": float(rate_match.group(1)) if rate_match else None,
                "delta_E_cm1": float(delta_match.group(1)) if delta_match else None,
                "temperature_K": float(temp_match.group(1)) if temp_match else None,
                "fc_percent": float(fc_ht_match.group(1)) if fc_ht_match else None,
                "ht_percent": float(fc_ht_match.group(2)) if fc_ht_match else None,
                "source_file": isc_file.name
            }

        except Exception as e:
            logger.error(f"Error parsing ISC file {isc_file}: {e}")

    # Calculate total rate
    if isc_data["ms_components"]:
        isc_data["total_rate_s1"] = _calculate_total_isc_rate(isc_data["ms_components"])
        return isc_data

    return None


def parse_ic_data(esd_dir: Path, state1: str, state2: str) -> Optional[Dict[str, Any]]:
    """Parse IC data from IC output files."""
    ic_file = esd_dir / f"{state1}_{state2}_IC.out"

    if not ic_file.exists():
        return None

    try:
        with open(ic_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        # Extract IC rate
        rate_match = re.search(
            r'(?:k_IC|internal conversion rate constant is)\s*=?\s*([\d.eE+-]+)\s*s-?1',
            content
        )
        delta_e_match = re.search(r'0-0 energy difference:\s*([-\d.]+)\s*cm-1', content)
        temp_match = re.search(r'Temperature used:\s*([-\d.]+)\s*K', content)

        return {
            "rate_s1": float(rate_match.group(1)) if rate_match else None,
            "delta_E_cm1": float(delta_e_match.group(1)) if delta_e_match else None,
            "temperature_K": float(temp_match.group(1)) if temp_match else None,
            "source_file": ic_file.name
        }

    except Exception as e:
        logger.error(f"Error parsing IC file {ic_file}: {e}")

    return None


def parse_fluor_data(esd_dir: Path, state1: str = "S1", state2: str = "S0") -> Optional[Dict[str, Any]]:
    """Parse fluorescence rate data from ESD(FLUOR) output."""
    fluor_file = esd_dir / f"{state1}_{state2}_FLUOR.out"
    if not fluor_file.exists():
        return None
    try:
        with fluor_file.open("r", encoding="utf-8", errors="ignore") as f:
            content = f.read()

        rate_match = re.search(
            r'(?:k[_\s-]*f|k[_\s-]*fluor|fluorescence\s+rate\s+constant\s+(?:is)?|calculated\s+fluorescence\s+rate\s+constant\s+is)\s*=?\s*([\d.eE+-]+)\s*s-?1',
            content,
            flags=re.IGNORECASE,
        )
        temp_match = re.search(r'Temperature used:\s*([-\d.]+)\s*K', content)
        delta_e_match = re.search(r'0-0 energy difference:\s*([-\d.]+)\s*cm-1', content)

        return {
            "rate_s1": float(rate_match.group(1)) if rate_match else None,
            "temperature_K": float(temp_match.group(1)) if temp_match else None,
            "delta_E_cm1": float(delta_e_match.group(1)) if delta_e_match else None,
            "source_file": fluor_file.name,
        }
    except Exception as e:
        logger.error(f"Error parsing FLUOR file {fluor_file}: {e}")
        return None


def parse_phosp_data(esd_dir: Path, state1: str = "T1", state2: str = "S0") -> Optional[Dict[str, Any]]:
    """Parse phosphorescence rates from ESD(PHOSP) output.

    The PHOSP workflow may run multiple IROOT subjobs in a single output via $new_job.
    We parse all phosphorescence rate constants and return k1/k2/k3 plus arithmetic mean.
    """
    phosp_file = esd_dir / f"{state1}_{state2}_PHOSP.out"
    if not phosp_file.exists():
        return None
    try:
        with phosp_file.open("r", encoding="utf-8", errors="ignore") as f:
            content = f.read()

        # Determine expected IROOTs from the corresponding input (if present)
        expected_iroots: list[str] = []
        inp_file = esd_dir / f"{state1}_{state2}_PHOSP.inp"
        if inp_file.exists():
            try:
                inp_text = inp_file.read_text(encoding="utf-8", errors="ignore")
                expected_iroots = re.findall(r"^\s*IROOT\s+(\d+)\s*$", inp_text, flags=re.IGNORECASE | re.MULTILINE)
            except Exception:
                expected_iroots = []

        # Collect all occurrences (often one per subjob / IROOT)
        matches = re.findall(
            r'(?:k[_\s-]*p|k[_\s-]*phosp|phosphorescence\s+rate\s+constant\s+(?:is)?|calculated\s+phosphorescence\s+rate\s+constant\s+is)\s*=?\s*([\d.eE+-]+)\s*s-?1',
            content,
            flags=re.IGNORECASE,
        )
        rates = [float(x) for x in matches] if matches else []

        temp_match = re.search(r'Temperature used:\s*([-\d.]+)\s*K', content)
        delta_e_match = re.search(r'0-0 energy difference:\s*([-\d.]+)\s*cm-1', content)

        # Map in order of appearance onto expected IROOTs if we know them,
        # otherwise default to 1..N.
        if expected_iroots:
            iroot_rates = {iroot: None for iroot in expected_iroots}
            for idx, rate in enumerate(rates):
                if idx >= len(expected_iroots):
                    break
                iroot_rates[expected_iroots[idx]] = rate
        else:
            iroot_rates = {str(i + 1): rates[i] for i in range(len(rates))}
        mean_rate = sum(rates) / len(rates) if rates else None

        return {
            "iroot_rates_s1": iroot_rates,
            "rate_mean_s1": mean_rate,
            "temperature_K": float(temp_match.group(1)) if temp_match else None,
            "delta_E_cm1": float(delta_e_match.group(1)) if delta_e_match else None,
            "source_file": phosp_file.name,
        }
    except Exception as e:
        logger.error(f"Error parsing PHOSP file {phosp_file}: {e}")
        return None


def parse_esd_summary(project_dir: Path) -> Dict[str, Any]:
    """Parse summary data from ESD.txt if present (rates, SOC, FC/HT)."""
    esd_txt = project_dir / "ESD.txt"
    if not esd_txt.exists():
        return {}

    summary_data: Dict[str, Any] = {"isc": {}, "ic": {}}

    try:
        with open(esd_txt, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()

        for line in lines:
            ms_match = re.match(r'^\s*([ST]\d)>([ST]\d)\(Ms=([+\-]?\d)\):\s*([-\deE+.]+)', line)
            if ms_match:
                s1, s2, ms_val, rate = ms_match.groups()
                temp = re.search(r'T=([-\d.]+)\s*K', line)
                delta = re.search(r'0-0=([-\d.]+)', line)
                soc = re.search(r'SOC=([^\s,]+)\+i([^\s,]+)', line)
                fc_ht = re.search(r'FC=([-\d.]+)%,\s*HT=([-\d.]+)%', line)
                key = f"{s1}_{s2}"
                ms_key = {"0": "ms_0", "+1": "ms_p1", "-1": "ms_m1"}.get(ms_val, f"ms_{ms_val}")
                entry = summary_data["isc"].setdefault(key, {"ms_components": {}})
                soc_re_f = float(soc.group(1)) if soc else None
                soc_im_f = float(soc.group(2)) if soc else None
                soc_abs = None
                if soc_re_f is not None and soc_im_f is not None:
                    soc_abs = (soc_re_f ** 2 + soc_im_f ** 2) ** 0.5
                entry["ms_components"][ms_key] = {
                    "soc_re_cm1": soc_re_f,
                    "soc_im_cm1": soc_im_f,
                    "soc_abs_cm1": soc_abs,
                    "rate_s1": float(rate),
                    "delta_E_cm1": float(delta.group(1)) if delta else None,
                    "temperature_K": float(temp.group(1)) if temp else None,
                    "fc_percent": float(fc_ht.group(1)) if fc_ht else None,
                    "ht_percent": float(fc_ht.group(2)) if fc_ht else None,
                    "source_file": "ESD.txt"
                }
                continue

            total_match = re.match(r'^\s*([ST]\d)>([ST]\d)\s*\(total\):\s*([-\deE+.]+)\s*s\^-?1', line)
            if total_match:
                s1, s2, rate = total_match.groups()
                key = f"{s1}_{s2}"
                entry = summary_data["isc"].setdefault(key, {"ms_components": {}})
                entry["total_rate_s1"] = float(rate)
                entry["source_file"] = "ESD.txt"
                continue

            ic_match = re.match(
                r'^\s*([ST]\d)>([ST]\d):\s*([-\deE+.]+)\s*s\^-?1\s*\(T=([-\d.]+)\s*K,\s*.*?0-0=([-\d.]+)',
                line
            )
            if ic_match:
                s1, s2, rate, temp, delta = ic_match.groups()
                key = f"{s1}_{s2}"
                summary_data["ic"][key] = {
                    "rate_s1": float(rate),
                    "temperature_K": float(temp),
                    "delta_E_cm1": float(delta),
                    "source_file": "ESD.txt"
                }

    except Exception as e:
        logger.error(f"Error parsing ESD.txt summary: {e}")

    return summary_data


def collect_esd_data(project_dir: Path) -> Dict[str, Any]:
    """
    Collect all ESD calculation data from a project directory.

    Args:
        project_dir: Path to project directory containing initial.inp, ESD/, OCCUPIER/, etc.

    Returns:
        Complete data structure as dictionary (ready for JSON)
    """
    from delfin.esd_input_generator import _resolve_state_filename

    logger.info(f"Collecting ESD data from {project_dir}")

    esd_dir = project_dir / "ESD"
    occupier_dir = project_dir / "OCCUPIER"
    summary_data = parse_esd_summary(project_dir)
    control_flags = parse_control_flags(project_dir)

    # Read ESD_modus from CONTROL.txt to determine correct filenames for hybrid1 mode
    esd_mode = 'tddft'  # default
    control_path = project_dir / "CONTROL.txt"
    if control_path.exists():
        try:
            content = control_path.read_text(encoding="utf-8", errors="ignore")
            mode_match = re.search(r'^\s*ESD_modus\s*=\s*(.+)', content, flags=re.IGNORECASE | re.MULTILINE)
            if mode_match:
                esd_mode = mode_match.group(1).strip().lower()
                if "|" in esd_mode:
                    esd_mode = esd_mode.split("|")[0].strip()
        except Exception as exc:  # noqa: BLE001
            logger.warning("Failed to parse ESD_modus from CONTROL.txt: %s", exc)

    # Parse DELFIN.txt for redox potentials and summary info
    delfin_summary = parse_delfin_txt(project_dir)

    # Get git commit info for reproducibility tracking
    git_commit = get_git_commit_info()

    data = {
        "git_commit": git_commit,
        "metadata": parse_initial_inp(project_dir),
        "ground_state_S0": {},
        "excited_states": {},
        "oxidized_state": {},
        "oxidized_states": {},
        "reduced_state": {},
        "reduced_states": {},
        "vibrational_frequencies": {},
        "emission": {},
        "intersystem_crossing": {},
        "internal_conversion": {},
        "fluorescence_rates": {},
        "phosphorescence_rates": {},
        "occupier": {},
        "photophysical_rates": {},
        "delfin_summary": delfin_summary,
        "control_flags": control_flags,
    }

    # Parse S0
    s0_out = esd_dir / "S0.out"
    if s0_out.exists():
        logger.info("Parsing S0 ground state")
        data["ground_state_S0"]["optimization"] = {
            "converged": True,
            "geometry_file": "S0.xyz"
        }

        # Parse charge and multiplicity from S0 input file in ESD directory
        s0_inp = esd_dir / "S0.inp"
        cm_data = parse_charge_multiplicity(s0_inp)
        if cm_data:
            data["ground_state_S0"]["optimization"].update(cm_data)

        scf_data = parse_scf_energy(s0_out)
        if scf_data:
            data["ground_state_S0"]["optimization"].update(scf_data)

        data["ground_state_S0"]["thermochemistry"] = parse_thermochemistry(s0_out)
        data["ground_state_S0"]["orbitals"] = parse_orbitals(s0_out)
        data["ground_state_S0"]["dipole_moment"] = parse_dipole_moment(s0_out)

        # Parse polarizability if available
        from delfin.parser import parse_polarizability
        polarizability = parse_polarizability(s0_out)
        if polarizability:
            data["ground_state_S0"]["polarizability"] = polarizability

        # Parse hyperpolarizability if available
        from delfin.parser import parse_hyperpolarizability, calculate_beta_properties
        beta_tensor = parse_hyperpolarizability(s0_out)
        if beta_tensor:
            dipole = data["ground_state_S0"]["dipole_moment"] or {}
            dipole_x = dipole.get("x_au", 0.0)
            dipole_y = dipole.get("y_au", 0.0)
            dipole_z = dipole.get("z_au", 0.0)
            beta_props = calculate_beta_properties(beta_tensor, dipole_x, dipole_y, dipole_z)
            data["ground_state_S0"]["hyperpolarizability"] = {
                "tensor_au": beta_tensor,
                **beta_props
            }

        # Parse absorption spectrum from S0
        # Try S0_TDDFT.out first, fall back to S0.out if it doesn't exist
        s0_tddft = esd_dir / "S0_TDDFT.out"
        if not s0_tddft.exists():
            s0_tddft = s0_out if s0_out.exists() else None

        if s0_tddft and s0_tddft.exists():
            transitions = parse_absorption_spectrum(s0_tddft)
            data["ground_state_S0"]["tddft_absorption"] = {
                "comment": f"Parsed from {s0_tddft.name}",
                "transitions": [
                    {
                        "from_state": t.from_state,
                        "to_state": t.to_state,
                        "energy_eV": t.energy_ev,
                        "wavelength_nm": t.wavelength_nm,
                        "oscillator_strength": t.fosc,
                        "excitations": t.excitations if t.excitations else [],
                        "homo_number": t.homo_number
                    }
                    for t in transitions
                ]
            }

    # Parse oxidized and reduced states
    def parse_charged_series(prefix: str, target_states: Dict[str, Any], legacy_key: str):
        # Parse up to three steps (1, 2, 3) if they exist
        for step in [1, 2, 3]:
            base = f"{prefix}_step_{step}"
            out_file = project_dir / f"{base}.out"
            if not out_file.exists():
                continue

            logger.info(f"Parsing {prefix} state step {step}")
            entry = {
                "optimization": {
                    "converged": True,
                    "geometry_file": f"{base}.xyz"
                },
                "thermochemistry": parse_thermochemistry(out_file),
                "orbitals": parse_orbitals(out_file),
                "dipole_moment": parse_dipole_moment(out_file)
            }

            cm_data = parse_charge_multiplicity(project_dir / f"{base}.inp")
            if cm_data:
                entry["optimization"].update(cm_data)

            scf_data = parse_scf_energy(out_file)
            if scf_data:
                entry["optimization"].update(scf_data)

            target_states[base] = entry

        # Legacy single-key exposure for step 1
        step1_key = f"{prefix}_step_1"
        if step1_key in target_states:
            data[legacy_key] = target_states[step1_key]

    parse_charged_series("ox", data["oxidized_states"], "oxidized_state")
    parse_charged_series("red", data["reduced_states"], "reduced_state")

    # Parse excited states (S1-S6, T1-T6)
    for state in ["S1", "S2", "S3", "S4", "S5", "S6", "T1", "T2", "T3", "T4", "T5", "T6"]:
        # Resolve filename for hybrid1 mode: S1_second.out, T2_second.out, etc.
        state_out_filename = _resolve_state_filename(state, 'out', esd_mode)
        state_out = esd_dir / state_out_filename
        if not state_out.exists():
            continue

        logger.info(f"Parsing excited state {state} (from {state_out_filename})")

        # Geometry file also uses _second_deltaSCF suffix in hybrid1
        state_xyz_filename = _resolve_state_filename(state, 'xyz', esd_mode)

        state_data = {
            "optimization": {
                "converged": True,
                "geometry_file": state_xyz_filename
            },
            "thermochemistry": parse_thermochemistry(state_out),
            "dipole_moment": parse_dipole_moment(state_out),
            "orbitals": parse_orbitals(state_out)
        }

        # Parse charge and multiplicity from input file
        state_inp = esd_dir / f"{state}.inp"
        cm_data = parse_charge_multiplicity(state_inp)
        if cm_data:
            state_data["optimization"].update(cm_data)

        # Parse TDDFT keywords from the corresponding input, if available
        tddft_inp = esd_dir / f"{state}_TDDFT.inp"
        tddft_settings = parse_tddft_settings(tddft_inp if tddft_inp.exists() else state_inp)
        if tddft_settings:
            state_data["tddft_settings"] = tddft_settings

        scf_data = parse_scf_energy(state_out)
        if scf_data:
            state_data["optimization"].update(scf_data)

        # Parse emission spectrum (transitions FROM this state)
        # For hybrid1 mode: TDDFT excitations come from first_TDDFT (step 1)
        # For deltaSCF mode: TDDFT excitations come from separate _TDDFT job
        # Priority: first_TDDFT.out > _TDDFT.out > state.out
        tddft_first_out = esd_dir / f"{state}_first_TDDFT.out"
        tddft_out = esd_dir / f"{state}_TDDFT.out"
        if tddft_first_out.exists():
            tddft_source = tddft_first_out
        elif tddft_out.exists():
            tddft_source = tddft_out
        else:
            tddft_source = state_out
        transitions = parse_absorption_spectrum(tddft_source)
        if transitions:
            state_data["tddft_from_geometry"] = {
                "comment": f"Vertical transitions FROM {state} geometry (parsed from {tddft_source.name})",
                "transitions": [
                    {
                        "from_state": t.from_state,
                        "to_state": t.to_state,
                        "energy_eV": t.energy_ev,
                        "wavelength_nm": t.wavelength_nm,
                        "oscillator_strength": t.fosc,
                        "excitations": t.excitations if t.excitations else [],
                        "homo_number": t.homo_number
                    }
                    for t in transitions
                ]
            }

        data["excited_states"][state] = state_data

    # Parse OCCUPIER data for all *_OCCUPIER folders
    occupier_results: Dict[str, Any] = {}
    for folder in sorted(project_dir.glob("*_OCCUPIER")):
        parsed = parse_occupier_folder(folder)
        if parsed:
            key = folder.name.replace("_OCCUPIER", "")
            # try to attach orbitals from the matching parent out (if exists)
            parent_out = project_dir / f"{key}.out"
            parsed["orbitals"] = parse_orbitals_from_occuper(parent_out)
            occupier_results[key] = parsed
    if occupier_results:
        data["occupier"] = occupier_results

    # Parse ISC data
    isc_pairs = [("S1", "T1"), ("S1", "T2"), ("S2", "T1")]
    for state1, state2 in isc_pairs:
        isc_data = parse_isc_data(esd_dir, state1, state2)
        key = f"{state1}_{state2}"
        summary_isc = summary_data.get("isc", {}).get(key)
        if summary_isc:
            # Ensure total rate present
            if "total_rate_s1" not in summary_isc and summary_isc.get("ms_components"):
                summary_isc["total_rate_s1"] = _calculate_total_isc_rate(summary_isc["ms_components"])
            data["intersystem_crossing"][key] = summary_isc
        elif isc_data:
            data["intersystem_crossing"][key] = isc_data

    # Add any additional ISC entries from ESD.txt not covered above
    for key, entry in summary_data.get("isc", {}).items():
        if key not in data["intersystem_crossing"]:
            if "total_rate_s1" not in entry and entry.get("ms_components"):
                entry["total_rate_s1"] = _calculate_total_isc_rate(entry["ms_components"])
            data["intersystem_crossing"][key] = entry

    # Parse IC data
    ic_pairs = [("S1", "S0"), ("S2", "S1"), ("T2", "T1")]
    for state1, state2 in ic_pairs:
        ic_data = parse_ic_data(esd_dir, state1, state2)
        key = f"{state1}_{state2}"
        summary_ic = summary_data.get("ic", {}).get(key)
        if summary_ic:
            data["internal_conversion"][key] = summary_ic
        elif ic_data:
            data["internal_conversion"][key] = ic_data

    # Add any additional IC entries from ESD.txt not covered above
    for key, entry in summary_data.get("ic", {}).items():
        if key not in data["internal_conversion"]:
            data["internal_conversion"][key] = entry

    # Parse fluorescence/phosphorescence rates if present
    fluor = parse_fluor_data(esd_dir, "S1", "S0")
    if fluor:
        data["fluorescence_rates"]["S1_S0"] = fluor

    phosp = parse_phosp_data(esd_dir, "T1", "S0")
    if phosp:
        data["phosphorescence_rates"]["T1_S0"] = phosp

    logger.info("ESD data collection complete")
    return data


def save_esd_data_json(project_dir: Path, output_json: Optional[Path] = None) -> Path:
    """
    Collect all ESD data and save to JSON file.

    Args:
        project_dir: Path to project directory
        output_json: Optional output JSON path (default: project_dir/DELFIN_Data.json)

    Returns:
        Path to saved JSON file
    """
    data = collect_esd_data(project_dir)

    if output_json is None:
        output_json = project_dir / "DELFIN_Data.json"

    with open(output_json, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    logger.info(f"Saved DELFIN data to {output_json}")
    return output_json
