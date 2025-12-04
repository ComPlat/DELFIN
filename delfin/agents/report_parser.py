"""
Report Parser for DELFIN

Extracts data from ORCA .out files, OCCUPIER.txt, and DELFIN.txt for AI report generation.
This module provides FACTUAL data extraction without interpretation.
"""

from __future__ import annotations
import re
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass, field


@dataclass
class ConformerData:
    """Data for a single conformer"""
    count: int = 0
    method: str = ""
    algorithm: str = ""
    goat_used: bool = False
    crest_used: bool = False


@dataclass
class GeometryData:
    """Geometry optimization data"""
    method: str = ""
    functional: str = ""
    basis_set: str = ""
    dispersion: str = ""
    solvation: str = ""
    final_energy_ev: Optional[float] = None
    final_energy_hartree: Optional[float] = None


@dataclass
class FrequencyData:
    """Vibrational frequency data"""
    has_imaginary: bool = False
    imaginary_frequencies: List[float] = field(default_factory=list)
    intense_modes: List[Tuple[float, str]] = field(default_factory=list)  # (freq, description)


@dataclass
class OrbitalData:
    """Molecular orbital data"""
    homo_ev: Optional[float] = None
    lumo_ev: Optional[float] = None
    gap_ev: Optional[float] = None
    preferred_multiplicity: Optional[int] = None
    preferred_brokensym: Optional[str] = None
    spin_contamination: Optional[float] = None


@dataclass
class ExcitedStateData:
    """Excited state calculation data"""
    num_states: int = 0
    method: str = ""
    singlet_count: int = 0
    triplet_count: int = 0
    intense_absorptions: List[Tuple[float, float]] = field(default_factory=list)  # (wavelength_nm, oscillator)
    s0_s1_ev: Optional[float] = None
    s0_s1_nm: Optional[float] = None
    s0_t1_ev: Optional[float] = None
    s0_t1_nm: Optional[float] = None
    s1_s0_emission_ev: Optional[float] = None
    s1_s0_emission_nm: Optional[float] = None
    t1_s0_phosphorescence_ev: Optional[float] = None
    t1_s0_phosphorescence_nm: Optional[float] = None
    e00_ev: Optional[float] = None


@dataclass
class RedoxData:
    """Redox potential data"""
    e_red: Optional[float] = None
    e_red_2: Optional[float] = None
    e_ox: Optional[float] = None
    e_ox_2: Optional[float] = None
    reference: str = ""


@dataclass
class CalculationStatus:
    """Status and errors from calculations"""
    has_errors: bool = False
    has_warnings: bool = False
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    scf_converged: Optional[bool] = None
    geometry_converged: Optional[bool] = None


@dataclass
class DELFINReportData:
    """Complete DELFIN calculation data for report generation"""
    # Basic info
    compound_name: str = ""
    charge: Optional[int] = None
    multiplicity: Optional[int] = None

    # Calculation details
    conformers: Optional[ConformerData] = None
    geometry: Optional[GeometryData] = None
    frequencies: Optional[FrequencyData] = None
    orbitals: Optional[OrbitalData] = None
    excited_states: Optional[ExcitedStateData] = None
    redox: Optional[RedoxData] = None
    status: Optional[CalculationStatus] = None

    # Software info
    software_packages: List[str] = field(default_factory=list)

    def __post_init__(self):
        """Initialize nested dataclasses if None"""
        if self.conformers is None:
            self.conformers = ConformerData()
        if self.geometry is None:
            self.geometry = GeometryData()
        if self.frequencies is None:
            self.frequencies = FrequencyData()
        if self.orbitals is None:
            self.orbitals = OrbitalData()
        if self.excited_states is None:
            self.excited_states = ExcitedStateData()
        if self.redox is None:
            self.redox = RedoxData()
        if self.status is None:
            self.status = CalculationStatus()


class ReportParser:
    """Parse DELFIN output files to extract data for report generation"""

    @staticmethod
    def parse_control_txt(file_path: Path) -> Dict[str, Any]:
        """Parse CONTROL.txt file for method and module information"""
        data = {}

        if not file_path.exists():
            return data

        content = file_path.read_text(encoding='utf-8')

        # Extract method line (contains full computational details)
        method_match = re.search(r'method\s*=\s*(.+?)(?:\n|$)', content, re.IGNORECASE)
        if method_match:
            data['method'] = method_match.group(1).strip()

        # Check which modules are enabled
        data['XTB_GOAT'] = bool(re.search(r'XTB_GOAT\s*=\s*yes', content, re.IGNORECASE))
        data['CREST'] = bool(re.search(r'CREST\s*=\s*yes', content, re.IGNORECASE))
        data['XTB_OPT'] = bool(re.search(r'XTB_OPT\s*=\s*yes', content, re.IGNORECASE))
        data['ESD_modul'] = bool(re.search(r'ESD_modul\s*=\s*yes', content, re.IGNORECASE))
        data['ESD_frequency'] = bool(re.search(r'ESD_frequency\s*=\s*yes', content, re.IGNORECASE))
        data['E_00'] = bool(re.search(r'E_00\s*=\s*yes', content, re.IGNORECASE))
        data['OCCUPIER'] = bool(re.search(r'used_method\s*=\s*OCCUPIER', content, re.IGNORECASE))

        # Extract GOAT settings if enabled
        if data['XTB_GOAT']:
            goat_iter_match = re.search(r'GOAT_iterations\s*=\s*(\d+)', content, re.IGNORECASE)
            if goat_iter_match:
                data['GOAT_iterations'] = int(goat_iter_match.group(1))

        # Extract solvation model
        solvent_match = re.search(r'solvent\s*=\s*(\w+)', content, re.IGNORECASE)
        if solvent_match:
            data['solvent'] = solvent_match.group(1)

        return data

    @staticmethod
    def parse_delfin_txt(file_path: Path) -> Dict[str, Any]:
        """Parse DELFIN.txt file for summary data"""
        data = {}

        if not file_path.exists():
            return data

        content = file_path.read_text(encoding='utf-8')

        # Extract method line
        method_match = re.search(r'Method freq:\s*(.+)', content)
        if method_match:
            data['method_line'] = method_match.group(1).strip()

        # Extract charge and multiplicity
        charge_match = re.search(r'Charge:\s*(-?\d+)', content)
        if charge_match:
            data['charge'] = int(charge_match.group(1))

        mult_match = re.search(r'Multiplicity:\s*(\d+)', content)
        if mult_match:
            data['multiplicity'] = int(mult_match.group(1))

        # Extract redox potentials - support both SCE and Fc+/Fc references
        redox_section = re.search(
            r'Calculated properties \(V vs\. (SCE|Fc\+/Fc)\):(.*?)(?:Experimental|Literature|TOTAL RUN TIME)',
            content,
            re.DOTALL
        )
        if redox_section:
            reference = redox_section.group(1)
            data['reference'] = reference
            redox_text = redox_section.group(2)

            e_red_match = re.search(r'E_red\s*=\s*([-\d.]+)', redox_text)
            if e_red_match:
                data['e_red'] = float(e_red_match.group(1))

            e_red2_match = re.search(r'E_red_2\s*=\s*([-\d.]+)', redox_text)
            if e_red2_match:
                data['e_red_2'] = float(e_red2_match.group(1))

            e_ox_match = re.search(r'E_ox\s*=\s*([-\d.]+)', redox_text)
            if e_ox_match:
                data['e_ox'] = float(e_ox_match.group(1))

            e_ox2_match = re.search(r'E_ox_2\s*=\s*([-\d.]+)', redox_text)
            if e_ox2_match:
                data['e_ox_2'] = float(e_ox2_match.group(1))

        return data

    @staticmethod
    def parse_occupier_txt(file_path: Path) -> Dict[str, Any]:
        """Parse OCCUPIER.txt file for orbital and energy data"""
        data = {}

        if not file_path.exists():
            return data

        content = file_path.read_text(encoding='utf-8')

        # Extract HOMO/LUMO energies (in eV)
        homo_match = re.search(r'HOMO.*?(-?\d+\.\d+)\s*eV', content)
        if homo_match:
            data['homo_ev'] = float(homo_match.group(1))

        lumo_match = re.search(r'LUMO.*?(-?\d+\.\d+)\s*eV', content)
        if lumo_match:
            data['lumo_ev'] = float(lumo_match.group(1))

        if 'homo_ev' in data and 'lumo_ev' in data:
            data['gap_ev'] = data['lumo_ev'] - data['homo_ev']

        # Extract preferred configuration (OCCUPIER selection)
        preferred_match = re.search(
            r'FINAL SINGLE POINT ENERGY \(1\)\s*=\s*[-\d.]+\s*\(H\)\s*<-- PREFERRED VALUE\s*\n\s*multiplicity\s*(\d+)(?:,\s*BrokenSym\s*([\d,]+))?',
            content
        )
        if preferred_match:
            data['preferred_multiplicity'] = int(preferred_match.group(1))
            if preferred_match.group(2):
                data['preferred_brokensym'] = preferred_match.group(2)

        # Extract spin contamination for preferred configuration
        spin_contam_match = re.search(
            r'FINAL SINGLE POINT ENERGY \(1\).*?<-- PREFERRED VALUE.*?Spin Contamination \(⟨S²⟩ - S\(S\+1\)\)\s*:\s*([-\d.]+|N/A)',
            content,
            re.DOTALL
        )
        if spin_contam_match and spin_contam_match.group(1) != 'N/A':
            try:
                data['spin_contamination'] = float(spin_contam_match.group(1))
            except ValueError:
                pass

        # Extract final energy - look for PREFERRED VALUE first
        energy_match = re.search(r'FINAL SINGLE POINT ENERGY \(1\)\s*=\s*([-\d.]+)\s*\(H\)\s*<-- PREFERRED VALUE', content)
        if not energy_match:
            # Fallback to any FINAL SINGLE POINT ENERGY
            energy_match = re.search(r'FINAL SINGLE POINT ENERGY\s*(?:\(\d+\))?\s*=?\s*([-\d.]+)', content)

        if energy_match:
            data['final_energy_hartree'] = float(energy_match.group(1))
            data['final_energy_ev'] = float(energy_match.group(1)) * 27.2114  # Hartree to eV

        return data

    @staticmethod
    def check_calculation_status(file_path: Path) -> Dict[str, Any]:
        """
        Check ORCA output for errors, warnings, and convergence status.

        Args:
            file_path: Path to ORCA .out file

        Returns:
            Dictionary with status information
        """
        status = {
            'has_errors': False,
            'has_warnings': False,
            'errors': [],
            'warnings': [],
            'scf_converged': None,
            'geometry_converged': None
        }

        if not file_path.exists():
            return status

        try:
            content = file_path.read_text(encoding='utf-8', errors='ignore')
        except Exception:
            return status

        # Check for ORCA errors
        if 'ORCA TERMINATED NORMALLY' not in content:
            status['has_errors'] = True
            status['errors'].append(f'Calculation did not terminate normally in {file_path.name}')

        # Check for specific error messages
        error_patterns = [
            (r'ERROR.*', 'ORCA error'),
            (r'FATAL ERROR.*', 'Fatal error'),
            (r'SCF NOT CONVERGED', 'SCF convergence failure'),
            (r'The optimization did not converge', 'Geometry optimization did not converge'),
            (r'ABORTING.*', 'Calculation aborted')
        ]

        for pattern, error_type in error_patterns:
            matches = re.finditer(pattern, content, re.IGNORECASE)
            for match in matches:
                error_msg = match.group(0).strip()
                if error_msg and error_msg not in status['errors']:
                    status['has_errors'] = True
                    status['errors'].append(f'{error_type}: {error_msg[:100]}')

        # Check for warnings
        warning_patterns = [
            (r'WARNING.*', 'Warning'),
            (r'ATTENTION.*', 'Attention'),
            (r'Note:.*convergence', 'Convergence note')
        ]

        for pattern, warning_type in warning_patterns:
            matches = re.finditer(pattern, content, re.IGNORECASE)
            for match in matches:
                warning_msg = match.group(0).strip()
                if warning_msg and warning_msg not in status['warnings']:
                    status['has_warnings'] = True
                    status['warnings'].append(f'{warning_type}: {warning_msg[:100]}')

        # Check SCF convergence
        if 'SCF CONVERGED' in content or 'The SCF converged' in content:
            status['scf_converged'] = True
        elif 'SCF NOT CONVERGED' in content:
            status['scf_converged'] = False

        # Check geometry convergence
        if 'THE OPTIMIZATION HAS CONVERGED' in content.upper():
            status['geometry_converged'] = True
        elif 'THE OPTIMIZATION DID NOT CONVERGE' in content.upper():
            status['geometry_converged'] = False

        # Limit number of errors/warnings to avoid overwhelming output
        status['errors'] = status['errors'][:10]
        status['warnings'] = status['warnings'][:10]

        return status

    @staticmethod
    def parse_orca_output(file_path: Path) -> Dict[str, Any]:
        """Parse ORCA .out file for detailed calculation data"""
        data = {}

        if not file_path.exists():
            return data

        content = file_path.read_text(encoding='utf-8', errors='ignore')

        # Extract method information
        method_match = re.search(r'!\s+(.+?)(?:\n|$)', content)
        if method_match:
            data['method_line'] = method_match.group(1).strip()

        # Extract vibrational frequencies
        freq_section = re.search(
            r'VIBRATIONAL FREQUENCIES(.*?)(?:\n\n|\Z)',
            content,
            re.DOTALL
        )
        if freq_section:
            freq_lines = freq_section.group(1).split('\n')
            frequencies = []
            imaginary = []

            for line in freq_lines:
                # Match frequency lines like "  0:         0.00 cm**-1"
                freq_match = re.search(r'\d+:\s+([-\d.]+)\s+cm\*\*-1', line)
                if freq_match:
                    freq = float(freq_match.group(1))
                    frequencies.append(freq)
                    if freq < 0:
                        imaginary.append(freq)

            data['frequencies'] = frequencies
            data['imaginary_frequencies'] = imaginary
            data['has_imaginary'] = len(imaginary) > 0

        # Extract excited states
        # Look for TD-DFT or CIS excited states
        excited_section = re.search(
            r'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS(.*?)(?:\n\n|\Z)',
            content,
            re.DOTALL
        )
        if excited_section:
            excited_text = excited_section.group(1)
            transitions = []

            # Parse transitions: State   Energy  Wavelength  fosc
            for line in excited_text.split('\n'):
                match = re.search(
                    r'(\d+)\s+([\d.]+)\s+\(([\d.]+)\s*nm\).*?f\s*=\s*([\d.]+)',
                    line
                )
                if match:
                    state_num = int(match.group(1))
                    energy_ev = float(match.group(2))
                    wavelength_nm = float(match.group(3))
                    fosc = float(match.group(4))
                    transitions.append({
                        'state': state_num,
                        'energy_ev': energy_ev,
                        'wavelength_nm': wavelength_nm,
                        'oscillator_strength': fosc
                    })

            data['transitions'] = transitions

        return data

    @staticmethod
    def parse_esd_directory(esd_dir: Path) -> Dict[str, Any]:
        """Parse ESD directory for excited state dynamics data"""
        data = {}

        if not esd_dir.exists() or not esd_dir.is_dir():
            return data

        # Look for state calculations (S0.out, S1.out, T1.out, etc.)
        state_files = list(esd_dir.glob('S*.out')) + list(esd_dir.glob('T*.out'))

        states = {}
        for state_file in state_files:
            state_name = state_file.stem  # S0, S1, T1, etc.
            state_data = ReportParser.parse_orca_output(state_file)
            if state_data:
                states[state_name] = state_data

        data['states'] = states

        # Look for transition calculations (ISC, IC)
        isc_files = list(esd_dir.glob('ISC*.out'))
        ic_files = list(esd_dir.glob('IC*.out'))

        data['isc_count'] = len(isc_files)
        data['ic_count'] = len(ic_files)

        return data

    @staticmethod
    def extract_calculation_summary(working_dir: Path) -> DELFINReportData:
        """
        Extract all relevant data from DELFIN calculation directory.

        Args:
            working_dir: Path to DELFIN calculation directory

        Returns:
            DELFINReportData object with all extracted information
        """
        report_data = DELFINReportData()

        # Parse CONTROL.txt first to get method and enabled modules
        control_txt = working_dir / 'CONTROL.txt'
        if control_txt.exists():
            control_data = ReportParser.parse_control_txt(control_txt)

            # Store module usage
            if control_data.get('XTB_GOAT'):
                report_data.conformers.goat_used = True
                report_data.conformers.algorithm = "GOAT"
            if control_data.get('CREST'):
                report_data.conformers.crest_used = True
            if control_data.get('XTB_OPT'):
                report_data.conformers.method = "GFN2-xTB"

            # Store method if not yet defined
            if not report_data.geometry.method and control_data.get('method'):
                report_data.geometry.method = control_data['method']

        # Parse DELFIN.txt
        delfin_txt = working_dir / 'DELFIN.txt'
        if delfin_txt.exists():
            delfin_data = ReportParser.parse_delfin_txt(delfin_txt)

            report_data.charge = delfin_data.get('charge')
            report_data.multiplicity = delfin_data.get('multiplicity')

            # Parse method line
            method_line = delfin_data.get('method_line', '')
            if method_line:
                report_data.geometry.method = method_line

                # Extract components
                if 'PBE0' in method_line:
                    report_data.geometry.functional = 'PBE0'
                elif 'B3LYP' in method_line:
                    report_data.geometry.functional = 'B3LYP'

                if 'def2-SVP' in method_line:
                    report_data.geometry.basis_set = 'def2-SVP'
                elif 'def2-TZVP' in method_line:
                    report_data.geometry.basis_set = 'def2-TZVP'

                if 'D4' in method_line:
                    report_data.geometry.dispersion = 'D4'
                elif 'D3' in method_line:
                    report_data.geometry.dispersion = 'D3'

                if 'CPCM' in method_line:
                    solv_match = re.search(r'CPCM\((\w+)\)', method_line)
                    if solv_match:
                        report_data.geometry.solvation = f'CPCM({solv_match.group(1)})'

            # Redox data
            if 'e_red' in delfin_data:
                report_data.redox.e_red = delfin_data['e_red']
            if 'e_red_2' in delfin_data:
                report_data.redox.e_red_2 = delfin_data['e_red_2']
            if 'e_ox' in delfin_data:
                report_data.redox.e_ox = delfin_data['e_ox']
            if 'e_ox_2' in delfin_data:
                report_data.redox.e_ox_2 = delfin_data['e_ox_2']
            if 'reference' in delfin_data:
                report_data.redox.reference = delfin_data['reference']
            else:
                report_data.redox.reference = "SCE"  # Default

        # Parse OCCUPIER.txt (initial or main)
        occupier_txt = working_dir / 'OCCUPIER.txt'
        if not occupier_txt.exists():
            occupier_txt = working_dir / 'initial_OCCUPIER' / 'OCCUPIER.txt'

        if occupier_txt.exists():
            occupier_data = ReportParser.parse_occupier_txt(occupier_txt)

            if 'homo_ev' in occupier_data:
                report_data.orbitals.homo_ev = occupier_data['homo_ev']
            if 'lumo_ev' in occupier_data:
                report_data.orbitals.lumo_ev = occupier_data['lumo_ev']
            if 'gap_ev' in occupier_data:
                report_data.orbitals.gap_ev = occupier_data['gap_ev']
            if 'preferred_multiplicity' in occupier_data:
                report_data.orbitals.preferred_multiplicity = occupier_data['preferred_multiplicity']
            if 'preferred_brokensym' in occupier_data:
                report_data.orbitals.preferred_brokensym = occupier_data['preferred_brokensym']
            if 'spin_contamination' in occupier_data:
                report_data.orbitals.spin_contamination = occupier_data['spin_contamination']
            if 'final_energy_ev' in occupier_data:
                report_data.geometry.final_energy_ev = occupier_data['final_energy_ev']
            if 'final_energy_hartree' in occupier_data:
                report_data.geometry.final_energy_hartree = occupier_data['final_energy_hartree']

        # Parse ESD directory if it exists
        esd_dir = working_dir / 'ESD'
        if esd_dir.exists():
            esd_data = ReportParser.parse_esd_directory(esd_dir)
            # TODO: Extract excited state data from ESD calculations

        # Detect software packages used
        report_data.software_packages = ['ORCA', 'DELFIN']
        if (working_dir / 'xtb.out').exists() or 'XTB' in str(working_dir):
            report_data.software_packages.append('XTB')

        # Check all .out files in directory and subdirectories for errors/warnings
        out_files = list(working_dir.rglob('*.out'))

        for out_file in out_files:
            file_status = ReportParser.check_calculation_status(out_file)

            # Aggregate errors and warnings
            for error in file_status['errors']:
                if error not in report_data.status.errors:
                    report_data.status.errors.append(error)
                    report_data.status.has_errors = True

            for warning in file_status['warnings']:
                if warning not in report_data.status.warnings:
                    report_data.status.warnings.append(warning)
                    report_data.status.has_warnings = True

            # Update convergence status (prioritize main calculation files)
            if 'initial' in out_file.name.lower():
                if file_status['scf_converged'] is not None:
                    report_data.status.scf_converged = file_status['scf_converged']
                if file_status['geometry_converged'] is not None:
                    report_data.status.geometry_converged = file_status['geometry_converged']

        return report_data
