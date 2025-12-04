"""
DELFIN Report Assistant

AI-powered scientific report generation based on extracted calculation data.
Generates factual, non-hallucinating reports in academic style.
"""

from __future__ import annotations
from pathlib import Path
from typing import Optional

from .base_provider import BaseLLMProvider, Message
from .report_parser import DELFINReportData, ReportParser


# System prompt for report generation - enforces factual reporting
REPORT_SYSTEM_PROMPT = """You are a scientific report generator for computational chemistry calculations performed with DELFIN.

Your task is to generate a DETAILED, publication-quality summary of DELFIN calculations in academic style.

CRITICAL RULES:
1. NEVER hallucinate or invent data
2. ONLY report values explicitly provided in the input data
3. If a value is None or missing, DO NOT mention it
4. Use precise scientific language with full technical detail
5. Report ALL energies, wavelengths, and values with appropriate precision
6. ALWAYS mention DELFIN as the automation software
7. Be comprehensive and thorough - this is for scientific publications

Required structure:
1. CONFORMER SEARCH (if available): "A total of X conformers were identified following global geometry optimization using the [METHOD] method in combination with the [ALGORITHM] algorithm."

2. SOFTWARE: "All calculations were carried out employing the ORCA, XTB, and DELFIN software packages..."

3. METHOD DETAILS: Full specification including:
   - Functional (e.g., PBE0)
   - Basis set (e.g., def2-SVP)
   - Auxiliary basis (e.g., def2/J)
   - Dispersion correction (e.g., D4)
   - Integration method (e.g., RIJCOSX)
   - Solvation model (e.g., CPCM(MeCN))
   - Any special basis sets for metals

4. ENERGIES: Report with full precision
   - Final single-point energy in both eV and Hartree
   - HOMO energy in eV
   - LUMO energy in eV
   - HOMO-LUMO gap in eV

5. VIBRATIONAL ANALYSIS:
   - List ALL most intense vibrational modes with frequencies in cm⁻¹
   - Describe the character of the modes (in-plane, out-of-plane, symmetric, asymmetric, stretching, bending)
   - State explicitly if no imaginary frequencies were observed
   - If imaginary frequencies exist, report them

6. EXCITED STATES (if available):
   - Total number of electronic transitions calculated
   - Number of singlet and triplet states
   - ALL intense absorption peaks with wavelengths in nm and oscillator strengths
   - S₀ → S₁ vertical excitation energy (eV and nm)
   - S₀ → T₁ vertical excitation energy (eV and nm)
   - S₁ → S₀ emission energy (eV and nm)
   - T₁ → S₀ phosphorescence energy (eV and nm)
   - E00 energy (eV)

7. REDOX POTENTIALS (if available):
   - Report all E_ox, E_ox_2, E_red, E_red_2 values
   - Include reference electrode (vs. Fc+/Fc or vs. SCE)

8. CALCULATION STATUS:
   - Note SCF and geometry convergence
   - If errors or warnings occurred, mention them

Style requirements:
- Use past tense throughout
- Use subscripts: S₀, S₁, T₁, etc.
- Include all units: eV, nm, cm⁻¹, V
- Write in paragraph form, not bullet points
- Length: 300-500 words for complete dataset
- Professional academic tone suitable for publication in J. Am. Chem. Soc., Inorg. Chem., etc.
"""


# Template for report structure
REPORT_TEMPLATE = """Generate a comprehensive, publication-quality scientific report based on the following calculation data.

CALCULATION DATA:
{data_section}

STRUCTURE YOUR REPORT AS FOLLOWS:

PARAGRAPH 1 - Method and Software:
- Mention conformer search if data available
- State: "All calculations were carried out employing the ORCA, XTB, and DELFIN software packages"
- Include full computational method with all technical details (functional, basis set, auxiliary basis, dispersion, integration method, solvation)
- Mention any special basis sets for specific atoms

PARAGRAPH 2 - Energies and Electronic Structure:
- Report final single-point energy in BOTH eV and Hartree with full precision
- Report HOMO energy, LUMO energy, and HOMO-LUMO gap
- Use proper formatting: "The final single-point energy of the optimized structure was determined to be −X.XX eV (−X.XXXXXX Eh)."

PARAGRAPH 3 - Vibrational Analysis:
- List ALL most intense vibrational modes with frequencies
- Describe the character of each mode (e.g., "in-plane symmetric vibration", "C=O stretching")
- Explicitly state: "No imaginary (negative) frequencies were observed, confirming that the optimized structure corresponds to a local minimum on the potential energy surface" OR report imaginary frequencies if present

PARAGRAPH 4 - Excited States (if available):
- Total number of transitions and state types (singlet/triplet)
- ALL intense absorption peaks with wavelengths in nm
- S₀ → S₁ and S₀ → T₁ vertical excitation energies with both eV and nm
- Emission energies (S₁ → S₀ fluorescence and T₁ → S₀ phosphorescence)
- E00 energy if available

PARAGRAPH 5 - Redox Properties (if available):
- All redox potentials with reference electrode
- Format: "E_red = −X.XXX V (vs. Fc⁺/Fc)"

PARAGRAPH 6 - Calculation Status (if errors/warnings present):
- Note convergence status
- Mention any significant errors or warnings

CRITICAL REQUIREMENTS:
- Use ALL available data from the input
- Report values with appropriate precision (energies: 2 decimals in eV, 6 decimals in Eh; potentials: 3 decimals)
- Use proper subscripts and superscripts (S₀, S₁, T₁, Fc⁺/Fc)
- Professional academic tone
- Past tense throughout
- 300-500 words for full dataset
- Do NOT invent data - only use provided values
"""


class ReportAssistant:
    """AI-powered assistant for generating scientific reports from DELFIN calculations"""

    def __init__(self, provider: BaseLLMProvider):
        """
        Initialize Report Assistant.

        Args:
            provider: LLM provider (Claude, OpenAI, etc.)
        """
        self.provider = provider

    def generate_report(
        self,
        report_data: DELFINReportData,
        output_path: Optional[Path] = None
    ) -> str:
        """
        Generate a scientific report from extracted calculation data.

        Args:
            report_data: Extracted calculation data
            output_path: Optional path to save report (default: REPORT.txt)

        Returns:
            Generated report text
        """
        print("╔══════════════════════════════════════════════════════════╗")
        print("║     DELFIN Report Generator (AI-powered)                 ║")
        print("║                                                          ║")
        print(f"║  Provider: {self.provider.provider_name:20s}                      ║")
        print(f"║  Model:    {self.provider.model:30s}          ║")
        print("╚══════════════════════════════════════════════════════════╝")
        print()

        # Format data for the AI
        data_section = self._format_data_for_ai(report_data)

        # Create conversation
        messages = [
            Message(role="system", content=REPORT_SYSTEM_PROMPT),
            Message(role="user", content=REPORT_TEMPLATE.format(data_section=data_section))
        ]

        # Generate report
        print("Generating scientific report...")
        try:
            response = self.provider.chat(
                messages=messages,
                temperature=0.0,  # Deterministic for factual reporting
            )

            # Try different keys depending on provider
            report_text = (
                response.data.get('content') or
                response.data.get('text') or
                ''
            ).strip()

            if not report_text:
                # Debug: show what we got
                print(f"DEBUG: Response data keys: {response.data.keys()}")
                print(f"DEBUG: Response data: {response.data}")
                report_text = "Error: No report generated."

            # Save to file if requested
            if output_path:
                output_path = Path(output_path)
                output_path.write_text(report_text, encoding='utf-8')
                print(f"\nReport saved to: {output_path}")

            print("\n" + "="*60)
            print("GENERATED REPORT")
            print("="*60)
            print(report_text)
            print("="*60)

            return report_text

        except Exception as e:
            error_msg = f"Error generating report: {e}"
            print(f"\n{error_msg}")
            return error_msg

    def _format_data_for_ai(self, data: DELFINReportData) -> str:
        """
        Format extracted data into a structured text for the AI.

        Args:
            data: Extracted calculation data

        Returns:
            Formatted data string
        """
        lines = []

        # Basic information
        lines.append("=== BASIC INFORMATION ===")
        if data.compound_name:
            lines.append(f"Compound: {data.compound_name}")
        if data.charge is not None:
            lines.append(f"Charge: {data.charge}")
        if data.multiplicity is not None:
            lines.append(f"Multiplicity: {data.multiplicity}")

        # Conformer data
        if data.conformers and data.conformers.count > 0:
            lines.append("\n=== CONFORMER SEARCH ===")
            lines.append(f"Number of conformers: {data.conformers.count}")
            if data.conformers.method:
                lines.append(f"Method: {data.conformers.method}")
            if data.conformers.algorithm:
                lines.append(f"Algorithm: {data.conformers.algorithm}")

        # Geometry optimization
        if data.geometry and data.geometry.method:
            lines.append("\n=== GEOMETRY OPTIMIZATION ===")
            lines.append(f"Method: {data.geometry.method}")
            if data.geometry.functional:
                lines.append(f"Functional: {data.geometry.functional}")
            if data.geometry.basis_set:
                lines.append(f"Basis set: {data.geometry.basis_set}")
            if data.geometry.dispersion:
                lines.append(f"Dispersion: {data.geometry.dispersion}")
            if data.geometry.solvation:
                lines.append(f"Solvation: {data.geometry.solvation}")

        # Software packages
        if data.software_packages:
            lines.append(f"Software: {', '.join(data.software_packages)}")

        # Final energy
        if data.geometry and data.geometry.final_energy_ev is not None:
            lines.append("\n=== ENERGIES ===")
            lines.append(f"Final single-point energy: {data.geometry.final_energy_ev:.2f} eV")
            if data.geometry.final_energy_hartree is not None:
                lines.append(f"  ({data.geometry.final_energy_hartree:.6f} Eh)")

        # Orbital data
        if data.orbitals and data.orbitals.homo_ev is not None:
            lines.append("\n=== MOLECULAR ORBITALS ===")
            lines.append(f"HOMO: {data.orbitals.homo_ev:.2f} eV")
            if data.orbitals.lumo_ev is not None:
                lines.append(f"LUMO: {data.orbitals.lumo_ev:.2f} eV")
            if data.orbitals.gap_ev is not None:
                lines.append(f"HOMO-LUMO gap: {data.orbitals.gap_ev:.2f} eV")

            # OCCUPIER preferred configuration
            if data.orbitals.preferred_multiplicity is not None:
                lines.append(f"\nPreferred electron configuration (OCCUPIER):")
                lines.append(f"  Multiplicity: {data.orbitals.preferred_multiplicity}")
                if data.orbitals.preferred_brokensym:
                    lines.append(f"  Broken symmetry: {data.orbitals.preferred_brokensym}")
                if data.orbitals.spin_contamination is not None:
                    lines.append(f"  Spin contamination: {data.orbitals.spin_contamination:.3f}")

        # Vibrational frequencies
        if data.frequencies:
            lines.append("\n=== VIBRATIONAL ANALYSIS ===")
            if data.frequencies.intense_modes:
                freqs = [f"{freq:.0f} cm⁻¹" for freq, _ in data.frequencies.intense_modes]
                lines.append(f"Most intense modes: {', '.join(freqs)}")
            if data.frequencies.has_imaginary:
                lines.append(f"Imaginary frequencies: {len(data.frequencies.imaginary_frequencies)}")
                freqs_str = ', '.join(f"{f:.1f}" for f in data.frequencies.imaginary_frequencies)
                lines.append(f"  Values: {freqs_str} cm⁻¹")
            else:
                lines.append("No imaginary frequencies observed")

        # Excited states
        if data.excited_states and data.excited_states.num_states > 0:
            lines.append("\n=== EXCITED STATE CALCULATIONS ===")
            lines.append(f"Total electronic transitions: {data.excited_states.num_states}")
            if data.excited_states.singlet_count > 0:
                lines.append(f"Singlet states: {data.excited_states.singlet_count}")
            if data.excited_states.triplet_count > 0:
                lines.append(f"Triplet states: {data.excited_states.triplet_count}")

            if data.excited_states.intense_absorptions:
                lines.append("\nIntense absorption peaks:")
                for wl, osc in data.excited_states.intense_absorptions:
                    lines.append(f"  {wl:.0f} nm (f = {osc:.3f})")

            if data.excited_states.s0_s1_ev is not None:
                lines.append(f"\nS₀ → S₁ excitation:")
                lines.append(f"  Energy: {data.excited_states.s0_s1_ev:.2f} eV")
                if data.excited_states.s0_s1_nm is not None:
                    lines.append(f"  Wavelength: {data.excited_states.s0_s1_nm:.0f} nm")

            if data.excited_states.s0_t1_ev is not None:
                lines.append(f"\nS₀ → T₁ excitation:")
                lines.append(f"  Energy: {data.excited_states.s0_t1_ev:.2f} eV")
                if data.excited_states.s0_t1_nm is not None:
                    lines.append(f"  Wavelength: {data.excited_states.s0_t1_nm:.0f} nm")

            if data.excited_states.s1_s0_emission_ev is not None:
                lines.append(f"\nS₁ → S₀ fluorescence:")
                lines.append(f"  Energy: {data.excited_states.s1_s0_emission_ev:.2f} eV")
                if data.excited_states.s1_s0_emission_nm is not None:
                    lines.append(f"  Wavelength: {data.excited_states.s1_s0_emission_nm:.0f} nm")

            if data.excited_states.t1_s0_phosphorescence_ev is not None:
                lines.append(f"\nT₁ → S₀ phosphorescence:")
                lines.append(f"  Energy: {data.excited_states.t1_s0_phosphorescence_ev:.2f} eV")
                if data.excited_states.t1_s0_phosphorescence_nm is not None:
                    lines.append(f"  Wavelength: {data.excited_states.t1_s0_phosphorescence_nm:.0f} nm")

            if data.excited_states.e00_ev is not None:
                lines.append(f"\nE00 energy: {data.excited_states.e00_ev:.2f} eV")

        # Redox potentials
        if data.redox and (data.redox.e_ox is not None or data.redox.e_red is not None):
            lines.append("\n=== REDOX POTENTIALS ===")
            ref = f" (vs. {data.redox.reference})" if data.redox.reference else ""
            if data.redox.e_ox is not None:
                lines.append(f"E_ox: {data.redox.e_ox:.3f} V{ref}")
            if data.redox.e_ox_2 is not None:
                lines.append(f"E_ox_2: {data.redox.e_ox_2:.3f} V{ref}")
            if data.redox.e_red is not None:
                lines.append(f"E_red: {data.redox.e_red:.3f} V{ref}")
            if data.redox.e_red_2 is not None:
                lines.append(f"E_red_2: {data.redox.e_red_2:.3f} V{ref}")

        # Calculation status (errors and warnings)
        if data.status:
            if data.status.has_errors or data.status.has_warnings:
                lines.append("\n=== CALCULATION STATUS ===")

                if data.status.scf_converged is not None:
                    lines.append(f"SCF converged: {'Yes' if data.status.scf_converged else 'No'}")
                if data.status.geometry_converged is not None:
                    lines.append(f"Geometry converged: {'Yes' if data.status.geometry_converged else 'No'}")

                if data.status.has_errors:
                    lines.append(f"\nErrors detected ({len(data.status.errors)}):")
                    for error in data.status.errors[:5]:  # Limit to 5 most important
                        lines.append(f"  - {error}")

                if data.status.has_warnings:
                    lines.append(f"\nWarnings detected ({len(data.status.warnings)}):")
                    for warning in data.status.warnings[:5]:  # Limit to 5 most important
                        lines.append(f"  - {warning}")

        return "\n".join(lines)


def generate_report_from_directory(
    working_dir: Path | str,
    provider: BaseLLMProvider,
    output_path: Optional[Path | str] = None
) -> str:
    """
    Convenience function to generate a report from a DELFIN calculation directory.

    Args:
        working_dir: Path to DELFIN calculation directory
        provider: LLM provider for report generation
        output_path: Optional path to save report (default: working_dir/REPORT.txt)

    Returns:
        Generated report text
    """
    working_dir = Path(working_dir)

    if output_path is None:
        output_path = working_dir / "REPORT.txt"
    else:
        output_path = Path(output_path)

    # Extract data
    print("Extracting calculation data...")
    report_data = ReportParser.extract_calculation_summary(working_dir)

    # Generate report
    assistant = ReportAssistant(provider)
    report = assistant.generate_report(report_data, output_path)

    return report
