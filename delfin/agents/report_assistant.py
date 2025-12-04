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
REPORT_SYSTEM_PROMPT = """You are an expert computational chemist writing a detailed scientific report based on DELFIN calculations.

Write as if YOU personally conducted this computational investigation and are now reporting your findings to scientific peers. Use natural scientific language with appropriate interpretation and discussion.

CRITICAL RULES - STRICT FACTUAL REQUIREMENTS:
1. NEVER hallucinate or invent data - ONLY use values explicitly provided in the data section
2. NEVER guess chemical interpretations (d-electron count, oxidation states, etc.) unless explicitly given
3. NEVER say "four unpaired electrons" if data says three
4. NEVER mention HOMO/LUMO if those values are not provided
5. NEVER describe vibrational modes if frequency data is missing
6. Write in expert style BUT stay 100% factual
7. Calculate energy differences in kcal/mol ONLY if multiple energies are provided
8. Interpret spin contamination values ONLY if they are provided
9. If data is missing, simply don't mention that aspect - DO NOT say "data was not provided"
10. Use transitional phrases but connect ONLY the facts you have

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

Writing style:
- Use past tense for what was done: "were calculated", "was determined", "proved to be"
- Use present tense for interpretations: "indicates", "suggests", "is consistent with"
- Write in flowing paragraphs with transitional phrases: "Furthermore,", "Notably,", "This finding suggests..."
- Use subscripts: S₀, S₁, T₁, etc.
- Include all units: eV, nm, cm⁻¹, V, kcal/mol
- Add brief chemical interpretations and significance
- Connect findings logically between paragraphs
- Length: 400-600 words for complete dataset
- Tone: Expert computational chemist explaining their findings to scientific peers
- Think: "How would an experienced researcher discuss these results in a group meeting?"
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

PARAGRAPH 2 - Electronic Configuration and Spin State Analysis:
Write this as your expert analysis of the spin state investigation:
- Start with: "To determine the ground electronic state, the OCCUPIER module systematically evaluated X spin multiplicities..."
- Report energies for ALL tested multiplicities
- Calculate and report energy differences in kcal/mol (1 Hartree = 627.5 kcal/mol)
- Discuss which state was selected: "The [multiplicity] state emerged as the energetically preferred configuration, lying X.X kcal/mol below the [next closest state]"
- Interpret spin contamination: values near 0 indicate "excellent spin purity" or "minimal multiconfigurational character"
- Report Boltzmann populations and interpret: "At room temperature, the [state] represents >99% of the ensemble, confirming its dominance"
- DO NOT guess metal identity or oxidation state unless explicitly provided
- Example: "To determine the ground electronic state, OCCUPIER systematically evaluated three spin multiplicities (M = 2, 4, and 6). The quartet state proved to be the most stable, with an energy of −3056.472 Eh. This configuration lies 6.9 kcal/mol below the doublet state (−3056.461 Eh, calculated as (−3056.461 − (−3056.472)) × 627.5 kcal/mol/Eh) and 55.2 kcal/mol below the sextet (−3056.384 Eh). The minimal spin contamination (⟨S²⟩ − S(S+1) = 0.012) confirms excellent spin purity for this state. Boltzmann analysis at 298.15 K indicates that the quartet state comprises 99.8% of the thermal ensemble, with the doublet contributing 0.2% and the sextet being essentially unpopulated. The three unpaired electrons in the quartet configuration indicate an open-shell electronic structure."

PARAGRAPH 3 - Energies and Orbitals:
- Report final single-point energy in BOTH eV and Hartree with full precision
- Report HOMO energy, LUMO energy, and HOMO-LUMO gap
- Use proper formatting: "The final single-point energy of the optimized structure was determined to be −X.XX eV (−X.XXXXXX Eh)."

PARAGRAPH 4 - Vibrational Analysis:
Write this with structural interpretation:
- Start with validation: "Frequency analysis confirmed the validity of the optimized geometry"
- If no imaginary modes: "with no imaginary frequencies detected, indicating a true minimum on the potential energy surface"
- List most intense modes with chemical interpretation
- Example: "The most intense vibrational modes were observed at 1486, 1634, and 1672 cm⁻¹, corresponding to aromatic C=C stretching vibrations, and at 3144 cm⁻¹, assigned to aromatic C−H stretching modes"
- Connect to structure when possible: "The absence of low-frequency modes below 100 cm⁻¹ suggests a relatively rigid molecular framework"

PARAGRAPH 5 - Excited States (if available):
- Total number of transitions and state types (singlet/triplet)
- ALL intense absorption peaks with wavelengths in nm
- S₀ → S₁ and S₀ → T₁ vertical excitation energies with both eV and nm
- Emission energies (S₁ → S₀ fluorescence and T₁ → S₀ phosphorescence)
- E00 energy if available

PARAGRAPH 6 - Redox-Induced Electronic Structure Changes (if OCCUPIER steps available):
This is CRITICAL - discuss how electronic configuration CHANGES with oxidation state:
- Compare initial vs red_step_1 vs red_step_2 (or ox_step_1, etc.)
- Report charge change: "Upon first reduction (initial +2 → red_step_1 +1)..."
- Report multiplicity changes: "The preferred multiplicity decreased from M=4 to M=3"
- Report unpaired electron changes: "The number of unpaired electrons changed from 3 to 2"
- Discuss energy differences between steps in kcal/mol
- Compare spin contamination across steps
- Interpret: "This multiplicity change suggests the added electron occupies a metal d-orbital..."
- Example: "Upon first reduction from the +2 to +1 oxidation state, OCCUPIER analysis reveals a striking change in electronic structure. The preferred multiplicity shifts from M = 4 (quartet, 3 unpaired electrons, E = −3056.472 Eh) to M = 3 (triplet, 2 unpaired electrons, E = −3056.XXX Eh). The reduction thus stabilizes a lower-spin configuration, with the added electron likely pairing with one of the metal d-electrons. The second reduction to the neutral species further reduces the multiplicity to M = 2 (doublet, 1 unpaired electron), indicating progressive spin-pairing upon sequential electron addition."

PARAGRAPH 7 - Redox Potentials (if available):
Report calculated potentials and connect to electronic structure changes:
- E_red, E_red_2, E_ox, E_ox_2 with reference electrode
- Connect to OCCUPIER analysis if available: "The calculated E_red = −1.330 V corresponds to the quartet→triplet transition discussed above"
- Discuss separation between successive redox events

PARAGRAPH 8 - Calculation Status (if errors/warnings present):
- Note convergence status
- Mention any significant errors or warnings

CRITICAL REQUIREMENTS:
- Use ALL available data from the input
- Report values with appropriate precision (energies: 2 decimals in eV, 6 decimals in Eh; potentials: 3 decimals)
- Use proper subscripts and superscripts (S₀, S₁, T₁, Fc⁺/Fc)
- Professional academic tone with expert interpretation
- Past tense for actions, present for interpretations
- 400-600 words for full dataset

ABSOLUTE PROHIBITIONS:
- Do NOT invent any numbers or values
- Do NOT guess electronic configurations or oxidation states
- Do NOT mention data that wasn't provided
- Do NOT contradict the provided data (e.g., if it says 3 unpaired electrons, never say 4)
- Do NOT add "data not available" disclaimers - just skip missing sections entirely
- If you calculate something (like kcal/mol), show your work using the provided Hartree values
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
        if data.conformers and (data.conformers.count > 0 or data.conformers.goat_used or data.conformers.crest_used):
            lines.append("\n=== CONFORMER SEARCH ===")
            if data.conformers.count > 0:
                lines.append(f"Number of conformers: {data.conformers.count}")
            if data.conformers.method:
                lines.append(f"Method: {data.conformers.method}")
            if data.conformers.goat_used:
                lines.append(f"Algorithm: GOAT (Global Optimization with Ansatz Trees)")
            if data.conformers.crest_used:
                lines.append("CREST conformer ensemble generation: enabled")

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

            # OCCUPIER preferred configuration and all tests
            if data.orbitals.preferred_multiplicity is not None:
                lines.append(f"\nOCCUPIER electron configuration analysis:")
                lines.append(f"  Preferred multiplicity: {data.orbitals.preferred_multiplicity}")
                if data.orbitals.unpaired_electrons is not None:
                    lines.append(f"  Unpaired electrons: {data.orbitals.unpaired_electrons}")
                if data.orbitals.preferred_brokensym:
                    lines.append(f"  Broken symmetry: {data.orbitals.preferred_brokensym}")
                if data.orbitals.spin_contamination is not None:
                    lines.append(f"  Spin contamination: {data.orbitals.spin_contamination:.3f}")

            # All multiplicity tests
            if data.orbitals.all_multiplicities_tested:
                lines.append(f"\nAll multiplicity states tested:")
                for test in data.orbitals.all_multiplicities_tested:
                    pref_marker = " ← PREFERRED" if test['is_preferred'] else ""
                    lines.append(f"  M={test['multiplicity']}: E = {test['energy_hartree']:.6f} Eh, S² contamination = {test['spin_contamination']:.3f if test['spin_contamination'] is not None else 'N/A'}{pref_marker}")

            # Boltzmann populations
            if data.orbitals.boltzmann_populations:
                lines.append(f"\nBoltzmann populations at 298.15 K:")
                for mult, pop in sorted(data.orbitals.boltzmann_populations.items()):
                    lines.append(f"  M={mult}: {pop:.1f}%")

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

        # OCCUPIER steps - evolution of electronic structure with oxidation state
        if data.occupier_steps:
            lines.append("\n=== OCCUPIER REDOX EVOLUTION ===")
            lines.append(f"Total redox steps analyzed: {len(data.occupier_steps)}")

            for step in data.occupier_steps:
                lines.append(f"\n{step.step_name.upper()}:")
                if step.charge is not None:
                    lines.append(f"  Charge: {step.charge:+d}")
                if step.orbitals and step.orbitals.preferred_multiplicity:
                    lines.append(f"  Preferred multiplicity: {step.orbitals.preferred_multiplicity}")
                    if step.orbitals.unpaired_electrons is not None:
                        lines.append(f"  Unpaired electrons: {step.orbitals.unpaired_electrons}")
                    if step.orbitals.spin_contamination is not None:
                        lines.append(f"  Spin contamination: {step.orbitals.spin_contamination:.3f}")
                if step.energy_hartree is not None:
                    lines.append(f"  Energy: {step.energy_hartree:.6f} Eh ({step.energy_ev:.2f} eV)")

                # Show all tested multiplicities for this step
                if step.orbitals and step.orbitals.all_multiplicities_tested:
                    lines.append(f"  Multiplicities tested:")
                    for test in step.orbitals.all_multiplicities_tested:
                        marker = " ← PREFERRED" if test['is_preferred'] else ""
                        lines.append(f"    M={test['multiplicity']}: {test['energy_hartree']:.6f} Eh{marker}")

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
