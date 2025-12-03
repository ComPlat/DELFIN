"""
DELFIN Control Assistant (Complete Rewrite)

Interactive LLM-powered assistant for creating CONTROL files with comprehensive workflow.
"""

from __future__ import annotations
import os
from pathlib import Path
from typing import Dict, Any, Optional, List

from .base_provider import BaseLLMProvider, Message
from .control_schemas import STEP_SCHEMAS
from .control_template import CONTROL_TEMPLATE, DEFAULT_VALUES
from .delfin_knowledge import DELFIN_KNOWLEDGE, SYSTEM_TYPE_RECOMMENDATIONS
from .control_examples import format_example_for_prompt


class ControlAssistantV2:
    """Interactive assistant for creating DELFIN CONTROL files (V2)"""

    def __init__(self, provider: BaseLLMProvider):
        """
        Initialize Control Assistant.

        Args:
            provider: LLM provider (Claude, Ollama, etc.)
        """
        self.provider = provider
        self.conversation_history: List[Message] = []
        self.collected_data: Dict[str, Any] = DEFAULT_VALUES.copy()

    def create_control_interactive(self, output_path: str = "CONTROL.txt") -> Dict[str, Any]:
        """
        Create CONTROL file through interactive conversation with new workflow.

        Args:
            output_path: Path to save CONTROL file

        Returns:
            Dictionary with collected parameters
        """
        print("╔══════════════════════════════════════════════════════════╗")
        print("║     DELFIN Control Assistant (AI-powered)                ║")
        print("║                                                          ║")
        print(f"║  Provider: {self.provider.provider_name:20s}                      ║")
        print(f"║  Model:    {self.provider.model:30s}          ║")
        print("╚══════════════════════════════════════════════════════════╝")
        print()

        # Step 0: Optional AI consultation
        print("━━━ AI Consultation (Optional) ━━━")
        print()
        print("DELFIN can calculate:")
        print("  • Redox potentials (oxidation/reduction)")
        print("  • Excited states (absorption/emission)")
        print("  • E_00 energies")
        print("  • Intersystem crossing (ISC) and internal conversion (IC)")
        print("  • Ground state optimizations")
        print()
        print("Would you like recommendations from the AI assistant?")
        print("Example: 'I need redox potentials for an iron complex'")
        print("         'Carbazole excited states and E_00'")
        print()
        user_goal = input("Describe your calculation goal (or press Enter to skip): ").strip()

        # Create system prompt with DELFIN knowledge
        system_prompt = self._create_system_prompt()
        if user_goal:
            # Add relevant CONTROL example based on user goal
            example_prompt = format_example_for_prompt(user_goal)
            system_prompt += "\n\n" + example_prompt

        self.conversation_history.append(Message(role="system", content=system_prompt))

        if user_goal:
            self.conversation_history.append(Message(
                role="user",
                content=f"""User goal: {user_goal}

Please provide:
1. BRIEF recommendation of key parameters (charge, solvent, method, functional)
2. DETAILED explanation of special requirements (e.g., if they mention ISC/IC, E_00, spin states)
3. What DELFIN modules to enable (ESD_modul, E_00, OCCUPIER, etc.) and WHY
4. Any important caveats or requirements

Format your answer in 2 sections:
**Quick Setup:**
- charge: ...
- method: ...
- functional: ...

**Detailed Explanation:**
[Explain the special requirements and how DELFIN handles them]

Keep total length under 15 lines."""
            ))
            # Get AI recommendations
            try:
                response = self.provider.chat(
                    messages=self.conversation_history,
                    temperature=0.0
                )
                print(f"\n💡 AI Recommendations:\n{response.data.get('text', 'No recommendations')}\n")
                print("━" * 60)
                print("Now let's set up your CONTROL file step by step.")
                print("━" * 60)
                self.conversation_history.append(Message(
                    role="assistant",
                    content=response.data.get('text', '')
                ))
            except Exception as e:
                print(f"Note: AI consultation failed ({e}), continuing with manual setup...")
                pass

        # === NEW WORKFLOW ===

        # Step 1: Charge (direct input, no AI suggestion)
        print("\n━━━ Step 1: Charge ━━━")
        print("What is the total charge of your molecule?")
        print("Examples: 0 (neutral), +1 (cation), +2 (dication), -1 (anion)")
        while True:
            charge_input = input("Charge: ").strip()
            try:
                charge = int(charge_input)
                self.collected_data["charge"] = str(charge)
                print(f"✓ Charge set to {charge}")
                break
            except ValueError:
                print(f"Error: '{charge_input}' is not a valid integer. Please enter a number (e.g., 0, +1, -1).")

        # Step 2: Solvent
        print("\n━━━ Step 2: Solvent ━━━")

        # Ask for solvation model
        print("Select solvation model:")
        print("  1) CPCM - Conductor-like Polarizable Continuum Model (general-purpose, reliable)")
        print("  2) SMD - Solvation Model based on Density (more accurate, slower)")
        print("  3) ALPB - Analytical Linearized Poisson-Boltzmann (fast, compatible with xTB)")
        print()

        solvation_choice = self._ask_user_direct_choice(
            prompt="Select solvation model [1-3]: ",
            choices={
                "1": "CPCM",
                "2": "SMD",
                "3": "ALPB",
            },
            explanations={
                "1": "CPCM: Best general-purpose choice, well-tested, reliable",
                "2": "SMD: Higher accuracy for polar solvents, but computationally more expensive",
                "3": "ALPB: Fast, works well with xTB pre-optimization",
            }
        )

        # Ask for solvent
        print("\nSelect solvent:")
        print("Common solvents:")
        print("  1) Water         2) MeCN (Acetonitrile)    3) DMF")
        print("  4) DMSO          5) Acetone                6) THF")
        print("  7) DCM           8) Chloroform             9) Toluene")
        print(" 10) Methanol     11) Ethanol               12) other")
        print()

        solvent_map = {
            "1": "Water", "2": "MeCN", "3": "DMF", "4": "DMSO",
            "5": "Acetone", "6": "THF", "7": "DCM", "8": "Chloroform",
            "9": "Toluene", "10": "Methanol", "11": "Ethanol",
        }

        solvent_input = input("Select solvent [1-12]: ").strip()

        if solvent_input == "12":
            # User wants to enter custom solvent
            print("\nAvailable solvents: Water, MeCN, DMF, DMSO, Acetone, THF, DCM,")
            print("Chloroform, Toluene, Benzene, Hexane, Methanol, Ethanol, DiethylEther")
            solvent = input("Enter solvent name: ").strip()
        elif solvent_input in solvent_map:
            solvent = solvent_map[solvent_input]
        else:
            print(f"Invalid choice, defaulting to Water")
            solvent = "Water"

        print(f"✓ Solvation: {solvation_choice} with {solvent}")

        solvent_data = {
            "implicit_solvation_model": solvation_choice,
            "solvent": solvent
        }
        self.collected_data.update(solvent_data)

        # Step 3: Structure quality / XTB optimization
        print("\n━━━ Step 3: Structure Quality ━━━")
        print("Where does the input structure come from?")
        print("  1) ChemDraw - 2D drawing converted to 3D")
        print("  2) SMILES - String notation converted to 3D")
        print("  3) GFN-xTB / GFN2-xTB - Semi-empirical optimization")
        print("  4) PM6 / PM7 - Semi-empirical optimization")
        print("  5) XRD - X-ray diffraction structure")
        print("  6) DFT_optimized - Already optimized with DFT")
        print("  7) other_semi_empirical - Other semi-empirical method")
        print()

        structure_data = self._ask_user_direct_choice(
            prompt="Select structure source [1-7]: ",
            choices={
                "1": {"structure_source": "ChemDraw", "XTB_OPT": "yes", "XTB_GOAT": "yes"},
                "2": {"structure_source": "SMILES", "XTB_OPT": "yes", "XTB_GOAT": "yes"},
                "3": {"structure_source": "GFN-xTB", "XTB_OPT": "yes", "XTB_GOAT": "yes"},
                "4": {"structure_source": "PM6", "XTB_OPT": "yes", "XTB_GOAT": "yes"},
                "5": {"structure_source": "XRD", "XTB_OPT": "no", "XTB_GOAT": "no"},
                "6": {"structure_source": "DFT_optimized", "XTB_OPT": "no", "XTB_GOAT": "no"},
                "7": {"structure_source": "other_semi_empirical", "XTB_OPT": "yes", "XTB_GOAT": "yes"},
            },
            explanations={
                "1": "ChemDraw structures need XTB pre-optimization (unrealistic bond lengths)",
                "2": "SMILES structures need XTB pre-optimization (converted from string)",
                "3": "GFN-xTB optimized but may benefit from GOAT for robustness",
                "4": "PM6/PM7 optimized but may benefit from XTB refinement",
                "5": "XRD structures are high-quality experimental data (no pre-opt needed)",
                "6": "DFT optimized structures are already high-quality (no pre-opt needed)",
                "7": "Other semi-empirical methods may need XTB refinement",
            }
        )
        self.collected_data.update(structure_data)

        # Store XTB options
        self.collected_data["xtb_opt"] = structure_data.get("XTB_OPT", "no")
        self.collected_data["xtb_goat"] = structure_data.get("XTB_GOAT", "no")

        # Step 4: Imaginary frequencies
        print("\n━━━ Step 4: Imaginary Frequencies ━━━")
        print("Do you want to check for and eliminate imaginary frequencies?")
        print()
        print("What are imaginary frequencies?")
        print("  Imaginary frequencies indicate that your structure is NOT a stable minimum,")
        print("  but rather a transition state (saddle point on the potential energy surface).")
        print("  This leads to INCORRECT energies and can cause calculation failures.")
        print()
        print("Why eliminate them?")
        print("  IMAG mode distorts the structure along imaginary modes and reoptimizes")
        print("  to find the nearest TRUE minimum energy structure.")
        print()
        print("Options:")
        print("  1) yes - Check and eliminate imaginary frequencies (RECOMMENDED)")
        print("  2) no  - Skip this step (only if you're sure structure is already optimized)")
        print()

        imag_choice = self._ask_user_direct_choice(
            prompt="Check for imaginary frequencies? [1-2]: ",
            choices={
                "1": "yes",
                "2": "no",
            },
            explanations={
                "1": "RECOMMENDED: Ensures your structure is a true minimum, prevents wrong energies",
                "2": "Skip if structure is already well-optimized (e.g., from previous DFT calc)",
            }
        )

        if imag_choice == "yes":
            print("\nSelect IMAG elimination method:")
            print("  1) IMAG_option=1 - Simple mode distortion")
            print("  2) IMAG_option=2 - Distortion + reoptimization (RECOMMENDED)")
            print("  3) IMAG_option=3 - Multiple distortions and reoptimizations")
            print()

            imag_option = self._ask_user_direct_choice(
                prompt="Select IMAG_option [1-3]: ",
                choices={
                    "1": "1",
                    "2": "2",
                    "3": "3",
                },
                explanations={
                    "1": "Fastest, but may not eliminate all imaginary frequencies",
                    "2": "RECOMMENDED: Good balance between speed and thoroughness",
                    "3": "Most thorough, but slowest (use for difficult cases)",
                }
            )

            self.collected_data["imag"] = "yes"
            self.collected_data["IMAG"] = "yes"
            self.collected_data["imag_option"] = imag_option
            self.collected_data["IMAG_option"] = imag_option
            print(f"✓ IMAG check enabled with option {imag_option}")
        else:
            self.collected_data["imag"] = "no"
            self.collected_data["IMAG"] = "no"
            print("✓ IMAG check disabled")

        # Step 5: Redox potentials setup
        print("\n━━━ Step 5: Redox Potentials Setup ━━━")
        redox_setup = self._ask_step("redox_setup", STEP_SCHEMAS["redox_setup"],
                                      "Which redox states do you want to calculate? "
                                      "For example: oxidation_steps='1,2,3' means three oxidation steps.")
        self.collected_data.update(redox_setup)

        # Step 6: Method selection
        print("\n━━━ Step 6: Method Selection ━━━")
        method_data = self._ask_step("method", STEP_SCHEMAS["method"],
                                      "What type of system are you calculating? "
                                      "For transition metal complexes → OCCUPIER (spin states are non-trivial). "
                                      "For TADF organic molecules → classic. "
                                      "For expert control → manually (define spin states and broken-symmetry yourself).")
        self.collected_data.update(method_data)

        # Step 7: MANUALLY configuration (if method=manually)
        if method_data.get("method") == "manually":
            print("\n━━━ Step 7: MANUALLY Configuration ━━━")
            print("You selected 'manually' method. Let's configure the electronic structure for each redox state.")
            print("You need to specify multiplicity and additional ORCA keywords (e.g., BrokenSym flags) for each state.")
            print()

            manually_config = self._configure_manually_interactive()
            self.collected_data.update(manually_config)

        # Step 8: Excited states / E_00
        print("\n━━━ Step 8: Excited State Calculations ━━━")
        excited_data = self._ask_step("excited_states", STEP_SCHEMAS["excited_states"],
                                       "Do you want to calculate E_00 energies or excited state redox potentials? "
                                       "NOTE: This is only available for closed-shell systems. "
                                       "You can calculate singlet E_00, triplet E_00, or both.")
        self.collected_data.update(excited_data)

        if excited_data.get("calculate_excited_states") == "yes":
            self.collected_data["e_00"] = excited_data.get("E_00", "yes")
            self.collected_data["excitation"] = excited_data.get("excitation", "s")
            self.collected_data["absorption_spec"] = excited_data.get("absorption_spec", "no")
            self.collected_data["emission_spec"] = excited_data.get("emission_spec", "no")

        # Step 9: Level of theory
        print("\n━━━ Step 9: Level of Theory ━━━")
        theory_data = self._ask_step("level_of_theory", STEP_SCHEMAS["level_of_theory"],
                                      f"What level of theory do you want? "
                                      f"System type: {method_data.get('system_type', 'unknown')}. "
                                      f"I can recommend functional, basis sets, dispersion correction, and relativistic treatment.")
        self.collected_data.update(theory_data)

        # Step 10: Additional analysis
        print("\n━━━ Step 10: Additional Analysis ━━━")
        analysis_data = self._ask_step("analysis", STEP_SCHEMAS["analysis"],
                                        "Do you want to plot MOs or perform Loewdin population analysis? "
                                        "These help understand electron distribution and bonding.")
        self.collected_data.update(analysis_data)

        # Step 11: Computational resources
        print("\n━━━ Step 11: Computational Resources ━━━")
        resources_data = self._ask_step("resources", STEP_SCHEMAS["resources"],
                                         "How many CPU cores (PAL) and how much memory per core (maxcore in MB) are available? "
                                         "Total memory = PAL × maxcore. "
                                         "You can also enable parallel workflows to run multiple jobs simultaneously.")
        self.collected_data.update(resources_data)

        # Step 12: Auto recovery
        print("\n━━━ Step 12: Automatic Error Recovery ━━━")
        recovery_data = self._ask_step("auto_recovery", STEP_SCHEMAS["auto_recovery"],
                                        "Do you want to enable automatic error recovery? "
                                        "This automatically retries failed calculations with adjusted settings. "
                                        "Job timeouts prevent stuck calculations from running forever.")
        self.collected_data.update(recovery_data)

        # Generate CONTROL file
        print("\n━━━ Generating CONTROL file... ━━━")
        self._write_control_file(output_path)

        print(f"\n✓ CONTROL file created: {output_path}")

        # AI Review of CONTROL file
        print("\n━━━ AI Review of CONTROL File ━━━")
        print("Let the AI check your CONTROL file for errors and optimization suggestions...")
        try:
            self._ai_review_control_file(output_path)
        except Exception as e:
            print(f"Note: AI review failed ({e}), but CONTROL file was created successfully.")

        print(f"\n✓ You can now run: delfin")

        return self.collected_data

    def _create_system_prompt(self) -> str:
        """Create comprehensive system prompt with DELFIN knowledge"""
        return f"""You are DELFIN Control Assistant, an expert in computational chemistry and the DELFIN software.

{DELFIN_KNOWLEDGE}

CRITICAL RULES:
1. LISTEN to what the user wants to calculate FIRST
2. You MUST use ONLY values from the provided schema
3. You MUST NOT invent or guess values
4. Provide intelligent recommendations based on system type
5. Explain options clearly and concisely
6. For transition metal complexes → recommend OCCUPIER (spin states non-trivial)
7. For TADF organic molecules → recommend classic method
8. For E_00 calculations → only available for closed-shell systems
9. Always explain WHY you recommend certain settings
10. Use structured output with the exact schema provided

Remember: DELFIN users are computational chemists. Be professional, accurate, and helpful.
"""

    def _ask_step(self, step_name: str, schema: Dict[str, Any], prompt: str) -> Dict[str, Any]:
        """
        Ask user for parameters in this step.

        Args:
            step_name: Name of the step
            schema: JSON schema for this step
            prompt: Human-readable prompt for this step

        Returns:
            Dictionary with collected parameters
        """
        # Create message for this step
        self.conversation_history.append(Message(role="user", content=prompt))

        # Get LLM response
        response = self.provider.chat(
            messages=self.conversation_history,
            schema=schema,
            temperature=0.0
        )

        # Display to user and confirm
        print("\nAI suggests:")
        for key, value in response.data.items():
            print(f"  {key:35s} = {value}")

        while True:
            confirm = input("\nAccept these values? [y/n/edit/?]: ").lower().strip()
            if confirm == 'y' or confirm == '':
                self.conversation_history.append(Message(
                    role="assistant",
                    content=f"Suggested: {response.data}"
                ))
                return response.data
            elif confirm == 'n':
                # Ask again
                user_feedback = input("What would you like to change? ")
                self.conversation_history.append(Message(
                    role="assistant",
                    content=f"Suggested: {response.data}"
                ))
                self.conversation_history.append(Message(
                    role="user",
                    content=f"Please adjust: {user_feedback}"
                ))
                return self._ask_step(step_name, schema, prompt)
            elif confirm == 'edit':
                # Manual editing
                return self._manual_edit(response.data, schema)
            elif (confirm == '?' or confirm == 'help' or
                  confirm.startswith('what') or confirm.startswith('explain') or
                  confirm.startswith('warum') or confirm.startswith('why') or
                  'why' in confirm.lower() or 'warum' in confirm.lower() or 'what' in confirm.lower()):
                # User wants explanation
                self.conversation_history.append(Message(
                    role="assistant",
                    content=f"Suggested: {response.data}"
                ))
                self.conversation_history.append(Message(
                    role="user",
                    content=f"Please explain: {confirm}"
                ))
                # Get explanation from AI
                explanation_response = self.provider.chat(
                    messages=self.conversation_history,
                    temperature=0.0
                )
                print(f"\n💡 AI: {explanation_response.data.get('text', 'No explanation available')}\n")
                self.conversation_history.append(Message(
                    role="assistant",
                    content=explanation_response.data.get('text', 'No explanation available')
                ))
                # Continue the loop - ask again if they accept
                continue
            else:
                # If user types anything else, treat it as feedback for changes
                self.conversation_history.append(Message(
                    role="assistant",
                    content=f"Suggested: {response.data}"
                ))
                self.conversation_history.append(Message(
                    role="user",
                    content=f"Please adjust: {confirm}"
                ))
                return self._ask_step(step_name, schema, prompt)

    def _ask_user_direct_choice(
        self,
        prompt: str,
        choices: Dict[str, Any],  # Can be Dict[str, str] or Dict[str, Dict]
        explanations: Dict[str, str] = None
    ):
        """
        Ask user to directly select from numbered choices (no AI inference).

        Args:
            prompt: Input prompt to display
            choices: Dict mapping choice number to result data (can be string or dict)
            explanations: Optional dict mapping choice number to explanation

        Returns:
            Selected choice data (string or dict)
        """
        while True:
            user_input = input(prompt).strip()

            # Check for help/explanation request
            if user_input == '?' or user_input.lower() == 'help':
                if explanations:
                    print("\nExplanations:")
                    for num, explanation in explanations.items():
                        print(f"  {num}) {explanation}")
                    print()
                else:
                    print("No detailed explanations available.\n")
                continue

            # Validate choice
            if user_input in choices:
                result = choices[user_input]
                if isinstance(result, dict):
                    print(f"✓ Selected: {result}")
                else:
                    print(f"✓ Selected: {result}")
                return result
            else:
                print(f"Invalid choice '{user_input}'. Please select from: {', '.join(choices.keys())}")
                print("Type '?' for explanations.\n")

    def _manual_edit(self, data: Dict[str, Any], schema: Dict[str, Any]) -> Dict[str, Any]:
        """Allow manual editing of parameters"""
        edited = data.copy()

        for key in schema["properties"]:
            current = edited.get(key, "")
            new_value = input(f"{key} [{current}]: ").strip()
            if new_value:
                # Try to parse to correct type
                prop_type = schema["properties"][key].get("type")
                if prop_type == "integer":
                    try:
                        edited[key] = int(new_value)
                    except ValueError:
                        print(f"Warning: {new_value} is not an integer, keeping {current}")
                else:
                    edited[key] = new_value

        return edited

    def _configure_manually_interactive(self) -> Dict[str, Any]:
        """
        Interactive configuration for MANUALLY method.
        Asks user for multiplicity and additions for each redox state.
        """
        manually_data = {}

        print("\nConfiguring MANUALLY settings...")
        print("For each state, provide:")
        print("  - multiplicity: Spin multiplicity (1=singlet, 2=doublet, 3=triplet, etc.)")
        print("  - additions: Additional ORCA keywords (e.g., '%scf BrokenSym 1,1 end')")
        print()

        # Initial state
        print("Initial state (charge={})".format(self.collected_data.get("charge", 0)))
        manually_data["multiplicity_0"] = input("  multiplicity_0: ").strip()
        manually_data["additions_0"] = input("  additions_0: ").strip()

        # TDDFT, T1, S1
        print("\nExcited state keywords:")
        manually_data["additions_TDDFT"] = input("  additions_TDDFT: ").strip()
        manually_data["additions_T1"] = input("  additions_T1: ").strip()
        manually_data["additions_S1"] = input("  additions_S1: ").strip()

        # Oxidation states
        ox_steps = self.collected_data.get("oxidation_steps", "")
        if ox_steps:
            ox_list = [x.strip() for x in ox_steps.split(",") if x.strip()]
            for i in ox_list:
                print(f"\nOxidation +{i} (charge={int(self.collected_data.get('charge', 0)) + int(i)})")
                manually_data[f"multiplicity_ox{i}"] = input(f"  multiplicity_ox{i}: ").strip()
                manually_data[f"additions_ox{i}"] = input(f"  additions_ox{i}: ").strip()

        # Reduction states
        red_steps = self.collected_data.get("reduction_steps", "")
        if red_steps:
            red_list = [x.strip() for x in red_steps.split(",") if x.strip()]
            for i in red_list:
                print(f"\nReduction -{i} (charge={int(self.collected_data.get('charge', 0)) - int(i)})")
                manually_data[f"multiplicity_red{i}"] = input(f"  multiplicity_red{i}: ").strip()
                manually_data[f"additions_red{i}"] = input(f"  additions_red{i}: ").strip()

        return manually_data

    def _ai_review_control_file(self, control_file_path: str):
        """
        AI reviews the generated CONTROL file for errors and suggests optimizations.
        """
        # Read the generated CONTROL file
        control_content = Path(control_file_path).read_text()

        # Ask AI to review
        review_prompt = f"""Review this DELFIN CONTROL file for errors and optimization opportunities.

CONTROL file:
```
{control_content}
```

Please check for:
1. **Critical Errors**: Missing parameters, incompatible settings, wrong syntax
2. **Warnings**: Suboptimal choices that may cause problems
3. **Optimization Suggestions**: Better settings for this calculation type

Provide a CONCISE review (max 10 lines). Format:
✓ OK: [things that look good]
⚠ Warning: [potential issues]
💡 Suggestion: [optimizations]
"""

        self.conversation_history.append(Message(
            role="user",
            content=review_prompt
        ))

        response = self.provider.chat(
            messages=self.conversation_history,
            temperature=0.0
        )

        review_text = response.data.get('text', 'No review available')
        print(f"\n{review_text}\n")

        self.conversation_history.append(Message(
            role="assistant",
            content=review_text
        ))

    def _write_control_file(self, output_path: str):
        """Write collected data to CONTROL file using template"""

        # Build MANUALLY section if method=manually
        if self.collected_data.get("method") == "manually":
            manually_lines = []
            for key in ["multiplicity_0", "additions_0", "additions_TDDFT", "additions_T1", "additions_S1",
                        "multiplicity_ox1", "additions_ox1", "multiplicity_ox2", "additions_ox2",
                        "multiplicity_ox3", "additions_ox3",
                        "multiplicity_red1", "additions_red1", "multiplicity_red2", "additions_red2",
                        "multiplicity_red3", "additions_red3"]:
                value = self.collected_data.get(key, "")
                manually_lines.append(f"{key}={value}")
            self.collected_data["manually_section"] = "\n".join(manually_lines)
        else:
            self.collected_data["manually_section"] = ""

        # Fill template
        control_content = CONTROL_TEMPLATE.format(**self.collected_data)

        # Write to file
        Path(output_path).write_text(control_content)
