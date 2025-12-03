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

        # Step 0: Ask what user wants to calculate
        print("━━━ What do you want to calculate? ━━━")
        print()
        print("DELFIN can calculate:")
        print("  • Redox potentials (oxidation/reduction)")
        print("  • Excited states (absorption/emission)")
        print("  • E_00 energies")
        print("  • Intersystem crossing (ISC) and internal conversion (IC)")
        print("  • Ground state optimizations")
        print()
        user_goal = input("Tell me what you want to calculate (or press Enter to skip): ").strip()

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
                content=f"I want to calculate: {user_goal}. Please help me set up DELFIN accordingly."
            ))
            # Let AI understand the goal
            try:
                response = self.provider.chat(
                    messages=self.conversation_history,
                    temperature=0.0
                )
                print(f"\n💡 AI: {response.data.get('text', 'Understood!')}\n")
                self.conversation_history.append(Message(
                    role="assistant",
                    content=response.data.get('text', 'Understood!')
                ))
            except Exception:
                pass  # Continue even if this fails

        # === NEW WORKFLOW ===

        # Step 1: Charge
        print("\n━━━ Step 1: Charge ━━━")
        charge_data = self._ask_step("charge", STEP_SCHEMAS["charge"],
                                      "What is the total charge of your molecule?")
        self.collected_data.update(charge_data)

        # Step 2: Solvent
        print("\n━━━ Step 2: Solvent ━━━")
        solvent_data = self._ask_step("solvent", STEP_SCHEMAS["solvent"],
                                       "What solvation model and solvent do you want to use?")
        self.collected_data.update(solvent_data)

        # Step 3: Structure quality / XTB optimization
        print("\n━━━ Step 3: Structure Quality ━━━")
        structure_data = self._ask_step("structure_quality", STEP_SCHEMAS["structure_quality"],
                                         "Is the input structure an approximation, or is it a reliable geometry? "
                                         "If approximation, I recommend XTB_OPT=yes and XTB_GOAT=yes for pre-optimization.")
        self.collected_data.update(structure_data)

        # Auto-set XTB options if structure is approximation
        if structure_data.get("is_approximation") == "yes":
            self.collected_data["xtb_opt"] = structure_data.get("XTB_OPT", "yes")
            self.collected_data["xtb_goat"] = structure_data.get("XTB_GOAT", "yes")

        # Step 4: Imaginary frequencies
        print("\n━━━ Step 4: Imaginary Frequencies ━━━")
        imag_data = self._ask_step("imag", STEP_SCHEMAS["imag"],
                                    "Do you want to check for and eliminate imaginary frequencies? "
                                    "Imaginary frequencies indicate transition states (not stable minima), "
                                    "which can cause incorrect energies. IMAG_option=2 is recommended.")
        self.collected_data.update(imag_data)
        if imag_data.get("IMAG") == "yes":
            self.collected_data["imag"] = "yes"
            self.collected_data["imag_option"] = imag_data.get("IMAG_option", "2")

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
            confirm = input("\nAccept these values? [y/n/edit]: ").lower().strip()
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
