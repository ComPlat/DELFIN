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

        # Step 1: Charge (direct input, with optional AI help)
        print("\n━━━ Step 1: Charge ━━━")
        print("What is the total charge of your molecule?")
        print("Examples: 0 (neutral), +1 (cation), +2 (dication), -1 (anion)")
        print("Type '?' for help or enter a number")
        print()

        charge = self._ask_with_ai_help(
            prompt="Charge: ",
            value_type="integer",
            help_context="molecular charge, oxidation states, and why it matters for calculations"
        )
        self.collected_data["charge"] = str(charge)
        print(f"✓ Charge set to {charge}")

        # Step 2: Solvent
        print("\n━━━ Step 2: Solvent ━━━")
        print("Solvation model: CPCM (default), SMD, or ALPB")
        print("Type '?' for help or enter model name")
        print()

        solvation_model = self._ask_with_ai_help(
            prompt="Solvation model [CPCM/SMD/ALPB]: ",
            value_type="choice",
            allowed_values=["CPCM", "SMD", "ALPB", "cpcm", "smd", "alpb"],
            help_context="solvation models (CPCM vs SMD vs ALPB), their accuracy and computational cost"
        ).upper()

        print(f"✓ Solvation model: {solvation_model}")
        print()
        print("Common solvents: Water, MeCN, DMF, DMSO, Acetone, THF, DCM, Toluene")
        print("Type '?' for help or enter solvent name")
        print()

        solvent = self._ask_with_ai_help(
            prompt="Solvent: ",
            value_type="string",
            help_context="solvent selection, polarity effects, and common solvents for your calculation"
        )

        print(f"✓ Solvent: {solvent}")

        solvent_data = {
            "implicit_solvation_model": solvation_model,
            "solvent": solvent
        }
        self.collected_data.update(solvent_data)

        # Step 3: Structure quality / XTB optimization
        print("\n━━━ Step 3: Structure Quality ━━━")
        print("Where does your structure come from?")
        print("Examples: ChemDraw, SMILES, XRD, DFT_optimized, GFN-xTB, PM6")
        print("Type '?' for help")
        print()

        structure_source = self._ask_with_ai_help(
            prompt="Structure source: ",
            value_type="string",
            help_context="structure sources (ChemDraw/SMILES/XRD/DFT), quality assessment, and when XTB pre-optimization is needed"
        )
        self.collected_data["structure_source"] = structure_source

        print()
        print("Do you want XTB pre-optimization?")
        print("Recommended: yes for ChemDraw/SMILES, no for XRD/DFT_optimized")
        print("Type '?' for help or enter yes/no")
        print()

        xtb_opt = self._ask_with_ai_help(
            prompt="XTB_OPT [yes/no]: ",
            value_type="choice",
            allowed_values=["yes", "no"],
            help_context="XTB pre-optimization, when it's needed, and XTB_GOAT optimizer"
        )
        self.collected_data["xtb_opt"] = xtb_opt
        self.collected_data["XTB_OPT"] = xtb_opt

        if xtb_opt == "yes":
            print()
            xtb_goat = self._ask_with_ai_help(
                prompt="XTB_GOAT (robust optimizer) [yes/no]: ",
                value_type="choice",
                allowed_values=["yes", "no"],
                help_context="XTB GOAT optimizer for difficult geometries"
            )
            self.collected_data["xtb_goat"] = xtb_goat
            self.collected_data["XTB_GOAT"] = xtb_goat
        else:
            self.collected_data["xtb_goat"] = "no"
            self.collected_data["XTB_GOAT"] = "no"

        print(f"✓ Structure: {structure_source}, XTB_OPT={xtb_opt}")

        # Step 4: Imaginary frequencies
        print("\n━━━ Step 4: Imaginary Frequencies ━━━")
        print("Check for imaginary frequencies (transition states)?")
        print("Imaginary frequencies = NOT a stable minimum → wrong energies")
        print("IMAG mode finds the nearest true minimum structure")
        print("Type '?' for detailed explanation")
        print()

        imag_check = self._ask_with_ai_help(
            prompt="Check for imaginary frequencies? [yes/no]: ",
            value_type="choice",
            allowed_values=["yes", "no"],
            help_context="imaginary frequencies, why they cause problems, and how IMAG fixes them"
        )

        if imag_check == "yes":
            print()
            print("IMAG elimination method:")
            print("  1 = Simple distortion (fastest)")
            print("  2 = Distortion + reoptimization (recommended)")
            print("  3 = Multiple attempts (most thorough)")
            print("Type '?' for help")
            print()

            imag_option = self._ask_with_ai_help(
                prompt="IMAG_option [1/2/3]: ",
                value_type="choice",
                allowed_values=["1", "2", "3"],
                help_context="IMAG elimination methods, their thoroughness and computational cost"
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
        print("Do you want to calculate redox potentials?")
        print("Type '?' for help")
        print()

        calc_redox = self._ask_with_ai_help(
            prompt="Calculate redox potentials? [yes/no]: ",
            value_type="choice",
            allowed_values=["yes", "no"],
            help_context="redox potentials, oxidation and reduction, and how they're calculated"
        )

        if calc_redox == "yes":
            self.collected_data["calc_initial"] = "yes"

            print()
            print("Oxidation steps (e.g., '1,2' for two steps, or leave empty)")
            print("Type '?' for help")
            print()

            ox_steps = self._ask_with_ai_help(
                prompt="Oxidation steps (comma-separated or empty): ",
                value_type="string",
                help_context="oxidation steps, how many to calculate, and notation"
            )
            self.collected_data["oxidation_steps"] = ox_steps.strip()

            print()
            print("Reduction steps (e.g., '1,2' for two steps, or leave empty)")
            print("Type '?' for help")
            print()

            red_steps = self._ask_with_ai_help(
                prompt="Reduction steps (comma-separated or empty): ",
                value_type="string",
                help_context="reduction steps, how many to calculate, and notation"
            )
            self.collected_data["reduction_steps"] = red_steps.strip()

            print(f"✓ Redox: ox={ox_steps if ox_steps else 'none'}, red={red_steps if red_steps else 'none'}")
        else:
            self.collected_data["calc_initial"] = "yes"
            self.collected_data["oxidation_steps"] = ""
            self.collected_data["reduction_steps"] = ""
            print("✓ No redox calculations")

        # Step 6: Method selection
        print("\n━━━ Step 6: Method Selection ━━━")
        print("Which method for spin state determination?")
        print("  classic = Known spin state (organic molecules)")
        print("  OCCUPIER = Automatic spin exploration (metal complexes)")
        print("  manually = Expert manual control")
        print("Type '?' for help")
        print()

        method = self._ask_with_ai_help(
            prompt="Method [classic/OCCUPIER/manually]: ",
            value_type="choice",
            allowed_values=["classic", "OCCUPIER", "manually"],
            help_context="calculation methods: classic vs OCCUPIER vs manually, when to use which"
        )
        self.collected_data["method"] = method
        print(f"✓ Method: {method}")

        # Step 7: MANUALLY configuration (if method=manually)
        if method == "manually":
            print("\n━━━ Step 7: MANUALLY Configuration ━━━")
            print("You selected 'manually' method. Let's configure the electronic structure for each redox state.")
            print("You need to specify multiplicity and additional ORCA keywords (e.g., BrokenSym flags) for each state.")
            print()

            manually_config = self._configure_manually_interactive()
            self.collected_data.update(manually_config)

        # Step 8: Excited states / E_00
        print("\n━━━ Step 8: Excited State Calculations ━━━")
        print("Calculate excited states (E_00, TDDFT)?")
        print("Note: Only for closed-shell ground states")
        print("Type '?' for help")
        print()

        calc_excited = self._ask_with_ai_help(
            prompt="Calculate excited states? [yes/no]: ",
            value_type="choice",
            allowed_values=["yes", "no"],
            help_context="excited states, E_00 energies, TDDFT, and when they're applicable"
        )

        if calc_excited == "yes":
            self.collected_data["calculate_excited_states"] = "yes"
            self.collected_data["e_00"] = "yes"
            self.collected_data["E_00"] = "yes"

            print()
            print("Excitation type: s (singlet), t (triplet), or s,t (both)")
            print("Type '?' for help")
            print()

            excitation = self._ask_with_ai_help(
                prompt="Excitation [s/t/s,t]: ",
                value_type="string",
                help_context="singlet vs triplet excitations, and when to calculate both"
            )
            self.collected_data["excitation"] = excitation

            print()
            print("Calculate ESD (ISC/IC rates)?")
            print()

            esd = self._ask_with_ai_help(
                prompt="ESD_modul [yes/no]: ",
                value_type="choice",
                allowed_values=["yes", "no"],
                help_context="ESD module for ISC and IC rate calculations"
            )
            self.collected_data["ESD_modul"] = esd

            print(f"✓ Excited states: excitation={excitation}, ESD={esd}")
        else:
            self.collected_data["calculate_excited_states"] = "no"
            self.collected_data["e_00"] = "no"
            self.collected_data["E_00"] = "no"
            print("✓ No excited state calculations")

        # Step 9: Level of theory
        print("\n━━━ Step 9: Level of Theory ━━━")
        print("DFT functional (e.g., PBE0, wB97X-V, B3LYP)")
        print("Type '?' for help")
        print()

        functional = self._ask_with_ai_help(
            prompt="Functional: ",
            value_type="string",
            help_context="DFT functionals, which to use for different systems (metal complexes, organics, excited states)"
        )
        self.collected_data["functional"] = functional

        print()
        print("Basis set (e.g., def2-TZVP, def2-SVP)")
        print()

        basisset = self._ask_with_ai_help(
            prompt="Basis set: ",
            value_type="string",
            help_context="basis sets, size vs accuracy tradeoff"
        )
        self.collected_data["main_basisset"] = basisset
        self.collected_data["metal_basisset"] = basisset

        print(f"✓ Theory: {functional}/{basisset}")

        # Step 10: Computational resources
        print("\n━━━ Step 10: Computational Resources ━━━")
        print("How many CPU cores to use?")
        print()

        pal = self._ask_with_ai_help(
            prompt="PAL (CPU cores): ",
            value_type="integer",
            help_context="parallel calculations, number of CPU cores"
        )
        self.collected_data["PAL"] = pal

        print()
        print("Memory per core in MB (total memory = PAL × maxcore)")
        print()

        maxcore = self._ask_with_ai_help(
            prompt="maxcore (MB per core): ",
            value_type="integer",
            help_context="memory allocation per CPU core"
        )
        self.collected_data["maxcore"] = maxcore

        print(f"✓ Resources: {pal} cores, {maxcore} MB/core (total: {pal*maxcore} MB)")

        # Step 11: Auto recovery (optional - skip for now, use defaults)
        print("\n━━━ Skipping optional advanced settings ━━━")
        print("(Auto-recovery, parallel workflows, etc. - using defaults)")
        self.collected_data["auto_recovery"] = "yes"
        self.collected_data["enable_job_timeouts"] = "yes"

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

    def _ask_with_ai_help(
        self,
        prompt: str,
        value_type: str = "string",
        help_context: str = "",
        allowed_values: List[str] = None
    ):
        """
        Ask user for input with optional AI help.

        User can:
        - Enter a value directly
        - Type '?' or a question to get AI explanation
        - After AI response, enter the value

        Args:
            prompt: Input prompt
            value_type: "integer", "string", "choice"
            help_context: Context for AI to explain
            allowed_values: For choice type, list of allowed values

        Returns:
            User's input (converted to appropriate type)
        """
        while True:
            user_input = input(prompt).strip()

            # Check if user wants help
            if (user_input == '?' or user_input.lower() == 'help' or
                'warum' in user_input.lower() or 'why' in user_input.lower() or
                'what' in user_input.lower() or 'wie' in user_input.lower() or
                'how' in user_input.lower() or user_input.endswith('?')):

                # Get AI explanation
                question = user_input if user_input not in ['?', 'help'] else f"Explain {help_context}"
                self.conversation_history.append(Message(
                    role="user",
                    content=f"User question about {help_context}: {question}. "
                            f"Provide a SHORT explanation (max 5 lines)."
                ))

                try:
                    response = self.provider.chat(
                        messages=self.conversation_history,
                        temperature=0.0
                    )
                    print(f"\n💡 AI: {response.data.get('text', 'No explanation available')}\n")
                    self.conversation_history.append(Message(
                        role="assistant",
                        content=response.data.get('text', '')
                    ))
                except Exception as e:
                    print(f"AI help not available: {e}\n")

                # Ask again for the actual value
                continue

            # Validate and return the input
            if value_type == "integer":
                try:
                    return int(user_input)
                except ValueError:
                    print(f"Error: '{user_input}' is not a valid integer. Please enter a number.\n")
                    continue

            elif value_type == "choice" and allowed_values:
                if user_input in allowed_values:
                    return user_input
                else:
                    print(f"Error: Please choose from: {', '.join(allowed_values)}\n")
                    continue

            else:  # string
                if user_input:
                    return user_input
                else:
                    print("Error: Please enter a value.\n")
                    continue

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
