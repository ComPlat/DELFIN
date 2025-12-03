"""
DELFIN Control Assistant

Interactive LLM-powered assistant for creating CONTROL files.
Supports multiple LLM providers (Claude, GPT, Ollama).
"""

from __future__ import annotations
import os
from pathlib import Path
from typing import Dict, Any, Optional, List

from .base_provider import BaseLLMProvider, Message
from .control_schemas import CONTROL_SCHEMA, STEP_SCHEMAS
from .control_examples import format_example_for_prompt
from ..define import create_control_file


class ControlAssistant:
    """Interactive assistant for creating DELFIN CONTROL files"""

    def __init__(self, provider: BaseLLMProvider):
        """
        Initialize Control Assistant.

        Args:
            provider: LLM provider (Claude, Ollama, etc.)
        """
        self.provider = provider
        self.conversation_history: List[Message] = []
        self.collected_data: Dict[str, Any] = {}

    def create_control_interactive(self, output_path: str = "CONTROL.txt") -> Dict[str, Any]:
        """
        Create CONTROL file through interactive conversation.

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

        # Create system prompt with relevant example
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
        else:
            # No specific goal - add generic system prompt
            self.conversation_history.append(Message(role="system", content=system_prompt))

        # Step 1: Basic information
        print("━━━ Step 1: Basic Information ━━━")
        basic_info = self._ask_step("basic_info", STEP_SCHEMAS["basic_info"])
        self.collected_data.update(basic_info)

        # Step 2: Level of theory
        print("\n━━━ Step 2: Level of Theory ━━━")
        theory = self._ask_step("level_of_theory", STEP_SCHEMAS["level_of_theory"])
        self.collected_data.update(theory)

        # Step 3: Solvation
        print("\n━━━ Step 3: Solvation ━━━")
        solvation = self._ask_step("solvation", STEP_SCHEMAS["solvation"])
        self.collected_data.update(solvation)

        # Step 4: Redox method
        print("\n━━━ Step 4: Redox Calculation ━━━")
        redox = self._ask_step("redox", STEP_SCHEMAS["redox"])
        self.collected_data.update(redox)

        # Step 5: Resources
        print("\n━━━ Step 5: Computational Resources ━━━")
        resources = self._ask_step("resources", STEP_SCHEMAS["resources"])
        self.collected_data.update(resources)

        # Generate CONTROL file
        print("\n━━━ Generating CONTROL file... ━━━")
        self._write_control_file(output_path)

        print(f"\n✓ CONTROL file created: {output_path}")
        print(f"\n✓ You can now run: delfin")

        return self.collected_data

    def _create_system_prompt(self) -> str:
        """Create system prompt with strict rules"""
        return """You are DELFIN Control Assistant, an expert in computational chemistry and the DELFIN software.

DELFIN is a workflow manager for ORCA, xTB, and CREST that can:
- Calculate redox potentials (E_ox, E_red) using OCCUPIER method
- Optimize geometries and find preferred electron configurations
- Calculate excited states (S0, S1, T1, T2) and E_00 energies
- Compute absorption and emission spectra
- Calculate intersystem crossing (ISC) and internal conversion (IC) rates
- Handle spin states automatically with the OCCUPIER method

CRITICAL RULES:
1. LISTEN to what the user wants to calculate FIRST
2. You MUST use ONLY values from the provided schema
3. You MUST NOT invent or guess values
4. If the user mentions ISC, IC, or excited states, you should enable:
   - ESD_modul=yes
   - absorption_spec and/or emission_spec
   - E_00=yes for E_00 energies
5. For redox calculations, recommend OCCUPIER method (automatic spin states)
6. Explain options clearly but concisely
7. Use structured output with the exact schema provided

When the user describes their calculation goal:
- Understand what DELFIN features they need
- Adjust your suggestions accordingly
- Explain why you're recommending certain settings

Remember: DELFIN users are computational chemists. Be professional, accurate, and helpful."""

    def _ask_step(self, step_name: str, schema: Dict[str, Any]) -> Dict[str, Any]:
        """
        Ask user for parameters in this step.

        Args:
            step_name: Name of the step
            schema: JSON schema for this step

        Returns:
            Dictionary with collected parameters
        """
        # Create prompt for this step
        user_message = self._create_step_prompt(step_name, schema)
        self.conversation_history.append(Message(role="user", content=user_message))

        # Get LLM response
        response = self.provider.chat(
            messages=self.conversation_history,
            schema=schema,
            temperature=0.0
        )

        # Display to user and confirm
        print("\nAI suggests:")
        for key, value in response.data.items():
            print(f"  {key:30s} = {value}")

        while True:
            confirm = input("\nAccept these values? [y/n/edit]: ").lower()
            if confirm == 'y':
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
                return self._ask_step(step_name, schema)
            elif confirm == 'edit':
                # Manual editing
                return self._manual_edit(response.data, schema)

    def _create_step_prompt(self, step_name: str, schema: Dict[str, Any]) -> str:
        """Create prompt for asking about this step"""
        prompts = {
            "basic_info": """Let's start with basic molecular properties.

I need to know:
1. The total charge of your system
2. The spin multiplicity (2S+1, where S is total spin)

For example:
- Neutral closed-shell molecule: charge=0, multiplicity=1
- Cation with unpaired electron: charge=+1, multiplicity=2
- Triplet state: multiplicity=3

What are the charge and multiplicity of your system?""",

            "level_of_theory": """Now let's choose the level of theory.

I need:
1. DFT functional (e.g., PBE0, wB97X-V, B3LYP)
2. Basis set for light atoms (e.g., def2-SVP, def2-TZVP)
3. Basis set for metals (e.g., def2-TZVP, SARC-ZORA-TZVP)

Common choices:
- Fast screening: PBE0 / def2-SVP / def2-TZVP
- Accurate: wB97X-V / def2-TZVP / def2-TZVP
- With heavy elements: include ZORA or DKH

What level of theory do you want?""",

            "solvation": """Let's set up solvation.

I need:
1. Solvation model (CPCM, SMD, or ALPB)
2. Solvent (Water, MeCN, DMF, etc.)

Note:
- CPCM: Good general-purpose model
- SMD: More accurate but slower
- ALPB: Fast, works with xTB

What solvation do you want?""",

            "redox": """Now for the redox calculation setup.

I need:
1. Method: classic, manually, or OCCUPIER
   - classic: Simple SCF calculations
   - manually: You specify spin states
   - OCCUPIER: Automatic electron configuration search
2. Which states to calculate (initial, oxidation, reduction)
3. How many oxidation/reduction steps

What redox setup do you want?""",

            "resources": """Finally, computational resources.

I need:
1. PAL: Number of CPU cores for ORCA
2. maxcore: Memory per core in MB

Example:
- Workstation: PAL=8, maxcore=4000 (32 GB total)
- Cluster node: PAL=32, maxcore=3800 (121 GB total)

What resources are available?""",
        }

        return prompts.get(step_name, f"Please provide values for: {list(schema['properties'].keys())}")

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

    def _write_control_file(self, output_path: str):
        """Write collected data to CONTROL file"""
        # Build CONTROL content
        lines = []

        # Basic info
        lines.append(f"input_file=input.txt")
        lines.append(f"charge={self.collected_data.get('charge', 0)}")

        # Solvation
        lines.append("------------------------------------")
        lines.append("Solvation:")
        lines.append(f"implicit_solvation_model={self.collected_data.get('implicit_solvation_model', 'CPCM')}")
        lines.append(f"solvent={self.collected_data.get('solvent', 'Water')}")

        # Level of theory
        lines.append("------------------------------------")
        lines.append("Level of Theory:")
        lines.append(f"functional={self.collected_data.get('functional', 'PBE0')}")
        lines.append(f"main_basisset={self.collected_data.get('main_basisset', 'def2-SVP')}")
        lines.append(f"metal_basisset={self.collected_data.get('metal_basisset', 'def2-TZVP')}")

        if self.collected_data.get('disp_corr'):
            lines.append(f"disp_corr={self.collected_data['disp_corr']}")
        if self.collected_data.get('ri_jkx'):
            lines.append(f"ri_jkx={self.collected_data['ri_jkx']}")
        if self.collected_data.get('relativity'):
            lines.append(f"relativity={self.collected_data['relativity']}")

        # Redox
        lines.append("------------------------------------")
        lines.append("Redox steps:")
        lines.append(f"calc_initial={self.collected_data.get('calc_initial', 'yes')}")
        lines.append(f"oxidation_steps={self.collected_data.get('oxidation_steps', '')}")
        lines.append(f"reduction_steps={self.collected_data.get('reduction_steps', '')}")
        lines.append(f"method={self.collected_data.get('method', 'OCCUPIER')}")

        # Resources
        lines.append("------------------------------------")
        lines.append(f"PAL={self.collected_data.get('PAL', 8)}")
        lines.append(f"maxcore={self.collected_data.get('maxcore', 4000)}")

        # Write to file
        Path(output_path).write_text("\n".join(lines))
