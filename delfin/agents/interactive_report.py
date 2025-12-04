"""
Interactive Report Generation

Allows users to guide the AI agent in creating customized reports by specifying
what additional information should be included.
"""

from __future__ import annotations
from pathlib import Path
from typing import List, Optional

from .base_provider import BaseLLMProvider, Message
from .report_parser import DELFINReportData, ReportParser


INTERACTIVE_SYSTEM_PROMPT = """You are an interactive scientific report assistant for computational chemistry.

Your task is to help users create customized reports by:
1. Understanding what additional information they want in the report
2. Searching through calculation output files to find that information
3. Incorporating it into the scientific report

You have access to various output files:
- DELFIN.txt (summary of calculations)
- OCCUPIER.txt (orbital and energy data)
- ORCA .out files (detailed calculation results)
- ESD directory (excited state calculations)

When a user requests specific information:
1. Acknowledge what they're looking for
2. Explain where you'll search for it
3. Report what you found (or didn't find)
4. Suggest the format for including it in the report

Be precise and factual. Never invent data.
"""


class InteractiveReportGenerator:
    """Interactive report generation with user guidance"""

    def __init__(self, provider: BaseLLMProvider, working_dir: Path):
        """
        Initialize interactive report generator.

        Args:
            provider: LLM provider
            working_dir: Directory containing calculation results
        """
        self.provider = provider
        self.working_dir = Path(working_dir)
        self.conversation_history: List[Message] = []
        self.report_data: Optional[DELFINReportData] = None
        self.custom_sections: List[str] = []

        # Initialize with system prompt
        self.conversation_history.append(
            Message(role="system", content=INTERACTIVE_SYSTEM_PROMPT)
        )

    def start_interactive_session(self) -> DELFINReportData:
        """
        Start interactive session to customize report content.

        Returns:
            DELFINReportData with additional custom sections
        """
        print("╔══════════════════════════════════════════════════════════╗")
        print("║   Interactive DELFIN Report Generator                    ║")
        print("╚══════════════════════════════════════════════════════════╝")
        print()

        # First, extract standard data
        print("Extracting standard calculation data...")
        self.report_data = ReportParser.extract_calculation_summary(self.working_dir)
        print("✓ Standard data extracted\n")

        # Show what was found
        print("Found the following standard information:")
        self._print_summary()
        print()

        # Interactive loop
        print("You can now request additional information for the report.")
        print("Examples:")
        print("  - 'Add information about spin contamination'")
        print("  - 'Include bond lengths'")
        print("  - 'Add NBO charges'")
        print("  - 'Include all warnings from the calculations'")
        print()
        print("Type 'done' when finished, or 'skip' to use standard report only.\n")

        while True:
            user_request = input("What else should be included? ").strip()

            if not user_request or user_request.lower() == 'done':
                print("\n✓ Interactive session complete!")
                break

            if user_request.lower() == 'skip':
                print("\n✓ Using standard report only.")
                break

            # Process user request with AI
            response = self._process_user_request(user_request)
            print(f"\nAgent: {response}\n")

        return self.report_data

    def _print_summary(self):
        """Print summary of extracted data"""
        data = self.report_data

        if data.charge is not None:
            print(f"  • Charge: {data.charge}, Multiplicity: {data.multiplicity}")

        if data.geometry and data.geometry.functional:
            print(f"  • Method: {data.geometry.functional}/{data.geometry.basis_set}")

        if data.geometry and data.geometry.final_energy_ev:
            print(f"  • Final Energy: {data.geometry.final_energy_ev:.2f} eV")

        if data.orbitals and data.orbitals.homo_ev:
            print(f"  • HOMO/LUMO: {data.orbitals.homo_ev:.2f} / {data.orbitals.lumo_ev:.2f} eV")

        if data.redox and (data.redox.e_ox is not None or data.redox.e_red is not None):
            redox_parts = []
            if data.redox.e_ox is not None:
                redox_parts.append(f"E_ox = {data.redox.e_ox:.3f} V")
            if data.redox.e_red is not None:
                redox_parts.append(f"E_red = {data.redox.e_red:.3f} V")
            if redox_parts:
                print(f"  • Redox: {', '.join(redox_parts)}")

        if data.excited_states and data.excited_states.num_states > 0:
            print(f"  • Excited states: {data.excited_states.num_states} transitions")

        if data.status and data.status.has_errors:
            print(f"  ⚠ Errors detected: {len(data.status.errors)}")

        if data.status and data.status.has_warnings:
            print(f"  ⚠ Warnings detected: {len(data.status.warnings)}")

    def _process_user_request(self, request: str) -> str:
        """
        Process user request with AI assistance.

        Args:
            request: User's request for additional information

        Returns:
            AI response
        """
        # Add user message to conversation
        self.conversation_history.append(
            Message(role="user", content=self._build_search_prompt(request))
        )

        # Get AI response
        try:
            response = self.provider.chat(
                messages=self.conversation_history,
                temperature=0.0
            )

            assistant_message = response.data.get('content', 'No response generated.')

            # Add assistant response to history
            self.conversation_history.append(
                Message(role="assistant", content=assistant_message)
            )

            # Store custom section
            self.custom_sections.append(f"User Request: {request}\nAgent Response: {assistant_message}")

            return assistant_message

        except Exception as e:
            return f"Error processing request: {e}"

    def _build_search_prompt(self, user_request: str) -> str:
        """
        Build prompt for AI to search for requested information.

        Args:
            user_request: What the user wants to add

        Returns:
            Formatted prompt
        """
        # Get list of available files
        available_files = []

        for pattern in ['*.out', '*.txt', 'ESD/*.out']:
            files = list(self.working_dir.glob(pattern))
            available_files.extend([str(f.relative_to(self.working_dir)) for f in files])

        prompt = f"""User Request: "{user_request}"

Available files in calculation directory ({self.working_dir.name}):
{chr(10).join('  - ' + f for f in available_files[:20])}

Your task:
1. Identify which files likely contain the requested information
2. Explain what to search for (keywords, patterns)
3. If you can find it in the standard extracted data, mention that
4. Suggest how to present this information in the report

Respond in a helpful, concise manner (2-4 sentences).
"""
        return prompt

    def get_custom_sections(self) -> str:
        """
        Get custom sections as formatted text.

        Returns:
            Formatted custom sections
        """
        if not self.custom_sections:
            return ""

        return "\n\n".join(self.custom_sections)


def generate_interactive_report(
    working_dir: Path | str,
    provider: BaseLLMProvider,
    output_path: Optional[Path | str] = None
) -> str:
    """
    Generate report with interactive user guidance.

    Args:
        working_dir: Calculation directory
        provider: LLM provider
        output_path: Output file path

    Returns:
        Generated report text
    """
    working_dir = Path(working_dir)

    if output_path is None:
        output_path = working_dir / "REPORT_INTERACTIVE.txt"
    else:
        output_path = Path(output_path)

    # Start interactive session
    generator = InteractiveReportGenerator(provider, working_dir)
    report_data = generator.start_interactive_session()

    # Generate report with custom sections
    from .report_assistant import ReportAssistant

    assistant = ReportAssistant(provider)

    # Add custom sections to conversation if any
    if generator.custom_sections:
        print("\nGenerating customized report with your additions...")
        custom_info = generator.get_custom_sections()

        # Create custom prompt
        custom_prompt = f"""Generate a scientific report that includes both standard data and the following custom sections requested by the user:

{custom_info}

Make sure to incorporate all requested information naturally into the report structure.
"""
        # TODO: Pass custom prompt to report assistant

    else:
        print("\nGenerating standard report...")

    report = assistant.generate_report(report_data, output_path)

    return report
