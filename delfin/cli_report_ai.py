"""
CLI command for AI-powered report generation

Usage: delfin report-ai [OPTIONS] [DIRECTORY]
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path

from delfin.common.logging import get_logger

logger = get_logger(__name__)


def run_report_ai_command(args: list[str]) -> int:
    """
    Run AI-powered report generation command.

    Args:
        args: Command-line arguments (after 'report-ai')

    Returns:
        Exit code (0 = success, non-zero = error)
    """
    parser = argparse.ArgumentParser(
        prog="delfin report-ai",
        description="Generate an AI-powered scientific report from DELFIN calculation results"
    )

    parser.add_argument(
        "directory",
        nargs="?",
        default=".",
        help="Directory containing DELFIN calculation results (default: current directory)"
    )

    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output file path (default: REPORT.txt in calculation directory)"
    )

    parser.add_argument(
        "--provider",
        default="auto",
        choices=["auto", "claude", "openai", "ollama"],
        help="LLM provider to use (default: auto-detect)"
    )

    parser.add_argument(
        "--model",
        default=None,
        help="Model to use (default: provider-specific default)"
    )

    parser.add_argument(
        "--interactive",
        "-i",
        action="store_true",
        help="Enable interactive mode to customize report content"
    )

    parser.add_argument(
        "--api-key",
        default=None,
        help="API key for the provider (overrides environment variable)"
    )

    parser.add_argument(
        "--base-url",
        default=None,
        help="Base URL for API (for custom OpenAI-compatible endpoints)"
    )

    parsed_args = parser.parse_args(args)

    # Resolve directory
    working_dir = Path(parsed_args.directory).resolve()

    if not working_dir.exists():
        logger.error(f"Directory not found: {working_dir}")
        return 1

    if not working_dir.is_dir():
        logger.error(f"Not a directory: {working_dir}")
        return 1

    # Check for DELFIN.txt or OCCUPIER.txt to verify it's a DELFIN directory
    has_delfin = (working_dir / "DELFIN.txt").exists()
    has_occupier = (working_dir / "OCCUPIER.txt").exists() or (working_dir / "initial_OCCUPIER" / "OCCUPIER.txt").exists()

    if not has_delfin and not has_occupier:
        logger.warning(
            f"Directory {working_dir} does not appear to contain DELFIN results.\n"
            f"Expected to find DELFIN.txt or OCCUPIER.txt.\n"
            f"Continuing anyway..."
        )

    # Determine output path
    if parsed_args.output:
        output_path = Path(parsed_args.output).resolve()
    else:
        output_path = working_dir / "REPORT.txt"

    # Import agent components
    try:
        from delfin.agents.provider_factory import create_provider
        from delfin.agents.report_assistant import generate_report_from_directory
    except ImportError as e:
        logger.error(f"Failed to import AI agent modules: {e}")
        logger.error("Make sure all dependencies are installed (anthropic, openai, or ollama)")
        return 1

    # Create LLM provider
    try:
        provider_kwargs = {}
        if parsed_args.api_key:
            provider_kwargs['api_key'] = parsed_args.api_key
        if parsed_args.base_url:
            provider_kwargs['base_url'] = parsed_args.base_url

        provider = create_provider(
            provider_name=parsed_args.provider,
            model=parsed_args.model,
            **provider_kwargs
        )
    except Exception as e:
        logger.error(f"Failed to create LLM provider: {e}")
        logger.error("\nTip: Set up your API key:")
        logger.error("  - For Claude: export ANTHROPIC_API_KEY=your_key")
        logger.error("  - For OpenAI: export OPENAI_API_KEY=your_key")
        logger.error("  - For Ollama: ensure ollama is running (ollama serve)")
        return 1

    # Generate report (interactive or automatic)
    try:
        if parsed_args.interactive:
            # Interactive mode
            from delfin.agents.interactive_report import generate_interactive_report

            report = generate_interactive_report(
                working_dir=working_dir,
                provider=provider,
                output_path=output_path
            )
        else:
            # Automatic mode
            report = generate_report_from_directory(
                working_dir=working_dir,
                provider=provider,
                output_path=output_path
            )

        if report:
            print(f"\nSuccess! Report generated at: {output_path}")
            return 0
        else:
            logger.error("Report generation returned empty result")
            return 1

    except Exception as e:
        logger.error(f"Report generation failed: {e}", exc_info=True)
        return 1
