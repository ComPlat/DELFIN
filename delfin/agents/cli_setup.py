"""
CLI command for interactive CONTROL setup with AI assistant
"""

import argparse
import sys

from .provider_factory import create_provider, list_available_providers
from .control_assistant import ControlAssistantV2


def run_setup_command(argv: list[str] = None) -> int:
    """
    Run `delfin setup` command.

    Usage:
        delfin setup --interactive
        delfin setup --provider=ollama --model=llama3.1
        delfin setup --list-providers
    """
    parser = argparse.ArgumentParser(
        prog="delfin setup",
        description="Interactive CONTROL file creation with AI assistant"
    )
    parser.add_argument(
        "--interactive", "-i",
        action="store_true",
        help="Start interactive setup with AI assistant"
    )
    parser.add_argument(
        "--provider",
        type=str,
        default="auto",
        choices=["auto", "claude", "openai", "ollama"],
        help="LLM provider to use (default: auto-detect)"
    )
    parser.add_argument(
        "--model",
        type=str,
        help="Model name (default: provider-specific default)"
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        default="CONTROL.txt",
        help="Output path for CONTROL file (default: CONTROL.txt)"
    )
    parser.add_argument(
        "--list-providers",
        action="store_true",
        help="List available LLM providers and their status"
    )

    args = parser.parse_args(argv)

    # List providers
    if args.list_providers:
        _print_providers()
        return 0

    # Interactive setup
    if args.interactive:
        try:
            provider = create_provider(
                provider_name=args.provider,
                model=args.model
            )

            assistant = ControlAssistantV2(provider)
            assistant.create_control_interactive(output_path=args.output)
            return 0

        except Exception as exc:
            print(f"\nError: {exc}")
            print("\nTroubleshooting:")
            print("  1. Check if your LLM provider is configured:")
            print("     delfin setup --list-providers")
            print("  2. For Claude: export ANTHROPIC_API_KEY=...")
            print("  3. For OpenAI: export OPENAI_API_KEY=...")
            print("  4. For Ollama: ollama serve")
            return 1

    # No action specified
    parser.print_help()
    return 0


def _print_providers():
    """Print available providers and their status"""
    providers = list_available_providers()

    print("╔══════════════════════════════════════════════════════════╗")
    print("║          Available LLM Providers                        ║")
    print("╚══════════════════════════════════════════════════════════╝")
    print()

    for name, info in providers.items():
        status = "✓" if info.get("configured") else ("○" if info.get("available") else "✗")
        print(f"{status} {name.upper():10s}")

        if info.get("available"):
            if info.get("models"):
                print(f"   Models: {', '.join(info['models'][:3])}")
                if len(info['models']) > 3:
                    print(f"           ... and {len(info['models']) - 3} more")
            print(f"   Status: {info.get('note', 'Unknown')}")
        else:
            print(f"   {info.get('note', 'Not available')}")
        print()

    print("Legend:")
    print("  ✓ = Configured and ready")
    print("  ○ = Installed but needs configuration")
    print("  ✗ = Not installed")
    print()
    print("Quick start:")
    print("  # Use Ollama (local, free)")
    print("  ollama pull llama3.1")
    print("  delfin setup --interactive --provider=ollama")
    print()
    print("  # Use Claude (cloud, API key required)")
    print("  export ANTHROPIC_API_KEY=sk-...")
    print("  delfin setup --interactive --provider=claude")
