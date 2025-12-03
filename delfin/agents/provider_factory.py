"""
LLM Provider Factory

Creates the appropriate provider based on configuration.
"""

from __future__ import annotations
import os
from typing import Optional

from .base_provider import BaseLLMProvider


def create_provider(
    provider_name: str = "auto",
    model: Optional[str] = None,
    api_key: Optional[str] = None,
    **kwargs
) -> BaseLLMProvider:
    """
    Create an LLM provider.

    Args:
        provider_name: Provider name ("claude", "openai", "ollama", or "auto")
        model: Model identifier (None = use provider default)
        api_key: API key (None = try environment variable)
        **kwargs: Provider-specific options

    Returns:
        Configured LLM provider

    Raises:
        ValueError: If provider is not supported or not configured
    """

    # Auto-detect provider
    if provider_name == "auto":
        provider_name = _auto_detect_provider()

    provider_name = provider_name.lower()

    # Create provider
    if provider_name == "claude":
        from .claude_provider import ClaudeProvider
        default_model = model or "claude-3-5-sonnet-20241022"
        return ClaudeProvider(model=default_model, api_key=api_key, **kwargs)

    elif provider_name == "ollama":
        from .ollama_provider import OllamaProvider
        default_model = model or "llama3.1"
        return OllamaProvider(model=default_model, api_key=api_key, **kwargs)

    elif provider_name == "openai" or provider_name == "gpt":
        from .openai_provider import OpenAIProvider
        default_model = model or "gpt-4"
        return OpenAIProvider(model=default_model, api_key=api_key, **kwargs)

    else:
        raise ValueError(
            f"Unknown provider: {provider_name}. "
            f"Supported: claude, openai, ollama"
        )


def _auto_detect_provider() -> str:
    """
    Auto-detect which provider to use based on environment.

    Priority:
    1. DELFIN_LLM_PROVIDER environment variable
    2. ANTHROPIC_API_KEY → claude
    3. OPENAI_API_KEY → openai
    4. Ollama running locally → ollama
    5. Fallback: ollama (user needs to install)

    Returns:
        Provider name
    """
    # Explicit environment variable
    env_provider = os.environ.get("DELFIN_LLM_PROVIDER")
    if env_provider:
        return env_provider.lower()

    # Check for API keys
    if os.environ.get("ANTHROPIC_API_KEY"):
        return "claude"
    if os.environ.get("OPENAI_API_KEY"):
        return "openai"

    # Check if Ollama is running
    try:
        import ollama
        # Try to list models - if this works, Ollama is running
        ollama.list()
        return "ollama"
    except Exception:
        pass

    # Default to Ollama (user will get helpful error message)
    return "ollama"


def list_available_providers() -> dict:
    """
    List all available providers and their status.

    Returns:
        Dictionary with provider info
    """
    providers = {}

    # Claude
    try:
        from .claude_provider import ClaudeProvider
        has_key = bool(os.environ.get("ANTHROPIC_API_KEY"))
        providers["claude"] = {
            "available": True,
            "configured": has_key,
            "models": ["claude-3-5-sonnet-20241022", "claude-3-opus-20240229"],
            "note": "Requires ANTHROPIC_API_KEY" if not has_key else "Ready"
        }
    except ImportError:
        providers["claude"] = {
            "available": False,
            "note": "Install with: pip install anthropic"
        }

    # OpenAI
    try:
        from .openai_provider import OpenAIProvider
        has_key = bool(os.environ.get("OPENAI_API_KEY"))
        providers["openai"] = {
            "available": True,
            "configured": has_key,
            "models": ["gpt-4", "gpt-4-turbo", "gpt-3.5-turbo"],
            "note": "Requires OPENAI_API_KEY" if not has_key else "Ready"
        }
    except ImportError:
        providers["openai"] = {
            "available": False,
            "note": "Install with: pip install openai"
        }

    # Ollama
    try:
        import ollama
        try:
            models = ollama.list()
            providers["ollama"] = {
                "available": True,
                "configured": True,
                "models": [m["name"] for m in models.get("models", [])],
                "note": f"Running locally ({len(models.get('models', []))} models installed)"
            }
        except Exception:
            providers["ollama"] = {
                "available": True,
                "configured": False,
                "note": "Installed but not running. Start with: ollama serve"
            }
    except ImportError:
        providers["ollama"] = {
            "available": False,
            "note": "Install with: pip install ollama && ollama pull llama3.1"
        }

    return providers
