"""
Base LLM Provider Interface

Abstract base class for all LLM providers (Claude, GPT, Ollama, etc.)
"""

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Union


@dataclass
class StructuredOutput:
    """Structured output from LLM with validation"""
    data: Dict[str, Any]
    raw_response: Any
    provider: str
    model: str


@dataclass
class Message:
    """Chat message"""
    role: str  # "user", "assistant", "system"
    content: str


class BaseLLMProvider(ABC):
    """Abstract base class for LLM providers"""

    def __init__(self, model: str, api_key: Optional[str] = None, **kwargs):
        """
        Initialize LLM provider.

        Args:
            model: Model identifier (e.g., "claude-3-5-sonnet-20241022", "gpt-4", "llama3")
            api_key: API key for cloud providers (None for local models)
            **kwargs: Provider-specific configuration
        """
        self.model = model
        self.api_key = api_key
        self.config = kwargs

    @abstractmethod
    def chat(
        self,
        messages: List[Message],
        schema: Optional[Dict[str, Any]] = None,
        tools: Optional[List[Dict[str, Any]]] = None,
        temperature: float = 0.0,
        **kwargs
    ) -> StructuredOutput:
        """
        Send chat messages and get structured response.

        Args:
            messages: List of chat messages
            schema: JSON schema for structured output (optional)
            tools: List of available tools/functions (optional)
            temperature: Sampling temperature (0.0 = deterministic)
            **kwargs: Provider-specific parameters

        Returns:
            StructuredOutput with validated data

        Raises:
            ValueError: If response doesn't match schema
            RuntimeError: If API call fails
        """
        pass

    @abstractmethod
    def validate_schema(self, data: Dict[str, Any], schema: Dict[str, Any]) -> bool:
        """
        Validate data against JSON schema.

        Args:
            data: Data to validate
            schema: JSON schema

        Returns:
            True if valid, False otherwise
        """
        pass

    @property
    @abstractmethod
    def provider_name(self) -> str:
        """Return provider name (e.g., 'claude', 'openai', 'ollama')"""
        pass

    @property
    @abstractmethod
    def supports_structured_output(self) -> bool:
        """Return True if provider natively supports structured output"""
        pass

    @property
    @abstractmethod
    def supports_tool_use(self) -> bool:
        """Return True if provider supports function/tool calling"""
        pass
