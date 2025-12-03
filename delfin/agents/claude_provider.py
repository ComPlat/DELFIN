"""
Claude (Anthropic) LLM Provider
"""

from __future__ import annotations
import json
import os
from typing import Any, Dict, List, Optional

from .base_provider import BaseLLMProvider, Message, StructuredOutput

try:
    from anthropic import Anthropic
    ANTHROPIC_AVAILABLE = True
except ImportError:
    ANTHROPIC_AVAILABLE = False


class ClaudeProvider(BaseLLMProvider):
    """Claude (Anthropic) provider with structured output support"""

    def __init__(self, model: str = "claude-3-5-sonnet-20241022", api_key: Optional[str] = None, **kwargs):
        if not ANTHROPIC_AVAILABLE:
            raise RuntimeError("Anthropic SDK not installed. Install with: pip install anthropic")

        super().__init__(model, api_key, **kwargs)
        self.api_key = api_key or os.environ.get("ANTHROPIC_API_KEY")
        if not self.api_key:
            raise ValueError("ANTHROPIC_API_KEY not found in environment or constructor")

        self.client = Anthropic(api_key=self.api_key)

    def chat(
        self,
        messages: List[Message],
        schema: Optional[Dict[str, Any]] = None,
        tools: Optional[List[Dict[str, Any]]] = None,
        temperature: float = 0.0,
        **kwargs
    ) -> StructuredOutput:
        """Send chat to Claude with optional structured output"""

        # Convert messages to Anthropic format
        system_msg = None
        formatted_messages = []
        for msg in messages:
            if msg.role == "system":
                system_msg = msg.content
            else:
                formatted_messages.append({"role": msg.role, "content": msg.content})

        # Build request parameters
        params = {
            "model": self.model,
            "messages": formatted_messages,
            "temperature": temperature,
            "max_tokens": kwargs.get("max_tokens", 4096),
        }

        if system_msg:
            params["system"] = system_msg

        # Add tools if schema is provided (Claude doesn't have native JSON mode yet)
        if schema:
            tool = {
                "name": "output_structured_data",
                "description": "Output structured data matching the required schema",
                "input_schema": schema
            }
            params["tools"] = [tool]
            params["tool_choice"] = {"type": "tool", "name": "output_structured_data"}

        # Make API call
        response = self.client.messages.create(**params)

        # Extract structured data
        if schema:
            # Tool use response
            for content_block in response.content:
                if content_block.type == "tool_use":
                    data = content_block.input
                    if not self.validate_schema(data, schema):
                        raise ValueError(f"Claude response doesn't match schema: {data}")
                    return StructuredOutput(
                        data=data,
                        raw_response=response,
                        provider=self.provider_name,
                        model=self.model
                    )
            raise RuntimeError("Claude didn't use tool despite tool_choice")
        else:
            # Regular text response
            text = response.content[0].text if response.content else ""
            return StructuredOutput(
                data={"text": text},
                raw_response=response,
                provider=self.provider_name,
                model=self.model
            )

    def validate_schema(self, data: Dict[str, Any], schema: Dict[str, Any]) -> bool:
        """Validate data against JSON schema"""
        try:
            import jsonschema
            jsonschema.validate(instance=data, schema=schema)
            return True
        except ImportError:
            # Fallback: basic type checking
            return self._basic_schema_validation(data, schema)
        except Exception:
            return False

    def _basic_schema_validation(self, data: Dict[str, Any], schema: Dict[str, Any]) -> bool:
        """Basic schema validation without jsonschema library"""
        if "properties" not in schema:
            return True

        for key, prop_schema in schema["properties"].items():
            if key not in data:
                if key in schema.get("required", []):
                    return False
                continue

            value = data[key]
            prop_type = prop_schema.get("type")

            # Check type
            if prop_type == "integer" and not isinstance(value, int):
                return False
            if prop_type == "string" and not isinstance(value, str):
                return False
            if prop_type == "boolean" and not isinstance(value, bool):
                return False

            # Check enum
            if "enum" in prop_schema and value not in prop_schema["enum"]:
                return False

            # Check min/max
            if "minimum" in prop_schema and value < prop_schema["minimum"]:
                return False
            if "maximum" in prop_schema and value > prop_schema["maximum"]:
                return False

        return True

    @property
    def provider_name(self) -> str:
        return "claude"

    @property
    def supports_structured_output(self) -> bool:
        return True  # Via tool use

    @property
    def supports_tool_use(self) -> bool:
        return True
