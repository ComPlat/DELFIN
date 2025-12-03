"""
OpenAI (GPT) Provider
"""

from __future__ import annotations
import json
import os
from typing import Any, Dict, List, Optional

from .base_provider import BaseLLMProvider, Message, StructuredOutput

try:
    from openai import OpenAI
    OPENAI_AVAILABLE = True
except ImportError:
    OPENAI_AVAILABLE = False


class OpenAIProvider(BaseLLMProvider):
    """OpenAI GPT provider with structured output support"""

    def __init__(self, model: str = "gpt-4", api_key: Optional[str] = None, **kwargs):
        if not OPENAI_AVAILABLE:
            raise RuntimeError("OpenAI SDK not installed. Install with: pip install openai")

        super().__init__(model, api_key, **kwargs)
        self.api_key = api_key or os.environ.get("OPENAI_API_KEY")
        if not self.api_key:
            raise ValueError("OPENAI_API_KEY not found in environment or constructor")

        self.client = OpenAI(api_key=self.api_key)

    def chat(
        self,
        messages: List[Message],
        schema: Optional[Dict[str, Any]] = None,
        tools: Optional[List[Dict[str, Any]]] = None,
        temperature: float = 0.0,
        **kwargs
    ) -> StructuredOutput:
        """Send chat to OpenAI with optional structured output"""

        # Convert messages to OpenAI format
        formatted_messages = [
            {"role": msg.role, "content": msg.content}
            for msg in messages
        ]

        # Build request parameters
        params = {
            "model": self.model,
            "messages": formatted_messages,
            "temperature": temperature,
        }

        # Add JSON mode or function calling for structured output
        if schema:
            # Use function calling for structured output
            function = {
                "name": "output_data",
                "description": "Output structured data",
                "parameters": schema
            }
            params["functions"] = [function]
            params["function_call"] = {"name": "output_data"}

        # Make API call
        response = self.client.chat.completions.create(**params)

        # Extract structured data
        if schema:
            # Function calling response
            message = response.choices[0].message
            if message.function_call:
                data = json.loads(message.function_call.arguments)
                if not self.validate_schema(data, schema):
                    raise ValueError(f"OpenAI response doesn't match schema: {data}")
                return StructuredOutput(
                    data=data,
                    raw_response=response,
                    provider=self.provider_name,
                    model=self.model
                )
            raise RuntimeError("OpenAI didn't call function despite function_call parameter")
        else:
            # Regular text response
            text = response.choices[0].message.content or ""
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
            # Fallback: basic validation
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
        return "openai"

    @property
    def supports_structured_output(self) -> bool:
        return True  # Via function calling

    @property
    def supports_tool_use(self) -> bool:
        return True
