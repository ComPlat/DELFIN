"""
Ollama (Local LLM) Provider
"""

from __future__ import annotations
import json
import re
from typing import Any, Dict, List, Optional

from .base_provider import BaseLLMProvider, Message, StructuredOutput

try:
    import ollama
    OLLAMA_AVAILABLE = True
except ImportError:
    OLLAMA_AVAILABLE = False


class OllamaProvider(BaseLLMProvider):
    """Ollama provider for local LLMs"""

    def __init__(self, model: str = "llama3.1", api_key: Optional[str] = None, **kwargs):
        if not OLLAMA_AVAILABLE:
            raise RuntimeError("Ollama SDK not installed. Install with: pip install ollama")

        super().__init__(model, api_key, **kwargs)
        self.host = kwargs.get("host", "http://localhost:11434")

        # Initialize Ollama client
        self.client = ollama.Client(host=self.host)

    def chat(
        self,
        messages: List[Message],
        schema: Optional[Dict[str, Any]] = None,
        tools: Optional[List[Dict[str, Any]]] = None,
        temperature: float = 0.0,
        **kwargs
    ) -> StructuredOutput:
        """Send chat to Ollama with optional structured output"""

        # Convert messages to Ollama format
        formatted_messages = [
            {"role": msg.role, "content": msg.content}
            for msg in messages
        ]

        # Add JSON schema instructions if provided
        if schema:
            schema_prompt = self._create_schema_prompt(schema)
            formatted_messages.append({
                "role": "system",
                "content": schema_prompt
            })

        # Make API call
        response = self.client.chat(
            model=self.model,
            messages=formatted_messages,
            options={
                "temperature": temperature,
            },
            format="json" if schema else None,  # Request JSON output
        )

        # Extract and parse response
        text = response["message"]["content"]

        if schema:
            # Try to parse JSON
            try:
                data = json.loads(text)
                if not self.validate_schema(data, schema):
                    # Attempt to fix common issues
                    data = self._fix_common_json_issues(data, schema)
                    if not self.validate_schema(data, schema):
                        raise ValueError(f"Ollama response doesn't match schema: {data}")

                return StructuredOutput(
                    data=data,
                    raw_response=response,
                    provider=self.provider_name,
                    model=self.model
                )
            except json.JSONDecodeError:
                # Try to extract JSON from markdown code blocks
                json_match = re.search(r'```json\s*(\{.*?\})\s*```', text, re.DOTALL)
                if json_match:
                    try:
                        data = json.loads(json_match.group(1))
                        if self.validate_schema(data, schema):
                            return StructuredOutput(
                                data=data,
                                raw_response=response,
                                provider=self.provider_name,
                                model=self.model
                            )
                    except json.JSONDecodeError:
                        pass

                raise ValueError(f"Ollama didn't return valid JSON: {text}")
        else:
            # Regular text response
            return StructuredOutput(
                data={"text": text},
                raw_response=response,
                provider=self.provider_name,
                model=self.model
            )

    def _create_schema_prompt(self, schema: Dict[str, Any]) -> str:
        """Create a prompt that explains the schema to the LLM"""
        prompt = "You MUST respond with valid JSON matching this exact schema:\n\n"
        prompt += json.dumps(schema, indent=2)
        prompt += "\n\nRULES:\n"
        prompt += "1. Output ONLY valid JSON, no markdown, no explanations\n"
        prompt += "2. All required fields must be present\n"
        prompt += "3. Use only allowed enum values if specified\n"
        prompt += "4. Respect min/max constraints\n"
        prompt += "5. Do NOT invent values - only use allowed options\n"

        # Extract constraints for emphasis
        if "properties" in schema:
            prompt += "\nALLOWED VALUES:\n"
            for key, prop in schema["properties"].items():
                if "enum" in prop:
                    prompt += f"  - {key}: {', '.join(map(str, prop['enum']))}\n"
                elif "minimum" in prop and "maximum" in prop:
                    prompt += f"  - {key}: {prop['minimum']} to {prop['maximum']}\n"

        return prompt

    def _fix_common_json_issues(self, data: Dict[str, Any], schema: Dict[str, Any]) -> Dict[str, Any]:
        """Try to fix common JSON issues from local LLMs"""
        fixed = data.copy()

        if "properties" not in schema:
            return fixed

        for key, prop_schema in schema["properties"].items():
            if key not in fixed:
                continue

            value = fixed[key]
            prop_type = prop_schema.get("type")

            # Try to coerce to correct type
            if prop_type == "integer" and isinstance(value, str):
                try:
                    fixed[key] = int(value)
                except ValueError:
                    pass
            elif prop_type == "integer" and isinstance(value, float):
                fixed[key] = int(value)
            elif prop_type == "string" and not isinstance(value, str):
                fixed[key] = str(value)

        return fixed

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
        return "ollama"

    @property
    def supports_structured_output(self) -> bool:
        return True  # Via format="json"

    @property
    def supports_tool_use(self) -> bool:
        return False  # Ollama doesn't support tool use yet
