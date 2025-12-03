# DELFIN AI Agent Framework

LLM-powered assistants for interactive CONTROL file generation and reporting.

## Features

✅ **Model-Agnostic**: Support for Claude, GPT-4, and Ollama (local LLMs)
✅ **Structured Output**: No hallucinations - LLMs must use predefined values
✅ **Interactive Setup**: Step-by-step CONTROL file creation
✅ **Validation**: Triple-checked against JSON schemas
✅ **Flexible**: Easy to add new providers

## Installation

```bash
# Basic installation (no AI features)
pip install delfin-complat

# With AI features
pip install delfin-complat[ai]

# Or install providers individually
pip install anthropic  # For Claude
pip install openai     # For GPT
pip install ollama     # For local LLMs
```

## Quick Start

### 1. Check Available Providers

```bash
delfin setup --list-providers
```

### 2. Interactive Setup with Ollama (Local, Free)

```bash
# Install and start Ollama
ollama pull llama3.1
ollama serve

# Run interactive setup
delfin setup --interactive --provider=ollama
```

### 3. Interactive Setup with Claude (Cloud, API Required)

```bash
# Set API key
export ANTHROPIC_API_KEY=sk-ant-...

# Run interactive setup
delfin setup --interactive --provider=claude
```

## Architecture

```
┌─────────────────────────────────────────────────────┐
│           DELFIN Agent Framework                    │
│                                                     │
│  ┌──────────────────────────────────────────────┐  │
│  │  Abstract LLM Provider Interface             │  │
│  │  - Structured output with JSON schema        │  │
│  │  - Validation layer (no hallucinations!)     │  │
│  └──────────────────────────────────────────────┘  │
│           ↓            ↓            ↓               │
│  ┌─────────────┐ ┌──────────┐ ┌─────────────┐     │
│  │   Claude    │ │  OpenAI  │ │   Ollama    │     │
│  │  Provider   │ │ Provider │ │  Provider   │     │
│  └─────────────┘ └──────────┘ └─────────────┘     │
└─────────────────────────────────────────────────────┘
```

## How It Prevents Hallucinations

### 1. JSON Schema Enforcement

```python
schema = {
    "type": "object",
    "properties": {
        "functional": {
            "type": "string",
            "enum": ["PBE0", "B3LYP", "wB97X-V", "TPSS"]  # Only these!
        },
        "charge": {
            "type": "integer",
            "minimum": -10,
            "maximum": 10  # Must be in range
        }
    }
}
```

**LLM CANNOT return values outside this schema!**

### 2. Validation Layer

```python
def validate_response(data, schema):
    # Triple-check LLM output
    if not matches_schema(data, schema):
        raise ValueError("LLM hallucinated!")
    return data
```

### 3. Tool Use / Function Calling

LLM can **only** call predefined functions with validated parameters.

## Python API

```python
from delfin.agents import create_provider, ControlAssistant

# Create provider (auto-detect)
provider = create_provider()

# Or specify explicitly
provider = create_provider(
    provider_name="ollama",
    model="llama3.1"
)

# Create assistant
assistant = ControlAssistant(provider)

# Interactive setup
config = assistant.create_control_interactive(
    output_path="CONTROL.txt"
)
```

## Adding a New Provider

```python
from delfin.agents.base_provider import BaseLLMProvider, Message, StructuredOutput

class MyProvider(BaseLLMProvider):
    def chat(self, messages, schema=None, **kwargs):
        # Your implementation
        response = self.client.chat(messages)

        # Validate against schema
        if schema and not self.validate_schema(response, schema):
            raise ValueError("Invalid response")

        return StructuredOutput(
            data=response,
            provider="my_provider",
            model=self.model
        )

    # Implement other abstract methods...
```

## Environment Variables

```bash
# Provider selection
export DELFIN_LLM_PROVIDER=ollama  # or claude, openai

# API keys
export ANTHROPIC_API_KEY=sk-ant-...
export OPENAI_API_KEY=sk-...

# Ollama host (if not localhost)
export OLLAMA_HOST=http://server:11434
```

## Report Generator (Coming Soon)

The report generator will use **pure parsing** (NO LLM) to extract data from DELFIN outputs:

```python
from delfin.agents import ReportGenerator

# Parse results (no LLM - no hallucinations!)
report = ReportGenerator.from_workspace("./")

# Get structured data
print(report.redox_potentials)  # From DELFIN.txt
print(report.geometries)         # From XYZ files
print(report.orbitals)           # From OCCUPIER results

# Generate markdown/PDF
report.to_markdown("report.md")
report.to_pdf("report.pdf")
```

**Key principle**: Report data comes **only** from files, never from LLM!

## Troubleshooting

### Ollama not found
```bash
# Install Ollama
curl -fsSL https://ollama.com/install.sh | sh

# Pull a model
ollama pull llama3.1

# Start server
ollama serve
```

### Claude API key not working
```bash
# Check key
echo $ANTHROPIC_API_KEY

# Test with curl
curl https://api.anthropic.com/v1/messages \
  -H "x-api-key: $ANTHROPIC_API_KEY" \
  -H "anthropic-version: 2023-06-01" \
  -H "content-type: application/json" \
  -d '{"model":"claude-3-5-sonnet-20241022","max_tokens":10,"messages":[{"role":"user","content":"Hi"}]}'
```

### LLM returns invalid values
This should **never** happen due to schema enforcement. If it does:
1. Check provider implementation
2. Update provider to use stricter validation
3. File a bug report

## License

Same as DELFIN (LGPL-3.0-or-later)
