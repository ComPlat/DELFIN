# DELFIN AI-Powered Report Generation

## Overview

The AI-powered report generator automatically creates scientific summaries of DELFIN calculations in academic style. It extracts data from output files and uses an LLM to generate a factual, well-structured report.

## Features

- **Automatic Data Extraction**: Parses DELFIN.txt, OCCUPIER.txt, and ORCA output files
- **Factual Reporting**: Generates reports based only on extracted data (no hallucination)
- **Multiple LLM Providers**: Supports Claude, OpenAI, and Ollama
- **Academic Style**: Produces publication-ready text in past tense with proper units
- **Comprehensive Coverage**: Includes geometry, energies, orbitals, frequencies, excited states, and redox potentials

## Usage

### Basic Usage

Generate a report from a DELFIN calculation directory:

```bash
delfin report-ai /path/to/calculation
```

This will:
1. Extract data from all output files
2. Generate a scientific report using AI
3. Save it as `REPORT.txt` in the calculation directory

### Advanced Options

```bash
# Specify output file
delfin report-ai /path/to/calculation -o my_report.txt

# Choose LLM provider
delfin report-ai /path/to/calculation --provider claude
delfin report-ai /path/to/calculation --provider openai
delfin report-ai /path/to/calculation --provider ollama

# Specify model
delfin report-ai /path/to/calculation --provider claude --model claude-3-5-sonnet-20241022
```

## Setup

### 1. Install Dependencies

The report generator requires at least one LLM provider:

```bash
# For Claude (recommended)
pip install anthropic
export ANTHROPIC_API_KEY="your_api_key"

# For OpenAI
pip install openai
export OPENAI_API_KEY="your_api_key"

# For Ollama (local, free)
pip install ollama
ollama serve  # Start Ollama server
ollama pull llama3.1  # Download model
```

### 2. Check Available Providers

```bash
delfin setup
```

This will show which LLM providers are available and configured.

## Example Report

Here's an example of a generated report:

```
A total of 4 conformers were identified following global geometry optimization
using the GFN2-xTB method in combination with the GOAT algorithm. Subsequent
calculations included geometry optimization, vibrational frequency analysis, and
excited-state computations. All calculations were carried out employing the ORCA,
XTB, and DELFIN software packages at the PBE0/def2-SVP with def2/J, D4, RIJCOSX,
and CPCM(MeCN) level of theory.

The final single-point energy of the optimized structure was determined to be
−77,827.38 eV. The energies of the highest occupied molecular orbital (HOMO)
and lowest unoccupied molecular orbital (LUMO) were calculated as -5.24 eV and
-1.51 eV, respectively, yielding a HOMO–LUMO gap of 3.73 eV.

The most intense vibrational modes were found at 1486, 1634, 1672, and 3144 cm⁻¹,
all corresponding to in-plane symmetric vibrations. Notably, no imaginary
(negative) frequencies were observed, confirming that the optimized structure
corresponds to a local minimum on the potential energy surface.

Excited-state calculations revealed 25 electronic transitions, including both
singlet and triplet states. The most intense absorption peaks were predicted at
259 nm and 410 nm. The S₀ → S₁ vertical excitation energy was calculated to be
3.03 eV (410 nm), while the S₀ → T₁ excitation occurred at 2.42 eV (513 nm).

Emission energies were also computed: the S₁ → S₀ fluorescence transition
corresponds to 2.84 eV (438 nm), and the T₁ → S₀ phosphorescence transition
to 2.20 eV (564 nm). The calculated E00 energy is 2.96 eV.
```

## Data Sources

The report generator extracts data from:

### DELFIN.txt
- Molecular charge and multiplicity
- Computational method details
- Redox potentials (E_ox, E_red)

### OCCUPIER.txt
- Final single-point energies
- HOMO/LUMO energies and gap
- Spin contamination data

### ORCA .out files
- Vibrational frequencies
- Imaginary modes
- Normal mode descriptions

### ESD directory
- Excited state energies
- Absorption spectra
- Emission energies
- ISC and IC data

## Customization

### Programmatic Usage

```python
from pathlib import Path
from delfin.agents import create_provider, generate_report_from_directory

# Create LLM provider
provider = create_provider(provider_name="claude")

# Generate report
report = generate_report_from_directory(
    working_dir=Path("/path/to/calculation"),
    provider=provider,
    output_path=Path("my_report.txt")
)

print(report)
```

### Custom Parsing

```python
from delfin.agents.report_parser import ReportParser

# Extract data manually
report_data = ReportParser.extract_calculation_summary(working_dir)

# Access extracted data
print(f"HOMO: {report_data.orbitals.homo_ev} eV")
print(f"LUMO: {report_data.orbitals.lumo_ev} eV")
print(f"E_ox: {report_data.redox.e_ox} V")
```

## Architecture

### report_parser.py
- Extracts factual data from output files
- No interpretation or analysis
- Returns structured `DELFINReportData` objects

### report_assistant.py
- Uses LLM to generate scientific text
- Enforces factual reporting (no hallucination)
- Formats data in academic style
- Temperature = 0.0 for deterministic output

### cli_report_ai.py
- Command-line interface
- Provider auto-detection
- Error handling and validation

## Best Practices

1. **Always verify generated reports**: While the AI is constrained to factual data, always review the output
2. **Use appropriate models**: Claude or GPT-4 recommended for best quality
3. **Provide complete data**: Ensure all relevant output files are present
4. **Check units**: The report includes units, but verify they match your expectations

## Troubleshooting

### "No LLM provider available"
- Install at least one provider: `pip install anthropic` or `pip install openai`
- Set API key: `export ANTHROPIC_API_KEY=your_key`
- Or use Ollama locally: `ollama serve && ollama pull llama3.1`

### "Directory does not contain DELFIN results"
- Ensure DELFIN.txt or OCCUPIER.txt exists in the directory
- Check that the calculation completed successfully

### "Report generation failed"
- Check your API key is valid
- Ensure you have internet connection (for cloud providers)
- Try a different provider: `--provider ollama` for local generation

## Future Enhancements

Planned features:
- Integration into main pipeline (automatic report at end)
- Support for multiple conformers
- Comparison between calculated and experimental values
- Export to LaTeX/Word formats
- Citation generation
