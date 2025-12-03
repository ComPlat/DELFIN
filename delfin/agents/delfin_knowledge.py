"""
DELFIN Knowledge Base for AI Assistant

This module contains comprehensive documentation about all DELFIN features,
parameters, and best practices to enable the AI assistant to provide
expert guidance.
"""

DELFIN_KNOWLEDGE = """
# DELFIN - Comprehensive Knowledge Base

## Overview
DELFIN is an advanced workflow manager for computational chemistry calculations using ORCA, xTB, and CREST.
It specializes in redox potential calculations, excited state calculations, and automated spin state exploration.

## Key Features

### 1. Redox Potential Calculations
- Oxidation potentials (E_ox)
- Reduction potentials (E_red)
- Excited state redox potentials (*E_ox, *E_red)
- Multiple oxidation/reduction steps
- Three calculation methods: classic, manually, OCCUPIER

### 2. Excited State Calculations
- E_00 energies (singlet and triplet)
- Absorption spectra
- Emission spectra
- Intersystem crossing (ISC) rates
- Internal conversion (IC) rates
- TDDFT and deltaSCF methods

### 3. Automatic Spin State Exploration (OCCUPIER)
- Automatically determines optimal spin states
- Handles complex transition metal systems
- Explores broken-symmetry solutions
- Configurable exploration depth (deep, deep2, deep3, flat)

## Workflow Steps (in order)

### Step 1: Charge
- Total molecular charge
- Examples: 0 (neutral), +1 (cation), +2 (dicationic complex), -1 (anion)

### Step 2: Solvent
- Implicit solvation models: CPCM, SMD, ALPB
- Common solvents: Water, MeCN, DMF, DMSO, Acetone, THF, DCM, etc.
- CPCM: General-purpose, reliable
- SMD: More accurate but slower
- ALPB: Fast, compatible with xTB

### Step 3: Structure Quality
- Ask if the structure is an approximation or reliable geometry
- If approximation → Recommend XTB_OPT=yes and XTB_GOAT=yes
- XTB_OPT: Pre-optimization with xTB (fast, improves geometry)
- XTB_GOAT: xTB GOAT optimizer (robust for difficult cases)

### Step 4: Imaginary Frequencies (IMAG)
- IMAG=yes: Check for and eliminate imaginary frequencies
- IMAG eliminates transition states and ensures true minima
- IMAG_option=2: Recommended (distorts along imaginary mode and reoptimizes)
- IMAG_scope=initial: Only check initial structure
- Why: Imaginary frequencies indicate transition states, not stable minima
  → This causes incorrect energies and failed calculations

### Step 5: Redox Potentials Setup
- calc_initial: Calculate initial state before oxidation/reduction
- oxidation_steps: e.g., "1,2,3" for three oxidation steps
- reduction_steps: e.g., "1,2" for two reduction steps

### Step 6: Method Selection (CRITICAL)
Three methods available:

#### A) OCCUPIER (recommended for transition metal complexes)
- Automatically determines spin states
- Explores broken-symmetry solutions
- Handles open-shell systems
- Use when: Transition metal complexes, spin state is not trivial
- Recommendation: "For your iron/cobalt/nickel complex, OCCUPIER is recommended
  because spin states are non-trivial and multiple configurations are possible."

#### B) classic (recommended for TADF organic molecules)
- Simple SCF calculation with single spin state
- User must know the correct multiplicity
- Use when: Closed-shell organic molecules, TADF emitters
- Recommendation: "For TADF organic molecules, classic method is recommended
  if you know the spin state. For closed-shell systems, use multiplicity=1."

#### C) manually (expert mode)
- User defines spin states and broken-symmetry flags manually
- Full control over electronic configuration
- Use when: Expert user wants specific configurations
- Requires detailed configuration for each redox state
- When selected → Enter interactive MANUALLY configuration mode

### Step 7: MANUALLY Configuration (if method=manually)
Interactive setup for each state:
- multiplicity_0: Initial state multiplicity
- additions_0: Additional ORCA keywords for initial state
- additions_TDDFT: Keywords for TDDFT calculations
- additions_T1: Keywords for T1 state
- additions_S1: Keywords for S1 state
- multiplicity_ox1/ox2/ox3: Multiplicities for oxidation states
- additions_ox1/ox2/ox3: Additional keywords for oxidation states
- multiplicity_red1/red2/red3: Multiplicities for reduction states
- additions_red1/red2/red3: Additional keywords for reduction states

Example MANUALLY configuration:
```
multiplicity_0=1
additions_0=%scf BrokenSym 1,1 end
multiplicity_ox1=2
additions_ox1=%scf BrokenSym 2,1 end
multiplicity_red1=2
additions_red1=%scf BrokenSym 2,1 end
```

### Step 8: Excited State Redox Potentials / E_00
- Ask: Do you want to calculate excited state redox potentials or E_00 energies?
- IMPORTANT: Only available for closed-shell systems
- E_00=yes: Calculate E_00 energy (singlet-triplet gap)
- excitation: "s" (singlet), "t" (triplet), or "s,t" (both)
- *E_ox / *E_red: Excited state redox potentials

Options when E_00=yes:
- Singlet E_00: Energy gap between S0 and S1
- Triplet E_00: Energy gap between S0 and T1
- Both: Calculate both S1 and T1 E_00 energies

### Step 9: Level of Theory
Key parameters:
- functional: DFT functional (PBE0, wB97X-V, B3LYP, M06-2X, etc.)
- main_basisset: Basis set for light atoms (def2-SVP, def2-TZVP)
- metal_basisset: Basis set for metals (def2-TZVP, SARC-ZORA-TZVP)
- disp_corr: Dispersion correction (D4, D3BJ, D3)
- ri_jkx: RI approximation (RIJCOSX, RIJK)
- relativity: Relativistic treatment (ZORA, DKH, X2C, none)

Recommendations based on system:
- Transition metals with 3d elements (Fe, Co, Ni): PBE0, def2-TZVP, ZORA
- Transition metals with 4d/5d elements: ZORA or DKH required
- Organic molecules: PBE0 or wB97X-V, def2-TZVP
- TADF emitters: wB97X-V or M06-2X, def2-TZVP
- Fast screening: PBE0, def2-SVP
- High accuracy: wB97X-V or DLPNO-CCSD(T), def2-TZVP or def2-QZVP

### Step 10: Additional Analysis
- print_MOs=yes: Generate molecular orbital diagrams
- print_Loewdin_population_analysis=yes: Mulliken/Loewdin population analysis
  → Useful for understanding electron distribution

### Step 11: Computational Resources
- PAL: Number of CPU cores (e.g., 8, 12, 16, 32)
- maxcore: Memory per core in MB (e.g., 4000 = 4 GB per core)
- Total memory = PAL × maxcore (e.g., 12 cores × 6000 MB = 72 GB)

Advanced parallelization:
- parallel_workflows=yes: Run multiple calculations in parallel
- pal_jobs: Number of parallel jobs (e.g., 4)
- Calculation: PAL_per_job = PAL / pal_jobs
  Example: PAL=12, pal_jobs=3 → Each job uses 4 cores

### Step 12: Automatic Recovery & Error Handling
- enable_auto_recovery=yes: Automatically retry failed calculations
- max_recovery_attempts: Number of retry attempts (default: 1)
- enable_adaptive_parallelism=yes: Adjust parallelism if memory errors occur
- enable_job_timeouts=yes: Kill stuck calculations
- job_timeout_hours: Maximum time for any job (default: 36h)

Why use auto recovery?
- SCF convergence failures → Retry with different algorithm
- Memory errors → Reduce parallelism
- Geometry optimization failures → Adjust optimizer settings

## OCCUPIER-Specific Settings

### OCCUPIER_method
- auto: DELFIN automatically determines the best approach
- manually: User provides custom spin state sequences

### OCCUPIER_tree
- deep: Standard exploration (9 configurations per redox state)
- deep2: Extended exploration (more broken-symmetry solutions)
- deep3: Exhaustive exploration (all reasonable configurations)
- flat: Minimal exploration (only high-spin and low-spin)
- own: Custom user-defined sequence

### OCCUPIER_sequence_profiles
Defines the sequence of spin states and broken-symmetry configurations:
- even electron number: Singlet (m=1), Triplet (m=3), Quintet (m=5), ...
- odd electron number: Doublet (m=2), Quartet (m=4), Sextet (m=6), ...
- BS flag: Broken-symmetry notation (e.g., "1,1" means 1 unpaired α, 1 unpaired β)

Example even_seq (for Fe²⁺, d6 system):
```
{"index": 1, "m": 1, "BS": "",    "from": 0}     # Singlet
{"index": 2, "m": 1, "BS": "1,1", "from": 1}     # Singlet broken-symmetry
{"index": 4, "m": 3, "BS": "",    "from": 1}     # Triplet
{"index": 7, "m": 5, "BS": "",    "from": 4}     # Quintet (often ground state for Fe²⁺)
```

### occupier_selection
- tolerance: Select all states within energy threshold (epsilon)
- truncation: Select top N states
- rounding: Round to nearest integer

### frequency_calculation_OCCUPIER
- yes: Calculate frequencies for all OCCUPIER-found states
- no: Skip frequency calculations (faster)

## ESD Module (Excited State Dynamics)

Enable with: ESD_modul=yes

Calculates:
- ISC rates: Intersystem crossing (e.g., S1→T1)
- IC rates: Internal conversion (e.g., S1→S0)
- Spin-orbit coupling matrix elements

Configuration:
- ESD_modus: "TDDFT" or "deltaSCF"
- states: e.g., "S0,S1,T1,T2"
- ISCs: e.g., "S1>T1,T1>S1"
- ICs: e.g., "S1>S0,T2>T1"
- DOSOC=TRUE: Calculate spin-orbit coupling

## Common Calculation Types

### 1. Redox Potentials for Transition Metal Complex
- method=OCCUPIER
- OCCUPIER_tree=deep
- functional=PBE0
- relativity=ZORA
- IMAG=yes

### 2. TADF Emitter (Organic Molecule)
- method=classic
- E_00=yes, excitation=s
- functional=wB97X-V or M06-2X
- IMAG=yes
- ESD_modul=yes (for ISC rates)

### 3. Excited State Redox Potentials
- E_00=yes
- *E_ox / *E_red experimental values
- Must be closed-shell system
- method=classic or OCCUPIER

### 4. High-Spin vs Low-Spin Energy Difference
- method=OCCUPIER
- OCCUPIER_tree=deep
- Compare energies of different spin states

## Troubleshooting & Best Practices

### SCF Convergence Issues
- Increase maxiter (default: 125 → try 250)
- Enable auto_recovery
- Use better initial guess (initial_guess=PModel)

### Memory Errors
- Reduce PAL
- Reduce maxcore
- Enable adaptive_parallelism

### Imaginary Frequencies
- IMAG=yes, IMAG_option=2
- Reoptimize structure
- Check if structure is transition state

### OCCUPIER Not Finding Good States
- Try OCCUPIER_tree=deep2 or deep3
- Adjust occupier_epsilon (default: 5e-4)
- Increase maxiter_occupier

### Geometry Optimization Not Converging
- Use XTB_OPT=yes first
- Increase opt_timeout_hours
- Check for symmetry breaking
"""

# System type recommendations
SYSTEM_TYPE_RECOMMENDATIONS = {
    "transition_metal_complex": {
        "method": "OCCUPIER",
        "functional": "PBE0",
        "relativity": "ZORA",
        "metal_basisset": "SARC-ZORA-TZVP",
        "reason": "Transition metal complexes have non-trivial spin states. OCCUPIER automatically explores all relevant configurations."
    },
    "tadf_organic": {
        "method": "classic",
        "functional": "wB97X-V",
        "main_basisset": "def2-TZVP",
        "E_00": "yes",
        "ESD_modul": "yes",
        "reason": "TADF emitters are typically closed-shell organic molecules. Use classic method with excited state calculations."
    },
    "organic_closed_shell": {
        "method": "classic",
        "functional": "PBE0",
        "main_basisset": "def2-TZVP",
        "reason": "Closed-shell organic molecules are straightforward. Classic method with singlet ground state."
    },
    "radical": {
        "method": "classic",
        "functional": "PBE0",
        "main_basisset": "def2-TZVP",
        "reason": "Organic radicals typically have well-defined doublet ground state."
    }
}
