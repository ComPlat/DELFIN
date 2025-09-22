# DELFIN Methodology Documentation

## Overview

DELFIN (DFT-based automated prediction of preferred spin states and associated redox potentials) is a comprehensive computational pipeline that automates density functional theory (DFT) calculations for predicting molecular spin states and redox potentials. The software interfaces with ORCA, xTB and CREST qantum chemistry package to perform systematic calculations across multiple spin states and redox conditions.

## Theoretical Foundation

### 1. Density Functional Theory Framework

DELFIN employs density functional theory (DFT) as implemented in ORCA 6.1.0 for all electronic structure calculations. The choice of exchange-correlation functional, basis sets, and computational parameters is systematically determined based on the molecular composition and user-defined requirements.

### 2. Spin State Prediction

The core methodology revolves around calculating energies for different spin multiplicities (2S+1) to identify the preferred electronic configuration:

- **Even electron systems**: Calculations typically performed for singlet (S=0), triplet (S=1), quintet (S=2) states
- **Odd electron systems**: Calculations for doublet (S=1/2), quartet (S=3/2), sextet (S=5/2) states

The ground state is determined by comparing total electronic energies or Gibbs free energies across all calculated spin states.

### 3. Redox Potential Calculations

Redox potentials are computed using the thermodynamic cycle approach:

```
E°(M^n+/M^(n-1)+) = [G(M^n+) - G(M^(n-1)+)] / nF - E_ref
```

Where:
- G(M^n+), G(M^(n-1)+): Gibbs free energies of oxidized and reduced species
- n: Number of electrons transferred
- F: Faraday constant
- E_ref: Reference electrode potential (solvent-dependent)

## Computational Methodology

### 1. Workflow Architecture

DELFIN implements three primary calculation modes:

#### OCCUPIER Method
- **Purpose**: Systematic exploration of spin states with broken-symmetry DFT
- **Implementation**: Sequential calculations with wavefunction passing between steps
- **Features**:
  - Approximate spin projection (APMethod)
  - Broken symmetry calculations for antiferromagnetic coupling
  - Per-atom basis set assignment for metal centers

#### Classic Method
- **Purpose**: Standard DFT calculations for each redox state
- **Implementation**: Independent calculations for each oxidation state
- **Features**:
  - Geometry optimization for each state
  - Frequency calculations for thermodynamic corrections
  - Excited state calculations (TD-DFT) when requested

#### Manual Method
- **Purpose**: User-defined calculation sequences
- **Implementation**: Flexible workflow based on CONTROL.txt specifications

### 2. Basis Set Selection Policy

DELFIN implements an intelligent basis set selection policy based on metal composition:

```python
# Metal classification from utils.py:61-78
def classify_tm_presence(found_metals):
    """
    Classify transition metals:
    - '3d': Only 3d metals (Sc-Zn)
    - '4d5d': Only 4d/5d metals
    - 'mixed': Both 3d and 4d/5d present
    - 'none': No transition metals
    """
```

**Policy Implementation**:
- **3d metals only**: Standard basis sets (def2-TZVP), no relativistic effects
- **4d/5d metals**: Extended basis sets (def2-TZVPP), relativistic corrections (ZORA/X2C)
- **Mixed systems**: Conservative approach favoring larger basis sets and relativity

### 3. Per-Atom Basis Set Assignment

For transition metal complexes, DELFIN employs sophisticated per-atom basis set assignment:

```python
# From occupier.py:174-296
def read_and_modify_file_OCCUPIER(...):
    """
    - Global method line with main_basisset
    - Per-atom NewGTO "metal_basisset" for metals
    - Optional per-atom NewGTO for first coordination sphere
    """
```

**First Coordination Sphere Detection**:
- Uses covalent radii from Pyykkö (2009) or Cordero (2008) datasets
- Bonds identified by: `d(M-X) ≤ scale × (r_cov(M) + r_cov(X))`
- Default scaling factor: 1.20 (configurable)

### 4. Relativity Treatment

Relativistic effects are automatically applied based on metal composition:

```python
# From utils.py:93-99
def _should_use_rel(found_metals):
    """Relativity for 4d/5d or mixed metal systems"""
    return classify_tm_presence(found_metals) in ('4d5d', 'mixed')
```

**Available Methods**:
- **ZORA**: Zeroth-order regular approximation (default for relativistic systems)
- **X2C**: Exact two-component method
- **DKH**: Douglas-Kroll-Hess transformation

### 5. Solvent Effects

Implicit solvation is handled through ORCA's continuum models:

```python
# From config.py:142-161
def get_E_ref(config):
    """Reference electrode potentials (V vs SHE)"""
    solvent_E_ref = {
        "dmf": 4.795, "dcm": 4.805, "acetonitrile": 4.745,
        "thf": 4.905, "dmso": 4.780, "acetone": 4.825
    }
```

**Implementation**:
- CPCM (Conductor-like Polarizable Continuum Model) default
- SMD (Solvation Model based on Density) available
- Solvent-specific reference electrode potentials for redox calculations

### 6. Geometry Optimization and Conformational Search

DELFIN integrates with xTB and CREST for pre-optimization and conformational analysis:

```python
# From xtb_crest.py:23-64
def XTB(multiplicity, charge, config):
    """Semi-empirical geometry optimization with xTB"""

def run_crest_workflow(...):
    """Conformational search with CREST"""
```

**Workflow**:
1. **XTB pre-optimization**: Fast geometry optimization using GFNn-xTB methods
2. **GOAT/CREST conformational search**: Systematic exploration of conformational space
3. **DFT refinement**: High-level optimization of selected conformers

### 7. Energy Extraction and Analysis

DELFIN implements robust energy parsing from ORCA output files:

```python
# From energies.py:30-51
def find_gibbs_energy(filename):
    """Extract Gibbs free energy from ORCA output"""

def find_electronic_energy(filename):
    """Extract final single point energy"""
```

**Energy Types**:
- **Electronic energy**: SCF convergence energy
- **Gibbs free energy**: Thermodynamically corrected energy (when frequencies calculated)
- **Zero-point energy**: Vibrational zero-point correction

### 8. Excited State Calculations

For photochemical applications, DELFIN supports time-dependent DFT calculations:

```python
# From energies.py:99-137
def find_state1_mit_SOC(filename):
    """Extract S1 state energy with spin-orbit coupling"""

def check_and_execute_SOC(filename, config):
    """Handle SOC-corrected excited states"""
```

**Features**:
- TD-DFT for singlet and triplet excited states
- Spin-orbit coupling corrections (optional)
- Emission spectra calculations from optimized excited states

## Validation and Quality Control

### 1. Convergence Criteria

- SCF convergence: Configurable maximum iterations (default: 125)
- Geometry optimization: ORCA default criteria
- Frequency calculations: Analytical Hessian when available

### 2. Error Handling

- Automatic detection of calculation failures
- Wavefunction passing between related calculations
- Recalculation mode for incomplete jobs

### 3. Output Validation

- Systematic checking for "ORCA TERMINATED NORMALLY" in output files
- Energy extraction with error handling for malformed outputs
- Logging of calculation progress and issues

## Implementation Details

### 1. File Management

DELFIN maintains organized directory structures:
- Input files: `input.txt`, `input2.txt`, etc.
- ORCA inputs: `.inp` files with systematic naming
- Results: `DELFIN.txt` and `OCCUPIER.txt` summary reports

### 2. Parallel Processing

- Configurable number of CPU cores (`PAL` parameter)
- Memory management (`maxcore` specification)
- Efficient parallelization within ORCA calculations

### 3. Reproducibility

- Complete parameter logging in output files
- Systematic file naming conventions
- Configuration file preservation

## References

The methodology implemented in DELFIN builds upon established computational chemistry protocols and leverages the capabilities of:

- **ORCA**: Frank Neese et al., *J. Chem. Phys.* **152**, 224108 (2020)
- **xTB methods**: Bannwarth et al., *WIREs Comput. Mol. Sci.* **11**, e1493 (2021)
- **CREST**: Pracht et al., *J. Chem. Phys.* **160**, 114110 (2024)

The specific implementation details for transition metal basis set policies, relativistic corrections, and redox potential calculations represent novel automation approaches developed for this software package.