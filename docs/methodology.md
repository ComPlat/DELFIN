# DELFIN Methodology Documentation

## Overview

DELFIN (DFT-based automated prediction of preferred spin states and associated redox potentials) is a fully automated workflow that determines preferred spin states and computes ground-state redox potentials for open- and closed-shell systems, encompassing up to three sequential reduction and oxidation steps. For closed-shell species, it further provides excited-state redox potentials and E₀₀ energies. Starting from 2D structures or SMILES strings, DELFIN integrates GOAT/CREST with xTB for global conformer exploration and refines geometries with high-level DFT in ORCA, including automated energy and wavefunction analyses to resolve preferred configurations and track orbital occupations across redox sequences.

## Theoretical Foundation

### 1. Density Functional Theory Framework

DELFIN employs density functional theory (DFT) as implemented in ORCA 6.1.0 for all electronic structure calculations. The choice of exchange-correlation functional, basis sets, and computational parameters is systematically determined based on the molecular composition and user-defined requirements.

### 2. Spin State Prediction

The core methodology revolves around calculating energies for different spin multiplicities (2S+1) to identify the preferred electronic configuration:

- **Even electron systems**: Calculations typically performed for singlet (S=0), triplet (S=1), quintet (S=2) states
- **Odd electron systems**: Calculations for doublet (S=1/2), quartet (S=3/2), sextet (S=5/2) states

The ground state is determined by comparing total electronic energies or Gibbs free energies across all calculated spin states.

### 3. Redox Potential Calculations

Redox potentials are computed using the thermodynamic cycle approach, following the methodology established by Marenich et al. (2014):

```
E°(M^n+/M^(n-1)+) = [G(M^n+) - G(M^(n-1)+)] / nF - E_ref
```

Where:
- G(M^n+), G(M^(n-1)+): Gibbs free energies of oxidized and reduced species
- n: Number of electrons transferred
- F: Faraday constant (96,485 C/mol)
- E_ref: Reference electrode potential (solvent-dependent)

**Reference Electrode Potentials (V vs SHE):**
- DMF: 4.795 V
- DCM: 4.805 V
- Acetonitrile: 4.745 V
- THF: 4.905 V
- DMSO: 4.780 V
- Acetone: 4.825 V

By default, DELFIN reports potentials versus the ferrocene/ferrocenium (Fc/Fc⁺) couple, applying solvent-specific corrections based on established experimental correlations.

## Computational Methodology

### 1. Workflow Architecture

DELFIN implements three primary calculation modes:

#### OCCUPIER Method
- **Purpose**: Systematic exploration of spin states with broken-symmetry DFT for transition metal complexes and systems with non-innocent ligands
- **Implementation**: Sequential calculations with wavefunction passing between steps
- **Theoretical Framework**: Based on broken-symmetry (BS) DFT approach widely used to describe antiferromagnetic coupling between metal-centered d electrons and ligand-centered π* electrons
- **Features**:
  - Approximate spin projection using Noodleman/Ruiz/Yamaguchi methods
  - Broken symmetry calculations for systems with ⟨S²⟩ ≈ 1 → BS(1,1), ⟨S²⟩ ≈ 2 → BS(1,2)
  - Per-atom basis set assignment for metal centers and first coordination sphere
  - Automatic electronic state assignment using energy and spin contamination analysis

#### Classic Method
- **Purpose**: Standard DFT calculations for each redox state
- **Implementation**: Independent calculations for each oxidation state
- **Features**:
  - Geometry optimization for each state
  - Frequency calculations for thermodynamic corrections
  - Excited state calculations (TD-DFT) when requested
  - Shared orchestration in `parallel_classic_manually.py` with automatic multiplicity detection from the current electron count

#### Manual Method
- **Purpose**: User-defined calculation sequences
- **Implementation**: Flexible workflow based on CONTROL.txt specifications
  - Uses the same parallel orchestrator but honours explicit `multiplicity_*` and `additions_*` entries from `CONTROL.txt`

### Input File Preparation

- `delfin.common.paths.resolve_path` normalises all filesystem paths (user expansion plus absolute paths) so CONTROL and geometry assets remain reproducible across runs.
- `convert_xyz_to_input_txt` and `create_control_file` emit paired console messages and log records so that batch experiments retain a human-readable audit trail without changing existing behaviour.
- `normalize_input_file` converts CONTROL-referenced `.xyz` inputs into ORCA-ready `.txt` geometries prior to starting the pipeline.

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
- **3d metals only**: Non-relativistic def2-SVP family as global main basis. Optional per-metal overrides are taken from the `metal_basisset` control key when provided.
- **4d/5d metals**: Relativistic ZORA-def2-SVP family as main basis together with any `metal_basisset_rel` overrides (e.g. SARC-ZORA-TZVP when present in the CONTROL file).
- **Mixed systems**: Treated like 4d/5d metals to ensure relativistic settings for all transition metals in the structure.

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

#### Mathematical Formulation for Spin Projection
When broken-symmetry solutions are obtained, approximate spin projection corrects for spin contamination:

```
E_projected = E_BS + (E_HS - E_BS) × [S(S+1) - ⟨S²⟩_BS] / [⟨S²⟩_HS - ⟨S²⟩_BS]
```

Where:
- E_BS: Broken-symmetry energy
- E_HS: High-spin energy
- S: Target spin quantum number
- ⟨S²⟩: Expectation value of spin operator

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

### 2. Parallel Processing and Global Resource Management

DELFIN implements sophisticated multi-level parallelization with a **global job manager singleton** to ensure consistent resource allocation across all computational workflows:

#### Global Job Manager Architecture
- **Singleton pattern**: A single `GlobalJobManager` instance coordinates all workflows throughout execution
- **Centralized PAL management**: CPU cores (`PAL`) are read once from `CONTROL.txt` at startup and managed globally
- **Thread-safe coordination**: All parallel workflows coordinate through a shared `DynamicCorePool` with proper locking
- **Subprocess coordination**: OCCUPIER subprocesses receive their allocated PAL via environment variables (`DELFIN_CHILD_GLOBAL_MANAGER`)
- **No double allocation**: When oxidation and reduction workflows run in parallel, cores are automatically split (e.g., PAL=12 → 6 cores per workflow)

#### Dynamic Core Pool Management
- **Intelligent resource allocation**: Dynamic distribution of CPU cores across concurrent jobs through the global pool
- **Memory-aware scheduling**: Jobs allocated based on both core and memory requirements
- **Priority-based execution**: High-priority jobs receive preferential resource allocation
- **Automatic PAL splitting**: When `parallel_workflows = yes` and both ox/red workflows are configured, PAL is divided between workflows
- **Sequential mode**: When `parallel_workflows = no`, each workflow uses the full PAL sequentially

#### Workflow-Level Parallelization
- **Oxidation/Reduction workflows**: Independent execution of redox workflows with automatic core splitting
- **OCCUPIER FoB parallelization**: Multiple Fragments of Bulk (FoB) calculations executed simultaneously within allocated resources
- **Dependency-aware scheduling**: Jobs with dependencies scheduled appropriately to maximize throughput
- **Thread-safe file operations**: `thread_safe_helpers.py` ensures proper handling of CONTROL.txt updates and geometry files

#### ORCA Integration
- **PAL process management**: Configurable number of CPU cores per ORCA calculation (`PAL` parameter)
- **Consistent PAL enforcement**: All ORCA jobs respect the globally managed PAL limits
- **Memory optimization**: Per-job memory allocation (`maxcore` specification)
- **Cluster environment detection**: Automatic resource detection on SLURM/PBS/LSF systems via `cluster_utils.py`

#### Implementation Details
- **Main process**: Initializes global manager with full PAL from CONTROL.txt
- **Parallel execution**: Calculates `reduced_pal = max(2, total_pal // 2)` when both workflows run
- **CONTROL.txt propagation**: Each subprocess folder gets updated CONTROL.txt with correct PAL value
- **Environment bootstrap**: Subprocesses initialize their global manager from `DELFIN_CHILD_GLOBAL_MANAGER` JSON payload

This architecture ensures that the total CPU allocation (PAL) is never exceeded, regardless of how many workflows or subprocesses are running simultaneously.

### 3. Reproducibility

- Complete parameter logging in output files
- Systematic file naming conventions
- Configuration file preservation

## Performance and Validation

### Benchmark Performance
DELFIN has been extensively validated across multiple datasets totaling over 600 experimentally characterized systems:

#### Ground-State Redox Potentials
- **Mean Absolute Deviation (MAD)**: 0.08-0.17 eV depending on system class
- **Adjusted R²**: 0.93-0.97 for most datasets
- **Coverage**: Organic molecules, transition metal complexes, photoredox catalysts

#### Excited-State Properties
- **E₀₀ energies**: MAD = 0.10 eV, R² = 0.76
- **Excited-state redox potentials**: MAD = 0.12-0.16 eV

#### Computational Efficiency
- **Global optimization**: GOAT/CREST integration for robust geometry prediction
- **Parallel execution**: Multi-level parallelization for oxidation/reduction workflows
- **Automated workflows**: Minimal user intervention required

### Accuracy Comparison
DELFIN demonstrates competitive or superior accuracy compared to manual DFT calculations while providing full automation for non-experts in quantum chemistry.

## Computational Parameters and Convergence Criteria

### OCCUPIER Selection Parameters
- **Tolerance threshold**: `occupier_epsilon = 5e-4` (default energy tolerance)
- **Maximum iterations**: `maxiter_occupier = 100`
- **Numerical precision**: `occupier_precision = 3-6` (configurable)
- **Selection method**: `tolerance | truncation | rounding`

### State Selection Criteria
- **Clean override window**: `0.003 Hartree` (energy window for clean state preference)
- **Quality improvement threshold**: `0.05 Hartree` (minimum improvement for state ranking)
- **Bias window**: `0.003 Hartree` (energy bias tolerance)
- **Spin contamination tolerance**: System-dependent ⟨S²⟩ analysis

### ORCA Convergence Settings
- **SCF convergence**: Default ORCA criteria with configurable `maxiter`
- **Geometry optimization**: Standard BFGS with automatic coordinate propagation
- **Frequency calculations**: Analytical Hessian when available

### Grid and Numerical Settings
- **DFT integration grid**: ORCA defaults (Grid4 or Grid5 for heavy atoms)
- **RI approximation**: RIJCOSX with def2/J or SARC/J auxiliary sets
- **Memory allocation**: `maxcore` parameter (auto-detected or user-specified)

### 9. Broken-Symmetry DFT and Electronic Structure Analysis

#### Theoretical Background
DELFIN employs broken-symmetry (BS) DFT to describe systems with antiferromagnetic coupling, particularly relevant for transition metal complexes with non-innocent ligands. This approach addresses the challenge of modeling systems where unpaired electrons are coupled between metal centers and ligand π* orbitals.

#### Spin Contamination and State Selection
- **Spin contamination analysis**: Uses ⟨S²⟩ deviations from ideal values to identify appropriate electronic descriptions
- **State ranking**: Based on electronic energies (or Gibbs free energies) with spin-contamination diagnostics as tie-breakers
- **Approximate spin projection**: When BS-DFT is used, projected energies correct for spin contamination effects

#### Exchange-Correlation Functional Selection
Recent benchmarks identify optimal functionals for spin-state energetics:
- **PBE0**: Robust hybrid functional for transition metal spin states
- **TPSSh**: Alternative hybrid showing good performance for spin-crossover systems
- **r²SCAN**: Semilocal meta-GGA with excellent spin-state accuracy

### 10. Time-Dependent DFT for Excited States

For closed-shell species, DELFIN computes S₁/T₁ optimizations and E₀₀ transition energies:

#### S₁ State Optimization
- **Method**: TD-DFT with multiple singlet roots (typically ~10)
- **State tracking**: FOLLOWROOT keyword ensures consistent state tracking during optimization
- **Vertical transitions**: Absorption/emission spectra calculated from optimized geometries

#### Spin-Orbit Coupling
- **Optional SOC corrections**: Account for relativistic effects in excited states
- **Implementation**: Perturbative treatment of spin-orbit coupling for improved accuracy

## Error Analysis and Quality Control

### Systematic Error Sources
- **Functional dependence**: Choice of exchange-correlation functional affects spin-state energetics
- **Basis set incompleteness**: Finite basis set effects on energy differences
- **Solvation model limitations**: Continuum models vs. explicit solvation effects
- **Conformational sampling**: Global optimization may not capture all relevant minima

### Statistical Error Analysis
- **Energy precision**: Typically 0.001-0.01 Hartree depending on system size and convergence
- **Spin contamination**: ⟨S²⟩ deviations monitored for broken-symmetry states
- **Geometry optimization**: RMS gradient convergence to 10⁻⁴ Hartree/Bohr

### Quality Metrics
- **State assignment confidence**: Based on energy gaps and spin contamination analysis
- **Wavefunction stability**: Automatic detection of wavefunction instabilities
- **Convergence monitoring**: Systematic tracking of SCF and geometry convergence

### Error Propagation
- **Redox potential uncertainty**: Propagated from individual state energies (typically ±0.05-0.15 eV)
- **Temperature effects**: Entropy contributions affect Gibbs free energy differences
- **Reference electrode**: Solvent-specific E_ref uncertainties

## Scope and Limitations

### System Requirements
- **Software dependencies**: ORCA 6.1.0+, optional xTB 6.7.1+ and CREST 3.0.2
- **Hardware requirements**: Minimum 4 GB RAM, multicore CPU recommended
- **Operating systems**: Linux, macOS, Windows (with appropriate ORCA installation)

### Methodological Scope
- **Redox steps**: Up to 3 sequential oxidation/reduction steps
- **Spin multiplicities**: Automatic detection from singlet to high-spin states
- **System types**: Closed-shell and open-shell molecules, transition metal complexes
- **Solvents**: Implicit solvation (CPCM/SMD), optional explicit solvation

### Current Limitations
- **System size**: Practical limit ~500 atoms for routine calculations
- **Heavy atoms**: Requires appropriate relativistic treatment (4d/5d metals)
- **Strongly correlated systems**: Single-reference DFT may be insufficient
- **Dynamic effects**: Static calculations do not capture dynamic solvation

### Accuracy Expectations
- **Ground-state redox potentials**: MAD = 0.08-0.17 eV (system dependent)
- **Excited-state properties**: MAD = 0.10-0.16 eV (higher uncertainty)
- **Spin-state assignment**: >95% success rate for well-defined systems
- **Geometry prediction**: RMSD < 1.0 Å for most transition metal complexes

### Best Practice Recommendations
- **Input preparation**: Use reasonable starting geometries, validate with XRD when available
- **Functional selection**: PBE0-D4 recommended for most systems, test alternatives for spin-crossover
- **Basis sets**: Default selection adequate for most cases, consider larger sets for accuracy-critical applications
- **Convergence**: Monitor SCF and geometry convergence, increase iterations if needed

## References

The methodology implemented in DELFIN builds upon established computational chemistry protocols and leverages the capabilities of:

### Core Software Components
- **ORCA**: Neese, F. et al. *J. Chem. Phys.* **152**, 224108 (2020). doi:10.1063/5.0004608
- **xTB methods**: Bannwarth, C. et al. *WIREs Comput. Mol. Sci.* **11**, e1493 (2021). doi:10.1002/wcms.1493
- **CREST**: Pracht, P. et al. *J. Chem. Phys.* **160**, 114110 (2024). doi:10.1063/5.0197592

### Theoretical Foundations
- **DFT Framework**: Kohn, W. & Sham, L. J. *Phys. Rev.* **140**, A1133 (1965). doi:10.1103/PhysRev.140.A1133
- **Broken-Symmetry DFT**: Noodleman, L. *J. Chem. Phys.* **74**, 5737 (1981). doi:10.1063/1.440939
- **Redox Potential Calculations**: Marenich, A. V. et al. *Phys. Chem. Chem. Phys.* **16**, 15068 (2014). doi:10.1039/c4cp01572j
- **Spin-State Energetics**: Radoń, M. et al. *Chem. Sci.* **15**, 20189 (2024). doi:10.1039/d4sc04903d
- **PBE0 Functional**: Adamo, C. & Barone, V. *J. Chem. Phys.* **110**, 6158 (1999). doi:10.1063/1.478522
- **ZORA Relativistic Method**: van Lenthe, E. et al. *J. Chem. Phys.* **99**, 4597 (1993). doi:10.1063/1.466059
- **TD-DFT**: Petrenko, T. et al. *J. Chem. Phys.* **134**, 054116 (2011). doi:10.1063/1.3533441

### Solvation Models
- **CPCM**: Garcia-Ratés, M. & Neese, F. *J. Comput. Chem.* **41**, 922 (2020). doi:10.1002/jcc.26139
- **SMD**: Marenich, A. V. et al. *J. Phys. Chem. B* **113**, 6378 (2009). doi:10.1021/jp810292n

The specific implementation details for transition metal basis set policies, relativistic corrections, automated spin-state assignment, and parallel processing represent novel automation approaches developed for this software package.
