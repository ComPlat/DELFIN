# delfin/define.py
# -*- coding: utf-8 -*-
from delfin.common.logging import get_logger
from delfin.common.paths import resolve_path

logger = get_logger(__name__)

TEMPLATE = """input_file=input.txt
NAME=
SMILES=
charge=[CHARGE]
------------------------------------
Solvation:
implicit_solvation_model=CPCM
solvent=[SOLVENT]
XTB_SOLVATOR=no
number_explicit_solv_molecules=2
------------------------------------
Global geometry optimisation:
xTB_method=XTB2
XTB_OPT=no
XTB_GOAT=no
CREST=no
multiplicity_global_opt=
------------------------------------
IMAG=yes
IMAG_scope=initial
IMAG_option=2
allow_imaginary_freq=0
IMAG_sp_energy_window=1e-3
IMAG_optimize_candidates=no
------------------------------------
calc_prop_of_interest=no
properties_of_interest=IP,EA
------------------------------------
Redox steps:
calc_initial=yes
oxidation_steps=1,2,3
reduction_steps=1,2,3
method=classic|manually|OCCUPIER
calc_potential_method=2
------------------------------------
ESD module (excited state dynamics):
ESD_modul=no
ESD_modus=TDDFT|deltaSCF|hybrid1
ESD_frequency=yes
states=S1,T1,S2,T2
ISCs=S1>T1,T1>S1
ICs=S2>S1
emission_rates=f,p
phosp_IROOT=1,2,3
phosp_keywords=
fluor_keywords=
TROOTSSL=-1,0,1
addition_S0=
DOHT=TRUE
ESD_LINES=LORENTZ
ESD_LINEW=50
ESD_INLINEW=250
ESD_NPOINTS=131072
ESD_MAXTIME=12000
hybrid1_geom_MaxIter=60
--------------------
Electrical Properties:
elprop_Dipole=no
elprop_Quadrupole=no
elprop_Hyperpol=no
elprop_Polar=no
elprop_PolarVelocity=no
elprop_PolarDipQuad=no
elprop_PolarQuadQuad=no
--------------------
deltaSCF Settings:
deltaSCF_DOMOM=true
deltaSCF_PMOM=true
deltaSCF_keepinitialref=true
deltaSCF_SOSCFHESSUP=LSR1
deltaSCF_keywords=FreezeAndRelease
deltaSCF_maxiter=300
deltaSCF_SOSCFConvFactor=500
deltaSCF_SOSCFMaxStep=0.1
--------------------
TDDFT Settings:
TDDFT_TDDFT_maxiter=500
TDDFT_nroots=15
TDDFT_maxdim=30
TDDFT_TDA=FALSE
TDDFT_followiroot=true
TDDFT_SOC=false
------------------------------------
MANUALLY:
multiplicity_0=
additions_0=
additions_TDDFT=
additions_T1=
additions_S1=
multiplicity_ox1=
additions_ox1=
multiplicity_ox2=
additions_ox2=
multiplicity_ox3=
additions_ox3=
multiplicity_red1=
additions_red1=
multiplicity_red2=
additions_red2=
multiplicity_red3=
additions_red3=
------------------------------------
Level of Theory:
functional=PBE0
disp_corr=D4
ri_jkx=RIJCOSX
relativity=ZORA
aux_jk=def2/J
aux_jk_rel=SARC/J
main_basisset=def2-SVP
main_basisset_rel=ZORA-def2-SVP
metal_basisset=def2-TZVP
metal_basisset_rel=SARC-ZORA-TZVP
first_coordination_sphere_metal_basisset=no
first_coordination_sphere_scale=1.3
geom_opt=OPT
freq_type=FREQ
initial_guess=PModel
temperature=298.15
maxiter=125
qmmm_option=QM/PBEH-3c
------------------------------------
Reference value:
E_ref=
------------------------------------
Literature_reference=
reference_CV=V Vs. Fc+/Fc
E_00_exp=
E_red_exp=
E_red_2_exp=
E_red_3_exp=
E_ox_exp=
E_ox_2_exp=
E_ox_3_exp=
*E_red_exp=
*E_ox_exp=
------------------------------------
Prints:
print_MOs=no
print_Loewdin_population_analysis=no
------------------------------------
Resource Settings:
PAL=12
maxcore=6000
parallel_workflows=yes
pal_jobs=3
enable_job_timeouts=no
job_timeout_hours=36
opt_timeout_hours=14
frequency_timeout_hours=36
sp_timeout_hours=3
------------------------------------
Automatic Error Recovery & Retry:
enable_auto_recovery=yes
max_recovery_attempts=1
enable_adaptive_parallelism=yes
enable_performance_metrics=yes
------------------------------------
OCCUPIER-Settings:
--------------------
OCCUPIER_method=auto|manually
OCCUPIER_tree=own
OWN_TREE_PURE_WINDOW=3
OWN_progressive_from=no
fob_equal_weights=yes
frequency_calculation_OCCUPIER=no
occupier_selection=tolerance|truncation|rounding
occupier_precision=3
occupier_epsilon=5e-4
clean_override_window_h=0.002
clean_quality_improvement=0.05
clean_quality_good=0.05
maxiter_occupier=125
geom_opt_OCCUPIER=OPT
pass_wavefunction=no
approximate_spin_projection_APMethod=2
--------------------
OCCUPIER_sequence_profiles:
-3,-2,-1,0,+1,+2,+3=[
even electron number:
even_seq = [
  {"index": 1, "m": 1, "BS": "",    "from": 0},
  {"index": 2, "m": 1, "BS": "1,1", "from": 1},
  {"index": 3, "m": 1, "BS": "2,2", "from": 2},
  {"index": 4, "m": 3, "BS": "",    "from": 1},
  {"index": 5, "m": 3, "BS": "3,1", "from": 4},
  {"index": 6, "m": 3, "BS": "4,2", "from": 5},
  {"index": 7, "m": 5, "BS": "",    "from": 4}
]
-------------------
odd electron number:
odd_seq = [
  {"index": 1, "m": 2, "BS": "",    "from": 0},
  {"index": 2, "m": 2, "BS": "2,1", "from": 1},
  {"index": 3, "m": 2, "BS": "3,2", "from": 2},
  {"index": 4, "m": 4, "BS": "",    "from": 1},
  {"index": 5, "m": 4, "BS": "4,1", "from": 4},
  {"index": 6, "m": 4, "BS": "5,2", "from": 5},
  {"index": 7, "m": 6, "BS": "",    "from": 4}
]
]

INFOS:
-------------------------------------------------
Available METHODS: classic, manually, OCCUPIER
Available OX_STEPS: 1 ; 1,2 ; 1,2,3 ; 2 ; 3 ; 2,3 ; 1,3
Available RED_STEPS: 1 ; 1,2 ; 1,2,3 ; 2 ; 3 ; 2,3 ; 1,3
Available IMPLICIT SOLVATION MODELS: CPCM ; CPCMC ; SMD
Available dispersion corrections DISP_CORR: D4 ; D3 ; D3BJ ; D3ZERO ; NONE
ESD MODULE: Set ESD_modul=yes and specify states=[S0,S1,T1] for excited state calculations
Available states: S0, S1, S2, T1, T2, T3
ISCs: S1>T1, T1>S0, etc. (intersystem crossing rates)
ICs: S1>S0, S2>S1, etc. (internal conversion rates)
Available qmmm_option: QM/XTB ; QM/PBEH-3C ; QM/HF-3C ; QM/r2SCAN-3C (for QM/MM calculations)
Available freq_type: FREQ (analytic, default) ; NUMFREQ (numerical, required for WB97X-V and other DFT-NL functionals)
EXPLICIT SOLVATION MODEL IS VERY EXPENSIVE!!!!!
IMAG_option:
  1 -> red/ox OCCUPIER continues immediately (IMAG and OCCUPIER run in parallel)
  2 -> red/ox OCCUPIER waits for IMAG to finish and uses the refined geometry
OCCUPIER_option:
# flat   -> legacy flat sequences with BS
# deep2  -> only pure states, no BS (simple testing)
# deep   -> adaptive BS evolution (reduction: BS(m-1,1) or BS(M±1,N); oxidation: pure only)
# own    -> rule-based mode using custom/auto sequences
OWN_TREE_PURE_WINDOW: Number of pure states to test around previous winner (default: 1)
OWN_progressive_from: Sequential pure states (yes) or all parallel (no, default)
  - no  -> All pure states start from "from": 0 (maximum parallelization)
  - yes -> Pure states build sequentially: m=1 -> m=3 -> m=5 (safer but slower)
-------------------------------------------------
OCCUPIER Selection Parameters (Fine-tuning):
clean_override_window_h: Energy window (Hartree) for considering cleaner candidates (default: 0.002)
  → If a candidate has much lower spin contamination and energy within this window, it may be preferred
  → Reduce (e.g., 0.001) to strictly prefer lowest energy regardless of spin contamination
  → Increase (e.g., 0.004) to allow cleaner solutions even with slightly higher energy
clean_quality_improvement: Minimum quality improvement to trigger clean-override (default: 0.05)
  → Quality = spin contamination deviation for non-BS, or BS pair mismatch for BS states
  → Increase (e.g., 0.10) to require larger quality differences before preferring cleaner states
clean_quality_good: Absolute quality threshold for clean-override (default: 0.05)
  → If quality < this value, candidate is considered very clean and may be preferred
-------------------------------------------------
ESD (Excited State Dynamics):
ESD_MODUL: yes/no - Enable ESD calculations in separate ESD/ directory
states: Comma-separated list of states to calculate (S0, S1, T1, T2)
ISCs: Comma-separated list of intersystem crossings (e.g., S1>T1, T1>S1)
ICs: Comma-separated list of internal conversions (e.g., S1>S0, T1>T2)
ESD_modus: TDDFT (default), deltaSCF, or hybrid1 - Method for excited state calculations
  → TDDFT: Standard TDDFT-based excited state optimization
  → deltaSCF: Direct deltaSCF optimization (may collapse to S0 for some systems)
  → hybrid1: Two-step approach - TDDFT OPT (no FREQ) → deltaSCF OPT FREQ
             Creates {state}_first_TDDFT.inp and {state}_second.inp
             Prevents collapse to S0 by using TDDFT pre-optimization
ESD_frequency: yes/no - Enable frequency calculations (default: yes)
  → When 'no': Only geometry optimization (OPT) is performed
  → ISC/IC calculations are disabled (require ZPE from frequency calculations)
  → Only Final single point energies are available
  → Useful for quick tests or when frequencies are not needed
ESD_followiroot: true/false - Follow TDDFT root during geometry optimization (TDDFT mode only, default: true)

-------------------------------------------------
Automatic Error Recovery & Retry System:

⭐ Intelligent Error Recovery (Merged Retry + Recovery):
enable_auto_recovery: yes/no - Enable intelligent error detection and recovery (default: yes)
  NOTE: Old parameter 'orca_retry_enabled' still supported for backward compatibility
max_recovery_attempts: Number of retries after initial failure (default: 1)
  NOTE: Old parameter 'orca_retry_max_attempts' still supported for backward compatibility
  → 1 = try once more (2 runs total), 2 = try twice more (3 runs total), etc.
  → Handles BOTH ORCA-specific errors AND transient system errors
  → Uses MOREAD to continue from last .gbw (no restart from scratch!)
  → Automatically updates geometry from latest .xyz file
  → Creates GBW backup to prevent ORCA deletion on failure

⏱️ Job Lifetime Management:
enable_job_timeouts: yes/no - Enable job timeout limits (default: yes)
  → When yes: Jobs are terminated after timeout (job_timeout_hours, opt_timeout_hours, etc.)
  → When no: Jobs run indefinitely until completion or error
  → Disable for difficult systems that need >24h runtime
  → Individual timeout settings still configurable via job_timeout_hours, opt_timeout_hours, etc.

📊 Performance Monitoring:
enable_adaptive_parallelism: yes/no - Dynamically adjust core count based on job requirements
enable_performance_metrics: yes/no - Track and log job performance statistics

How auto_recovery works:
1. ORCA runs with your input file
2. If it fails, DELFIN detects the specific error type:
   a) ORCA-specific errors (SCF, TRAH, geometry, etc.)
   b) Transient system errors (disk full, network timeout, etc.)
3. Automatically modifies input file with appropriate fixes (MOREAD, NoAutoTRAH, etc.)
4. For transient errors: Applies exponential backoff (2s, 4s, 8s, ...) before retry
5. Continues from last .gbw state (no restart from scratch!)
6. Progressively escalates fixes if problem persists
7. Tracks state to prevent infinite loops

Supported error types:
- SCF convergence failures → SlowConv/VerySlowConv, KDIIS, SOSCF, high damping
- TRAH segmentation faults → NoAutoTRAH
- DIIS errors → KDIIS, SOSCF (second-order SCF)
- Geometry convergence → Fixed trust radius (negative), Hessian recalculation
- MPI crashes → Reduced parallelism, different MPI transport
- Memory errors → Reduced maxcore/PAL
- Frequency failures → NumFreq or skip
- Transient system errors → Exponential backoff retry (disk full, network, I/O)

See docs/RECOVERY_STRATEGIES_DETAILED.txt for all settings and ORCA Manual references
-------------------------------------------------
"""
# -------------------------------------------------------------------------------------------------------
def convert_xyz_to_input_txt(src_xyz: str, dst_txt: str = "input.txt") -> str:
    """Convert an XYZ file to input.txt by dropping the first two lines."""
    src_path = resolve_path(src_xyz)
    dst_path = resolve_path(dst_txt)

    if not src_path.exists():
        message = f"XYZ source '{src_xyz}' not found. Creating empty {dst_txt} instead."
        print(message)
        logger.warning(message)
        dst_path.touch(exist_ok=True)
        return dst_txt

    lines = src_path.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)
    content = "".join(lines[2:]) if len(lines) >= 2 else ""
    if content and not content.endswith("\n"):
        content += "\n"

    dst_path.write_text(content, encoding="utf-8")
    message = f"Converted '{src_xyz}' → '{dst_txt}' (dropped first two lines)."
    print(message)
    logger.info(message)
    return dst_txt
# -------------------------------------------------------------------------------------------------------
def create_control_file(filename: str = "CONTROL.txt",
                        input_file: str = "input.txt",
                        overwrite: bool = False) -> None:
    """
    Create a CONTROL.txt and create an input file.
    If input_file ends with '.xyz', convert it to 'input.txt' by dropping the first two lines.
    """
    # If user passed an .xyz, convert to input.txt and use that in CONTROL.txt
    target_input = input_file
    if str(input_file).lower().endswith(".xyz"):
        target_input = convert_xyz_to_input_txt(input_file, "input.txt")
    else:
        # Ensure empty input file exists
        target_path = resolve_path(target_input)
        if not target_path.exists():
            target_path.touch()
            message = f"{target_input} has been created (empty)."
            print(message)
            logger.info(message)

    control_path = resolve_path(filename)

    if control_path.exists() and not overwrite:
        message = f"{filename} already exists. Use --overwrite to replace it."
        print(message)
        logger.warning(message)
        return

    content = TEMPLATE.replace("{INPUT_FILE}", target_input)
    control_path.write_text(content, encoding="utf-8")
    message = f"{filename} has been written (input_file={target_input})."
    print(message)
    logger.info(message)
# -------------------------------------------------------------------------------------------------------
