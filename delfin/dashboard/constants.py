"""Constants, templates, and style definitions for the DELFIN Dashboard."""

import ipywidgets as widgets

# ---------------------------------------------------------------------------
# Job time limits
# ---------------------------------------------------------------------------
JOB_TIME_LIMITS = {
    '12h': '12:00:00',
    '24h': '24:00:00',
    '36h': '36:00:00',
    '48h': '48:00:00',
    '60h': '60:00:00',
}

# ---------------------------------------------------------------------------
# Status colours (used in Job Status tab)
# ---------------------------------------------------------------------------
STATUS_COLORS = {
    'PENDING':   '#9E9E9E',
    'RUNNING':   '#4CAF50',
    'COMPLETED': '#2196F3',
    'FAILED':    '#f44336',
    'CANCELLED': '#FF9800',
    'TIMEOUT':   '#b71c1c',
}

# ---------------------------------------------------------------------------
# Common widget styling
# ---------------------------------------------------------------------------
COMMON_LAYOUT = widgets.Layout(width='500px')
COMMON_STYLE = {'description_width': 'initial'}

# ---------------------------------------------------------------------------
# Search suggestions for Calculations Browser
# ---------------------------------------------------------------------------
CALC_SEARCH_OPTIONS = [
    'ABSORPTION SPECTRUM',
    'ERROR',
    'FINAL SINGLE POINT ENERGY',
    'Final Gibbs free energy',
    'HURRAY',
    'JOB NUMBER ',
    'LOEWDIN POPULATION ANALYSIS',
    'LOEWDIN REDUCED ORBITAL CHARGES',
    'MULLIKEN POPULATION ANALYSIS',
    'MULLIKEN REDUCED ORBITAL CHARGES',
    'ORBITAL ENERGIES',
    'SCF CONVERGED AFTER',
    'TD-DFT EXCITED STATES',
    'Total Enthalpy',
    'WARNING',
]

# ---------------------------------------------------------------------------
# Job status table CSS
# ---------------------------------------------------------------------------
JOB_TABLE_CSS = """
<style>
    .job-table { font-family: monospace; font-size: 12px; border-collapse: collapse; width: 100%; }
    .job-table th, .job-table td { padding: 6px 10px; text-align: left; border-bottom: 1px solid #ddd; }
    .job-table th { background-color: #2196F3; color: white; }
    .job-table tr:hover { background-color: #f5f5f5; }
</style>
"""

# ---------------------------------------------------------------------------
# DEFAULT CONTROL.txt template
# ---------------------------------------------------------------------------
DEFAULT_CONTROL = """input_file=input.txt
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
deltaSCF_PMOM=false
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
BrokenSym_0=
multiplicity_ox1=
BrokenSym_ox1=
multiplicity_ox2=
BrokenSym_ox2=
multiplicity_ox3=
BrokenSym_ox3=
multiplicity_red1=
BrokenSym_red1=
multiplicity_red2=
BrokenSym_red2=
multiplicity_red3=
BrokenSym_red3=
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
print_MOs=yes
print_Loewdin_population_analysis=no
------------------------------------
Resource Settings:
PAL=40
maxcore=6000
parallel_workflows=yes
pal_jobs=4
orca_parallel_strategy=auto
enable_job_timeouts=no
job_timeout_hours=36
opt_timeout_hours=14
frequency_timeout_hours=36
sp_timeout_hours=3
------------------------------------
Automatic Error Recovery & Retry:
enable_auto_recovery=yes
max_recovery_attempts=3
enable_adaptive_parallelism=yes
enable_performance_metrics=yes
------------------------------------
OCCUPIER-Settings:
--------------------
OCCUPIER_method=auto
OCCUPIER_tree=own
OWN_TREE_PURE_WINDOW=3
OWN_progressive_from=no
fob_equal_weights=yes
frequency_calculation_OCCUPIER=no
occupier_selection=tolerance
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
--------------------
ORCA base overrides (optional):
keyword:basename=[]
additions:basename=[]
"""

# ---------------------------------------------------------------------------
# Only GOAT template (XTB_GOAT=yes, no redox steps)
# ---------------------------------------------------------------------------
ONLY_GOAT_TEMPLATE = """input_file=input.txt
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
XTB_GOAT=yes
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
calc_initial=no
oxidation_steps=
reduction_steps=
method=classic
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
deltaSCF_PMOM=false
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
BrokenSym_0=
multiplicity_ox1=
BrokenSym_ox1=
multiplicity_ox2=
BrokenSym_ox2=
multiplicity_ox3=
BrokenSym_ox3=
multiplicity_red1=
BrokenSym_red1=
multiplicity_red2=
BrokenSym_red2=
multiplicity_red3=
BrokenSym_red3=
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
print_MOs=yes
print_Loewdin_population_analysis=no
------------------------------------
Resource Settings:
PAL=40
maxcore=6000
parallel_workflows=yes
pal_jobs=4
orca_parallel_strategy=auto
enable_job_timeouts=no
job_timeout_hours=36
opt_timeout_hours=14
frequency_timeout_hours=36
sp_timeout_hours=3
------------------------------------
Automatic Error Recovery & Retry:
enable_auto_recovery=yes
max_recovery_attempts=3
enable_adaptive_parallelism=yes
enable_performance_metrics=yes
------------------------------------
OCCUPIER-Settings:
--------------------
OCCUPIER_method=auto
OCCUPIER_tree=own
OWN_TREE_PURE_WINDOW=3
OWN_progressive_from=no
fob_equal_weights=yes
frequency_calculation_OCCUPIER=no
occupier_selection=tolerance
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
--------------------
ORCA base overrides (optional):
keyword:basename=[]
additions:basename=[]
"""
