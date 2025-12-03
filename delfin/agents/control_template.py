"""
DELFIN CONTROL file template with all available options
"""

CONTROL_TEMPLATE = """input_file=input.txt
NAME={name}
SMILES={smiles}
charge={charge}
------------------------------------
Solvation:
implicit_solvation_model={solvation_model}
solvent={solvent}
XTB_SOLVATOR={xtb_solvator}
number_explicit_solv_molecules={explicit_solv_molecules}
------------------------------------
Global geometry optimisation:
xTB_method={xtb_method}
XTB_OPT={xtb_opt}
XTB_GOAT={xtb_goat}
CREST={crest}
multiplicity_global_opt={multiplicity_global_opt}
------------------------------------
IMAG={imag}
IMAG_scope={imag_scope}
IMAG_option={imag_option}
allow_imaginary_freq={allow_imaginary_freq}
IMAG_sp_energy_window={imag_sp_energy_window}
IMAG_optimize_candidates={imag_optimize_candidates}
------------------------------------
Redox steps:
calc_initial={calc_initial}
oxidation_steps={oxidation_steps}
reduction_steps={reduction_steps}
method={method}
calc_potential_method={calc_potential_method}
------------------------------------
E_00={e_00}
excitation={excitation}
S1_opt={s1_opt}
triplet_flag={triplet_flag}
absorption_spec={absorption_spec}
emission_spec={emission_spec}
NROOTS={nroots}
TDA={tda}
NACME={nacme}
ETF={etf}
DONTO={donto}
DOSOC={dosoc}
singlet exitation:
IROOT={iroot}
FOLLOWIROOT={followiroot}
mcore_E00={mcore_e00}
------------------------------------
MANUALLY:
{manually_section}
------------------------------------
Level of Theory:
functional={functional}
disp_corr={disp_corr}
ri_jkx={ri_jkx}
ri_soc={ri_soc}
relativity={relativity}
aux_jk={aux_jk}
aux_jk_rel={aux_jk_rel}
main_basisset={main_basisset}
main_basisset_rel={main_basisset_rel}
metal_basisset={metal_basisset}
metal_basisset_rel={metal_basisset_rel}
first_coordination_sphere_metal_basisset={first_coordination_sphere_metal_basisset}
first_coordination_sphere_scale={first_coordination_sphere_scale}
geom_opt={geom_opt}
freq_type={freq_type}
initial_guess={initial_guess}
temperature={temperature}
maxiter={maxiter}
qmmm_option={qmmm_option}
----------------
deltaSCF:
deltaSCF_DOMOM={deltascf_domom}
deltaSCF_PMOM={deltascf_pmom}
deltaSCF_keepinitialref={deltascf_keepinitialref}
deltaSCF_SOSCFHESSUP={deltascf_soscfhessup}
----------------
ESD_modul={esd_modul}
ESD_modus={esd_modus}
ESD_TDDFT_maxiter={esd_tddft_maxiter}
ESD_nroots={esd_nroots}
ESD_maxdim={esd_maxdim}
states={esd_states}
ISCs={esd_iscs}
ICs={esd_ics}
------------------------------------
Reference value:
E_ref={e_ref}
------------------------------------
Literature_reference={literature_reference}
reference_CV={reference_cv}
E_00_exp={e_00_exp}
E_red_exp={e_red_exp}
E_red_2_exp={e_red_2_exp}
E_red_3_exp={e_red_3_exp}
E_ox_exp={e_ox_exp}
E_ox_2_exp={e_ox_2_exp}
E_ox_3_exp={e_ox_3_exp}
*E_red_exp={excited_e_red_exp}
*E_ox_exp={excited_e_ox_exp}
------------------------------------
Prints:
print_MOs={print_mos}
print_Loewdin_population_analysis={print_loewdin}
------------------------------------
Resource Settings:
PAL={pal}
maxcore={maxcore}
parallel_workflows={parallel_workflows}
pal_jobs={pal_jobs}
enable_job_timeouts={enable_job_timeouts}
job_timeout_hours={job_timeout_hours}
opt_timeout_hours={opt_timeout_hours}
frequency_timeout_hours={frequency_timeout_hours}
sp_timeout_hours={sp_timeout_hours}
------------------------------------
Automatic Error Recovery & Retry:
enable_auto_recovery={enable_auto_recovery}
max_recovery_attempts={max_recovery_attempts}
enable_adaptive_parallelism={enable_adaptive_parallelism}
enable_performance_metrics={enable_performance_metrics}
------------------------------------
OCCUPIER-Settings:
--------------------
OCCUPIER_method={occupier_method}
OCCUPIER_tree={occupier_tree}
OWN_TREE_PURE_WINDOW={own_tree_pure_window}
OWN_progressive_from={own_progressive_from}
frequency_calculation_OCCUPIER={frequency_calculation_occupier}
occupier_selection={occupier_selection}
occupier_precision={occupier_precision}
occupier_epsilon={occupier_epsilon}
maxiter_occupier={maxiter_occupier}
geom_opt_OCCUPIER={geom_opt_occupier}
pass_wavefunction={pass_wavefunction}
approximate_spin_projection_APMethod={approximate_spin_projection_apmethod}
--------------------
OCCUPIER_sequence_profiles:
-3,-2,-1,0,+1,+2,+3={occupier_sequence_profiles}
"""

# Default values for all parameters
DEFAULT_VALUES = {
    "name": "",
    "smiles": "",
    "charge": "0",
    "solvation_model": "CPCM",
    "solvent": "Water",
    "xtb_solvator": "no",
    "explicit_solv_molecules": "2",
    "xtb_method": "XTB2",
    "xtb_opt": "no",
    "xtb_goat": "no",
    "crest": "no",
    "multiplicity_global_opt": "",
    "imag": "yes",
    "imag_scope": "initial",
    "imag_option": "2",
    "allow_imaginary_freq": "0",
    "imag_sp_energy_window": "1e-3",
    "imag_optimize_candidates": "no",
    "calc_initial": "yes",
    "oxidation_steps": "1,2,3",
    "reduction_steps": "1,2,3",
    "method": "classic",
    "calc_potential_method": "2",
    "e_00": "no",
    "excitation": "s",
    "s1_opt": "TDDFT",
    "triplet_flag": "FALSE",
    "absorption_spec": "no",
    "emission_spec": "no",
    "nroots": "15",
    "tda": "FALSE",
    "nacme": "TRUE",
    "etf": "TRUE",
    "donto": "FALSE",
    "dosoc": "TRUE",
    "iroot": "1",
    "followiroot": "TRUE",
    "mcore_e00": "10000",
    "manually_section": "",
    "functional": "PBE0",
    "disp_corr": "D4",
    "ri_jkx": "RIJCOSX",
    "ri_soc": "RI-SOMF(1X)",
    "relativity": "ZORA",
    "aux_jk": "def2/J",
    "aux_jk_rel": "SARC/J",
    "main_basisset": "def2-SVP",
    "main_basisset_rel": "ZORA-def2-SVP",
    "metal_basisset": "def2-TZVP",
    "metal_basisset_rel": "SARC-ZORA-TZVP",
    "first_coordination_sphere_metal_basisset": "no",
    "first_coordination_sphere_scale": "1.3",
    "geom_opt": "OPT",
    "freq_type": "FREQ",
    "initial_guess": "PModel",
    "temperature": "298.15",
    "maxiter": "125",
    "qmmm_option": "QM/PBEH-3c",
    "deltascf_domom": "true",
    "deltascf_pmom": "true",
    "deltascf_keepinitialref": "true",
    "deltascf_soscfhessup": "LBFGS",
    "esd_modul": "no",
    "esd_modus": "TDDFT",
    "esd_tddft_maxiter": "1000",
    "esd_nroots": "15",
    "esd_maxdim": "30",
    "esd_states": "S0,S1,T1,T2",
    "esd_iscs": "S1>T1,T1>S1,S1>T2,T2>S1",
    "esd_ics": "S1>S0,T2>T1",
    "e_ref": "",
    "literature_reference": "",
    "reference_cv": "V Vs. Fc+/Fc",
    "e_00_exp": "",
    "e_red_exp": "",
    "e_red_2_exp": "",
    "e_red_3_exp": "",
    "e_ox_exp": "",
    "e_ox_2_exp": "",
    "e_ox_3_exp": "",
    "excited_e_red_exp": "",
    "excited_e_ox_exp": "",
    "print_mos": "no",
    "print_loewdin": "no",
    "pal": "12",
    "maxcore": "6000",
    "parallel_workflows": "yes",
    "pal_jobs": "4",
    "enable_job_timeouts": "yes",
    "job_timeout_hours": "36",
    "opt_timeout_hours": "14",
    "frequency_timeout_hours": "36",
    "sp_timeout_hours": "3",
    "enable_auto_recovery": "yes",
    "max_recovery_attempts": "1",
    "enable_adaptive_parallelism": "yes",
    "enable_performance_metrics": "yes",
    "occupier_method": "auto",
    "occupier_tree": "deep",
    "own_tree_pure_window": "3",
    "own_progressive_from": "no",
    "frequency_calculation_occupier": "no",
    "occupier_selection": "tolerance",
    "occupier_precision": "3",
    "occupier_epsilon": "5e-4",
    "maxiter_occupier": "125",
    "geom_opt_occupier": "OPT",
    "pass_wavefunction": "no",
    "approximate_spin_projection_apmethod": "2",
    "occupier_sequence_profiles": "[\neven electron number:\neven_seq = [\n  {\"index\": 1, \"m\": 1, \"BS\": \"\",    \"from\": 0},\n  {\"index\": 2, \"m\": 1, \"BS\": \"1,1\", \"from\": 1},\n  {\"index\": 3, \"m\": 1, \"BS\": \"2,2\", \"from\": 2},\n  {\"index\": 4, \"m\": 3, \"BS\": \"\",    \"from\": 1},\n  {\"index\": 5, \"m\": 3, \"BS\": \"3,1\", \"from\": 4},\n  {\"index\": 6, \"m\": 3, \"BS\": \"4,2\", \"from\": 5},\n  {\"index\": 7, \"m\": 5, \"BS\": \"\",    \"from\": 4},\n  {\"index\": 8, \"m\": 5, \"BS\": \"5,1\", \"from\": 7},\n  {\"index\": 9, \"m\": 5, \"BS\": \"6,2\", \"from\": 8}\n]\n-------------------\nodd electron number:\nodd_seq = [\n  {\"index\": 1, \"m\": 2, \"BS\": \"\",    \"from\": 0},\n  {\"index\": 2, \"m\": 2, \"BS\": \"2,1\", \"from\": 1},\n  {\"index\": 3, \"m\": 2, \"BS\": \"3,2\", \"from\": 2},\n  {\"index\": 4, \"m\": 4, \"BS\": \"\",    \"from\": 1},\n  {\"index\": 5, \"m\": 4, \"BS\": \"4,1\", \"from\": 4},\n  {\"index\": 6, \"m\": 4, \"BS\": \"5,2\", \"from\": 5},\n  {\"index\": 7, \"m\": 6, \"BS\": \"\",    \"from\": 4},\n  {\"index\": 8, \"m\": 6, \"BS\": \"6,1\", \"from\": 7},\n  {\"index\": 9, \"m\": 6, \"BS\": \"7,2\", \"from\": 8}\n]\n]"
}
