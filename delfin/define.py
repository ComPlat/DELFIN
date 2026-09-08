# delfin/define.py
# -*- coding: utf-8 -*-
import re

from delfin.common.logging import get_logger
from delfin.common.paths import resolve_path

logger = get_logger(__name__)

TEMPLATE = """input_file=input.txt
NAME=
SMILES=
charge=[CHARGE]
------------------------------------------------------
Solvation:
----------
implicit_solvation_model=CPCM
solvent=[SOLVENT]
XTB_SOLVATOR=no
number_explicit_solv_molecules=2
------------------------------------------------------
SMILES conversion:
----------
smiles_converter=[QUICK|NORMAL|GUPPY|ARCHITECTOR]
------------------------------------------------------
Global geometry optimisation:
----------
xTB_method=XTB2
XTB_preOPT=no
global_optimizer=[GOAT|CREST]
multiplicity_global_opt=
------------------------------------------------------
Imaginary mode elimination:
----------
IMAG=yes
IMAG_scope=initial
IMAG_option=2
allow_imaginary_freq=0
IMAG_sp_energy_window=1e-3
IMAG_optimize_candidates=no
------------------------------------------------------
Redox steps:
------------------
method=[classic|manually|OCCUPIER]
calc_initial=yes
oxidation_steps=
reduction_steps=
calc_potential_method=2
E_ref=
---------------------------------
MANUALLY:
----------
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
---------------------------------
OCCUPIER-Settings:
----------
OCCUPIER_method=auto
OWN_TREE_PURE_WINDOW=3
OWN_progressive_from=no
fob_equal_weights=yes
OCCUPIER_compare=FSPE
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
---------------------------------
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
---------------------------------
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
---------------------------------
ORCA base overrides (optional):
----------
keyword:basename=[]
additions:basename=[]
---------------------------------
CO2 Coordination:
----------
co2_coordination=off
co2_species_delta=0
------------------------------------------------------
calc_prop_of_interest=no
properties_of_interest=IP,EA
reorganisation_energy=lambda_p,lambda_m
------------------------------------------------------
ESD module (excited state dynamics):
------------------
ESD_modul=no
ESD_modus=[TDDFT|deltaSCF|hybrid1]
ESD_T1_opt=[uks|tddft]
ESD_frequency=yes
states=[S1,T1,S2,T2]
ISCs=[S1>T1,T1>S1]
ICs=[S2>S1]
emission_rates=[f,p]
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
---------------------------------
Electrical Properties:
----------
elprop_Dipole=no
elprop_Quadrupole=no
elprop_Hyperpol=no
elprop_Polar=no
elprop_PolarVelocity=no
elprop_PolarDipQuad=no
elprop_PolarQuadQuad=no
---------------------------------
deltaSCF Settings:
----------
deltaSCF_DOMOM=true
deltaSCF_PMOM=false
deltaSCF_keepinitialref=true
deltaSCF_SOSCFHESSUP=LSR1
deltaSCF_keywords=FreezeAndRelease
deltaSCF_maxiter=300
deltaSCF_SOSCFConvFactor=500
deltaSCF_SOSCFMaxStep=0.1
---------------------------------
TDDFT Settings:
----------
TDDFT_TDDFT_maxiter=500
TDDFT_nroots=15
TDDFT_maxdim=30
TDDFT_TDA=TRUE
TDDFT_followiroot=true
TDDFT_SOC=false
------------------------------------------------------
Level of Theory:
----------
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
------------------------------------------------------
Prints:
----------
print_MOs=no
print_Loewdin_population_analysis=no
------------------------------------------------------
Resource Settings:
----------
PAL=48
maxcore=6000
parallel_workflows=yes
pal_jobs=4
orca_parallel_strategy=auto
enable_job_timeouts=no
job_timeout_hours=36
opt_timeout_hours=14
frequency_timeout_hours=36
sp_timeout_hours=3
------------------------------------------------------
Automatic Error Recovery & Retry:
----------
enable_auto_recovery=yes
max_recovery_attempts=3
enable_adaptive_parallelism=yes
enable_performance_metrics=yes
------------------------------------------------------
GUPPY_settings:
----------
GUPPY_RUNS=20
GUPPY_GOAT=0
GUPPY_PARALLEL_JOBS=4
GUPPY_SEED=31
------------------------------------------------------
xTB Hyperpolarizability (sTD-DFT-xTB):
----------
hyperpol_xTB=no
hyperpol_xTB_xyz=start.txt
hyperpol_xTB_preopt=none
hyperpol_xTB_engine=std2
hyperpol_xTB_bfw=no
hyperpol_xTB_wavelengths=
hyperpol_xTB_energy_window=15.0
------------------------------------------------------
xTB TADF Screening:
----------
tadf_xTB=no
tadf_xTB_xyz=start.txt
tadf_xTB_preopt=none
tadf_xTB_excited_method=stda
tadf_xTB_bfw=no
tadf_xTB_energy_window=10.0
tadf_xTB_run_t1_opt=yes
------------------------------------------------------
Thermodynamics:
----------
thermodynamics=no
thermodynamics_mode=[auto|reaction]
thermodynamics_reaction=a*{SMILES}+b*{SMILES}...>>>c*{SMILES}+d*{SMILES}...
n_explicit_solvent=6
logK_exp=
thdy_smiles_converter=[QUICK|NORMAL|GUPPY|ARCHITECTOR]
thdy_preopt=[none|xtb|crest|goat]
"""
# -------------------------------------------------------------------------------------------------------
def _is_coordinate_line(line: str) -> bool:
    """True if `line` looks like an XYZ atom row: 'El x y z' (≥3 numeric cols)."""
    ls = line.strip()
    if not ls or ls == "*":
        return False
    parts = ls.split()
    if len(parts) < 4:
        return False
    if not re.match(r"^[A-Za-z]{1,2}", parts[0]):
        return False
    try:
        float(parts[1]); float(parts[2]); float(parts[3])
    except ValueError:
        return False
    return True
# -------------------------------------------------------------------------------------------------------
def _coordinate_lines(text: str) -> list:
    """Return the atom rows from a coordinate block (also works on a full XYZ:
    the natoms/comment header lines are not coordinate rows and are dropped)."""
    return [ln for ln in text.splitlines() if _is_coordinate_line(ln)]
# -------------------------------------------------------------------------------------------------------
def _has_coordinates(path) -> bool:
    """True if `path` exists and holds at least one XYZ-style atom row."""
    try:
        return bool(_coordinate_lines(path.read_text(encoding="utf-8", errors="ignore")))
    except Exception:
        return False
# -------------------------------------------------------------------------------------------------------
def convert_input_txt_to_xyz(src_txt: str, dst_xyz: str = "input.xyz") -> str:
    """Build an XYZ file from a DELFIN coordinate block by prepending the 2-line
    XYZ header (atom count + comment). Inverse of ``convert_xyz_to_input_txt``;
    used to recover a missing ``input.xyz`` on a recalc whose geometry survived
    only as ``input.txt`` / ``start.txt``."""
    src_path = resolve_path(src_txt)
    dst_path = resolve_path(dst_xyz)

    atom_lines = _coordinate_lines(src_path.read_text(encoding="utf-8", errors="ignore"))
    if not atom_lines:
        raise ValueError(f"No XYZ coordinates found in '{src_txt}'; cannot build '{dst_xyz}'.")

    body = "\n".join(line.strip() for line in atom_lines)
    content = (
        f"{len(atom_lines)}\n"
        f"Regenerated by DELFIN from {src_path.name} (recalc recovery)\n"
        f"{body}\n"
    )
    dst_path.write_text(content, encoding="utf-8")
    message = f"Built '{dst_xyz}' from '{src_txt}' ({len(atom_lines)} atoms, recalc recovery)."
    print(message)
    logger.info(message)
    return dst_xyz
# -------------------------------------------------------------------------------------------------------
def convert_xyz_to_input_txt(src_xyz: str, dst_txt: str = "input.txt") -> str:
    """Convert an XYZ file to input.txt by dropping the first two lines."""
    src_path = resolve_path(src_xyz)
    dst_path = resolve_path(dst_txt)

    if not src_path.exists():
        # Recalc recovery: the .xyz named in CONTROL.txt was not staged to the
        # run directory, but a real coordinate block may still survive as the
        # target .txt (or start.txt). Rebuild the missing .xyz from it instead of
        # clobbering the geometry with an empty file — otherwise the downstream
        # electron-count / OCCUPIER steps abort on an empty geometry.
        recovery_src = None
        if _has_coordinates(dst_path):
            recovery_src = dst_path
        else:
            start_path = src_path.parent / "start.txt"
            if start_path != dst_path and _has_coordinates(start_path):
                recovery_src = start_path
        if recovery_src is not None:
            convert_input_txt_to_xyz(str(recovery_src), str(src_path))
            # Ensure dst_txt also carries the bare coordinate block.
            if recovery_src != dst_path:
                convert_xyz_to_input_txt(str(src_path), str(dst_path))
            logger.warning(
                "Regenerated missing '%s' from '%s' (recalc recovery); "
                "geometry preserved instead of aborting.", src_xyz, recovery_src.name,
            )
            return dst_txt

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

    content = TEMPLATE
    control_path.write_text(content, encoding="utf-8")
    message = f"{filename} has been written (input_file={target_input})."
    print(message)
    logger.info(message)
# -------------------------------------------------------------------------------------------------------
