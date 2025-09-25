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
allow_imaginary_freq=0
------------------------------------
Redox steps:
calc_initial=yes
oxidation_steps=1,2,3
reduction_steps=1,2,3
method=classic|manually|OCCUPIER
calc_potential_method=2
------------------------------------
E_00=no
excitation=s|t
triplet_flag=FALSE
absorption_spec=no
emission_spec=no
NROOTS=15
TDA=FALSE
DONTO=FALSE
DOSOC=FALSE
singlet exitation:
IROOT=1
FOLLOWIROOT=TRUE
mcore_E00=10000
------------------------------------
MANUALLY:
multiplicity_0=
additions_0=m,n
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
ri_soc=RI-SOMF(1X)
relativity=ZORA
aux_jk=def2/J
aux_jk_rel=SARC/J
main_basisset=def2-SVP
main_basisset_rel=ZORA-def2-SVP
metal_basisset=def2-TZVP
metal_basisset_rel=SARC-ZORA-TZVP
first_coordination_sphere_metal_basisset=no
first_coordination_sphere_scale=
geom_opt=OPT
initial_guess=PModel
temperature=298.15
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
maxcore=3800
parallel_workflows=no
pal_jobs=4
------------------------------------
OCCUPIER-Settings:
--------------------
frequency_calculation_OCCUPIER=no
occupier_selection=tolerance|truncation|rounding
occupier_precision=3
occupier_epsilon=5e-4
maxiter_occupier=100
geom_opt_OCCUPIER=OPT
pass_wavefunction=no
approximate_spin_projection_APMethod=2
--------------------
even electron number:
even_seq = [
  {"index": 1, "m": 1, "BS": "",    "from": 0},
  {"index": 2, "m": 1, "BS": "1,1", "from": 1},
  {"index": 3, "m": 1, "BS": "2,2", "from": 2},
  {"index": 4, "m": 3, "BS": "",    "from": 1},
  {"index": 5, "m": 3, "BS": "3,1", "from": 4},
  {"index": 6, "m": 3, "BS": "4,2", "from": 5},
  {"index": 7, "m": 5, "BS": "",    "from": 4},
  {"index": 8, "m": 5, "BS": "5,1", "from": 7},
  {"index": 9, "m": 5, "BS": "6,2", "from": 8}
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
  {"index": 7, "m": 6, "BS": "",    "from": 4},
  {"index": 8, "m": 6, "BS": "6,1", "from": 7},
  {"index": 9, "m": 6, "BS": "7,2", "from": 8}
]

INFOS:
-------------------------------------------------
Available METHODS: classic, manually, OCCUPIER
Available ADDITIONS: %SCF BrokenSym m,n END
Available OX_STEPS: 1 ; 1,2 ; 1,2,3 ; 2 ; 3 ; 2,3 ; 1,3
Available RED_STEPS: 1 ; 1,2 ; 1,2,3 ; 2 ; 3 ; 2,3 ; 1,3
Available IMPLICIT SOLVATION MODELS: CPCM ; CPCMC ; SMD
Available dispersion corrections DISP_CORR: D4 ; D3 ; D3BJ ; D3ZERO ; NONE
Available EXCITATIONS: s (singulet) ; t (triplet) (s is more difficult to converge, there may be no convergence).
E_00 can only be calculated for closed shell systems (use classic or manually!)
EXPLICIT SOLVATION MODEL IS VERY EXPENSIVE, especially in combination with OCCUPIER!!!!!

use yes/no not Yes/No !!!!!!
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
