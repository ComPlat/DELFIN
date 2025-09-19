import os, time, re, sys, logging, argparse
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
from pathlib import Path
from .define import create_control_file
from .cleanup import cleanup_all
from .config import read_control_file, get_E_ref
from .utils import search_transition_metals, set_main_basisset, calculate_total_electrons_txt
from .orca import run_orca
from .imag import run_IMAG
from .xyz_io import (
    read_and_modify_file_1,
    read_xyz_and_create_input2,
    read_xyz_and_create_input3,
    read_xyz_and_create_input4,
)
from .xtb_crest import XTB, XTB_GOAT, run_crest_workflow, XTB_SOLVATOR
from .energies import find_gibbs_energy, find_ZPE, find_electronic_energy
from .report import generate_summary_report_DELFIN as generate_summary_report
from .copy_helpers import read_occupier_file, prepare_occ_folder, prepare_occ_folder_2, copy_if_exists


def _build_parser() -> argparse.ArgumentParser:
   
    description = (
        "DELFIN – DFT-based automated prediction of preferred spin states and associated redox potentials pipeline\n\n"
        "Prerequisites:\n"
        "  • ORCA 6.1.0 installed and available in PATH (academic license required)\n"
        "  • Recommended for some workflows: XTB and CREST available in PATH\n"
        "  • Create and edit CONTROL.txt (or run `delfin --define`) before running calculations\n"
        "  • Input geometry should be in XYZ format (atom count + comment line + coordinates)\n\n"
        "Default behavior:\n"
        "  • If no options are provided, DELFIN runs the calculation pipeline using CONTROL.txt\n"
        "    and the referenced input file.\n\n"
        "Notes on --define:\n"
        "  • If you pass an .xyz file to --define (e.g. --define=foo.xyz), DELFIN will convert it\n"
        "    to 'input.txt' by removing the first two lines and will set input_file=input.txt in CONTROL.txt.\n"
        "  • If you pass a non-.xyz name (e.g. --define=mycoords.txt), an empty file with that name\n"
        "    is created and referenced in CONTROL.txt.\n"
        "  • If you omit a value (just --define), 'input.txt' is created by default.\n\n"
        "Notes on --recalc:\n"
        "  • Only (re)runs external jobs whose output (.out) files are missing or appear incomplete.\n"
        "  • Existing results are preserved; parsing/aggregation is redone from what is on disk.\n"
        "  • A job is considered complete if its .out contains typical ORCA end markers such as\n"
        "    'ORCA TERMINATED NORMALLY'.\n"
    )
    epilog = (
        "Examples:\n"
        "  delfin\n"
        "      Run the calculation pipeline using CONTROL.txt and the referenced input file.\n\n"
        "  delfin --define\n"
        "      Generate CONTROL.txt and an empty input.txt (default) and exit.\n\n"
        "  delfin --define=input.xyz\n"
        "      Convert input.xyz → input.txt (drop first two lines), write CONTROL.txt with\n"
        "      input_file=input.txt, then exit.\n\n"
        "  delfin --cleanup\n"
        "      Remove intermediate files/folders from previous runs and exit.\n\n"
        "  delfin --recalc\n"
        "      Re-parse existing outputs and (re)run only external jobs with missing/incomplete .out files.\n"
    )
    p = argparse.ArgumentParser(
        prog="delfin",
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True,
    )
    p.add_argument(
        "-D", "--define",
        nargs="?", const="input.txt", metavar="INPUTFILE",
        help=("Generate CONTROL.txt and create an input file.\n"
              "If INPUTFILE ends with '.xyz', it will be converted to 'input.txt' by dropping the first two lines.\n"
              "If INPUTFILE is omitted, 'input.txt' is created.")
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite CONTROL.txt and the input file if they already exist."
    )
    p.add_argument(
        "-C", "--cleanup",
        action="store_true",
        help="Clean up intermediate files/folders and exit."
    )
    p.add_argument(
        "-V", "--version",
        action="version",
        version="DELFIN 1.0.0",
        help="Show version and exit."
    )
    p.add_argument(
        "--recalc",
        action="store_true",
        help="Only (re)run external jobs whose .out files are missing or incomplete."
    )
    return p



def main():
    # ---- Parse flags first; --help/--version handled by argparse automatically ----
    parser = _build_parser()
    args, _unknown = parser.parse_known_args(sys.argv[1:])
    RECALC_MODE = bool(getattr(args, "recalc", False))
    os.environ["DELFIN_RECALC"] = "1" if RECALC_MODE else "0"

    if RECALC_MODE:
        # IMPORTANT: override the global bindings so all call sites use the wrappers
        global run_orca, XTB, XTB_GOAT, run_crest_workflow, XTB_SOLVATOR

        _run_orca_real = run_orca
        _XTB_real = XTB
        _XTB_GOAT_real = XTB_GOAT
        _CREST_real = run_crest_workflow
        _SOLV_real = XTB_SOLVATOR


        OK_MARKER = "ORCA TERMINATED NORMALLY"

        def _run_orca_wrapper(inp_file, out_file):
            need = True
            if os.path.exists(out_file):
                try:
                    with open(out_file, "r", errors="ignore") as f:
                        need = (OK_MARKER not in f.read())
                    if not need:
                        logging.info("[recalc] skipping ORCA; %s appears complete.", out_file)
                        return None
                except Exception as e:
                    logging.debug("[recalc] could not check %s (%s) -> will run", out_file, e)
            logging.info("[recalc] (re)running ORCA for %s", out_file)
            return _run_orca_real(inp_file, out_file)


        def _xtb_wrapper(multiplicity, charge, config):
            # Skip if typical XTB artifacts or a marker exist
            artifacts = ("xtbopt.xyz", "xtb.trj", "xtbopt.log")
            marker = Path(".delfin_done_xtb")
            if marker.exists() or any(os.path.exists(a) for a in artifacts):
                logging.info("[recalc] skipping XTB; artifacts/marker found.")
                return None
            res = _XTB_real(multiplicity, charge, config)
            try:
                marker.touch()
            except Exception:
                pass
            return res

        def _goat_wrapper(multiplicity, charge, config):
            artifacts = ("GOAT.txt", "goat.out", "goat.log")
            marker = Path(".delfin_done_goat")
            if marker.exists() or any(os.path.exists(a) for a in artifacts):
                logging.info("[recalc] skipping XTB_GOAT; artifacts/marker found.")
                return None
            res = _XTB_GOAT_real(multiplicity, charge, config)
            try:
                marker.touch()
            except Exception:
                pass
            return res

        def _crest_wrapper(PAL, solvent, charge, multiplicity):
            artifacts = ("crest_conformers.xyz", "crest_best.xyz", "crest.energies", "crest.out")
            marker = Path(".delfin_done_crest")
            if marker.exists() or any(os.path.exists(a) for a in artifacts):
                logging.info("[recalc] skipping CREST; artifacts/marker found.")
                return None
            res = _CREST_real(PAL, solvent, charge, multiplicity)
            try:
                marker.touch()
            except Exception:
                pass
            return res

        def _solv_wrapper(input_path, multiplicity, charge, solvent, n_solv, config):
            # If your implementation produces a specific, stable output, prefer checking for it.
            marker = Path(".delfin_done_xtb_solvator")
            if marker.exists():
                logging.info("[recalc] skipping XTB_SOLVATOR; marker found.")
                return None
            res = _SOLV_real(input_path, multiplicity, charge, solvent, n_solv, config)
            try:
                marker.touch()
            except Exception:
                pass
            return res

        # Swap in wrappers in THIS module
        run_orca = _run_orca_wrapper
        XTB = _xtb_wrapper
        XTB_GOAT = _goat_wrapper
        run_crest_workflow = _crest_wrapper
        XTB_SOLVATOR = _solv_wrapper

        # --- CRUCIAL: also patch the modules that captured their own references ---
        from . import orca as _orca_mod
        _orca_mod.run_orca = _run_orca_wrapper

        from . import occupier as _occupier_mod
        _occupier_mod.run_orca = _run_orca_wrapper


    # Only define template and exit
    if args.define:
        create_control_file(filename="CONTROL.txt",
                            input_file=args.define,
                            overwrite=args.overwrite)
        return 0

    # Only cleanup and exit
    if args.cleanup:
        cleanup_all()
        print("Cleanup done.")
        return 0


    # --------------------- From here: normal pipeline run with banner --------------------
    print("""
                          ******************
                          *     DELFIN     *
                          ******************

    #############################################################
    #                           -***-                           #
    #                        M. Hartmann                        #
    #   Automates ORCA 6.1.0, xTB 6.7.1 and CREST 3.0.2 runs    #
    #                       Version 1.0.0                       #
    #                           -***-                           #
    #############################################################
          
With DELFIN, it is possible to automatically identify preferred electron configurations of any hypothetical 
or known system, track the changes in orbital occupations upon reduction and oxidation, calculate 
redox potentials, and, for closed-shell species, calculate E_00 energies and excited-state redox potentials.
DELFIN does not address the fundamental question of whether a hypothetical system is chemically viable or 
synthetically accessible.

To use DELFIN, install ORCA 6.1.0 and add it to your PATH. CREST 3.0.2 is optional.

ORCA 6.1.0 is available free of charge for academic use.
CREST 3.0 is released under the GNU General Public License (GPL).
""")

    # ---- Friendly checks for missing CONTROL.txt / input file ----
    control_file_path = os.path.join(os.getcwd(), "CONTROL.txt")
    if not os.path.exists(control_file_path):
        print("CONTROL.txt was not found in the current directory.")
        print("Tip: run `delfin --define` to generate a template, or see `delfin --help` for usage.")
        return 2

    # Read CONTROL.txt once and derive all settings from it
    config = read_control_file(control_file_path)
    E_ref = get_E_ref(config)

    # Default to input.txt if not set
    input_file = (config.get('input_file') or 'input.txt').strip() or 'input.txt'
    if not os.path.exists(input_file):
        print(f"Input file '{input_file}' was not found.")
        print("Tip: run `delfin --define=your.xyz` to convert an XYZ into input.txt, "
              "or create the file manually and update CONTROL.txt (input_file=...).")
        return 2 

    NAME = (config.get('NAME') or '').strip()

    # Examples: filenames, parsing, conversions with sensible defaults, etc.
    xyz_file = "initial.xyz"
    xyz_file2 = "red_step_1.xyz"
    xyz_file3 = "red_step_2.xyz"
    xyz_file4 = "ox_step_1.xyz"
    xyz_file8 = "ox_step_2.xyz"

    output_file = "initial.inp"
    output_file3 = "absorption_td.inp"
    output_file4 = "e_state_opt.inp"
    output_file5 = "ox_step_1.inp"
    output_file6 = "red_step_1.inp"
    output_file7 = "red_step_2.inp"
    output_file8 = "red_step_3.inp"
    output_file9 = "ox_step_2.inp"
    output_file10 = "ox_step_3.inp"
    output_file11 = "emission_td.inp"

    try:
        charge = int(str(config.get('charge', 0)).strip())
    except ValueError:
        logging.error("Invalid 'charge' in CONTROL.txt; falling back to 0.")
        charge = 0
    try:
        PAL = int(str(config.get('PAL', 6)).strip())
    except ValueError:
        logging.error("Invalid 'PAL' in CONTROL.txt; falling back to 6.")
        PAL = 6
    try:
        number_explicit_solv_molecules = int(str(config.get('number_explicit_solv_molecules', 0)).strip())
    except ValueError:
        logging.error("Invalid 'number_explicit_solv_molecules'; falling back to 0.")
        number_explicit_solv_molecules = 0

    min_fspe_index = None
    solvent = (config.get('solvent') or '').strip()
    start_time = time.time()

    print(f"used Method: {config.get('method', 'UNDEFINED')}\n")

    metals = search_transition_metals(input_file)
    if metals:
        logging.info("Found transition metals:")
        for metal in metals:
            logging.info(metal)
    else:
        logging.info("No transition metals found in the file.")

    main_basisset, metal_basisset = set_main_basisset(metals, config)

    D45_SET = {
        'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
        'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn'
    }
    use_rel = any(m in D45_SET for m in metals)
    if not use_rel:
        if str(config.get('relativity', '')).lower() != 'none':
            logging.info("3d-only system detected → relativity=none (ZORA/X2C/DKH is deactivated).")
        config['relativity'] = 'none' 


    total_electrons_txt, multiplicity_guess = calculate_total_electrons_txt(control_file_path)
    try:
        total_electrons_txt = int(total_electrons_txt)
    except (TypeError, ValueError):
        logging.error("Could not parse total electrons from CONTROL.txt; assuming 0.")
        total_electrons_txt = 0

    total_electrons = total_electrons_txt - charge
    is_even = (total_electrons % 2 == 0)

    cfg_mult = None
    try:
        cfg_mult_raw = config.get('multiplicity_global_opt') if config is not None else None
        cfg_mult = int(cfg_mult_raw) if cfg_mult_raw not in (None, "") else None
        if cfg_mult is not None and cfg_mult <= 0:
            cfg_mult = None
    except (TypeError, ValueError):
        cfg_mult = None

    try:
        ctl_mult_raw = multiplicity_guess.strip() if isinstance(multiplicity_guess, str) else multiplicity_guess
        ctl_mult = int(ctl_mult_raw) if ctl_mult_raw not in (None, "") else None
        if ctl_mult is not None and ctl_mult <= 0:
            ctl_mult = None
    except (TypeError, ValueError):
        ctl_mult = None

    multiplicity = cfg_mult if cfg_mult is not None else (ctl_mult if ctl_mult is not None else (1 if is_even else 2))






    # ------------------- OCCUPIER --------------------
    if config['method'] == "OCCUPIER":

        if config['XTB_OPT'] == "yes":
            XTB(multiplicity, charge, config)

        if config['XTB_GOAT'] == "yes":
            XTB_GOAT(multiplicity, charge, config)

        if config['CREST'] == "yes":
            run_crest_workflow(PAL, solvent, charge, multiplicity)

        if config['XTB_SOLVATOR'] == "no":

            
            if "yes" in config.get("calc_initial", ""):
                print("\nOCCUPIER for the initial system:\n")
                prepare_occ_folder("initial_OCCUPIER", charge_delta=0)
                        
            if "1" in config.get("oxidation_steps", ""):
                print("\nOCCUPIER for the first oxidation step:\n")
                prepare_occ_folder_2("ox_step_1_OCCUPIER", source_occ_folder="initial_OCCUPIER", charge_delta=+1, config=config)
             
            if "2" in config.get("oxidation_steps", ""):
                print("\nOCCUPIER for the second oxidation step:\n")
                prepare_occ_folder_2("ox_step_2_OCCUPIER", source_occ_folder="ox_step_1_OCCUPIER", charge_delta=+2, config=config)
            
            if "3" in config.get("oxidation_steps", ""):
                print("\nOCCUPIER for the third oxidation step:\n") 
                prepare_occ_folder_2("ox_step_3_OCCUPIER", source_occ_folder="ox_step_2_OCCUPIER", charge_delta=+3, config=config)
                

            if "1" in config.get("reduction_steps", ""):
                print("\nOCCUPIER for the first reduction step:\n")
                prepare_occ_folder_2("red_step_1_OCCUPIER", source_occ_folder="initial_OCCUPIER", charge_delta=-1, config=config)
                 
            if "2" in config.get("reduction_steps", ""):
                print("\nOCCUPIER for the second reduction step:\n")
                prepare_occ_folder_2("red_step_2_OCCUPIER", source_occ_folder="red_step_1_OCCUPIER", charge_delta=-2, config=config)
            
            if "3" in config.get("reduction_steps", ""):
                print("\nOCCUPIER for the third reduction step:\n") 
                prepare_occ_folder_2("red_step_3_OCCUPIER", source_occ_folder="red_step_2_OCCUPIER", charge_delta=-3, config=config)
                

            xyz_file_initial_OCCUPIER = f"input_initial_OCCUPIER.xyz"
            xyz_file_ox_step_1_OCCUPIER = f"input_ox_step_1_OCCUPIER.xyz"
            xyz_file_ox_step_2_OCCUPIER = f"input_ox_step_2_OCCUPIER.xyz"
            xyz_file_ox_step_3_OCCUPIER = f"input_ox_step_3_OCCUPIER.xyz"
            xyz_file_red_step_1_OCCUPIER = f"input_red_step_1_OCCUPIER.xyz"
            xyz_file_red_step_2_OCCUPIER = f"input_red_step_2_OCCUPIER.xyz"
            xyz_file_red_step_3_OCCUPIER = f"input_red_step_3_OCCUPIER.xyz"

            if config['frequency_calculation_OCCUPIER'] == "no":
                if config['calc_initial'] == "yes":
                    multiplicity_0, additions_0, min_fspe_index = read_occupier_file("initial_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                    read_xyz_and_create_input3(xyz_file_initial_OCCUPIER, output_file, charge, multiplicity_0, solvent, metals, metal_basisset, main_basisset, config, additions_0)
                    run_orca(output_file, "initial.out")
                    run_IMAG("initial.out", "initial", charge, multiplicity_0, solvent, metals, config, main_basisset, metal_basisset, additions_0)
                    logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization of the initial system complete!")

                if config['absorption_spec'] == "yes":
                    multiplicity = 1
                    additions=config['additions_TDDFT']
                    read_xyz_and_create_input2(xyz_file_initial_OCCUPIER, output_file3, charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
                    run_orca(output_file3, "absorption_spec.out")
                    logging.info("TD-DFT absorption spectra calculation complete!")   

                if config['E_00'] == "yes":

                    if "t" in config.get("excitation", ""):
                        multiplicity = 3
                        additions=config['additions_TDDFT']
                        read_xyz_and_create_input3(xyz_file, "t1_state_opt.inp", charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
                        run_orca("t1_state_opt.inp", "t1_state_opt.out")
                        logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization of T_1 complete!")
                        if config['emission_spec'] == "yes":
                            multiplicity = 1
                            additions=config['additions_TDDFT']
                            read_xyz_and_create_input2("t1_state_opt.xyz", "emission_t1.inp", charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
                            run_orca("emission_t1.inp", "emission_t1.out")
                            logging.info("TD-DFT T1 emission spectra calculation complete!") 

                    if "s" in config.get("excitation", ""):
                        multiplicity = 1
                        additions=config['additions_TDDFT']
                        read_xyz_and_create_input4(xyz_file, "s1_state_opt.inp", charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
                        run_orca("s1_state_opt.inp", "s1_state_opt.out")
                        logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization of S_1 complete!")
                        if config['emission_spec'] == "yes":
                            multiplicity = 1
                            additions=config['additions_TDDFT']
                            read_xyz_and_create_input2("s1_state_opt.xyz", "emission_s1.inp", charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
                            run_orca("emission_s1.inp", "emission_s1.out")
                            logging.info("TD-DFT S1 emission spectra calculation complete!") 

                if "1" in config['oxidation_steps']:
                    multiplicity_ox1, additions_ox1, min_fspe_index = read_occupier_file("ox_step_1_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                    charge = int(config['charge']) + 1
                    read_xyz_and_create_input3(xyz_file_ox_step_1_OCCUPIER, output_file5, charge, multiplicity_ox1, solvent, metals, metal_basisset, main_basisset, config, additions_ox1)
                    run_orca(output_file5, "ox_step_1.out")
                    logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization cation!")

                if "2" in config['oxidation_steps']:
                    multiplicity_ox2, additions_ox2, min_fspe_index = read_occupier_file("ox_step_2_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                    charge = int(config['charge']) + 2
                    read_xyz_and_create_input3(xyz_file_ox_step_2_OCCUPIER, output_file9, charge, multiplicity_ox2, solvent, metals, metal_basisset, main_basisset, config, additions_ox2)
                    run_orca(output_file9, "ox_step_2.out")
                    logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization cation!")

                if "3" in config['oxidation_steps']:
                    multiplicity_ox3, additions_ox3, min_fspe_index = read_occupier_file("ox_step_3_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                    charge = int(config['charge']) + 3
                    read_xyz_and_create_input3(xyz_file_ox_step_3_OCCUPIER, output_file10, charge, multiplicity_ox3, solvent, metals, metal_basisset, main_basisset, config, additions_ox3)
                    run_orca(output_file10, "ox_step_3.out")
                    logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization cation!")




                if "1" in config['reduction_steps']:
                    multiplicity_red1, additions_red1, min_fspe_index = read_occupier_file("red_step_1_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                    charge = int(config['charge']) - 1
                    read_xyz_and_create_input3(xyz_file_red_step_1_OCCUPIER, output_file6, charge, multiplicity_red1, solvent, metals, metal_basisset, main_basisset, config, additions_red1)
                    run_orca(output_file6, "red_step_1.out")
                    logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization anion!")

                if "2" in config['reduction_steps']:
                    multiplicity_red2, additions_red2, min_fspe_index = read_occupier_file("red_step_2_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                    charge = int(config['charge']) - 2
                    read_xyz_and_create_input3(xyz_file_red_step_2_OCCUPIER, output_file7, charge, multiplicity_red2, solvent, metals, metal_basisset, main_basisset, config, additions_red2)
                    run_orca(output_file7, "red_step_2.out")
                    logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization dianion!")

                if "3" in config['reduction_steps']:
                    multiplicity_red3, additions_red3, min_fspe_index = read_occupier_file("red_step_3_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                    charge = int(config['charge']) - 3
                    read_xyz_and_create_input3(xyz_file_red_step_3_OCCUPIER, output_file8, charge, multiplicity_red3, solvent, metals, metal_basisset, main_basisset, config, additions_red3)
                    run_orca(output_file8, "red_step_3.out")
                    logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization trianion!")

            if str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() == "yes":
                # read OCCUPIER info (unchanged)
                multiplicity_0, additions_0, min_fspe_index = read_occupier_file(
                    "initial_OCCUPIER", "OCCUPIER.txt", None, None, None, config
                )

                # initial
                copy_if_exists("./initial_OCCUPIER", "initial.out", "initial.xyz")

                # oxidation steps (corrected folder names)
                copy_if_exists("./ox_step_1_OCCUPIER", "ox_step_1.out", "ox_step_1.xyz")
                copy_if_exists("./ox_step_2_OCCUPIER", "ox_step_2.out", "ox_step_2.xyz")
                copy_if_exists("./ox_step_3_OCCUPIER", "ox_step_3.out", "ox_step_3.xyz")

                # reduction steps (use the right step folders)
                copy_if_exists("./red_step_1_OCCUPIER", "red_step_1.out", "red_step_1.xyz")
                copy_if_exists("./red_step_2_OCCUPIER", "red_step_2.out", "red_step_2.xyz")
                copy_if_exists("./red_step_3_OCCUPIER", "red_step_3.out", "red_step_3.xyz")

        if config['XTB_SOLVATOR'] == "yes":

            if "yes" in config.get("calc_initial", ""):
                print("\nOCCUPIER for the initial system:\n")
                prepare_occ_folder("initial_OCCUPIER", charge_delta=0)
               
            if "1" in config.get("oxidation_steps", ""):
                print("\nOCCUPIER for the first oxidation step:\n")
                prepare_occ_folder_2("ox_step_1_OCCUPIER", source_occ_folder="initial_OCCUPIER", charge_delta=+1, config=config)
 
            if "2" in config.get("oxidation_steps", ""):
                print("\nOCCUPIER for the second oxidation step:\n")
                prepare_occ_folder_2("ox_step_2_OCCUPIER", source_occ_folder="ox_step_1_OCCUPIER", charge_delta=+2, config=config)

            if "3" in config.get("oxidation_steps", ""):
                print("\nOCCUPIER for the third oxidation step:\n") 
                prepare_occ_folder_2("ox_step_3_OCCUPIER", source_occ_folder="ox_step_2_OCCUPIER", charge_delta=+3, config=config)
                

            if "1" in config.get("reduction_steps", ""):
                print("\nOCCUPIER for the first reduction step:\n")
                prepare_occ_folder_2("red_step_1_OCCUPIER", source_occ_folder="initial_OCCUPIER", charge_delta=-1, config=config)

            if "2" in config.get("reduction_steps", ""):
                print("\nOCCUPIER for the second reduction step:\n") 
                prepare_occ_folder_2("red_step_2_OCCUPIER", source_occ_folder="red_step_1_OCCUPIER", charge_delta=-2, config=config)

            if "3" in config.get("reduction_steps", ""):
                print("\nOCCUPIER for the third reduction step:\n") 
                prepare_occ_folder_2("red_step_3_OCCUPIER", source_occ_folder="red_step_2_OCCUPIER", charge_delta=-3, config=config)
                


            multiplicity_0, additions_0, min_fspe_index = read_occupier_file("initial_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
            XTB_SOLVATOR(str(Path("input_initial_OCCUPIER.xyz").resolve()), multiplicity_0, charge, solvent, number_explicit_solv_molecules, config)


            if config['calc_initial'] == "yes":
                read_and_modify_file_1(input_file, output_file, charge, multiplicity_0, solvent, metals, metal_basisset, main_basisset, config, additions_0)
                run_orca(output_file, "initial.out")
                run_IMAG("initial.out", "initial", charge, multiplicity_0, solvent, metals, config, main_basisset, metal_basisset, additions_0)
                logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization of the initial system complete!")

            if config['absorption_spec'] == "yes":
                multiplicity = 1
                additions=config['additions_TDDFT']
                read_xyz_and_create_input2(xyz_file, output_file3, charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
                run_orca(output_file3, "absorption_spec.out")
                logging.info("TD-DFT absorption spectra calculation complete!")


            if config['E_00'] == "yes":
                additions = config.get('additions_TDDFT', '')
                if "t" in config.get("excitation", ""):
                    multiplicity = 3
                    read_xyz_and_create_input3(xyz_file, "t1_state_opt.inp", charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
                    run_orca("t1_state_opt.inp", "t1_state_opt.out")
                    logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization of T_1 complete!")
                    if config['emission_spec'] == "yes":
                        multiplicity = 1
                        additions=config['additions_TDDFT']
                        read_xyz_and_create_input2("t1_state_opt.xyz", "emission_t1.inp", charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
                        run_orca("emission_t1.inp", "emission_t1.out")
                        logging.info("TD-DFT emission spectra calculation complete!") 

                if "s" in config.get("excitation", ""):
                    multiplicity = 1
                    read_xyz_and_create_input4(xyz_file, "s1_state_opt.inp", charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
                    run_orca("s1_state_opt.inp", "s1_state_opt.out")
                    logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization of S_1 complete!")
                    if config['emission_spec'] == "yes":
                        multiplicity = 1
                        additions=config['additions_TDDFT']
                        read_xyz_and_create_input2("s1_state_opt.xyz", "emission_s1.inp", charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
                        run_orca("emission_s1.inp", "emission_s1.out")
                        logging.info("TD-DFT emission spectra calculation complete!") 



            if "1" in config['oxidation_steps']:
                multiplicity_ox1, additions_ox1, min_fspe_index = read_occupier_file("ox_step_1_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                charge = int(config['charge']) + 1
                read_xyz_and_create_input3(xyz_file, output_file5, charge, multiplicity_ox1, solvent, metals, metal_basisset, main_basisset, config, additions_ox1)
                run_orca(output_file5, "ox_step_1.out")
                logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization cation!")

            if "2" in config['oxidation_steps']:
                multiplicity_ox2, additions_ox2, min_fspe_index = read_occupier_file("ox_step_2_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                charge = int(config['charge']) + 2
                read_xyz_and_create_input3(xyz_file4, output_file9, charge, multiplicity_ox2, solvent, metals, metal_basisset, main_basisset, config, additions_ox2)
                run_orca(output_file9, "ox_step_2.out")
                logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization cation!")

            if "3" in config['oxidation_steps']:
                multiplicity_ox3, additions_ox3, min_fspe_index = read_occupier_file("ox_step_3_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                charge = int(config['charge']) + 3
                read_xyz_and_create_input3(xyz_file8, output_file10, charge, multiplicity_ox3, solvent, metals, metal_basisset, main_basisset, config, additions_ox3)
                run_orca(output_file10, "ox_step_3.out")
                logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization cation!")



            if "1" in config['reduction_steps']:
                multiplicity_red1, additions_red1, min_fspe_index = read_occupier_file("red_step_1_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                charge = int(config['charge']) - 1
                read_xyz_and_create_input3(xyz_file, output_file6, charge, multiplicity_red1, solvent, metals, metal_basisset, main_basisset, config, additions_red1)
                run_orca(output_file6, "red_step_1.out")
                logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization anion!")

            if "2" in config['reduction_steps']:
                multiplicity_red2, additions_red2, min_fspe_index = read_occupier_file("red_step_2_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                charge = int(config['charge']) - 2
                read_xyz_and_create_input3(xyz_file2, output_file7, charge, multiplicity_red2, solvent, metals, metal_basisset, main_basisset, config, additions_red2)
                run_orca(output_file7, "red_step_2.out")
                logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization dianion!")

            if "3" in config['reduction_steps']:
                multiplicity_red3, additions_red3, min_fspe_index = read_occupier_file("red_step_3_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                charge = int(config['charge']) - 3
                read_xyz_and_create_input3(xyz_file3, output_file8, charge, multiplicity_red3, solvent, metals, metal_basisset, main_basisset, config, additions_red3)
                run_orca(output_file8, "red_step_3.out")
                logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization trianion!")




    # ------------------- classic --------------------
    if config['method'] == "classic":

        if config['XTB_OPT'] == "yes":
            XTB(multiplicity, charge, config)

        if config['XTB_GOAT'] == "yes":
            XTB_GOAT(multiplicity, charge, config)

        if config['CREST'] == "yes":
            run_crest_workflow(PAL, solvent, charge, multiplicity)

        if config['XTB_SOLVATOR'] == "yes":
            XTB_SOLVATOR(str(Path("input.txt").resolve()), multiplicity, charge, solvent, number_explicit_solv_molecules, config)



        additions= ""
        if config['calc_initial'] == "yes":
            read_and_modify_file_1(input_file, output_file, charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
            run_orca(output_file, "initial.out")
            run_IMAG("initial.out", "initial", charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
            logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization of the initial system complete!")

        if config['absorption_spec'] == "yes":
            multiplicity = 1
            additions=config['additions_TDDFT']
            read_xyz_and_create_input2(xyz_file, output_file3, charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
            run_orca(output_file3, "absorption_spec.out")
            logging.info("TD-DFT absorption spectra calculation complete!")

        if config['E_00'] == "yes":

            if "t" in config.get("excitation", ""):
                multiplicity = 3
                read_xyz_and_create_input3(xyz_file, "t1_state_opt.inp", charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
                run_orca("t1_state_opt.inp", "t1_state_opt.out")
                logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization of T_1 complete!")
                if config['emission_spec'] == "yes":
                    multiplicity = 1
                    additions=config['additions_TDDFT']
                    read_xyz_and_create_input2("t1_state_opt.xyz", "emission_t1.inp", charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
                    run_orca("emission_t1.inp", "emission_t1.out")
                    logging.info("TD-DFT emission spectra calculation complete!") 

            if "s" in config.get("excitation", ""):
                multiplicity = 1
                read_xyz_and_create_input4(xyz_file, "s1_state_opt.inp", charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
                run_orca("s1_state_opt.inp", "s1_state_opt.out")
                logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization of S_1 complete!")
                if config['emission_spec'] == "yes":
                    multiplicity = 1
                    additions=config['additions_TDDFT']
                    read_xyz_and_create_input2("s1_state_opt.xyz", "emission_s1.inp", charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
                    run_orca("emission_s1.inp", "emission_s1.out")
                    logging.info("TD-DFT emission spectra calculation complete!") 


        if "1" in config['oxidation_steps']:
            charge = int(config['charge']) + 1
            total_electrons = total_electrons_txt - charge
            is_even = total_electrons % 2 == 0
            if is_even:
                multiplicity = 1
            else:
                multiplicity = 2
            read_xyz_and_create_input3(xyz_file, output_file5, charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
            run_orca(output_file5, "ox_step_1.out")
            logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization cation!")


        if "2" in config['oxidation_steps']:
            charge = int(config['charge']) + 2
            total_electrons = total_electrons_txt - charge
            is_even = total_electrons % 2 == 0
            if is_even:
                multiplicity = 1
            else:
                multiplicity = 2
            read_xyz_and_create_input3(xyz_file4, output_file9, charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
            run_orca(output_file9, "ox_step_2.out")
            logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization dication!")

        if "3" in config['oxidation_steps']:
            charge = int(config['charge']) + 3
            total_electrons = total_electrons_txt - charge
            is_even = total_electrons % 2 == 0
            if is_even:
                multiplicity = 1
            else:
                multiplicity = 2
            read_xyz_and_create_input3(xyz_file8, output_file10, charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
            run_orca(output_file10, "ox_step_3.out")
            logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization trication!")



        if "1" in config['reduction_steps']:
            charge = int(config['charge']) - 1
            total_electrons = total_electrons_txt - charge
            is_even = total_electrons % 2 == 0
            if is_even:
                multiplicity = 1
            else:
                multiplicity = 2
            read_xyz_and_create_input3(xyz_file, output_file6, charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
            run_orca(output_file6, "red_step_1.out")
            logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization anion!")

        if "2" in config['reduction_steps']:
            charge = int(config['charge']) - 2
            total_electrons = total_electrons_txt - charge
            is_even = total_electrons % 2 == 0
            if is_even:
                multiplicity = 1
            else:
                multiplicity = 2
            read_xyz_and_create_input3(xyz_file2, output_file7, charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
            run_orca(output_file7, "red_step_2.out")
            logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization dianion!")

        if "3" in config['reduction_steps']:
            charge = int(config['charge']) - 3
            total_electrons = total_electrons_txt - charge
            is_even = total_electrons % 2 == 0
            if is_even:
                multiplicity = 1
            else:
                multiplicity = 2
            read_xyz_and_create_input3(xyz_file3, output_file8, charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
            run_orca(output_file8, "red_step_3.out")
            logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization trianion!")



    # ------------------- manually --------------------
    if config['method'] == "manually":

        multiplicity = config['multiplicity_0']

        if config['XTB_OPT'] == "yes":
            XTB(multiplicity, charge, config)

        if config['XTB_GOAT'] == "yes":
            XTB_GOAT(multiplicity, charge, config)

        if config['CREST'] == "yes":
            run_crest_workflow(PAL, solvent, charge, multiplicity)

        if config['XTB_SOLVATOR'] == "yes":
            XTB_SOLVATOR(str(Path("input.txt").resolve()), multiplicity, charge, solvent, number_explicit_solv_molecules, config)

        wert = config.get('additions_0', "")
        if isinstance(wert, str):
            if re.fullmatch(r"\d+,\d+", wert):
                additions = f"%scf BrokenSym {wert} end"
            else:
                additions = wert
        elif isinstance(wert, list):
            additions = f"%scf BrokenSym {','.join(map(str, wert))} end"
        else:
            additions = ""
        if config['calc_initial'] == "yes":
            read_and_modify_file_1(input_file, output_file, charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
            run_orca(output_file, "initial.out")
            run_IMAG("initial.out", "initial", charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
            logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization of the initial system complete!")

        if config['absorption_spec'] == "yes":
            multiplicity = 1
            additions=config['additions_TDDFT']
            read_xyz_and_create_input2(xyz_file, output_file3, charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
            run_orca(output_file3, "absorption_spec.out")
            logging.info("TD-DFT absorption spectra calculation complete!")   

        if config['E_00'] == "yes":

            if "t" in config.get("excitation", ""):
                multiplicity = 3
                read_xyz_and_create_input3(xyz_file, "t1_state_opt.inp", charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
                run_orca("t1_state_opt.inp", "t1_state_opt.out")
                logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization of T_1 complete!")
                if config['emission_spec'] == "yes":
                    multiplicity = 1
                    additions=config['additions_TDDFT']
                    read_xyz_and_create_input2("t1_state_opt.xyz", "emission_t1.inp", charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
                    run_orca("emission_t1.inp", "emission_t1.out")
                    logging.info("TD-DFT emission spectra calculation complete!") 

            if "s" in config.get("excitation", ""):
                multiplicity = 1
                read_xyz_and_create_input4(xyz_file, "s1_state_opt.inp", charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
                run_orca("s1_state_opt.inp", "s1_state_opt.out")
                logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization of S_1 complete!")
                if config['emission_spec'] == "yes":
                    multiplicity = 1
                    additions=config['additions_TDDFT']
                    read_xyz_and_create_input2("s1_state_opt.xyz", "emission_s1.inp", charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
                    run_orca("emission_s1.inp", "emission_s1.out")
                    logging.info("TD-DFT emission spectra calculation complete!") 


        if "1" in config['oxidation_steps']:
            charge = int(config['charge']) + 1
            multiplicity = config['multiplicity_ox1']
            wert = config.get('additions_ox1', "")

            if isinstance(wert, str):
                if re.fullmatch(r"\d+,\d+", wert):
                    additions = f"%scf BrokenSym {wert} end"
                else:
                    additions = wert
            elif isinstance(wert, list):
                additions = f"%scf BrokenSym {','.join(map(str, wert))} end"
            else:
                additions = ""
            read_xyz_and_create_input3(xyz_file, output_file5, charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
            run_orca(output_file5, "ox_step_1.out")
            logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization cation!")

        if "2" in config['oxidation_steps']:
            charge = int(config['charge']) + 2
            multiplicity = config['multiplicity_ox2']
            wert = config.get('additions_ox2', "")

            if isinstance(wert, str):
                if re.fullmatch(r"\d+,\d+", wert):
                    additions = f"%scf BrokenSym {wert} end"
                else:
                    additions = wert
            elif isinstance(wert, list):
                additions = f"%scf BrokenSym {','.join(map(str, wert))} end"
            else:
                additions = ""
            read_xyz_and_create_input3(xyz_file4, output_file9, charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
            run_orca(output_file9, "ox_step_2.out")
            logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization cation!")

        if "3" in config['oxidation_steps']:
            charge = int(config['charge']) + 3
            multiplicity = config['multiplicity_ox3']
            wert = config.get('additions_ox3', "")

            if isinstance(wert, str):
                if re.fullmatch(r"\d+,\d+", wert):
                    additions = f"%scf BrokenSym {wert} end"
                else:
                    additions = wert
            elif isinstance(wert, list):
                additions = f"%scf BrokenSym {','.join(map(str, wert))} end"
            else:
                additions = ""
            read_xyz_and_create_input3(xyz_file8, output_file10, charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
            run_orca(output_file10, "ox_step_3.out")
            logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization cation!")





        if "1" in config['reduction_steps']:
            charge = int(config['charge']) - 1
            multiplicity = config['multiplicity_red1']
            wert = config.get('additions_red1', "")

            if isinstance(wert, str):
                if re.fullmatch(r"\d+,\d+", wert):
                    additions = f"%scf BrokenSym {wert} end"
                else:
                    additions = wert
            elif isinstance(wert, list):
                additions = f"%scf BrokenSym {','.join(map(str, wert))} end"
            else:
                additions = ""
            read_xyz_and_create_input3(xyz_file, output_file6, charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
            run_orca(output_file6, "red_step_1.out")
            logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization anion!")

        if "2" in config['reduction_steps']:
            charge = int(config['charge']) - 2
            multiplicity = config['multiplicity_red2']
            wert = config.get('additions_red2', "")

            if isinstance(wert, str):
                if re.fullmatch(r"\d+,\d+", wert):
                    additions = f"%scf BrokenSym {wert} end"
                else:
                    additions = wert
            elif isinstance(wert, list):
                additions = f"%scf BrokenSym {','.join(map(str, wert))} end"
            else:
                additions = ""
            read_xyz_and_create_input3(xyz_file2, output_file7, charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
            run_orca(output_file7, "red_step_2.out")
            logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization dianion!")

        if "3" in config['reduction_steps']:
            charge = int(config['charge']) - 3
            multiplicity = config['multiplicity_red3']
            wert = config.get('additions_red3', "")

            if isinstance(wert, str):
                if re.fullmatch(r"\d+,\d+", wert):
                    additions = f"%scf BrokenSym {wert} end"
                else:
                    additions = wert
            elif isinstance(wert, list):
                additions = f"%scf BrokenSym {','.join(map(str, wert))} end"
            else:
                additions = ""
            read_xyz_and_create_input3(xyz_file3, output_file8, charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
            run_orca(output_file8, "red_step_3.out")
            logging.info(f"{config['functional']} {main_basisset} freq & geometry optimization trianion!")






    # Initialize variables to None to avoid UnboundLocalError
    free_gibbs_minus_1 = None
    free_gibbs_minus_2 = None
    free_gibbs_minus_3 = None
    free_gibbs_plus_1 = None
    free_gibbs_plus_2 = None
    free_gibbs_plus_3 = None

    filename0 = 'initial.out'
    filename_s1 = 's1_state_opt.out'
    filename_t1 = 't1_state_opt.out'
    filename_minus_1 = 'red_step_1.out'
    filename_minus_2 = 'red_step_2.out'
    filename_minus_3 = 'red_step_3.out'
    filename_plus_1 = 'ox_step_1.out'
    filemane_plus_2 = 'ox_step_2.out'
    filemane_plus_3 = 'ox_step_3.out'
    filename3 = 'absorption_spec.out'
    filename4 = 'emission_td.out'

    ZPE_S0 = find_ZPE(filename0)
    if ZPE_S0 is not None:
        logging.info(f"ZPE S_0 (eV): {ZPE_S0}")
    #else:
        #logging.error("ZPE S_0 not found in 'initial.out' or conversion failed.")

    ZPE_T1 = None
    ZPE_S1 = None

    if config['E_00'] == "yes":

        if "t" in config.get("excitation", ""):
                ZPE_T1 = find_ZPE(filename_t1)
                if ZPE_T1 is not None:
                    logging.info(f"ZPE T_1 (eV): {ZPE_T1}")
                #else:
                    #logging.error("ZPE T_1 not found in 't1_state_opt.out' or conversion failed.")

        if "s" in config.get("excitation", ""):
                ZPE_S1 = find_ZPE(filename_s1)
                if ZPE_S1 is not None:
                    logging.info(f"ZPE S_1 (eV): {ZPE_S1}")
                #else:
                    #logging.error("ZPE S_1 not found in 's1_state_opt.out' or conversion failed.")

        E_0 = None
        E_S1 = None
        E_T1 = None
        if config['E_00'] == "yes":
            E_0 = find_electronic_energy(filename0)
            if E_0 is not None:
                logging.info(f"Electronic energy S_0 (Eh): {E_0}")
            #else:
                #logging.error("Electronic energy S_0 not found in 'initial.out' or conversion failed.")

        if "t" in config.get("excitation", ""):
            E_T1 = find_electronic_energy(filename_t1)
            if E_T1 is not None:
                logging.info(f"Electronic energy T_1 (Eh): {E_T1}")
            #else:
                #logging.error("Electronic energy S_1 not found in 't1_state_opt.out' or conversion failed.")
        if "s" in config.get("excitation", ""):
            E_S1 = find_electronic_energy(filename_s1)
            if E_S1 is not None:
                logging.info(f"Electronic energy S_1 (Eh): {E_S1}")
            #else:
                #logging.error("Electronic energy S_1 not found in 's1_state_opt.out' or conversion failed.")







    free_gibbs_0 = find_gibbs_energy(filename0)
    if free_gibbs_0 is not None:
        logging.info(f"Free Gibbs Free Energy 0 (H): {free_gibbs_0}")
    #else:
        #logging.error("Final Gibbs free energy not found in 'initial.out' or conversion failed.")



    if "1" in config['oxidation_steps']:
        free_gibbs_plus_1 = find_gibbs_energy(filename_plus_1)
        if free_gibbs_plus_1 is not None:
            logging.info(f"Free Gibbs Free Energy +1 (H): {free_gibbs_plus_1}")
        #else:
            #logging.error("Final Gibbs free energy not found in 'ox_step_1.out' or conversion failed.")

    if "2" in config['oxidation_steps']:
        free_gibbs_plus_2 = find_gibbs_energy(filemane_plus_2)
        if free_gibbs_plus_2 is not None:
            logging.info(f"Free Gibbs Free Energy +2 (H): {free_gibbs_plus_2}")
        #else:
            #logging.error("Final Gibbs free energy not found in 'ox_step_2.out' or conversion failed.")

    if "3" in config['oxidation_steps']:
        free_gibbs_plus_3 = find_gibbs_energy(filemane_plus_3)
        if free_gibbs_plus_3 is not None:
            logging.info(f"Free Gibbs Free Energy +3 (H): {free_gibbs_plus_3}")
        #else:
            #logging.error("Final Gibbs free energy not found in 'ox_step_3.out' or conversion failed.")

    if "1" in config['reduction_steps']:
        free_gibbs_minus_1 = find_gibbs_energy(filename_minus_1)
        if free_gibbs_minus_1 is not None:
            logging.info(f"Free Gibbs Free Energy -1 (H): {free_gibbs_minus_1}")
        #else:
            #logging.error("Final Gibbs free energy not found in 'red_step_1.out' or conversion failed.")

    if "2" in config['reduction_steps']:
        free_gibbs_minus_2 = find_gibbs_energy(filename_minus_2)
        if free_gibbs_minus_2 is not None:
            logging.info(f"Free Gibbs Free Energy -2 (H): {free_gibbs_minus_2}")
        #else:
            #logging.error("Final Gibbs free energy not found in 'red_step_2.out' or conversion failed.")

    if "3" in config['reduction_steps']:
        free_gibbs_minus_3 = find_gibbs_energy(filename_minus_3)   
        if free_gibbs_minus_3 is not None:
            logging.info(f"Free Gibbs Free Energy -3 (H): {free_gibbs_minus_3}")
        #else:
            #logging.error("Final Gibbs free energy not found in 'red_step_3.out' or conversion failed.")



    # ----------------- Calculations ------------------------
    
    E_ox = None
    E_ox_2 = None
    E_ox_3 = None
    E_red = None
    E_red_2 = None
    E_red_3 = None
    E_00_t1 = None
    E_00_s1 = None

    # --- read selection (default -> '2'), no new parser module needed ---
    _sel_raw = config.get('calc_potential_method', config.get('calc_method', 2))
    if isinstance(_sel_raw, (list, tuple, set)):
        _sel_tokens = [str(x) for x in _sel_raw]
    else:
        _sel_tokens = re.split(r'[\s,]+', str(_sel_raw).strip())

    use_m1 = '1' in _sel_tokens
    use_m2 = '2' in _sel_tokens
    use_m3 = '3' in _sel_tokens
    if not (use_m1 or use_m2 or use_m3):
        use_m1 = True  # default

    # ---------- constants ----------
    conv = 2625.499639479947            # kJ/mol per Hartree
    F    = 96.4853321233100184          # kJ/(V·mol)

    # ---------- containers ----------
    m1_avg  = {}  # Method 1: averages (multi-e as average)
    m2_step = {}  # Method 2: step-wise (1e steps)
    m3_mix  = {}  # Method 3: (M1 + M2)/2 as requested

    # ---------- METHOD 1 (averages) ----------
    # Oxidations
    if free_gibbs_0 is not None and free_gibbs_plus_1 is not None:
        m1_avg['E_ox']   = -((free_gibbs_0 - free_gibbs_plus_1) * conv) / F - E_ref
    if free_gibbs_0 is not None and free_gibbs_plus_2 is not None:
        m1_avg['E_ox_2'] = -((free_gibbs_0 - free_gibbs_plus_2) * conv) / (2 * F) - E_ref
    if free_gibbs_0 is not None and free_gibbs_plus_3 is not None:
        m1_avg['E_ox_3'] = -((free_gibbs_0 - free_gibbs_plus_3) * conv) / (3 * F) - E_ref
    # Reductions
    if free_gibbs_0 is not None and free_gibbs_minus_1 is not None:
        m1_avg['E_red']   = ((free_gibbs_0 - free_gibbs_minus_1) * conv) / F - E_ref
    if free_gibbs_0 is not None and free_gibbs_minus_2 is not None:
        m1_avg['E_red_2'] = ((free_gibbs_0 - free_gibbs_minus_2) * conv) / (2 * F) - E_ref
    if free_gibbs_0 is not None and free_gibbs_minus_3 is not None:
        m1_avg['E_red_3'] = ((free_gibbs_0 - free_gibbs_minus_3) * conv) / (3 * F) - E_ref

    # ---------- METHOD 2 (step-wise) ----------
    # Oxidations
    if free_gibbs_0 is not None and free_gibbs_plus_1 is not None:
        m2_step['E_ox']   = -((free_gibbs_0 - free_gibbs_plus_1) * conv) / F - E_ref
    if free_gibbs_plus_1 is not None and free_gibbs_plus_2 is not None:
        m2_step['E_ox_2'] = -((free_gibbs_plus_1 - free_gibbs_plus_2) * conv) / F - E_ref
    if free_gibbs_plus_2 is not None and free_gibbs_plus_3 is not None:
        m2_step['E_ox_3'] = -((free_gibbs_plus_2 - free_gibbs_plus_3) * conv) / F - E_ref
    # Reductions
    if free_gibbs_0 is not None and free_gibbs_minus_1 is not None:
        m2_step['E_red']   = ((free_gibbs_0 - free_gibbs_minus_1) * conv) / F - E_ref
    if free_gibbs_minus_1 is not None and free_gibbs_minus_2 is not None:
        m2_step['E_red_2'] = ((free_gibbs_minus_1 - free_gibbs_minus_2) * conv) / F - E_ref
    if free_gibbs_minus_2 is not None and free_gibbs_minus_3 is not None:
        m2_step['E_red_3'] = ((free_gibbs_minus_2 - free_gibbs_minus_3) * conv) / F - E_ref

    # ---------- METHOD 3 (M1+M2)/2 on the published outputs ----------
    def _avg_or_none(a, b):
        return (a + b) / 2 if (a is not None and b is not None) else None

    # 1e steps: both methods define the same thing; mean equals them
    m3_mix['E_ox']   = _avg_or_none(m1_avg.get('E_ox'),   m2_step.get('E_ox'))
    m3_mix['E_red']  = _avg_or_none(m1_avg.get('E_red'),  m2_step.get('E_red'))
    # 2e/3e: mean of M1-average and M2-step
    m3_mix['E_ox_2']  = _avg_or_none(m1_avg.get('E_ox_2'),  m2_step.get('E_ox_2'))
    m3_mix['E_ox_3']  = _avg_or_none(m1_avg.get('E_ox_3'),  m2_step.get('E_ox_3'))
    m3_mix['E_red_2'] = _avg_or_none(m1_avg.get('E_red_2'), m2_step.get('E_red_2'))
    m3_mix['E_red_3'] = _avg_or_none(m1_avg.get('E_red_3'), m2_step.get('E_red_3'))

    # ---------- FINAL SELECTION: assign legacy variable names for the report ----------
    # Priority: if multiple selected, prefer 3 > 2 > 1. Fallback if missing.
    def _pick(key):
        if use_m3 and m3_mix.get(key) is not None:
            src = 'M3'; val = m3_mix[key]
        elif use_m2 and m2_step.get(key) is not None:
            src = 'M2'; val = m2_step[key]
        elif use_m1 and m1_avg.get(key) is not None:
            src = 'M1'; val = m1_avg[key]
        else:
            # secondary fallbacks (ensure something is set if available)
            for src, bag in (('M3', m3_mix), ('M2', m2_step), ('M1', m1_avg)):
                if bag.get(key) is not None:
                    val = bag[key]
                    logging.warning(f"Using fallback from {src} for {key} (selected method missing).")
                    return val
            val = None
            #logging.warning(f"{key} unavailable in all methods.")
            return val
        logging.info(f"[{src}] {key} = {val}")
        return val

    # Assign the names your report expects:
    E_ox    = _pick('E_ox')
    E_ox_2  = _pick('E_ox_2')
    E_ox_3  = _pick('E_ox_3')
    E_red   = _pick('E_red')
    E_red_2 = _pick('E_red_2')
    E_red_3 = _pick('E_red_3')



    if config['E_00'] == "yes":

        if "t" in config.get("excitation", ""):
            if ZPE_S0 is not None and ZPE_T1 is not None:
                if E_0 is not None and E_T1 is not None:
                    E_00_t1 = ((E_T1 - E_0) + (ZPE_T1 - ZPE_S0)) * 27.211386245988
                    logging.info(f"E_00_t (eV): {E_00_t1}")
                else:
                    logging.error("E_00_t calculation cannot be performed due to missing state1 data.")
            else:
                logging.error("E_00_t calculation cannot be performed due to missing ZPE data.")
        
        if "s" in config.get("excitation", ""):
            if ZPE_S0 is not None and ZPE_S1 is not None:
                if E_0 is not None and E_S1 is not None:
                    E_00_s1 = ((E_S1 - E_0) + (ZPE_S1 - ZPE_S0)) * 27.211386245988
                    logging.info(f"E_00_s (eV): {E_00_s1}")
                else:
                    logging.error("E_00_s calculation cannot be performed due to missing state1 data.")
            else:
                logging.error("E_00_s calculation cannot be performed due to missing ZPE data.")


    charge = int(config['charge'])

    if config['method'] == "OCCUPIER":
        multiplicity = multiplicity_0

    if config['method'] == "manually":
        multiplicity = int(config['multiplicity_0'])

    if config['method'] == "classic":
        total_electrons_txt, multiplicity = calculate_total_electrons_txt(control_file_path)
        total_electrons_txt = int(total_electrons_txt)  # Ensure integer

        total_electrons = total_electrons_txt - charge
        is_even = total_electrons % 2 == 0
        multiplicity = 1 if is_even else 2  # Set correct multiplicity

    E_ref = get_E_ref(config)
    end_time = time.time()
    duration = end_time - start_time
    generate_summary_report(charge, multiplicity, solvent, E_ox, E_ox_2, E_ox_3, E_red, E_red_2, E_red_3, E_00_t1, E_00_s1, metals, metal_basisset, NAME, main_basisset, config, duration, E_ref)

    cleanup_all(".")

if __name__ == "__main__":
    main()


