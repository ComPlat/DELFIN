from __future__ import annotations
import os, time, re, sys, argparse, difflib
from pathlib import Path

from delfin.common.logging import configure_logging, get_logger
from delfin.common.paths import get_runtime_dir, resolve_path
from delfin.cluster_utils import auto_configure_resources, detect_cluster_environment
from delfin.define import convert_xyz_to_input_txt
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
from .reporting import generate_summary_report_DELFIN as generate_summary_report
from .copy_helpers import read_occupier_file, prepare_occ_folder, prepare_occ_folder_2, copy_if_exists
from .thread_safe_helpers import prepare_occ_folder_2_threadsafe
from .cli_helpers import _avg_or_none, _build_parser
from .cli_recalc import setup_recalc_mode, patch_modules_for_recalc
from .cli_banner import print_delfin_banner, validate_required_files, get_file_paths
from .cli_calculations import calculate_redox_potentials, select_final_potentials
from .parallel_classic import (
    execute_classic_parallel_workflows,
    execute_classic_sequential_workflows,
    execute_manually_parallel_workflows,
    normalize_parallel_token,
)
from .parallel_occupier_integration import (
    OccupierExecutionContext,
    run_occupier_orca_jobs,
)

logger = get_logger(__name__)


def _execute_oxidation_workflow(config):
    """Execute oxidation steps workflow."""

    logger.info("Starting oxidation workflow")

    # Remember original working directory
    original_cwd = Path.cwd()

    try:
        if "1" in config.get("oxidation_steps", ""):
            print("\nOCCUPIER for the first oxidation step:\n")
            success = prepare_occ_folder_2_threadsafe("ox_step_1_OCCUPIER", source_occ_folder="initial_OCCUPIER", charge_delta=+1, config=config, original_cwd=original_cwd)
            if not success:
                logger.error("Failed to prepare ox_step_1_OCCUPIER")
                return False

        if "2" in config.get("oxidation_steps", ""):
            print("\nOCCUPIER for the second oxidation step:\n")
            success = prepare_occ_folder_2_threadsafe("ox_step_2_OCCUPIER", source_occ_folder="ox_step_1_OCCUPIER", charge_delta=+2, config=config, original_cwd=original_cwd)
            if not success:
                logger.error("Failed to prepare ox_step_2_OCCUPIER")
                return False

        if "3" in config.get("oxidation_steps", ""):
            print("\nOCCUPIER for the third oxidation step:\n")
            success = prepare_occ_folder_2_threadsafe("ox_step_3_OCCUPIER", source_occ_folder="ox_step_2_OCCUPIER", charge_delta=+3, config=config, original_cwd=original_cwd)
            if not success:
                logger.error("Failed to prepare ox_step_3_OCCUPIER")
                return False

        logger.info("Oxidation workflow completed")
        return True

    except Exception as e:
        logger.error(f"Oxidation workflow failed: {e}")
        return False


def _execute_reduction_workflow(config):
    """Execute reduction steps workflow."""

    logger.info("Starting reduction workflow")

    # Remember original working directory
    original_cwd = Path.cwd()

    try:
        if "1" in config.get("reduction_steps", ""):
            print("\nOCCUPIER for the first reduction step:\n")
            success = prepare_occ_folder_2_threadsafe("red_step_1_OCCUPIER", source_occ_folder="initial_OCCUPIER", charge_delta=-1, config=config, original_cwd=original_cwd)
            if not success:
                logger.error("Failed to prepare red_step_1_OCCUPIER")
                return False

        if "2" in config.get("reduction_steps", ""):
            print("\nOCCUPIER for the second reduction step:\n")
            success = prepare_occ_folder_2_threadsafe("red_step_2_OCCUPIER", source_occ_folder="red_step_1_OCCUPIER", charge_delta=-2, config=config, original_cwd=original_cwd)
            if not success:
                logger.error("Failed to prepare red_step_2_OCCUPIER")
                return False

        if "3" in config.get("reduction_steps", ""):
            print("\nOCCUPIER for the third reduction step:\n")
            success = prepare_occ_folder_2_threadsafe("red_step_3_OCCUPIER", source_occ_folder="red_step_2_OCCUPIER", charge_delta=-3, config=config, original_cwd=original_cwd)
            if not success:
                logger.error("Failed to prepare red_step_3_OCCUPIER")
                return False

        logger.info("Reduction workflow completed")
        return True

    except Exception as e:
        logger.error(f"Reduction workflow failed: {e}")
        return False


def _execute_parallel_workflows(config):
    """Execute oxidation and reduction workflows in parallel."""
    from concurrent.futures import ThreadPoolExecutor, as_completed

    workflows = []

    # Add oxidation workflow if configured
    if config.get("oxidation_steps", ""):
        workflows.append(("oxidation", _execute_oxidation_workflow, config))

    # Add reduction workflow if configured
    if config.get("reduction_steps", ""):
        workflows.append(("reduction", _execute_reduction_workflow, config))

    if not workflows:
        logger.warning("No oxidation or reduction steps configured")
        return True

    logger.info(f"Starting {len(workflows)} workflows in parallel")

    # Execute workflows in parallel
    with ThreadPoolExecutor(max_workers=len(workflows)) as executor:
        # Submit all workflows
        future_to_workflow = {}
        for workflow_name, workflow_func, workflow_config in workflows:
            future = executor.submit(workflow_func, workflow_config)
            future_to_workflow[future] = workflow_name

        # Wait for completion and collect results
        all_success = True
        for future in as_completed(future_to_workflow):
            workflow_name = future_to_workflow[future]
            try:
                success = future.result()
                if success:
                    logger.info(f"{workflow_name} workflow completed successfully")
                else:
                    logger.error(f"{workflow_name} workflow failed")
                    all_success = False
            except Exception as e:
                logger.error(f"{workflow_name} workflow raised exception: {e}")
                all_success = False

    if all_success:
        logger.info("All parallel workflows completed successfully")
    else:
        logger.error("Some parallel workflows failed")

    return all_success


def _execute_sequential_workflows(config):
    """Execute oxidation and reduction workflows sequentially."""
    all_success = True

    # Execute oxidation workflow first
    if config.get("oxidation_steps", ""):
        success = _execute_oxidation_workflow(config)
        if not success:
            all_success = False

    # Execute reduction workflow second
    if config.get("reduction_steps", ""):
        success = _execute_reduction_workflow(config)
        if not success:
            all_success = False

    return all_success


def _normalize_input_file(config, control_path: Path) -> str:
    input_entry = (config.get('input_file') or 'input.txt').strip() or 'input.txt'
    entry_path = Path(input_entry)
    if entry_path.is_absolute():
        input_path = resolve_path(entry_path)
    else:
        input_path = resolve_path(control_path.parent / entry_path)
    if input_path.suffix.lower() == '.xyz':
        target = input_path.with_suffix('.txt')
        convert_xyz_to_input_txt(str(input_path), str(target))
        config['input_file'] = str(target)
        return str(target)
    config['input_file'] = str(input_path)
    return str(input_path)






def main(argv: list[str] | None = None) -> int:
    configure_logging()
    # ---- Parse flags first; --help/--version handled by argparse automatically ----
    parser = _build_parser()
    args, _ = parser.parse_known_args(argv if argv is not None else sys.argv[1:])
    RECALC_MODE = bool(getattr(args, "recalc", False))
    os.environ["DELFIN_RECALC"] = "1" if RECALC_MODE else "0"

    if RECALC_MODE:
        # IMPORTANT: override the global bindings so all call sites use the wrappers
        global run_orca, XTB, XTB_GOAT, run_crest_workflow, XTB_SOLVATOR

        wrappers, reals = setup_recalc_mode()

        # Swap in wrappers in THIS module
        run_orca = wrappers['run_orca']
        XTB = wrappers['XTB']
        XTB_GOAT = wrappers['XTB_GOAT']
        run_crest_workflow = wrappers['run_crest_workflow']
        XTB_SOLVATOR = wrappers['XTB_SOLVATOR']

        # Patch other modules that captured their own references
        patch_modules_for_recalc(wrappers)


    # Only define template and exit
    if args.define:
        create_control_file(filename=str(resolve_path(args.control)),
                            input_file=args.define,
                            overwrite=args.overwrite)
        return 0

    # Only cleanup and exit
    if args.cleanup:
        cleanup_all(str(get_runtime_dir()))
        print("Cleanup done.")
        return 0


    # --------------------- From here: normal pipeline run with banner --------------------
    print_delfin_banner()

    # ---- Friendly checks for missing CONTROL.txt / input file ----
    # Read CONTROL.txt once and derive all settings from it
    control_file_path = resolve_path(args.control)
    try:
        config = read_control_file(str(control_file_path))
    except ValueError as exc:
        logger.error("Invalid CONTROL configuration: %s", exc)
        return 2

    # Auto-configure cluster resources if not explicitly set
    config = auto_configure_resources(config)

    # Populate optional flags with safe defaults so reduced CONTROL files remain usable
    default_config = {
        'XTB_OPT': 'no',
        'XTB_GOAT': 'no',
        'CREST': 'no',
        'XTB_SOLVATOR': 'no',
        'calc_initial': 'yes',
        'oxidation_steps': '',
        'reduction_steps': '',
        'parallel_workflows': 'yes',  # parallel ox/red workflows
        'pal_jobs': None,
        'absorption_spec': 'no',
        'emission_spec': 'no',
        'E_00': 'no',
        'additions_TDDFT': '',
        'DONTO': 'FALSE',
        'DOSOC': 'FALSE',
        'FOLLOWIROOT': 'TRUE',
        'IROOT': '1',
        'NROOTS': '15',
        'TDA': 'FALSE',
        'implicit_solvation_model': 'CPCM',
        'maxcore': 3800,
        'maxiter_occupier': 100,
        'mcore_E00': 10000,
        'multiplicity_0': None,
        'multiplicity_ox1': None,
        'multiplicity_ox2': None,
        'multiplicity_ox3': None,
        'multiplicity_red1': None,
        'multiplicity_red2': None,
        'multiplicity_red3': None,
        'out_files': None,
        'inp_files': None,
    }
    for key, value in default_config.items():
        config.setdefault(key, value)

    # Validate required files
    normalized_input = _normalize_input_file(config, control_file_path)
    success, error_code, _ = validate_required_files(config, control_file_path)
    input_file = normalized_input
    if not success:
        return error_code

    E_ref = get_E_ref(config) 

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
        logger.error("Invalid 'charge' in CONTROL.txt; falling back to 0.")
        charge = 0
    try:
        PAL = int(str(config.get('PAL', 6)).strip())
    except ValueError:
        logger.error("Invalid 'PAL' in CONTROL.txt; falling back to 6.")
        PAL = 6
    try:
        number_explicit_solv_molecules = int(str(config.get('number_explicit_solv_molecules', 0)).strip())
    except ValueError:
        logger.error("Invalid 'number_explicit_solv_molecules'; falling back to 0.")
        number_explicit_solv_molecules = 0

    solvent = (config.get('solvent') or '').strip()
    start_time = time.time()

    print(f"used Method: {config.get('method', 'UNDEFINED')}\n")

    metals = search_transition_metals(input_file)
    if metals:
        logger.info("Found transition metals:")
        for metal in metals:
            logger.info(metal)
    else:
        logger.info("No transition metals found in the file.")

    main_basisset, metal_basisset = set_main_basisset(metals, config)

    D45_SET = {
        'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
        'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn'
    }
    use_rel = any(m in D45_SET for m in metals)
    if not use_rel:
        if str(config.get('relativity', '')).lower() != 'none':
            logger.info("3d-only system detected → relativity=none (ZORA/X2C/DKH is deactivated).")
        config['relativity'] = 'none' 


    total_electrons_txt, multiplicity_guess = calculate_total_electrons_txt(str(control_file_path))
    try:
        total_electrons_txt = int(total_electrons_txt)
    except (TypeError, ValueError):
        logger.error("Could not parse total electrons from CONTROL.txt; assuming 0.")
        total_electrons_txt = 0

    total_electrons = total_electrons_txt - charge
    is_even = (total_electrons % 2 == 0)

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

    # Ensure optional multiplicity slots share the detected ground-state multiplicity by default
    for mult_key in (
        'multiplicity_0',
        'multiplicity_ox1',
        'multiplicity_ox2',
        'multiplicity_ox3',
        'multiplicity_red1',
        'multiplicity_red2',
        'multiplicity_red3',
    ):
        if config.get(mult_key) in (None, ''):
            config[mult_key] = multiplicity

    try:
        multiplicity_0 = int(config.get('multiplicity_0'))
    except (TypeError, ValueError):
        multiplicity_0 = config.get('multiplicity_0')

    raw_method = str(config.get('method', '')).strip()
    typo_corrections = {
        'calssic': 'classic',
    }
    normalized_method = typo_corrections.get(raw_method.lower(), raw_method.lower())

    method_aliases = {
        'classic': 'classic',
        'manually': 'manually',
        'occupier': 'OCCUPIER',
    }

    canonical_method = method_aliases.get(normalized_method)
    if canonical_method is None:
        logger.error(
            "Unknown method '%s'. Supported methods: %s.",
            raw_method,
            ", ".join(sorted(method_aliases.values())),
        )
        if raw_method:
            suggestions = difflib.get_close_matches(raw_method, method_aliases.keys(), n=1)
            if suggestions:
                logger.error("Did you mean '%s'?", suggestions[0])
        return 2

    if canonical_method != raw_method and canonical_method.lower() != raw_method.lower():
        logger.warning("CONTROL.txt method '%s' interpreted as '%s'.", raw_method, canonical_method)

    config['method'] = canonical_method






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
                logger.info("Initial OCCUPIER completed successfully")

            # IMPORTANT: Only start ox/red workflows AFTER initial_OCCUPIER is completely finished
            # Check if we need to run ox/red workflows
            needs_ox_red = config.get("oxidation_steps", "").strip() or config.get("reduction_steps", "").strip()

            if needs_ox_red:
                parallel_mode = normalize_parallel_token(config.get('parallel_workflows', 'auto'))
                has_ox = bool(config.get("oxidation_steps", "").strip())
                has_red = bool(config.get("reduction_steps", "").strip())
                can_parallel = has_ox and has_red

                if parallel_mode == 'enable' or (parallel_mode == 'auto' and can_parallel):
                    if can_parallel:
                        logger.info("Executing oxidation and reduction workflows in parallel")
                        _execute_parallel_workflows(config)
                    else:
                        logger.info("Only one workflow configured → executing sequentially")
                        _execute_sequential_workflows(config)
                else:
                    logger.info("Executing workflows sequentially")
                    _execute_sequential_workflows(config)
            else:
                logger.info("No oxidation or reduction steps configured")
                

            parallel_mode = normalize_parallel_token(config.get('parallel_workflows', 'auto'))
            parallel_workflows_enabled = parallel_mode != 'disable'
            metals_list = list(metals) if isinstance(metals, (list, tuple, set)) else ([metals] if metals else [])

            if str(config.get('frequency_calculation_OCCUPIER', 'no')).lower() != "yes":
                occ_context = OccupierExecutionContext(
                    charge=charge,
                    solvent=solvent,
                    metals=metals_list,
                    main_basisset=main_basisset,
                    metal_basisset=metal_basisset,
                    config=config,
                )
                if not run_occupier_orca_jobs(occ_context, parallel_workflows_enabled):
                    logger.error("OCCUPIER post-processing failed; aborting run")
                    return 1

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
                logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization of the initial system complete!")

            if config['absorption_spec'] == "yes":
                multiplicity = 1
                additions=config['additions_TDDFT']
                read_xyz_and_create_input2(xyz_file, output_file3, charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
                run_orca(output_file3, "absorption_spec.out")
                logger.info("TD-DFT absorption spectra calculation complete!")


            if config['E_00'] == "yes":
                additions = config.get('additions_TDDFT', '')
                if "t" in config.get("excitation", ""):
                    multiplicity = 3
                    read_xyz_and_create_input3(xyz_file, "t1_state_opt.inp", charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
                    run_orca("t1_state_opt.inp", "t1_state_opt.out")
                    logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization of T_1 complete!")
                    if config['emission_spec'] == "yes":
                        multiplicity = 1
                        additions=config['additions_TDDFT']
                        read_xyz_and_create_input2("t1_state_opt.xyz", "emission_t1.inp", charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
                        run_orca("emission_t1.inp", "emission_t1.out")
                        logger.info("TD-DFT emission spectra calculation complete!") 

                if "s" in config.get("excitation", ""):
                    multiplicity = 1
                    read_xyz_and_create_input4(xyz_file, "s1_state_opt.inp", charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
                    run_orca("s1_state_opt.inp", "s1_state_opt.out")
                    logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization of S_1 complete!")
                    if config['emission_spec'] == "yes":
                        multiplicity = 1
                        additions=config['additions_TDDFT']
                        read_xyz_and_create_input2("s1_state_opt.xyz", "emission_s1.inp", charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
                        run_orca("emission_s1.inp", "emission_s1.out")
                        logger.info("TD-DFT emission spectra calculation complete!") 



            if "1" in config['oxidation_steps']:
                multiplicity_ox1, additions_ox1, min_fspe_index = read_occupier_file("ox_step_1_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                charge = int(config['charge']) + 1
                read_xyz_and_create_input3(xyz_file, output_file5, charge, multiplicity_ox1, solvent, metals, metal_basisset, main_basisset, config, additions_ox1)
                run_orca(output_file5, "ox_step_1.out")
                logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization cation!")

            if "2" in config['oxidation_steps']:
                multiplicity_ox2, additions_ox2, min_fspe_index = read_occupier_file("ox_step_2_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                charge = int(config['charge']) + 2
                read_xyz_and_create_input3(xyz_file4, output_file9, charge, multiplicity_ox2, solvent, metals, metal_basisset, main_basisset, config, additions_ox2)
                run_orca(output_file9, "ox_step_2.out")
                logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization cation!")

            if "3" in config['oxidation_steps']:
                multiplicity_ox3, additions_ox3, min_fspe_index = read_occupier_file("ox_step_3_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                charge = int(config['charge']) + 3
                read_xyz_and_create_input3(xyz_file8, output_file10, charge, multiplicity_ox3, solvent, metals, metal_basisset, main_basisset, config, additions_ox3)
                run_orca(output_file10, "ox_step_3.out")
                logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization cation!")



            if "1" in config['reduction_steps']:
                multiplicity_red1, additions_red1, min_fspe_index = read_occupier_file("red_step_1_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                charge = int(config['charge']) - 1
                read_xyz_and_create_input3(xyz_file, output_file6, charge, multiplicity_red1, solvent, metals, metal_basisset, main_basisset, config, additions_red1)
                run_orca(output_file6, "red_step_1.out")
                logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization anion!")

            if "2" in config['reduction_steps']:
                multiplicity_red2, additions_red2, min_fspe_index = read_occupier_file("red_step_2_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                charge = int(config['charge']) - 2
                read_xyz_and_create_input3(xyz_file2, output_file7, charge, multiplicity_red2, solvent, metals, metal_basisset, main_basisset, config, additions_red2)
                run_orca(output_file7, "red_step_2.out")
                logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization dianion!")

            if "3" in config['reduction_steps']:
                multiplicity_red3, additions_red3, min_fspe_index = read_occupier_file("red_step_3_OCCUPIER", "OCCUPIER.txt", None, None, None, config)
                charge = int(config['charge']) - 3
                read_xyz_and_create_input3(xyz_file3, output_file8, charge, multiplicity_red3, solvent, metals, metal_basisset, main_basisset, config, additions_red3)
                run_orca(output_file8, "red_step_3.out")
                logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization trianion!")




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

        ground_multiplicity = multiplicity
        additions_ground = ""
        parallel_cli_mode = normalize_parallel_token(config.get('parallel_workflows', 'auto'))
        parallel_cli_enabled = parallel_cli_mode != 'disable'

        if not parallel_cli_enabled:
            if config['calc_initial'] == "yes":
                read_and_modify_file_1(input_file, output_file, charge, ground_multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions_ground)
                run_orca(output_file, "initial.out")
                run_IMAG("initial.out", "initial", charge, ground_multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions_ground)
                logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization of the initial system complete!")

            if config['absorption_spec'] == "yes":
                additions_td = config.get('additions_TDDFT', '')
                read_xyz_and_create_input2(xyz_file, output_file3, charge, 1, solvent, metals, config, main_basisset, metal_basisset, additions_td)
                run_orca(output_file3, "absorption_spec.out")
                logger.info("TD-DFT absorption spectra calculation complete!")

            if config['E_00'] == "yes":

                if "t" in config.get("excitation", ""):
                    read_xyz_and_create_input3(xyz_file, "t1_state_opt.inp", charge, 3, solvent, metals, metal_basisset, main_basisset, config, additions_ground)
                    run_orca("t1_state_opt.inp", "t1_state_opt.out")
                    logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization of T_1 complete!")
                    if config['emission_spec'] == "yes":
                        additions_td = config.get('additions_TDDFT', '')
                        read_xyz_and_create_input2("t1_state_opt.xyz", "emission_t1.inp", charge, 1, solvent, metals, config, main_basisset, metal_basisset, additions_td)
                        run_orca("emission_t1.inp", "emission_t1.out")
                        logger.info("TD-DFT emission spectra calculation complete!") 

                if "s" in config.get("excitation", ""):
                    read_xyz_and_create_input4(xyz_file, "s1_state_opt.inp", charge, 1, solvent, metals, metal_basisset, main_basisset, config, additions_ground)
                    run_orca("s1_state_opt.inp", "s1_state_opt.out")
                    logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization of S_1 complete!")
                    if config['emission_spec'] == "yes":
                        additions_td = config.get('additions_TDDFT', '')
                        read_xyz_and_create_input2("s1_state_opt.xyz", "emission_s1.inp", charge, 1, solvent, metals, config, main_basisset, metal_basisset, additions_td)
                        run_orca("emission_s1.inp", "emission_s1.out")
                        logger.info("TD-DFT emission spectra calculation complete!") 
        else:
            logger.info("[classic] Initial and excited-state ORCA jobs delegated to parallel scheduler")

        multiplicity = ground_multiplicity
        classic_kwargs = {
            'total_electrons_txt': total_electrons_txt,
            'xyz_file': xyz_file,
            'xyz_file2': xyz_file2,
            'xyz_file3': xyz_file3,
            'xyz_file4': xyz_file4,
            'xyz_file8': xyz_file8,
            'output_file5': output_file5,
            'output_file9': output_file9,
            'output_file10': output_file10,
            'output_file6': output_file6,
            'output_file7': output_file7,
            'output_file8': output_file8,
            'solvent': solvent,
            'metals': metals,
            'metal_basisset': metal_basisset,
            'main_basisset': main_basisset,
            'additions': additions_ground,
            'input_file_path': input_file,
            'output_initial': output_file,
            'ground_multiplicity': ground_multiplicity,
            'include_excited_jobs': parallel_cli_enabled,
        }

        oxidation_requested = bool(str(config.get('oxidation_steps', '')).strip())
        reduction_requested = bool(str(config.get('reduction_steps', '')).strip())

        if oxidation_requested or reduction_requested:
            if parallel_cli_enabled:
                logger.info("[classic] Dispatching oxidation/reduction workflows to parallel scheduler")
                if not execute_classic_parallel_workflows(config, **classic_kwargs):
                    logger.error("Classic parallel execution failed; falling back to sequential mode")
                    execute_classic_sequential_workflows(config, **classic_kwargs)
            else:
                logger.info("[classic] Parallel workflows disabled in CONTROL.txt; running sequentially")
                execute_classic_sequential_workflows(config, **classic_kwargs)
        else:
            logger.info("[classic] No oxidation or reduction steps configured; skipping workflows")


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
        ground_additions = additions
        try:
            ground_multiplicity = int(str(config.get('multiplicity_0', 1)).strip())
        except Exception:
            ground_multiplicity = 1

        parallel_cli_mode = normalize_parallel_token(config.get('parallel_workflows', 'auto'))
        parallel_cli_enabled = parallel_cli_mode != 'disable'

        if not parallel_cli_enabled:
            if config['calc_initial'] == "yes":
                read_and_modify_file_1(input_file, output_file, charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
                run_orca(output_file, "initial.out")
                run_IMAG("initial.out", "initial", charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
                logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization of the initial system complete!")

            if config['absorption_spec'] == "yes":
                multiplicity = 1
                additions=config['additions_TDDFT']
                read_xyz_and_create_input2(xyz_file, output_file3, charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
                run_orca(output_file3, "absorption_spec.out")
                logger.info("TD-DFT absorption spectra calculation complete!")   

            if config['E_00'] == "yes":

                if "t" in config.get("excitation", ""):
                    multiplicity = 3
                    read_xyz_and_create_input3(xyz_file, "t1_state_opt.inp", charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
                    run_orca("t1_state_opt.inp", "t1_state_opt.out")
                    logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization of T_1 complete!")
                    if config['emission_spec'] == "yes":
                        multiplicity = 1
                        additions=config['additions_TDDFT']
                        read_xyz_and_create_input2("t1_state_opt.xyz", "emission_t1.inp", charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
                        run_orca("emission_t1.inp", "emission_t1.out")
                        logger.info("TD-DFT emission spectra calculation complete!") 

                if "s" in config.get("excitation", ""):
                    multiplicity = 1
                    read_xyz_and_create_input4(xyz_file, "s1_state_opt.inp", charge, multiplicity, solvent, metals, metal_basisset, main_basisset, config, additions)
                    run_orca("s1_state_opt.inp", "s1_state_opt.out")
                    logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization of S_1 complete!")
                    if config['emission_spec'] == "yes":
                        multiplicity = 1
                        additions=config['additions_TDDFT']
                        read_xyz_and_create_input2("s1_state_opt.xyz", "emission_s1.inp", charge, multiplicity, solvent, metals, config, main_basisset, metal_basisset, additions)
                        run_orca("emission_s1.inp", "emission_s1.out")
                        logger.info("TD-DFT emission spectra calculation complete!") 
        else:
            logger.info("[manually] Initial and excited-state ORCA jobs delegated to parallel scheduler")

        oxidation_requested = bool(str(config.get('oxidation_steps', '')).strip())
        reduction_requested = bool(str(config.get('reduction_steps', '')).strip())

        manual_kwargs = {
            'total_electrons_txt': total_electrons_txt,
            'xyz_file': xyz_file,
            'xyz_file2': xyz_file2,
            'xyz_file3': xyz_file3,
            'xyz_file4': xyz_file4,
            'xyz_file8': xyz_file8,
            'output_file5': output_file5,
            'output_file9': output_file9,
            'output_file10': output_file10,
            'output_file6': output_file6,
            'output_file7': output_file7,
            'output_file8': output_file8,
            'solvent': solvent,
            'metals': metals,
            'metal_basisset': metal_basisset,
            'main_basisset': main_basisset,
            'additions': ground_additions,
            'input_file_path': input_file,
            'output_initial': output_file,
            'ground_multiplicity': ground_multiplicity,
            'ground_additions': ground_additions,
            'include_excited_jobs': parallel_cli_enabled,
        }

        run_manual_sequential = True
        if oxidation_requested or reduction_requested:
            if parallel_cli_enabled:
                logger.info("[manually] Dispatching oxidation/reduction workflows to parallel scheduler")
                if execute_manually_parallel_workflows(config, **manual_kwargs):
                    run_manual_sequential = False
                else:
                    logger.error("Manual parallel execution failed; falling back to sequential mode")
            else:
                logger.info("[manually] Parallel workflows disabled in CONTROL.txt; running sequentially")

        if run_manual_sequential:
            if not (oxidation_requested or reduction_requested):
                logger.info("[manually] No oxidation or reduction steps configured; skipping workflows")

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
                logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization cation!")

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
                logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization cation!")

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
                logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization cation!")

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
                logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization anion!")

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
                logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization dianion!")

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
                logger.info(f"{config['functional']} {main_basisset} freq & geometry optimization trianion!")






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
        logger.info(f"ZPE S_0 (eV): {ZPE_S0}")
    #else:
        #logger.error("ZPE S_0 not found in 'initial.out' or conversion failed.")

    ZPE_T1 = None
    ZPE_S1 = None

    E_0 = None
    E_S1 = None
    E_T1 = None

    if config['E_00'] == "yes":

        if "t" in config.get("excitation", ""):
                ZPE_T1 = find_ZPE(filename_t1)
                if ZPE_T1 is not None:
                    logger.info(f"ZPE T_1 (eV): {ZPE_T1}")
                #else:
                    #logger.error("ZPE T_1 not found in 't1_state_opt.out' or conversion failed.")

        if "s" in config.get("excitation", ""):
                ZPE_S1 = find_ZPE(filename_s1)
                if ZPE_S1 is not None:
                    logger.info(f"ZPE S_1 (eV): {ZPE_S1}")
                #else:
                    #logger.error("ZPE S_1 not found in 's1_state_opt.out' or conversion failed.")

        if config['E_00'] == "yes":
            E_0 = find_electronic_energy(filename0)
            if E_0 is not None:
                logger.info(f"Electronic energy S_0 (Eh): {E_0}")
            #else:
                #logger.error("Electronic energy S_0 not found in 'initial.out' or conversion failed.")

        if "t" in config.get("excitation", ""):
            E_T1 = find_electronic_energy(filename_t1)
            if E_T1 is not None:
                logger.info(f"Electronic energy T_1 (Eh): {E_T1}")
            #else:
                #logger.error("Electronic energy S_1 not found in 't1_state_opt.out' or conversion failed.")
        if "s" in config.get("excitation", ""):
            E_S1 = find_electronic_energy(filename_s1)
            if E_S1 is not None:
                logger.info(f"Electronic energy S_1 (Eh): {E_S1}")
            #else:
                #logger.error("Electronic energy S_1 not found in 's1_state_opt.out' or conversion failed.")







    free_gibbs_0 = find_gibbs_energy(filename0)
    if free_gibbs_0 is not None:
        logger.info(f"Free Gibbs Free Energy 0 (H): {free_gibbs_0}")
    #else:
        #logger.error("Final Gibbs free energy not found in 'initial.out' or conversion failed.")



    if "1" in config['oxidation_steps']:
        free_gibbs_plus_1 = find_gibbs_energy(filename_plus_1)
        if free_gibbs_plus_1 is not None:
            logger.info(f"Free Gibbs Free Energy +1 (H): {free_gibbs_plus_1}")
        #else:
            #logger.error("Final Gibbs free energy not found in 'ox_step_1.out' or conversion failed.")

    if "2" in config['oxidation_steps']:
        free_gibbs_plus_2 = find_gibbs_energy(filemane_plus_2)
        if free_gibbs_plus_2 is not None:
            logger.info(f"Free Gibbs Free Energy +2 (H): {free_gibbs_plus_2}")
        #else:
            #logger.error("Final Gibbs free energy not found in 'ox_step_2.out' or conversion failed.")

    if "3" in config['oxidation_steps']:
        free_gibbs_plus_3 = find_gibbs_energy(filemane_plus_3)
        if free_gibbs_plus_3 is not None:
            logger.info(f"Free Gibbs Free Energy +3 (H): {free_gibbs_plus_3}")
        #else:
            #logger.error("Final Gibbs free energy not found in 'ox_step_3.out' or conversion failed.")

    if "1" in config['reduction_steps']:
        free_gibbs_minus_1 = find_gibbs_energy(filename_minus_1)
        if free_gibbs_minus_1 is not None:
            logger.info(f"Free Gibbs Free Energy -1 (H): {free_gibbs_minus_1}")
        #else:
            #logger.error("Final Gibbs free energy not found in 'red_step_1.out' or conversion failed.")

    if "2" in config['reduction_steps']:
        free_gibbs_minus_2 = find_gibbs_energy(filename_minus_2)
        if free_gibbs_minus_2 is not None:
            logger.info(f"Free Gibbs Free Energy -2 (H): {free_gibbs_minus_2}")
        #else:
            #logger.error("Final Gibbs free energy not found in 'red_step_2.out' or conversion failed.")

    if "3" in config['reduction_steps']:
        free_gibbs_minus_3 = find_gibbs_energy(filename_minus_3)   
        if free_gibbs_minus_3 is not None:
            logger.info(f"Free Gibbs Free Energy -3 (H): {free_gibbs_minus_3}")
        #else:
            #logger.error("Final Gibbs free energy not found in 'red_step_3.out' or conversion failed.")



    # ----------------- Calculations ------------------------

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
    def _pick(_key):
        if use_m3 and m3_mix.get(_key) is not None:
            src = 'M3'; val = m3_mix[_key]
        elif use_m2 and m2_step.get(_key) is not None:
            src = 'M2'; val = m2_step[_key]
        elif use_m1 and m1_avg.get(_key) is not None:
            src = 'M1'; val = m1_avg[_key]
        else:
            # secondary fallbacks (ensure something is set if available)
            for src, bag in (('M3', m3_mix), ('M2', m2_step), ('M1', m1_avg)):
                if bag.get(_key) is not None:
                    val = bag[_key]
                    logger.warning(f"Using fallback from {src} for {_key} (selected method missing).")
                    return val
            val = None
            #logger.warning(f"{_key} unavailable in all methods.")
            return val
        logger.info(f"[{src}] {_key} = {val}")
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
                    logger.info(f"E_00_t (eV): {E_00_t1}")
                else:
                    logger.error("E_00_t calculation cannot be performed due to missing state1 data.")
            else:
                logger.error("E_00_t calculation cannot be performed due to missing ZPE data.")
        
        if "s" in config.get("excitation", ""):
            if ZPE_S0 is not None and ZPE_S1 is not None:
                if E_0 is not None and E_S1 is not None:
                    E_00_s1 = ((E_S1 - E_0) + (ZPE_S1 - ZPE_S0)) * 27.211386245988
                    logger.info(f"E_00_s (eV): {E_00_s1}")
                else:
                    logger.error("E_00_s calculation cannot be performed due to missing state1 data.")
            else:
                logger.error("E_00_s calculation cannot be performed due to missing ZPE data.")


    charge = int(config['charge'])

    if config['method'] == "OCCUPIER":
        multiplicity = multiplicity_0

    if config['method'] == "manually":
        multiplicity = int(config['multiplicity_0'])

    if config['method'] == "classic":
        total_electrons_txt, multiplicity = calculate_total_electrons_txt(str(control_file_path))
        total_electrons_txt = int(total_electrons_txt)  # Ensure integer

        total_electrons = total_electrons_txt - charge
        is_even = total_electrons % 2 == 0
        multiplicity = 1 if is_even else 2  # Set correct multiplicity

    E_ref = get_E_ref(config)
    end_time = time.time()
    duration = end_time - start_time
    generate_summary_report(charge, multiplicity, solvent, E_ox, E_ox_2, E_ox_3, E_red, E_red_2, E_red_3, E_00_t1, E_00_s1, metals, metal_basisset, NAME, main_basisset, config, duration, E_ref)

    if not args.no_cleanup:
        cleanup_all(str(get_runtime_dir()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
