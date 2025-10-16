from __future__ import annotations
import os, time, sys, argparse
from pathlib import Path

from delfin.common.logging import configure_logging, get_logger
from delfin.common.paths import get_runtime_dir, resolve_path
from delfin.cluster_utils import auto_configure_resources, detect_cluster_environment
from delfin.global_manager import get_global_manager
from .define import create_control_file
from .cleanup import cleanup_all
from .config import read_control_file, get_E_ref
from .utils import search_transition_metals, set_main_basisset, calculate_total_electrons_txt
from .orca import run_orca
from .imag import run_IMAG
from .xtb_crest import XTB, XTB_GOAT, run_crest_workflow, XTB_SOLVATOR
from .energies import find_gibbs_energy, find_ZPE, find_electronic_energy
from .reporting import generate_summary_report_DELFIN as generate_summary_report
from .cli_helpers import _avg_or_none, _build_parser
from .cli_recalc import setup_recalc_mode, patch_modules_for_recalc
from .cli_banner import print_delfin_banner, validate_required_files, get_file_paths
from .pipeline import (
    FileBundle,
    PipelineContext,
    compute_summary,
    interpret_method_alias,
    normalize_input_file,
    run_classic_phase,
    run_manual_phase,
    run_occuper_phase,
)

logger = get_logger(__name__)


def _normalize_input_file(config, control_path: Path) -> str:
    return normalize_input_file(config, control_path)






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

    # Initialize global job manager with configuration
    global_mgr = get_global_manager()
    global_mgr.initialize(config)
    logger.info("Global job manager initialized")

    def _finalize(exit_code: int) -> int:
        """Shutdown global resources and perform optional cleanup before exiting."""
        try:
            global_mgr.shutdown()
            logger.info("Global job manager shutdown complete")
        except Exception as exc:  # noqa: BLE001
            logger.warning(f"Global manager shutdown raised: {exc}")

        if not args.no_cleanup:
            cleanup_all(str(get_runtime_dir()))
        return exit_code

    # Populate optional flags with safe defaults so reduced CONTROL files remain usable
    default_config = {
        'XTB_OPT': 'no',
        'XTB_GOAT': 'no',
        'CREST': 'no',
        'XTB_SOLVATOR': 'no',
        'calc_initial': 'yes',
        'oxidation_steps': '',
        'reduction_steps': '',
        'parallel_workflows': 'yes',
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
        'maxiter': 125,
        'maxiter_occupier': 125,
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

    # Canonical file names used throughout workflows
    file_bundle = FileBundle(
        xyz_initial="initial.xyz",
        xyz_red1="red_step_1.xyz",
        xyz_red2="red_step_2.xyz",
        xyz_red3="red_step_3.xyz",
        xyz_ox1="ox_step_1.xyz",
        xyz_ox2="ox_step_2.xyz",
        xyz_ox3="ox_step_3.xyz",
        output_initial="initial.inp",
        output_absorption="absorption_td.inp",
        output_t1="t1_state_opt.inp",
        output_s1="s1_state_opt.inp",
        output_emission="emission_td.inp",
        output_ox1="ox_step_1.inp",
        output_ox2="ox_step_2.inp",
        output_ox3="ox_step_3.inp",
        output_red1="red_step_1.inp",
        output_red2="red_step_2.inp",
        output_red3="red_step_3.inp",
    )

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

    metals_list = list(metals) if metals else []

    pipeline_ctx = PipelineContext(
        config=config,
        control_file_path=control_file_path,
        input_file=input_file,
        charge=charge,
        PAL=PAL,
        multiplicity=multiplicity,
        solvent=solvent,
        metals=metals_list,
        main_basisset=main_basisset,
        metal_basisset=metal_basisset,
        number_explicit_solv_molecules=number_explicit_solv_molecules,
        total_electrons_txt=total_electrons_txt,
        start_time=start_time,
        name=NAME,
        file_bundle=file_bundle,
    )

    raw_method = str(config.get('method', '')).strip()
    canonical_method, suggestion = interpret_method_alias(raw_method)

    if canonical_method not in {"classic", "manually", "OCCUPIER"}:
        if suggestion is not None:
            logger.error(
                "Unknown method '%s'. Did you mean '%s'?",
                raw_method,
                suggestion,
            )
        else:
            logger.error("Unknown method '%s'. Supported methods: classic, manually, OCCUPIER", raw_method)
        return 2

    if canonical_method != raw_method and canonical_method.lower() != raw_method.lower():
        logger.warning("CONTROL.txt method '%s' interpreted as '%s'.", raw_method, canonical_method)

    config['method'] = canonical_method






    # ------------------- OCCUPIER --------------------
    if config['method'] == "OCCUPIER":
        if not run_occuper_phase(pipeline_ctx):
            return _finalize(1)

    # ------------------- classic --------------------
    if config['method'] == "classic":
        run_classic_phase(pipeline_ctx)


    # ------------------- manually --------------------
    if config['method'] == "manually":
        run_manual_phase(pipeline_ctx)


    # Finalize redox and emission summary
    if config['method'] == "OCCUPIER":
        mul0 = pipeline_ctx.extra.get('multiplicity_0')
        if mul0 is not None:
            try:
                pipeline_ctx.multiplicity = int(mul0)
            except Exception:
                pipeline_ctx.multiplicity = mul0  # fallback to raw value
    elif config['method'] == "manually":
        try:
            pipeline_ctx.multiplicity = int(config.get('multiplicity_0', pipeline_ctx.multiplicity))
        except Exception:
            pass
    elif config['method'] == "classic":
        try:
            total_electrons_txt, mult_guess = calculate_total_electrons_txt(str(control_file_path))
            total_electrons_txt = int(total_electrons_txt)
            total_electrons = total_electrons_txt - pipeline_ctx.charge
            pipeline_ctx.multiplicity = 1 if total_electrons % 2 == 0 else 2
        except Exception:
            pass

    summary = compute_summary(pipeline_ctx, E_ref)

    charge = pipeline_ctx.charge
    solvent = pipeline_ctx.solvent
    metals_list = pipeline_ctx.metals
    main_basis = pipeline_ctx.main_basisset
    metal_basis = pipeline_ctx.metal_basisset

    generate_summary_report(
        charge,
        summary.multiplicity,
        solvent,
        summary.E_ox,
        summary.E_ox_2,
        summary.E_ox_3,
        summary.E_red,
        summary.E_red_2,
        summary.E_red_3,
        summary.E_00_t1,
        summary.E_00_s1,
        metals_list,
        metal_basis,
        NAME,
        main_basis,
        config,
        summary.duration,
        E_ref,
    )

    return _finalize(0)


if __name__ == "__main__":
    sys.exit(main())
