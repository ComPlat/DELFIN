# cli_helpers.py
# Helper functions extracted from cli.py for better code organization

import argparse
from typing import Optional

from delfin import __version__


def _avg_or_none(a: Optional[float], b: Optional[float]) -> Optional[float]:
    """Calculate average of two values if both are not None."""
    return (a + b) / 2 if (a is not None and b is not None) else None


def _build_parser() -> argparse.ArgumentParser:
    """Build and configure the DELFIN command line argument parser."""
    description = (
        "DELFIN – automated ORCA/xTB/CREST workflow for preferred spin states, redox potentials,\n"
        "and optional excited-state analysis.\n\n"
        "Typical usage:\n"
        "  delfin [WORKSPACE]\n"
        "  delfin [WORKSPACE] --recalc\n"
        "  delfin [WORKSPACE] --report {text,docx}\n\n"
        "Subcommands:\n"
        "  cleanup   Remove DELFIN artifacts / stop ORCA jobs (see: delfin cleanup --help)\n"
        "  stop      Gracefully stop running DELFIN processes in a workspace\n"
        "  co2       CO2 Coordinator workflow (template + scan runner)\n\n"
        "Prerequisites:\n"
        "  • ORCA 6.1.1 installed and available in PATH (academic license required)\n"
        "  • Recommended for some workflows: XTB and CREST available in PATH\n"
        "  • Create and edit CONTROL.txt (or run `delfin --define`) before running calculations\n"
        "  • Input geometry should be in XYZ format (atom count + comment line + coordinates)\n"
        "  • Optional: insert a line containing only '$' to split QM/MM regions – DELFIN keeps the\n"
        "    resulting QM/XTB setup for all follow-up jobs automatically\n\n"
        "Default behavior:\n"
        "  • If no options are provided, DELFIN runs the calculation pipeline using CONTROL.txt\n"
        "    and the referenced input file.\n"
        "  • The global scheduler attaches `delfin_run.log` in the working directory and OCCUPIER\n"
        "    subprocesses emit their own `occupier.log` alongside the standard ORCA outputs.\n\n"
        "Project/workspace handling:\n"
        "  • WORKSPACE is an optional directory argument. If given, DELFIN looks for WORKSPACE/CONTROL.txt.\n"
        "  • Use --control to point to a specific CONTROL file anywhere.\n\n"
        "Notes on --define:\n"
        "  • If you pass an .xyz file to --define (e.g. --define=foo.xyz), DELFIN will convert it\n"
        "    to 'input.txt' by removing the first two lines and will set input_file=input.txt in CONTROL.txt.\n"
        "  • If you pass a non-.xyz name (e.g. --define=mycoords.txt), an empty file with that name\n"
        "    is created and referenced in CONTROL.txt.\n"
        "  • If you omit a value (just --define), 'input.txt' is created by default.\n\n"
        "Notes on --control:\n"
        "  • Use --control=/path/to/CONTROL.txt to run with a CONTROL file outside the current directory.\n"
        "  • If the referenced CONTROL file points to an .xyz geometry, DELFIN automatically converts it\n"
        "    to a matching .txt input before the run starts (atom count/comment lines are removed).\n"
        "  • This is particularly useful for staged HPC jobs or batch workflows.\n\n"
        "Parallel workflows:\n"
        "  • Set `parallel_workflows=yes` (default) to let all jobs (oxidation/reduction steps, initial,\n"
        "    TDDFT, OCCUPIER, XTB/CREST) share PAL dynamically and run in parallel.\n"
        "  • Set `parallel_workflows=no` to force sequential execution; the scheduler still enforces PAL\n"
        "    globally but runs one job at a time.\n\n"
        "Notes on cleanup:\n"
        "  • `--cleanup` removes intermediate files from the runtime directory and exits.\n"
        "  • `delfin cleanup --orca` additionally terminates hanging ORCA jobs and purges OCCUPIER folders.\n\n"
        "Notes on --no-cleanup:\n"
        "  • Skips removal of intermediate files at the end of a run (they reside in DELFIN_SCRATCH when set).\n"
        "  • Handy when debugging or inspecting intermediates after automated runs.\n\n"
        "Notes on --purge:\n"
        "  • Removes only DELFIN-generated artifacts (OCCUPIER folders, ORCA inputs/outputs, logs).\n"
        "  • CONTROL.txt and the referenced input file are always preserved; unknown files are skipped.\n"
        "  • Always asks for confirmation before deleting files.\n\n"
        "Notes on --recalc:\n"
        "  • Only (re)runs external jobs whose output (.out) files are missing or appear incomplete.\n"
        "  • Existing results are preserved; parsing/aggregation is redone from what is on disk.\n"
        "  • A job is considered complete if its .out contains typical ORCA end markers such as\n"
        "    'ORCA TERMINATED NORMALLY'.\n\n"
        "Notes on --imag:\n"
        "  • Only runs IMAG elimination on existing .out/.hess files (skips OCCUPIER/classic/manually workflows).\n"
        "  • Respects IMAG_scope setting in CONTROL.txt (initial or all).\n"
        "  • After IMAG elimination completes, automatically generates a report with updated potentials.\n"
        "  • Useful when you already have converged structures but need to eliminate imaginary frequencies.\n"
    )
    epilog = (
        "Examples:\n"
        "  delfin\n"
        "      Run the pipeline in the current directory (expects ./CONTROL.txt).\n\n"
        "  delfin /path/to/project\n"
        "      Run the pipeline in /path/to/project (expects /path/to/project/CONTROL.txt).\n\n"
        "  delfin --define\n"
        "      Generate CONTROL.txt and an empty input.txt (default) and exit.\n\n"
        "  delfin /path/to/project --define\n"
        "      Generate /path/to/project/CONTROL.txt and /path/to/project/input.txt.\n\n"
        "  delfin --define=input.xyz\n"
        "      Convert input.xyz → input.txt (drop first two lines), write CONTROL.txt with\n"
        "      input_file=input.txt, then exit.\n\n"
        "  delfin --cleanup\n"
        "      Remove intermediate files/folders from previous runs and exit.\n\n"
        "  delfin cleanup --orca\n"
        "      Terminate running ORCA jobs, purge OCCUPIER folders, and clean up temporary files.\n\n"
        "  delfin stop --workspace /path/to/project\n"
        "      Gracefully stop running DELFIN processes associated with the workspace.\n\n"
        "  delfin co2 --define\n"
        "      Write CO2 coordinator templates (CONTROL + co2.xyz) and exit.\n\n"
        "  delfin --purge\n"
        "      Delete DELFIN artifacts in the working directory (keeps CONTROL.txt/input, skips unfamiliar files).\n\n"
        "  delfin --recalc\n"
        "      Re-parse existing outputs and (re)run only missing/incomplete jobs.\n\n"
        "  delfin /path/to/project --recalc\n"
        "      Same as above, but for a different workspace.\n\n"
        "  delfin --json\n"
        "      Build DELFIN_Data.json from on-disk results and exit.\n\n"
        "  delfin --afp\n"
        "      Generate AFP_spectrum.png from existing outputs and exit.\n\n"
        "  delfin --imag\n"
        "      Run IMAG elimination on existing .out/.hess files, then generate report with updated potentials.\n"
    )
    p = argparse.ArgumentParser(
        prog="delfin",
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True,
    )
    p.add_argument(
        "workspace",
        nargs="?",
        default=".",
        help="Workspace directory (default: current directory). If set, CONTROL.txt is read from there.",
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
        "--control",
        default="CONTROL.txt",
        metavar="FILE",
        help="Path to CONTROL file (default: CONTROL.txt). Use this to point to a CONTROL file outside WORKSPACE.",
    )
    p.add_argument(
        "-C", "--cleanup",
        action="store_true",
        help="Clean up intermediate files/folders and exit (same as `delfin cleanup`)."
    )
    p.add_argument(
        "--purge",
        action="store_true",
        help="Remove DELFIN-generated artifacts (OCCUPIER folders, ORCA inputs/outputs, logs) while keeping CONTROL.txt and the referenced input file."
    )
    p.add_argument(
        "--no-cleanup",
        action="store_true",
        help="Keep intermediate files instead of removing them at the end of the run."
    )
    p.add_argument(
        "-V", "--version",
        action="version",
        version=f"DELFIN {__version__}",
        help="Show version and exit."
    )
    p.add_argument(
        "--recalc",
        action="store_true",
        help="Only (re)run external jobs whose .out files are missing or incomplete."
    )
    p.add_argument(
        "--occupier-override",
        action="append",
        default=[],
        metavar="STAGE=INDEX",
        help="In --recalc mode, force OCCUPIER Preferred Index for a stage (e.g., red_step_2_OCCUPIER=2). Can be passed multiple times."
    )
    p.add_argument(
        "--report",
        nargs="?",
        const="text",
        choices=["text", "docx"],
        help="Recompute redox potentials from existing output files (text) or build DELFIN.docx (docx)."
    )
    p.add_argument(
        "--imag",
        action="store_true",
        help="Only run IMAG elimination on existing .out/.hess files, then generate report."
    )
    p.add_argument(
        "--json",
        action="store_true",
        help="Generate DELFIN_Data.json for the project (if used alone, only JSON is built and then exit)."
    )
    p.add_argument(
        "--json-output",
        metavar="FILE",
        help="Optional custom path for DELFIN_Data.json (default: <workspace>/DELFIN_Data.json)."
    )
    p.add_argument(
        "--afp",
        action="store_true",
        help="Generate combined Absorption-Fluorescence-Phosphorescence spectrum plot (AFP_spectrum.png)."
    )
    p.add_argument(
        "--afp-fwhm",
        type=float,
        default=50.0,
        metavar="NM",
        help="FWHM for Gaussian broadening in AFP plot (default: 50.0 nm)."
    )
    return p
