# cli_helpers.py
# Helper functions extracted from cli.py for better code organization

import argparse
from typing import Optional, Dict, Any


def _avg_or_none(a: Optional[float], b: Optional[float]) -> Optional[float]:
    """Calculate average of two values if both are not None."""
    return (a + b) / 2 if (a is not None and b is not None) else None


def _build_parser() -> argparse.ArgumentParser:
    """Build and configure the DELFIN command line argument parser."""
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
        "Notes on --control:\n"
        "  • Use --control=/path/to/CONTROL.txt to run with a CONTROL file outside the current directory.\n"
        "  • If the referenced CONTROL file points to an .xyz geometry, DELFIN automatically converts it\n"
        "    to a matching .txt input before the run starts (atom count/comment lines are removed).\n"
        "  • This is particularly useful for staged HPC jobs or batch workflows.\n\n"
        "Notes on --no-cleanup:\n"
        "  • Skips removal of intermediate files at the end of a run (they reside in DELFIN_SCRATCH when set).\n"
        "  • Handy when debugging or inspecting intermediates after automated runs.\n\n"
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
        "--control",
        default="CONTROL.txt",
        metavar="FILE",
        help="Path to CONTROL file (default: CONTROL.txt)."
    )
    p.add_argument(
        "-C", "--cleanup",
        action="store_true",
        help="Clean up intermediate files/folders and exit."
    )
    p.add_argument(
        "--no-cleanup",
        action="store_true",
        help="Keep intermediate files instead of removing them at the end of the run."
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
