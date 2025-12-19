# cli_banner.py
# Banner and initial setup utilities for DELFIN CLI

from pathlib import Path

from delfin.common.banners import build_standard_banner


def print_delfin_banner():
    """Print the DELFIN application banner."""
    banner = build_standard_banner(header_indent=4, info_indent=4)
    info = """
With DELFIN, it is possible to automatically identify preferred electron configurations of any hypothetical
or known system, track the changes in orbital occupations upon reduction and oxidation, calculate
redox potentials, and, for closed-shell species, determine S0, S1, T1, and T2 states including automated
intersystem crossings (ISC) and internal conversions (IC), alongside E_00 energies and excited-state redox potentials.
DELFIN does not address the fundamental question of whether a hypothetical system is chemically viable or
synthetically accessible.

To use DELFIN, install ORCA 6.1.1 and add it to your PATH. CREST 3.0.2 is optional.

ORCA 6.1.1 is available free of charge for academic use under the ORCA License Agreement.
It is developed and copyrighted by the Max Planck Society; all non-academic use
requires a separate commercial license from the ORCA team. The full licensing terms
and registration steps are provided via the official ORCA portal:
https://orcaforum.kofo.mpg.de/app.php/dlext/?cat=3. Users are responsible for
ensuring compliance with the ORCA license and applicable export regulations.

xTB 6.7.1+ is distributed under the GNU Lesser General Public License (LGPL v3)
by the Grimme group. Source code and binaries are available from
https://github.com/grimme-lab/xtb. Please observe the LGPL requirements if you
modify or redistribute xTB.

CREST 3.0.2 is released under the GNU General Public License (GPL v3) and
maintained by the Grimme group. See https://crest-lab.github.io/ for downloads
and license details. DELFIN does not bundle CREST or xTB; users must obtain them
separately and comply with their licenses.
""".strip()
    print(f"\n{banner}\n\n{info}\n")

def validate_required_files(config, control_path):
    """Validate that required control and input files exist.

    Args:
        config: Configuration dictionary
        control_path: Path to CONTROL.txt file

    Returns:
        tuple: (success: bool, error_code: int, input_file: str)
    """
    # Check CONTROL.txt
    if not control_path.exists():
        print(f"CONTROL file was not found at {control_path}")
        print("Tip: run `delfin --define` to generate a template, or see `delfin --help` for usage.")
        return False, 2, None

    # Check input file
    input_file_entry = (config.get('input_file') or 'input.txt').strip() or 'input.txt'
    input_path = Path(input_file_entry)
    if not input_path.is_absolute():
        input_path = control_path.parent / input_path
    input_path = input_path.expanduser()
    if not input_path.exists():
        print(f"Input file '{input_path}' was not found.")
        print("Tip: run `delfin --define=your.xyz` to convert an XYZ into input.txt, "
              "or create the file manually and update CONTROL.txt (input_file=...).")
        return False, 2, None

    return True, 0, str(input_path)


def get_file_paths():
    """Get standard file paths used throughout the workflow.

    Returns:
        dict: Dictionary of standard file paths
    """
    return {
        'xyz_files': {
            'initial': "initial.xyz",
            'red_step_1': "red_step_1.xyz",
            'red_step_2': "red_step_2.xyz",
            'ox_step_1': "ox_step_1.xyz",
            'ox_step_2': "ox_step_2.xyz"
        },
        'input_files': {
            'initial': "initial.inp",
            'absorption_td': "absorption_td.inp",
            'e_state_opt': "e_state_opt.inp"
        },
        'output_files': {
            'initial': "initial.out",
            'absorption_td': "absorption_td.out",
            'e_state_opt': "e_state_opt.out",
            'red_step_1': "red_step_1.out",
            'red_step_2': "red_step_2.out",
            'red_step_3': "red_step_3.out",
            'ox_step_1': "ox_step_1.out",
            'ox_step_2': "ox_step_2.out",
            'ox_step_3': "ox_step_3.out"
        }
    }
