# cli_banner.py
# Banner and initial setup utilities for DELFIN CLI

from pathlib import Path


def print_delfin_banner():
    """Print the DELFIN application banner."""
    print("""
                          ******************
                          *     DELFIN     *
                          ******************

    #############################################################
    #                           -***-                           #
    #                          ComPlat                          #
    #                        M. Hartmann                        #
    #           Karlsruhe Institute of Technology (KIT)         #
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
