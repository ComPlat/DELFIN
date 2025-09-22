# cli_banner.py
# Banner and initial setup utilities for DELFIN CLI

import os


def print_delfin_banner():
    """Print the DELFIN application banner."""
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


def validate_required_files(config):
    """Validate that required control and input files exist.

    Args:
        config: Configuration dictionary

    Returns:
        tuple: (success: bool, error_code: int, input_file: str)
    """
    # Check CONTROL.txt
    control_file_path = os.path.join(os.getcwd(), "CONTROL.txt")
    if not os.path.exists(control_file_path):
        print("CONTROL.txt was not found in the current directory.")
        print("Tip: run `delfin --define` to generate a template, or see `delfin --help` for usage.")
        return False, 2, None

    # Check input file
    input_file = (config.get('input_file') or 'input.txt').strip() or 'input.txt'
    if not os.path.exists(input_file):
        print(f"Input file '{input_file}' was not found.")
        print("Tip: run `delfin --define=your.xyz` to convert an XYZ into input.txt, "
              "or create the file manually and update CONTROL.txt (input_file=...).")
        return False, 2, None

    return True, 0, input_file


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