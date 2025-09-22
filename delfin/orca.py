import logging
import subprocess
import sys
from shutil import which
from typing import Optional

ORCA_PLOT_INPUT_TEMPLATE = (
    "1\n"
    "1\n"
    "4\n"
    "100\n"
    "5\n"
    "7\n"
    "2\n"
    "{index}\n"
    "10\n"
    "11\n"
)

def find_orca_executable() -> Optional[str]:
    """Locate ORCA executable in system PATH.

    Returns:
        Path to ORCA executable, or None if not found
    """
    orca_path = which("orca")
    if not orca_path:
        logging.error("ORCA executable not found. Please ensure ORCA is installed and in your PATH.")
        return None
    return orca_path


def _run_orca_subprocess(orca_path: str, input_file_path: str, output_log: str) -> bool:
    """Run ORCA subprocess and capture output. Returns True when successful."""
    with open(output_log, "w") as output_file:
        try:
            subprocess.run([orca_path, input_file_path], check=True, stdout=output_file, stderr=output_file)
        except subprocess.CalledProcessError as error:
            logging.error(f"Error running ORCA: {error}")
            return False
    return True

def run_orca(input_file_path: str, output_log: str) -> None:
    """Execute ORCA calculation with specified input file.

    Runs ORCA subprocess with input file and captures output to log file.
    Logs success/failure and handles subprocess errors.

    Args:
        input_file_path: Path to ORCA input file (.inp)
        output_log: Path for ORCA output file (.out)
    """
    orca_path = find_orca_executable()
    if not orca_path:
        return

    if _run_orca_subprocess(orca_path, input_file_path, output_log):
        logging.info(f"ORCA run successful for '{input_file_path}'")

def run_orca_IMAG(input_file_path: str, iteration: int) -> None:
    """Execute ORCA calculation for imaginary frequency workflow.

    Specialized ORCA runner for IMAG workflow with iteration-specific
    output naming and enhanced error handling.

    Args:
        input_file_path: Path to ORCA input file
        iteration: Iteration number for output file naming
    """
    orca_path = find_orca_executable()
    if not orca_path:
        logging.error("Cannot run ORCA IMAG calculation because the ORCA executable was not found in PATH.")
        sys.exit(1)

    output_log = f"output_{iteration}.out"
    if _run_orca_subprocess(orca_path, input_file_path, output_log):
        logging.info(f"ORCA run successful for '{input_file_path}', output saved to '{output_log}'")
    else:
        sys.exit(1)

def run_orca_plot(homo_index: int) -> None:
    """Generate molecular orbital plots around HOMO using orca_plot.

    Creates orbital plots for orbitals from HOMO-10 to HOMO+10
    using ORCA's orca_plot utility with automated input.

    Args:
        homo_index: Index of the HOMO orbital
    """
    for index in range(homo_index - 10, homo_index + 11):
        success, stderr_output = _run_orca_plot_for_index(index)
        if success:
            logging.info(f"orca_plot ran successfully for index {index}")
        else:
            logging.error(f"orca_plot encountered an error for index {index}: {stderr_output}")


def _run_orca_plot_for_index(index: int) -> tuple[bool, str]:
    """Run orca_plot for a single orbital index and return success flag and stderr."""
    process = subprocess.Popen(
        ["orca_plot", "input.gbw", "-i"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    _, stderr = process.communicate(input=_prepare_orca_plot_input(index))
    return process.returncode == 0, stderr.decode()


def _prepare_orca_plot_input(index: int) -> bytes:
    """Build the scripted user input for orca_plot."""
    return ORCA_PLOT_INPUT_TEMPLATE.format(index=index).encode()
