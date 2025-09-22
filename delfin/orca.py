import platform
import subprocess
import logging
import sys
from typing import Optional

def find_orca_executable() -> Optional[str]:
    """Locate ORCA executable in system PATH.

    Uses platform-appropriate command ('which' on Unix, 'where' on Windows)
    to find the ORCA executable.

    Returns:
        Path to ORCA executable, or None if not found
    """
    command = "where" if platform.system() == "Windows" else "which"
    result = subprocess.run([command, "orca"], capture_output=True, text=True)
    orca_path = result.stdout.strip()
    if not orca_path:
        logging.error("ORCA executable not found. Please ensure ORCA is installed and in your PATH.")
        return None
    return orca_path

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
    with open(output_log, "w") as output_file:
        try:
            subprocess.run([orca_path, input_file_path], check=True, stdout=output_file, stderr=output_file)
            logging.info(f"ORCA run successful for '{input_file_path}'")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running ORCA: {e}")

def run_orca_IMAG(input_file_path: str, iteration: int) -> None:
    """Execute ORCA calculation for imaginary frequency workflow.

    Specialized ORCA runner for IMAG workflow with iteration-specific
    output naming and enhanced error handling.

    Args:
        input_file_path: Path to ORCA input file
        iteration: Iteration number for output file naming
    """
    orca_path = find_orca_executable()
    output_log = f"output_{iteration}.out"
    with open(output_log, "w") as output_file:
        try:
            subprocess.run([orca_path, input_file_path], check=True, stdout=output_file, stderr=output_file)
            logging.info(f"ORCA run successful for '{input_file_path}', output saved to '{output_log}'")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running ORCA: {e}")
            sys.exit(1)

def run_orca_plot(homo_index: int) -> None:
    """Generate molecular orbital plots around HOMO using orca_plot.

    Creates orbital plots for orbitals from HOMO-10 to HOMO+10
    using ORCA's orca_plot utility with automated input.

    Args:
        homo_index: Index of the HOMO orbital
    """
    for index in range(homo_index - 10, homo_index + 11):
        process = subprocess.Popen(['orca_plot', 'input.gbw', '-i'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        inputs = f'1\n1\n4\n100\n5\n7\n2\n{index}\n10\n11\n'.encode()
        stdout, stderr = process.communicate(input=inputs)
        if process.returncode == 0:
            logging.info(f"orca_plot ran successfully for index {index}")
        else:
            logging.error(f"orca_plot encountered an error for index {index}: {stderr.decode()}")
