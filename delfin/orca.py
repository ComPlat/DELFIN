import os
import subprocess
import sys
from pathlib import Path
from shutil import which
from typing import Iterable, Optional

from delfin.common.logging import get_logger

logger = get_logger(__name__)

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

def _validate_candidate(candidate: str) -> Optional[str]:
    """Return a usable executable path when candidate points to a file."""
    if not candidate:
        return None

    expanded = Path(candidate.strip()).expanduser()
    if not expanded.is_file():
        return None

    if not os.access(expanded, os.X_OK):
        return None

    return str(expanded.resolve())


def _iter_orca_candidates() -> Iterable[str]:
    """Yield potential ORCA paths from environment and helper tools."""
    env_keys = ("ORCA_BINARY", "ORCA_PATH")
    for key in env_keys:
        value = os.environ.get(key)
        if value:
            yield value

    which_targets = ["orca"]
    if sys.platform.startswith("win"):
        which_targets.append("orca.exe")

    for target in which_targets:
        located = which(target)
        if located:
            yield located

    locator = which("orca_locate")
    if locator:
        try:
            result = subprocess.run([locator], check=False, capture_output=True, text=True)
        except Exception as exc:
            logger.debug(f"Failed to query orca_locate: {exc}")
        else:
            if result.returncode != 0:
                logger.debug(
                    "orca_locate returned non-zero exit status %s with stderr: %s",
                    result.returncode,
                    result.stderr.strip(),
                )
            else:
                for line in result.stdout.splitlines():
                    stripped = line.strip()
                    if stripped:
                        yield stripped


def find_orca_executable() -> Optional[str]:
    """Locate a valid ORCA executable by validating several candidate sources."""
    for candidate in _iter_orca_candidates():
        valid_path = _validate_candidate(candidate)
        if valid_path:
            return valid_path

        logger.debug(f"Discarding invalid ORCA candidate path: {candidate!r}")

    logger.error("ORCA executable not found. Please ensure ORCA is installed and in your PATH.")
    return None


def _run_orca_subprocess(orca_path: str, input_file_path: str, output_log: str, timeout: Optional[int] = None) -> bool:
    """Run ORCA subprocess and capture output. Returns True when successful."""
    with open(output_log, "w") as output_file:
        try:
            subprocess.run([orca_path, input_file_path], check=True, stdout=output_file, stderr=output_file, timeout=timeout)
        except subprocess.TimeoutExpired:
            logger.error(f"ORCA calculation timed out after {timeout} seconds for '{input_file_path}'")
            return False
        except subprocess.CalledProcessError as error:
            logger.error(f"Error running ORCA: {error}")
            return False
    return True

def run_orca(input_file_path: str, output_log: str, timeout: Optional[int] = None) -> None:
    """Execute ORCA calculation with specified input file.

    Runs ORCA subprocess with input file and captures output to log file.
    Logs success/failure and handles subprocess errors.

    Args:
        input_file_path: Path to ORCA input file (.inp)
        output_log: Path for ORCA output file (.out)
        timeout: Optional timeout in seconds for ORCA calculation
    """
    orca_path = find_orca_executable()
    if not orca_path:
        return

    if _run_orca_subprocess(orca_path, input_file_path, output_log, timeout):
        logger.info(f"ORCA run successful for '{input_file_path}'")

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
        logger.error("Cannot run ORCA IMAG calculation because the ORCA executable was not found in PATH.")
        sys.exit(1)

    output_log = f"output_{iteration}.out"
    if _run_orca_subprocess(orca_path, input_file_path, output_log):
        logger.info(f"ORCA run successful for '{input_file_path}', output saved to '{output_log}'")
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
            logger.info(f"orca_plot ran successfully for index {index}")
        else:
            logger.error(f"orca_plot encountered an error for index {index}: {stderr_output}")


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
