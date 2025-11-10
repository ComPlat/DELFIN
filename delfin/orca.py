import os
import signal
import subprocess
import sys
import threading
import time
import socket
from pathlib import Path
from shutil import which
from typing import Dict, Iterable, Optional

from delfin.common.logging import get_logger
from delfin.common.paths import get_runtime_dir
from delfin.global_manager import get_global_manager

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

_RUN_SCRATCH_DIR: Optional[Path] = None


def _ensure_openmpi_subdir(base_path: Path) -> None:
    """Create the OpenMPI-specific scratch subdirectory expected by ORCA."""
    if base_path is None:
        return
    try:
        hostname = os.environ.get("ORCA_MPI_HOST") or socket.gethostname()
    except Exception:
        hostname = "localhost"
    uid = None
    if hasattr(os, "getuid"):
        try:
            uid = os.getuid()
        except Exception:
            uid = None
    if uid is None:
        uid = os.environ.get("ORCA_MPI_UID") or os.environ.get("UID")
    if uid is None:
        uid = os.getpid()
    subdir = base_path / f"ompi.{hostname}.{uid}"
    try:
        subdir.mkdir(parents=True, exist_ok=True)
    except Exception:
        logger.debug("Could not create OpenMPI scratch subdir %s", subdir, exc_info=True)


def _ensure_orca_scratch_dir() -> Path:
    """Create (once) and return a run-specific scratch directory for ORCA."""
    global _RUN_SCRATCH_DIR
    if _RUN_SCRATCH_DIR is not None:
        return _RUN_SCRATCH_DIR

    base_candidates = [
        os.environ.get("ORCA_SCRDIR"),
        os.environ.get("ORCA_TMPDIR"),
        os.environ.get("DELFIN_SCRATCH"),
    ]

    for candidate in base_candidates:
        if candidate:
            base_path = Path(candidate).expanduser()
            break
    else:
        base_path = get_runtime_dir().joinpath(".orca_scratch")

    base_path.mkdir(parents=True, exist_ok=True)

    run_label = os.environ.get("DELFIN_RUN_TOKEN")
    if not run_label:
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        cwd_name = Path.cwd().name or "delfin"
        run_label = f"{cwd_name}_{os.getpid()}_{timestamp}"

    scratch_dir = base_path / run_label
    scratch_dir.mkdir(parents=True, exist_ok=True)
    _ensure_openmpi_subdir(scratch_dir)
    _RUN_SCRATCH_DIR = scratch_dir
    return scratch_dir


def _prepare_orca_environment(extra_scratch: Optional[Path] = None) -> Dict[str, str]:
    """Return a subprocess environment with isolated ORCA scratch settings.

    Args:
        extra_scratch: Optional subdirectory appended to the run scratch path.
    """
    env = os.environ.copy()
    scratch_dir = _ensure_orca_scratch_dir()
    if extra_scratch is not None:
        scratch_dir = scratch_dir / extra_scratch
        scratch_dir.mkdir(parents=True, exist_ok=True)
        _ensure_openmpi_subdir(scratch_dir)

    scratch_str = str(scratch_dir)
    env["ORCA_SCRDIR"] = scratch_str
    env["ORCA_TMPDIR"] = scratch_str
    env.setdefault("TMPDIR", scratch_str)
    env.setdefault("DELFIN_RUN_TOKEN", scratch_dir.name)
    return env

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


def _run_orca_subprocess(
    orca_path: str,
    input_file_path: str,
    output_log: str,
    timeout: Optional[int] = None,
    scratch_subdir: Optional[Path] = None,
    working_dir: Optional[Path] = None,
) -> bool:
    """Run ORCA subprocess and capture output. Returns True when successful."""
    process = None
    monitor_thread = None
    stop_event = threading.Event()
    manager = None
    registration_token: Optional[str] = None

    try:
        with open(output_log, "w") as output_file:
            # Use Popen with process group to ensure all child processes can be killed
            # start_new_session creates a new process group, making cleanup easier
            process = subprocess.Popen(
                [orca_path, input_file_path],
                stdout=output_file,
                stderr=output_file,
                env=_prepare_orca_environment(scratch_subdir),
                start_new_session=True,  # Create new process group
                cwd=str(working_dir) if working_dir is not None else None,
            )

            try:
                manager = get_global_manager()
                try:
                    cwd_hint = Path(input_file_path).resolve().parent
                except Exception:
                    cwd_hint = Path.cwd()
                registration_token = manager.register_subprocess(
                    process,
                    label=input_file_path,
                    cwd=str(cwd_hint),
                )
            except Exception:
                logger.debug("Failed to register ORCA subprocess for tracking", exc_info=True)

            # Start progress monitoring thread
            input_name = Path(input_file_path).stem
            monitor_thread = threading.Thread(
                target=_monitor_orca_progress,
                args=(output_log, stop_event, input_name),
                daemon=True
            )
            monitor_thread.start()

            # Wait for completion
            return_code = process.wait(timeout=timeout)

            # Stop progress monitor
            stop_event.set()
            if monitor_thread:
                monitor_thread.join(timeout=2)

            if return_code != 0:
                logger.error(f"ORCA failed with return code {return_code} for {input_file_path}")
                logger.error(f"Check {output_log} for details")
                return False

            # Check if ORCA actually terminated normally
            success_marker = _check_orca_success(output_log)
            if not success_marker:
                logger.error(f"ORCA did not terminate normally for {input_file_path}")
                logger.error(f"Check {output_log} for error messages")
                return False

            return True

    except subprocess.TimeoutExpired:
        logger.error(f"ORCA timeout after {timeout}s")
        stop_event.set()
        if process:
            _kill_process_group(process)
        return False
    except KeyboardInterrupt:
        logger.warning("ORCA interrupted by user (Ctrl+C)")
        stop_event.set()
        if process:
            _kill_process_group(process)
        raise
    except Exception as e:
        logger.error(f"ORCA subprocess error: {e}")
        stop_event.set()
        if process:
            _kill_process_group(process)
        return False
    finally:
        # Ensure monitor thread is stopped
        stop_event.set()
        if monitor_thread and monitor_thread.is_alive():
            monitor_thread.join(timeout=1)
        if manager and registration_token:
            try:
                manager.unregister_subprocess(registration_token)
            except Exception:
                logger.debug("Failed to unregister ORCA subprocess %s", registration_token, exc_info=True)


def _check_orca_success(output_file: str) -> bool:
    """Check if ORCA terminated normally by looking for success marker."""
    try:
        with open(output_file, 'r') as f:
            content = f.read()
            return 'ORCA TERMINATED NORMALLY' in content
    except Exception as e:
        logger.debug(f"Could not check ORCA success marker: {e}")
        return False


def _monitor_orca_progress(output_file: str, stop_event: threading.Event, input_name: str):
    """Monitor ORCA output file and log progress updates.

    Runs in a background thread and logs interesting progress markers:
    - SCF iterations
    - Geometry optimization steps
    - Numerical frequency displacements
    """
    last_size = 0
    last_log_time = time.time()
    log_interval = 60  # Log every 60 seconds

    # Patterns to detect in output
    patterns = {
        'scf': r'ITER\s+Energy\s+Delta-E',
        'opt': r'OPTIMIZATION\s+RUN',
        'freq_disp': r'Calculating gradient on displaced geometry\s+(\d+)\s+\(of\s+(\d+)\)',
        'terminating': r'ORCA TERMINATED',
    }

    while not stop_event.is_set():
        try:
            if not os.path.exists(output_file):
                time.sleep(2)
                continue

            current_size = os.path.getsize(output_file)
            current_time = time.time()

            # Check if file is growing and enough time has passed
            if current_size > last_size and (current_time - last_log_time) >= log_interval:
                with open(output_file, 'r') as f:
                    # Read last 50 lines for recent activity
                    f.seek(max(0, current_size - 5000))
                    recent_lines = f.read().split('\n')[-50:]

                    for line in recent_lines:
                        # Check for frequency displacement progress
                        if 'displaced geometry' in line:
                            import re
                            match = re.search(r'(\d+)\s+\(of\s+(\d+)\)', line)
                            if match:
                                current, total = match.groups()
                                percent = (int(current) / int(total)) * 100
                                logger.info(f"[{input_name}] Frequency calculation: {current}/{total} ({percent:.1f}%)")
                                last_log_time = current_time
                                break

                        # Check for geometry optimization
                        elif 'GEOMETRY OPTIMIZATION CYCLE' in line:
                            logger.info(f"[{input_name}] Geometry optimization running...")
                            last_log_time = current_time
                            break

                        # Check for SCF convergence
                        elif 'SCF CONVERGED' in line:
                            logger.info(f"[{input_name}] SCF converged")
                            last_log_time = current_time
                            break

                last_size = current_size

            time.sleep(5)  # Check every 5 seconds

        except Exception as e:
            logger.debug(f"Progress monitor error: {e}")
            time.sleep(5)


def _kill_process_group(process: subprocess.Popen) -> None:
    """Kill entire process group including all child processes (like mpirun)."""
    if process.poll() is None:  # Process still running
        try:
            # Send SIGTERM to entire process group
            os.killpg(os.getpgid(process.pid), signal.SIGTERM)
            logger.info(f"Sent SIGTERM to process group {os.getpgid(process.pid)}")

            # Wait a bit for graceful shutdown
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                # Force kill if still running
                os.killpg(os.getpgid(process.pid), signal.SIGKILL)
                logger.warning(f"Sent SIGKILL to process group {os.getpgid(process.pid)}")
                process.wait()
        except ProcessLookupError:
            # Process already terminated
            pass
        except Exception as e:
            logger.error(f"Error killing process group: {e}")

def run_orca(
    input_file_path: str,
    output_log: str,
    timeout: Optional[int] = None,
    *,
    scratch_subdir: Optional[Path] = None,
    working_dir: Optional[Path] = None,
) -> bool:
    """Execute ORCA calculation with specified input file.

    Runs ORCA subprocess with input file and captures output to log file.
    Logs success/failure and handles subprocess errors.

    Args:
        input_file_path: Path to ORCA input file (.inp)
        output_log: Path for ORCA output file (.out)
        timeout: Optional timeout in seconds for ORCA calculation

    Returns:
        bool: True if ORCA completed successfully, False otherwise
    """
    orca_path = find_orca_executable()
    if not orca_path:
        return False

    input_path = Path(input_file_path)
    output_path = Path(output_log)

    if working_dir is not None:
        working_dir = Path(working_dir)
        if not input_path.is_absolute():
            input_path = (Path.cwd() / input_path).resolve()
        if not output_path.is_absolute():
            output_path = (Path.cwd() / output_path).resolve()
    else:
        if not input_path.is_absolute():
            input_path = input_path.resolve()
        if not output_path.is_absolute():
            output_path = output_path.resolve()

    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
    except Exception:
        # Directory might already exist or creation could fail due to permissions;
        # defer handling to the subprocess call / open().
        pass

    if _run_orca_subprocess(
        orca_path,
        str(input_path),
        str(output_path),
        timeout,
        scratch_subdir=scratch_subdir,
        working_dir=working_dir,
    ):
        logger.info(f"ORCA run successful for '{input_file_path}'")
        return True
    return False

def run_orca_IMAG(input_file_path: str, iteration: int, *, working_dir: Optional[Path] = None) -> bool:
    """Execute ORCA calculation for imaginary frequency workflow.

    Specialized ORCA runner for IMAG workflow with iteration-specific
    output naming and enhanced error handling.

    Args:
        input_file_path: Path to ORCA input file
        iteration: Iteration number for output file naming
        working_dir: Directory in which ORCA should be executed
    """
    orca_path = find_orca_executable()
    if not orca_path:
        logger.error("Cannot run ORCA IMAG calculation because the ORCA executable was not found in PATH.")
        sys.exit(1)

    input_path = Path(input_file_path)
    if working_dir is not None:
        working_dir = Path(working_dir)
        output_log_path = working_dir / f"output_{iteration}.out"
        if not input_path.is_absolute():
            # Provide ORCA with an absolute path when running inside working_dir
            input_path = (Path.cwd() / input_path).resolve()
    else:
        output_log_path = Path(f"output_{iteration}.out")
        if not input_path.is_absolute():
            input_path = input_path.resolve()

    if _run_orca_subprocess(
        orca_path,
        str(input_path),
        str(output_log_path),
        working_dir=working_dir,
    ):
        logger.info(f"ORCA run successful for '{input_file_path}', output saved to '{output_log_path}'")
        return True

    logger.error(f"ORCA IMAG calculation failed for '{input_file_path}'. See '{output_log_path}' for details.")
    return False

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
