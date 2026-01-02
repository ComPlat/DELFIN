import os
import re
import signal
import subprocess
import sys
import threading
import time
import socket
import shutil
from pathlib import Path
from shutil import which
import tempfile
from typing import Dict, Iterable, Optional

try:
    import psutil
except ImportError:
    psutil = None  # type: ignore

from delfin.common.logging import get_logger
from delfin.common.progress_display import ProgressDisplay
from delfin.global_manager import get_global_manager
from delfin.orca_recovery import (
    OrcaErrorDetector,
    OrcaErrorType,
    OrcaInputModifier,
    RecoveryStrategy,
    RetryStateTracker,
)

logger = get_logger(__name__)

# Compile regex patterns at module level for performance
_CPHF_PATTERN = re.compile(r'ITERATION\s+(\d+):.*?\(\s*[\d.]+\s+sec\s+(\d+)/(\d+)\s+done\)')
_FREQ_PATTERN = re.compile(r'(\d+)\s+\(of\s+(\d+)\)')

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
        base_path = Path(tempfile.gettempdir()).joinpath("delfin_orca_scratch")

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


def get_orca_scratch_dir() -> Optional[Path]:
    """Return the active ORCA scratch directory for this run, if any."""
    return _RUN_SCRATCH_DIR


def cleanup_orca_scratch_dir() -> Optional[Path]:
    """Remove the active ORCA scratch directory created for this run."""
    global _RUN_SCRATCH_DIR
    target = _RUN_SCRATCH_DIR
    _RUN_SCRATCH_DIR = None
    if target is None:
        return None

    try:
        if target.exists():
            shutil.rmtree(target)
            logger.info("Removed ORCA scratch directory %s", target)
    except Exception:
        logger.warning("Failed to remove ORCA scratch directory %s", target, exc_info=True)
    return target


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

    # Get manager early so we can check shutdown status
    try:
        manager = get_global_manager()
    except Exception:
        logger.debug("Failed to get global manager", exc_info=True)

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

            # Register subprocess for signal-based cleanup
            if manager:
                try:
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

            # Wait for completion (check for shutdown request periodically)
            # This allows Ctrl+C to interrupt even when waiting for ORCA
            shutdown_check_interval = 0.5  # Check every 500ms
            remaining_timeout = timeout if timeout else float('inf')
            start_time = time.time()

            while True:
                try:
                    # Check if shutdown was requested (e.g., Ctrl+C)
                    if manager and manager._shutdown_requested.is_set():
                        logger.info("Shutdown requested - terminating ORCA process")
                        _kill_process_group(process)
                        raise KeyboardInterrupt("Shutdown requested")

                    # Wait for process with short timeout
                    wait_timeout = min(shutdown_check_interval, remaining_timeout) if timeout else shutdown_check_interval
                    return_code = process.wait(timeout=wait_timeout)
                    break  # Process finished
                except subprocess.TimeoutExpired:
                    # Process still running, check if we should keep waiting
                    if timeout:
                        elapsed = time.time() - start_time
                        remaining_timeout = timeout - elapsed
                        if remaining_timeout <= 0:
                            # Overall timeout expired
                            raise
                    # Otherwise, loop and check again
                    continue

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

            # BUGFIX: Ensure all MPI child processes are terminated
            # ORCA sometimes leaves MPI workers running even after main process exits
            # Give them a short grace period, then force cleanup
            time.sleep(0.5)  # Brief grace period for clean MPI shutdown
            try:
                _ensure_process_group_terminated(process, grace_timeout=2.0)
            except Exception as e:
                logger.debug(f"Process cleanup after successful ORCA run: {e}")

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
    """Monitor ORCA output file and display progress updates.

    Runs in a background thread and displays/logs interesting progress markers:
    - CPHF/POPLE solver iterations (numerical frequencies)
    - Numerical frequency displacements
    - Convergence markers

    Automatically detects TTY and uses live updates when available,
    falls back to logging in batch environments (SLURM, PBS, etc.).

    Args:
        output_file: Path to ORCA output file to monitor
        stop_event: Threading event to signal monitoring should stop
        input_name: Name of input file (for display purposes)
    """
    # Initialize progress display (auto-detects TTY)
    progress = ProgressDisplay()

    last_size = 0
    last_update_time = time.time()
    update_interval = 5  # Update every 5 seconds (reduced from 60s)

    # Track last displayed progress to avoid redundant updates
    last_progress_state = None

    while not stop_event.is_set():
        try:
            # Wait for output file to exist
            if not os.path.exists(output_file):
                time.sleep(2)
                continue

            current_size = os.path.getsize(output_file)
            current_time = time.time()

            # Check if file is growing
            if current_size > last_size:
                # Only check for updates if enough time has passed
                if (current_time - last_update_time) >= update_interval:
                    try:
                        with open(output_file, 'r', errors='replace') as f:
                            # Read last 5000 bytes for recent activity
                            f.seek(max(0, current_size - 5000))
                            recent_content = f.read()
                            recent_lines = recent_content.split('\n')[-50:]

                            progress_message = None
                            progress_key = None  # To detect duplicates
                            is_persistent = False

                            # Check patterns in priority order (most specific first)
                            # Process lines in reverse to get most recent match
                            for line in reversed(recent_lines):

                                # Priority 1: CPHF/POPLE solver progress (numerical frequencies)
                                # Example: "ITERATION   5:  (  3.47 sec  162/162 done)"
                                if 'ITERATION' in line and 'done)' in line:
                                    match = _CPHF_PATTERN.search(line)
                                    if match:
                                        iteration, current, total = match.groups()
                                        percent = (int(current) / int(total)) * 100
                                        progress_message = (
                                            f"[{input_name}] CPHF/POPLE iteration {iteration}: "
                                            f"{current}/{total} ({percent:.1f}%)"
                                        )
                                        progress_key = f"cphf_{iteration}_{current}"
                                        break

                                # Priority 2: Frequency displacement progress
                                # Example: "Calculating gradient on displaced geometry 24 (of 162)"
                                elif 'displaced geometry' in line:
                                    match = _FREQ_PATTERN.search(line)
                                    if match:
                                        current, total = match.groups()
                                        percent = (int(current) / int(total)) * 100
                                        progress_message = (
                                            f"[{input_name}] Frequency calculation: "
                                            f"displaced geometry {current}/{total} ({percent:.1f}%)"
                                        )
                                        progress_key = f"freq_{current}"
                                        break

                                # Priority 3: Optimization convergence (persistent message)
                                elif 'THE OPTIMIZATION HAS CONVERGED' in line:
                                    progress_message = f"[{input_name}] Geometry optimization converged"
                                    progress_key = "opt_converged"
                                    is_persistent = True
                                    break

                                # Priority 4: SCF convergence (persistent message)
                                elif 'SCF CONVERGED' in line:
                                    progress_message = f"[{input_name}] SCF converged"
                                    progress_key = "scf_converged"
                                    is_persistent = True
                                    break

                            # Only update display if we found NEW progress (different from last state)
                            # This ensures no output when nothing has changed
                            if progress_message and progress_key and progress_key != last_progress_state:
                                progress.update(progress_message, persistent=is_persistent)
                                last_update_time = current_time
                                last_progress_state = progress_key

                    except OSError as e:
                        logger.debug(f"Could not read ORCA output file {output_file}: {e}")
                    except Exception as e:
                        logger.debug(f"Error processing ORCA output: {e}")

                # Always update last_size when file grows
                last_size = current_size

            # Sleep between checks
            time.sleep(5)

        except Exception as e:
            logger.debug(f"Progress monitor error: {e}")
            time.sleep(5)

    # Finalize progress display when stopping
    # This ensures clean output (newline in TTY mode)
    progress.finalize()


def _ensure_process_group_terminated(process: subprocess.Popen, grace_timeout: float = 2.0) -> None:
    """Ensure all processes in the process group are terminated, even orphaned children.

    ORCA sometimes leaves MPI worker processes running even after the main process exits.
    This function checks for such orphans and cleans them up.
    """
    try:
        pgid = os.getpgid(process.pid)
    except (ProcessLookupError, OSError):
        # Process group already gone
        return

    if psutil is None:
        # psutil not available, use fallback
        logger.debug("psutil not available, using process group kill fallback")
        try:
            os.killpg(pgid, signal.SIGTERM)
            time.sleep(grace_timeout)
            os.killpg(pgid, signal.SIGKILL)
        except (ProcessLookupError, OSError):
            pass
        return

    try:
        # Find all processes in this process group
        orphans = []
        for proc in psutil.process_iter(['pid', 'name', 'status']):
            try:
                if os.getpgid(proc.pid) == pgid and proc.pid != process.pid:
                    orphans.append(proc)
            except (psutil.NoSuchProcess, ProcessLookupError, OSError):
                continue

        if not orphans:
            return  # All clean

        logger.debug(f"Found {len(orphans)} orphaned processes in group {pgid}")

        # Try graceful termination first
        for proc in orphans:
            try:
                proc.terminate()
            except (psutil.NoSuchProcess, ProcessLookupError):
                pass

        # Wait for grace period
        gone, alive = psutil.wait_procs(orphans, timeout=grace_timeout)

        # Force kill any survivors
        for proc in alive:
            try:
                logger.warning(f"Force killing orphaned ORCA process {proc.pid} ({proc.name()})")
                proc.kill()
            except (psutil.NoSuchProcess, ProcessLookupError):
                pass

    except Exception as e:
        # If psutil fails, fall back to process group kill
        logger.debug(f"Orphan cleanup via psutil failed: {e}, using fallback")
        try:
            os.killpg(pgid, signal.SIGTERM)
            import time
            time.sleep(grace_timeout)
            os.killpg(pgid, signal.SIGKILL)
        except (ProcessLookupError, OSError):
            pass


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




def run_orca_with_intelligent_recovery(
    input_file_path: str,
    output_log: str,
    timeout: Optional[int] = None,
    *,
    scratch_subdir: Optional[Path] = None,
    working_dir: Optional[Path] = None,
    config: Optional[Dict] = None,
) -> bool:
    """Execute ORCA with intelligent error detection and automatic recovery.

    This function combines transient error retry logic with intelligent
    error-specific recovery strategies. It detects specific ORCA failures
    (SCF convergence, TRAH crashes, etc.) and automatically modifies the
    input file with appropriate fixes, using MOREAD to continue from the
    last successful state.

    Recovery features:
    - Automatic error classification (SCF, TRAH, geometry, MPI, etc.)
    - Progressive escalation of fixes across retry attempts
    - MOREAD-based continuation from last .gbw state
    - Preservation of basis sets and geometry
    - State tracking to prevent infinite loops

    Args:
        input_file_path: Path to ORCA input file (.inp)
        output_log: Path for ORCA output file (.out)
        timeout: Optional timeout in seconds
        scratch_subdir: Optional scratch subdirectory
        working_dir: Optional working directory
        config: Configuration dict with recovery settings

    Returns:
        bool: True if ORCA completed successfully, False otherwise

    Configuration options (in CONTROL.txt):
        enable_auto_recovery: Enable intelligent error recovery (yes/no, default: yes)
        max_recovery_attempts: Maximum recovery attempts per error type (default: 2)
    """
    config = config or {}

    # Check if intelligent recovery is enabled
    # Support both new and old config parameter names for backward compatibility
    recovery_enabled = (
        _parse_bool_config(config.get('enable_auto_recovery', 'no'))
        or _parse_bool_config(config.get('orca_retry_enabled', 'no'))  # Old parameter
    )

    manager = get_global_manager()

    def _shutdown_requested() -> bool:
        try:
            return bool(getattr(manager, "_shutdown_requested", None) and manager._shutdown_requested.is_set())
        except Exception:
            return False

    if not recovery_enabled:
        # Fall back to direct ORCA execution (no retry, no recovery)
        return run_orca(
            input_file_path,
            output_log,
            timeout,
            scratch_subdir=scratch_subdir,
            working_dir=working_dir,
        )

    # Intelligent recovery is enabled
    inp_path = Path(input_file_path)
    out_path = Path(output_log)
    work_dir = working_dir or inp_path.parent

    # Initialize recovery components
    detector = OrcaErrorDetector()
    state_file = work_dir / ".delfin_recovery_state.json"
    tracker = RetryStateTracker(state_file)

    # Support both new and old config parameter names
    max_recovery_attempts = int(
        config.get('max_recovery_attempts')
        or config.get('orca_retry_max_attempts', 3)  # Old parameter
    )

    job_name = inp_path.stem
    current_inp = inp_path

    logger.info(f"Starting ORCA with intelligent recovery enabled for {job_name}")

    # Track all attempted error types to prevent loops
    attempted_errors = set()

    for overall_attempt in range(1, max_recovery_attempts + 2):  # +1 for initial attempt
        if _shutdown_requested():
            logger.warning("Shutdown requested; aborting recovery for %s", job_name)
            raise KeyboardInterrupt

        # Run ORCA
        success = run_orca(
            str(current_inp),
            output_log,
            timeout,
            scratch_subdir=scratch_subdir,
            working_dir=working_dir,
        )

        if success:
            if overall_attempt > 1:
                logger.info(
                    f"✓ ORCA succeeded after {overall_attempt - 1} recovery attempt(s) for {job_name}"
                )
            return True

        # Job failed - analyze error
        if _shutdown_requested():
            logger.warning("Shutdown requested after failure; aborting recovery for %s", job_name)
            raise KeyboardInterrupt

        error_type = detector.analyze_output(out_path)

        if not error_type or error_type == OrcaErrorType.UNKNOWN:
            logger.error(f"ORCA failed with unrecoverable error for {job_name}")
            return False

        # Check if we've already tried to recover from this error type
        error_key = f"{job_name}_{error_type.value}"
        if error_key in attempted_errors:
            logger.error(
                f"Error type {error_type.value} persists after recovery attempt for {job_name}"
            )
            return False

        # Check if we should attempt recovery
        if not tracker.should_retry(job_name, error_type, max_recovery_attempts):
            logger.error(
                f"Max recovery attempts ({max_recovery_attempts}) reached "
                f"for {error_type.value} in {job_name}"
            )
            return False

        # Get recovery strategy
        attempt = tracker.get_attempt(job_name, error_type) + 1
        strategy = RecoveryStrategy(error_type, attempt, config)

        # Get modifications to check for backoff delay
        mods = strategy.get_modifications()

        logger.warning(
            f"⚠ ORCA failed with {error_type.value} for {job_name}. "
            f"Applying recovery strategy (attempt {attempt}/{max_recovery_attempts})..."
        )

        # Handle exponential backoff for transient errors
        if "backoff_delay" in mods:
            delay = mods["backoff_delay"]
            logger.warning(
                f"Transient system error detected. Waiting {delay}s before retry (exponential backoff)..."
            )
            time.sleep(delay)

        # Modify input file
        modifier = OrcaInputModifier(current_inp, config)
        new_inp = modifier.apply_recovery(strategy)

        if new_inp == current_inp:
            logger.error(f"Failed to create recovery input for {job_name}")
            return False

        # Prepare the new input for continuation (update geometry from xyz, ensure GBW backup)
        # This is already done by OrcaInputModifier, but we ensure it's applied
        # Note: The modifier already calls _update_geometry_from_xyz and _add_moread internally

        # Update tracking
        tracker.increment_attempt(job_name, error_type)
        attempted_errors.add(error_key)

        # Backup the failed output file before retry (for debugging)
        if out_path.exists():
            backup_num = attempt
            backup_path = out_path.with_suffix(f'.old{backup_num}.out')
            try:
                shutil.copy2(out_path, backup_path)
                logger.info(f"Backed up failed output: {backup_path.name}")
            except Exception as e:
                logger.warning(f"Failed to backup output file: {e}")

        # Use modified input for next attempt
        current_inp = new_inp

        logger.info(f"Retrying with recovery input: {new_inp.name}")

    # All recovery attempts exhausted
    logger.error(
        f"ORCA failed after {max_recovery_attempts} recovery attempts for {job_name}"
    )
    return False




def _parse_bool_config(value) -> bool:
    """Parse boolean value from config."""
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.lower() in ('yes', 'true', '1', 'on')
    return bool(value)


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
