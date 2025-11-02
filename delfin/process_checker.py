"""Process conflict detection for ORCA jobs."""

import os
from pathlib import Path
from typing import List, Tuple

from delfin.common.logging import get_logger

logger = get_logger(__name__)


def find_competing_orca_processes(inp_file: str) -> List[Tuple[int, str]]:
    """Find ORCA processes working on the same input/gbw files.

    Returns list of (pid, cmdline) tuples for competing processes.
    """
    inp_path = Path(inp_file).resolve()
    base_name = inp_path.stem  # e.g., "input2" from "input2.inp"
    work_dir = inp_path.parent
    try:
        work_dir = work_dir.resolve()
    except OSError:
        work_dir = work_dir

    competing: List[Tuple[int, str]] = []

    try:
        # Find all ORCA-related processes
        for pid_dir in Path("/proc").iterdir():
            if not pid_dir.name.isdigit():
                continue

            try:
                pid = int(pid_dir.name)
                if pid == os.getpid():
                    continue  # skip our own process

                exe_path = pid_dir / "exe"
                try:
                    target = os.readlink(exe_path)
                except (OSError, RuntimeError):
                    continue

                if not any(keyword in target for keyword in ["orca", "mpirun"]):
                    continue

                try:
                    raw_cmd = (pid_dir / "cmdline").read_bytes()
                    if not raw_cmd:
                        continue
                    tokens = [tok for tok in raw_cmd.split(b"\0") if tok]
                    cmdline = " ".join(tok.decode(errors="ignore") for tok in tokens)
                except (OSError, RuntimeError, UnicodeDecodeError):
                    tokens = []
                    cmdline = target

                # Require exact match of both basename and working directory
                if not any(Path(tok.decode(errors="ignore")).stem == base_name for tok in tokens if tok):
                    continue

                try:
                    proc_cwd = Path(os.readlink(pid_dir / "cwd")).resolve()
                except (OSError, RuntimeError):
                    continue

                if proc_cwd != work_dir:
                    continue

                competing.append((pid, cmdline))
            except ValueError:
                continue

    except Exception as e:
        logger.debug(f"Error checking for competing processes: {e}")

    return competing


def check_and_warn_competing_processes(inp_file: str) -> bool:
    """Check for competing ORCA processes and warn if found.

    Returns True if competing processes were found, False otherwise.
    """
    competing = find_competing_orca_processes(inp_file)

    if competing:
        logger.warning(
            f"Found {len(competing)} competing ORCA process(es) for {Path(inp_file).name}:"
        )
        for pid, cmdline in competing[:5]:  # Show max 5
            logger.warning(f"  PID {pid}: {cmdline[:100]}")

        logger.warning(
            "These processes may block ORCA execution. "
            "Consider stopping them or waiting for completion."
        )
        return True

    return False


def wait_for_process_cleanup(inp_file: str, max_checks: int = 3, interval: int = 2) -> bool:
    """Wait for competing processes to finish.

    Args:
        inp_file: Input file to check
        max_checks: Maximum number of checks
        interval: Seconds between checks

    Returns:
        True if processes cleared, False if still running
    """
    import time

    for i in range(max_checks):
        competing = find_competing_orca_processes(inp_file)
        if not competing:
            return True

        if i < max_checks - 1:
            logger.info(f"Waiting {interval}s for {len(competing)} competing process(es) to finish...")
            time.sleep(interval)

    return False
