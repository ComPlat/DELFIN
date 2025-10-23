"""Process conflict detection for ORCA jobs."""

import os
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple

from delfin.common.logging import get_logger

logger = get_logger(__name__)


def find_competing_orca_processes(inp_file: str) -> List[Tuple[int, str]]:
    """Find ORCA processes working on the same input/gbw files.

    Returns list of (pid, cmdline) tuples for competing processes.
    """
    inp_path = Path(inp_file).resolve()
    base_name = inp_path.stem  # e.g., "input2" from "input2.inp"
    work_dir = inp_path.parent

    competing = []

    try:
        # Find all ORCA-related processes
        result = subprocess.run(
            ["ps", "aux"],
            capture_output=True,
            text=True,
            timeout=5
        )

        for line in result.stdout.splitlines():
            # Look for orca, orca_leanscf_mpi, mpirun processes
            if not any(keyword in line for keyword in ['orca', 'mpirun']):
                continue

            # Check if this process is working on our files
            if base_name in line:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        pid = int(parts[1])
                        cmdline = ' '.join(parts[10:])

                        # Don't report our own process
                        if pid != os.getpid():
                            competing.append((pid, cmdline))
                    except (ValueError, IndexError):
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
