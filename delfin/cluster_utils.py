"""Cluster utilities for automatic resource detection and configuration."""

import os
import re
import subprocess
import time
from typing import Optional, Dict, Any
from delfin.common.logging import get_logger

logger = get_logger(__name__)

try:
    import psutil  # type: ignore
except ImportError:
    psutil = None  # type: ignore


def detect_cluster_environment() -> Dict[str, Any]:
    """Detect cluster scheduler and available resources.

    Returns:
        Dictionary with cluster information and resource limits
    """
    cluster_info = {
        'scheduler': None,
        'job_id': None,
        'cpus_available': None,
        'memory_available_mb': None,
        'walltime_remaining': None
    }

    # Detect SLURM
    if os.getenv('SLURM_JOB_ID'):
        cluster_info.update(_detect_slurm_resources())

    # Detect PBS/Torque
    elif os.getenv('PBS_JOBID'):
        cluster_info.update(_detect_pbs_resources())

    # Detect LSF
    elif os.getenv('LSB_JOBID'):
        cluster_info.update(_detect_lsf_resources())

    # Fallback to system detection
    else:
        cluster_info.update(_detect_system_resources())

    return cluster_info


def _detect_slurm_resources() -> Dict[str, Any]:
    """Detect SLURM-specific resources."""
    return {
        'scheduler': 'SLURM',
        'job_id': os.getenv('SLURM_JOB_ID'),
        'cpus_available': _get_int_env('SLURM_CPUS_PER_TASK') or _get_int_env('SLURM_NPROCS'),
        'memory_available_mb': _get_slurm_memory(),
        'partition': os.getenv('SLURM_JOB_PARTITION'),
        'node_list': os.getenv('SLURM_JOB_NODELIST'),
    }


def _detect_pbs_resources() -> Dict[str, Any]:
    """Detect PBS/Torque-specific resources."""
    return {
        'scheduler': 'PBS',
        'job_id': os.getenv('PBS_JOBID'),
        'cpus_available': _get_int_env('PBS_NUM_PPN') or _get_int_env('NCPUS'),
        'memory_available_mb': _get_pbs_memory(),
        'queue': os.getenv('PBS_QUEUE'),
    }


def _detect_lsf_resources() -> Dict[str, Any]:
    """Detect LSF-specific resources."""
    return {
        'scheduler': 'LSF',
        'job_id': os.getenv('LSB_JOBID'),
        'cpus_available': _get_int_env('LSB_DJOB_NUMPROC'),
        'memory_available_mb': _get_lsf_memory(),
        'queue': os.getenv('LSB_QUEUE'),
    }


def _detect_system_resources() -> Dict[str, Any]:
    """Fallback system resource detection."""
    if psutil is None:
        return {
            'scheduler': 'system',
            'cpus_available': os.cpu_count(),
            'memory_available_mb': None,
        }
    return {
        'scheduler': 'system',
        'cpus_available': psutil.cpu_count(),
        'memory_available_mb': psutil.virtual_memory().total // (1024 * 1024),
    }


def _get_int_env(var_name: str) -> Optional[int]:
    """Get integer value from environment variable."""
    value = os.getenv(var_name)
    if value:
        try:
            return int(value)
        except ValueError:
            logger.warning(f"Could not parse {var_name}={value} as integer")
    return None


def _get_slurm_memory() -> Optional[int]:
    """Parse SLURM memory allocation (handles MB/GB suffixes)."""
    mem_per_cpu = os.getenv('SLURM_MEM_PER_CPU')
    mem_per_node = os.getenv('SLURM_MEM_PER_NODE')

    if mem_per_node:
        return _parse_memory_string(mem_per_node)
    elif mem_per_cpu:
        cpus = _get_int_env('SLURM_CPUS_PER_TASK') or 1
        mem_mb = _parse_memory_string(mem_per_cpu)
        return mem_mb * cpus if mem_mb else None

    return None


def _get_pbs_memory() -> Optional[int]:
    """Parse PBS memory allocation."""
    # PBS memory can be in various formats
    mem_env = os.getenv('PBS_VMEM') or os.getenv('PBS_MEM')
    if mem_env:
        return _parse_memory_string(mem_env)
    return None


def _get_lsf_memory() -> Optional[int]:
    """Parse LSF memory allocation."""
    mem_env = os.getenv('LSB_DJOB_RUSAGE')
    if mem_env and 'mem=' in mem_env:
        # Extract memory from rusage string like "mem=8000"
        try:
            mem_part = [part for part in mem_env.split() if part.startswith('mem=')][0]
            mem_value = mem_part.split('=')[1]
            return _parse_memory_string(mem_value)
        except (IndexError, ValueError):
            pass
    return None


def _parse_memory_string(mem_str: str) -> Optional[int]:
    """Parse memory string with units (MB, GB, etc.) to MB."""
    if not mem_str:
        return None

    mem_str = mem_str.strip().upper()

    # Extract numeric part and unit
    import re
    match = re.match(r'(\d+(?:\.\d+)?)\s*([KMGT]?B?)', mem_str)
    if not match:
        try:
            return int(mem_str)  # Assume MB if no unit
        except ValueError:
            return None

    value, unit = match.groups()
    value = float(value)

    # Convert to MB
    if unit in ('', 'B', 'MB', 'M'):
        return int(value)
    elif unit in ('KB', 'K'):
        return int(value / 1024)
    elif unit in ('GB', 'G'):
        return int(value * 1024)
    elif unit in ('TB', 'T'):
        return int(value * 1024 * 1024)

    return int(value)  # Default to MB


def auto_configure_resources(config: Dict[str, Any]) -> Dict[str, Any]:
    """Automatically configure PAL and maxcore based on detected resources.

    Args:
        config: Current DELFIN configuration

    Returns:
        Updated configuration with auto-detected resources
    """
    cluster_info = detect_cluster_environment()
    updated_config = config.copy()

    # Auto-configure PAL (CPU cores)
    if cluster_info['cpus_available'] and not config.get('PAL'):
        suggested_pal = cluster_info['cpus_available']
        updated_config['PAL'] = suggested_pal
        logger.info(f"Auto-detected {cluster_info['cpus_available']} CPUs, setting PAL={suggested_pal}")

    # Auto-configure maxcore (memory per core in MB)
    if cluster_info['memory_available_mb'] and not config.get('maxcore'):
        cpus = updated_config.get('PAL', cluster_info['cpus_available'] or 1)
        # Reserve 1GB for system, distribute rest among cores
        available_mb = cluster_info['memory_available_mb'] - 1024
        if available_mb > 0:
            suggested_maxcore = max(500, min(available_mb // cpus, 4000))  # 500MB-4GB per core
            updated_config['maxcore'] = suggested_maxcore
            logger.info(f"Auto-detected {cluster_info['memory_available_mb']}MB RAM, setting maxcore={suggested_maxcore}MB per core")

    # Log cluster environment
    if cluster_info['scheduler'] != 'system':
        logger.info(f"Running on {cluster_info['scheduler']} cluster (Job ID: {cluster_info['job_id']})")

    return updated_config


def _slurm_remaining_seconds() -> Optional[int]:
    """Remaining walltime of the current SLURM job, or None if unknowable.

    ``squeue -h -o %L`` reports time left directly, in the same D-HH:MM:SS
    family of formats as TimeLimit. It is asked first because it needs no
    arithmetic; ``scontrol`` (TimeLimit minus RunTime) is the fallback for
    configurations where squeue is restricted on compute nodes.
    """
    job_id = os.getenv('SLURM_JOB_ID')
    if not job_id:
        return None

    try:
        result = subprocess.run(
            ['squeue', '-h', '-j', job_id, '-o', '%L'],
            capture_output=True, text=True, timeout=10,
        )
        if result.returncode == 0:
            remaining = _parse_walltime(result.stdout.strip())
            if remaining is not None:
                return remaining
    except (OSError, subprocess.SubprocessError):
        pass

    try:
        result = subprocess.run(
            ['scontrol', 'show', 'job', job_id],
            capture_output=True, text=True, timeout=10,
        )
        if result.returncode == 0:
            fields = dict(
                re.findall(r'\b(TimeLimit|RunTime)=(\S+)', result.stdout)
            )
            limit = _parse_walltime(fields.get('TimeLimit', ''))
            runtime = _parse_walltime(fields.get('RunTime', ''))
            if limit is not None and runtime is not None:
                return max(0, limit - runtime)
    except (OSError, subprocess.SubprocessError):
        pass

    return None


def get_walltime_limit() -> Optional[int]:
    """Get remaining walltime in seconds, or None when there is no limit.

    None means "unknown or unlimited" and every caller must read it that way:
    a wrong number here would make the scheduler refuse jobs that do in fact
    fit, which is worse than not knowing.
    """
    # SLURM
    if os.getenv('SLURM_JOB_ID'):
        remaining = _slurm_remaining_seconds()
        if remaining is not None:
            return remaining
        # The job wrapper stamps the deadline at startup, which still holds
        # when squeue/scontrol are unreachable from the compute node.
        end_epoch = os.getenv('DELFIN_JOB_END_EPOCH')
        if end_epoch:
            try:
                return max(0, int(float(end_epoch)) - int(time.time()))
            except (TypeError, ValueError):
                pass
        return None

    # PBS
    elif os.getenv('PBS_JOBID'):
        walltime = os.getenv('PBS_WALLTIME')
        if walltime:
            return _parse_walltime(walltime)

    # LSF
    elif os.getenv('LSB_JOBID'):
        # LSF typically doesn't expose remaining walltime easily
        pass

    return None


def _parse_walltime(walltime_str: str) -> Optional[int]:
    """Parse a walltime string to seconds, or None when it carries no number.

    Covers what SLURM emits for TimeLimit, RunTime and squeue's %L — that is
    ``D-HH:MM:SS``, ``D-HH:MM``, ``HH:MM:SS`` and ``MM:SS`` — plus the bare
    hour count used by PBS. The placeholders SLURM uses for "no limit"
    (UNLIMITED, NOT_SET, INVALID) map to None rather than to zero, so an
    unlimited job is never mistaken for one that is out of time.
    """
    if not walltime_str:
        return None

    walltime_str = walltime_str.strip()
    if not walltime_str or walltime_str.upper() in {
        'UNLIMITED', 'NOT_SET', 'INVALID', 'N/A', 'NONE'
    }:
        return None

    days = 0
    if '-' in walltime_str:
        day_part, _, walltime_str = walltime_str.partition('-')
        try:
            days = int(day_part)
        except ValueError:
            return None

    # Formats: HH:MM:SS, MM:SS (and D-HH:MM once the day part is split off)
    if ':' in walltime_str:
        parts = walltime_str.split(':')
        try:
            values = [int(p) for p in parts]
        except ValueError:
            return None
        if len(values) == 3:
            hours, minutes, seconds = values
        elif len(values) == 2:
            # With a day prefix SLURM writes D-HH:MM; without one it is MM:SS.
            hours, minutes, seconds = (values[0], values[1], 0) if days else (0, values[0], values[1])
        else:
            return None
        return days * 86400 + hours * 3600 + minutes * 60 + seconds

    # Format: hours (integer or float)
    try:
        hours = float(walltime_str)
        return days * 86400 + int(hours * 3600)
    except ValueError:
        pass

    return None
