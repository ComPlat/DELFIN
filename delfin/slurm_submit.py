"""Which SLURM partitions a job may start in, and time limits sbatch accepts.

Kept free of the dashboard on purpose: ``delfin-step --slurm``,
``delfin-pipeline --slurm`` and the tools runtime submit jobs as well, and
importing the dashboard backend just to ask costs four seconds and loads
ipywidgets. The dashboard's backend calls these same functions, so a job
submitted from a terminal and one submitted from the dashboard are placed
alike.
"""

from __future__ import annotations

import os
import re
import socket
import subprocess
from pathlib import Path
from typing import Mapping, Optional, Sequence

#: Partitions a site's CPU jobs may start in. SLURM starts a job listed for
#: several partitions in whichever can run it first, "with no regard given to
#: the partition name ordering" (sbatch(1)), so a second partition of suitable
#: nodes is a second queue at no cost to the request. Measured on bwUniCluster
#: 3.0 with sbatch --test-only for 40 cores, 240 GB, 2 days: cpu (80 nodes)
#: expected a start twelve days after cpu_il (264 nodes).
PROFILE_PARTITIONS: dict[str, tuple[str, ...]] = {
    'bwunicluster3': ('cpu', 'cpu_il'),
}


def detect_site_profile() -> str:
    """The site profile of this machine from its host name, or ''."""
    try:
        fqdn = socket.getfqdn().lower()
    except Exception:
        fqdn = ''
    hostname = os.environ.get('HOSTNAME', '').lower() or fqdn.split('.')[0]
    if 'scc.kit.edu' in fqdn or hostname.startswith('uc3'):
        return 'bwunicluster3'
    return ''


def partition_list(value) -> tuple:
    """'cpu, cpu_il' or ['cpu', 'cpu_il'] as ('cpu', 'cpu_il'), each once."""
    if not value:
        return ()
    items = value.split(',') if isinstance(value, str) else list(value)
    out: list[str] = []
    for item in items:
        name = str(item).strip()
        if name and name not in out:
            out.append(name)
    return tuple(out)


def configured_partitions(explicit=None, *, profile: Optional[str] = None) -> tuple:
    """The candidates: given, else ``runtime.slurm.partitions`` in the settings,
    else ``DELFIN_SLURM_PARTITIONS``, else the site profile's."""
    found = partition_list(explicit)
    if not found:
        try:
            from delfin.user_settings import load_settings

            slurm = (load_settings().get('runtime') or {}).get('slurm') or {}
            found = partition_list(slurm.get('partitions', ''))
        except Exception:
            found = ()
    if not found:
        found = partition_list(os.environ.get('DELFIN_SLURM_PARTITIONS', ''))
    if not found:
        found = PROFILE_PARTITIONS.get(detect_site_profile() if profile is None else profile, ())
    return found


def partition_accepts(name, job_dir, submit_script, *, resource_args: Sequence[str] = (),
                      submit_env: Optional[Mapping[str, str]] = None,
                      sbatch: str = 'sbatch') -> Optional[bool]:
    """True or False as ``sbatch --test-only`` answered, None when SLURM could
    not be asked. Nothing is submitted."""
    try:
        done = subprocess.run(
            [sbatch, '--test-only', f'--partition={name}', *resource_args, str(submit_script)],
            cwd=str(job_dir), capture_output=True, text=True,
            env=dict(submit_env) if submit_env is not None else None, timeout=30,
        )
    except Exception:
        return None
    said = f'{done.stdout or ""}\n{done.stderr or ""}'
    if done.returncode == 0 and 'to start at' in said:
        return True
    if 'allocation failure' in said or 'error' in said.lower():
        return False
    return None


def choose_partitions(submit_script, job_dir, *, partitions=None, resource_args: Sequence[str] = (),
                      submit_env: Optional[Mapping[str, str]] = None, sbatch: str = 'sbatch',
                      profile: Optional[str] = None) -> str:
    """The candidate partitions a job fits, as one ``--partition`` value.

    A partition whose nodes cannot hold the request makes SLURM reject the
    whole job when it is listed, so each candidate is asked first. '' means
    nothing to choose: the script's own partition, or the cluster default.
    """
    candidates = configured_partitions(partitions, profile=profile)
    if len(candidates) <= 1:
        return candidates[0] if candidates else ''
    return ','.join(
        name for name in candidates
        if partition_accepts(name, job_dir, submit_script, resource_args=resource_args,
                             submit_env=submit_env, sbatch=sbatch)
    )


def script_sets_partition(script_text: str) -> bool:
    return bool(re.search(r'^#SBATCH[ \t]+(?:-p\b|--partition\b)', script_text, re.MULTILINE))


def sbatch_command(sbatch: str, submit_script, job_dir=None) -> list:
    """``sbatch`` for a script that states its own resources, with the
    partitions it fits unless the script names one.

    Never raises: when nothing can be chosen the command is plain
    ``sbatch script``, which is what these callers ran before.
    """
    path = Path(submit_script)
    cmd = [sbatch]
    try:
        if not script_sets_partition(path.read_text(encoding='utf-8', errors='replace')):
            chosen = choose_partitions(path, Path(job_dir) if job_dir else path.parent, sbatch=sbatch)
            if chosen:
                cmd.append(f'--partition={chosen}')
    except Exception:
        pass
    cmd.append(str(path))
    return cmd


def time_limit_seconds(time_limit) -> int:
    """SLURM's time formats in seconds: D-HH:MM:SS, D-HH:MM, D-HH,
    HH:MM:SS, MM:SS, or bare minutes."""
    text = str(time_limit).strip()
    if '-' in text:
        day_part, text = text.split('-', 1)
        parts = text.split(':')
        # After a day count the first field is hours, not minutes.
        while len(parts) < 3:
            parts.append('0')
        hours, minutes, seconds = (int(x) for x in parts[:3])
        return int(day_part) * 86400 + hours * 3600 + minutes * 60 + seconds
    parts = text.split(':')
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + int(parts[2])
    if len(parts) == 2:
        return int(parts[0]) * 60 + int(parts[1])
    return int(parts[0]) * 60  # SLURM treats a bare number as minutes


def normalize_time_limit(value) -> str:
    """A time limit sbatch accepts, from what was typed into a time field.

    SLURM's own forms pass through (``48:00:00``, ``2-00:00:00``, ``90``).
    Durations as people write them are converted: ``48h`` and ``2d`` become
    ``2-00:00:00``, ``90min`` becomes ``01:30:00``. Anything else raises
    ValueError with the reason.
    """
    text = str(value if value is not None else '').strip()
    if not text:
        raise ValueError('the time limit is empty')
    if text.upper() in ('UNLIMITED', 'INFINITE'):
        return text.upper()
    compact = text.replace(' ', '').lower()
    found = re.fullmatch(r'(?:(\d+)d)?(?:(\d+)h)?(?:(\d+)(?:min|m))?(?:(\d+)s)?', compact)
    if found and any(found.groups()):
        days, hours, minutes, seconds = (int(part or 0) for part in found.groups())
        total = days * 86400 + hours * 3600 + minutes * 60 + seconds
        if total <= 0:
            raise ValueError('the time limit is zero')
        whole_days, rest = divmod(total, 86400)
        hours, rest = divmod(rest, 3600)
        minutes, seconds = divmod(rest, 60)
        clock = f'{hours:02d}:{minutes:02d}:{seconds:02d}'
        return f'{whole_days}-{clock}' if whole_days else clock
    if not re.fullmatch(r'\d+-\d+(?::\d+){0,2}|\d+(?::\d+){0,2}', text):
        raise ValueError('it is not a time SLURM understands')
    if time_limit_seconds(text) <= 0:
        raise ValueError('the time limit is zero')
    return text
