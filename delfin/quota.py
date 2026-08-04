"""Report HOME disk usage without ever walking the directory tree.

DELFIN runs on shared cluster file systems. Measuring usage with ``du`` there
means one ``stat()`` per file against a central metadata server -- on a HOME
with a few hundred thousand files that is hundreds of thousands of RPCs per
call, which degrades the file system for every user on the machine. Cluster
operations treat it as abuse and block accounts for it.

Every provider below instead asks the quota accounting the file system already
maintains, in a single call. If no provider matches, usage is reported as
unknown; falling back to a tree walk is never correct.
"""

from __future__ import annotations

import getpass
import os
import subprocess
from pathlib import Path
from typing import Callable, List, NamedTuple, Optional, Tuple

_TIMEOUT_SECONDS = 20

# Filesystem magic names as reported by ``stat -f -c %T`` / statvfs, mapped to
# the provider that understands them. Used only to order the attempts.
_FS_HINTS = {
    'lustre': 'lustre',
    'gpfs': 'gpfs',
    'beegfs': 'beegfs',
    'fhgfs': 'beegfs',
}


class QuotaUsage(NamedTuple):
    """Bytes used, and the soft limit if the file system reports one."""

    used_bytes: int
    soft_limit_bytes: Optional[int]
    source: str


def _run(cmd: List[str]) -> Optional[str]:
    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            check=False,
            timeout=_TIMEOUT_SECONDS,
        )
    except Exception:
        return None
    if result.returncode != 0:
        return None
    return result.stdout or None


def _int_or_none(token: str) -> Optional[int]:
    # Over-quota values carry a trailing '*'; unset limits show as '-' or '0'.
    token = token.strip().rstrip('*')
    if not token or not token.lstrip('-').isdigit():
        return None
    try:
        return int(token)
    except ValueError:
        return None


def parse_lfs_quota(output: str) -> Optional[Tuple[int, Optional[int]]]:
    """Parse ``lfs quota -q -u <user> <path>`` (Lustre).

    Fields: ``<mountpoint> <kbytes> <bquota> <blimit> <bgrace> <files> ...``.
    The mountpoint may wrap onto its own line, so numbers are collected from
    the first numeric run rather than by column index.
    """
    numbers: List[int] = []
    for token in output.split():
        value = _int_or_none(token)
        if value is not None:
            numbers.append(value)
        elif numbers:
            break
    if not numbers:
        return None
    soft = numbers[1] if len(numbers) > 1 and numbers[1] > 0 else None
    return numbers[0] * 1024, (soft * 1024 if soft else None)


def parse_posix_quota(output: str) -> Optional[Tuple[int, Optional[int]]]:
    """Parse ``quota -u <user>`` (POSIX / NFS quotas, blocks in 1K units).

    Layout::

        Disk quotas for user foo (uid 1000):
             Filesystem  blocks   quota   limit   grace   files   quota ...
              /dev/sda1  123456  200000  250000            1234      0 ...

    Multiple file systems may be listed; the first with non-zero usage wins,
    which is the common single-quota case and a safe default otherwise.
    """
    for line in output.splitlines():
        stripped = line.strip()
        if not stripped or stripped.lower().startswith('disk quotas'):
            continue
        fields = stripped.split()
        if len(fields) < 3 or not fields[0].startswith('/'):
            continue
        blocks = _int_or_none(fields[1])
        if blocks is None:
            continue
        soft = _int_or_none(fields[2])
        return blocks * 1024, (soft * 1024 if soft else None)
    return None


def parse_mmlsquota(output: str) -> Optional[Tuple[int, Optional[int]]]:
    """Parse ``mmlsquota -u <user> --block-size 1K <fs>`` (IBM Spectrum Scale).

    Layout::

        Filesystem type  KB      quota    limit  in_doubt  grace | files ...
        gpfs1      USR   123456  200000   250000        0   none |  1234 ...
    """
    for line in output.splitlines():
        fields = line.split()
        if len(fields) < 5:
            continue
        if fields[1].upper() not in {'USR', 'GRP', 'FILESET'}:
            continue
        blocks = _int_or_none(fields[2])
        if blocks is None:
            continue
        soft = _int_or_none(fields[3])
        return blocks * 1024, (soft * 1024 if soft else None)
    return None


def parse_beegfs_quota(output: str) -> Optional[Tuple[int, Optional[int]]]:
    """Parse ``beegfs-ctl --getquota --uid <user> --csv`` (BeeGFS, bytes).

    Layout: ``name,id,size,hard,files,hard`` with a header line.
    """
    for line in output.splitlines():
        fields = [f.strip() for f in line.split(',')]
        if len(fields) < 4:
            continue
        used = _int_or_none(fields[2])
        if used is None:
            continue
        hard = _int_or_none(fields[3])
        return used, (hard if hard else None)
    return None


def _filesystem_hint(path: Path) -> Optional[str]:
    """Best-effort file system family for *path*, or None if undetectable."""
    output = _run(['stat', '-f', '-c', '%T', str(path)])
    if not output:
        return None
    name = output.strip().lower()
    for needle, family in _FS_HINTS.items():
        if needle in name:
            return family
    return None


def _mount_point(path: Path) -> Path:
    """Deepest mount point containing *path* (falls back to *path* itself)."""
    try:
        current = path.resolve()
    except OSError:
        current = path
    try:
        while current != current.parent and not os.path.ismount(str(current)):
            current = current.parent
        return current
    except OSError:
        return path


def _provider_lustre(path: Path, user: str) -> Optional[Tuple[int, Optional[int]]]:
    output = _run(['lfs', 'quota', '-q', '-u', user, str(path)])
    return parse_lfs_quota(output) if output else None


def _provider_gpfs(path: Path, user: str) -> Optional[Tuple[int, Optional[int]]]:
    output = _run(
        ['mmlsquota', '-u', user, '--block-size', '1K', str(_mount_point(path))]
    )
    if not output:
        output = _run(['mmlsquota', '-u', user, '--block-size', '1K'])
    return parse_mmlsquota(output) if output else None


def _provider_beegfs(path: Path, user: str) -> Optional[Tuple[int, Optional[int]]]:
    output = _run(
        [
            'beegfs-ctl',
            '--getquota',
            '--uid',
            user,
            '--csv',
            f'--mount={_mount_point(path)}',
        ]
    )
    if not output:
        output = _run(['beegfs-ctl', '--getquota', '--uid', user, '--csv'])
    return parse_beegfs_quota(output) if output else None


def _provider_posix(path: Path, user: str) -> Optional[Tuple[int, Optional[int]]]:
    output = _run(['quota', '-u', user])
    return parse_posix_quota(output) if output else None


_PROVIDERS: List[Tuple[str, Callable[[Path, str], Optional[Tuple[int, Optional[int]]]]]] = [
    ('lustre', _provider_lustre),
    ('gpfs', _provider_gpfs),
    ('beegfs', _provider_beegfs),
    ('posix', _provider_posix),
]


def home_usage(home_dir, user: Optional[str] = None) -> Optional[QuotaUsage]:
    """Return quota-reported usage for *home_dir*, or None if unavailable.

    Providers are tried cheapest-first: the family matching the detected file
    system, then the remaining ones, since a site may export quotas through a
    tool that does not match the reported file system name. None means the
    caller should display "unknown" -- never a tree walk.
    """
    path = Path(home_dir)
    try:
        user = user or getpass.getuser()
    except Exception:
        return None

    hint = _filesystem_hint(path)
    ordered = sorted(_PROVIDERS, key=lambda item: item[0] != hint)

    for name, provider in ordered:
        try:
            result = provider(path, user)
        except Exception:
            continue
        if result and result[0] > 0:
            return QuotaUsage(result[0], result[1], name)
    return None
