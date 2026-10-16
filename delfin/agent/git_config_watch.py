"""A shell command does not leave a program in the repository's git config.

Git runs programs named in a repository's config (fsmonitor, filters,
pager, textconv, ssh command, credential helper, aliases with ``!``,
include paths that pull in more config, ...). A command in the sandbox
can write the workspace; on a host where the sandbox cannot keep the
config read-only (Landlock only, no sandbox, a layer that failed), it can
also write ``.git/config``. What it puts there would run in the next git
command outside the sandbox -- DELFIN's own reads go through ``safe_git``,
the user's own git in a terminal does not.

``take`` reads the config files of the workspace's git directories before
a command; ``check`` reads them after. A program-naming entry
(``safe_git.PROGRAM_KEY_RE``) that was not there before -- a new key, or
a known key with a new value -- puts the previous bytes back (atomically,
same mode) and returns one line for the model and the user. Every other
config change stays: ``git config user.name``, a new remote, an upstream.

Entries are read with ``git config --file <path> --null --list``, which
parses the file and runs nothing. Never raises.

Rejected: refusing the command up front (a script builds the path at run
time; the gate cannot read it); reverting every config change (costs
ordinary version control); checking only on hosts without bubblewrap (a
layer can fail at run time, and the check is cheap: two small files).
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .safe_git import PROGRAM_KEY_RE

#: Inside a git directory, the files whose entries name programs.
_WATCHED = ("config", "config.worktree")
_MAX_BYTES = 1 << 20


@dataclass
class Snapshot:
    files: dict = field(default_factory=dict)   # path -> bytes | None


def _git_dirs(workspace) -> list[str]:
    try:
        from .api_client import _git_metadata_roots
        return list(_git_metadata_roots(workspace))
    except Exception:
        return []


def _read(path: Path) -> Optional[bytes]:
    try:
        st = path.lstat()
    except OSError:
        return None
    try:
        import stat as _stat
        if not _stat.S_ISREG(st.st_mode) or st.st_size > _MAX_BYTES:
            return b""
        return path.read_bytes()
    except OSError:
        return None


def take(workspace) -> Snapshot:
    snap = Snapshot()
    for g in _git_dirs(workspace):
        for name in _WATCHED:
            p = Path(g) / name
            snap.files[str(p)] = _read(p)
    return snap


def _program_entries(data: Optional[bytes]) -> set[tuple[str, str]]:
    """(key, value) pairs that name a program, parsed by git from *data*."""
    if not data:
        return set()
    try:
        with tempfile.NamedTemporaryFile(prefix="delfin-gitcfg-", delete=False) as fh:
            fh.write(data)
            tmp = fh.name
        try:
            env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
            proc = subprocess.run(
                ["git", "config", "--file", tmp, "--null", "--list"],
                capture_output=True, timeout=10, env=env, check=False)
        finally:
            os.unlink(tmp)
    except Exception:
        return {("<unreadable>", "")}
    if proc.returncode != 0:
        # A file git cannot parse is a file git will refuse to use; it
        # names nothing that runs. Record it as such.
        return set()
    out: set[tuple[str, str]] = set()
    for entry in proc.stdout.split(b"\0"):
        if not entry:
            continue
        text = entry.decode("utf-8", "replace")
        key, _, value = text.partition("\n")
        key = key.strip()
        if key and PROGRAM_KEY_RE.fullmatch(key.lower()):
            out.add((key.lower(), value))
    return out


def _restore(path: Path, data: Optional[bytes]) -> bool:
    try:
        if data is None:
            path.unlink(missing_ok=True)
            return True
        # The file's own mode is kept; a file that cannot be read for its
        # mode gets owner-only.
        mode = 0o600
        try:
            mode = path.stat().st_mode & 0o777
        except OSError:
            pass
        fd, tmp = tempfile.mkstemp(prefix=".delfin-restore-", dir=str(path.parent))
        try:
            with os.fdopen(fd, "wb") as fh:
                fh.write(data)
                fh.flush()
                os.fsync(fh.fileno())
            os.chmod(tmp, mode)
            os.replace(tmp, path)
        except Exception:
            try:
                os.unlink(tmp)
            except OSError:
                pass
            raise
        return True
    except Exception:
        return False


def check(snap: Optional[Snapshot]) -> list[str]:
    """Undo a new program entry in a watched file; one note per file."""
    notes: list[str] = []
    if snap is None:
        return notes
    for path_s, before in snap.files.items():
        path = Path(path_s)
        after = _read(path)
        if after == before:
            continue
        added = _program_entries(after) - _program_entries(before)
        if not added:
            continue
        keys = sorted({k for k, _ in added})
        ok = _restore(path, before)
        try:
            from .api_client import _record_security_event
            _record_security_event(
                "git_config_program", "bash",
                f"{path} gained {', '.join(keys)[:200]}; "
                + ("restored" if ok else "could not be restored"),
                blocked=ok)
        except Exception:
            pass
        if ok:
            notes.append(
                f"{path}: the command added {', '.join(keys)}, which names a "
                "program git would run outside the sandbox. The file is "
                "back as it was before the command (with it, any other "
                "change the command made to this file). A config change "
                "that names no program is kept; a wanted program key is "
                "for the user to set in a terminal.")
        else:
            notes.append(
                f"{path}: the command added {', '.join(keys)}, which names a "
                "program git would run outside the sandbox, and the file "
                "could not be restored. Check it before the next git command.")
    return notes
