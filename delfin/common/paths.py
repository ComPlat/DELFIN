"""Path utilities for DELFIN."""
from __future__ import annotations

import os
import threading
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Union

SCRATCH_ENV = "DELFIN_SCRATCH"

# Process-global CWD guard.
#
# ``os.chdir`` mutates the working directory of the WHOLE process, not just the
# calling thread. DELFIN runs OCCUPIER/FoB jobs concurrently in a single-process
# ThreadPoolExecutor, so any thread that ``chdir``s while another thread does
# CWD-relative file I/O can make that I/O land in the wrong directory. This is a
# real defect: it leaked OCCUPIER ``CONTROL.txt`` writes from a sub-folder into
# the main run directory. Every thread that must ``chdir`` has to serialize on
# THIS single re-entrant lock (RLock — nested chdir helpers occur) and hold it
# for the entire chdir → work → restore span. Use ``pushd`` for that.
GLOBAL_CWD_LOCK = threading.RLock()


@contextmanager
def pushd(path: Union[str, Path]) -> Iterator[Path]:
    """Temporarily ``chdir`` into *path*, serialized on :data:`GLOBAL_CWD_LOCK`.

    The process-global lock is held for the whole ``with`` block, so no other
    thread can change the working directory underneath the caller. The previous
    working directory is always restored on exit.
    """
    target = resolve_path(path)
    with GLOBAL_CWD_LOCK:
        prev = os.getcwd()
        os.chdir(target)
        try:
            yield target
        finally:
            try:
                os.chdir(prev)
            except OSError:
                pass


def resolve_path(path: Union[str, Path]) -> Path:
    """Return absolute, user-expanded path without altering unresolved behaviour."""
    candidate = Path(path).expanduser()
    try:
        return candidate.resolve()
    except FileNotFoundError:
        return candidate


def get_runtime_dir() -> Path:
    """Return the working directory for DELFIN runs (respects DELFIN_SCRATCH)."""
    env_path = os.getenv(SCRATCH_ENV)
    if env_path:
        candidate = Path(env_path).expanduser()
        candidate.mkdir(parents=True, exist_ok=True)
        try:
            return candidate.resolve()
        except FileNotFoundError:
            return candidate
    return Path.cwd()


def scratch_path(*parts: Union[str, Path], create: bool = False) -> Path:
    """Build a path relative to the runtime directory (scratch-aware)."""
    base = get_runtime_dir()
    if not parts:
        path = base
    else:
        path = base.joinpath(*(Path(p) for p in parts))
    if create:
        path.parent.mkdir(parents=True, exist_ok=True)
    return path


def ensure_relative_link(src: Union[str, Path], dest: Union[str, Path]) -> bool:
    """Create a relative symlink (or hardlink fallback) at dest pointing to src.

    Returns True if dest already exists or the link was created successfully.
    """
    src_path = Path(src)
    dest_path = Path(dest)

    if dest_path.exists() or dest_path.is_symlink():
        return True
    if not src_path.exists():
        return False

    try:
        rel = os.path.relpath(src_path, start=dest_path.parent)
        dest_path.symlink_to(rel)
        return True
    except OSError:
        pass

    try:
        os.link(src_path, dest_path)
        return True
    except OSError:
        return False

__all__ = [
    "SCRATCH_ENV",
    "GLOBAL_CWD_LOCK",
    "pushd",
    "resolve_path",
    "get_runtime_dir",
    "scratch_path",
    "ensure_relative_link",
]
