"""Path utilities for DELFIN."""
from __future__ import annotations

import os
from pathlib import Path
from typing import Union

SCRATCH_ENV = "DELFIN_SCRATCH"


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

__all__ = [
    "SCRATCH_ENV",
    "resolve_path",
    "get_runtime_dir",
    "scratch_path",
]
