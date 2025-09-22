"""Path utilities for DELFIN."""
from pathlib import Path
from typing import Union


def resolve_path(path: Union[str, Path]) -> Path:
    """Return absolute, user-expanded path without altering unresolved behaviour."""
    candidate = Path(path).expanduser()
    try:
        return candidate.resolve()
    except FileNotFoundError:
        return candidate
