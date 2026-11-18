from __future__ import annotations

from pathlib import Path
from typing import Any


def report_file_name(value: Any, *, default: str) -> str:
    name = str(value or "").strip() or default
    if (name in (".", "..") or "/" in name or "\\" in name
            or Path(name).name != name):
        raise ValueError(
            f"'report_name' must be a file name inside work_dir, not a path: {name!r}"
        )
    return name
