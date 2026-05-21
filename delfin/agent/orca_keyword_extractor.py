"""ORCA-specific wrapper around the generic ``manual_extractor``.

Kept as a thin backward-compat layer so older imports
``from delfin.agent.orca_keyword_extractor import …`` still work.
New code should use ``manual_extractor`` directly with a
``ProgramConfig``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .manual_extractor import (
    ORCA_CONFIG,
    extract_namespace,
    is_real_keyword as _is_real_keyword_generic,
)


def extract_keyword_namespace(
    index_path: Path | None = None,
    *,
    sanitize_section_ids: bool = True,
) -> dict[str, dict[str, Any]]:
    """Extract the ORCA keyword namespace.  Wrapper over the generic
    pipeline with the built-in ``ORCA_CONFIG``."""
    return extract_namespace(
        ORCA_CONFIG,
        index_path=index_path,
        sanitize_section_ids=sanitize_section_ids,
    )


def is_real_keyword(
    block: str, keyword: str,
    namespace: dict[str, dict[str, Any]] | None = None,
    path: Path | None = None,
) -> bool:
    """Backward-compatible ORCA-only wrapper."""
    return _is_real_keyword_generic(
        program="orca", block=block, keyword=keyword,
        namespace=namespace, path=path,
    )


__all__ = ["extract_keyword_namespace", "is_real_keyword"]
