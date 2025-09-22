"""Shared logging helpers for DELFIN."""
from __future__ import annotations

import logging
import sys
from typing import Optional

DEFAULT_LOG_LEVEL = logging.INFO
DEFAULT_LOG_FORMAT = "%(levelname)s: %(message)s"


def configure_logging(
    level: Optional[int] = None,
    fmt: Optional[str] = None,
    stream = None,
    force: bool = False,
) -> None:
    """Configure root logging once with DELFIN defaults.

    - `level`: defaults to INFO if not provided.
    - `fmt`: default format keeps CLI output compact.
    - `stream`: defaults to stderr (ideal for SLURM job output).
    - `force`: pass True to reconfigure even if handlers exist.
    """

    root_logger = logging.getLogger()
    if root_logger.handlers and not force:
        return

    logging.basicConfig(
        level=level or DEFAULT_LOG_LEVEL,
        format=fmt or DEFAULT_LOG_FORMAT,
        stream=stream or sys.stderr,
        force=force,
    )


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """Return a module-level logger; mirrors logging.getLogger without side effects."""
    if name is None:
        return logging.getLogger()
    return logging.getLogger(name)
