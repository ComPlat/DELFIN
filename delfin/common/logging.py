"""Shared logging helpers for DELFIN."""
from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional, Dict

from delfin.common.progress_display import ProgressAwareStreamHandler

DEFAULT_LOG_LEVEL = logging.INFO
DEFAULT_LOG_FORMAT = "%(levelname)s: %(message)s"
FILE_LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
FILE_LOG_DATEFMT = "%Y-%m-%d %H:%M:%S"

_FILE_HANDLERS: Dict[str, logging.Handler] = {}


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

    # Remove existing handlers if forcing reconfiguration
    if force:
        for handler in root_logger.handlers[:]:
            root_logger.removeHandler(handler)
            handler.close()

    # Create a progress-aware stream handler instead of using basicConfig
    handler = ProgressAwareStreamHandler(stream or sys.stderr)
    handler.setLevel(level or DEFAULT_LOG_LEVEL)
    formatter = logging.Formatter(fmt or DEFAULT_LOG_FORMAT)
    handler.setFormatter(formatter)

    root_logger.addHandler(handler)
    root_logger.setLevel(level or DEFAULT_LOG_LEVEL)


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """Return a module-level logger; mirrors logging.getLogger without side effects."""
    if name is None:
        return logging.getLogger()
    return logging.getLogger(name)


def add_file_handler(
    path,
    level: Optional[int] = None,
    mode: str = "a",
    fmt: str = FILE_LOG_FORMAT,
    datefmt: str = FILE_LOG_DATEFMT,
) -> logging.Handler:
    """Attach a single shared file handler to the root logger."""
    target = Path(path).expanduser()
    target.parent.mkdir(parents=True, exist_ok=True)
    resolved = str(target.resolve())

    if resolved in _FILE_HANDLERS:
        handler = _FILE_HANDLERS[resolved]
        if level is not None:
            handler.setLevel(level)
        return handler

    handler = logging.FileHandler(resolved, mode=mode, encoding="utf-8")
    handler.setLevel(level or DEFAULT_LOG_LEVEL)
    handler.setFormatter(logging.Formatter(fmt=fmt, datefmt=datefmt))

    root_logger = logging.getLogger()
    root_logger.addHandler(handler)

    _FILE_HANDLERS[resolved] = handler
    return handler
