"""Shared logging helpers for DELFIN."""
import logging
from typing import Optional


def get_logger(name: Optional[str] = None) -> logging.Logger:
    """Return a module-level logger; mirrors logging.getLogger without side effects."""
    if name is None:
        return logging.getLogger()
    return logging.getLogger(name)
