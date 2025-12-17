"""
File I/O utility functions for DELFIN.

Provides consistent file reading/writing helpers with error handling.
"""

from pathlib import Path
from typing import Optional, Union


def safe_read_text(
    path: Union[str, Path],
    encoding: str = 'utf-8',
    errors: str = 'ignore'
) -> Optional[str]:
    """Safely read text file with consistent error handling.

    Args:
        path: File path to read
        encoding: Text encoding to use
        errors: How to handle encoding errors

    Returns:
        File content as string, or None if file not found
    """
    try:
        with open(path, 'r', encoding=encoding, errors=errors) as f:
            return f.read()
    except FileNotFoundError:
        return None


def safe_read_lines(
    path: Union[str, Path],
    encoding: str = 'utf-8',
    errors: str = 'ignore'
) -> Optional[list[str]]:
    """Safely read text file as list of lines.

    Args:
        path: File path to read
        encoding: Text encoding to use
        errors: How to handle encoding errors

    Returns:
        List of lines, or None if file not found
    """
    try:
        with open(path, 'r', encoding=encoding, errors=errors) as f:
            return f.readlines()
    except FileNotFoundError:
        return None
