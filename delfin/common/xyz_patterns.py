"""
Shared regex patterns and utilities for XYZ file parsing.
"""

import re


# XYZ coordinate line pattern (element + X Y Z coordinates)
XYZ_COORD_LINE_RE = re.compile(
    r"^\s*[A-Za-z]{1,2}[A-Za-z0-9()]*\s+"      # Atom label, optional index
    r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?\s+"      # X coordinate
    r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?\s+"      # Y coordinate
    r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?"         # Z coordinate
)


def count_xyz_coord_lines(lines) -> int:
    """Return number of lines that look like XYZ coordinates.

    Args:
        lines: List of text lines to check

    Returns:
        Count of lines matching XYZ coordinate pattern
    """
    return sum(1 for line in lines if XYZ_COORD_LINE_RE.match(line))
