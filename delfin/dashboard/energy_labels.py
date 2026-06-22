"""Helpers for formatting energy-bearing XYZ frame labels."""

from __future__ import annotations

import html
import re

HARTREE_TO_KCAL_PER_MOL = 627.509474
_COMMENT_ENERGY_RE = re.compile(
    r'(?<![A-Za-z0-9_])([+-]?(?:(?:\d+\.\d*)|(?:\d*\.\d+)|(?:\d+(?:[Ee][+-]?\d+))))(?![A-Za-z0-9_])'
)


def extract_comment_energy(comment):
    """Extract a Hartree energy from an XYZ comment line when present."""
    text = str(comment or "").strip()
    if not text:
        return None
    epreopt_match = re.search(
        r'Epreopt\s*=\s*([+-]?(?:(?:\d+\.\d*)|(?:\d*\.\d+)|(?:\d+(?:[Ee][+-]?\d+))))',
        text,
        flags=re.IGNORECASE,
    )
    candidate = epreopt_match.group(1) if epreopt_match else None
    if candidate is None:
        generic_match = _COMMENT_ENERGY_RE.search(text)
        candidate = generic_match.group(1) if generic_match else None
    if candidate is None:
        return None
    try:
        return float(candidate)
    except Exception:
        return None


def format_xyz_comment_label(comment, frames=None, max_chars=100):
    """Render frame comment text and append absolute energy information."""
    raw_comment = str(comment or "")
    clipped_comment = html.escape(raw_comment[:max_chars])
    if len(raw_comment) > max_chars:
        clipped_comment += "..."

    energy_hartree = extract_comment_energy(raw_comment)
    if energy_hartree is None:
        return clipped_comment

    abs_kcal = energy_hartree * HARTREE_TO_KCAL_PER_MOL
    delta_line = ""
    frame_list = list(frames or [])
    if len(frame_list) > 1:
        try:
            reference_energy = extract_comment_energy(frame_list[0][0])
        except Exception:
            reference_energy = None
        if reference_energy is not None:
            delta_kcal = (energy_hartree - reference_energy) * HARTREE_TO_KCAL_PER_MOL
            delta_line = f"<br>ΔE = {delta_kcal:.2f} kcal/mol"
    return (
        "<span style='color:#666;font-size:0.9em;'>"
        f"E = {energy_hartree:.12f} Hartree; "
        f"{abs_kcal:.2f} kcal/mol"
        f"{delta_line}"
        "</span>"
    )
