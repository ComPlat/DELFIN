"""CE1: reviewer s21 findings F1/F2 on crest_ensemble (0b4e85e6).

F1: the bare-float fallback in _energy_from_comment takes the FIRST
parsable float, so a comment like "2 -12.34 Eh" (CREST variants write a
frame index first) parses energy 2.0. Convention pinned here: the
fallback takes the LAST parsable float in the comment.
F2: a NaN energy parses and propagates silently into relative_kcal's
min()-based ranking. Pinned here: non-finite energies are refused at
parse time.

Both cases RED on 0b4e85e6, GREEN after the fix commit.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from delfin.tools.adapters.crest_ensemble import (
    parse_ensemble,
    relative_kcal,
)


def _write(tmp_path: Path, text: str) -> Path:
    p = tmp_path / "ensemble.xyz"
    p.write_text(textwrap.dedent(text))
    return p


def test_first_frame_index_is_not_read_as_energy(tmp_path: Path):
    p = _write(tmp_path, """\
        3
        1 -12.34 Eh
        C 0.0 0.0 0.0
        H 1.09 0.0 0.0
        H -0.55 0.95 0.0
        """)
    frames = parse_ensemble(p)
    assert frames[0]["energy_Eh"] == pytest.approx(-12.34)


def test_nan_energy_is_refused_not_propagated(tmp_path: Path):
    p = _write(tmp_path, """\
        3
        energy: nan Eh
        C 0.0 0.0 0.0
        H 1.09 0.0 0.0
        H -0.55 0.95 0.0
        """)
    with pytest.raises(ValueError, match="finite"):
        parse_ensemble(p)


def test_relative_kcal_never_sees_a_non_finite_value():
    """Belt and braces: the ranking helper itself refuses non-finite."""
    frames = [{"energy_Eh": float("nan"), "degeneracy": 1}]
    with pytest.raises(ValueError, match="finite"):
        relative_kcal(frames)
