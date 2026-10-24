"""CE1: level-tag convention and degeneracy policy of crest_ensemble.

Reviewer s21 findings, pinned as tests (written against the phase-1 commit
9b3a7eda; the level-tag storing and the partial-tag refusal are NEW rules —
these cases must be RED before the fix commit).

Convention being pinned:
- A comment line that carries a recognizable level token (gfn*, xtb*, or a
  named functional such as r2scan-3c/b3lyp/wb97x/pbe0) tags the frame's
  level of theory.
- If SOME frames are tagged and others are not, the level is undetectable
  for the ensemble as a whole: weighting refuses loudly (ValueError), it
  never silently assumes one level.
- If NO frame carries a tag (plain CREST comments), the energies are on the
  search's single surface; weighting accepts and the caller states the
  level in the report.
- Degeneracy: parse_ensemble always reports degeneracy 1 per frame; a
  mirror-image pair (same energy) is grouped only by an explicit
  collapse_degeneracies call, never silently.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from delfin.tools.adapters.crest_ensemble import (
    parse_ensemble,
    boltzmann_weights,
)


def _write(tmp_path: Path, text: str) -> Path:
    p = tmp_path / "ensemble.xyz"
    p.write_text(textwrap.dedent(text))
    return p


def test_parse_stores_the_level_tag_from_the_comment(tmp_path: Path):
    p = _write(tmp_path, """\
        3
        energy: -12.34 Eh gfn2
        C 0.0 0.0 0.0
        H 1.09 0.0 0.0
        H -0.55 0.95 0.0
        """)
    frames = parse_ensemble(p)
    assert frames[0]["method"] == "gfn2"


def test_a_functional_tag_is_recognized_as_the_level(tmp_path: Path):
    p = _write(tmp_path, """\
        3
        energy: -500.1 Eh r2scan-3c
        C 0.0 0.0 0.0
        H 1.09 0.0 0.0
        H -0.55 0.95 0.0
        """)
    frames = parse_ensemble(p)
    assert frames[0]["method"] == "r2scan-3c"


def test_partial_tags_refuse_weighting_loudly(tmp_path: Path):
    """Some frames tagged, one not: the level of the ensemble is
    undetectable — weighting must refuse, not assume."""
    p = _write(tmp_path, """\
        3
        energy: -12.34 Eh gfn2
        C 0.0 0.0 0.0
        H 1.09 0.0 0.0
        H -0.55 0.95 0.0
        3
        energy: -12.30 Eh
        C 0.0 0.0 0.0
        H 1.09 0.0 0.0
        H 0.55 0.95 0.0
        """)
    frames = parse_ensemble(p)
    with pytest.raises(ValueError, match="level"):
        boltzmann_weights(frames, temperature_K=298.15)


def test_fully_untagged_ensemble_is_weighted_at_the_search_surface(tmp_path: Path):
    """No tags anywhere: the plain CREST case — accept, single surface."""
    p = _write(tmp_path, """\
        3
        energy: -12.34 Eh
        C 0.0 0.0 0.0
        H 1.09 0.0 0.0
        H -0.55 0.95 0.0
        """)
    frames = parse_ensemble(p)
    weights = boltzmann_weights(frames, temperature_K=298.15)
    assert weights[0] == pytest.approx(1.0)


def test_parse_reports_degeneracy_one_per_frame(tmp_path: Path):
    """parse_ensemble never invents a degeneracy; grouping is explicit."""
    p = _write(tmp_path, """\
        3
        energy: -12.340000 Eh
        C 0.0 0.0 0.0
        H 1.09 0.0 0.0
        H -0.55 0.95 0.0
        3
        energy: -12.340000 Eh
        C 0.0 0.0 0.0
        H 1.09 0.0 0.0
        H 0.55 -0.95 0.0
        """)
    frames = parse_ensemble(p)
    assert all(int(f["degeneracy"]) == 1 for f in frames)
