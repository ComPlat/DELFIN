"""A requested conformer search has a skill that gets charge and spin right.

Input: a request for a conformer search or a global minimum (ORCA GOAT,
CREST). Before: no skill covered it, so the model wrote GOAT and CREST
inputs from memory. The errors this skill exists for are each a silent
wrong result, not a crash:
- xtb and CREST take the number of unpaired electrons (``--uhf``), ORCA
  takes the multiplicity, and passing one as the other runs the wrong
  spin state;
- xtb reads ``.CHRG`` / ``.UHF`` files in the working directory, which
  override a missing flag;
- GOAT gives every conformer degeneracy 1 and no factor for mirror
  images, so populations need the degeneracy counted explicitly.
The skill runs what was asked and offers refinement; it is not a
standard protocol the model applies on its own.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from delfin.agent import skills

NAME = "conformer-search"


@pytest.fixture
def skill(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    sk = skills.get_skill(NAME, domain="code")
    assert sk is not None, "no conformer-search skill in the code domain"
    return sk


def test_it_is_listed_with_a_description(skill):
    assert "GOAT" in skill.description and "CREST" in skill.description


def test_it_tells_unpaired_electrons_from_multiplicity(skill):
    body = skill.body
    assert "--uhf" in body and "2S+1" in body
    assert re.search(r"unpaired electrons\s*\(2S\)", body)
    assert ".CHRG" in body and ".UHF" in body


def test_it_counts_degeneracy_for_populations(skill):
    assert "g_i exp(-E_i/RT)" in skill.body
    assert "mirror images" in skill.body.lower()


def test_it_runs_what_was_asked_and_offers_the_rest(skill):
    body = " ".join(skill.body.split())
    assert "run it only when it was asked for" in body
    assert "only as an offer" in body


def test_it_names_no_path_or_credential():
    from tests.test_c1_skill_text_lint import _violations
    text = (Path(__file__).resolve().parents[1]
            / "delfin/agent/pack/skills" / f"{NAME}.md").read_text()
    assert _violations(text) == []


def test_it_says_where_the_results_are(skill):
    # File names as ORCA 6.1 (manual section 4.10) and CREST write them.
    for name in (".globalminimum.xyz", ".finalensemble.xyz",
                 "# Final ensemble info #", "crest_conformers.xyz",
                 "crest_best.xyz"):
        assert name in skill.body, name


def test_a_rotation_scan_is_checked_against_the_ensemble(skill):
    body = " ".join(skill.body.split())
    assert "every minimum of the scan must be a conformer of the ensemble" \
        in body
