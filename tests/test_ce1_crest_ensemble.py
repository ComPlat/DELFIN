"""CE1: ensemble parsing and Boltzmann weighting for CREST/CENSO ensembles.

Red control for delfin/tools/adapters/crest_ensemble.py. Written against the
unchanged tree: every case here must FAIL before the module exists.

Rules under test (wave-20 priority 1, scientific correctness):
- Boltzmann weights at a stated temperature sum to 1 and are computed from
  ONE level of theory only; a mixture is refused, not silently ranked.
- Degenerate conformers (equal energy) keep their degeneracy explicitly
  instead of being collapsed silently.
- Energies carry their unit (Eh); relative energies are derived, never
  assumed.
"""

from __future__ import annotations

import math
import textwrap
from pathlib import Path

import pytest

from delfin.tools.adapters.crest_ensemble import (
    parse_ensemble,
    boltzmann_weights,
    relative_kcal,
    collapse_degeneracies,
)

TWO_FRAME_ENSEMBLE = textwrap.dedent("""\
    3
    energy: -12.345670 Eh
    C 0.0 0.0 0.0
    H 1.09 0.0 0.0
    H -0.55 0.95 0.0
    3
    energy: -12.340000 Eh
    C 0.0 0.0 0.0
    H 1.09 0.0 0.0
    H 0.55 0.95 0.0
    """)


def test_parse_reads_frames_with_energy_in_eh(tmp_path: Path):
    p = tmp_path / "crest_conformers.xyz"
    p.write_text(TWO_FRAME_ENSEMBLE)
    frames = parse_ensemble(p)
    assert len(frames) == 2
    assert frames[0]["energy_Eh"] == pytest.approx(-12.345670)
    assert frames[1]["energy_Eh"] == pytest.approx(-12.340000)
    assert frames[0]["symbols"] == ["C", "H", "H"]


def test_relative_energies_in_kcal_are_derived_from_the_lowest(tmp_path: Path):
    p = tmp_path / "crest_conformers.xyz"
    p.write_text(TWO_FRAME_ENSEMBLE)
    frames = parse_ensemble(p)
    rel = relative_kcal(frames)
    assert rel[0] == pytest.approx(0.0, abs=1e-9)
    # (−12.340000 + 12.345670) Eh * 627.5094740631 kcal/mol per Eh
    expected = (-12.340000 + 12.345670) * 627.5094740631
    assert rel[1] == pytest.approx(expected, rel=1e-6)


def test_boltzmann_weights_sum_to_one_at_the_stated_temperature(tmp_path: Path):
    p = tmp_path / "crest_conformers.xyz"
    p.write_text(TWO_FRAME_ENSEMBLE)
    frames = parse_ensemble(p)
    weights = boltzmann_weights(frames, temperature_K=298.15)
    assert sum(weights) == pytest.approx(1.0, abs=1e-12)
    assert weights[0] > weights[1]  # the lower-energy frame is more populated


def test_boltzmann_weights_are_temperature_dependent():
    energies = [{"energy_Eh": e} for e in (-12.0, -11.9)]
    hot = boltzmann_weights(energies, temperature_K=2000.0)
    cold = boltzmann_weights(energies, temperature_K=100.0)
    assert hot[1] > cold[1]  # high T evens the populations out


def test_boltzmann_refuses_energies_from_two_levels():
    """The trap: an ensemble ranked with energies from two different levels."""
    mixed = [
        {"energy_Eh": -12.0, "method": "gfn2"},
        {"energy_Eh": -11.9, "method": "r2scan-3c"},
    ]
    with pytest.raises(ValueError, match="level of theory"):
        boltzmann_weights(mixed, temperature_K=298.15)


def test_degenerate_conformers_are_counted_not_collapsed(tmp_path: Path):
    frames_txt = TWO_FRAME_ENSEMBLE + textwrap.dedent("""\
        3
        energy: -12.340000 Eh
        C 0.0 0.0 0.0
        H 1.09 0.0 0.0
        H 0.55 -0.95 0.0
        """)
    p = tmp_path / "crest_conformers.xyz"
    p.write_text(frames_txt)
    frames = parse_ensemble(p)
    groups = collapse_degeneracies(frames, tolerance_kcal=0.01)
    assert len(groups) == 2
    lowest = groups[0]
    assert lowest["energy_Eh"] == pytest.approx(-12.345670)
    assert lowest["degeneracy"] == 1
    assert groups[1]["degeneracy"] == 2


def test_boltzmann_weights_respect_degeneracy(tmp_path: Path):
    frames_txt = TWO_FRAME_ENSEMBLE + textwrap.dedent("""\
        3
        energy: -12.340000 Eh
        C 0.0 0.0 0.0
        H 1.09 0.0 0.0
        H 0.55 -0.95 0.0
        """)
    p = tmp_path / "crest_conformers.xyz"
    p.write_text(frames_txt)
    frames = parse_ensemble(p)
    groups = collapse_degeneracies(frames, tolerance_kcal=0.01)
    weights = boltzmann_weights(groups, temperature_K=298.15)
    assert sum(weights) == pytest.approx(1.0, abs=1e-12)
    # the twofold-degenerate second group carries twice its bare weight
    bare = boltzmann_weights(
        [{"energy_Eh": g["energy_Eh"], "degeneracy": 1} for g in groups],
        temperature_K=298.15,
    )
    assert weights[1] == pytest.approx(2 * bare[1] / (bare[0] + 2 * bare[1]),
                                       rel=1e-9)
