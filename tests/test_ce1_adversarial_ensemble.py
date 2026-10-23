"""CE1 adversarial tests for delfin/tools/adapters/crest_ensemble.py.

Reviewer nacht-s21. Each test pins a behaviour the parser must have on
hostile or edge inputs. Failures are findings against the builder's
commits 9b3a7eda / 0b4e85e6.
"""

from pathlib import Path

import pytest

from delfin.tools.adapters.crest_ensemble import (
    boltzmann_weights,
    collapse_degeneracies,
    parse_ensemble,
    relative_kcal,
)


def _write_xyz(tmp_path: Path, body: str) -> Path:
    p = tmp_path / "ensemble.xyz"
    p.write_text(body)
    return p


# --- edge: single frame ------------------------------------------------------

def test_single_frame_weights_are_exactly_one(tmp_path):
    p = _write_xyz(
        tmp_path,
        "1\nenergy: -12.34 Eh\nH 0.0 0.0 0.0\n",
    )
    frames = parse_ensemble(p)
    w = boltzmann_weights(frames, temperature_K=298.15)
    assert len(w) == 1
    assert w[0] == pytest.approx(1.0, abs=1e-12)


# --- edge: duplicate energies ------------------------------------------------

def test_duplicate_energies_get_equal_weights(tmp_path):
    body = "".join(
        "1\nenergy: -12.34 Eh\nH 0.0 0.0 0.0\n" for _ in range(3)
    )
    frames = parse_ensemble(_write_xyz(tmp_path, body))
    w = boltzmann_weights(frames, temperature_K=298.15)
    assert w[0] == pytest.approx(1.0 / 3.0, abs=1e-12)
    assert w[1] == pytest.approx(1.0 / 3.0, abs=1e-12)
    assert w[2] == pytest.approx(1.0 / 3.0, abs=1e-12)


def test_weight_sum_is_one_within_1e12(tmp_path):
    body = (
        "1\nenergy: -12.34 Eh gfn2\nH 0.0 0.0 0.0\n"
        "1\nenergy: -12.30 Eh gfn2\nH 0.0 1.0 0.0\n"
        "1\nenergy: -11.00 Eh gfn2\nH 1.0 0.0 0.0\n"
    )
    frames = parse_ensemble(_write_xyz(tmp_path, body))
    w = boltzmann_weights(frames, temperature_K=500.0)
    assert abs(sum(w) - 1.0) < 1e-12


# --- attack: comment-line number ambiguity -----------------------------------

def test_first_frame_index_is_not_read_as_energy(tmp_path):
    # CREST writes a frame index before the energy in some variants;
    # the parser must not mistake the index for the energy (Eh).
    p = _write_xyz(
        tmp_path,
        "1\n2 -12.34 Eh\nH 0.0 0.0 0.0\n",
    )
    frames = parse_ensemble(p)
    assert frames[0]["energy_Eh"] == pytest.approx(-12.34, abs=1e-9)


# --- attack: NaN propagates silently through relative_kcal -------------------

def test_nan_energy_is_refused_not_propagated(tmp_path):
    p = _write_xyz(
        tmp_path,
        "1\nenergy: nan Eh\nH 0.0 0.0 0.0\n"
        "1\nenergy: -12.34 Eh\nH 0.0 0.0 0.0\n",
    )
    frames = parse_ensemble(p)
    with pytest.raises(ValueError):
        relative_kcal(frames)


# --- attack: collapse_degeneracies groups across methods ---------------------

def test_collapse_does_not_merge_frames_of_different_methods(tmp_path):
    body = (
        "1\nenergy: -12.3400 Eh gfn2\nH 0.0 0.0 0.0\n"
        "1\nenergy: -12.3399 Eh r2scan-3c\nH 0.0 0.0 0.0\n"
    )
    frames = parse_ensemble(_write_xyz(tmp_path, body))
    groups = collapse_degeneracies(frames)
    # 0.0001 Eh apart = 0.0628 kcal/mol > the 0.01 kcal/mol tolerance, so no
    # merge regardless; but a mixed-method group must never exist either way:
    for g in groups:
        assert g.get("method") == "gfn2" or g.get("method") == "r2scan-3c"
    # The mixed-level refusal of boltzmann_weights must survive collapsing:
    with pytest.raises(ValueError):
        boltzmann_weights(groups, temperature_K=298.15)


# --- attack: unicode look-alike in a comment ---------------------------------

def test_unicode_minus_sign_is_not_silently_ignored(tmp_path):
    # U+2212 MINUS SIGN instead of ASCII hyphen: must not read as +12.34.
    p = _write_xyz(
        tmp_path,
        "1\nenergy: −12.34 Eh\nH 0.0 0.0 0.0\n",
    )
    with pytest.raises(ValueError):
        parse_ensemble(p)
