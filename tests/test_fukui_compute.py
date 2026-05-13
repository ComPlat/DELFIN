"""Unit tests for delfin.fukui — Fukui formula + charge parsing."""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

import pytest

from delfin import fukui


# ---------------------------------------------------------------------------
# compute_fukui_from_charges
# ---------------------------------------------------------------------------

def test_compute_fukui_sentinel_values():
    """Verify formula against hand-computed sentinel values."""
    q_neutral = [0.10, -0.20, 0.05, 0.05]
    q_anion = [-0.05, -0.30, -0.02, -0.03]
    q_cation = [0.25, -0.10, 0.15, 0.20]

    result = fukui.compute_fukui_from_charges(q_neutral, q_anion, q_cation)

    # f_plus = q_neutral - q_anion
    assert result["f_plus"] == pytest.approx([0.15, 0.10, 0.07, 0.08])
    # f_minus = q_cation - q_neutral
    assert result["f_minus"] == pytest.approx([0.15, 0.10, 0.10, 0.15])
    # f_zero = (q_cation - q_anion) / 2
    assert result["f_zero"] == pytest.approx([0.15, 0.10, 0.085, 0.115])


def test_compute_fukui_formaldehyde_qualitative():
    """Replicate OPI-notebook qualitative finding (HCHO):
    f_plus maximum on C, f_minus maximum on O, f_zero maximum on O.

    Numbers are synthetic but ordered to match the chemistry the OPI
    notebook ends with.
    """
    # Atom order: C, O, H, H
    # LUMO is C-centric (anion deposits charge there) → big f_plus on C.
    # HOMO is O-centric (cation removes charge there) → big f_minus on O,
    # and because the cation swing on O dominates the f_zero average, O also
    # carries the largest radical Fukui index.
    q_neutral = [0.10, -0.30, 0.10, 0.10]
    q_anion = [-0.20, -0.40, 0.05, 0.05]
    q_cation = [0.20, 0.50, 0.15, 0.15]

    res = fukui.compute_fukui_from_charges(q_neutral, q_anion, q_cation)

    assert res["f_plus"][0] == max(res["f_plus"]), "f_plus max should be on C"
    assert res["f_minus"][1] == max(res["f_minus"]), "f_minus max should be on O"
    assert res["f_zero"][1] == max(res["f_zero"]), "f_zero max should be on O"


def test_compute_fukui_length_mismatch_raises():
    with pytest.raises(ValueError, match="length mismatch"):
        fukui.compute_fukui_from_charges([0.1, 0.2], [0.0], [0.3, 0.4])


# ---------------------------------------------------------------------------
# read_charges — parsing real ORCA-formatted Mulliken/Loewdin blocks
# ---------------------------------------------------------------------------

_FAKE_ORCA_MULLIKEN_CLOSED = textwrap.dedent("""\
    Some ORCA preamble blah blah
    ...

    -----------------------
    MULLIKEN ATOMIC CHARGES
    -----------------------
       0 C    :   -0.123456
       1 O    :   -0.234567
       2 H    :    0.178900
       3 H    :    0.179123
    Sum of atomic charges:   -0.0000000

    -----------------------
    Some other section
    """)

_FAKE_ORCA_MULLIKEN_OPEN = textwrap.dedent("""\
    -----------------------
    MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS
    -----------------------
       0 C    :   -0.300000    0.500000
       1 O    :   -0.100000    0.400000
       2 H    :    0.050000    0.050000
       3 H    :    0.050000    0.050000
    Sum of atomic charges:   -0.3000000

    """)

_FAKE_ORCA_LOEWDIN = textwrap.dedent("""\
    ----------------------
    LOEWDIN ATOMIC CHARGES
    ----------------------
       0 C    :    0.111000
       1 O    :   -0.222000
       2 H    :    0.055500
       3 H    :    0.055500

    """)


def _write_fake_out(tmp_path: Path, contents: str) -> Path:
    p = tmp_path / "fake.out"
    p.write_text(contents, encoding="utf-8")
    return p


def test_read_charges_mulliken_closed_shell(tmp_path):
    p = _write_fake_out(tmp_path, _FAKE_ORCA_MULLIKEN_CLOSED)
    charges = fukui.read_charges(p, scheme="mulliken")
    assert charges == pytest.approx([-0.123456, -0.234567, 0.178900, 0.179123])


def test_read_charges_mulliken_open_shell(tmp_path):
    p = _write_fake_out(tmp_path, _FAKE_ORCA_MULLIKEN_OPEN)
    charges = fukui.read_charges(p, scheme="mulliken")
    # Spin populations must NOT bleed into charges
    assert charges == pytest.approx([-0.3, -0.1, 0.05, 0.05])


def test_read_charges_loewdin(tmp_path):
    p = _write_fake_out(tmp_path, _FAKE_ORCA_LOEWDIN)
    charges = fukui.read_charges(p, scheme="loewdin")
    assert charges == pytest.approx([0.111, -0.222, 0.0555, 0.0555])


def test_read_atoms_and_charges_returns_symbols(tmp_path):
    p = _write_fake_out(tmp_path, _FAKE_ORCA_MULLIKEN_CLOSED)
    symbols, charges = fukui.read_atoms_and_charges(p, scheme="mulliken")
    assert symbols == ["C", "O", "H", "H"]
    assert charges == pytest.approx([-0.123456, -0.234567, 0.178900, 0.179123])


def test_read_charges_unknown_scheme_raises(tmp_path):
    p = _write_fake_out(tmp_path, _FAKE_ORCA_MULLIKEN_CLOSED)
    with pytest.raises(ValueError, match="scheme must be"):
        fukui.read_charges(p, scheme="mayer")


def test_read_charges_missing_block_raises(tmp_path):
    p = _write_fake_out(tmp_path, "ORCA preamble without any charge block\n")
    with pytest.raises(ValueError, match="failed to parse"):
        fukui.read_charges(p, scheme="mulliken")


# ---------------------------------------------------------------------------
# JSON serialization + marker
# ---------------------------------------------------------------------------

def test_write_and_load_fukui_result_roundtrip(tmp_path):
    fukui.write_fukui_result_json(
        tmp_path,
        atoms=["C", "O", "H", "H"],
        scheme="mulliken",
        q_neutral=[0.1, -0.2, 0.05, 0.05],
        q_anion=[-0.05, -0.3, -0.02, -0.03],
        q_cation=[0.25, -0.1, 0.15, 0.2],
        fukui={
            "f_plus":  [0.15, 0.10, 0.07, 0.08],
            "f_minus": [0.15, 0.10, 0.10, 0.15],
            "f_zero":  [0.15, 0.10, 0.085, 0.115],
        },
        orca_settings={"functional": "B3LYP", "basis": "def2-SVP"},
        geometry_origin="user_xyz",
    )

    loaded = fukui.load_fukui_result(tmp_path)
    assert loaded["atoms"] == ["C", "O", "H", "H"]
    assert loaded["scheme"] == "mulliken"
    assert loaded["f_plus"] == pytest.approx([0.15, 0.10, 0.07, 0.08])
    assert loaded["orca_settings"]["functional"] == "B3LYP"
    assert loaded["geometry_origin"] == "user_xyz"
    assert "timestamp_utc" in loaded


def test_write_fukui_result_includes_extra(tmp_path):
    fukui.write_fukui_result_json(
        tmp_path,
        atoms=["C"],
        scheme="loewdin",
        q_neutral=[0.0], q_anion=[0.0], q_cation=[0.0],
        fukui={"f_plus": [0.0], "f_minus": [0.0], "f_zero": [0.0]},
        orca_settings={},
        geometry_origin="opt",
        extra={"job_label": "fukui_test_xyz"},
    )
    loaded = fukui.load_fukui_result(tmp_path)
    assert loaded["job_label"] == "fukui_test_xyz"


def test_marker_lifecycle(tmp_path):
    assert not fukui.is_fukui_dir(tmp_path)
    fukui.write_marker(tmp_path)
    assert fukui.is_fukui_dir(tmp_path)
    assert (tmp_path / fukui.FUKUI_MARKER).exists()
