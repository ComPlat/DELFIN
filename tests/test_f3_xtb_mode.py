"""MISSION F3 -- xtb post-processing as 5th selectable fallback mode.

Covers the ``DELFIN_FFFREE_FALLBACK_MODE`` dispatch contract extensions:

  * ``"xtb"``  -- ETKDG embed + GFN2-xTB relaxation per conformer.
  * ``"all"``  -- emit grip + uff + xtb variants per conformer (A/B/C).
  * existing modes (``uff``/``grip``/``none``/``both``) still parse the same.

The xtb branch is gracefully skipped when xtb is not installed -- the
output then carries the ``-raw`` suffix instead of ``-xtb`` so callers
can detect partial-polish without crashing.

Determinism: ETKDGv3 seed=42 + xtb ``--norestart`` + ``OMP_NUM_THREADS=1``
-> two consecutive calls return byte-identical XYZ blocks.

Default-OFF byte-identity: with no env flag set ``resolve_fallback_mode``
still returns ``"uff"`` so the legacy pipeline is unchanged.
"""
from __future__ import annotations

import os
from typing import List, Tuple

import numpy as np
import pytest

from delfin import smiles_converter as sc
from delfin.fffree import embed_fallback as ef


_METAL_SMI = "[OH2+][Cd-4]([Cl])([Cl])[OH2+]"
_ORGANIC_SMI = "CCO"
_SIMPLE_TMC = "N[Pt](N)(Cl)Cl"


# ---------------------------------------------------------------------------
# resolve_fallback_mode -- F3 mode extensions
# ---------------------------------------------------------------------------


def test_resolve_default_still_uff_after_f3():
    """F3 must not change the default: unset env still maps to ``uff``."""
    assert ef.resolve_fallback_mode({}) == "uff"


@pytest.mark.parametrize("val,expected", [
    ("xtb", "xtb"),
    ("XTB", "xtb"),
    ("  xtb  ", "xtb"),
    ("all", "all"),
    ("ALL", "all"),
    ("  all  ", "all"),
])
def test_resolve_new_modes(val, expected):
    assert ef.resolve_fallback_mode(
        {"DELFIN_FFFREE_FALLBACK_MODE": val}
    ) == expected


def test_valid_modes_constant_contains_new_modes():
    """The published constant must include the F3 additions."""
    assert "xtb" in ef._VALID_MODES
    assert "all" in ef._VALID_MODES
    # F2 modes still present
    for m in ("grip", "uff", "none", "both"):
        assert m in ef._VALID_MODES


# ---------------------------------------------------------------------------
# xtb binary locator
# ---------------------------------------------------------------------------


def test_find_xtb_binary_or_none():
    """Returns a path-like string or None; never raises."""
    res = ef._find_xtb_binary()
    assert res is None or (isinstance(res, str) and os.path.isfile(res))


# ---------------------------------------------------------------------------
# xtb XYZ parser
# ---------------------------------------------------------------------------


def test_parse_xtb_xyz_basic():
    text = (
        "3\n"
        "comment\n"
        "C 0.000000 0.000000 0.000000\n"
        "H 1.000000 0.000000 0.000000\n"
        "H 0.000000 1.000000 0.000000\n"
    )
    coords = ef._parse_xtb_xyz(text, 3)
    assert coords is not None
    assert coords.shape == (3, 3)
    np.testing.assert_allclose(coords[1], [1.0, 0.0, 0.0])
    np.testing.assert_allclose(coords[2], [0.0, 1.0, 0.0])


def test_parse_xtb_xyz_count_mismatch_returns_none():
    text = "5\ncomment\nC 0 0 0\nH 1 0 0\n"
    assert ef._parse_xtb_xyz(text, 5) is None


def test_parse_xtb_xyz_truncated_returns_none():
    text = "2\ncomment\nC 0 0 0\n"  # only 1 atom line
    assert ef._parse_xtb_xyz(text, 2) is None


def test_parse_xtb_xyz_garbage_returns_none():
    assert ef._parse_xtb_xyz("not an xyz", 3) is None
    assert ef._parse_xtb_xyz("", 3) is None


# ---------------------------------------------------------------------------
# embed_isomers -- xtb mode
# ---------------------------------------------------------------------------


def test_embed_isomers_xtb_organic_returns_outputs():
    """Organic CCO via xtb-mode returns 2 conformers labelled -xtb or -raw."""
    r = ef.embed_isomers(_ORGANIC_SMI, max_isomers=2, polish="xtb")
    assert r is not None and len(r) == 2
    for xyz, lab in r:
        assert lab.startswith("embed-conf")
        assert lab.endswith("-xtb") or lab.endswith("-raw")
        assert xyz.count("\n") >= 2  # at least 3 atom lines


def test_embed_isomers_xtb_graceful_skip_when_missing(monkeypatch):
    """When xtb binary cannot be found, output gets ``-raw`` suffix (not ``-xtb``)."""
    monkeypatch.setattr(ef, "_find_xtb_binary", lambda: None)
    r = ef.embed_isomers(_ORGANIC_SMI, max_isomers=1, polish="xtb")
    assert r is not None and len(r) == 1
    _, lab = r[0]
    assert lab.endswith("-raw"), (
        f"Expected ``-raw`` when xtb unavailable, got {lab!r}"
    )


def test_embed_isomers_uff_mode_emits_uff_or_raw():
    r = ef.embed_isomers(_ORGANIC_SMI, max_isomers=2, polish="uff")
    assert r is not None and len(r) == 2
    for _, lab in r:
        assert lab.endswith("-uff") or lab.endswith("-raw")


# ---------------------------------------------------------------------------
# embed_isomers -- all mode
# ---------------------------------------------------------------------------


def test_embed_isomers_all_mode_emits_three_per_conformer():
    """``all`` mode emits grip + uff + xtb per conformer (3 outputs each)."""
    n_conf = 2
    r = ef.embed_isomers(_ORGANIC_SMI, max_isomers=n_conf, polish="all")
    assert r is not None
    assert len(r) == 3 * n_conf, (
        f"Expected {3*n_conf} outputs, got {len(r)}"
    )
    # The labels should be ordered grip, uff, xtb per conformer.
    for k in range(n_conf):
        slot_grip = r[3 * k][1]
        slot_uff = r[3 * k + 1][1]
        slot_xtb = r[3 * k + 2][1]
        assert slot_grip.endswith("-grip") or slot_grip.endswith("-raw")
        assert slot_uff.endswith("-uff") or slot_uff.endswith("-raw")
        assert slot_xtb.endswith("-xtb") or slot_xtb.endswith("-raw")


def test_embed_isomers_all_mode_metal_smi_includes_polish_branches():
    """A simple TMC under ``all`` mode produces grip-branch with -grip suffix."""
    r = ef.embed_isomers(_SIMPLE_TMC, max_isomers=1, polish="all")
    assert r is not None and len(r) == 3
    # grip branch should succeed on a real TMC with a metal -- detection
    # only works when metal+donors exist.
    suffixes = [lab.rsplit("-", 1)[-1] for _, lab in r]
    # At least one of {grip, uff, raw} per slot.
    assert any(s == "grip" or s == "raw" for s in suffixes[:1])
    assert any(s == "uff" or s == "raw" for s in suffixes[1:2])
    assert any(s == "xtb" or s == "raw" for s in suffixes[2:3])


# ---------------------------------------------------------------------------
# Determinism: 2-run byte-identical
# ---------------------------------------------------------------------------


def test_embed_isomers_2run_byte_identical_xtb():
    """xtb mode 2-run on the same SMILES is byte-identical."""
    smi = _ORGANIC_SMI
    r1 = ef.embed_isomers(smi, max_isomers=2, polish="xtb")
    r2 = ef.embed_isomers(smi, max_isomers=2, polish="xtb")
    assert r1 is not None and r2 is not None
    assert r1 == r2, "xtb 2-run output must be byte-identical"


def test_embed_isomers_2run_byte_identical_all():
    """``all`` mode 2-run on the same SMILES is byte-identical."""
    smi = _ORGANIC_SMI
    r1 = ef.embed_isomers(smi, max_isomers=2, polish="all")
    r2 = ef.embed_isomers(smi, max_isomers=2, polish="all")
    assert r1 is not None and r2 is not None
    assert r1 == r2, "all 2-run output must be byte-identical"


# ---------------------------------------------------------------------------
# Convenience wrappers
# ---------------------------------------------------------------------------


def test_xtb_embed_fallback_wrapper_delegates():
    r = ef.xtb_embed_fallback(_ORGANIC_SMI, max_isomers=1)
    assert r is not None and len(r) == 1
    _, lab = r[0]
    assert lab.endswith("-xtb") or lab.endswith("-raw")


def test_all_embed_fallback_wrapper_delegates():
    r = ef.all_embed_fallback(_ORGANIC_SMI, max_isomers=1)
    assert r is not None and len(r) == 3
    suffixes = [lab.rsplit("-", 1)[-1] for _, lab in r]
    assert suffixes[0] in ("grip", "raw")
    assert suffixes[1] in ("uff", "raw")
    assert suffixes[2] in ("xtb", "raw")


# ---------------------------------------------------------------------------
# Dispatch contract -- smiles_to_xyz_isomers integration
# ---------------------------------------------------------------------------


@pytest.fixture
def clean_env(monkeypatch):
    """Strip every fffree env flag for a known starting point."""
    for var in (
        "DELFIN_FFFREE_BUILDER",
        "DELFIN_FFFREE_NO_FALLBACK",
        "DELFIN_FFFREE_PURE_TRACK3",
        "DELFIN_FFFREE_FALLBACK_MODE",
        "DELFIN_FFFREE_FORENSIK_LOG",
    ):
        monkeypatch.delenv(var, raising=False)
    yield


def test_mode_xtb_implies_dispatch(monkeypatch, clean_env):
    """``mode=xtb`` alone must dispatch fffree even without BUILDER."""
    calls = []
    sentinel = [("X 0.0 0.0 0.0", "lab")]

    def _spy(smiles, max_isomers=50):
        calls.append(smiles)
        return sentinel

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy,
    )
    monkeypatch.setenv("DELFIN_FFFREE_FALLBACK_MODE", "xtb")
    res, err = sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert calls == [_METAL_SMI]
    assert res == sentinel and err is None


def test_mode_xtb_runs_embed_fallback_on_native_fail(monkeypatch, clean_env):
    """``mode=xtb`` + fffree-native fail -> embed_isomers(..., polish='xtb')."""

    def _spy_fail(smiles, max_isomers=50):
        return None

    embed_calls = []
    sentinel = [("X 0.0 0.0 0.0", "embed-conf0-xtb")]

    def _embed_spy(smiles, *, max_isomers=16, polish="grip"):
        embed_calls.append((smiles, polish))
        return sentinel

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy_fail,
    )
    monkeypatch.setattr(
        "delfin.fffree.embed_fallback.embed_isomers", _embed_spy,
    )
    monkeypatch.setenv("DELFIN_FFFREE_FALLBACK_MODE", "xtb")
    res, err = sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert embed_calls == [(_METAL_SMI, "xtb")]
    assert res == sentinel and err is None


def test_mode_all_runs_embed_fallback_with_all_polish(monkeypatch, clean_env):
    """``mode=all`` + fffree-native fail -> embed_isomers(..., polish='all')."""

    def _spy_fail(smiles, max_isomers=50):
        return None

    embed_calls = []
    sentinel = [
        ("X 0.0 0.0 0.0", "embed-conf0-grip"),
        ("X 0.0 0.0 0.0", "embed-conf0-uff"),
        ("X 0.0 0.0 0.0", "embed-conf0-xtb"),
    ]

    def _embed_spy(smiles, *, max_isomers=16, polish="grip"):
        embed_calls.append((smiles, polish))
        return sentinel

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy_fail,
    )
    monkeypatch.setattr(
        "delfin.fffree.embed_fallback.embed_isomers", _embed_spy,
    )
    monkeypatch.setenv("DELFIN_FFFREE_FALLBACK_MODE", "all")
    res, err = sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert embed_calls == [(_METAL_SMI, "all")]
    assert res == sentinel and err is None
    # The all-mode dispatch returns the embed result directly (no UFF
    # concat) since the embed output already contains a UFF branch.
    labels = [lab for _, lab in res]
    assert labels == [
        "embed-conf0-grip", "embed-conf0-uff", "embed-conf0-xtb",
    ]


def test_default_off_byte_identical_no_dispatch_after_f3(monkeypatch, clean_env):
    """F3 must not change default-OFF byte-identity: no env flags = no dispatch."""
    calls = []

    def _spy(smiles, max_isomers=50):
        calls.append(smiles)
        return None

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy,
    )
    monkeypatch.setattr(sc, "RDKIT_AVAILABLE", False)
    res, err = sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert res == [] and err == "RDKit is not installed"
    assert calls == [], "fffree dispatch must NOT fire when all flags OFF"


def test_forensik_log_records_embed_xtb_status(monkeypatch, clean_env, tmp_path):
    """``mode=xtb`` + fffree-native fail -> log status ``embed-xtb``."""

    def _spy_fail(smiles, max_isomers=50):
        return None

    def _embed_spy(smiles, *, max_isomers=16, polish="grip"):
        return [("X 0.0 0.0 0.0", "embed-conf0-xtb")]

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy_fail,
    )
    monkeypatch.setattr(
        "delfin.fffree.embed_fallback.embed_isomers", _embed_spy,
    )
    log_path = tmp_path / "forensik.tsv"
    monkeypatch.setenv("DELFIN_FFFREE_FALLBACK_MODE", "xtb")
    monkeypatch.setenv("DELFIN_FFFREE_FORENSIK_LOG", str(log_path))
    res, err = sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert res and err is None
    assert log_path.exists()
    text = log_path.read_text().strip().splitlines()
    assert len(text) == 1
    parts = text[0].split("\t")
    assert parts[1] == "embed-xtb"


def test_forensik_log_records_embed_all_status(monkeypatch, clean_env, tmp_path):
    """``mode=all`` + fffree-native fail -> log status ``embed-all``."""

    def _spy_fail(smiles, max_isomers=50):
        return None

    def _embed_spy(smiles, *, max_isomers=16, polish="grip"):
        return [
            ("X 0.0 0.0 0.0", "embed-conf0-grip"),
            ("X 0.0 0.0 0.0", "embed-conf0-uff"),
            ("X 0.0 0.0 0.0", "embed-conf0-xtb"),
        ]

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy_fail,
    )
    monkeypatch.setattr(
        "delfin.fffree.embed_fallback.embed_isomers", _embed_spy,
    )
    log_path = tmp_path / "forensik.tsv"
    monkeypatch.setenv("DELFIN_FFFREE_FALLBACK_MODE", "all")
    monkeypatch.setenv("DELFIN_FFFREE_FORENSIK_LOG", str(log_path))
    res, err = sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert res and err is None
    parts = log_path.read_text().strip().split("\t")
    assert parts[1] == "embed-all"
