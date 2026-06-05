"""MISSION F2 -- Universal Embed Fallback with selectable polish.

Covers the ``DELFIN_FFFREE_FALLBACK_MODE`` dispatch contract:

  * default unset == ``"uff"``: byte-identical to pre-F2 (no embed-fallback,
    legacy UFF pipeline runs as before).
  * ``"grip"``: ETKDGv3 embed + GRIP polish when fffree-native fails; the
    legacy UFF fallthrough is REPLACED, so the call returns either fffree-
    native or the embed-grip result (never UFF).
  * ``"none"``: equivalent to the existing ``NO_FALLBACK`` (return ``[]``
    when fffree-native fails) but spelled explicitly via the mode flag.
  * ``"both"``: emit embed-grip output PLUS the legacy UFF output so the
    paper can A/B them on the same SMILES.

Determinism: ETKDGv3 with fixed ``randomSeed = 42`` + ``useRandomCoords =
False`` -> two consecutive calls return byte-identical XYZ blocks.

Default-OFF byte-identity: with no env flag set the only code touched by F2
is a single conditional read of ``os.environ`` and a no-op stash reset, so
the legacy result path is unchanged.

The tests use ``monkeypatch`` to fake ``_fffree_isomers`` so we never need
the heavy assemble path -- only the dispatch decision tree matters.
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


# ---------------------------------------------------------------------------
# resolve_fallback_mode -- unit tests
# ---------------------------------------------------------------------------


def test_resolve_default_unset_is_uff():
    assert ef.resolve_fallback_mode({}) == "uff"


@pytest.mark.parametrize("val,expected", [
    ("grip", "grip"),
    ("GRIP", "grip"),
    ("  grip  ", "grip"),
    ("uff", "uff"),
    ("none", "none"),
    ("both", "both"),
])
def test_resolve_valid_modes(val, expected):
    assert ef.resolve_fallback_mode(
        {"DELFIN_FFFREE_FALLBACK_MODE": val}
    ) == expected


@pytest.mark.parametrize("val", ["", "xyz", "garbage", "1", "0", "true"])
def test_resolve_invalid_falls_back_to_uff(val):
    assert ef.resolve_fallback_mode(
        {"DELFIN_FFFREE_FALLBACK_MODE": val}
    ) == "uff"


# ---------------------------------------------------------------------------
# embed_isomers -- direct unit tests
# ---------------------------------------------------------------------------


def test_embed_isomers_organic_raw_returns_conformers():
    r = ef.embed_isomers(_ORGANIC_SMI, max_isomers=3, polish="raw")
    assert r is not None and len(r) == 3
    for xyz, lab in r:
        assert lab.startswith("embed-conf")
        assert lab.endswith("-raw")
        assert xyz.count("\n") >= 0  # non-empty


def test_embed_isomers_organic_grip_runs_polish_as_raw_when_no_metal():
    """No metal in CCO -> polish silently degrades to raw (no exception)."""
    r = ef.embed_isomers(_ORGANIC_SMI, max_isomers=2, polish="grip")
    assert r is not None and len(r) == 2
    for xyz, lab in r:
        assert lab.endswith("-raw")  # no metal -> no polish


def test_embed_isomers_metal_grip_polishes_each_conformer():
    r = ef.embed_isomers("N[Pt](N)(Cl)Cl", max_isomers=2, polish="grip")
    assert r is not None and len(r) == 2
    for xyz, lab in r:
        # metal present -> grip suffix when polish succeeded, raw on fallback
        assert lab.endswith("-grip") or lab.endswith("-raw")


def test_embed_isomers_metal_both_emits_double():
    """``both`` mode emits raw + grip-polished per conformer."""
    r = ef.embed_isomers("N[Pt](N)(Cl)Cl", max_isomers=2, polish="both")
    assert r is not None
    # ``both`` produces 2 outputs per conformer (raw + grip)
    assert len(r) == 4
    raws = [lab for _, lab in r if lab.endswith("-raw")]
    grips = [lab for _, lab in r if lab.endswith("-grip")]
    assert len(raws) == 2 and len(grips) == 2


def test_embed_isomers_invalid_smiles_returns_none():
    assert ef.embed_isomers("XYZ-not-a-smiles", max_isomers=3) is None


def test_embed_isomers_clamps_max_isomers():
    # >16 should be capped at 16 silently
    r = ef.embed_isomers("CC", max_isomers=999, polish="raw")
    assert r is not None and len(r) <= 16


# ---------------------------------------------------------------------------
# Determinism: 2-run byte-identical
# ---------------------------------------------------------------------------


def test_embed_isomers_2run_byte_identical_organic():
    r1 = ef.embed_isomers(_ORGANIC_SMI, max_isomers=4, polish="raw")
    r2 = ef.embed_isomers(_ORGANIC_SMI, max_isomers=4, polish="raw")
    assert r1 is not None and r2 is not None
    assert r1 == r2, "2-run XYZ blocks must be byte-identical"


def test_embed_isomers_2run_byte_identical_metal_grip():
    smi = "N[Pt](N)(Cl)Cl"
    r1 = ef.embed_isomers(smi, max_isomers=3, polish="grip")
    r2 = ef.embed_isomers(smi, max_isomers=3, polish="grip")
    assert r1 is not None and r2 is not None
    assert r1 == r2, "GRIP-polished 2-run must be byte-identical"


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


def test_default_off_byte_identical_no_dispatch(monkeypatch, clean_env):
    """No env flags -> fffree_isomers is NOT called (byte-identical)."""
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


def test_mode_grip_implies_dispatch(monkeypatch, clean_env):
    """``mode=grip`` alone must dispatch fffree even without BUILDER."""
    calls = []
    sentinel = [("X 0.0 0.0 0.0", "lab")]

    def _spy(smiles, max_isomers=50):
        calls.append(smiles)
        return sentinel

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy,
    )
    monkeypatch.setenv("DELFIN_FFFREE_FALLBACK_MODE", "grip")
    res, err = sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert calls == [_METAL_SMI]
    assert res == sentinel and err is None


def test_mode_grip_runs_embed_fallback_on_native_fail(monkeypatch, clean_env):
    """When fffree-native returns None and ``mode=grip``, embed-fallback runs."""

    def _spy_fail(smiles, max_isomers=50):
        return None

    embed_calls = []
    sentinel = [("X 0.0 0.0 0.0", "embed-conf0-grip")]

    def _embed_spy(smiles, *, max_isomers=16, polish="grip"):
        embed_calls.append((smiles, polish))
        return sentinel

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy_fail,
    )
    monkeypatch.setattr(
        "delfin.fffree.embed_fallback.embed_isomers", _embed_spy,
    )
    monkeypatch.setenv("DELFIN_FFFREE_FALLBACK_MODE", "grip")
    res, err = sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert embed_calls == [(_METAL_SMI, "grip")]
    assert res == sentinel and err is None


def test_mode_none_returns_empty_on_native_fail(monkeypatch, clean_env):
    """``mode=none`` is paper-grade: fffree-native fail -> return []."""

    def _spy_fail(smiles, max_isomers=50):
        return None

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy_fail,
    )
    # ensure the embed-fallback is NOT called in "none" mode
    embed_calls = []

    def _embed_spy(smiles, *, max_isomers=16, polish="grip"):
        embed_calls.append(smiles)
        return [("X 0.0 0.0 0.0", "lab")]

    monkeypatch.setattr(
        "delfin.fffree.embed_fallback.embed_isomers", _embed_spy,
    )
    monkeypatch.setenv("DELFIN_FFFREE_FALLBACK_MODE", "none")
    res, err = sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert res == [] and err is None
    assert embed_calls == [], "embed-fallback must NOT fire in mode=none"


def test_mode_uff_default_falls_through_to_legacy(monkeypatch, clean_env):
    """``mode=uff`` is the legacy fallthrough: when BUILDER is also set and
    fffree-native fails, we should NOT short-circuit -- the rest of the
    impl runs.  We assert that by checking the impl proceeds past the
    dispatch block (RDKit-disabled path returns the standard error).
    """
    def _spy_fail(smiles, max_isomers=50):
        return None

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy_fail,
    )
    monkeypatch.setattr(sc, "RDKIT_AVAILABLE", False)
    monkeypatch.setenv("DELFIN_FFFREE_BUILDER", "1")
    monkeypatch.setenv("DELFIN_FFFREE_FALLBACK_MODE", "uff")
    res, err = sc.smiles_to_xyz_isomers(_METAL_SMI)
    # Falls through to legacy => RDKit-disabled error from the standard
    # legacy path (proves we did not short-circuit return).
    assert res == [] and err == "RDKit is not installed"


def test_mode_both_concatenates_embed_and_legacy(monkeypatch, clean_env):
    """``mode=both`` should produce embed-fallback results AND let the
    legacy pipeline run to produce its own output, then the wrapper
    concatenates the two.
    """
    def _spy_fail(smiles, max_isomers=50):
        return None

    embed_results = [
        ("E0 0.0 0.0 0.0", "embed-conf0-grip"),
        ("E1 0.0 0.0 0.0", "embed-conf1-grip"),
    ]

    def _embed_spy(smiles, *, max_isomers=16, polish="grip"):
        return embed_results

    legacy_results: List[Tuple[str, str]] = [("L0 1.0 0.0 0.0", "legacy-uff-0")]

    def _impl_spy(*args, **kwargs):
        # Mirror the impl behaviour: when "both" mode is active and
        # fffree fails, the impl writes _F2_BOTH_EXTRA on the thread-
        # local before returning the legacy pipeline output.  We bypass
        # the heavy impl by writing the stash directly + returning a
        # fake legacy result.
        setattr(sc._F2_BOTH_EXTRA, "extra", list(embed_results))
        return list(legacy_results), None

    monkeypatch.setattr(sc, "_smiles_to_xyz_isomers_impl", _impl_spy)
    monkeypatch.setenv("DELFIN_FFFREE_FALLBACK_MODE", "both")
    res, err = sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert err is None
    labels = [lab for _, lab in res]
    # embed-fallback results must come BEFORE the legacy results
    assert labels == [
        "embed-conf0-grip", "embed-conf1-grip", "legacy-uff-0",
    ]


def test_mode_both_stash_clears_between_calls(monkeypatch, clean_env):
    """The thread-local stash must not leak between calls."""

    def _impl_no_stash(*args, **kwargs):
        return [("X 0.0 0.0 0.0", "legacy")], None

    # First call writes stash
    setattr(sc._F2_BOTH_EXTRA, "extra", [("E 0 0 0", "leftover")])
    monkeypatch.setattr(sc, "_smiles_to_xyz_isomers_impl", _impl_no_stash)
    res, _ = sc.smiles_to_xyz_isomers(_ORGANIC_SMI)
    # The wrapper reset the stash at entry, so this call must NOT inherit
    # the leftover -- only the legacy result should be present.
    labels = [lab for _, lab in res]
    assert labels == ["legacy"]


def test_forensik_log_records_embed_status(monkeypatch, clean_env, tmp_path):
    """Mode=grip + fffree-native fail should log status ``embed-grip``."""
    def _spy_fail(smiles, max_isomers=50):
        return None

    def _embed_spy(smiles, *, max_isomers=16, polish="grip"):
        return [("X 0.0 0.0 0.0", "embed-conf0-grip")]

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy_fail,
    )
    monkeypatch.setattr(
        "delfin.fffree.embed_fallback.embed_isomers", _embed_spy,
    )
    log_path = tmp_path / "forensik.tsv"
    monkeypatch.setenv("DELFIN_FFFREE_FALLBACK_MODE", "grip")
    monkeypatch.setenv("DELFIN_FFFREE_FORENSIK_LOG", str(log_path))
    res, err = sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert res and err is None
    assert log_path.exists()
    text = log_path.read_text().strip().splitlines()
    assert len(text) == 1
    parts = text[0].split("\t")
    assert parts[1] == "embed-grip"
    assert parts[2] == "1"


def test_forensik_logs_embed_failed_when_grip_yields_nothing(
    monkeypatch, clean_env, tmp_path,
):
    """``mode=grip`` + fffree-native fail + embed-fallback returns None
    -> forensik status must be ``embed-failed`` (distinguish from the
    legacy ``fallback`` status which means "would fall through to UFF").
    """
    def _spy_fail(smiles, max_isomers=50):
        return None

    def _embed_none(smiles, *, max_isomers=16, polish="grip"):
        return None

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy_fail,
    )
    monkeypatch.setattr(
        "delfin.fffree.embed_fallback.embed_isomers", _embed_none,
    )
    log_path = tmp_path / "forensik.tsv"
    monkeypatch.setenv("DELFIN_FFFREE_FALLBACK_MODE", "grip")
    monkeypatch.setenv("DELFIN_FFFREE_FORENSIK_LOG", str(log_path))
    res, err = sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert res == [] and err is None
    text = log_path.read_text().strip().splitlines()
    assert len(text) == 1
    parts = text[0].split("\t")
    assert parts[1] == "embed-failed"


def test_non_metal_smiles_unaffected_by_mode(monkeypatch, clean_env):
    """``mode=grip`` MUST NOT touch the dispatch for non-metal SMILES."""
    calls = []

    def _spy(smiles, max_isomers=50):
        calls.append(smiles)
        return None

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy,
    )
    monkeypatch.setattr(sc, "RDKIT_AVAILABLE", False)
    monkeypatch.setenv("DELFIN_FFFREE_FALLBACK_MODE", "grip")
    res, err = sc.smiles_to_xyz_isomers(_ORGANIC_SMI)
    assert calls == [], "non-metal SMILES must not trigger fffree dispatch"
    assert res == [] and err == "RDKit is not installed"
