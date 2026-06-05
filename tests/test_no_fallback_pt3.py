"""MISSION A5 — NO_FALLBACK + PT3: kill the legacy-UFF fallback for metals.

Covers the dispatch behaviour of ``smiles_to_xyz_isomers`` when the new
``DELFIN_FFFREE_NO_FALLBACK`` env-flag (Mission D1 2026-06-05: PT3 no longer
auto-implies NO_FALLBACK -- legacy 2792332-style fallback is restored under
PT3 alone; NO_FALLBACK must be set explicitly for the "honest fffree" mode)
is set.  Verifies:

  * default-OFF byte-identity (no env, no behavioural change)
  * NO_FALLBACK=1 alone implicitly enables fffree dispatch (even without BUILDER)
  * fffree success returns the native isomer list unchanged
  * fffree failure returns ``([], None)`` (no legacy fallthrough)
  * forensik log file receives one ``smi\\tstatus\\tcount`` line per call
  * 2-run determinism with all flags ON
  * non-metal SMILES are not touched by the dispatch (no flag consultation)

Module-level mock-friendly: only the converter dispatch path is exercised, no
RDKit/OB heavy lifting.
"""
from __future__ import annotations

import os

import pytest

from delfin import smiles_converter as sc


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_METAL_SMI = "[OH2+][Cd-4]([Cl])([Cl])[OH2+]"   # any explicit-metal SMILES


@pytest.fixture
def clean_env(monkeypatch):
    """Strip every NO_FALLBACK / PT3 / BUILDER flag for a known starting point."""
    for var in (
        "DELFIN_FFFREE_BUILDER",
        "DELFIN_FFFREE_NO_FALLBACK",
        "DELFIN_FFFREE_PURE_TRACK3",
        "DELFIN_FFFREE_FORENSIK_LOG",
    ):
        monkeypatch.delenv(var, raising=False)
    yield


# ---------------------------------------------------------------------------
# Determinism + byte-identity (default OFF)
# ---------------------------------------------------------------------------


def test_no_flags_set_does_not_dispatch_fffree(monkeypatch, clean_env):
    """With NO env flags set, the fffree dispatch must be skipped entirely.

    We patch ``_fffree_isomers`` so that, if it WERE called, the test
    would observe a side-effect (counter increment).  No call ⇒ off path.
    """
    calls = []

    def _spy(smiles, max_isomers=50):
        calls.append(smiles)
        return None

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy,
    )
    # Short-circuit RDKit availability so the rest of the routine bails fast.
    monkeypatch.setattr(sc, "RDKIT_AVAILABLE", False)
    res, err = sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert res == [] and err == "RDKit is not installed"
    assert calls == [], "fffree dispatch must NOT fire when all flags OFF"


def test_no_fallback_alone_implicitly_enables_dispatch(monkeypatch, clean_env):
    """``DELFIN_FFFREE_NO_FALLBACK=1`` must dispatch fffree even when
    ``FFFREE_BUILDER`` is unset (the auto-implication added by MISSION A5).
    """
    sentinel = [("X 0.0 0.0 0.0", "lab")]
    calls = []

    def _spy(smiles, max_isomers=50):
        calls.append(smiles)
        return list(sentinel)

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy,
    )
    monkeypatch.setattr(sc, "RDKIT_AVAILABLE", True)
    monkeypatch.setattr(sc, "contains_metal", lambda _smi: True)
    monkeypatch.setenv("DELFIN_FFFREE_NO_FALLBACK", "1")

    res, err = sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert err is None
    assert res == sentinel
    assert calls == [_METAL_SMI], "NO_FALLBACK=1 must implicitly enable dispatch"


def test_pure_track3_alone_does_not_dispatch_fffree(monkeypatch, clean_env):
    """``DELFIN_FFFREE_PURE_TRACK3=1`` alone (without BUILDER or NO_FALLBACK)
    MUST NOT auto-dispatch fffree.  Mission D1 (2026-06-05): the A5
    auto-implication broke the voll-pool quality on Internals axes by
    forcing fffree-native dispatch on PT3-only runs.  The pre-A5 contract
    is that PT3 is a downstream flag (sort/dedup/conformer) and only
    BUILDER or NO_FALLBACK actually dispatches the fffree builder.
    """
    calls = []

    def _spy(smiles, max_isomers=50):
        calls.append(smiles)
        return [("X 0.0 0.0 0.0", "lab")]

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy,
    )
    monkeypatch.setattr(sc, "RDKIT_AVAILABLE", True)
    monkeypatch.setattr(sc, "contains_metal", lambda _smi: True)
    monkeypatch.setattr(sc, "_probe_hapto_groups_from_smiles", lambda _smi: [])
    monkeypatch.setattr(sc, "_prepare_mol_for_embedding", lambda *_a, **_k: None)
    monkeypatch.setenv("DELFIN_FFFREE_PURE_TRACK3", "1")

    try:
        sc.smiles_to_xyz_isomers(_METAL_SMI)
    except Exception:
        pass
    assert calls == [], (
        "PT3=1 alone must NOT dispatch fffree -- legacy pipeline owns "
        "the metal-SMILES, matching 2792332-style behaviour"
    )


def test_no_fallback_blocks_legacy_on_fffree_failure(monkeypatch, clean_env):
    """When fffree returns None, NO_FALLBACK=1 must return ([], None) and
    NEVER hand control to the legacy pipeline.
    """

    def _spy(smiles, max_isomers=50):
        return None     # fffree fails to build

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy,
    )
    monkeypatch.setattr(sc, "RDKIT_AVAILABLE", True)
    monkeypatch.setattr(sc, "contains_metal", lambda _smi: True)
    monkeypatch.setenv("DELFIN_FFFREE_NO_FALLBACK", "1")

    # If the no-fallback gate ever lets execution flow past line 26370 the
    # legacy pipeline would call _probe_hapto_groups_from_smiles -- patch it
    # to a sentinel that fails loudly so the test would crash, not silently
    # pass.
    def _legacy_marker(_smi):
        raise AssertionError("LEGACY PIPELINE INVOKED -- NO_FALLBACK violated")

    monkeypatch.setattr(sc, "_probe_hapto_groups_from_smiles", _legacy_marker)

    res, err = sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert res == []
    assert err is None


def test_fallback_allowed_when_neither_flag_set(monkeypatch, clean_env):
    """With BUILDER=1 but no NO_FALLBACK/PT3, fffree-failure MUST fall through
    to the legacy pipeline (existing behaviour preserved).
    """

    def _spy(smiles, max_isomers=50):
        return None     # fffree fails

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy,
    )
    monkeypatch.setattr(sc, "RDKIT_AVAILABLE", True)
    monkeypatch.setattr(sc, "contains_metal", lambda _smi: True)
    monkeypatch.setenv("DELFIN_FFFREE_BUILDER", "1")

    legacy_called = []

    def _legacy_marker(_smi):
        legacy_called.append(_smi)
        return []   # no hapto

    monkeypatch.setattr(sc, "_probe_hapto_groups_from_smiles", _legacy_marker)
    # Short-circuit the rest of the legacy pipeline (we only care that the
    # legacy entry point was reached).
    monkeypatch.setattr(sc, "_prepare_mol_for_embedding",
                        lambda *_a, **_k: None)

    # The full legacy routine will subsequently bail -- we only assert legacy
    # was entered.
    try:
        sc.smiles_to_xyz_isomers(_METAL_SMI)
    except Exception:
        pass
    assert legacy_called and legacy_called[0] == _METAL_SMI, (
        "BUILDER=1 alone must let fffree-failure fall through to legacy"
    )


# ---------------------------------------------------------------------------
# Forensik breadcrumb
# ---------------------------------------------------------------------------


def test_forensik_log_records_native(monkeypatch, clean_env, tmp_path):
    log = tmp_path / "forensik.tsv"
    monkeypatch.setenv("DELFIN_FFFREE_FORENSIK_LOG", str(log))
    monkeypatch.setenv("DELFIN_FFFREE_NO_FALLBACK", "1")
    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers",
        lambda smi, max_isomers=50: [("X 0 0 0", "lab1"), ("X 1 1 1", "lab2")],
    )
    monkeypatch.setattr(sc, "RDKIT_AVAILABLE", True)
    monkeypatch.setattr(sc, "contains_metal", lambda _smi: True)

    sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert log.exists()
    parts = log.read_text().strip().split("\t")
    assert parts[1] == "native"
    assert parts[2] == "2"


def test_forensik_log_records_blocked(monkeypatch, clean_env, tmp_path):
    log = tmp_path / "forensik.tsv"
    monkeypatch.setenv("DELFIN_FFFREE_FORENSIK_LOG", str(log))
    monkeypatch.setenv("DELFIN_FFFREE_NO_FALLBACK", "1")
    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers",
        lambda smi, max_isomers=50: None,
    )
    monkeypatch.setattr(sc, "RDKIT_AVAILABLE", True)
    monkeypatch.setattr(sc, "contains_metal", lambda _smi: True)
    monkeypatch.setattr(sc, "_probe_hapto_groups_from_smiles", lambda _smi: [])

    sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert log.exists()
    parts = log.read_text().strip().split("\t")
    assert parts[1] == "blocked"
    assert parts[2] == "0"


def test_forensik_log_records_fallback(monkeypatch, clean_env, tmp_path):
    log = tmp_path / "forensik.tsv"
    monkeypatch.setenv("DELFIN_FFFREE_FORENSIK_LOG", str(log))
    monkeypatch.setenv("DELFIN_FFFREE_BUILDER", "1")     # NO_FALLBACK off
    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers",
        lambda smi, max_isomers=50: None,
    )
    monkeypatch.setattr(sc, "RDKIT_AVAILABLE", True)
    monkeypatch.setattr(sc, "contains_metal", lambda _smi: True)
    monkeypatch.setattr(sc, "_probe_hapto_groups_from_smiles", lambda _smi: [])
    monkeypatch.setattr(sc, "_prepare_mol_for_embedding",
                        lambda *_a, **_k: None)
    try:
        sc.smiles_to_xyz_isomers(_METAL_SMI)
    except Exception:
        pass
    assert log.exists()
    parts = log.read_text().strip().split("\t")
    assert parts[1] == "fallback"
    assert parts[2] == "0"


def test_forensik_log_not_touched_when_unset(monkeypatch, clean_env, tmp_path):
    """With no DELFIN_FFFREE_FORENSIK_LOG env-var, no file must be created
    even when the dispatch fires.
    """
    log = tmp_path / "would_be_log.tsv"
    monkeypatch.setenv("DELFIN_FFFREE_NO_FALLBACK", "1")
    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers",
        lambda smi, max_isomers=50: [("X 0 0 0", "lab")],
    )
    monkeypatch.setattr(sc, "RDKIT_AVAILABLE", True)
    monkeypatch.setattr(sc, "contains_metal", lambda _smi: True)

    sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert not log.exists()


# ---------------------------------------------------------------------------
# Determinism (2-run): identical input + identical env => identical output
# ---------------------------------------------------------------------------


def test_two_run_determinism_with_no_fallback(monkeypatch, clean_env):
    """Two consecutive calls with identical env must produce the same
    isomer list (the fffree backend itself is the only source of non-
    determinism here; we feed a fixed sentinel)."""
    monkeypatch.setenv("DELFIN_FFFREE_NO_FALLBACK", "1")
    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers",
        lambda smi, max_isomers=50: [("X 0 0 0", "lab")],
    )
    monkeypatch.setattr(sc, "RDKIT_AVAILABLE", True)
    monkeypatch.setattr(sc, "contains_metal", lambda _smi: True)

    r1, _ = sc.smiles_to_xyz_isomers(_METAL_SMI)
    r2, _ = sc.smiles_to_xyz_isomers(_METAL_SMI)
    assert r1 == r2


# ---------------------------------------------------------------------------
# Non-metal SMILES — dispatch must remain inert
# ---------------------------------------------------------------------------


def test_non_metal_smiles_skips_fffree_dispatch(monkeypatch, clean_env):
    """A purely organic SMILES (no metal atom, no contains_metal hit) must
    skip the fffree dispatch entirely, even with NO_FALLBACK on.  This
    guarantees PT3 mode does not break organic conversion.
    """
    calls = []

    def _spy(smiles, max_isomers=50):
        calls.append(smiles)
        return [("X 0 0 0", "lab")]

    monkeypatch.setattr(
        "delfin.fffree.converter_backend._fffree_isomers", _spy,
    )
    monkeypatch.setattr(sc, "RDKIT_AVAILABLE", False)   # short-circuit
    monkeypatch.setattr(sc, "contains_metal", lambda _smi: False)
    monkeypatch.setenv("DELFIN_FFFREE_NO_FALLBACK", "1")
    monkeypatch.setenv("DELFIN_FFFREE_PURE_TRACK3", "1")

    res, _ = sc.smiles_to_xyz_isomers("CCO")    # ethanol, no metal
    assert calls == [], (
        "non-metal SMILES must not trigger fffree dispatch under any flag"
    )
