"""Tests for Welle-5l Track-2 xTB-cascade wire-in (default-OFF infrastructure).

Validates that ``delfin.smiles_converter._apply_xtb_cascade_if_enabled``
is a bit-exact no-op when ``DELFIN_CASCADE_REFINER`` is unset (the
default), and that the env-flag + class-conditional gate plumbing
behaves as documented.

Hard contract: when both env-flags are unset (default) the helper
returns the input ``results`` list unchanged so that callers can wire
it into the production pipeline without changing default behaviour.

All tests skip cleanly if the module or RDKit is missing.
"""
from __future__ import annotations

import os

import pytest

_smiles_converter = pytest.importorskip(
    "delfin.smiles_converter",
    reason="smiles_converter module required",
)
Chem = pytest.importorskip("rdkit.Chem", reason="RDKit required for tests")

_apply = getattr(_smiles_converter, "_apply_xtb_cascade_if_enabled", None)
if _apply is None:
    pytest.skip("xTB cascade wire-in not present", allow_module_level=True)


_SAMPLE_XYZ = (
    "Fe     0.000000      0.000000      0.000000\n"
    "Cl     2.300000      0.000000      0.000000\n"
    "Cl    -2.300000      0.000000      0.000000\n"
    "Cl     0.000000      2.300000      0.000000\n"
    "Cl     0.000000     -2.300000      0.000000\n"
)


@pytest.fixture(autouse=True)
def _clear_env(monkeypatch):
    """Ensure each test starts with all cascade env-flags unset."""
    monkeypatch.delenv("DELFIN_CASCADE_REFINER", raising=False)
    monkeypatch.delenv("DELFIN_CASCADE_REFINER_CLASSES", raising=False)
    monkeypatch.delenv("DELFIN_CASCADE_REFINER_TIMEOUT_S", raising=False)
    monkeypatch.delenv("DELFIN_CASCADE_REFINER_MAX_ITER", raising=False)
    monkeypatch.delenv("DELFIN_CASCADE_REFINER_GFN", raising=False)


def _mol_for_sigma_fe_cl4():
    """Tetrahedral Fe(II)Cl4(2-) — RDKit-parsable, sigma class, charge -2."""
    mol = Chem.MolFromSmiles("[Cl-].[Cl-].[Cl-].[Cl-].[Fe+2]")
    return mol


def test_default_off_no_op_bit_identical():
    """No env-flag set → helper returns the input list unchanged."""
    mol = _mol_for_sigma_fe_cl4()
    assert mol is not None
    results = [(_SAMPLE_XYZ, "frame_a"), (_SAMPLE_XYZ, "frame_b")]
    out = _apply(mol, results, dual_parse_done=False)
    assert out is results or out == results
    # XYZ strings must be byte-identical
    assert out[0][0] == _SAMPLE_XYZ
    assert out[1][0] == _SAMPLE_XYZ


def test_dual_parse_skipped():
    """Inner dual-parse call → helper is a no-op even if flag is set."""
    mol = _mol_for_sigma_fe_cl4()
    os.environ["DELFIN_CASCADE_REFINER"] = "1"
    try:
        results = [(_SAMPLE_XYZ, "frame_a")]
        out = _apply(mol, results, dual_parse_done=True)
        # bit-identical regardless of flag because dual_parse_done=True
        assert out[0][0] == _SAMPLE_XYZ
    finally:
        os.environ.pop("DELFIN_CASCADE_REFINER", None)


def test_empty_results_returned_unchanged():
    """Empty list short-circuit."""
    mol = _mol_for_sigma_fe_cl4()
    assert _apply(mol, [], dual_parse_done=False) == []


def test_mol_none_skipped():
    """``mol=None`` skips the cascade (no class classification possible)."""
    results = [(_SAMPLE_XYZ, "frame_a")]
    out = _apply(None, results, dual_parse_done=False)
    assert out[0][0] == _SAMPLE_XYZ


def test_class_conditional_off_class_no_op():
    """Class allow-list excludes our class → no-op."""
    mol = _mol_for_sigma_fe_cl4()
    os.environ["DELFIN_CASCADE_REFINER_CLASSES"] = "multi_hapto"
    try:
        results = [(_SAMPLE_XYZ, "frame_a")]
        out = _apply(mol, results, dual_parse_done=False)
        # Bit-identical because the molecule's class is not in the allow-list
        assert out[0][0] == _SAMPLE_XYZ
    finally:
        os.environ.pop("DELFIN_CASCADE_REFINER_CLASSES", None)


def test_organic_only_short_circuit():
    """Non-metal molecule → cascade is a no-op even when enabled."""
    mol = Chem.MolFromSmiles("CCO")  # ethanol, no metal
    assert mol is not None
    os.environ["DELFIN_CASCADE_REFINER"] = "1"
    try:
        results = [("C     0.0 0.0 0.0\nC    1.5 0.0 0.0\nO    3.0 0.0 0.0\n",
                    "frame_a")]
        out = _apply(mol, results, dual_parse_done=False)
        # Bit-identical because no metal triggers the metal-only gate
        assert out[0][0] == results[0][0]
    finally:
        os.environ.pop("DELFIN_CASCADE_REFINER", None)


def test_helper_signature_matches_dispatch_pattern():
    """Helper must accept (mol, results, dual_parse_done) per pipeline contract."""
    import inspect
    sig = inspect.signature(_apply)
    params = list(sig.parameters)
    assert params == ["mol", "results", "dual_parse_done"], (
        f"Signature drift: expected (mol, results, dual_parse_done), got {params}"
    )
