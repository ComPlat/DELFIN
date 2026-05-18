"""Tests for the Welle-5m-Z Treatment-Matrix Dispatcher.

Covers:
- Master gate OFF returns empty dict (bit-exact OFF)
- All 6 patterns (A..F) activate the correct env-flags via the stub classifier
- The :func:`treatment_scope` context manager applies and restores os.environ
- Real-world SMILES (YIRQIC, YIVROM) trigger the expected patterns when the
  dedicated system classifier (Agent Y) is available
"""

from __future__ import annotations

import os

import pytest

from delfin import _treatment_matrix as tm


# ---------------------------------------------------------------------------
# Synthetic classifier outputs — exercise predicates without RDKit dependency
# ---------------------------------------------------------------------------
_BASE = {
    "class": "sigma",
    "CN": 4,
    "donor_heterogeneity": "homo",
    "chelate_pattern": "all-mono",
    "metal_block": "3d",
    "q_magnitude": 0,
    "hapticity_max": 0,
    "bulky_flag": False,
    "cell_key": "sigma_CN4_homo_all-mono_3d_q0_eta0_nobulky",
}


def _features(**overrides):
    """Helper: synthetic feature dict with overrides."""
    out = dict(_BASE)
    out.update(overrides)
    return out


# ---------------------------------------------------------------------------
# Master gate
# ---------------------------------------------------------------------------
def test_master_gate_off_returns_empty(monkeypatch):
    """When the master gate is OFF the dispatcher must return ``{}`` for
    every input — bit-exact OFF behaviour."""
    monkeypatch.delenv("DELFIN_5M_TREATMENT_MATRIX", raising=False)

    class _DummyMol:
        pass

    assert tm.dispatch_treatment(_DummyMol()) == {}
    assert tm.dispatch_treatment(None) == {}


def test_master_gate_on_with_none_mol_returns_empty(monkeypatch):
    """Defensive: ``None`` molecule still returns empty even when ON."""
    monkeypatch.setenv("DELFIN_5M_TREATMENT_MATRIX", "1")
    assert tm.dispatch_treatment(None) == {}


# ---------------------------------------------------------------------------
# Per-pattern unit tests via the predicates directly
# ---------------------------------------------------------------------------
def test_pattern_a_sigma_cn6_maximal_asym_with_chelate():
    """Pattern A — sigma+CN6+maximal-asym+(any chelate)."""
    f = _features(
        CN=6,
        donor_heterogeneity="maximal-asym",
        chelate_pattern="mono-plus-bid",
    )
    assert tm._pattern_a(f) is True
    # Negative: missing maximal-asym
    assert tm._pattern_a(_features(CN=6, chelate_pattern="mono-plus-bid")) is False
    # Negative: wrong CN
    assert (
        tm._pattern_a(
            _features(
                CN=5,
                donor_heterogeneity="maximal-asym",
                chelate_pattern="mono-plus-bid",
            )
        )
        is False
    )
    # Negative: all-mono (no chelate)
    assert (
        tm._pattern_a(
            _features(
                CN=6,
                donor_heterogeneity="maximal-asym",
                chelate_pattern="all-mono",
            )
        )
        is False
    )


def test_pattern_b_tris_bidentate():
    """Pattern B — tris-bidentate."""
    assert tm._pattern_b(_features(chelate_pattern="tris-bid")) is True
    assert tm._pattern_b(_features(chelate_pattern="bis-bid")) is False
    assert tm._pattern_b(_features(chelate_pattern="all-mono")) is False


def test_pattern_c_extreme_charge():
    """Pattern C — |q| >= 4."""
    assert tm._pattern_c(_features(q_magnitude=4)) is True
    assert tm._pattern_c(_features(q_magnitude=5)) is True
    assert tm._pattern_c(_features(q_magnitude=3)) is False
    assert tm._pattern_c(_features(q_magnitude=0)) is False


def test_pattern_d_hapto_eta5():
    """Pattern D — hapto + eta-5 piano-stool."""
    assert (
        tm._pattern_d(_features(**{"class": "hapto"}, hapticity_max=5)) is True
    )
    assert (
        tm._pattern_d(_features(**{"class": "hapto"}, hapticity_max=6)) is False
    )
    assert (
        tm._pattern_d(_features(**{"class": "sigma"}, hapticity_max=5)) is False
    )


def test_pattern_e_bulky_alkyl():
    """Pattern E — bulky flag set."""
    assert tm._pattern_e(_features(bulky_flag=True)) is True
    assert tm._pattern_e(_features(bulky_flag=False)) is False


def test_pattern_f_multi_hapto():
    """Pattern F — multi-hapto cluster."""
    assert tm._pattern_f(_features(**{"class": "multi-hapto"})) is True
    assert tm._pattern_f(_features(**{"class": "hapto"})) is False
    assert tm._pattern_f(_features(**{"class": "sigma"})) is False


# ---------------------------------------------------------------------------
# End-to-end dispatch via monkeypatched classifier
# ---------------------------------------------------------------------------
def test_dispatch_pattern_a_activates_burnside_and_theorem_d(monkeypatch):
    monkeypatch.setenv("DELFIN_5M_TREATMENT_MATRIX", "1")
    monkeypatch.setattr(
        tm,
        "_classify",
        lambda mol: _features(
            CN=6,
            donor_heterogeneity="maximal-asym",
            chelate_pattern="mono-plus-bid",
        ),
    )
    out = tm.dispatch_treatment(object())
    assert out["DELFIN_BURNSIDE_FULL"] == "1"
    assert out["DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE"] == "1"


def test_dispatch_pattern_b_activates_theorem_d(monkeypatch):
    monkeypatch.setenv("DELFIN_5M_TREATMENT_MATRIX", "1")
    monkeypatch.setattr(
        tm,
        "_classify",
        lambda mol: _features(chelate_pattern="tris-bid"),
    )
    out = tm.dispatch_treatment(object())
    assert out == {"DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE": "1"}


def test_dispatch_pattern_c_activates_extreme_charge(monkeypatch):
    monkeypatch.setenv("DELFIN_5M_TREATMENT_MATRIX", "1")
    monkeypatch.setattr(tm, "_classify", lambda mol: _features(q_magnitude=4))
    out = tm.dispatch_treatment(object())
    assert out == {"DELFIN_EXTREME_CHARGE_FALLBACK": "1"}


def test_dispatch_pattern_d_activates_cp_piano_stool(monkeypatch):
    monkeypatch.setenv("DELFIN_5M_TREATMENT_MATRIX", "1")
    monkeypatch.setattr(
        tm,
        "_classify",
        lambda mol: _features(**{"class": "hapto"}, hapticity_max=5),
    )
    out = tm.dispatch_treatment(object())
    assert out == {"DELFIN_5J_A_CP_PIANO_STOOL": "1"}


def test_dispatch_pattern_e_activates_rotamer_diversity(monkeypatch):
    monkeypatch.setenv("DELFIN_5M_TREATMENT_MATRIX", "1")
    monkeypatch.setattr(tm, "_classify", lambda mol: _features(bulky_flag=True))
    out = tm.dispatch_treatment(object())
    assert out["DELFIN_5L_T6_ROTAMER_DIVERSITY"] == "1"
    assert out["DELFIN_5L_T6_ROTAMER_K"] == "3"


def test_dispatch_pattern_f_activates_mm_rigid_drag(monkeypatch):
    monkeypatch.setenv("DELFIN_5M_TREATMENT_MATRIX", "1")
    monkeypatch.setattr(
        tm,
        "_classify",
        lambda mol: _features(**{"class": "multi-hapto"}),
    )
    out = tm.dispatch_treatment(object())
    assert out == {"DELFIN_5J_G_MM_RIGID_DRAG": "1"}


def test_dispatch_multi_pattern_merge(monkeypatch):
    """Multiple matching patterns merge — YIRQIC-style sigma+CN6+maximal-asym
    + extreme charge fires both Pattern A and Pattern C."""
    monkeypatch.setenv("DELFIN_5M_TREATMENT_MATRIX", "1")
    monkeypatch.setattr(
        tm,
        "_classify",
        lambda mol: _features(
            CN=6,
            donor_heterogeneity="maximal-asym",
            chelate_pattern="mono-plus-bid",
            q_magnitude=5,
        ),
    )
    out = tm.dispatch_treatment(object())
    # Pattern A flags
    assert out["DELFIN_BURNSIDE_FULL"] == "1"
    assert out["DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE"] == "1"
    # Pattern C flag
    assert out["DELFIN_EXTREME_CHARGE_FALLBACK"] == "1"


# ---------------------------------------------------------------------------
# treatment_scope context manager
# ---------------------------------------------------------------------------
def test_treatment_scope_no_op_when_master_gate_off(monkeypatch):
    """Scope must NOT touch os.environ when master gate is OFF."""
    monkeypatch.delenv("DELFIN_5M_TREATMENT_MATRIX", raising=False)
    monkeypatch.delenv("DELFIN_BURNSIDE_FULL", raising=False)

    before = dict(os.environ)
    with tm.treatment_scope(object()) as flags:
        assert flags == {}
        # No flags appear in env
        assert "DELFIN_BURNSIDE_FULL" not in os.environ
    after = dict(os.environ)
    assert before == after


def test_treatment_scope_applies_and_restores(monkeypatch):
    """Scope sets the flags on entry and restores prior env on exit."""
    monkeypatch.setenv("DELFIN_5M_TREATMENT_MATRIX", "1")
    monkeypatch.delenv("DELFIN_BURNSIDE_FULL", raising=False)
    monkeypatch.delenv("DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE", raising=False)

    monkeypatch.setattr(
        tm,
        "_classify",
        lambda mol: _features(
            CN=6,
            donor_heterogeneity="maximal-asym",
            chelate_pattern="mono-plus-bid",
        ),
    )

    assert "DELFIN_BURNSIDE_FULL" not in os.environ
    with tm.treatment_scope(object()) as flags:
        assert os.environ["DELFIN_BURNSIDE_FULL"] == "1"
        assert os.environ["DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE"] == "1"
        assert flags["DELFIN_BURNSIDE_FULL"] == "1"
    # Restored
    assert "DELFIN_BURNSIDE_FULL" not in os.environ
    assert "DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE" not in os.environ


def test_treatment_scope_restores_pre_existing_value(monkeypatch):
    """If an overridden env-flag already had a value, it must be restored
    (not deleted) on scope exit."""
    monkeypatch.setenv("DELFIN_5M_TREATMENT_MATRIX", "1")
    monkeypatch.setenv("DELFIN_BURNSIDE_FULL", "preexisting")
    monkeypatch.setattr(
        tm,
        "_classify",
        lambda mol: _features(
            CN=6,
            donor_heterogeneity="maximal-asym",
            chelate_pattern="mono-plus-bid",
        ),
    )
    with tm.treatment_scope(object()):
        assert os.environ["DELFIN_BURNSIDE_FULL"] == "1"
    assert os.environ["DELFIN_BURNSIDE_FULL"] == "preexisting"


# ---------------------------------------------------------------------------
# Real-SMILES validation (only runs when RDKit + Agent Y's classifier are
# importable; otherwise skipped).
# ---------------------------------------------------------------------------
_YIRQIC_SMILES = (
    "[O+]#[C][Re-5]1([Br])([C]#[O+])([P+](C2=CC=CC=C2)(C2=CC=CC=C2)"
    "C2=CC=CC=C2)[C+]2N(C3=CC=CC=C3)C=CN2C2=CC=CC=[N+]21"
)
_YIVROM_SMILES = (
    "CCN(CC)C1=NC(C2=CC=CC=C2)=[O+][Fe-3]23([O+]=C(C4=CC=CC=C4)"
    "N=C(N(CC)CC)[S]2)([O+]=C(C2=CC=CC=C2)N=C(N(CC)CC)[S]3)[S]1"
)


def _parse_or_skip(smiles):
    try:
        from rdkit import Chem
    except Exception:
        pytest.skip("RDKit not available")
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        pytest.skip(f"RDKit failed to parse SMILES: {smiles[:40]}")
    return mol


def test_real_yirqic_matches_pattern_a_and_c(monkeypatch):
    """X10-YIRQIC: sigma + CN6 + maximal-asym + chelate + |q|=5
    expected → Pattern A (Burnside+TheoremD) AND Pattern C (extreme charge)."""
    mol = _parse_or_skip(_YIRQIC_SMILES)
    if tm.classifier_source() != "system_classifier":
        pytest.skip(
            "Real-SMILES test requires delfin._system_classifier (Agent Y)"
        )
    monkeypatch.setenv("DELFIN_5M_TREATMENT_MATRIX", "1")
    out = tm.dispatch_treatment(mol)
    assert out.get("DELFIN_BURNSIDE_FULL") == "1"
    assert out.get("DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE") == "1"
    assert out.get("DELFIN_EXTREME_CHARGE_FALLBACK") == "1"


def test_real_yivrom_matches_pattern_b(monkeypatch):
    """X10-YIVROM: tris-bid chelate-pattern → Pattern B (Theorem-D)."""
    mol = _parse_or_skip(_YIVROM_SMILES)
    if tm.classifier_source() != "system_classifier":
        pytest.skip(
            "Real-SMILES test requires delfin._system_classifier (Agent Y)"
        )
    monkeypatch.setenv("DELFIN_5M_TREATMENT_MATRIX", "1")
    out = tm.dispatch_treatment(mol)
    assert out.get("DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE") == "1"
    # YIVROM is NOT charge-extreme (|q|=3)
    assert "DELFIN_EXTREME_CHARGE_FALLBACK" not in out


# ---------------------------------------------------------------------------
# Introspection
# ---------------------------------------------------------------------------
def test_matrix_summary_has_six_patterns():
    summary = tm.matrix_summary()
    assert len(summary) == 6
    names = [entry["predicate_name"] for entry in summary]
    assert names == [
        "_pattern_a",
        "_pattern_b",
        "_pattern_c",
        "_pattern_d",
        "_pattern_e",
        "_pattern_f",
    ]
    # All flags are strings — round-trippable through os.environ
    for entry in summary:
        for k, v in entry["flags"].items():
            assert isinstance(k, str) and isinstance(v, str)


def test_classifier_source_is_string():
    assert tm.classifier_source() in {"system_classifier", "stub"}
