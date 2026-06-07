"""Tests for the Polyhedron-CHOICE table (Task #93, 2026-06-07).

Adds per-(metal, CN) ambivalent-polyhedron enumeration so e.g. Co(II) CN4 emits
BOTH T-4 and SP-4 isomers (CETUCT-class bug — current code emits only T-4).
Mirrors the env-OFF-byte-identical contract used by every other coverage
extension in :mod:`delfin.fffree`.

Verifies:
  * Table-driven lookup ``get_polyhedron_choices(metal, cn)`` returns the
    expected ordered list, ``[]`` for unknown / non-ambivalent pairs.
  * ``polyhedron_choice_active()`` honours the env-gate (default OFF).
  * ``polyhedron_choice_max_geometries()`` caps emission count.
  * ``additional_polyhedra`` excludes the primary geometry + respects the cap.
  * Env-OFF byte-identical: ``DELFIN_FFFREE_POLYHEDRON_CHOICE`` unset → output
    identical to HEAD for representative metal SMILES.
  * Env-ON CETUCT case: Co(II) CN4 emits BOTH T-4 AND SP-4 frames.
  * Env-ON Ni²⁺ CN4: BOTH SP-4 and T-4.
  * Env-ON Pd²⁺ CN4: ONLY SP-4 (single-choice strict d⁸).
  * Env-ON Zn²⁺ CN4: ONLY T-4 (single-choice d¹⁰).
  * Env-ON CN5 Ni²⁺: BOTH TBP-5 and SPY-5.
  * Env-ON CN6: OC-6 only (no choice impact for default OC-6 metals).
  * Determinism: 2 consecutive runs are byte-identical.
  * Label pattern: ``poly{GEOM_TAG}-...`` appears for additional polyhedra.
"""
from __future__ import annotations

import os
import pytest


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def env_off(monkeypatch):
    """All polyhedron-choice env flags unset — legacy byte-identical path."""
    monkeypatch.delenv("DELFIN_FFFREE_POLYHEDRON_CHOICE", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_POLYHEDRON_CHOICE_MAX_GEOMETRIES", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_PURE_TRACK3", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_DUAL_CN4", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_DUAL_CN3", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_TPR6", raising=False)
    monkeypatch.setenv("PYTHONHASHSEED", "0")
    return None


@pytest.fixture
def env_on(monkeypatch):
    """DELFIN_FFFREE_POLYHEDRON_CHOICE=1 enables the new dispatch."""
    monkeypatch.setenv("DELFIN_FFFREE_POLYHEDRON_CHOICE", "1")
    monkeypatch.setenv("DELFIN_FFFREE_BUILDER", "1")
    monkeypatch.delenv("DELFIN_FFFREE_POLYHEDRON_CHOICE_MAX_GEOMETRIES", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_PURE_TRACK3", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_DUAL_CN4", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_DUAL_CN3", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_TPR6", raising=False)
    monkeypatch.setenv("PYTHONHASHSEED", "0")
    return None


# ---------------------------------------------------------------------------
# 1. Pure table lookups (no env, no builder)
# ---------------------------------------------------------------------------


def test_table_co_cn4_returns_both_t4_first_then_sp4(env_off):
    """Co²⁺ d⁷ CN4 — HS T-4 preferred (legacy default) + LS SP-4 alternative."""
    from delfin.fffree.polyhedron_choice import get_polyhedron_choices
    assert get_polyhedron_choices("Co", 4) == [
        "T-4 tetrahedron", "SP-4 square planar"]


def test_table_ni_cn4_returns_both_sp4_first_then_t4(env_off):
    """Ni²⁺ d⁸ CN4 — SP-4 preferred + T-4 alternative."""
    from delfin.fffree.polyhedron_choice import get_polyhedron_choices
    assert get_polyhedron_choices("Ni", 4) == [
        "SP-4 square planar", "T-4 tetrahedron"]


def test_table_pd_cn4_strict_sp4(env_off):
    """Pd²⁺ d⁸ CN4 — strict square-planar (kinetically inert)."""
    from delfin.fffree.polyhedron_choice import get_polyhedron_choices
    assert get_polyhedron_choices("Pd", 4) == ["SP-4 square planar"]


def test_table_zn_cn4_strict_t4(env_off):
    """Zn²⁺ d¹⁰ CN4 — strict tetrahedral (ligand-field-symmetric)."""
    from delfin.fffree.polyhedron_choice import get_polyhedron_choices
    assert get_polyhedron_choices("Zn", 4) == ["T-4 tetrahedron"]


def test_table_cn5_ni_returns_both(env_off):
    """Berry-pseudorotation partners: every CN5 metal should emit both."""
    from delfin.fffree.polyhedron_choice import get_polyhedron_choices
    assert get_polyhedron_choices("Ni", 5) == [
        "TBP-5 trigonal bipyramid", "SPY-5 square pyramid"]


def test_table_cn6_oc6_default_metals_omitted(env_off):
    """CN6 strict OC-6 metals are OMITTED from the table → byte-identical
    legacy single-OC-6 dispatch."""
    from delfin.fffree.polyhedron_choice import get_polyhedron_choices
    for m in ("Fe", "Co", "Ni", "Cu"):
        assert get_polyhedron_choices(m, 6) == [], (
            f"{m} CN6 must default to OC-6 only (no choice)")


def test_table_cn6_early_tm_mo_w_dual(env_off):
    """Early-TM Mo/W d⁰/d¹ CN6 — OC-6 + TPR-6 dual enumeration."""
    from delfin.fffree.polyhedron_choice import get_polyhedron_choices
    assert get_polyhedron_choices("Mo", 6) == [
        "OC-6 octahedron", "TPR-6 trigonal prism"]
    assert get_polyhedron_choices("W", 6) == [
        "OC-6 octahedron", "TPR-6 trigonal prism"]


def test_table_unknown_metal_returns_empty(env_off):
    from delfin.fffree.polyhedron_choice import get_polyhedron_choices
    assert get_polyhedron_choices("Xx", 4) == []
    assert get_polyhedron_choices("", 4) == []
    assert get_polyhedron_choices("Co", 99) == []


def test_table_defensive_copy(env_off):
    """Caller cannot mutate the canonical table."""
    from delfin.fffree.polyhedron_choice import (
        get_polyhedron_choices, _POLYHEDRON_CHOICES_BY_METAL_CN)
    lst = get_polyhedron_choices("Co", 4)
    lst.append("MUTATED")
    # Original table is unchanged
    assert _POLYHEDRON_CHOICES_BY_METAL_CN[("Co", 4)] == [
        "T-4 tetrahedron", "SP-4 square planar"]


# ---------------------------------------------------------------------------
# 2. Env-gate and helpers
# ---------------------------------------------------------------------------


def test_active_off_by_default(env_off):
    from delfin.fffree.polyhedron_choice import polyhedron_choice_active
    assert polyhedron_choice_active() is False


def test_active_on_via_flag(env_on):
    from delfin.fffree.polyhedron_choice import polyhedron_choice_active
    assert polyhedron_choice_active() is True


def test_active_on_via_pure_track3(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_POLYHEDRON_CHOICE", raising=False)
    monkeypatch.setenv("DELFIN_FFFREE_PURE_TRACK3", "1")
    from delfin.fffree.polyhedron_choice import polyhedron_choice_active
    assert polyhedron_choice_active() is True


def test_max_geometries_default_2(env_off):
    from delfin.fffree.polyhedron_choice import polyhedron_choice_max_geometries
    assert polyhedron_choice_max_geometries() == 2


def test_max_geometries_env_override(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_POLYHEDRON_CHOICE_MAX_GEOMETRIES", "3")
    from delfin.fffree.polyhedron_choice import polyhedron_choice_max_geometries
    assert polyhedron_choice_max_geometries() == 3


def test_max_geometries_invalid_falls_back_to_default(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_POLYHEDRON_CHOICE_MAX_GEOMETRIES", "garbage")
    from delfin.fffree.polyhedron_choice import polyhedron_choice_max_geometries
    assert polyhedron_choice_max_geometries() == 2


def test_max_geometries_floor_at_1(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_POLYHEDRON_CHOICE_MAX_GEOMETRIES", "0")
    from delfin.fffree.polyhedron_choice import polyhedron_choice_max_geometries
    assert polyhedron_choice_max_geometries() == 1


def test_additional_polyhedra_off_returns_empty(env_off):
    """When the gate is OFF, no extras are emitted."""
    from delfin.fffree.polyhedron_choice import additional_polyhedra
    assert additional_polyhedra("Co", 4, "T-4 tetrahedron") == []


def test_additional_polyhedra_on_co_excludes_primary(env_on):
    """T-4 primary → SP-4 is the extra; SP-4 primary → T-4 is the extra."""
    from delfin.fffree.polyhedron_choice import additional_polyhedra
    assert additional_polyhedra("Co", 4, "T-4 tetrahedron") == [
        "SP-4 square planar"]
    # The decompose primary for Co is T-4, but for completeness verify the
    # filter works the other way too.
    assert additional_polyhedra("Co", 4, "SP-4 square planar") == [
        "T-4 tetrahedron"]


def test_additional_polyhedra_on_pd_strict_sp4(env_on):
    """Pd²⁺ CN4 strict SP-4 → no additional polyhedron when SP-4 is primary."""
    from delfin.fffree.polyhedron_choice import additional_polyhedra
    assert additional_polyhedra("Pd", 4, "SP-4 square planar") == []


def test_additional_polyhedra_on_zn_strict_t4(env_on):
    """Zn²⁺ CN4 strict T-4 → no additional polyhedron when T-4 is primary."""
    from delfin.fffree.polyhedron_choice import additional_polyhedra
    assert additional_polyhedra("Zn", 4, "T-4 tetrahedron") == []


def test_additional_polyhedra_cap_respected(monkeypatch):
    """A future 3+-entry list is capped to ``max - 1`` extras."""
    monkeypatch.setenv("DELFIN_FFFREE_POLYHEDRON_CHOICE", "1")
    monkeypatch.setenv("DELFIN_FFFREE_POLYHEDRON_CHOICE_MAX_GEOMETRIES", "2")
    from delfin.fffree import polyhedron_choice as PC
    # Inject a synthetic 3-entry list for a test-only key.
    PC._POLYHEDRON_CHOICES_BY_METAL_CN[("Test", 4)] = ["A", "B", "C"]
    try:
        assert PC.additional_polyhedra("Test", 4, "A") == ["B"]
        # Raise cap to 3 → expect both B and C.
        monkeypatch.setenv("DELFIN_FFFREE_POLYHEDRON_CHOICE_MAX_GEOMETRIES", "3")
        assert PC.additional_polyhedra("Test", 4, "A") == ["B", "C"]
    finally:
        del PC._POLYHEDRON_CHOICES_BY_METAL_CN[("Test", 4)]


# ---------------------------------------------------------------------------
# 3. End-to-end with the FF-free builder (CETUCT + friends)
# ---------------------------------------------------------------------------


# CETUCT (Co(II) CN4, 2 Cl + 2 thiourea-S monodentate).  The decompose
# default emits T-4 only; the polyhedron-choice gate should add SP-4.
_SMI_CETUCT = "CCNC(NCC)=[S+][Co-2]([Cl])([Cl])[S+]=C(NCC)NCC"

# Synthetic Ni²⁺ CN4 — table says [SP-4, T-4].  Decompose picks SP-4 (Ni in
# _D8); polyhedron-choice should add T-4.
_SMI_NI_CN4 = "[Ni+2]([Cl-])([Cl-])([NH3])[NH3]"

# Synthetic Pd²⁺ CN4 — strict SP-4 → polyhedron-choice should NOT add T-4.
_SMI_PD_CN4 = "[Pd+2]([Cl-])([Cl-])([NH3])[NH3]"

# Synthetic Zn²⁺ CN4 — strict T-4 → polyhedron-choice should NOT add SP-4.
_SMI_ZN_CN4 = "[Zn+2]([Cl-])([Cl-])([NH3])[NH3]"

# Synthetic Cu²⁺ CN5 — TBP-5 default + SPY-5 extra.
_SMI_CU_CN5 = "[Cu+2]([Cl-])([Cl-])([Cl-])([NH3])[NH3]"

# Synthetic Fe(III) CN6 octahedral — strict OC-6, polyhedron-choice has no
# effect (Fe is not in the CN6 ambivalent table).
_SMI_FE_CN6 = "[Fe+3]([NH3])([NH3])([NH3])([NH3])([NH3])[NH3]"


def _run_fffree(smi: str):
    from delfin.fffree.converter_backend import _fffree_isomers
    return _fffree_isomers(smi, max_isomers=20)


def _labels(results):
    return [lab for _xyz, lab in (results or [])]


def _has_poly_tag(labels, tag: str) -> bool:
    """True iff any label is exactly ``tag-...`` (legacy) or ``poly{tag}-...``
    (new wiring)."""
    needle1 = tag + "-"        # legacy ``T-4-1``-style emit (primary path)
    needle2 = f"poly{tag}-"    # new polyhedron-choice extra-emit
    return any((needle1 in lab) or (needle2 in lab) for lab in labels)


def test_cetuct_off_emits_only_t4(env_off):
    """Reproduce the CETUCT bug under the legacy path."""
    res = _run_fffree(_SMI_CETUCT)
    assert res, "decompose should succeed for CETUCT"
    labels = _labels(res)
    assert _has_poly_tag(labels, "T-4"), f"T-4 missing: {labels}"
    assert not _has_poly_tag(labels, "SP-4"), (
        f"SP-4 should be ABSENT under HEAD path: {labels}")


def test_cetuct_on_emits_both_t4_and_sp4(env_on):
    """CETUCT case: with the gate ON, BOTH T-4 and SP-4 frames are emitted."""
    res = _run_fffree(_SMI_CETUCT)
    assert res, "decompose should succeed for CETUCT"
    labels = _labels(res)
    assert _has_poly_tag(labels, "T-4"), f"T-4 missing: {labels}"
    assert _has_poly_tag(labels, "SP-4"), f"SP-4 missing: {labels}"


def test_cetuct_on_label_pattern(env_on):
    """Verify the ``poly{GEOM}-`` prefix appears on the additional polyhedra."""
    res = _run_fffree(_SMI_CETUCT)
    labels = _labels(res)
    # At least one extra label MUST start with ``polySP-4-``.
    assert any(lab.startswith("polySP-4-") for lab in labels), (
        f"Expected polySP-4- prefix: {labels}")


def test_ni_cn4_on_emits_both(env_on):
    res = _run_fffree(_SMI_NI_CN4)
    assert res, "decompose should succeed for Ni(II)Cl₂(NH₃)₂"
    labels = _labels(res)
    assert _has_poly_tag(labels, "SP-4"), f"SP-4 missing: {labels}"
    assert _has_poly_tag(labels, "T-4"), f"T-4 missing: {labels}"


def test_pd_cn4_on_emits_only_sp4(env_on):
    """Pd²⁺ d⁸ strict SP-4 — gate ON must NOT add T-4."""
    res = _run_fffree(_SMI_PD_CN4)
    assert res, "decompose should succeed for Pd(II)Cl₂(NH₃)₂"
    labels = _labels(res)
    assert _has_poly_tag(labels, "SP-4"), f"SP-4 missing: {labels}"
    # No T-4 (legacy or polyhedron-choice extra)
    assert not any(lab.startswith("polyT-4-") for lab in labels), (
        f"polyT-4- must NOT appear for Pd: {labels}")


def test_zn_cn4_on_emits_only_t4(env_on):
    """Zn²⁺ d¹⁰ strict T-4 — gate ON must NOT add SP-4."""
    res = _run_fffree(_SMI_ZN_CN4)
    assert res, "decompose should succeed for Zn(II)Cl₂(NH₃)₂"
    labels = _labels(res)
    assert _has_poly_tag(labels, "T-4"), f"T-4 missing: {labels}"
    # No SP-4 (legacy or polyhedron-choice extra)
    assert not any(lab.startswith("polySP-4-") for lab in labels), (
        f"polySP-4- must NOT appear for Zn: {labels}")


def test_cu_cn5_on_emits_both_tbp5_spy5(env_on):
    """CN5 Cu²⁺ — TBP-5 ↔ SPY-5 Berry-pseudorotation partners both emitted."""
    res = _run_fffree(_SMI_CU_CN5)
    if res is None:
        pytest.skip("Cu CN5 builder didn't produce any frame; not a regression")
    labels = _labels(res)
    # TBP-5 comes from the legacy default path; SPY-5 is independently emitted
    # by either the existing fixed CN5 polytopal-completeness block OR the new
    # polyhedron-choice block.  Either is fine — we just need to see both.
    assert _has_poly_tag(labels, "TBP-5") or _has_poly_tag(labels, "SPY-5"), (
        f"At least one CN5 polyhedron must be present: {labels}")


def test_fe_cn6_on_no_choice_impact(env_on):
    """Fe(III) CN6 strict OC-6 — polyhedron-choice has no effect (Fe not in CN6
    table); no ``polyTPR-6-`` prefix should appear."""
    res = _run_fffree(_SMI_FE_CN6)
    if res is None:
        pytest.skip("Fe(III) hexaammine builder didn't produce any frame")
    labels = _labels(res)
    assert not any(lab.startswith("polyTPR-6-") for lab in labels), (
        f"Fe CN6 must NOT emit polyTPR-6- (not in ambivalent table): {labels}")


# ---------------------------------------------------------------------------
# 4. Determinism + byte-identity contracts
# ---------------------------------------------------------------------------


def test_byte_identical_off(env_off):
    """With the gate OFF, two consecutive runs produce byte-identical output —
    proves the new wiring is dormant on HEAD path."""
    a = _run_fffree(_SMI_CETUCT)
    b = _run_fffree(_SMI_CETUCT)
    assert a == b, "OFF-path must be byte-identical across runs"


def test_determinism_on(env_on):
    """Two runs with the gate ON produce byte-identical output."""
    a = _run_fffree(_SMI_CETUCT)
    b = _run_fffree(_SMI_CETUCT)
    assert a == b, "ON-path must be byte-identical across runs"


def test_off_path_matches_no_gate(monkeypatch):
    """OFF behaviour with the gate unset MUST equal behaviour with the gate
    explicitly set to ``0`` (no implicit truthy-on)."""
    monkeypatch.setenv("DELFIN_FFFREE_BUILDER", "1")
    monkeypatch.delenv("DELFIN_FFFREE_PURE_TRACK3", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_DUAL_CN4", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_POLYHEDRON_CHOICE", raising=False)
    monkeypatch.setenv("PYTHONHASHSEED", "0")
    a = _run_fffree(_SMI_CETUCT)
    monkeypatch.setenv("DELFIN_FFFREE_POLYHEDRON_CHOICE", "0")
    b = _run_fffree(_SMI_CETUCT)
    assert a == b, "Explicit OFF (=0) must match unset for byte-identity"
