"""Tests for the CN=5 polyhedron classifier (DELFIN_CN5_GEOM_AWARE).

Covers both APIs that the scope-extension adds to ``delfin.smiles_converter``:

  * ``_classify_cn5_geometry_from_labels`` — label-only facade (no RDKit
    required), consumed when the rich-mol kwargs are absent.
  * ``_classify_cn5_geometry``             — graph-based variant taking the
    RDKit mol + metal_idx + donor_idxs.

The acceptance contract for the scope-extension is six positive assertions
(see ``ITER-cn5_geom_audit.md``):

  1. d8 Ni(II) CN=5             → SP
  2. d0 Ti(IV) CN=5             → TBP
  3. d10 Cu(I)  CN=5            → TBP
  4. ≥ 2 π-acceptor donors      → SP override (Fe(CO)5 path)
  5. Default-OFF env-flag       → bit-exact legacy emission
  6. Malformed input            → safe-fallback ``'TBP'``

A seventh test ensures tridentate-chelate detection forces SP for d-counts
that would otherwise prefer TBP.

The file skips cleanly when RDKit is missing so it is harmless in lean CI
environments.
"""
from __future__ import annotations

import os
from typing import List

import pytest

pytest.importorskip("rdkit", reason="RDKit required for CN=5 classifier tests")
from rdkit import Chem  # noqa: E402 — after importorskip guard

from delfin.smiles_converter import (  # noqa: E402
    _METAL_GROUP_NUMBER,
    _PREFERRED_CN5_GEOMETRY,
    _classify_cn5_geometry,
    _classify_cn5_geometry_from_labels,
    _cn5_count_pi_acceptor_donors,
    _cn5_d_electron_count,
    _cn5_has_tridentate_chelate,
    _enumerate_topological_isomers,
)


# ----------------------------------------------------------------------
# d-electron count
# ----------------------------------------------------------------------


def test_d_count_ni_ii_is_d8() -> None:
    assert _cn5_d_electron_count("Ni", 2) == 8


def test_d_count_pd_ii_is_d8() -> None:
    assert _cn5_d_electron_count("Pd", 2) == 8


def test_d_count_pt_ii_is_d8() -> None:
    assert _cn5_d_electron_count("Pt", 2) == 8


def test_d_count_co_i_is_d8() -> None:
    assert _cn5_d_electron_count("Co", 1) == 8


def test_d_count_ti_iv_is_d0() -> None:
    assert _cn5_d_electron_count("Ti", 4) == 0


def test_d_count_zn_ii_is_d10() -> None:
    assert _cn5_d_electron_count("Zn", 2) == 10


def test_d_count_cu_i_is_d10() -> None:
    assert _cn5_d_electron_count("Cu", 1) == 10


def test_d_count_clamps_negative_to_zero() -> None:
    # Ti(VI) impossible -> clamp to 0
    assert _cn5_d_electron_count("Ti", 6) == 0


def test_d_count_clamps_above_ten() -> None:
    # Zn(-4) absurd -> clamp to 10
    assert _cn5_d_electron_count("Zn", -4) == 10


def test_d_count_main_group_returns_none() -> None:
    assert _cn5_d_electron_count("Ca", 2) is None
    assert _cn5_d_electron_count("Al", 3) is None
    assert _cn5_d_electron_count("K", 1) is None


def test_d_count_unknown_element_returns_none() -> None:
    assert _cn5_d_electron_count("XX", 0) is None


# ----------------------------------------------------------------------
# π-acceptor donor counter
# ----------------------------------------------------------------------


def test_pi_acceptor_count_phosphine() -> None:
    assert _cn5_count_pi_acceptor_donors(["P0", "P1", "N0", "N1", "O0"]) == 2


def test_pi_acceptor_count_zero_for_n_o() -> None:
    assert _cn5_count_pi_acceptor_donors(["N0", "N1", "O0", "O1", "Cl0"]) == 0


def test_pi_acceptor_count_arsine_and_stibine() -> None:
    assert _cn5_count_pi_acceptor_donors(["As0", "Sb0", "P0"]) == 3


def test_pi_acceptor_count_strips_class_digits() -> None:
    # "P12" must still be detected
    assert _cn5_count_pi_acceptor_donors(["P12", "N7", "P3"]) == 2


def test_pi_acceptor_count_empty() -> None:
    assert _cn5_count_pi_acceptor_donors([]) == 0


# ----------------------------------------------------------------------
# Tridentate chelate detector
# ----------------------------------------------------------------------


def test_tridentate_chelate_chain() -> None:
    # donors 0-1-2 chain → single connected component of 3
    assert _cn5_has_tridentate_chelate(
        [frozenset([0, 1]), frozenset([1, 2])], n_coord=5
    ) is True


def test_tridentate_chelate_two_separate_bidentate() -> None:
    # 0-1 and 2-3 separate → biggest component = 2 → False
    assert _cn5_has_tridentate_chelate(
        [frozenset([0, 1]), frozenset([2, 3])], n_coord=5
    ) is False


def test_tridentate_chelate_empty() -> None:
    assert _cn5_has_tridentate_chelate([], n_coord=5) is False


def test_tridentate_chelate_wrong_cn() -> None:
    # only fires for n_coord=5
    assert _cn5_has_tridentate_chelate(
        [frozenset([0, 1]), frozenset([1, 2])], n_coord=6
    ) is False


# ----------------------------------------------------------------------
# Label-only facade — required positive cases
# ----------------------------------------------------------------------


def test_facade_d8_ni_ii_is_sp() -> None:
    """d8 Ni(II) CN=5 → SP (acceptance #1)."""
    out = _classify_cn5_geometry_from_labels(
        metal_symbol="Ni",
        formal_charge=2,
        donor_labels=["N0", "N1", "N2", "N3", "N4"],
        chelate_pairs=[],
    )
    assert out == "SP"


def test_facade_d0_ti_iv_is_tbp() -> None:
    """d0 Ti(IV) CN=5 → TBP (acceptance #2)."""
    out = _classify_cn5_geometry_from_labels(
        metal_symbol="Ti",
        formal_charge=4,
        donor_labels=["Cl0", "Cl1", "Cl2", "Cl3", "Cl4"],
        chelate_pairs=[],
    )
    assert out == "TBP"


def test_facade_d10_cu_i_is_tbp() -> None:
    """d10 Cu(I) CN=5 → TBP (acceptance #3)."""
    out = _classify_cn5_geometry_from_labels(
        metal_symbol="Cu",
        formal_charge=1,
        donor_labels=["O0", "O1", "O2", "O3", "O4"],
        chelate_pairs=[],
    )
    assert out == "TBP"


def test_facade_bulky_two_phosphines_override_to_sp() -> None:
    """≥ 2 P/As/Sb donors flip TBP → SP (acceptance #4)."""
    # Fe(II) base preference is TBP (d6 actually -> SP via low-spin rule);
    # pick a d-count that would normally be TBP — d3 Cr(III).
    out = _classify_cn5_geometry_from_labels(
        metal_symbol="Cr",
        formal_charge=3,
        donor_labels=["P0", "P1", "Cl0", "Cl1", "Cl2"],
        chelate_pairs=[],
    )
    assert out == "SP"


def test_facade_tridentate_override_to_sp() -> None:
    """Tridentate (or larger) chelate flips TBP → SP."""
    chel = [frozenset([0, 1]), frozenset([1, 2])]
    out = _classify_cn5_geometry_from_labels(
        metal_symbol="Mn",
        formal_charge=2,
        donor_labels=["N0", "N1", "N2", "Cl0", "Cl1"],
        chelate_pairs=chel,
    )
    assert out == "SP"


def test_facade_safe_fallback_unknown_element() -> None:
    """Unknown metal element → safe TBP fallback (acceptance #6)."""
    out = _classify_cn5_geometry_from_labels(
        metal_symbol="Xx",
        formal_charge=0,
        donor_labels=["N0"] * 5,
        chelate_pairs=[],
    )
    assert out == "TBP"


def test_facade_safe_fallback_empty_args() -> None:
    """Empty donor list + empty metal symbol → TBP."""
    out = _classify_cn5_geometry_from_labels(
        metal_symbol="",
        formal_charge=0,
        donor_labels=[],
        chelate_pairs=None,
    )
    assert out == "TBP"


# ----------------------------------------------------------------------
# Rich graph-based variant — acceptance cases requiring an RDKit mol
# ----------------------------------------------------------------------


def _make_mol(smi: str):
    """Parse a SMILES without sanitising so coordination centres survive."""
    return Chem.MolFromSmiles(smi, sanitize=False)


def test_rich_ni_ii_penta_cyano_is_sp() -> None:
    """[Ni(NC)5]³⁻ — d8 Ni(II) → SP."""
    mol = _make_mol("[Ni+2](#[N-])(#[N-])(#[N-])(#[N-])#[N-]")
    assert mol is not None
    metal = next(a for a in mol.GetAtoms() if a.GetSymbol() == "Ni")
    donors = [n.GetIdx() for n in metal.GetNeighbors()]
    assert _classify_cn5_geometry(mol, metal.GetIdx(), donors) == "SP"


def test_rich_fe_co_5_picks_sp_via_pi_acceptor_override() -> None:
    """Fe(CO)5 — five carbonyl pi-acceptors → SP override."""
    mol = _make_mol("[Fe](=C=O)(=C=O)(=C=O)(=C=O)=C=O")
    assert mol is not None
    metal = next(a for a in mol.GetAtoms() if a.GetSymbol() == "Fe")
    donors = [n.GetIdx() for n in metal.GetNeighbors()]
    assert _classify_cn5_geometry(mol, metal.GetIdx(), donors) == "SP"


def test_rich_cu_i_aqua_is_tbp() -> None:
    """[Cu(OH2)5]+ — d10 Cu(I) → TBP."""
    mol = _make_mol("[Cu+](O)(O)(O)(O)O")
    assert mol is not None
    metal = next(a for a in mol.GetAtoms() if a.GetSymbol() == "Cu")
    donors = [n.GetIdx() for n in metal.GetNeighbors()]
    assert _classify_cn5_geometry(mol, metal.GetIdx(), donors) == "TBP"


def test_rich_pt_ii_penta_phosphine_is_sp() -> None:
    """[Pt(PH3)5]²⁺ — d8 Pt(II) + 5 P donors → SP."""
    mol = _make_mol("[Pt+2](P)(P)(P)(P)P")
    assert mol is not None
    metal = next(a for a in mol.GetAtoms() if a.GetSymbol() == "Pt")
    donors = [n.GetIdx() for n in metal.GetNeighbors()]
    assert _classify_cn5_geometry(mol, metal.GetIdx(), donors) == "SP"


def test_rich_safe_fallback_mol_none() -> None:
    """mol=None → safe TBP fallback."""
    assert _classify_cn5_geometry(None, 0, [0, 1, 2, 3, 4]) == "TBP"


def test_rich_safe_fallback_wrong_donor_count() -> None:
    """len(donor_idxs) != 5 → safe TBP fallback."""
    mol = _make_mol("[Ni+2](N)(N)(N)N")  # CN=4, not 5
    assert mol is not None
    metal = next(a for a in mol.GetAtoms() if a.GetSymbol() == "Ni")
    donors = [n.GetIdx() for n in metal.GetNeighbors()]
    assert _classify_cn5_geometry(mol, metal.GetIdx(), donors) == "TBP"


# ----------------------------------------------------------------------
# Env-flag bit-exact OFF (acceptance #5)
# ----------------------------------------------------------------------


def _strip_env(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("DELFIN_CN5_GEOM_AWARE", raising=False)


def test_default_off_legacy_emission_mn(monkeypatch: pytest.MonkeyPatch) -> None:
    """Mn (legacy preference: TBP) — first emitted geometry must be TBP."""
    _strip_env(monkeypatch)
    out = _enumerate_topological_isomers(
        ["N0", "N1", "Cl0", "Cl1", "O0"], 5, [], metal_symbol="Mn"
    )
    assert out, "enumerator must emit at least one isomer"
    assert out[0][0][0] == "TBP"


def test_default_off_legacy_emission_ni(monkeypatch: pytest.MonkeyPatch) -> None:
    """Ni (legacy preference: SP) — first emitted geometry must be SP."""
    _strip_env(monkeypatch)
    out = _enumerate_topological_isomers(
        ["N0", "N1", "Cl0", "Cl1", "O0"], 5, [], metal_symbol="Ni"
    )
    assert out, "enumerator must emit at least one isomer"
    assert out[0][0][0] == "SP"


def test_default_off_isomer_count_unchanged(monkeypatch: pytest.MonkeyPatch) -> None:
    """Default OFF: the isomer set is the same with or without
    ``metal_formal_charge`` / ``mol`` kwargs (back-compat)."""
    _strip_env(monkeypatch)
    base = _enumerate_topological_isomers(
        ["N0", "N1", "Cl0", "Cl1", "O0"], 5, [], metal_symbol="Mn"
    )
    rich = _enumerate_topological_isomers(
        ["N0", "N1", "Cl0", "Cl1", "O0"], 5, [],
        metal_symbol="Mn", metal_formal_charge=2,
    )
    assert base == rich


def test_env_on_tridentate_flips_first_geom(monkeypatch: pytest.MonkeyPatch) -> None:
    """Env ON + tridentate chelate on a TBP-preferred metal flips
    the FIRST emitted geometry from TBP → SP."""
    _strip_env(monkeypatch)
    out_off = _enumerate_topological_isomers(
        ["N0", "N1", "N2", "Cl0", "Cl1"], 5, [], metal_symbol="Mn"
    )
    monkeypatch.setenv("DELFIN_CN5_GEOM_AWARE", "1")
    chel = [frozenset([0, 1]), frozenset([1, 2])]
    out_on = _enumerate_topological_isomers(
        ["N0", "N1", "N2", "Cl0", "Cl1"], 5, chel,
        metal_symbol="Mn", metal_formal_charge=2,
    )
    assert out_off and out_on
    assert out_off[0][0][0] == "TBP"
    assert out_on[0][0][0] == "SP"


# ----------------------------------------------------------------------
# Group-number table integrity
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "metal,expected_group",
    [
        ("Ni", 10), ("Pd", 10), ("Pt", 10),
        ("Cu", 11), ("Ag", 11), ("Au", 11),
        ("Zn", 12), ("Cd", 12), ("Hg", 12),
        ("Ti", 4), ("Zr", 4), ("Hf", 4),
        ("Mo", 6), ("W", 6), ("Cr", 6),
    ],
)
def test_group_number_table(metal: str, expected_group: int) -> None:
    assert _METAL_GROUP_NUMBER[metal] == expected_group


def test_legacy_preferred_cn5_unchanged() -> None:
    """The legacy CN=5 preference table is not touched by the new helper."""
    assert _PREFERRED_CN5_GEOMETRY["Ni"] == "SP"
    assert _PREFERRED_CN5_GEOMETRY["Fe"] == "TBP"
    assert _PREFERRED_CN5_GEOMETRY["Mn"] == "TBP"
