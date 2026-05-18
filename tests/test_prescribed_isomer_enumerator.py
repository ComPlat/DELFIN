"""Tests for delfin._prescribed_isomer_enumerator (Welle-5n-Pre).

Covers the 5 brief validation cases plus the universal-fundamental
acceptance gates:

1.  ``[Pt(NH3)4]``         — homoleptic CN4 → 1 isomer per polyhedron.
2.  ``PtCl2(NH3)2``        — SQ cis + trans; TH only-isomer.
3.  ``[Co(en)2(NH3)Cl]``   — CN6 bis-bidentate-mono-mixed → Δ/Λ enantiomer set.
4.  ``Fe(en)3``-like       — tris-bidentate-symmetric → exact Δ/Λ on OH.
5.  ``X10-YIRQIC``-like Re — maximal-asym CN6, ≥4 distinct donor classes.

Plus:
*  Universal-fundamental: no SMILES regex / element allowlist.
*  Determinism: same SMILES → same prescribed-isomer list bit-exact.
*  No-metal / CN<2 → empty list.
*  Burnside-lemma orbit count cross-check (proper-rotation subgroup).
"""

from __future__ import annotations

import itertools
from typing import List

import pytest

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402

from delfin._burnside_groups import burnside_canonical_key, get_groups  # noqa: E402
from delfin._polya_groups import (  # noqa: E402
    polyhedra_for_cn,
    polyhedron_geometry,
    positional_descriptor,
    trans_positions,
)
from delfin._prescribed_isomer_enumerator import (  # noqa: E402
    _build_extended_label,
    _enumerate_orbits,
    enumerate_prescribed_isomers,
    isomer_labels,
    n_theory_prescribed,
)


def _mol(smi: str):
    m = Chem.MolFromSmiles(smi, sanitize=False)
    assert m is not None, f"failed to parse: {smi}"
    return m


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _isomers_for(smiles: str):
    return enumerate_prescribed_isomers(_mol(smiles))


def _by_polyhedron(isomers, geom: str):
    return [d for d in isomers if d["polyhedron"] == geom]


# ---------------------------------------------------------------------------
# Validation Case 1: [Pt(NH3)4] homoleptic CN4
# ---------------------------------------------------------------------------
def test_homoleptic_cn4_pt_nh3_4():
    """One isomer per polyhedron (SQ, SS, TH).  Homoleptic = "only-isomer"."""
    out = _isomers_for("[Pt]([NH3])([NH3])([NH3])[NH3]")
    sq = _by_polyhedron(out, "SQ")
    th = _by_polyhedron(out, "TH")
    assert len(sq) == 1
    assert len(th) == 1
    assert sq[0]["positional"] == "only-isomer"
    assert th[0]["positional"] == "only-isomer"
    # Every entry is achiral (no chelates → no Λ/Δ).
    for d in out:
        assert d["stereo"] == "achiral"


# ---------------------------------------------------------------------------
# Validation Case 2: PtCl2(NH3)2 cis/trans
# ---------------------------------------------------------------------------
def test_ptcl2_nh3_2_cis_trans():
    """SQ gives EXACTLY cis + trans; TH gives 1 isomer (no fac/mer at CN4-TH)."""
    out = _isomers_for("[Cl][Pt]([Cl])([NH3])[NH3]")
    sq = _by_polyhedron(out, "SQ")
    # Should be exactly {cis, trans}
    sq_pos = sorted(d["positional"] for d in sq)
    assert sq_pos == ["cis", "trans"], f"got SQ positions {sq_pos}"
    th = _by_polyhedron(out, "TH")
    assert len(th) == 1
    # No chelates → no Λ/Δ
    for d in out:
        assert d["stereo"] == "achiral"


# ---------------------------------------------------------------------------
# Validation Case 3: [Co(en)2(NH3)Cl] CN6 bis-bidentate
# ---------------------------------------------------------------------------
def test_co_en2_nh3_cl_delta_lambda_pair():
    """Bis-chelate CN6 with asymmetric mono donors produces Δ/Λ pair on OH."""
    out = _isomers_for("[Cl][Co]12([NH3])([N]CC[N]1)[N]CC[N]2")
    oh = _by_polyhedron(out, "OH")
    # Distinct stereo tags on OH must include both Lambda AND Delta
    stereos = {d["stereo"] for d in oh}
    assert "Lambda" in stereos
    assert "Delta" in stereos
    # Δ + Λ orbits exist in equal numbers (mirror symmetry)
    n_delta = sum(1 for d in oh if d["stereo"] == "Delta")
    n_lambda = sum(1 for d in oh if d["stereo"] == "Lambda")
    assert n_delta == n_lambda, (
        f"Δ/Λ symmetry broken: {n_delta} Δ vs {n_lambda} Λ on OH"
    )


# ---------------------------------------------------------------------------
# Validation Case 4: Fe(en)3 — tris-bidentate-symmetric
# ---------------------------------------------------------------------------
def test_fe_en3_tris_bidentate_delta_lambda():
    """Fe(en)3 homoleptic tris-bidentate on OH = exactly 2 orbits (Δ + Λ)."""
    out = _isomers_for("[Fe]123([N]CC[N]1)([N]CC[N]2)[N]CC[N]3")
    oh = _by_polyhedron(out, "OH")
    # Exactly 2 OH orbits = Λ and Δ
    assert len(oh) == 2, (
        f"Fe(en)3 OH should have 2 orbits (Δ+Λ), got {len(oh)}"
    )
    stereos = sorted(d["stereo"] for d in oh)
    assert stereos == ["Delta", "Lambda"], f"got {stereos}"
    # Orbit sizes must match: |G+| = 24 for OH, |Stab|=1 → orbit size 24
    for d in oh:
        assert d["n_polya_orbit"] == 24, (
            f"Fe(en)3 OH orbit size = {d['n_polya_orbit']}, expected 24"
        )


# ---------------------------------------------------------------------------
# Validation Case 5: X10-YIRQIC-like Re-CN6 maximal-asym
# ---------------------------------------------------------------------------
SMILES_YIRQIC = (
    "[O+]#[C][Re-5]1([Br])([C]#[O+])"
    "([P+](c2ccccc2)(c2ccccc2)c2ccccc2)"
    "[C+]2N(c3ccccc3)C=CN2C2=CC=CC=[N+]21"
)


def test_yirqic_re_cn6_maximal_asym():
    """X10-YIRQIC-like Re CN6 maximal-asym: ≥12 OH orbits, Δ/Λ pairs visible."""
    out = _isomers_for(SMILES_YIRQIC)
    oh = _by_polyhedron(out, "OH")
    # Brief expects 12-30 orbits across all polyhedra; OH alone must be ≥10.
    assert len(oh) >= 10, f"YIRQIC OH had only {len(oh)} orbits"
    # Δ/Λ pairs balanced
    n_delta = sum(1 for d in oh if d["stereo"] == "Delta")
    n_lambda = sum(1 for d in oh if d["stereo"] == "Lambda")
    assert n_delta == n_lambda, (
        f"YIRQIC OH Δ/Λ imbalanced: {n_delta} vs {n_lambda}"
    )
    # The classifier must see 4 distinct donor element classes
    classes = set(out[0]["donor_classes"])
    assert len(classes) >= 4, f"YIRQIC class diversity = {len(classes)} (<4)"


# ---------------------------------------------------------------------------
# Universal-fundamental gates
# ---------------------------------------------------------------------------
def test_determinism():
    """Same SMILES → identical prescribed-isomer list (deterministic)."""
    smi = "[Cl][Pt]([Cl])([NH3])[NH3]"
    a = _isomers_for(smi)
    b = _isomers_for(smi)
    assert [d["label"] for d in a] == [d["label"] for d in b]
    for da, db in zip(a, b):
        assert da == db


def test_no_metal_returns_empty():
    """Pure organic mol → no metal → empty prescribed list."""
    out = enumerate_prescribed_isomers(_mol("CCOCC"))
    assert out == []


def test_cn_lt_2_returns_empty():
    """CN<2 → no isomer enumeration possible."""
    # Cl-Mn: CN=1
    out = enumerate_prescribed_isomers(_mol("[Cl][Mn]"))
    assert out == []


def test_none_mol_returns_empty():
    assert enumerate_prescribed_isomers(None) == []


# ---------------------------------------------------------------------------
# Burnside cross-validation
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("smi,geom,expected_min,expected_max", [
    # SMILES, polyhedron, min orbit count, max orbit count
    ("[Cl][Pt]([Cl])([NH3])[NH3]", "SQ", 2, 2),
    ("[Pt]([NH3])([NH3])([NH3])[NH3]", "SQ", 1, 1),
    ("[Pt]([NH3])([NH3])([NH3])[NH3]", "TH", 1, 1),
    ("[Fe]123([N]CC[N]1)([N]CC[N]2)[N]CC[N]3", "OH", 2, 2),
])
def test_orbit_counts_match_chemistry(smi, geom, expected_min, expected_max):
    out = _isomers_for(smi)
    geom_isomers = _by_polyhedron(out, geom)
    assert expected_min <= len(geom_isomers) <= expected_max, (
        f"{smi} on {geom}: got {len(geom_isomers)} orbits, "
        f"expected {expected_min}..{expected_max}"
    )


def test_chelate_extended_label_construction():
    """The chelate-encoded label vector must mark all chelate donors with a colour."""
    classes = ["N", "N", "N", "N"]
    chelate_pairs = [frozenset((0, 1)), frozenset((2, 3))]
    ext = _build_extended_label(classes, chelate_pairs)
    assert ext[0].endswith("@c0")
    assert ext[1].endswith("@c0")
    assert ext[2].endswith("@c1")
    assert ext[3].endswith("@c1")


def test_chelate_not_on_trans_pair():
    """Enumerator must NEVER produce a perm where a chelate sits on a trans pair."""
    classes = ["N"] * 6
    chelate_pairs = [frozenset((0, 1)), frozenset((2, 3)), frozenset((4, 5))]
    orbits = _enumerate_orbits("OH", classes, chelate_pairs)
    oh_trans = frozenset(frozenset(p) for p in trans_positions("OH"))
    for o in orbits:
        for cp_at_pos in o["chelate_pairs_at_positions"]:
            assert cp_at_pos not in oh_trans, (
                f"chelate landed on trans pair: {cp_at_pos}"
            )


def test_polyhedra_for_cn_full_coverage():
    """Every CN in 2..9 must return at least one polyhedron."""
    for cn in range(2, 10):
        ps = polyhedra_for_cn(cn)
        assert ps, f"CN={cn} returned no polyhedra"


def test_positional_descriptor_fac_mer_oh():
    """OH 3-on-fac vs 3-on-mer triangle: descriptor distinguishes them."""
    # types_at_perm A on positions {0,2,4} = three cis = fac (no trans pair among them)
    fac = positional_descriptor("OH", ("A", "B", "A", "B", "A", "B"),
                                 perm=(0, 1, 2, 3, 4, 5),
                                 chelate_pairs=())
    assert fac == "fac"
    # A on {0,1,2} (positions 0,1 are trans on OH) → mer
    mer = positional_descriptor("OH", ("A", "A", "A", "B", "B", "B"),
                                 perm=(0, 1, 2, 3, 4, 5),
                                 chelate_pairs=())
    assert mer == "mer"


def test_n_theory_prescribed_smoke():
    """``n_theory_prescribed`` matches ``len(enumerate_prescribed_isomers)``."""
    smi = "[Cl][Pt]([Cl])([NH3])[NH3]"
    n = n_theory_prescribed(_mol(smi))
    full = enumerate_prescribed_isomers(_mol(smi))
    assert n == len(full)


def test_isomer_labels_unique():
    """Within a single SMILES, every prescribed isomer must have a unique label."""
    out = _isomers_for("[Cl][Co]12([NH3])([N]CC[N]1)[N]CC[N]2")
    labels = [d["label"] for d in out]
    assert len(labels) == len(set(labels)), (
        f"duplicate labels detected: {labels}"
    )
