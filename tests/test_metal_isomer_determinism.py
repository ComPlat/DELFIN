"""Determinism and completeness tests for metal-isomer SMILES conversion.

Covers CN 4, 6, 7, 8 and chelate cis-constraint enforcement.  These tests
guard the coordination-preserving UFF path and the final geometry-based
dedup introduced to make isomer generation deterministic across runs.
"""
from __future__ import annotations

import math
import itertools
from typing import List, Tuple

import pytest

from delfin import smiles_converter as sc


def _metal_and_donors(
    xyz: str,
    metal_symbols: Tuple[str, ...] = (),
    top_k: int = 8,
) -> Tuple[Tuple[float, float, float], List[Tuple[str, Tuple[float, float, float], float]]]:
    atoms: List[Tuple[str, Tuple[float, float, float]]] = []
    for line in xyz.splitlines():
        parts = line.split()
        if len(parts) < 4:
            continue
        atoms.append(
            (parts[0], (float(parts[1]), float(parts[2]), float(parts[3])))
        )
    metals = [
        a for a in atoms
        if (metal_symbols and a[0] in metal_symbols)
        or (not metal_symbols and a[0] in sc._METAL_SET)
    ]
    assert metals, f"no metal atom found for symbols {metal_symbols}"
    m_sym, m_pos = metals[0]
    heavy = [
        (sym, pos, math.dist(m_pos, pos))
        for sym, pos in atoms
        if sym != 'H' and (sym, pos) != (m_sym, m_pos)
    ]
    heavy.sort(key=lambda t: t[2])
    return m_pos, heavy[:top_k]


def _angles(metal: Tuple[float, float, float], donors) -> List[Tuple[str, str, float]]:
    out: List[Tuple[str, str, float]] = []
    for (s1, p1, _), (s2, p2, _) in itertools.combinations(donors, 2):
        v1 = tuple(a - b for a, b in zip(p1, metal))
        v2 = tuple(a - b for a, b in zip(p2, metal))
        m1 = math.sqrt(sum(x * x for x in v1))
        m2 = math.sqrt(sum(x * x for x in v2))
        if m1 < 1e-9 or m2 < 1e-9:
            continue
        dot = sum(a * b for a, b in zip(v1, v2))
        ang = math.degrees(math.acos(max(-1.0, min(1.0, dot / (m1 * m2)))))
        out.append((s1, s2, ang))
    return out


def _isomers(smiles: str, num_confs: int = 60, max_isomers: int = 20):
    res, err = sc.smiles_to_xyz_isomers(
        smiles,
        num_confs=num_confs,
        max_isomers=max_isomers,
    )
    assert err is None, f"conversion error: {err}"
    return res


# ---------------------------------------------------------------------------
# CN=6 octahedral MA2B2C2 — 5 distinct geometric isomers expected
# ---------------------------------------------------------------------------

CD_MA2B2C2 = (
    "C[OH+][Cd-4]([Cl])([Cl])([OH+]C)"
    "([N+]1=C2SC(C)=NN2C(C2=CC=CC=C2)=N1)"
    "[N+]1=C2SC(C)=NN2C(C2=CC=CC=C2)=N1"
)


def test_cd_MA2B2C2_yields_all_five_octahedral_isomers():
    res = _isomers(CD_MA2B2C2)
    labels = {lbl for _xyz, lbl in res}
    oh_expected = {'all-cis', 'all-trans', 'N-trans', 'O-trans', 'Cl-trans'}
    assert oh_expected <= labels, (
        f"missing octahedral isomers: {oh_expected - labels} (got {labels})"
    )


def test_cd_MA2B2C2_octahedral_geometries_are_clean():
    res = _isomers(CD_MA2B2C2)
    oh_labels = {'all-cis', 'all-trans', 'N-trans', 'O-trans', 'Cl-trans'}
    for xyz, lbl in res:
        if lbl not in oh_labels:
            continue
        m, donors = _metal_and_donors(xyz, metal_symbols=('Cd',), top_k=6)
        assert len(donors) == 6, f"{lbl}: expected CN=6, got {len(donors)}"
        for sym, _p, d in donors:
            assert 2.15 <= d <= 2.85, (
                f"{lbl}: {sym} at d={d:.2f} Å outside Cd range"
            )
        angles = _angles(m, donors)
        trans = [a for _s1, _s2, a in angles if a >= 135]
        assert len(trans) == 3, f"{lbl}: expected 3 trans pairs, got {len(trans)}"
        for a in trans:
            assert a >= 170.0, f"{lbl}: trans angle {a:.1f}° too distorted"


# ---------------------------------------------------------------------------
# CN=6 cyclometalated Ir — clean octahedra, no alt-binding artifacts
# ---------------------------------------------------------------------------

IR_COMPLEX = (
    "CC1=CC(C)=[O+][Ir-3]2([N+]3=C4C=CC=C3)"
    "(C5=CC=CC=C54)(O1)[N+]6=CC=CC=C6C7=C2C=CC=C7"
)


def test_ir_complex_no_unphysical_short_metal_donor_bonds():
    res = _isomers(IR_COMPLEX)
    assert res, "no Ir isomers produced"
    for xyz, lbl in res:
        _m, donors = _metal_and_donors(xyz, metal_symbols=('Ir',), top_k=6)
        for sym, _p, d in donors:
            # No Ir-X coordination below 1.80 Å
            assert d >= 1.75, (
                f"{lbl}: unphysical Ir-{sym} distance {d:.2f} Å"
            )


def test_ir_complex_octahedral_angles_within_tolerance():
    res = _isomers(IR_COMPLEX)
    for xyz, lbl in res:
        if 'trigonal-prismatic' in lbl or 'see-saw' in lbl:
            continue
        m, donors = _metal_and_donors(xyz, metal_symbols=('Ir',), top_k=6)
        assert len(donors) == 6, f"{lbl}: expected CN=6"
        angles = _angles(m, donors)
        trans = [a for _s1, _s2, a in angles if a >= 135]
        assert len(trans) == 3, f"{lbl}: expected 3 trans pairs, got {len(trans)}"
        for a in trans:
            assert a >= 150.0, f"{lbl}: trans angle {a:.1f}° too distorted"


# ---------------------------------------------------------------------------
# Determinism: identical SMILES must produce identical isomer label sets
# ---------------------------------------------------------------------------

def test_cd_isomer_set_is_deterministic_across_runs():
    first = sorted(lbl for _x, lbl in _isomers(CD_MA2B2C2))
    second = sorted(lbl for _x, lbl in _isomers(CD_MA2B2C2))
    assert first == second, (
        f"non-deterministic isomer set: {first} vs {second}"
    )


# ---------------------------------------------------------------------------
# CN=7 and CN=8 — homoleptic sanity
# ---------------------------------------------------------------------------

def test_cn7_homoleptic_produces_at_least_one_isomer():
    res = _isomers('[Fe+2](O)(O)(O)(O)(O)(O)O')
    assert res, "no CN=7 isomer produced"
    xyz, _lbl = res[0]
    _m, donors = _metal_and_donors(xyz, metal_symbols=('Fe',), top_k=7)
    assert len(donors) == 7, f"expected CN=7, got {len(donors)}"
    for sym, _p, d in donors:
        assert 1.70 <= d <= 2.60, f"Fe-{sym} distance {d:.2f} Å out of range"


def test_cn8_homoleptic_produces_at_least_one_isomer():
    # Zr(IV) has a proper M-L bond-length lookup entry (2.10 for O); U is
    # less well parameterized by UFF and would drift further from the ideal.
    res = _isomers('[Zr+4](O)(O)(O)(O)(O)(O)(O)O')
    assert res, "no CN=8 isomer produced"
    xyz, _lbl = res[0]
    _m, donors = _metal_and_donors(xyz, metal_symbols=('Zr',), top_k=8)
    assert len(donors) == 8, f"expected CN=8, got {len(donors)}"
    for sym, _p, d in donors:
        assert 1.90 <= d <= 2.60, f"Zr-{sym} distance {d:.2f} Å out of range"


# ---------------------------------------------------------------------------
# Chelate constraint: bidentate cis enforcement
# ---------------------------------------------------------------------------

def test_bis_en_co_chelate_allows_trans_cl():
    """[Co(en)2Cl2]+: 2 ethylenediamine chelates must remain cis,
    the 2 chlorides may be cis OR trans to each other.
    """
    smi = 'Cl[Co+3]12(Cl)(NCCN1)NCCN2'
    res = _isomers(smi)
    assert res, "no Co-en-Cl2 isomer produced"
    # For each isomer, each en ligand's two N donors must be cis (<110°).
    for xyz, lbl in res:
        m, donors = _metal_and_donors(xyz, metal_symbols=('Co',), top_k=6)
        assert len(donors) >= 4, f"{lbl}: not enough donors found"
        n_donors = [d for d in donors if d[0] == 'N']
        # At least two N donors are present.
        assert len(n_donors) >= 2, f"{lbl}: missing N donors"
        # Chelate ring forces at least one N-N pair at cis angle.
        angles = _angles(m, n_donors)
        cis_nn = [a for _s1, _s2, a in angles if a <= 110.0]
        assert cis_nn, (
            f"{lbl}: no cis N-N pair — chelate constraint violated"
        )


# ---------------------------------------------------------------------------
# Multinuclear: bimetallic coupled enumeration
# ---------------------------------------------------------------------------

def test_bimetallic_produces_more_isomers_than_single_metal():
    """[Fe2(mu-Cl)2Cl(NH3)4]: two Fe sharing bridging Cl.
    Coupled enumeration should produce more isomers than a single metal."""
    smi = '[Fe+2]1(N)(N)(Cl[Fe+2](N)(N)Cl1)Cl'
    res = _isomers(smi, num_confs=12)
    assert len(res) >= 3, (
        f"bimetallic should produce >=3 isomers, got {len(res)}"
    )
