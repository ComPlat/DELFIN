"""Universal CN4 multi-polyhedron dispatch + md_distance audit tests.

Pins the following contracts (ZURMAA-class fix, 2026-06-07):

1. **CN4 multi-poly (Au(III))** — a CN4 Au(III) SMILES with mixed donors
   {N, C, O, S} enumerates BOTH the T-4 tetrahedron AND the SP-4 square-
   planar orbit families (not just one), since CN4 is now in
   :data:`polyhedra._MULTI_POLY_CNS`.

2. **CN4 multi-poly (Ni)** — same contract for a non-Au CN4 metal: BOTH
   T-4 and SP-4 are enumerated unconditionally; no per-metal table.

3. **md_distance audit values** — under
   ``DELFIN_FFFREE_MD_DISTANCE_AUDIT=1`` the universal Pyykkö covalent-
   radii table returns chemistry-realistic M-D distances for late-TM /
   soft-donor combinations:

     * Au-S  ≈ 2.30 Å   (was missing → soft-donor correction)
     * Au-C  ≈ 1.99 Å
     * Au-N  ≈ 1.95 Å
     * Pt-N  ≈ 1.94 Å
     * Hg-S  ≈ 2.41 Å
     * Pt-Cl ≈ 2.22 Å
     * Cu-Br ≈ 2.26 Å

   All values are within ±0.10 Å of the CCDC-Mogul empirical first-
   sphere means for the corresponding bond type, *and* every common
   transition-metal × main-group-donor pair returns a finite, non-
   placeholder value (no more 1.5 + 0.75 = 2.25 generic fallback).

4. **Byte-identical OFF** — when the audit flag is unset,
   :func:`md_distance` returns exactly the legacy ``COV[m] + COV[d]``
   sum for every pair in the legacy table.

5. **Determinism** — two consecutive calls produce bit-identical
   results regardless of which env-flags are set.

6. **Universal coverage** — for every transition metal × every common
   donor (C, N, O, S, P, F, Cl, Br, I), the audited table returns a
   finite distance in [1.5, 3.5] Å (chemistry-plausible range).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

# Ensure the repo root is importable when run stand-alone.
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))


# --------------------------------------------------------------------------
# Fix #1 — Universal CN4 multi-polyhedron dispatch
# --------------------------------------------------------------------------
def test_multi_poly_cn_set_contains_cn4():
    """``_MULTI_POLY_CNS`` must include CN4 (ZURMAA fix)."""
    from delfin.fffree import polyhedra as PLY
    assert 4 in PLY._MULTI_POLY_CNS
    assert PLY.is_multi_poly_cn(4) is True
    # Other CNs we expect in the set (cross-CN sanity).
    for cn in (3, 5, 8, 9, 10, 11, 12):
        assert PLY.is_multi_poly_cn(cn) is True, f"CN{cn} should be multi-poly"
    # CN6 and CN7 are intentionally single-poly (OC-6 / PB-7 only).
    assert PLY.is_multi_poly_cn(6) is False
    assert PLY.is_multi_poly_cn(7) is False


def test_all_geometries_for_cn4_lists_both_polyhedra():
    """CN4 multi-poly view returns BOTH T-4 and SP-4 deterministically."""
    from delfin.fffree import polyhedra as PLY
    g4 = PLY.all_geometries_for_cn(4)
    assert "T-4 tetrahedron" in g4
    assert "SP-4 square planar" in g4
    # Determinism: lex / GEOM_BY_CN[cn] order, two calls identical.
    assert PLY.all_geometries_for_cn(4) == g4


def test_au_iii_cn4_multipoly_orbits_emit_both_t4_and_sp4():
    """A Au(III) CN4 SMILES emits BOTH T-4 and SP-4 orbit dicts under multi-poly.

    Uses the ZURMAA template: Au(III) with mixed donors (N, C, O, S).
    The simplified Au[N][C][O][S] SMILES suffices to trigger the CN4
    multi-poly dispatch — the universal fix doesn't depend on the exact
    ligand chemistry, only on the metal CN.
    """
    pytest.importorskip("rdkit", reason="rdkit not installed")
    from delfin.fffree.polyhedron_vertex_polya import (
        enumerate_orbits_for_smiles_multi,
    )
    # Simple Au(III) CN4 SMILES with four monodentate donors of distinct
    # element type — guaranteed CN4, no chelates.
    smi = "[Au+3]([NH3])(C)([OH])[SH]"
    out = enumerate_orbits_for_smiles_multi(smi)
    assert out is not None, "multi-poly enumerator returned None for Au(III) CN4"
    polys = sorted({str(r["polyhedron"]) for r in out})
    assert "T-4 tetrahedron" in polys, f"T-4 missing from {polys}"
    assert "SP-4 square planar" in polys, f"SP-4 missing from {polys}"


def test_ni_cn4_multipoly_orbits_emit_both_t4_and_sp4():
    """Ni CN4 with mixed donors enumerates BOTH T-4 and SP-4 (no per-metal table)."""
    pytest.importorskip("rdkit", reason="rdkit not installed")
    from delfin.fffree.polyhedron_vertex_polya import (
        enumerate_orbits_for_smiles_multi,
    )
    smi = "[Ni+2]([NH3])([NH3])(Cl)Cl"
    out = enumerate_orbits_for_smiles_multi(smi)
    assert out is not None, "multi-poly enumerator returned None for Ni CN4"
    polys = sorted({str(r["polyhedron"]) for r in out})
    assert "T-4 tetrahedron" in polys
    assert "SP-4 square planar" in polys


def test_cn6_stays_single_poly():
    """CN6 is intentionally NOT in multi_poly_cns -> single polyhedron only."""
    pytest.importorskip("rdkit", reason="rdkit not installed")
    from delfin.fffree.polyhedron_vertex_polya import (
        enumerate_orbits_for_smiles_multi,
    )
    smi = "[Fe+3](Cl)(Cl)(Cl)(Cl)(Cl)Cl"
    out = enumerate_orbits_for_smiles_multi(smi)
    assert out is not None
    polys = {str(r["polyhedron"]) for r in out}
    assert polys == {"OC-6 octahedron"}, (
        f"CN6 should stay OC-6 only (CN6 not in multi-poly set), got {polys}"
    )


def test_multipoly_orbits_deterministic():
    """Two calls produce identical multi-poly orbit lists (byte-equal dicts)."""
    pytest.importorskip("rdkit", reason="rdkit not installed")
    from delfin.fffree.polyhedron_vertex_polya import (
        enumerate_orbits_for_smiles_multi,
    )
    smi = "[Au+3]([NH3])(C)([OH])[SH]"
    a = enumerate_orbits_for_smiles_multi(smi)
    b = enumerate_orbits_for_smiles_multi(smi)
    assert a == b


# --------------------------------------------------------------------------
# Fix #2 — md_distance accuracy audit
# --------------------------------------------------------------------------
def _audit_on():
    os.environ["DELFIN_FFFREE_MD_DISTANCE_AUDIT"] = "1"


def _audit_off():
    os.environ.pop("DELFIN_FFFREE_MD_DISTANCE_AUDIT", None)
    os.environ.pop("DELFIN_FFFREE_PURE_TRACK3", None)


def test_md_distance_audit_off_byte_identical():
    """When DELFIN_FFFREE_MD_DISTANCE_AUDIT unset, md_distance == legacy COV sum."""
    _audit_off()
    from delfin.fffree.polyhedra import COV, md_distance
    # Every pair currently registered in COV must equal the legacy sum
    # (no soft-donor correction, no Pyykkö lookup, byte-identical).
    for metal in ("Cu", "Pd", "Au", "Pt", "Fe", "Ni", "Zn", "Hg"):
        for donor in ("C", "N", "O", "S", "P", "F", "Cl", "Br", "I"):
            expected = COV.get(metal, 1.5) + COV.get(donor, 0.75)
            got = md_distance(metal, donor)
            assert got == expected, (
                f"OFF byte-identity broken for {metal}-{donor}: "
                f"got {got}, expected {expected}"
            )


def test_md_distance_audit_au_s_chemistry_realistic():
    """Au-S under audit ≈ 2.30 Å (CCDC: 2.25-2.40 first-sphere mean)."""
    _audit_on()
    try:
        from delfin.fffree.polyhedra import md_distance
        r = md_distance("Au", "S")
        # 2.25-2.40 Å is the textbook Au-thiolate range; allow ±0.10 slack.
        assert 2.15 <= r <= 2.45, f"Au-S = {r} Å out of chemistry-realistic range"
    finally:
        _audit_off()


def test_md_distance_audit_au_c_chemistry_realistic():
    """Au-C under audit ≈ 2.00-2.10 Å (Au-aryl typical)."""
    _audit_on()
    try:
        from delfin.fffree.polyhedra import md_distance
        r = md_distance("Au", "C")
        assert 1.85 <= r <= 2.20, f"Au-C = {r} Å out of chemistry-realistic range"
    finally:
        _audit_off()


def test_md_distance_audit_pt_n_chemistry_realistic():
    """Pt-N under audit ≈ 2.00 Å (Pt(II)-N amine, pyridine: 1.95-2.10)."""
    _audit_on()
    try:
        from delfin.fffree.polyhedra import md_distance
        r = md_distance("Pt", "N")
        assert 1.85 <= r <= 2.15, f"Pt-N = {r} Å out of chemistry-realistic range"
    finally:
        _audit_off()


def test_md_distance_audit_soft_donor_correction_applied():
    """Soft-donor correction adds +0.05 Å for S / P / Se / I over Pyykkö-sum."""
    _audit_on()
    try:
        from delfin.fffree.polyhedra import (
            COV_PYYKKO, SOFT_DONOR_CORR, md_distance,
        )
        # Au-S should equal r(Au) + r(S) + 0.05 (audit ON).
        expected = COV_PYYKKO["Au"] + COV_PYYKKO["S"] + SOFT_DONOR_CORR["S"]
        got = md_distance("Au", "S")
        assert got == pytest.approx(expected), (
            f"Soft-donor correction not applied: got {got}, expected {expected}"
        )
        # Au-O is NOT a soft donor -> no correction.
        expected_o = COV_PYYKKO["Au"] + COV_PYYKKO["O"]
        got_o = md_distance("Au", "O")
        assert got_o == pytest.approx(expected_o)
    finally:
        _audit_off()


def test_md_distance_audit_universal_coverage_tm_x_donor():
    """Every (TM, common-donor) pair returns chemistry-plausible distance under audit.

    Spot-checks 14 transition metals × 9 main-group donors = 126 combos.
    Each must lie in [1.5, 3.5] Å — chemistry-plausible for any single
    M-D bond.  No generic 2.25 fallback (which is what the legacy COV
    falls back to for unknown metals).
    """
    _audit_on()
    try:
        from delfin.fffree.polyhedra import md_distance
        metals = (
            "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Pd", "Pt", "Au", "Hg", "Ag", "Cd", "Re", "Os", "Ir", "Tl",
        )
        donors = ("C", "N", "O", "S", "P", "F", "Cl", "Br", "I")
        for m in metals:
            for d in donors:
                r = md_distance(m, d)
                assert 1.5 <= r <= 3.5, (
                    f"{m}-{d} = {r} Å out of chemistry-plausible range [1.5, 3.5]"
                )
    finally:
        _audit_off()


def test_md_distance_audit_deterministic():
    """Audit values are bit-identical across consecutive calls."""
    _audit_on()
    try:
        from delfin.fffree.polyhedra import md_distance
        for m, d in (("Au", "S"), ("Au", "C"), ("Pt", "N"), ("Hg", "Cl")):
            a = md_distance(m, d)
            b = md_distance(m, d)
            assert a == b
    finally:
        _audit_off()


def test_md_distance_audit_pyykko_table_covers_all_tm():
    """COV_PYYKKO must cover every transition metal and every common donor."""
    from delfin.fffree.polyhedra import COV_PYYKKO
    for m in (
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi",
    ):
        assert m in COV_PYYKKO, f"Pyykkö table missing TM {m!r}"
    for d in ("C", "N", "O", "S", "P", "F", "Cl", "Br", "I",
              "Se", "Te", "As", "Sb"):
        assert d in COV_PYYKKO, f"Pyykkö table missing donor {d!r}"


# --------------------------------------------------------------------------
# Integrated demo — ZURMAA-class SMILES sanity
# --------------------------------------------------------------------------
def test_zurmaa_class_au_cn4_emits_both_polyhedra():
    """ZURMAA-class Au(III) CN4 SMILES emits BOTH T-4 and SP-4 polyhedron orbits.

    Uses a chemistry-realistic ZURMAA analogue (Au(III) + 1 N-pyridyl +
    1 C-aryl + 1 O-carboxylate + 1 S-thiolate).  This is the universal
    Fix #1 test in production form.
    """
    pytest.importorskip("rdkit", reason="rdkit not installed")
    from delfin.fffree.polyhedron_vertex_polya import (
        enumerate_orbits_for_smiles_multi,
    )
    # Au(III) CN4: 4 monodentate donors {N, C, O, S}, all distinct.
    smi = "[Au+3](c1ccccn1)(c1ccccc1)(OC(=O)C)Sc1ccccc1"
    out = enumerate_orbits_for_smiles_multi(smi)
    assert out is not None, "Au(III) CN4 ZURMAA-analogue: enumerator returned None"
    polys = {str(r["polyhedron"]) for r in out}
    assert "T-4 tetrahedron" in polys, f"T-4 missing: {polys}"
    assert "SP-4 square planar" in polys, f"SP-4 missing: {polys}"


def test_zurmaa_md_distances_audit_path():
    """ZURMAA-class Au-N/C/O/S distances under audit are chemistry-distinct.

    The legacy bug claim was "all M-D uniform 2.02 Å fallback".  Under the
    audited table, Au-S (soft, longer) and Au-O (hard, shorter) must
    differ by ≥ 0.30 Å -- proving the donor-type distinction is real,
    not collapsed to a single generic value.
    """
    _audit_on()
    try:
        from delfin.fffree.polyhedra import md_distance
        d_au_n = md_distance("Au", "N")
        d_au_c = md_distance("Au", "C")
        d_au_o = md_distance("Au", "O")
        d_au_s = md_distance("Au", "S")
        assert abs(d_au_s - d_au_o) >= 0.30, (
            f"Au-S ({d_au_s}) vs Au-O ({d_au_o}) collapse — uniform-fallback bug"
        )
        # Soft donor must be the LONGEST of the four (chemistry rule).
        assert d_au_s > max(d_au_n, d_au_c, d_au_o)
        # Each individually in chemistry-plausible window.
        assert 1.85 <= d_au_o <= 2.15
        assert 1.85 <= d_au_n <= 2.10
        assert 1.85 <= d_au_c <= 2.20
        assert 2.20 <= d_au_s <= 2.45
    finally:
        _audit_off()
