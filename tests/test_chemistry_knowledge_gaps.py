"""Tests for chemistry-knowledge layer gaps (Welle-3 Agent T6.2 audit).

Validates the *current* coverage of the chemistry tables in
``delfin.smiles_converter`` and pins the specific gaps identified in
``agent_workspace/quality_framework/iters/welle3_T6.2_chemistry_gaps_20260515.md``
so any future expansion is exercised by the unit suite.

These tests are designed to be **failing-by-design** for the documented
gaps (``xfail``) and **passing** for the parts of the knowledge layer
that already work.  When a gap is closed by a small data-patch (e.g. a
new ``_METAL_GROUP_NUMBER`` entry, a new ``_METAL_METAL_BOND_LENGTHS``
key) the corresponding ``xfail`` flips to a passing case.

All assertions reference master_v3 statistics:

  * top-150 (M, donor) pairs              — graph-level edge count
  * 9 unique direct M-M edges             — Ir-Sn, Rh-Sn, Os-Sn, Ru-Sn,
                                              Co-Sn, Fe-Sn, In-Ru, In-Ir,
                                              Co-In  (smiles_master_v3)
  * 14 lanthanides + 6 actinides          — Pu/Th/U/Np density 0.x %

No worktree side-effects: pure read-only inspection of module globals.
"""
from __future__ import annotations

import pytest

_sc = pytest.importorskip(
    "delfin.smiles_converter",
    reason="smiles_converter not importable",
)

_METAL_LIGAND_BOND_LENGTHS = getattr(_sc, "_METAL_LIGAND_BOND_LENGTHS", None)
_METAL_METAL_BOND_LENGTHS = getattr(_sc, "_METAL_METAL_BOND_LENGTHS", None)
_METAL_GROUP_NUMBER = getattr(_sc, "_METAL_GROUP_NUMBER", None)
_PI_ACCEPTOR_DONOR_ELEMS = getattr(_sc, "_PI_ACCEPTOR_DONOR_ELEMS", None)
_HAPTO_CENTROID_DISTANCES = getattr(_sc, "_HAPTO_CENTROID_DISTANCES", None)
_cn5_d_electron_count = getattr(_sc, "_cn5_d_electron_count", None)

if any(
    obj is None
    for obj in (
        _METAL_LIGAND_BOND_LENGTHS,
        _METAL_METAL_BOND_LENGTHS,
        _METAL_GROUP_NUMBER,
        _PI_ACCEPTOR_DONOR_ELEMS,
        _HAPTO_CENTROID_DISTANCES,
        _cn5_d_electron_count,
    )
):
    pytest.skip(
        "chemistry-knowledge surface not exported",
        allow_module_level=True,
    )


# --------------------------------------------------------------------------
# 1. M-L bond-length coverage of master_v3 top-150 pairs
# --------------------------------------------------------------------------
# Empirical: zero misses out of the 150 highest-frequency (M, donor)
# graph edges in pools/smiles_master_v3.txt (11 363 SMILES).  This guards
# against accidental deletion of entries that downstream code (UFF / OB
# pin) depends on.
TOP_PAIRS_MASTER_V3 = [
    ("Ru", "C"), ("Ru", "N"), ("Co", "N"), ("Cu", "N"), ("Ni", "N"),
    ("Fe", "N"), ("Ir", "C"), ("Ru", "P"), ("Re", "C"), ("Rh", "C"),
    ("Zn", "N"), ("Mn", "N"), ("Re", "N"), ("Mo", "C"), ("W", "C"),
    ("Ru", "Cl"), ("Pt", "N"), ("Mo", "N"), ("Cd", "N"), ("Rh", "N"),
    ("Os", "C"), ("Fe", "C"), ("Y", "N"), ("Cu", "O"), ("Ir", "P"),
    ("Zr", "N"), ("Rh", "P"), ("Ir", "N"), ("Pt", "C"), ("Pd", "N"),
    ("Y", "O"), ("Co", "O"), ("Cd", "O"), ("Ti", "N"), ("Mo", "O"),
    ("Cr", "N"), ("Os", "N"), ("Ni", "O"), ("Re", "O"), ("Zn", "O"),
    ("Mn", "C"), ("Fe", "O"), ("Ru", "O"), ("W", "N"), ("Cr", "C"),
    ("Ni", "P"), ("Re", "P"), ("Mn", "O"), ("Fe", "P"), ("Pd", "C"),
]


@pytest.mark.parametrize("metal,donor", TOP_PAIRS_MASTER_V3)
def test_metal_ligand_top50_present(metal: str, donor: str) -> None:
    assert (metal, donor) in _METAL_LIGAND_BOND_LENGTHS, (
        f"({metal!r}, {donor!r}) missing from _METAL_LIGAND_BOND_LENGTHS"
    )


# --------------------------------------------------------------------------
# 2. M-M bond-length coverage gap
# --------------------------------------------------------------------------
# The 24-key table covers homometallic Cr-Cr / Mo-Mo (multiple bonds),
# common homometallic dimers, and 4 heterometallic pairs.  master_v3
# contains direct M-M edges to Sn (TM-Sn stannyl) and In (TM-In indyl)
# that are *not* in the table — the fallback path produces 3.0+ Å which
# breaks downstream M-M validation.
TM_SN_PAIRS = [("Ir", "Sn"), ("Rh", "Sn"), ("Os", "Sn"), ("Ru", "Sn"),
               ("Co", "Sn"), ("Fe", "Sn")]
TM_IN_PAIRS = [("In", "Ru"), ("In", "Ir"), ("Co", "In")]


@pytest.mark.xfail(
    strict=True,
    reason="Gap G2: TM-Sn / TM-In heteronuclear M-M pairs absent "
    "(see welle3_T6.2 report).  Proposed patch adds 9 frozenset keys.",
)
@pytest.mark.parametrize("a,b", TM_SN_PAIRS + TM_IN_PAIRS)
def test_metal_metal_tm_main_group_pairs(a: str, b: str) -> None:
    assert frozenset({a, b}) in _METAL_METAL_BOND_LENGTHS


def test_metal_metal_homometallic_count() -> None:
    homo_keys = [k for k in _METAL_METAL_BOND_LENGTHS if len(k) == 1]
    assert len(homo_keys) >= 20, (
        f"expected >=20 homometallic entries, got {len(homo_keys)}"
    )


# --------------------------------------------------------------------------
# 3. d-electron count helper — lanthanide & actinide gaps
# --------------------------------------------------------------------------
LANTHANIDES_MAIN = ["Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb",
                    "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
ACTINIDES_MAIN = ["Th", "Pa", "U", "Np", "Pu"]


@pytest.mark.xfail(
    strict=True,
    reason="Gap G3a: lanthanides Ce-Lu return d=None (group entry "
    "missing).  Patch: add group-3 entries with note d=group-ox, where "
    "the f-block contribution is intentionally folded into a constant.",
)
@pytest.mark.parametrize("sym", LANTHANIDES_MAIN)
def test_d_electron_count_lanthanide_returns_int(sym: str) -> None:
    d = _cn5_d_electron_count(sym, 3)
    assert d is not None and 0 <= d <= 10


@pytest.mark.xfail(
    strict=True,
    reason="Gap G3b: actinides Th, Pa, U, Np, Pu — no group entry; "
    "Patch: add explicit group numbers from periodic-table convention.",
)
@pytest.mark.parametrize("sym", ACTINIDES_MAIN)
def test_d_electron_count_actinide_returns_int(sym: str) -> None:
    d = _cn5_d_electron_count(sym, 4)
    assert d is not None and 0 <= d <= 10


def test_d_electron_count_typical_tm_cases() -> None:
    # Pt(II) d8, Ni(II) d8, Cu(II) d9, V(V) d0, Cr(III) d3, Fe(III) d5,
    # Co(III) d6, Mn(II) d5, Zn(II) d10
    assert _cn5_d_electron_count("Pt", 2) == 8
    assert _cn5_d_electron_count("Ni", 2) == 8
    assert _cn5_d_electron_count("Cu", 2) == 9
    assert _cn5_d_electron_count("V", 5) == 0
    assert _cn5_d_electron_count("Cr", 3) == 3
    assert _cn5_d_electron_count("Fe", 3) == 5
    assert _cn5_d_electron_count("Co", 3) == 6
    assert _cn5_d_electron_count("Mn", 2) == 5
    assert _cn5_d_electron_count("Zn", 2) == 10


def test_d_electron_count_main_group_is_none() -> None:
    """Main-group metals correctly return None (no d-block chemistry)."""
    for sym in ("Al", "Ga", "In", "Tl", "Sn", "Pb", "Bi",
                "Li", "Na", "K", "Mg", "Ca", "Sr", "Ba", "Be"):
        assert _cn5_d_electron_count(sym, 2) is None, (
            f"{sym} unexpectedly returned a d-count"
        )


# --------------------------------------------------------------------------
# 4. Donor-element π-acceptor classification
# --------------------------------------------------------------------------
def test_pi_acceptor_includes_p_as_sb() -> None:
    for e in ("P", "As", "Sb"):
        assert e in _PI_ACCEPTOR_DONOR_ELEMS


def test_sigma_donors_excluded_from_pi_acceptor_label_table() -> None:
    """Pure σ donors N / O / S / halides must NOT be flagged as π-acceptors.

    NHC carbene C and CO are detected graph-based in
    :func:`_classify_cn5_geometry` so they are *intentionally* absent
    from the static element list.
    """
    for e in ("N", "O", "S", "F", "Cl", "Br", "I", "C"):
        assert e not in _PI_ACCEPTOR_DONOR_ELEMS


# --------------------------------------------------------------------------
# 5. Hapto centroid-distance coverage
# --------------------------------------------------------------------------
# Common Cp* / arene anchors (eta-5, eta-6) for d-block metals must be
# present; eta-3 (allyl) needed for Co, Ni, Pd, Pt, Rh, Ir, Mo.
HAPTO_CRITICAL = [
    ("Fe", 5), ("Ru", 5), ("Os", 5), ("Cr", 5), ("Mo", 5),
    ("W", 5), ("Re", 5), ("Mn", 5), ("Ti", 5), ("Zr", 5),
    ("Hf", 5), ("V", 5), ("Co", 5), ("Rh", 5), ("Ir", 5), ("Ni", 5),
    ("Cr", 6), ("Mo", 6), ("W", 6), ("Fe", 6), ("Ru", 6),
    ("Os", 6), ("Mn", 6), ("V", 6),
]


@pytest.mark.parametrize("metal,eta", HAPTO_CRITICAL)
def test_hapto_centroid_distance_known(metal: str, eta: int) -> None:
    assert (metal, eta) in _HAPTO_CENTROID_DISTANCES, (
        f"({metal!r}, {eta}) hapto centroid distance missing"
    )


# eta-2 (alkene) coverage gap for several common metals.
HAPTO_ETA2_GAP = ["Cu", "Ag", "Au", "Mo", "W", "Mn", "Os", "Ti", "V"]


@pytest.mark.xfail(
    strict=False,
    reason="Gap G4: η^2 (alkene π) centroid distances absent for many "
    "TMs that form Zeise-type Cu(I)/Ag(I)/Au(I) alkene adducts.",
)
@pytest.mark.parametrize("sym", HAPTO_ETA2_GAP)
def test_hapto_eta2_centroid_known(sym: str) -> None:
    assert (sym, 2) in _HAPTO_CENTROID_DISTANCES
