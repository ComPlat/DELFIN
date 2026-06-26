"""Tests for the FF-free donor-aware metal-centre resolution (delfin.manta.decompose).

Root cause locked here (eye-flagged ATOQUV, a Pd bis(distibine)): ``decompose``
identifies the coordination CENTRE with ``_bond_decollapse._is_metal``, which also
flags the metalloid / heavier-pnictogen-chalcogen / post-transition elements
(Sb / As / Bi / Te / Se / Ge / Sn / Pb).  A Pd coordinated by four Sb donors is
therefore counted as FIVE "metals" -> ``len(metals) != 1`` -> ``decompose`` bails to
legacy and the donors are never seated (the ligands float free at M-Sb ~3.2 A, CN 0).
The SAME two distibine ligands on a Pt centre are seated correctly by a metal-dependent
legacy branch (M-Sb 3.05 A, CN 4), so the classification was INCONSISTENT across the
central metal.

Fix (DELFIN_FFFREE_METALLOID_DONOR, default OFF -> byte-identical): when the historic
metal count is not 1, resolve the unique true (non-metalloid) transition/lanthanide
centre and treat the metalloid neighbours as DONORS via the existing Werner cleave, so
the complex builds the SAME way regardless of the central metal.

Contracts:
  * graph/element-rule resolver: a single true metal with metalloid donors bonded to it
    resolves to that metal; a metalloid cluster / a real bimetallic does NOT,
  * with the flag ON the Pd bis(distibine) decomposes natively (Pd CN4 SP-4, two Sb
    chelates) exactly like its Pt twin -> consistent across the central metal,
  * env-gate: flag unset/0 -> byte-identical (Sb-donor complex stays legacy, normal
    Werner complexes unchanged even with the flag ON),
  * determinism: same input -> byte-identical output,
  * the heavy metalloid M-D distance is realistic (Pyykko radius, not the 0.75 A
    placeholder) ONLY under the flag.
"""
import os

import pytest

from rdkit import Chem

from delfin.manta import decompose as DEC
from delfin.manta import polyhedra as PLY
from delfin.manta.converter_backend import _fffree_isomers

# ATOQUV: Pd(II) coordinated by two o-phenylene-bis(dimethylstibine) chelates (4 Sb).
ATOQUV_PD = ("[CH3][Sb+]1([CH3])[CH2]C2=CC=CC=C2[CH2][Sb+]([CH3])([CH3])[Pd-2]12"
             "[Sb+]([CH3])([CH3])[CH2]C1=CC=CC=C1[CH2][Sb+]2([CH3])[CH3]")
# ATOQOP: the SAME two distibine ligands on a Pt centre (the "perfect twin").
ATOQOP_PT = ("[CH3][Sb+]1([CH3])[CH2]C2=CC=CC=C2[CH2][Sb+]([CH3])([CH3])[Pt-2]12"
             "[Sb+]([CH3])([CH3])[CH2]C1=CC=CC=C1[CH2][Sb+]2([CH3])[CH3]")
# AKEVIX: an Sb cluster around Y (metalloid-only ring, no clean single donor) -> legacy.
AKEVIX = ("[Sb]12[Sb+]34[Sb]5[Sb+]16[Sb+]17[Sb+]8[Sb+]9%10[Sb+]1[Y-12]"
          "25836791%11%12[Sb]2[Sb+]41[Sb+]%11[Sb+]2%10%12")
# Non-metalloid Werner controls (must be untouched / never enter the resolver branch).
CISPLATIN = "N[Pt](N)(Cl)Cl"
COCL3NH3 = "[NH3][Co]([NH3])([NH3])([Cl])([Cl])[Cl]"
# A genuine bimetallic (two true TM) -> must stay legacy.
FE_FE = "[NH3][Fe]([NH3])([NH3])([Fe]([NH3])([NH3])[NH3])([NH3])[NH3]"

# To EMIT a coordinated frame the Sb-methyl ligand also needs the existing X-ray-short
# C-H collapse calibration (orthogonal coverage gate, not part of this fix): one methyl
# H lands at a short distance the historic 0.82*ideal floor false-flags.
_BUILD_FLAGS = {
    "DELFIN_FFFREE_BUILDER": "1",
    "DELFIN_FFFREE_METALLOID_DONOR": "1",
    "DELFIN_FFFREE_XH_COLLAPSE": "1",
    "PYTHONHASHSEED": "0",
}


@pytest.fixture
def clean_env(monkeypatch):
    """Remove every FF-free flag so each test sets exactly what it needs."""
    for k in list(os.environ):
        if k.startswith("DELFIN_FFFREE_"):
            monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("PYTHONHASHSEED", "0")
    yield monkeypatch


def _prep(smiles):
    from delfin.smiles_converter import _prepare_mol_for_embedding
    mol = _prepare_mol_for_embedding(smiles)
    if mol is None:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
    return mol


# --------------------------------------------------------------------------- #
# graph/element-rule resolver
# --------------------------------------------------------------------------- #
def test_resolver_picks_true_metal_over_metalloid(clean_env):
    clean_env.setenv("DELFIN_FFFREE_METALLOID_DONOR", "1")
    mol = _prep(ATOQUV_PD)
    c = DEC._resolve_metal_center(mol)
    assert c is not None
    assert mol.GetAtomWithIdx(c).GetSymbol() == "Pd"


def test_resolver_rejects_metalloid_cluster(clean_env):
    # AKEVIX: Y wrapped in an Sb cluster — the Sb are not clean single donors -> legacy.
    clean_env.setenv("DELFIN_FFFREE_METALLOID_DONOR", "1")
    mol = _prep(AKEVIX)
    # either no clean single resolution, or decompose declines it -> legacy
    assert DEC.decompose(AKEVIX) is None


def test_resolver_rejects_real_bimetallic(clean_env):
    clean_env.setenv("DELFIN_FFFREE_METALLOID_DONOR", "1")
    mol = _prep(FE_FE)
    assert DEC._resolve_metal_center(mol) is None
    assert DEC.decompose(FE_FE) is None


# --------------------------------------------------------------------------- #
# consistent native decomposition across the central metal
# --------------------------------------------------------------------------- #
def test_distibine_decomposes_consistently_with_flag(clean_env):
    clean_env.setenv("DELFIN_FFFREE_METALLOID_DONOR", "1")
    d_pd = DEC.decompose(ATOQUV_PD)
    d_pt = DEC.decompose(ATOQOP_PT)
    assert d_pd is not None and d_pt is not None
    # SAME coordination decomposition regardless of the central metal
    for d in (d_pd, d_pt):
        assert d["cn"] == 4
        assert d["geometry"] == "SP-4 square planar"
        assert d["has_chelate"] is True
        assert d["donor_elems"] == ["Sb", "Sb"]
    assert d_pd["metal"] == "Pd"
    assert d_pt["metal"] == "Pt"


def test_distibine_builds_coordinated_frame(clean_env):
    for k, v in _BUILD_FLAGS.items():
        clean_env.setenv(k, v)
    import numpy as np
    iso = _fffree_isomers(ATOQUV_PD)
    assert iso is not None and len(iso) >= 1
    xyz = iso[0][0]
    rows = [ln.split() for ln in xyz.strip().splitlines() if len(ln.split()) == 4]
    syms = [r[0] for r in rows]
    P = np.array([[float(x) for x in r[1:4]] for r in rows], float)
    mi = syms.index("Pd")
    sb = [i for i, s in enumerate(syms) if s == "Sb"]
    dists = sorted(float(np.linalg.norm(P[i] - P[mi])) for i in sb)
    # all four Sb donors coordinated at a realistic Pd-Sb distance (~2.6-2.9 A),
    # NOT the uncoordinated legacy float (~3.2 A) the OFF path produces.
    assert sum(1 for x in dists if x < 3.0) == 4
    assert all(2.4 < x < 3.0 for x in dists)


# --------------------------------------------------------------------------- #
# heavy-metalloid M-D distance is realistic only under the flag
# --------------------------------------------------------------------------- #
def test_metalloid_md_distance_flag_gated(clean_env):
    # OFF: Sb missing from COV -> 0.75 A placeholder (Pd-Sb ~ 2.14 A, too short).
    assert PLY.md_distance("Pd", "Sb") == pytest.approx(
        PLY.COV["Pd"] + 0.75, abs=1e-6)
    # ON: published Pyykko Sb radius -> realistic Pd-Sb ~ 2.79 A.
    clean_env.setenv("DELFIN_FFFREE_METALLOID_DONOR", "1")
    on = PLY.md_distance("Pd", "Sb")
    assert on == pytest.approx(PLY.COV["Pd"] + PLY._METALLOID_COV["Sb"], abs=1e-6)
    assert on > 2.6


# --------------------------------------------------------------------------- #
# env-gate: flag unset/0 -> byte-identical
# --------------------------------------------------------------------------- #
def test_legacy_when_flag_unset(clean_env):
    # No DELFIN_FFFREE_METALLOID_DONOR -> the distibine complexes stay legacy (None).
    clean_env.setenv("DELFIN_FFFREE_BUILDER", "1")
    assert DEC.decompose(ATOQUV_PD) is None
    assert DEC.decompose(ATOQOP_PT) is None


def test_flag_off_byte_identical_on_distibine(clean_env):
    # unset vs explicit 0 must give the identical (None) decompose result.
    unset = DEC.decompose(ATOQUV_PD)
    clean_env.setenv("DELFIN_FFFREE_METALLOID_DONOR", "0")
    explicit0 = DEC.decompose(ATOQUV_PD)
    assert unset is None and explicit0 is None


def test_flag_on_inert_on_non_metalloid(clean_env):
    clean_env.setenv("DELFIN_FFFREE_BUILDER", "1")
    # baseline (flag off) on normal Werner complexes
    base = _fffree_isomers(CISPLATIN)
    base_co = _fffree_isomers(COCL3NH3)
    assert base is not None and base_co is not None
    # flag on must not change a non-metalloid build (resolver branch never entered)
    clean_env.setenv("DELFIN_FFFREE_METALLOID_DONOR", "1")
    on = _fffree_isomers(CISPLATIN)
    on_co = _fffree_isomers(COCL3NH3)
    assert [x[0] for x in on] == [x[0] for x in base]
    assert [x[0] for x in on_co] == [x[0] for x in base_co]


# --------------------------------------------------------------------------- #
# determinism
# --------------------------------------------------------------------------- #
def test_distibine_build_is_deterministic(clean_env):
    for k, v in _BUILD_FLAGS.items():
        clean_env.setenv(k, v)
    a = _fffree_isomers(ATOQUV_PD)
    b = _fffree_isomers(ATOQUV_PD)
    assert a is not None
    assert [x[0] for x in a] == [x[0] for x in b]
