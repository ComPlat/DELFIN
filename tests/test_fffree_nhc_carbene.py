"""Tests for the FF-free NHC / carbene-C donor repair (delfin.manta.decompose).

Root cause locked here: a metal-bound carbene carbon (NHC = N-C-N imidazol-2-
yliden and kin) is encoded in the dative-bond SMILES as ``[C+]`` (charge +1,
valence 3).  After the Werner M-C cleave it drops to degree 2 but keeps the
``[C+]``/aromatic encoding artifact, so RDKit cannot kekulise the surrounding ring
(KekulizeException) -> the whole fragment split fails -> the system falls back to
legacy.  The fix repairs the cleaved carbene-C to the RDKit-valid representation of
a FREE singlet carbene ``:C:`` (neutral, valence 2, two non-bonding radical
electrons) and retries the split.

Contracts:
  * graph-only carbene detection: an NHC carbene-C is flagged, a lone sigma-C
    (methyl / carbonyl C) is NOT,
  * the repaired free carbene sanitizes (no KekulizeException) and is neutral
    divalent (valence 2),
  * with the flag ON an NHC complex decomposes natively (carbene-C kept as a
    C sigma-donor, denticity 1) and builds >=1 isomer,
  * env-gate: flag unset/0 -> byte-identical (NHC stays legacy, non-NHC unchanged),
  * determinism: same input -> byte-identical output.
"""
import os

import pytest

from rdkit import Chem

from delfin.manta import decompose as DEC
from delfin.manta.converter_backend import _fffree_isomers

# A small single-NHC silver(I) complex (FUFYOW): Ag bound to an acetate-O and one
# benzannulated NHC carbene-C.  CN2 (d10 Ag(I) linear).
FUFYOW = ("CC(=O)[O][Ag-][C+]1N(C)C2=C(C(=O)N(C)C(=O)N2C)N1CC1=CC=CC=C1")
# A bis-NHC silver(I) (ACIPAG): two benzimidazol-2-yliden carbenes, CN2.
ACIPAG = ("CCCCCCCCCCN1C2=CC=CC=C2N(CCCCCCCCCC)[C+]1[Ag-][C+]1N(CCCCCCCCCC)"
          "C2=CC=CC=C2N1CCCCCCCCCC")
# Non-NHC controls (must be untouched / never enter the repair branch).
CISPLATIN = "N[Pt](N)(Cl)Cl"
PTME2CL2 = "C[Pt](C)(Cl)Cl"            # sigma methyl-C, not a carbene
COCL3NH3 = "[NH3][Co]([NH3])([NH3])([Cl])([Cl])[Cl]"

# Flags needed for these particular NHC to *build* natively (orthogonal coverage
# gates, not part of the carbene fix): Ag-NHC is CN2 (linear) with bulky ligands.
_BUILD_FLAGS = {
    "DELFIN_FFFREE_NHC_CARBENE": "1",
    "DELFIN_FFFREE_CN_EXTEND": "1",
    "DELFIN_FFFREE_JOINT_DECLASH": "1",
    "DELFIN_FFFREE_MONO_HEAVY_CAP": "30",
    "PYTHONHASHSEED": "0",
}


@pytest.fixture
def clean_env(monkeypatch):
    """Remove every FF-free flag so each test sets exactly what it needs."""
    for k in list(os.environ):
        if k.startswith("DELFIN_FFFREE_"):
            monkeypatch.delenv(k, raising=False)
    yield monkeypatch


def _prep(smiles):
    from delfin.smiles_converter import _prepare_mol_for_embedding
    mol = _prepare_mol_for_embedding(smiles)
    if mol is None:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
    import delfin._bond_decollapse as bd
    metals = [a.GetIdx() for a in mol.GetAtoms() if bd._is_metal(a.GetSymbol())]
    return mol, metals[0]


# --------------------------------------------------------------------------- #
# graph-only detection
# --------------------------------------------------------------------------- #
def test_detect_flags_nhc_carbene(clean_env):
    mol, m = _prep(FUFYOW)
    matom = mol.GetAtomWithIdx(m)
    cands = DEC._carbene_carbon_idxs(mol, m, matom)
    assert len(cands) == 1
    a = mol.GetAtomWithIdx(cands[0])
    assert a.GetAtomicNum() == 6
    # carbene-C sits between two ring N (N-C-N)
    n_ring_het = sum(1 for nb in a.GetNeighbors()
                     if nb.GetAtomicNum() in (7, 8, 15, 16) and nb.IsInRing())
    assert n_ring_het >= 2


def test_detect_ignores_sigma_carbon(clean_env):
    # PtMe2Cl2: the two metal-bound carbons are methyls (no ring, no heteroatoms).
    mol, m = _prep(PTME2CL2)
    matom = mol.GetAtomWithIdx(m)
    assert DEC._carbene_carbon_idxs(mol, m, matom) == []


# --------------------------------------------------------------------------- #
# repaired representation = neutral divalent singlet carbene, sanitizes
# --------------------------------------------------------------------------- #
def test_repair_yields_neutral_divalent_carbene(clean_env):
    mol, m = _prep(FUFYOW)
    matom = mol.GetAtomWithIdx(m)
    cands = DEC._carbene_carbon_idxs(mol, m, matom)
    em = Chem.RWMol(mol)
    for d in [n.GetIdx() for n in matom.GetNeighbors()]:
        if em.GetBondBetweenAtoms(m, d) is not None:
            em.RemoveBond(m, d)
    # without repair the split raises a KekulizeException
    with pytest.raises(Exception):
        Chem.GetMolFrags(em, asMols=True, sanitizeFrags=True)
    # repair makes it sanitize
    DEC._repair_carbene_carbons(em, cands)
    frags = Chem.GetMolFrags(em, asMols=True, sanitizeFrags=True,
                             fragsMolAtomMapping=[])
    assert len(frags) >= 2
    # the carbene-C is now neutral divalent (free singlet carbene :C:)
    a = em.GetAtomWithIdx(cands[0])
    assert a.GetFormalCharge() == 0
    assert a.GetNumRadicalElectrons() == 2
    assert a.GetNoImplicit() is True


# --------------------------------------------------------------------------- #
# native decomposition + build with the flag on
# --------------------------------------------------------------------------- #
def test_nhc_decomposes_natively_with_flag(clean_env):
    for k, v in _BUILD_FLAGS.items():
        clean_env.setenv(k, v)
    d = DEC.decompose(FUFYOW)
    assert d is not None, "FUFYOW must decompose natively with NHC repair on"
    # carbene-C is classified as a normal C sigma-donor, denticity 1
    assert "C" in d["donor_elems"]
    assert all(lg["denticity"] == 1 for lg in d["ligands"])


def test_nhc_builds_isomers_with_flag(clean_env):
    for k, v in _BUILD_FLAGS.items():
        clean_env.setenv(k, v)
    iso = _fffree_isomers(FUFYOW)
    assert iso is not None and len(iso) >= 1
    # every emitted frame is a non-empty, finite XYZ
    for xyz, _label in iso:
        assert "Ag" in xyz
        for tok in xyz.split():
            try:
                v = float(tok)
                assert v == v and abs(v) != float("inf")
            except ValueError:
                pass


# --------------------------------------------------------------------------- #
# env-gate: flag unset/0 -> byte-identical
# --------------------------------------------------------------------------- #
def test_nhc_legacy_when_flag_unset(clean_env):
    # No DELFIN_FFFREE_* set at all -> NHC must stay legacy (None).
    assert _fffree_isomers(FUFYOW) is None
    assert _fffree_isomers(ACIPAG) is None


def test_flag_off_byte_identical_on_nhc(clean_env):
    clean_env.setenv("PYTHONHASHSEED", "0")
    # unset vs explicit 0 must give the identical (None) result on the NHC path.
    unset = _fffree_isomers(FUFYOW)
    clean_env.setenv("DELFIN_FFFREE_NHC_CARBENE", "0")
    explicit0 = _fffree_isomers(FUFYOW)
    assert unset == explicit0 is None


def test_flag_on_inert_on_non_nhc(clean_env):
    clean_env.setenv("PYTHONHASHSEED", "0")
    # baseline (flag off) on a normal Werner complex
    base = _fffree_isomers(CISPLATIN)
    base_co = _fffree_isomers(COCL3NH3)
    assert base is not None and base_co is not None
    # flag on must not change a non-NHC build (repair branch never entered)
    clean_env.setenv("DELFIN_FFFREE_NHC_CARBENE", "1")
    on = _fffree_isomers(CISPLATIN)
    on_co = _fffree_isomers(COCL3NH3)
    assert [x[0] for x in on] == [x[0] for x in base]
    assert [x[0] for x in on_co] == [x[0] for x in base_co]


# --------------------------------------------------------------------------- #
# determinism
# --------------------------------------------------------------------------- #
def test_nhc_build_is_deterministic(clean_env):
    for k, v in _BUILD_FLAGS.items():
        clean_env.setenv(k, v)
    a = _fffree_isomers(FUFYOW)
    b = _fffree_isomers(FUFYOW)
    assert a is not None
    assert [x[0] for x in a] == [x[0] for x in b]
