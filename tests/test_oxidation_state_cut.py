"""The oxidation state from the SMILES as written: cut every metal-donor bond (#528).

The SMILES DELFIN builds from write every metal-donor contact as a bond and mark an L-type donor
with a formal charge: a carbonyl O on the metal is [O+]=C, a pyridine N [N+], an ammine [NH3+],
the metal carries the balance ([Cu-2]).  The donor table of _oxidation_state reads [O+]=C as an
anion and gives Cu(acac)2 as d7; the cut (DELFIN_OXSTATE_CUT=1) gives d9.  Every expected value
below is fixed by the chemistry, not by an earlier run.
"""
import pytest

Chem = pytest.importorskip("rdkit.Chem")

from delfin.manta import _oxidation_state as OX  # noqa: E402

CASES = [
    # SMILES (DELFIN encoding), expected (OS, d) or None, note
    ("CC1=CC(C)=[O+][Cu-2]2([O]1)[O+]=C(C)C=C(C)[O]2", (2, 9), "Cu(acac)2"),
    ("[Cl][Pt-2]([Cl])([NH3+])[NH3+]", (2, 8), "cisplatin"),
    ("[Cl][Cd-2]([Cl])([N+]1=CC=CC=C1)[N+]1=CC=CC=C1", (2, 10), "CdCl2(py)2"),
    ("[O+]#[C][Cr-6]([C]#[O+])([C]#[O+])([C]#[O+])([C]#[O+])[C]#[O+]", (0, 6), "Cr(CO)6"),
    ("C1C[NH+]2CCC[NH+]3CC[NH+]4CCC[NH+](C1)[Ni-2]234", (2, 8), "[Ni(cyclam)]2+"),
    ("CC1=[O+][Zn-2]2([O]1)[O+]=C(C)[O]2", (2, 10), "Zn(OAc)2, chelating carboxylate"),
    ("CC(=O)[O][Cu-2]([OH2+])([OH2+])[O]C(C)=O", (2, 9), "Cu(OAc)2(H2O)2, monodentate"),
    ("CC(C)=[O+][Zn-2]([Cl])([Cl])[O+]=C(C)C", (2, 10), "ZnCl2(acetone)2, ketone O"),
    ("[Cl][Zn-2]([Cl])([OH+]C)[OH+]C", (2, 10), "ZnCl2(MeOH)2, alcohol"),
    ("CC#[N+][Cu-3]([N+]#CC)([N+]#CC)[N+]#CC", (1, 10), "[Cu(MeCN)4]+, nitrile"),
    ("[NH3+][Cu-2]([NH3+])([NH3+])[NH3+].[Cl-].[Cl-]", (2, 9), "[Cu(NH3)4]Cl2, counter-ions"),
    ("O=[Re]([Cl])([Cl])[Cl]", (5, 2), "Re(V) oxo"),
    ("[Cl][Fe-2]([Cl])([Cl])[N]=O", None, "nitrosyl: NO+ or NO- not in the SMILES"),
    ("[Cl][Pt-]1([Cl])[CH2][CH2]1", None, "hapto: adjacent donors"),
]


@pytest.mark.parametrize("smi,expected,note", CASES, ids=[c[2] for c in CASES])
@pytest.mark.parametrize("sanitize", [True, False])
def test_the_cut_reads_the_smiles_as_written(monkeypatch, smi, expected, note, sanitize):
    monkeypatch.setenv(OX.OXSTATE_CUT_FLAG, "1")
    mol = Chem.MolFromSmiles(smi, sanitize=sanitize)
    if mol is None:
        pytest.skip("RDKit cannot sanitise %s" % note)
    if not sanitize:
        mol.UpdatePropertyCache(strict=False)
    mi = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in OX._OX_GROUP)
    assert OX.oxidation_state(mol, mi) == expected, note


def test_switch_off_keeps_the_donor_table(monkeypatch):
    """Off is the old rule, byte for byte -- including its acac error (Cu(IV) d7)."""
    monkeypatch.delenv(OX.OXSTATE_CUT_FLAG, raising=False)
    mol = Chem.MolFromSmiles(CASES[0][0])
    mi = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "Cu")
    assert OX.oxidation_state(mol, mi) == (4, 7)
    monkeypatch.setenv(OX.OXSTATE_CUT_FLAG, "0")
    assert OX.oxidation_state(mol, mi) == (4, 7)


def test_several_metals_stay_undetermined(monkeypatch):
    monkeypatch.setenv(OX.OXSTATE_CUT_FLAG, "1")
    mol = Chem.MolFromSmiles("[Cl][Cu-]([Cl])[Cl+][Cu-]([Cl])[Cl]", sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    assert OX.oxidation_state(mol, 1) is None
