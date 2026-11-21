import pytest
import json
from rdkit import Chem

from delfin.tools.bio import prepare_ligand, write_sdf


@pytest.mark.parametrize("smiles", ["CC[NH3+]", "CC(=O)[O-]", "N[C@@H](C)C(=O)O", "[13CH3]CO", "F/C=C/F"])
def test_formal_charge_isotopes_and_stereo_survive_preparation_and_file_roundtrip(tmp_path, smiles):
    original = Chem.MolFromSmiles(smiles)
    expected = Chem.MolToSmiles(original, isomericSmiles=True)
    first = tmp_path / "first.sdf"
    write_sdf(prepare_ligand(smiles=smiles), first)
    reread = next(iter(Chem.SDMolSupplier(str(first), removeHs=True)))
    assert reread is not None
    assert Chem.MolToSmiles(reread, isomericSmiles=True) == expected
    second = tmp_path / "second.sdf"
    write_sdf(prepare_ligand(ligand_file=first), second)
    final = next(iter(Chem.SDMolSupplier(str(second), removeHs=True)))
    assert Chem.MolToSmiles(final, isomericSmiles=True) == expected
    assert Chem.GetFormalCharge(final) == Chem.GetFormalCharge(original)


def test_identity_distinguishes_protonation_and_stereochemistry_and_maps_added_hydrogens():
    values = [prepare_ligand(smiles=s) for s in ("CCN", "CC[NH3+]", "N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O")]
    identities = [json.loads(v.molecule.properties["DELFIN_IDENTITY"]) for v in values]
    assert len({i["chemical_state_id"] for i in identities}) == 4
    assert identities[1]["formal_charge"] == 1
    mapping = identities[0]["atom_mapping"]
    assert [row["input_index"] for row in mapping[:3]] == [0, 1, 2]
    assert all(row["origin"] == "added_hydrogen" and row["input_index"] is None for row in mapping[3:])


def test_state_identity_is_independent_of_smiles_atom_order():
    first = json.loads(prepare_ligand(smiles="CCO").molecule.properties["DELFIN_IDENTITY"])
    second = json.loads(prepare_ligand(smiles="OCC").molecule.properties["DELFIN_IDENTITY"])
    assert first["chemical_state_id"] == second["chemical_state_id"]


def test_removed_fragments_do_not_renumber_original_atom_mapping():
    identity = json.loads(prepare_ligand(smiles="[Na+].CCO", keep_largest=True).molecule.properties["DELFIN_IDENTITY"])
    heavy = [a for a in identity["atom_mapping"] if a["element"] != "H"]
    assert [a["input_index"] for a in heavy] == [1, 2, 3]


@pytest.mark.parametrize("first,second", [("CCN", "CC[NH3+]"), ("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O")])
def test_pose_comparison_refuses_a_different_charge_or_stereoisomer(first, second):
    from delfin.tools.bio import pose_rmsd, StructureFormatError
    with pytest.raises(StructureFormatError, match="not the same molecule"):
        pose_rmsd(prepare_ligand(smiles=first).molecule, prepare_ligand(smiles=second).molecule)
