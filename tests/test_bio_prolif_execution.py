import json

import pytest
from rdkit import Chem

from delfin.tools import list_steps, run_step


def molecule(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, (i * 1.2, (i % 3) * 0.7, (i % 2) * 0.3))
    mol.AddConformer(conf)
    return mol


@pytest.mark.parametrize('mode', ['success', 'reordered', 'lost_pose', 'bad_mapping', 'moved', 'changed_chemistry', 'bad_contact', 'lost_receptor'])
def test_prolif_validates_all_outputs(tmp_path, monkeypatch, mode):
    list_steps()
    import delfin.tools.adapters.bio_prolif as module
    protein = molecule('NCC(=O)O')
    for atom in protein.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f'{atom.GetSymbol()}{atom.GetIdx():<3}'[:4])
        info.SetResidueName('GLY')
        info.SetResidueNumber(1)
        info.SetChainId('A')
        info.SetIsHeteroAtom(False)
        atom.SetMonomerInfo(info)
    receptor = tmp_path / 'protein.pdb'
    Chem.MolToPDBFile(protein, str(receptor))
    pose = molecule('CO')
    poses = tmp_path / 'poses.sdf'
    with Chem.SDWriter(str(poses)) as writer:
        writer.write(pose)

    def worker(script, options, versions, **kwargs):
        compile(script, 'prolif_worker.py', 'exec')
        root = kwargs['work_dir']
        output_protein = Chem.RemoveHs(protein) if mode == 'lost_receptor' else protein
        Chem.MolToMolFile(output_protein, str(root / 'prolif_receptor.mol'), forceV3000=True)
        order = list(reversed(range(pose.GetNumAtoms()))) if mode == 'reordered' else list(range(pose.GetNumAtoms()))
        exported = Chem.RenumberAtoms(pose, order)
        if mode == 'moved':
            exported.GetConformer().SetAtomPosition(0, (10, 10, 10))
        if mode == 'changed_chemistry':
            exported.GetAtomWithIdx(1).SetIsotope(18)
        Chem.MolToMolFile(exported, str(root / 'pose.mol'), forceV3000=True)
        mapping = [{'parent_index': i, 'source_atom_index': source} for i, source in enumerate(order)]
        if mode == 'bad_mapping':
            mapping[0]['source_atom_index'] = True
        contacts = [{'interaction': 'Hydrophobic', 'metadata': {'parent_indices': {'ligand': [999 if mode == 'bad_contact' else 0], 'protein': [0]}, 'distance': 3.0}}]
        rows = [] if mode == 'lost_pose' else [{'record_id': 1, 'molecule_file': 'pose.mol', 'atom_mapping': mapping, 'interactions': contacts}]
        (root / 'interactions.json').write_text(json.dumps({'records': rows, 'receptor_molecule': 'prolif_receptor.mol'}))
        return {'distributions': versions}

    monkeypatch.setattr(module, 'run_worker', worker)
    result = run_step('prolif_analyze', poses_file=poses, receptor_file=receptor,
                      interactions=['Hydrophobic'], expected_version='fixture', bond_policy='rdkit_proximity',
                      shared_coordinate_frame=True, work_dir=tmp_path / 'run')
    assert result.ok == (mode in {'success', 'reordered'}), result.error
    if result.ok:
        assert result.data['records'][0]['interactions'][0]['metadata']['distance'] == 3.0


def test_wrong_permutation_of_identical_atoms_is_refused():
    from delfin.tools.adapters.bio_prolif import _validate_mapped_pose
    pose = molecule('CCO')
    mapping = [{'parent_index': i, 'source_atom_index': i} for i in range(pose.GetNumAtoms())]
    mapping[0]['source_atom_index'], mapping[1]['source_atom_index'] = 1, 0
    with pytest.raises(ValueError, match='coordinates|connectivity'):
        _validate_mapped_pose(pose, pose, mapping)
