from pathlib import Path

import pytest

from delfin.tools import run_step
from delfin.tools.bio import read_structure
from delfin.tools.bio.processing import load_ligands


@pytest.fixture
def inputs(tmp_path):
    fixtures = Path(__file__).resolve().parents[1] / 'examples' / 'biomolecular_foundation'
    receptor = tmp_path / 'receptor.pdb'
    receptor.write_bytes((fixtures / 'receptor.pdb').read_bytes())
    poses = tmp_path / 'poses.sdf'
    poses.write_bytes((fixtures / 'poses.sdf').read_bytes())
    return dict(receptor_file=receptor, poses_file=poses, record_id=2, shared_coordinate_frame=True)


def test_assembly_preserves_frame_and_records_all_atoms(inputs, tmp_path):
    result = run_step('bio_assemble_complex', **inputs, work_dir=tmp_path / 'run')
    assert result.ok, result.error
    atoms = read_structure(result.artifacts['complex_structure']).first_model.atoms
    original = read_structure(inputs['receptor_file']).first_model.atoms
    assert [a.position for a in atoms[:len(original)]] == [a.position for a in original]
    ligand = load_ligands(inputs['poses_file'])[1]
    assert len(atoms) == len(original) + ligand.GetNumAtoms()
    assert [a.serial for a in atoms] == list(range(1, len(atoms) + 1))
    for actual, expected in zip(atoms[len(original):], ligand.GetConformer().GetPositions()):
        assert actual.position == pytest.approx(expected, abs=0.00051)
    assert result.data['selection']['source_record'] == 2
    assert len(result.data['ligand_mapping']) == ligand.GetNumAtoms()
    assert result.artifacts['ligand'].is_file()


@pytest.mark.parametrize('changes', [{'record_id': 0}, {'record_id': 4}, {'shared_coordinate_frame': False}, {'ligand_chain': 'A'}])
def test_assembly_refuses_ambiguous_or_invalid_selection(inputs, tmp_path, changes):
    if changes.get('record_id') == 0 or changes.get('shared_coordinate_frame') is False:
        with pytest.raises(ValueError):
            run_step('bio_assemble_complex', **{**inputs, **changes}, work_dir=tmp_path / 'bad')
        return
    result = run_step('bio_assemble_complex', **{**inputs, **changes}, work_dir=tmp_path / 'bad')
    assert not result.ok


def test_receptor_serial_references_are_remapped(inputs, tmp_path):
    path = inputs['receptor_file']
    lines = path.read_text().splitlines()
    lines = [line[:6] + f'{int(line[6:11]) + 10:5d}' + line[11:] if line.startswith('ATOM') else line for line in lines]
    path.write_text('\n'.join(lines[:-1] + ['CONECT   11   12', 'END']) + '\n')
    result = run_step('bio_assemble_complex', **inputs, work_dir=tmp_path / 'run')
    assert result.ok, result.error
    assert 'CONECT    1    2' in result.artifacts['complex_structure'].read_text()
    assert result.data['receptor_mapping'][0]['source_serial'] == 11


def test_unsupported_serial_metadata_is_not_silently_dropped(inputs, tmp_path):
    path = inputs['receptor_file']
    path.write_text(path.read_text().replace('END', 'ANISOU    1\nEND'))
    result = run_step('bio_assemble_complex', **inputs, work_dir=tmp_path / 'bad')
    assert not result.ok and 'anisotropic' in result.error
