import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest

from delfin.tools import list_steps, run_step
from delfin.tools.bio import read_structure
from delfin.tools.bio.artifacts import validate_artifact


@pytest.fixture
def inputs(tmp_path):
    list_steps()
    source = Path(__file__).resolve().parents[1] / 'examples/biomolecular_foundation/receptor.pdb'
    lines = source.read_text().splitlines()[:-1]
    lines.append('ATOM      6  H   MET A   1      10.104  13.207  -2.004  1.00  0.00           H')
    path = tmp_path / 'protein.pdb'
    path.write_text('\n'.join(lines + ['END']) + '\n')
    return dict(structure_file=path, expected_version='fixture', forcefield='amber14/protein.ff14SB.xml',
                environment='vacuum', constraints='none', expected_charge=0)


def install_mock(monkeypatch, mode='success'):
    import delfin.tools.adapters.bio_openmm_system as module
    def worker(script, options, versions, **kwargs):
        compile(script, 'system_worker.py', 'exec')
        root = kwargs['work_dir']
        structure = read_structure(root / 'input.pdb')
        pairs = [(r, a) for r in structure.first_model.residues for a in r.atoms]
        atoms = [dict(index=i, name=a.name, element=a.element, chain=r.chain_id, residue_id=str(r.seq_id),
                      residue_name=r.name, insertion_code=r.insertion_code) for i, (r, a) in enumerate(pairs)]
        positions = dict(unit='nanometer', values=[[v / 10 for v in a.position] for r, a in pairs])
        if mode == 'moved':
            positions['values'][0][0] += 1
        if mode == 'wrong_units':
            positions['unit'] = 'angstrom'
        if mode == 'changed_atom':
            atoms[0]['name'] = 'WRONG'
        system = ET.Element('System')
        particles = ET.SubElement(system, 'Particles')
        for i in range(len(atoms) - (mode == 'lost_particle')):
            ET.SubElement(particles, 'Particle', mass='nan' if mode == 'bad_mass' else '1')
        forces = ET.SubElement(system, 'Forces')
        nonbonded = ET.SubElement(forces, 'Force', type='NonbondedForce', method='4' if mode == 'periodic' else '0')
        parameters = ET.SubElement(nonbonded, 'Particles')
        for atom in atoms:
            ET.SubElement(parameters, 'Particle', q='1' if mode == 'wrong_charge' else '0', sig='1', eps='0')
        if options['constraints'] == 'HBonds' and mode != 'lost_constraint':
            constraints = ET.SubElement(system, 'Constraints')
            ET.SubElement(constraints, 'Constraint', p1='0', p2='3', d='0.1')
        ET.ElementTree(system).write(root / 'system.xml')
        (root / 'forcefield.xml').write_text('<ForceField/>')
        bonds = [[0, 3]] if mode != 'bad_bond' else [[0, 999]]
        (root / 'topology.json').write_text(json.dumps(dict(atoms=atoms, bonds=bonds)))
        if mode != 'missing_positions':
            (root / 'positions.json').write_text(json.dumps(positions))
        return {'distributions': versions}
    monkeypatch.setattr(module, 'run_worker', worker)


@pytest.mark.parametrize('mode', ['success', 'moved', 'wrong_units', 'changed_atom', 'lost_particle', 'bad_mass', 'periodic', 'wrong_charge', 'bad_bond', 'missing_positions'])
def test_system_bundle_checks_engine_outputs(inputs, tmp_path, monkeypatch, mode):
    install_mock(monkeypatch, mode)
    result = run_step('bio_openmm_parameterize', **inputs, work_dir=tmp_path / 'run')
    assert result.ok == (mode == 'success'), result.error
    if result.ok:
        details = validate_artifact(result.artifacts['simulation_system'], 'simulation_system')
        assert details['particles'] == 6 and details['charge_e'] == 0


def test_nested_system_changes_invalidate_bundle(inputs, tmp_path, monkeypatch):
    install_mock(monkeypatch)
    result = run_step('bio_openmm_parameterize', **inputs, work_dir=tmp_path / 'run')
    assert result.ok, result.error
    path = result.artifacts['simulation_system']
    (path.parent / 'forcefield.xml').write_text('changed')
    with pytest.raises(ValueError, match='checksum'):
        validate_artifact(path, 'simulation_system')


def test_unsupported_complex_refused_before_worker(inputs, tmp_path, monkeypatch):
    import delfin.tools.adapters.bio_openmm_system as module
    path = inputs['structure_file']
    path.write_text(path.read_text().replace('MET', 'LIG'))
    monkeypatch.setattr(module, 'run_worker', lambda *a, **kw: pytest.fail('must not run'))
    result = run_step('bio_openmm_parameterize', **inputs, work_dir=tmp_path / 'run')
    assert not result.ok and 'canonical protein' in result.error


def test_simulation_bundle_cannot_escape_its_directory(inputs, tmp_path, monkeypatch):
    install_mock(monkeypatch)
    result = run_step('bio_openmm_parameterize', **inputs, work_dir=tmp_path / 'run')
    path = result.artifacts['simulation_system']
    value = json.loads(path.read_text())
    value['files']['system']['path'] = '../escape.xml'
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match='escapes'):
        validate_artifact(path, 'simulation_system')


@pytest.mark.parametrize('mode', ['success', 'lost_constraint'])
def test_hydrogen_constraints_must_match_topology(inputs, tmp_path, monkeypatch, mode):
    install_mock(monkeypatch, mode)
    inputs['constraints'] = 'HBonds'
    result = run_step('bio_openmm_parameterize', **inputs, work_dir=tmp_path / 'run')
    assert result.ok == (mode == 'success'), result.error
