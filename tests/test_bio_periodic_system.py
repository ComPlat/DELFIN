import json
import xml.etree.ElementTree as ET

import pytest

from delfin.tools import run_step
from delfin.tools.bio.artifacts import read_json, validate_artifact
from test_bio_simulation_system import inputs, install_mock


def solvation():
    return dict(box_size_nm=[4., 4., 4.], ionic_strength_m=0., neutralize=True,
                positive_ion='Na+', negative_ion='Cl-', seed=7, solute_charge=0,
                cutoff_nm=0.9, ewald_tolerance=0.0005)


def mock_solvation(monkeypatch, mode='success'):
    import delfin.tools.adapters.bio_openmm_system as module
    install_mock(monkeypatch)
    base = module.run_worker
    def worker(script, options, versions, **kwargs):
        runtime = base(script, options, versions, **kwargs)
        root = kwargs['work_dir']
        lines = (root / 'input.pdb').read_text().splitlines()[:-1]
        if mode == 'moved_protein':
            lines = [line.replace('11.104', '12.104') for line in lines]
        for index, (name, element, x, y, z) in enumerate([('O', 'O', 20., 20., 20.), ('H1', 'H', 20.9, 20., 20.), ('H2', 'H', 20., 20.9, 20.)], 7):
            lines.append(f'HETATM{index:5d} {name:^4s} HOH W   1    {x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{0.:6.2f}          {element:>2s}')
        (root / 'solvated.pdb').write_text('\n'.join(lines + ['END']) + '\n')
        topology = read_json(root / 'topology.json')
        positions = read_json(root / 'positions.json')
        for index, (name, element) in enumerate([('O', 'O'), ('H1', 'H'), ('H2', 'H')], 6):
            topology['atoms'].append(dict(index=index, name=name, element=element, chain='W', residue_id='1', residue_name='HOH', insertion_code=''))
        positions['values'] += [[2., 2., 2.], [2.09, 2., 2.], [2., 2.09, 2.]]
        if mode == 'moved_protein':
            positions['values'][0][0] += .1
        topology['bonds'] += [[6, 7], [6, 8]]
        (root / 'topology.json').write_text(json.dumps(topology))
        (root / 'positions.json').write_text(json.dumps(positions))
        xml = ET.parse(root / 'system.xml').getroot()
        vectors = ET.SubElement(xml, 'PeriodicBoxVectors')
        for i, name in enumerate('ABC'):
            ET.SubElement(vectors, name, **{key: str(4. if i == j else 0.) for j, key in enumerate('xyz')})
        if mode == 'wrong_box':
            vectors.find('A').set('x', '5')
        force = xml.find('./Forces/Force')
        force.set('method', '4'); force.set('cutoff', '1.2' if mode == 'wrong_cutoff' else '0.9'); force.set('ewaldTolerance', '0.0005')
        for mass, charge in [(16., -.834), (1., .417), (1., .417)]:
            ET.SubElement(xml.find('Particles'), 'Particle', mass=str(mass))
            ET.SubElement(force.find('Particles'), 'Particle', q=str(charge), sig='1', eps='0')
        if mode == 'charged_water':
            force.find('Particles')[-1].set('q', '0.5')
        ET.ElementTree(xml).write(root / 'system.xml')
        (root / 'solvent_forcefield.xml').write_text('<ForceField/>')
        return runtime
    monkeypatch.setattr(module, 'run_worker', worker)


@pytest.mark.parametrize('mode', ['success', 'moved_protein', 'wrong_box', 'wrong_cutoff', 'charged_water'])
def test_explicit_solvation_preserves_source_and_validates_periodic_system(inputs, tmp_path, monkeypatch, mode):
    mock_solvation(monkeypatch, mode)
    inputs.update(environment='tip3p_pme', solvation=solvation())
    result = run_step('bio_openmm_parameterize', **inputs, work_dir=tmp_path / 'solvate')
    assert result.ok == (mode == 'success'), result.error
    if result.ok:
        details = validate_artifact(result.artifacts['simulation_system'], 'simulation_system')
        assert details['periodic'] is True and details['particles'] == 9


@pytest.mark.parametrize('change', [{'cutoff_nm': 2.1}, {'seed': 0}, {'positive_ion': 'Ca2+'}, {'box_size_nm': [4, 4]}])
def test_invalid_solvation_policy_fails_before_launch(inputs, tmp_path, monkeypatch, change):
    import delfin.tools.adapters.bio_openmm_system as module
    monkeypatch.setattr(module, 'run_worker', lambda *a, **kw: pytest.fail('must not run'))
    inputs.update(environment='tip3p_pme', solvation={**solvation(), **change})
    with pytest.raises(ValueError):
        run_step('bio_openmm_parameterize', **inputs, work_dir=tmp_path / 'bad')


@pytest.mark.parametrize('mode', ['success', 'bad_box'])
def test_periodic_dynamics_requires_consistent_boxes(inputs, tmp_path, monkeypatch, mode):
    from test_bio_openmm_dynamics import mock_dynamics, settings
    mock_solvation(monkeypatch)
    inputs.update(environment='tip3p_pme', solvation=solvation())
    prepared = run_step('bio_openmm_parameterize', **inputs, work_dir=tmp_path / 'prepare')
    assert prepared.ok, prepared.error
    mock_dynamics(monkeypatch, mode)
    result = run_step('bio_openmm_run', **settings(prepared.artifacts['simulation_system']), work_dir=tmp_path / 'run')
    assert result.ok == (mode == 'success'), result.error
    if result.ok:
        assert validate_artifact(result.artifacts['trajectory'], 'trajectory')['periodic']
        assert validate_artifact(result.artifacts['checkpoint'], 'checkpoint')['step'] == 5


def test_periodic_analysis_needs_explicit_coordinate_policy(inputs, tmp_path, monkeypatch):
    from test_bio_openmm_dynamics import mock_dynamics, settings
    import delfin.tools.adapters.bio_trajectory_analysis as module
    mock_solvation(monkeypatch)
    inputs.update(environment='tip3p_pme', solvation=solvation())
    prepared = run_step('bio_openmm_parameterize', **inputs, work_dir=tmp_path / 'prepare')
    mock_dynamics(monkeypatch)
    dynamics = run_step('bio_openmm_run', **settings(prepared.artifacts['simulation_system']), work_dir=tmp_path / 'run')
    assert dynamics.ok, dynamics.error
    monkeypatch.setattr(module, 'run_worker', lambda *a, **kw: pytest.fail('must not launch'))
    result = run_step('bio_mdanalysis', system_file=prepared.artifacts['simulation_system'], trajectory_file=dynamics.artifacts['trajectory'],
                      expected_version='fixture', atom_indices=[0, 1, 2], alignment='fit', weighting='mass', discard_frames=0, stride=1,
                      work_dir=tmp_path / 'analysis')
    assert not result.ok and 'periodic_policy' in result.error
