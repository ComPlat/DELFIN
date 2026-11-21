"""Real RDKit identity checks around mocked OpenFF/OpenMM engine boundaries."""
import json
from itertools import combinations
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from delfin.tools import run_step
from delfin.tools.bio.artifacts import file_hash, read_json, validate_artifact
from delfin.tools.bio.ligand_parameters import ligand_input
from test_bio_simulation_system import inputs as inputs, install_mock


@pytest.fixture
def ligand(tmp_path):
    from delfin.tools import list_steps
    list_steps()
    mol = Chem.AddHs(Chem.MolFromSmiles('CO'))
    assert AllChem.EmbedMolecule(mol, randomSeed=7) == 0
    path = tmp_path / 'ligand.sdf'
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)
    charges = tmp_path / 'charges.json'
    charges.write_text(json.dumps(dict(schema_version=1, ligand_sha256=file_hash(path), unit='elementary_charge',
                                       method='synthetic test charges', source='test fixture; not a scientific model',
                                       values=[-.2, -.4] + [.15] * 4)))
    offxml = tmp_path / 'forcefield.offxml'
    offxml.write_text('<SMIRNOFF version="0.3"><vdW/><Electrostatics/></SMIRNOFF>')
    return dict(ligand_file=path, charges_file=charges, forcefield_file=offxml,
                expected_version='0.15.0', toolkit_version='0.17.0', openmm_version='8.4.0')


def mock_parameters(monkeypatch, mode='success'):
    import delfin.tools.adapters.bio_ligand_parameterize as module
    def worker(script, options, versions, **kwargs):
        compile(script, 'ligand_worker.py', 'exec')
        root = kwargs['work_dir']
        mol = ligand_input(root / 'ligand.sdf')
        charges = read_json(root / 'charges.json')['values']
        xml = ET.Element('ForceField')
        types = ET.SubElement(xml, 'AtomTypes')
        residue = ET.SubElement(ET.SubElement(xml, 'Residues'), 'Residue', name='fixture')
        nonbonded = ET.SubElement(xml, 'NonbondedForce', coulomb14scale=str(1/1.2), lj14scale='0.5')
        ET.SubElement(nonbonded, 'UseAttributeFromResidue', name='charge')
        for i, atom in enumerate(mol.GetAtoms()):
            ET.SubElement(types, 'Type', name=f't{i}', element=atom.GetSymbol(), mass=str(atom.GetMass()), **{'class': f't{i}'})
            ET.SubElement(residue, 'Atom', name=f'L{i+1:03d}', type=f't{i}', charge=str(charges[i] + (1 if mode == 'wrong_charge' and i == 0 else 0)))
            ET.SubElement(nonbonded, 'Atom', sigma='0.3', epsilon='0.1', **{'class': f't{i}'})
        bonds = ET.SubElement(xml, 'HarmonicBondForce')
        for b in mol.GetBonds():
            a, c = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            ET.SubElement(residue, 'Bond', atomName1=f'L{a+1:03d}', atomName2=f'L{c+1:03d}')
            ET.SubElement(bonds, 'Bond', class1=f't{a}', class2=f't{c}', length='0.15', k='100')
        angles = ET.SubElement(xml, 'HarmonicAngleForce')
        for center in mol.GetAtoms():
            for a, c in combinations(center.GetNeighbors(), 2):
                ET.SubElement(angles, 'Angle', class1=f't{a.GetIdx()}', class2=f't{center.GetIdx()}', class3=f't{c.GetIdx()}', angle='1.9', k='10')
        ET.SubElement(xml, 'PeriodicTorsionForce', ordering='smirnoff')
        if mode == 'missing_bond':
            bonds.remove(bonds[0])
        if mode == 'include':
            ET.SubElement(xml, 'Include', file='other.xml')
        ET.ElementTree(xml).write(root / 'parameters.xml')
        return {'distributions': versions}
    monkeypatch.setattr(module, 'run_worker', worker)


@pytest.mark.parametrize('mode', ['success', 'wrong_charge', 'missing_bond', 'include'])
def test_openff_bundle_checks_engine_outputs(ligand, tmp_path, monkeypatch, mode):
    mock_parameters(monkeypatch, mode)
    result = run_step('bio_ligand_parameterize', **ligand, work_dir=tmp_path / 'parameters')
    assert result.ok == (mode == 'success'), result.error
    if result.ok:
        assert validate_artifact(result.artifacts['ligand_parameters'], 'ligand_parameters')['atoms'] == 6
        (tmp_path / 'parameters/charges.json').write_text('{}')
        with pytest.raises(ValueError, match='changed'):
            validate_artifact(result.artifacts['ligand_parameters'], 'ligand_parameters')


@pytest.mark.parametrize('mode', ['charge_sum', 'wrong_hash', 'zero_charges', 'wrong_units', 'bad_offxml'])
def test_ligand_policy_refused_before_worker(ligand, tmp_path, monkeypatch, mode):
    import delfin.tools.adapters.bio_ligand_parameterize as module
    monkeypatch.setattr(module, 'run_worker', lambda *a, **k: pytest.fail('invalid policy must not launch'))
    charges = read_json(ligand['charges_file'])
    if mode == 'charge_sum': charges['values'][0] += 1
    if mode == 'wrong_hash': charges['ligand_sha256'] = 'wrong'
    if mode == 'zero_charges': charges['values'] = [0] * 6
    if mode == 'wrong_units': charges['unit'] = 'coulomb'
    if mode == 'bad_offxml': ligand['forcefield_file'].write_text('<broken')
    ligand['charges_file'].write_text(json.dumps(charges))
    result = run_step('bio_ligand_parameterize', **ligand, work_dir=tmp_path / 'bad')
    assert not result.ok and not (tmp_path / 'bad/worker.py').exists()


def prepared_complex(inputs, ligand, tmp_path, monkeypatch):
    from delfin.tools.bio.assembly import assemble_complex
    mock_parameters(monkeypatch)
    parameters = run_step('bio_ligand_parameterize', **ligand, work_dir=tmp_path / 'parameters')
    assert parameters.ok, parameters.error
    assembly = tmp_path / 'assembly'
    assembly.mkdir()
    assemble_complex(inputs['structure_file'], ligand['ligand_file'], assembly, record_id=1, shared_coordinate_frame=True)
    return {**inputs, 'structure_file': assembly / 'complex.pdb', 'ligand_parameters_file': parameters.artifacts['ligand_parameters']}


def mock_complex_system(monkeypatch, mode='success'):
    import delfin.tools.adapters.bio_openmm_system as module
    install_mock(monkeypatch)
    base = module.run_worker
    def worker(script, options, versions, **kwargs):
        result = base(script, options, versions, **kwargs)
        root = kwargs['work_dir']
        assert options['ligand_forcefield'] == 'ligand_parameters/parameters.xml'
        mol = ligand_input(root / 'ligand_parameters/ligand.sdf')
        values = read_json(root / 'ligand_parameters/charges.json')['values']
        xml = ET.parse(root / 'system.xml').getroot()
        params = xml.findall('./Forces/Force/Particles/Particle')
        for i, q in enumerate(values, 6):
            params[i].set('q', str(q))
            params[i].set('sig', '0.3')
            params[i].set('eps', '0.1')
            xml.findall('./Particles/Particle')[i].set('mass', str(mol.GetAtomWithIdx(i-6).GetMass()))
        if mode == 'wrong_lj': params[6].set('sig', '0.4')
        if mode == 'wrong_atom_charges':
            params[6].set('q', str(values[0] + .1)); params[7].set('q', str(values[1] - .1))
        ET.ElementTree(xml).write(root / 'system.xml')
        topology = read_json(root / 'topology.json')
        topology['bonds'] += [[b.GetBeginAtomIdx() + 6, b.GetEndAtomIdx() + 6] for b in mol.GetBonds()]
        if mode == 'covalent': topology['bonds'].append([0, 6])
        if mode == 'missing_ligand_bond': topology['bonds'].pop()
        (root / 'topology.json').write_text(json.dumps(topology))
        return result
    monkeypatch.setattr(module, 'run_worker', worker)


@pytest.mark.parametrize('mode', ['success', 'wrong_atom_charges', 'wrong_lj', 'covalent', 'missing_ligand_bond'])
def test_complex_system_uses_authoritative_ligand_graph(inputs, ligand, tmp_path, monkeypatch, mode):
    options = prepared_complex(inputs, ligand, tmp_path, monkeypatch)
    mock_complex_system(monkeypatch, mode)
    result = run_step('bio_openmm_complex_parameterize', **options, work_dir=tmp_path / 'system')
    assert result.ok == (mode == 'success'), result.error
    if result.ok:
        from delfin.tools.bio.dynamics import copy_system, system_identity
        copied = copy_system(result.artifacts['simulation_system'], tmp_path / 'relocated')
        assert system_identity(copied) == system_identity(result.artifacts['simulation_system'])
        assert validate_artifact(copied, 'simulation_system')['ligand']
        (copied.parent / 'ligand_parameters/parameters.xml').write_text('modified')
        with pytest.raises(ValueError, match='changed'):
            validate_artifact(copied, 'simulation_system')


def test_complex_coordinate_mismatch_prevents_worker(inputs, ligand, tmp_path, monkeypatch):
    import delfin.tools.adapters.bio_openmm_system as module
    options = prepared_complex(inputs, ligand, tmp_path, monkeypatch)
    path = options['structure_file']
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.startswith('HETATM'):
            lines[i] = line[:30] + f'{999.:8.3f}' + line[38:]
            break
    path.write_text('\n'.join(lines) + '\n')
    monkeypatch.setattr(module, 'run_worker', lambda *a, **k: pytest.fail('changed ligand must not launch'))
    result = run_step('bio_openmm_complex_parameterize', **options, work_dir=tmp_path / 'bad')
    assert not result.ok and 'coordinates' in result.error


@pytest.mark.parametrize('periodic', [False, True])
def test_four_node_complex_dag_retains_bundle_through_resume(inputs, ligand, tmp_path, monkeypatch, periodic):
    from delfin.tools import list_steps
    from delfin.tools.bio.execution import execute_plan
    from delfin.tools.bio.workflows import workflow_template
    from test_bio_openmm_dynamics import mock_dynamics, settings
    from test_bio_trajectory_analysis import execute_analysis_worker
    from test_bio_periodic_system import mock_solvation, solvation
    options = prepared_complex(inputs, ligand, tmp_path, monkeypatch)
    if periodic:
        mock_solvation(monkeypatch, setup=mock_complex_system)
    else:
        mock_complex_system(monkeypatch)
    mock_dynamics(monkeypatch)
    execute_analysis_worker(monkeypatch)
    adapters = list_steps()
    for name in ('bio_ligand_parameterize', 'bio_openmm_complex_parameterize', 'bio_openmm_run', 'bio_mdanalysis'):
        monkeypatch.setattr(adapters[name], 'runtime_preflight', lambda **kwargs: {'version': 'fixture'})
    plan = workflow_template('protein_ligand_md')
    plan['inputs']['complex']['path'] = str(options['structure_file'])
    plan['inputs']['ligand']['path'] = str(ligand['ligand_file'])
    nodes = {node['id']: node for node in plan['nodes']}
    nodes['ligand_parameters']['settings'] = {k: str(v) if isinstance(v, Path) else v for k, v in ligand.items() if k != 'ligand_file'}
    nodes['parameters']['settings'] = {k: v for k, v in inputs.items() if k != 'structure_file'}
    if periodic:
        nodes['parameters']['settings'].update(environment='tip3p_pme', solvation=solvation())
    nodes['simulate']['settings'] = {k: v for k, v in settings(None).items() if k != 'system_file'}
    nodes['analyze']['settings'] = dict(expected_version='fixture', atom_indices=[6, 7, 8], alignment='fit',
                                       weighting='mass', discard_frames=0, stride=1, periodic_policy='as_stored' if periodic else 'reject')
    output = tmp_path / 'dag'
    result = execute_plan(plan, output)
    assert result['status'] == 'success', result
    assert len(result['steps']) == 4
    assert len(execute_plan(plan, output, resume=True)['steps']) == 4
    (output / 'ligand_parameters_1/charges.json').write_text('{}')
    with pytest.raises(ValueError, match='output changed'):
        execute_plan(plan, output, resume=True)


@pytest.mark.parametrize('complete', [False, True])
def test_plan_builder_requires_complete_ligand_policy_without_running_engines(tmp_path, monkeypatch, complete):
    from examples.biomolecular_foundation.build_md_plan import main
    import delfin.tools.bio.optional_runtime as runtime
    monkeypatch.setattr(runtime, 'run_worker', lambda *a, **k: pytest.fail('plan construction cannot run engines'))
    output = tmp_path / 'plan.json'
    argv = ['--structure', '/planned/complex.pdb', '--output', str(output), '--charge', '0',
            '--steps', '5', '--report-interval', '2', '--temperature-k', '300', '--friction-per-ps', '1',
            '--timestep-ps', '0.001', '--seed', '7', '--constraint-tolerance', '0.00001', '--constraints', 'none',
            '--platform', 'Reference', '--minimize', 'yes', '--atom-indices', '6', '7', '8', '--alignment', 'fit',
            '--weighting', 'mass', '--discard-frames', '0', '--stride', '1', '--ligand', '/planned/ligand.sdf']
    if not complete:
        with pytest.raises(SystemExit):
            main(argv)
        assert not output.exists()
        return
    argv += ['--charges-file', '/planned/charges.json', '--ligand-forcefield', '/planned/chosen.offxml']
    assert main(argv) == 0
    plan = read_json(output)['plan']
    assert len(plan['nodes']) == 4
    assert plan['nodes'][0]['settings']['charges_file'] == '/planned/charges.json'
    assert plan['nodes'][1]['backend'] == 'parameterize_complex'
