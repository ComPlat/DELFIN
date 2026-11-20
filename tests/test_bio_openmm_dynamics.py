import json
import xml.etree.ElementTree as ET

import pytest

from delfin.tools import run_step
from delfin.tools.bio.artifacts import read_json, validate_artifact
from delfin.tools.bio.dynamics import frame_steps
from test_bio_simulation_system import inputs, install_mock


@pytest.fixture
def system(inputs, tmp_path, monkeypatch):
    install_mock(monkeypatch)
    result = run_step('bio_openmm_parameterize', **inputs, work_dir=tmp_path / 'parameterize')
    assert result.ok, result.error
    return result.artifacts['simulation_system']


def settings(system):
    return dict(system_file=system, expected_version='fixture', platform='Reference', temperature_k=300., friction_per_ps=1.,
                timestep_ps=0.001, seed=7, constraint_tolerance=1e-5, steps=5, report_interval=2, minimize=True)


def mock_dynamics(monkeypatch, mode='success'):
    import delfin.tools.adapters.bio_openmm_run as module
    calls = []
    def worker(script, options, versions, **kwargs):
        compile(script, 'dynamics_worker.py', 'exec')
        calls.append(options)
        root = kwargs['work_dir']; protocol = options['protocol']
        bundle = read_json(root / 'system/simulation_system.json')
        positions = read_json(root / 'system' / bundle['files']['positions']['path'])['values']
        frames = [dict(step=step, time_ps=options['start_time_ps'] + (step-options['start_step'])*protocol['timestep_ps'],
                       potential_kj_mol=-10., kinetic_kj_mol=1., positions_nm=positions)
                  for step in frame_steps(options['start_step'], options['steps'], options['report_interval'])]
        box = None
        if bundle['environment'] == 'tip3p_pme':
            from delfin.tools.bio.dynamics import system_box
            box = system_box(root / 'system/simulation_system.json')
            for frame in frames:
                frame['box_vectors_nm'] = box
        final = frames[-1]
        state = ET.Element('State', time=str(final['time_ps']), stepCount=str(final['step']))
        for group, tag in [('Positions', 'Position'), ('Velocities', 'Velocity')]:
            collection = ET.SubElement(state, group)
            for row in positions:
                ET.SubElement(collection, tag, **dict(zip('xyz', map(str, row))))
        if box is not None:
            vectors = ET.SubElement(state, 'PeriodicBoxVectors')
            for name, row in zip('ABC', box):
                ET.SubElement(vectors, name, **dict(zip('xyz', map(str, row))))
        if mode == 'bad_state':
            state.set('stepCount', '999')
        ET.ElementTree(state).write(root / 'state.xml')
        integrator = ET.Element('Integrator', type='LangevinMiddleIntegrator')
        for attr, key in [('stepSize', 'timestep_ps'), ('temperature', 'temperature_k'), ('friction', 'friction_per_ps'), ('constraintTolerance', 'constraint_tolerance'), ('randomSeed', 'seed')]:
            integrator.set(attr, str(protocol[key]))
        if mode == 'bad_integrator':
            integrator.set('stepSize', '0.9')
        ET.ElementTree(integrator).write(root / 'integrator.xml')
        (root / 'checkpoint.chk').write_bytes(b'synthetic opaque checkpoint; not loadable by OpenMM')
        minimized = dict(performed=options['minimize'])
        if options['minimize']:
            minimized.update(before_kj_mol=-5., after_kj_mol=5. if mode == 'energy_increase' else -10.)
        (root / 'worker_report.json').write_text(json.dumps(dict(final=final, platform=protocol['platform'], minimization=minimized)))
        if mode == 'bad_box':
            frames[1]['box_vectors_nm'] = [[9., 0., 0.], [0., 9., 0.], [0., 0., 9.]]
        if mode == 'missing_frame':
            frames.pop(1)
        if mode == 'bad_time':
            frames[1]['time_ps'] += 0.1
        if mode == 'bad_atoms':
            frames[1]['positions_nm'] = positions[:-1]
        (root / 'frames.jsonl').write_text(''.join(json.dumps(frame) + '\n' for frame in frames))
        return {'distributions': versions}
    monkeypatch.setattr(module, 'run_worker', worker)
    return calls


@pytest.mark.parametrize('mode', ['success', 'missing_frame', 'bad_time', 'bad_atoms', 'bad_state', 'bad_integrator', 'energy_increase'])
def test_dynamics_checks_trajectory_and_restart_outputs(system, tmp_path, monkeypatch, mode):
    mock_dynamics(monkeypatch, mode)
    result = run_step('bio_openmm_run', **settings(system), work_dir=tmp_path / 'dynamics')
    assert result.ok == (mode == 'success'), result.error
    if result.ok:
        assert validate_artifact(result.artifacts['trajectory'], 'trajectory')['frames'] == 4
        assert validate_artifact(result.artifacts['checkpoint'], 'checkpoint')['step'] == 5


def test_restart_keeps_absolute_time_and_refuses_protocol_change(system, tmp_path, monkeypatch):
    calls = mock_dynamics(monkeypatch)
    options = settings(system)
    first = run_step('bio_openmm_run', **options, work_dir=tmp_path / 'first')
    assert first.ok, first.error
    options.update(checkpoint_file=first.artifacts['checkpoint'], minimize=False)
    second = run_step('bio_openmm_run', **options, work_dir=tmp_path / 'second')
    assert second.ok, second.error
    assert second.data['validation']['final_step'] == 10
    assert calls[-1]['start_time_ps'] == 0.005
    options['temperature_k'] = 310.
    failed = run_step('bio_openmm_run', **options, work_dir=tmp_path / 'changed')
    assert not failed.ok and 'protocol changed' in failed.error
    assert len(calls) == 2


def test_minimization_only_emits_one_frame(system, tmp_path, monkeypatch):
    mock_dynamics(monkeypatch)
    options = settings(system); options['steps'] = 0
    result = run_step('bio_openmm_run', **options, work_dir=tmp_path / 'minimum')
    assert result.ok, result.error
    assert result.data['validation']['frames'] == 1


@pytest.mark.parametrize('key,value', [('seed', 0), ('timestep_ps', float('nan')), ('max_frames', 1), ('constraint_tolerance', 1)])
def test_invalid_protocol_never_launches(system, tmp_path, monkeypatch, key, value):
    calls = mock_dynamics(monkeypatch)
    options = settings(system); options[key] = value
    with pytest.raises(ValueError):
        run_step('bio_openmm_run', **options, work_dir=tmp_path / 'invalid')
    assert not calls
