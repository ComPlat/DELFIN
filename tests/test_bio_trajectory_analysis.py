import json
import sys
import types

import numpy as np
import pytest

from delfin.tools import run_step
from test_bio_openmm_dynamics import system as system, inputs as inputs, settings, mock_dynamics


@pytest.fixture
def trajectory(system, tmp_path, monkeypatch):
    mock_dynamics(monkeypatch)
    result = run_step('bio_openmm_run', **settings(system), work_dir=tmp_path / 'dynamics')
    assert result.ok, result.error
    return result.artifacts['trajectory']


def execute_analysis_worker(monkeypatch, mode='success'):
    import delfin.tools.adapters.bio_trajectory_analysis as module
    # Only the MDAnalysis fitting boundary is mocked; the streaming statistics run.
    def rotation(a, b, weights=None):
        u, _, vt = np.linalg.svd((a * weights[:, None]).T @ b)
        matrix = vt.T @ u.T
        if np.linalg.det(matrix) < 0:
            vt[-1] *= -1
            matrix = vt.T @ u.T
        return matrix, 0.
    def worker(script, options, versions, **kwargs):
        root = kwargs['work_dir']
        (root / 'options.json').write_text(json.dumps(options))
        with monkeypatch.context() as patch:
            patch.chdir(root)
            patch.setattr(sys, 'argv', ['worker.py', 'options.json'])
            patch.setitem(sys.modules, 'MDAnalysis', types.ModuleType('MDAnalysis'))
            patch.setitem(sys.modules, 'MDAnalysis.analysis', types.ModuleType('MDAnalysis.analysis'))
            patch.setitem(sys.modules, 'MDAnalysis.analysis.align', types.SimpleNamespace(rotation_matrix=rotation))
            exec(compile(script, 'analysis_worker.py', 'exec'), {})
        path = root / 'analysis.json'
        data = json.loads(path.read_text())
        if mode == 'lost_frame':
            data['frames'].pop()
        if mode == 'lost_atom':
            data['atoms'].pop()
        if mode == 'negative_statistic':
            data['frames'][0]['rmsd_nm'] = -1
        path.write_text(json.dumps(data))
        return {'distributions': versions}
    monkeypatch.setattr(module, 'run_worker', worker)


@pytest.mark.parametrize('mode', ['success', 'lost_frame', 'lost_atom', 'negative_statistic'])
def test_analysis_retains_selected_atoms_frames_and_units(system, trajectory, tmp_path, monkeypatch, mode):
    execute_analysis_worker(monkeypatch, mode)
    result = run_step('bio_mdanalysis', system_file=system, trajectory_file=trajectory, expected_version='fixture',
                      atom_indices=[0, 1, 2], alignment='fit', weighting='mass', discard_frames=1, stride=2,
                      work_dir=tmp_path / 'analysis')
    assert result.ok == (mode == 'success'), result.error
    if result.ok:
        assert result.data['frame_count'] == 2
        assert [row['step'] for row in result.data['frames']] == [2, 5]
        assert all(row['rmsd_nm'] < 1e-10 for row in result.data['frames'])
        assert all(row['rmsf_nm'] == 0 for row in result.data['atoms'])
        assert result.data['units']['distance'] == 'nanometer'


def test_no_analysis_for_a_changed_system(system, trajectory, tmp_path, monkeypatch):
    import delfin.tools.adapters.bio_trajectory_analysis as module
    monkeypatch.setattr(module, 'run_worker', lambda *a, **kw: pytest.fail('must not launch'))
    data = json.loads(system.read_text())
    (system.parent / data['files']['forcefield']['path']).write_text('changed')
    result = run_step('bio_mdanalysis', system_file=system, trajectory_file=trajectory, expected_version='fixture',
                      atom_indices=[0, 1, 2], alignment='none', weighting='uniform', discard_frames=0, stride=1,
                      work_dir=tmp_path / 'analysis')
    assert not result.ok and 'checksum' in result.error


def test_three_node_md_dag_executes_and_checks_resume(inputs, tmp_path, monkeypatch):
    from test_bio_simulation_system import install_mock
    from delfin.tools.bio.workflows import workflow_template
    from delfin.tools.bio.execution import execute_plan
    from delfin.tools import list_steps
    install_mock(monkeypatch)
    mock_dynamics(monkeypatch)
    execute_analysis_worker(monkeypatch)
    adapters = list_steps()
    for name in ('bio_openmm_parameterize', 'bio_openmm_run', 'bio_mdanalysis'):
        monkeypatch.setattr(adapters[name], 'runtime_preflight', lambda **kwargs: {'version': 'fixture'})
    plan = workflow_template('md')
    plan['inputs']['complex']['path'] = str(inputs['structure_file'])
    plan['nodes'][0]['settings'] = {k: v for k, v in inputs.items() if k != 'structure_file'}
    plan['nodes'][1]['settings'] = {k: v for k, v in settings(None).items() if k != 'system_file'}
    plan['nodes'][2]['settings'] = dict(expected_version='fixture', atom_indices=[0, 1, 2], alignment='fit',
                                       weighting='mass', discard_frames=0, stride=1)
    output = tmp_path / 'dag'
    first = execute_plan(plan, output)
    assert first['status'] == 'success', first
    assert len(first['steps']) == 3
    resumed = execute_plan(plan, output, resume=True)
    assert len(resumed['steps']) == 3
    path = output / 'simulate_1/frames.jsonl'
    path.write_text('changed')
    with pytest.raises(ValueError, match='output changed'):
        execute_plan(plan, output, resume=True)
