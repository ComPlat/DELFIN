import json
from pathlib import Path

import pytest

from delfin.tools import run_step


@pytest.fixture
def inputs():
    from delfin.tools import list_steps
    list_steps()
    source = Path(__file__).resolve().parents[1] / 'examples' / 'biomolecular_foundation' / 'receptor.pdb'
    return dict(structure_file=source, pdbfixer_version='fixture', openmm_version='fixture',
                ph=7.4, seed=7, waters='drop', missing_atoms='reject')


@pytest.mark.parametrize('mode', ['success', 'moved', 'lost', 'no_hydrogens', 'new_residue', 'heavy_added'])
def test_preparation_verifies_engine_outputs_and_maps_atoms(inputs, tmp_path, monkeypatch, mode):
    import delfin.tools.adapters.bio_receptor as module
    def fake_worker(script, options, versions, **kwargs):
        compile(script, 'worker.py', 'exec')
        assert versions == {'pdbfixer': 'fixture', 'openmm': 'fixture'}
        assert options['missing_atoms'] == 'reject'
        root = kwargs['work_dir']
        lines = (root / 'input.pdb').read_text().splitlines()[:-1]
        if mode == 'moved':
            lines = [line.replace('11.104', '21.104') for line in lines]
        if mode == 'lost':
            lines = [line for line in lines if not line.startswith('ATOM      2')]
        if mode != 'no_hydrogens':
            lines.append('ATOM      6  H   MET A   1      10.104  13.207  -2.004  1.00  0.00           H')
        if mode == 'new_residue':
            lines.append('ATOM      7  H   ALA Z  99      10.104  13.207  -2.004  1.00  0.00           H')
        if mode == 'heavy_added':
            lines.append('ATOM      7  CB  MET A   1      10.104  13.207  -2.004  1.00  0.00           C')
        (root / 'receptor.pdb').write_text('\n'.join(lines + ['END']) + '\n')
        (root / 'worker_report.json').write_text(json.dumps({'ph': 7.4, 'platform': 'Reference'}))
        (root / 'engine.log').write_text('mocked preparation')
        return {'distributions': versions}
    monkeypatch.setattr(module, 'run_worker', fake_worker)
    result = run_step('bio_receptor_prepare', **inputs, work_dir=tmp_path / 'run')
    assert result.ok == (mode == 'success'), result.error
    if result.ok:
        assert len(result.data['atom_mapping']) == 5
        assert len(result.data['added_atoms']) == 1


def test_unknown_chemistry_is_refused_before_launch(inputs, tmp_path, monkeypatch):
    import delfin.tools.adapters.bio_receptor as module
    source = tmp_path / 'modified.pdb'
    source.write_text(inputs['structure_file'].read_text().replace('MET', 'MSE'))
    inputs['structure_file'] = source
    monkeypatch.setattr(module, 'run_worker', lambda *a, **k: pytest.fail('must not launch'))
    result = run_step('bio_receptor_prepare', **inputs, work_dir=tmp_path / 'run')
    assert not result.ok
    assert not list((tmp_path / 'run').iterdir())


def test_missing_protonation_policy_is_a_validation_error(inputs, tmp_path):
    del inputs['ph']
    with pytest.raises(ValueError, match='ph'):
        run_step('bio_receptor_prepare', **inputs, work_dir=tmp_path / 'run')


def test_isolated_worker_uses_versioned_environment_and_reports_errors(tmp_path):
    import importlib.metadata
    import sys
    from delfin.tools.bio.optional_runtime import run_worker
    versions = {'rdkit': importlib.metadata.version('rdkit')}
    good = tmp_path / 'good'
    good.mkdir()
    script = "import json,sys; from pathlib import Path; Path('output.json').write_text(json.dumps(json.load(open(sys.argv[1]))))"
    result = run_worker(script, {'value': 7}, versions, python=sys.executable, work_dir=good, timeout=10)
    assert json.loads((good / 'output.json').read_text()) == {'value': 7}
    assert result['distributions']['rdkit']['version'] == versions['rdkit']
    bad = tmp_path / 'bad'
    bad.mkdir()
    with pytest.raises(RuntimeError, match='worker exited'):
        run_worker("raise ValueError('fixture failure')", {}, versions, python=sys.executable, work_dir=bad, timeout=10)
    assert 'fixture failure' in (bad / 'engine.log').read_text()
