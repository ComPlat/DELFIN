import json
from pathlib import Path

import pytest

from delfin.tools import list_steps, run_step
from delfin.tools.bio.artifacts import file_hash, validate_artifact
from delfin.tools.bio.assets import validate_assets


@pytest.fixture
def inputs(tmp_path):
    list_steps()
    from delfin.tools.adapters.bio_esmfold import _REQUIRED
    root = tmp_path / 'assets'
    root.mkdir()
    files = []
    for relative in _REQUIRED:
        path = root / relative
        path.parent.mkdir(exist_ok=True)
        path.write_text('synthetic weight identity; never loaded by a model')
        files.append({'path': relative, 'sha256': file_hash(path)})
    manifest = root / 'assets.json'
    manifest.write_text(json.dumps({'schema_version': 1, 'engine': 'esmfold', 'files': files}))
    fasta = tmp_path / 'input.fasta'
    fasta.write_text('>synthetic\nME\n')
    return dict(request_file=fasta, assets_file=manifest, expected_version='fixture', torch_version='fixture', device='cpu', seed=7)


def test_asset_changes_or_escape_are_refused(inputs):
    manifest = inputs['assets_file']
    value = json.loads(manifest.read_text())
    validate_assets(manifest, engine='esmfold')
    path = manifest.parent / value['files'][0]['path']
    path.write_text('changed')
    with pytest.raises(ValueError, match='changed'):
        validate_assets(manifest, engine='esmfold')
    value['files'][0]['path'] = '../escape.pt'
    manifest.write_text(json.dumps(value))
    with pytest.raises(ValueError, match='escapes'):
        validate_assets(manifest, engine='esmfold')


@pytest.mark.parametrize('mode', ['success', 'missing_model', 'bad_sequence', 'incomplete_backbone', 'bad_confidence'])
def test_prediction_checks_model_completeness_and_confidence(inputs, tmp_path, monkeypatch, mode):
    import delfin.tools.adapters.bio_esmfold as module
    def worker(script, options, versions, **kwargs):
        compile(script, 'worker.py', 'exec')
        assert options['sequence'] == 'ME'
        assert versions == {'fair-esm': 'fixture', 'torch': 'fixture'}
        root = kwargs['work_dir']
        if mode != 'missing_model':
            rows = []
            for residue, name in enumerate(('MET', 'GLY' if mode == 'bad_sequence' else 'GLU'), 1):
                for atom in ('N', 'CA', 'C', 'O'):
                    if mode == 'incomplete_backbone' and atom == 'O':
                        continue
                    serial = len(rows) + 1
                    rows.append(f'ATOM  {serial:5d} {atom:^4s} {name} A{residue:4d}    {float(serial):8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{80.:6.2f}          {atom[0]:>2s}')
            (root / 'model.pdb').write_text('\n'.join(rows + ['END']) + '\n')
        (root / 'confidence.json').write_text(json.dumps({'mean_plddt': 180 if mode == 'bad_confidence' else 80, 'ptm': 0.7}))
        (root / 'engine.log').write_text('mocked inference boundary')
        return {'distributions': versions}
    monkeypatch.setattr(module, 'run_worker', worker)
    result = run_step('esmfold_predict', **inputs, work_dir=tmp_path / 'run')
    assert result.ok == (mode == 'success'), result.error
    if result.ok:
        validate_artifact(result.artifacts['prediction_bundle'], 'prediction_bundle')
        bundle = json.loads(result.artifacts['prediction_bundle'].read_text())
        assert bundle['provenance_origin'] == 'engine_run'
        assert bundle['models'][0]['confidence']['sha256']


def test_wrong_assets_block_before_worker(inputs, tmp_path, monkeypatch):
    import delfin.tools.adapters.bio_esmfold as module
    monkeypatch.setattr(module, 'run_worker', lambda *a, **k: pytest.fail('must not run'))
    value = json.loads(inputs['assets_file'].read_text())
    value['files'].pop()
    inputs['assets_file'].write_text(json.dumps(value))
    result = run_step('esmfold_predict', **inputs, work_dir=tmp_path / 'run')
    assert not result.ok and 'required model assets' in result.error
