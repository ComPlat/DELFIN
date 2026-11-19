"""Boundary tests: a successful process alone is not a successful calculation."""
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from delfin import qm_runtime
from delfin.tools import run_step
from delfin.tools.bio import BindingBox
from delfin.tools.bio.optional_runtime import probe_distribution


@pytest.fixture
def docking(tmp_path, monkeypatch):
    atom = 'ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\n'
    for name in ('receptor.pdbqt', 'ligand.pdbqt'):
        (tmp_path / name).write_text(atom)
    (tmp_path / 'box.json').write_text(json.dumps(BindingBox(center=(0, 0, 0), size=(10, 10, 10)).to_dict()))
    binary = tmp_path / 'vina'
    binary.write_text('mock engine, never executed')
    monkeypatch.setattr(qm_runtime, 'resolve_tool', lambda name: SimpleNamespace(path=str(binary)))
    options = dict(receptor_file=tmp_path / 'receptor.pdbqt', ligand_file=tmp_path / 'ligand.pdbqt',
                   site_file=tmp_path / 'box.json', expected_version='fixture version 1', seed=7)
    return options, atom


def test_docking_records_actual_command_and_score(docking, tmp_path, monkeypatch):
    options, atom = docking
    calls = []
    def fake_run(name, args, **kw):
        calls.append((args, kw))
        assert kw['allow_install'] is False and kw['track_process'] is True
        if args == ['--version']:
            return subprocess.CompletedProcess(args, 0, 'fixture version 1\n', '')
        (Path(kw['cwd']) / 'poses.pdbqt').write_text('MODEL 1\nREMARK VINA RESULT: -1.5 0 0\n' + atom + 'ENDMDL\n')
        return subprocess.CompletedProcess(args, 0)
    monkeypatch.setattr(qm_runtime, 'run_tool', fake_run)
    result = run_step('vina_dock', work_dir=tmp_path / 'run', cores=2, **options)
    assert result.ok, result.error
    assert result.data['scores'][0]['value'] == -1.5
    argv = calls[-1][0]
    assert argv[argv.index('--cpu') + 1] == '2'
    assert argv[argv.index('--seed') + 1] == '7'
    assert result.data['scientific_validation'] == 'not_performed'


@pytest.mark.parametrize('failure', ['version', 'missing_output', 'incomplete_output', 'timeout', 'exit'])
def test_docking_refuses_unsuccessful_boundaries(docking, tmp_path, monkeypatch, failure):
    options, atom = docking
    def fake_run(name, args, **kw):
        if args == ['--version']:
            return subprocess.CompletedProcess(args, 0, 'wrong' if failure == 'version' else 'fixture version 1', '')
        if failure == 'timeout':
            raise subprocess.TimeoutExpired(args, 1)
        if failure == 'incomplete_output':
            (Path(kw['cwd']) / 'poses.pdbqt').write_text('MODEL 1\nREMARK VINA RESULT: -1 0 0\n' + atom)
        return subprocess.CompletedProcess(args, 2 if failure == 'exit' else 0)
    monkeypatch.setattr(qm_runtime, 'run_tool', fake_run)
    result = run_step('vina_dock', work_dir=tmp_path / 'run', **options)
    assert not result.ok
    if failure == 'timeout':
        assert result.error_kind.value == 'timeout'


def test_no_install_means_no_health_installer(monkeypatch):
    import delfin.qm_health
    monkeypatch.setattr(qm_runtime, 'resolve_tool', lambda name: None)
    monkeypatch.setattr(delfin.qm_health, 'ensure_tool', lambda name: pytest.fail('must not install'))
    with pytest.raises(FileNotFoundError):
        qm_runtime.run_tool('absent-fixture', [], allow_install=False)


def test_tracked_timeout_reaps_process_and_unregisters(tmp_path, monkeypatch):
    import delfin.global_manager
    processes, removed = [], []
    manager = SimpleNamespace(register_subprocess=lambda p, **kw: processes.append(p) or 'token',
                              unregister_subprocess=lambda token: removed.append(token))
    monkeypatch.setattr(delfin.global_manager, 'get_global_manager', lambda: manager)
    with pytest.raises(subprocess.TimeoutExpired):
        qm_runtime.run_tool(sys.executable, ['-I', '-c', 'import time; time.sleep(30)'],
                            cwd=tmp_path, timeout=0.1, track_process=True, allow_install=False)
    assert len(processes) == 1 and processes[0].poll() is not None
    assert removed == ['token']


def test_missing_optional_distribution_is_explicit(tmp_path):
    with pytest.raises(ImportError, match='unavailable'):
        probe_distribution('delfin-nonexistent-test-distribution', python=sys.executable,
                           expected_version='1', work_dir=tmp_path)


@pytest.mark.parametrize('count', [0, 2, 3])
def test_posebusters_requires_complete_report_and_preserves_failed_checks(tmp_path, monkeypatch, count):
    import delfin.tools.adapters.bio_posebusters as module
    fixtures = Path(__file__).resolve().parents[1] / 'examples' / 'biomolecular_foundation'
    def fake_console(distribution, entry, args, **kwargs):
        assert distribution == 'posebusters' and entry == 'bust'
        (kwargs['work_dir'] / 'checks.csv').write_text('file,molecule,position,sanitization\n' +
                                                       ''.join(f'poses.sdf,ligand,{i},False\n' for i in range(count)))
        return {'version': 'fixture'}
    monkeypatch.setattr(module, 'run_console', fake_console)
    result = run_step('posebusters_validate', poses_file=fixtures / 'poses.sdf', reference_file=fixtures / 'reference.sdf',
                      receptor_file=fixtures / 'receptor.pdb', expected_version='fixture', work_dir=tmp_path / 'run')
    assert result.ok == (count == 3), result.error
    if result.ok:
        assert all(row['sanitization'] == 'False' for row in result.data['checks'])


@pytest.mark.parametrize('operation', ['meeko_ligand_prepare', 'meeko_receptor_prepare', 'meeko_pose_export', 'plip_analyze'])
def test_optional_processing_uses_native_entry_and_validates_outputs(tmp_path, monkeypatch, operation):
    from rdkit import Chem
    import delfin.tools.adapters.bio_optional_processing as module
    fixtures = Path(__file__).resolve().parents[1] / 'examples' / 'biomolecular_foundation'
    atom = 'ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\n'
    source = fixtures / 'receptor.pdb'
    if operation == 'meeko_ligand_prepare':
        source = tmp_path / 'explicit_h.sdf'
        mol = Chem.SDMolSupplier(str(fixtures / 'reference.sdf'), removeHs=False)[0]
        with Chem.SDWriter(str(source)) as writer:
            writer.write(Chem.AddHs(mol, addCoords=True))
    elif operation == 'meeko_pose_export':
        source = tmp_path / 'input.pdbqt'
        source.write_text(atom)
    expected = {'meeko_ligand_prepare': ('mk_prepare_ligand.py', 'ligand.pdbqt'),
                'meeko_receptor_prepare': ('mk_prepare_receptor.py', 'receptor.pdbqt'),
                'meeko_pose_export': ('mk_export.py', 'poses.sdf'),
                'plip_analyze': ('plip', 'plip/report.xml')}
    entry, filename = expected[operation]
    def fake_console(distribution, actual_entry, args, **kwargs):
        assert actual_entry == entry
        root = kwargs['work_dir']
        target = root / filename
        target.parent.mkdir(exist_ok=True)
        if operation == 'meeko_pose_export':
            target.write_bytes((fixtures / 'poses.sdf').read_bytes())
        elif operation == 'plip_analyze':
            target.write_text('<report><bindingsite id="1"/></report>')
        else:
            target.write_text(atom)
        (root / 'engine.log').write_text('mock console boundary')
        return {'entry_point': entry, 'version': 'fixture'}
    monkeypatch.setattr(module, 'run_console', fake_console)
    result = run_step(operation, input_file=source, expected_version='fixture', work_dir=tmp_path / 'run')
    assert result.ok, result.error
    assert result.data['runtime']['entry_point'] == entry
    if operation == 'plip_analyze':
        assert result.data['binding_sites'] == 1
