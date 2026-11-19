import sys
from types import SimpleNamespace

import pytest

from delfin import qm_runtime
from delfin.tools.bio import environments as envs, optional_runtime as runtime


def absent(*args, **kwargs):
    raise ImportError('not available outside module')


def test_module_found_before_install_and_preserved_for_version_probe(tmp_path, monkeypatch):
    attempts = []
    def ensure(*args, **kwargs):
        attempts.append(kwargs['allow_install'])
        return absent()
    monkeypatch.setattr(envs, 'ensure_python', ensure)
    candidates = [SimpleNamespace(path=sys.executable, module_name='bio/prolif/old'),
                  SimpleNamespace(path=sys.executable, module_name='bio/prolif/2.2.3')]
    monkeypatch.setattr(envs, 'module_python_candidates', lambda item: candidates)
    def check(item, python, **kwargs):
        if kwargs['environment_module'].endswith('old'):
            raise ValueError('wrong version')
    monkeypatch.setattr(envs, 'check_python', check)
    calls = []
    monkeypatch.setattr(qm_runtime, 'run_tool', lambda executable, args, **kw: calls.append((args, kw)) or SimpleNamespace(returncode=0, stdout='2.2.3', stderr=''))
    result = runtime.probe_distribution('prolif', python=sys.executable, expected_version='2.2.3', work_dir=tmp_path, allow_install=True)
    assert attempts == [False]
    assert result['environment_module'] == 'bio/prolif/2.2.3'
    assert calls[0][0][:2] == ['-s', '-P']
    assert calls[0][1]['environment_module'] == result['environment_module']
    assert calls[0][1]['allow_install'] is False


def test_readonly_failed_modules_do_not_provision(tmp_path, monkeypatch):
    attempts = []
    monkeypatch.setattr(envs, 'ensure_python', lambda *a, **kw: attempts.append(kw['allow_install']) or absent())
    monkeypatch.setattr(envs, 'module_python_candidates', lambda item: [SimpleNamespace(path=sys.executable, module_name='bio/prolif/broken')])
    monkeypatch.setattr(envs, 'check_python', lambda *a, **kw: (_ for _ in ()).throw(ValueError('missing library')))
    with pytest.raises(ImportError, match='missing library'):
        runtime.probe_distribution('prolif', python=sys.executable, expected_version='2.2.3', work_dir=tmp_path)
    assert attempts == [False, False]


def test_worker_keeps_module_for_every_dependency_and_execution(tmp_path, monkeypatch):
    probes, executions = [], []
    def probe(name, **kwargs):
        probes.append((name, kwargs))
        return {'python': sys.executable, 'environment_module': 'bio/pdbfixer/1.12'}
    monkeypatch.setattr(runtime, 'probe_distribution', probe)
    monkeypatch.setattr(qm_runtime, 'run_tool', lambda executable, args, **kw: executions.append((args, kw)) or SimpleNamespace(returncode=0))
    runtime.run_worker('pass', {}, {'pdbfixer': '1.12', 'openmm': '8.4.0'}, python=sys.executable, work_dir=tmp_path, timeout=10)
    assert probes[1][1]['discover'] is False
    assert probes[1][1]['environment_module'] == 'bio/pdbfixer/1.12'
    assert executions[0][0][:2] == ['-s', '-P']
    assert executions[0][1]['environment_module'] == 'bio/pdbfixer/1.12'
    assert executions[0][1]['track_process'] is True


def test_same_interpreter_in_different_modules_is_not_deduplicated(monkeypatch):
    from delfin import system_tools
    monkeypatch.setattr(envs, '_module_python_on_path', lambda module: None)
    names = ['bio/prolif/2.2.3', 'bio/prolif/2.1.0']
    monkeypatch.setattr(system_tools, 'matching_modules', lambda patterns: names)
    monkeypatch.setattr(system_tools, 'discover_system_tool_candidates',
                        lambda *a, **kw: [SimpleNamespace(path=sys.executable, module_name=kw['module_patterns'][0])])
    assert [c.module_name for c in envs.module_python_candidates(envs.recipe('prolif'))] == names


def test_runtime_wraps_explicit_module_and_quotes_arguments(monkeypatch):
    import shlex
    calls = []
    monkeypatch.setattr(qm_runtime.subprocess, 'run', lambda cmd, **kw: calls.append(cmd) or SimpleNamespace(returncode=0))
    qm_runtime.run_tool(sys.executable, ['-c', 'print("$HOME; literal")'], environment_module='bio/prolif/2.2.3', allow_install=False)
    command = calls[0]
    assert command[:2] == ['bash', '-lc']
    assert command[2].startswith('module load bio/prolif/2.2.3 ')
    assert shlex.split(command[2].split(' && exec ', 1)[1])[-1] == 'print("$HOME; literal")'


def test_explicit_venv_python_retains_its_environment(tmp_path):
    import json
    import venv
    directory = tmp_path / 'isolated'
    venv.EnvBuilder(with_pip=False, symlinks=True).create(directory)
    python = directory / 'bin/python'
    result = qm_runtime.run_tool(str(python), ['-I', '-c', 'import sys,json;print(json.dumps(sys.prefix))'],
                                 capture_output=True, timeout=10, allow_install=False)
    assert result.returncode == 0
    assert json.loads(result.stdout) == str(directory)


def test_module_can_supply_python_through_dependency(monkeypatch):
    calls = []
    monkeypatch.setattr(qm_runtime, 'run_tool', lambda command, args, **kw: calls.append((command, args, kw)) or SimpleNamespace(returncode=0, stdout=sys.executable + '\n'))
    found = envs._module_python_on_path('bio/prolif/2.2.3')
    assert found.path == sys.executable and found.module_name == 'bio/prolif/2.2.3'
    assert calls[0][2]['track_process'] and not calls[0][2]['allow_install']
    assert calls[0][2]['timeout'] == 20
