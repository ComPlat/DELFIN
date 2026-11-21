import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from delfin.tools.bio import environments as envs
from conftest import child_env


@pytest.fixture(autouse=True)
def isolated(tmp_path, monkeypatch):
    monkeypatch.setenv('DELFIN_BIO_TOOLS_ROOT', str(tmp_path / 'bio'))
    envs._TRIED.clear()


def test_central_installer_and_resolver_share_recipes():
    from delfin import installer, qm_runtime
    assert set(installer.names('bio')) == set(envs.RECIPES)
    for item in envs.RECIPES.values():
        for command in item.commands:
            assert installer.find(command).name == item.name
            assert qm_runtime._SPEC_BY_CANONICAL[command].module_patterns


def test_existing_compatible_environment_prevents_install(monkeypatch):
    monkeypatch.setattr(envs, 'python_candidates', lambda item: iter(['/cluster/env/bin/python']))
    monkeypatch.setattr(envs, 'check_python', lambda item, python: {'version': item.version})
    monkeypatch.setattr(envs, 'install', lambda *a, **k: pytest.fail('existing installation must win'))
    assert envs.ensure_python('meeko', '0.7.1', allow_install=True) == '/cluster/env/bin/python'


def test_readonly_check_never_installs(monkeypatch):
    monkeypatch.setattr(envs, 'python_candidates', lambda item: iter([]))
    monkeypatch.setattr(envs, 'install', lambda *a, **k: pytest.fail('read-only check installed'))
    with pytest.raises(ImportError):
        envs.ensure_python('meeko', '0.7.1')


def test_execution_installs_once_and_uses_managed_python(monkeypatch):
    monkeypatch.setattr(envs, 'python_candidates', lambda item: iter([]))
    monkeypatch.setenv('DELFIN_AUTO_INSTALL_QM_TOOLS', '1')
    calls = []
    def install(names, **kwargs):
        calls.append((names, kwargs))
        return {'ok': True, 'results': []}
    monkeypatch.setattr(envs, 'install', install)
    selected = envs.ensure_python('meeko', '0.7.1', allow_install=True)
    assert selected == str(envs.prefix(envs.recipe('meeko', '0.7.1')) / 'bin/python')
    assert len(calls) == 1
    with pytest.raises(ImportError):
        envs.ensure_python('meeko', '0.7.1', allow_install=True)
    assert len(calls) == 1


def test_global_qm_optout_applies_to_bio(monkeypatch):
    monkeypatch.setenv('DELFIN_AUTO_INSTALL_QM_TOOLS', '0')
    monkeypatch.setattr(envs, 'python_candidates', lambda item: iter([]))
    monkeypatch.setattr(envs, 'install', lambda *a, **k: pytest.fail('installation disabled'))
    with pytest.raises(ImportError):
        envs.ensure_python('meeko', '0.7.1', allow_install=True)


def test_commands_are_fixed_and_versions_cannot_inject_solver_options():
    first, second = envs.recipe('meeko', '0.7.1'), envs.recipe('meeko', '0.8.0')
    assert envs.prefix(first) != envs.prefix(second)
    command = envs.installation_command(first, '/bin/micromamba')
    assert 'meeko=0.7.1' in command and '--override-channels' in command
    assert str(Path(sys.prefix)) != command[command.index('--prefix') + 1]
    for invalid in ('>=0.7', '0.7 --channel evil', 'https://example.org/package', ''):
        with pytest.raises(ValueError):
            envs.recipe('meeko', invalid)


def test_failed_install_is_not_reported_ready(tmp_path, monkeypatch):
    from delfin import installer
    monkeypatch.setattr(envs.shutil, 'which', lambda name: '/fake/micromamba')
    monkeypatch.setattr(installer, '_run', lambda *a, **k: (False, ['solver conflict']))
    result = envs.install(['meeko'])
    assert not result['ok']
    assert result['results'][0]['lines'] == ['solver conflict']
    assert Path(result['results'][0]['manifest']).is_file()


def test_installer_routes_bio_to_isolated_family(monkeypatch):
    from delfin import installer
    calls = []
    monkeypatch.setattr(envs, 'install', lambda names, **kw: calls.append(names) or {'ok': True, 'results': [{'lines': []}]})
    monkeypatch.setattr(envs, 'present', lambda name: True)
    result = installer.install(['meeko'])
    assert result['ok'] and calls == [['meeko']]


def test_repair_never_mutates_a_prefix_without_ownership_manifest(monkeypatch):
    from delfin import installer
    item = envs.recipe('meeko')
    envs.prefix(item).mkdir(parents=True)
    monkeypatch.setattr(envs.shutil, 'which', lambda name: '/fake/micromamba')
    monkeypatch.setattr(envs, 'check_python', lambda *a: (_ for _ in ()).throw(ValueError('broken')))
    monkeypatch.setattr(installer, '_run', lambda *a, **k: pytest.fail('unowned prefix must not be modified'))
    result = envs.install(['meeko'], repair=True)
    assert not result['ok'] and 'ownership' in result['results'][0]['lines'][0]


def test_existing_managed_binary_is_reusable_after_install_attempt(monkeypatch):
    item = envs.recipe('vina')
    binary = envs.prefix(item) / 'bin/vina'
    binary.parent.mkdir(parents=True)
    binary.write_text('fixture')
    envs._TRIED.add((item.name, item.version, str(envs.root())))
    monkeypatch.setattr(envs, 'check_binary', lambda *a: {'version': item.version})
    monkeypatch.setattr(envs, 'install', lambda *a, **k: pytest.fail('must reuse'))
    assert envs.ensure_binary('vina', item.version) == str(binary)


def test_installer_timeout_stops_children_holding_output_pipe(tmp_path):
    from delfin import installer
    import time
    start = time.monotonic()
    command = [sys.executable, '-c', 'import subprocess,sys,time; subprocess.Popen([sys.executable,"-c","import time;time.sleep(30)"]);time.sleep(30)']
    ok, lines = installer._run(command, cwd=tmp_path, env=child_env(tmp_path), on_line=None, timeout=0.2)
    assert not ok and time.monotonic() - start < 5
    assert any('stopped after' in line for line in lines)


@pytest.mark.parametrize('content', ['{broken', '[]', 'null'])
def test_corrupt_ownership_refuses_repair_without_crashing(monkeypatch, content):
    from delfin import installer
    item = envs.recipe('meeko')
    target = envs.prefix(item)
    target.mkdir(parents=True)
    manifest = envs.root() / (target.name + '.json')
    manifest.write_text(content)
    monkeypatch.setattr(envs.shutil, 'which', lambda name: '/fake/micromamba')
    monkeypatch.setattr(envs, 'check_python', lambda *a: (_ for _ in ()).throw(ValueError('broken')))
    monkeypatch.setattr(installer, '_run', lambda *a, **k: pytest.fail('must not mutate'))
    result = envs.install(['meeko'], repair=True)
    assert not result['ok']
    assert manifest.read_text() == content


def test_broken_package_provenance_prevents_ready_state(monkeypatch):
    from delfin import installer
    item = envs.recipe('meeko')
    def install(*args, **kwargs):
        metadata = envs.prefix(item) / 'conda-meta'
        metadata.mkdir(parents=True)
        (metadata / 'broken.json').write_text('{')
        return True, []
    monkeypatch.setattr(envs.shutil, 'which', lambda name: '/fake/micromamba')
    monkeypatch.setattr(installer, '_run', install)
    monkeypatch.setattr(envs, 'check_python', lambda *a: {'version': item.version})
    result = envs.install(['meeko'])
    assert not result['ok']
    assert 'Invalid package provenance' in result['results'][0]['lines'][0]


def test_posebusters_uses_shared_isolated_recipe_and_missing_tool_provisioning(monkeypatch):
    from delfin import installer
    item = envs.recipe('posebusters')
    assert installer.find('bust').name == item.name
    command = envs.installation_command(item, '/fake/micromamba')
    assert 'posebusters=0.5.1' in command
    assert 'rdkit>=2024.9' in command
    assert 'posebusters.cli' in item.modules
    monkeypatch.setattr(envs, 'python_candidates', lambda item: iter([]))
    monkeypatch.setenv('DELFIN_AUTO_INSTALL_QM_TOOLS', '1')
    calls = []
    monkeypatch.setattr(envs, 'install', lambda names, **kw: calls.append((names, kw)) or {'ok': True})
    assert envs.ensure_python(item.name, item.version, allow_install=True) == str(envs.prefix(item) / 'bin/python')
    assert calls == [(['posebusters'], {'versions': {'posebusters': '0.5.1'}})]
