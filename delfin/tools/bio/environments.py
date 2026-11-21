"""Shared bio tool recipes, discovery and isolated installation.

Listings never install. Execution may ensure a tool using the existing DELFIN
opt-out; managed prefixes are recipe-specific and never alter a cluster env.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import shlex
import subprocess
import sys


@dataclass(frozen=True)
class Recipe:
    name: str
    version: str
    modules: tuple[str, ...] = ()
    commands: tuple[str, ...] = ()
    dependencies: tuple[str, ...] = ()


RECIPES = {
    r.name: r for r in (
        Recipe('vina', '1.2.5', commands=('vina',)),
        Recipe('meeko', '0.7.1', modules=('meeko',), commands=('mk_prepare_ligand.py', 'mk_prepare_receptor.py', 'mk_export.py'), dependencies=('rdkit', 'gemmi', 'scipy')),
        Recipe('pdbfixer', '1.12', modules=('pdbfixer', 'openmm'), dependencies=('openmm=8.4.0',)),
        Recipe('openmm', '8.4.0', modules=('openmm',)),
        Recipe('openmmforcefields', '0.15.0', modules=('openmmforcefields.generators', 'openff.toolkit', 'openmm', 'rdkit'),
               dependencies=('openff-toolkit=0.17.0', 'openmm=8.4.0', 'rdkit')),
        Recipe('mdanalysis', '2.9.0', modules=('MDAnalysis', 'MDAnalysis.analysis.align', 'numpy'), dependencies=('biopython',)),
        Recipe('prolif', '2.2.3', modules=('prolif', 'MDAnalysis', 'rdkit')),
        Recipe('plip', '3.0.1', modules=('plip', 'openbabel'), commands=('plip',)),
        Recipe('posebusters', '0.5.1', modules=('posebusters', 'posebusters.cli', 'rdkit', 'pandas', 'numpy', 'yaml'),
               commands=('bust',), dependencies=('rdkit>=2024.9', 'pandas', 'numpy', 'pyyaml')),
    )
}
_TRIED = set()


def root() -> Path:
    return Path(os.environ.get('DELFIN_BIO_TOOLS_ROOT', str(Path.home() / '.local/share/delfin/bio_tools'))).expanduser().resolve()


def recipe(name: str, version: str | None = None) -> Recipe:
    if name not in RECIPES:
        raise ValueError(f'No managed bio installation recipe for {name}')
    value = RECIPES[name]
    if version is not None:
        if not isinstance(version, str) or not re.fullmatch(r'[0-9][A-Za-z0-9.!+_-]*', version):
            raise ValueError('version must be one exact release, not a solver expression or URL')
        value = replace(value, version=version)
    return value


def specification(item: Recipe) -> list[str]:
    return ['python=3.11', item.name + '=' + item.version, *item.dependencies]


def prefix(item: Recipe) -> Path:
    digest = hashlib.sha256(json.dumps(specification(item)).encode()).hexdigest()[:16]
    return root() / 'envs' / (item.name + '-' + digest)


def python_candidates(item: Recipe):
    explicit = os.environ.get('DELFIN_' + item.name.upper() + '_PYTHON')
    candidates = [explicit] if explicit else []
    candidates.extend([sys.executable, shutil.which('python')])
    from delfin import qm_runtime
    for command in item.commands:
        found = qm_runtime.resolve_tool(command)
        if found:
            candidates.append(str(Path(found.path).parent / 'python'))
    candidates.append(str(prefix(item) / 'bin/python'))
    seen = set()
    for candidate in candidates:
        if candidate and candidate not in seen and Path(candidate).is_file() and os.access(candidate, os.X_OK):
            seen.add(candidate)
            yield candidate


def present(name: str, version: str | None = None) -> bool:
    """Presence listing only; shared executable discovery may probe cluster modules."""
    item = recipe(name, version)
    if item.modules:
        for python in python_candidates(item):
            sites = list((Path(python).parent.parent / 'lib').glob('python*/site-packages'))
            if any(all((site.joinpath(*module.split('.')) / '__init__.py').is_file() or site.joinpath(*module.split('.')).with_suffix('.py').is_file()
                       for module in item.modules) for site in sites):
                return True
        return False
    from delfin import qm_runtime
    return any(qm_runtime.resolve_tool(command) is not None for command in item.commands)


def _module_python_on_path(module: str):
    from delfin import qm_runtime
    from delfin.system_tools import SystemToolCandidate
    script = f'module load {shlex.quote(module)} >/dev/null 2>&1 && (command -v python || command -v python3)'
    try:
        result = qm_runtime.run_tool('bash', ['-lc', script], capture_output=True,
                                     timeout=20, track_process=True, allow_install=False)
        lines = (result.stdout or '').strip().splitlines()
        if result.returncode == 0 and lines:
            path = Path(lines[-1])
            if path.is_absolute() and path.is_file() and os.access(path, os.X_OK):
                return SystemToolCandidate(str(path), 'module', module)
    except (OSError, subprocess.TimeoutExpired):
        pass
    return None


def module_python_candidates(item: Recipe):
    """Reuse QM discovery, preserving module identity even for a shared Python."""
    from delfin.system_tools import discover_system_tool_candidates, matching_modules
    patterns = (item.name, 'bio/' + item.name, 'chem/' + item.name)
    # Search one module at a time: different modules may expose the same Python.
    for module in matching_modules(patterns):
        primary = _module_python_on_path(module)
        if primary:
            yield primary
        for candidate in discover_system_tool_candidates(
                ('python', 'python3'), module_patterns=(module,),
                module_env_hints=('EBROOT' + item.name.upper(), item.name.upper() + '_HOME')):
            if candidate.module_name == module and (not primary or candidate.path != primary.path):
                yield candidate


def check_python(item: Recipe, python: str, *, environment_module: str | None = None) -> dict:
    from delfin import qm_runtime
    code = ('import importlib,importlib.metadata as m,json,sys; '
            'mods=json.loads(sys.argv[2]); [importlib.import_module(x) for x in mods]; '
            'print(json.dumps({"version":m.version(sys.argv[1])}))')
    flags = ['-s', '-P'] if environment_module else ['-I']
    done = qm_runtime.run_tool(python, [*flags, '-c', code, item.name, json.dumps(item.modules)],
                               capture_output=True, timeout=60, track_process=True, allow_install=False,
                               environment_module=environment_module)
    if done.returncode:
        raise ValueError((done.stderr or done.stdout or 'import check failed')[-2000:])
    result = json.loads(done.stdout)
    if result['version'] != item.version:
        raise ValueError(f'{item.name}: expected {item.version}, found {result["version"]}')
    return {'name': item.name, 'version': result['version'], 'python': python,
            'check': 'distribution version and imports; not a scientific reference run'}


def installation_command(item: Recipe, manager: str) -> list[str]:
    return [manager, 'create', '--yes', '--prefix', str(prefix(item)), '--override-channels',
            '--channel', 'conda-forge', '--strict-channel-priority', *specification(item)]


def check_binary(item: Recipe, binary: str) -> dict:
    from delfin import qm_runtime
    done = qm_runtime.run_tool(binary, ['--version'], capture_output=True, timeout=30,
                               track_process=True, allow_install=False)
    banner = ((done.stdout or '') + (done.stderr or '')).strip()
    if done.returncode or not re.search(r'(?<![0-9.])' + re.escape(item.version) + r'(?![0-9.])', banner):
        raise ValueError('installed binary version check failed: ' + banner)
    return {'binary': binary, 'version': item.version, 'banner': banner}


def install(names, *, on_line=None, timeout=1800, versions=None, repair=False) -> dict:
    from delfin import installer, qm_runtime
    import fcntl
    outcomes = []
    for name in names:
        item = recipe(name, (versions or {}).get(name))
        manager = next((shutil.which(n) for n in ('micromamba', 'mamba', 'conda') if shutil.which(n)), None)
        if not manager:
            outcomes.append({'name': name, 'ok': False, 'lines': ['No micromamba, mamba or conda available; use the existing DELFIN bootstrap first.']})
            continue
        directory = root()
        directory.mkdir(parents=True, exist_ok=True)
        target = prefix(item)
        with (directory / (target.name + '.lock')).open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            manifest = directory / (target.name + '.json')
            command = installation_command(item, manager)
            if target.exists():
                try:
                    if item.modules:
                        check_python(item, str(target / 'bin/python'))
                    else:
                        check_binary(item, str(target / 'bin' / item.commands[0]))
                    outcomes.append({'name': name, 'ok': True, 'lines': ['Existing managed environment retained.']})
                    continue
                except (ValueError, OSError, subprocess.TimeoutExpired):
                    pass
                try:
                    owned = json.loads(manifest.read_text()) if manifest.is_file() else {}
                except (ValueError, OSError):
                    owned = {}
                if not isinstance(owned, dict):
                    owned = {}
                if not repair or owned.get('prefix') != str(target) or owned.get('specification') != specification(item):
                    outcomes.append({'name': name, 'ok': False, 'lines': [f'Unhealthy prefix: {target}; explicit repair requires its matching ownership manifest.']})
                    continue
                command[1] = 'install'
                command.insert(2, '--force-reinstall')
            ownership = {'name': name, 'prefix': str(target), 'specification': specification(item), 'ok': False}
            manifest.write_text(json.dumps(ownership, indent=2))
            if on_line:
                on_line(f'{name} {item.version}: installing into isolated environment {target}')
            ok, lines = installer._run(command, cwd=directory, env=os.environ.copy(), on_line=on_line, timeout=timeout)
            qm_runtime.clear_resolver_cache()
            health = None
            if ok:
                try:
                    if item.modules:
                        health = check_python(item, str(target / 'bin/python'))
                    else:
                        health = check_binary(item, str(target / 'bin' / item.commands[0]))
                except (ValueError, OSError, subprocess.TimeoutExpired) as exc:
                    ok = False
                    lines.append(str(exc))
            # Conda metadata retains exact dependency versions, builds and hashes.
            packages = []
            for metadata in sorted((target / 'conda-meta').glob('*.json')):
                try:
                    value = json.loads(metadata.read_text())
                    if not isinstance(value, dict) or any(not value.get(key) for key in ('name', 'version', 'build')):
                        raise ValueError('missing package identity')
                    packages.append({key: value.get(key) for key in ('name', 'version', 'build', 'channel', 'sha256', 'md5')})
                except (ValueError, OSError) as exc:
                    ok = False
                    lines.append(f'Invalid package provenance {metadata.name}: {exc}')
            manifest.write_text(json.dumps({**ownership, 'ok': ok, 'argv': command, 'health': health,
                                            'packages': packages, 'lines': lines}, indent=2))
            outcomes.append({'name': name, 'ok': ok, 'lines': lines, 'manifest': str(manifest)})
    return {'ok': all(row['ok'] for row in outcomes), 'results': outcomes}


def ensure_python(name: str, version: str, *, allow_install: bool = False) -> str:
    item = recipe(name, version)
    failures = []
    for python in python_candidates(item):
        try:
            check_python(item, python)
            return python
        except (ValueError, OSError, subprocess.TimeoutExpired) as exc:
            failures.append(str(exc))
    from delfin.dashboard.gfn_optimize import auto_install_allowed
    key = (name, version, str(root()))
    if allow_install and auto_install_allowed() and key not in _TRIED:
        _TRIED.add(key)
        result = install([name], versions={name: version})
        if result['ok']:
            return str(prefix(item) / 'bin/python')
        failures.extend(line for row in result['results'] for line in row['lines'])
    raise ImportError(f'{name} {version} unavailable: ' + '; '.join(failures or ['no compatible environment; automatic installation not performed']))


def ensure_binary(name: str, version: str) -> str:
    item = recipe(name, version)
    binary = str(prefix(item) / 'bin' / item.commands[0])
    if Path(binary).is_file():
        check_binary(item, binary)
        return binary
    from delfin.dashboard.gfn_optimize import auto_install_allowed
    if not auto_install_allowed():
        raise ImportError(f'{name}: automatic installation is disabled')
    key = (name, version, str(root()))
    if key in _TRIED:
        raise ImportError(f'{name}: installation already attempted this session; inspect the installation manifest')
    _TRIED.add(key)
    outcome = install([name], versions={name: version})
    if not outcome['ok']:
        raise ImportError('; '.join(line for row in outcome['results'] for line in row['lines']))
    return str(prefix(item) / 'bin' / item.commands[0])
