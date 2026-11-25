"""Install.cmd can set up DELFIN on this PC, in a WSL distribution of its own.

The Windows half (Setup-Wsl.ps1) is checked by Test-Launcher.ps1 on
Windows; here: the Linux half (wsl_setup.sh) run under bash with stand-in
programs, and the rules the Windows half has to keep.
"""
import os
import re
import stat
import subprocess
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1] / 'tools' / 'windows-launcher'
_SETUP = _ROOT / 'wsl_setup.sh'
_REPO = 'https://github.com/ComPlat/DELFIN.git'
_COMMIT = re.search(r"\$script:Commit = '([0-9a-f]{40})'", (_ROOT / 'Setup-Wsl.ps1').read_text()).group(1)

pytestmark = pytest.mark.skipif(os.name == 'nt', reason='bash half; Windows half is Test-Launcher.ps1')


def _stub(folder, name, body):
    path = folder / name
    path.write_text('#!/bin/sh\n' + body + '\n')
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def _run(tmp_path, *args, root=True, extra=None):
    """wsl_setup.sh fed on stdin, as Setup-Wsl.ps1 feeds it through wsl.exe."""
    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir(exist_ok=True)
    log = tmp_path / 'calls.log'
    uid = '0' if root else '1000'
    # id -u answers root or not; "id <user>" fails until useradd made the user.
    _stub(bin_dir, 'id', f'''if [ "$1" = -u ]; then echo {uid}; exit 0; fi
[ -f "{tmp_path}/user-$1" ]''')
    # apt-get reads stdin on purpose: it must not eat the rest of the script.
    _stub(bin_dir, 'apt-get', f'cat >/dev/null; echo "apt-get $*" >> "{log}"')
    _stub(bin_dir, 'useradd', f'echo "useradd $*" >> "{log}"; for a; do last=$a; done; touch "{tmp_path}/user-$last"')
    for name, body in (extra or {}).items():
        _stub(bin_dir, name, body)
    home = tmp_path / 'home'
    home.mkdir(exist_ok=True)
    env = {'PATH': f'{bin_dir}:/usr/bin:/bin', 'HOME': str(home),
           'DELFIN_SETUP_WSL_CONF': str(tmp_path / 'wsl.conf')}
    result = subprocess.run(['bash', '-s', '--', *args], input=_SETUP.read_bytes(),
                            capture_output=True, env=env, timeout=60)
    calls = log.read_text() if log.exists() else ''
    return result, calls


def test_the_system_step_installs_tmux_and_makes_the_user(tmp_path):
    result, calls = _run(tmp_path, 'system', 'delfin')
    assert result.returncode == 0, result.stderr
    assert b'system ready' in result.stdout
    assert 'apt-get install -y -qq tmux git curl ca-certificates bzip2 build-essential gfortran meson ninja-build pkg-config libopenblas-dev' in calls
    assert 'useradd -m -s /bin/bash delfin' in calls
    assert (tmp_path / 'wsl.conf').read_text() == '[user]\ndefault=delfin\n'


def test_the_system_step_runs_again_without_changing_what_is_there(tmp_path):
    _run(tmp_path, 'system', 'delfin')
    (tmp_path / 'calls.log').unlink()
    result, calls = _run(tmp_path, 'system', 'delfin')
    assert result.returncode == 0, result.stderr
    assert 'useradd' not in calls
    assert (tmp_path / 'wsl.conf').read_text() == '[user]\ndefault=delfin\n'


def test_a_distribution_that_names_its_user_keeps_it(tmp_path):
    (tmp_path / 'wsl.conf').write_text('[user]\ndefault=someone\n')
    result, _ = _run(tmp_path, 'system', 'delfin')
    assert result.returncode == 0, result.stderr
    assert (tmp_path / 'wsl.conf').read_text() == '[user]\ndefault=someone\n'


@pytest.mark.parametrize('user', ['root', 'Delfin', '-rf', 'a b', 'x;id', ''])
def test_the_system_step_refuses_a_bad_user_name(tmp_path, user):
    result, calls = _run(tmp_path, 'system', user)
    assert result.returncode == 2
    assert calls == ''


def test_the_system_step_needs_root(tmp_path):
    result, calls = _run(tmp_path, 'system', 'delfin', root=False)
    assert result.returncode == 2
    assert b'needs root' in result.stdout
    assert calls == ''


def _fake_git(tmp_path, *, writes_env=True):
    """git clone puts an install.sh there that records how it was called."""
    log = tmp_path / 'calls.log'
    env_line = 'touch "$HOME/.delfin_env.sh"; mkdir -p "$here/.venv/bin"; ' \
               'printf "#!/bin/sh\\n" > "$here/.venv/bin/delfin"; chmod +x "$here/.venv/bin/delfin"' \
        if writes_env else 'true'
    install = f'here=$(cd "$(dirname "$0")" && pwd); echo "install.sh $*" >> "{log}"; {env_line}'
    return {'git': f'''echo "git $*" >> "{log}"
if [ "$1" = -C ]; then
  case "$3" in
    checkout) echo "$7" > "$2/.git/HEAD_COMMIT" ;;
    rev-parse) cat "$2/.git/HEAD_COMMIT" ;;
  esac
  exit 0
fi
dest=$6; mkdir -p "$dest/.git"
cat > "$dest/install.sh" <<'INSTALL'
{install}
INSTALL'''}


def test_the_delfin_step_clones_and_runs_the_installer(tmp_path):
    result, calls = _run(tmp_path, 'delfin', _REPO, _COMMIT, root=False,
                         extra=_fake_git(tmp_path))
    assert result.returncode == 0, result.stderr
    home = tmp_path / 'home'
    assert f'git clone -q --branch main {_REPO} {home}/DELFIN.partial' in calls
    assert f'git -C {home}/DELFIN.partial checkout -q -B main {_COMMIT}' in calls
    assert (home / 'DELFIN' / 'install.sh').exists() and not (home / 'DELFIN.partial').exists()
    assert 'install.sh --profile standard --no-orca --no-openmpi' in calls
    assert result.stdout.rstrip().endswith(b'DELFIN_WSL_READY')


def test_the_delfin_step_runs_again_without_cloning_again(tmp_path):
    """A second run keeps the checkout and lets install.sh retry what failed."""
    extra = _fake_git(tmp_path)
    _run(tmp_path, 'delfin', _REPO, _COMMIT, root=False, extra=extra)
    (tmp_path / 'calls.log').unlink()
    result, calls = _run(tmp_path, 'delfin', _REPO, _COMMIT, root=False, extra=extra)
    assert result.returncode == 0, result.stderr
    assert 'git clone' not in calls and 'checkout' not in calls
    assert calls.splitlines()[-1] == 'install.sh --profile standard --no-orca --no-openmpi'
    assert b'DELFIN_WSL_READY' in result.stdout


def test_a_failed_install_is_not_reported_ready(tmp_path):
    result, _ = _run(tmp_path, 'delfin', _REPO, _COMMIT, root=False,
                     extra=_fake_git(tmp_path, writes_env=False))
    assert result.returncode == 4
    assert b'DELFIN_WSL_READY' not in result.stdout


def test_the_delfin_step_runs_as_the_user_not_root(tmp_path):
    result, calls = _run(tmp_path, 'delfin', _REPO, _COMMIT, extra=_fake_git(tmp_path))
    assert result.returncode == 2
    assert calls == ''


@pytest.mark.parametrize('repo', ['http://example.org/x.git', 'https://x.org/a b', 'https://x.org/$(id)', '--upload-pack=x', ''])
def test_the_delfin_step_refuses_a_bad_repository(tmp_path, repo):
    result, calls = _run(tmp_path, 'delfin', repo, _COMMIT, root=False, extra=_fake_git(tmp_path))
    assert result.returncode == 2
    assert calls == ''


def _function(text, name):
    return re.search(r'^function ' + name + r'\b.*?^\}', text, re.M | re.S).group(0)


def test_setup_passes_only_fixed_values_to_wsl():
    setup = (_ROOT / 'Setup-Wsl.ps1').read_text(encoding='utf-8')
    arguments = _function(setup, 'Get-SetupArguments')
    assert "'-d ' + $script:Distribution + ' -u root --cd / -- bash -s -- system ' + $script:LinuxUser" in arguments
    assert "' --cd ~ -- bash -s -- delfin ' + $script:Repository + ' ' + $script:Commit" in arguments
    assert "$script:Distribution = 'DELFIN'" in setup
    assert "$script:LinuxUser = 'delfin'" in setup
    assert "$script:Repository = 'https://github.com/ComPlat/DELFIN.git'" in setup
    # Its own distribution: the person's other distributions stay untouched.
    assert "'--install ' + $script:Image + ' --name ' + $script:Distribution + ' --no-launch'" in setup
    # The Linux half goes on stdin, with LF line ends.
    assert "-replace \"`r`n\", \"`n\"" in _function(setup, 'Get-SetupScriptBytes')


def test_only_installing_wsl_itself_asks_for_administrator_rights():
    setup = (_ROOT / 'Setup-Wsl.ps1').read_text(encoding='utf-8')
    assert setup.count('-Verb RunAs') == 1
    feature = _function(setup, 'Install-WslFeature')
    assert "-ArgumentList '--install --no-distribution' -Verb RunAs -Wait -PassThru" in feature
    # Declining the prompt is said in words.
    assert 'NativeErrorCode -eq 1223' in feature
    for name in (n for n in os.listdir(_ROOT) if n.endswith(('.ps1', '.cmd')) and n != 'Setup-Wsl.ps1'):
        assert 'RunAs' not in (_ROOT / name).read_text(encoding='utf-8'), name


def test_the_restart_resume_runs_once_and_is_cleaned_up():
    """RunOnce under HKCU: Windows deletes it before it runs; no autostart stays."""
    setup = (_ROOT / 'Setup-Wsl.ps1').read_text(encoding='utf-8')
    assert "$script:ResumeKey = 'HKCU:\\Software\\Microsoft\\Windows\\CurrentVersion\\RunOnce'" in setup
    assert 'Register-Resume (Join-Path $PSScriptRoot' in _function(setup, 'Invoke-Setup')
    # Finished: a resume still pending is removed.
    assert 'Unregister-Resume' in setup.split('if ($FunctionsOnly) { return }')[1]
    uninstall = (_ROOT / 'Uninstall.ps1').read_text(encoding='utf-8')
    assert "Remove-ItemProperty -LiteralPath $resume -Name 'DELFINWslSetup'" in uninstall
    # Nothing in the package writes the lasting Run key or the Startup folder.
    for name in (n for n in os.listdir(_ROOT) if n.endswith(('.ps1', '.cmd'))):
        text = (_ROOT / name).read_text(encoding='utf-8')
        assert not re.search(r'(?i)CurrentVersion\\Run(?!Once)\b|\\Startup\\|HKLM:', text), name


def test_install_asks_before_setting_up_wsl():
    install = (_ROOT / 'Install.cmd').read_text(encoding='utf-8')
    ask = install.index('choice /C YN /N /M "Set up DELFIN on this PC too? [Y/N]: "')
    run = install.index('"%LOCALAPPDATA%\\DELFIN Launcher\\Setup-Wsl.ps1"')
    assert ask < run
    assert install[ask:run].count('if errorlevel 2 goto startapp') == 1
    # The installed copy runs, so the restart resume points at a lasting path.
    for name in ('Install.ps1', 'Uninstall.ps1'):
        text = (_ROOT / name).read_text(encoding='utf-8')
        assert "'Setup-Wsl.ps1','wsl_setup.sh'" in text, name


def test_known_wsl_failures_are_explained_in_words():
    setup = (_ROOT / 'Setup-Wsl.ps1').read_text(encoding='utf-8')
    errors = _function(setup, 'Get-FriendlyWslError')
    for code, words in (('0x80370102', 'BIOS/UEFI'), ('0x8007019e', 'Restart Windows'),
                        ('0x80072ee7', 'internet'), ('0x80070005', 'administrator rights')):
        assert code in errors and words in errors, code
    assert 'too old for WSL 2' in setup
    assert '$script:MinimumBuild = 19041' in setup


def test_a_download_that_is_not_the_pinned_commit_is_thrown_away(tmp_path):
    extra = _fake_git(tmp_path)
    # rev-parse answers another commit than the one asked for.
    extra['git'] = extra['git'].replace('rev-parse) cat "$2/.git/HEAD_COMMIT"', 'rev-parse) echo ' + 'f' * 40)
    result, calls = _run(tmp_path, 'delfin', _REPO, _COMMIT, root=False, extra=extra)
    assert result.returncode == 5
    assert 'install.sh' not in calls
    assert not (tmp_path / 'home' / 'DELFIN').exists() and not (tmp_path / 'home' / 'DELFIN.partial').exists()


@pytest.mark.parametrize('commit', ['', 'main', '6903045ff', 'A' * 40, 'g' * 40, '-' + 'a' * 39])
def test_the_delfin_step_needs_a_full_commit_id(tmp_path, commit):
    result, calls = _run(tmp_path, 'delfin', _REPO, commit, root=False, extra=_fake_git(tmp_path))
    assert result.returncode == 2
    assert calls == ''


def test_the_pinned_commit_is_on_main():
    """When the checkout has the history: the pin is a commit main contains."""
    repo = _ROOT.parents[1]
    known = subprocess.run(['git', '-C', str(repo), 'cat-file', '-e', _COMMIT + '^{commit}'], capture_output=True)
    if known.returncode != 0:
        pytest.skip('shallow checkout: the pinned commit is not in it')
    ancestor = subprocess.run(['git', '-C', str(repo), 'merge-base', '--is-ancestor', _COMMIT, 'HEAD'], capture_output=True)
    assert ancestor.returncode == 0


def test_uninstall_removes_the_distribution_only_when_asked():
    uninstall = (_ROOT / 'Uninstall.ps1').read_text(encoding='utf-8')
    assert re.search(r'^\s*\[switch\]\$RemoveWslDistribution,\s*$', uninstall, re.M)
    block = uninstall[uninstall.index('if ($RemoveWslDistribution -and $linux) {'):]
    assert "& $wsl --unregister DELFIN" in block
    assert uninstall.count('& $wsl --unregister') == 1
    # Asked only when the distribution is there; Enter keeps it.
    ask = uninstall.index("$answer = Read-Host \"Also delete DELFIN's Linux part on this PC")
    assert uninstall.index('if ($linux -and -not $RemoveWslDistribution) {') < ask < uninstall.index('& $wsl --unregister')
    assert "$RemoveWslDistribution = $answer -eq 'y' -or $answer -eq 'Y'" in uninstall
    assert "$linux = $names -contains 'DELFIN'" in uninstall


def test_nothing_turns_off_a_windows_protection():
    for name in (n for n in os.listdir(_ROOT) if n.endswith(('.ps1', '.cmd'))):
        text = (_ROOT / name).read_text(encoding='utf-8')
        assert not re.search(r'(?i)ExecutionPolicy\s+(Bypass|Unrestricted)|Set-MpPreference|Add-MpPreference|schtasks|Register-ScheduledTask', text), name
