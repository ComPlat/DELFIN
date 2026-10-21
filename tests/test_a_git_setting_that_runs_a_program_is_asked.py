"""A git setting whose value git later runs as a program is asked.

The bash gate asks before ``git config`` writes ``core.hooksPath`` or
``core.fsmonitor`` in every mode, and refuses with nobody to ask. These
tests pin the same rule for the other keys whose value git executes:
pager, editor, ssh command, askpass, sequence editor, gpg program,
credential helper, diff/filter/merge drivers, ``!`` aliases, the include
path, and the same keys set for one command with ``git -c``,
``--config-env`` or the ``GIT_CONFIG_*`` environment. Driven through the
real gate (``_run_permission_gate``); nothing is executed.

Control: reading a key, unsetting it, and harmless keys are not asked.
"""

from __future__ import annotations

import subprocess

import pytest

from conftest import child_env
from delfin.agent import job_monitor as jm
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

pytestmark = pytest.mark.skipif(
    subprocess.run(["git", "--version"], capture_output=True).returncode != 0,
    reason="git is required to build the repository")


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")
    cfg = tmp_path / "gitconfig"
    cfg.write_text("[user]\n\tname = t\n\temail = t@t\n")
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", str(cfg))
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    root = tmp_path / "repo"
    root.mkdir()
    env = child_env(tmp_path)
    env["GIT_CONFIG_GLOBAL"] = str(cfg)
    env["GIT_CONFIG_NOSYSTEM"] = "1"
    subprocess.run(["git", "init", "-q", str(root)], check=True,
                   capture_output=True, env=env)
    return root


def _perms(ws, mode, answer=None):
    perms = KitToolPermissions(workspace=ws, mode=mode)
    perms.asked = []
    if answer is not None:
        def _cb(*a, **k):
            perms.asked.append(a)
            return answer
        perms.confirm_callback = _cb
    return perms


def _gate(perms, cmd):
    return _DocToolExecutor()._run_permission_gate(
        "bash", {"command": cmd}, perms)


PERSISTENT = [
    "git config core.sshCommand 'ssh -i k'",
    "git config --global core.pager less",
    "git config core.PAGER 'sh x'",
    "git config set core.editor vim",
    "git config core.askPass ./p.sh",
    "git config sequence.editor ./e.sh",
    "git config gpg.program ./g.sh",
    "git config gpg.ssh.program ./g.sh",
    "git config GPG.X509.Program ./g.sh",
    "git config credential.helper store",
    "git config --add credential.https://example.org.helper ./h.sh",
    "git config diff.pdf.textconv ./t.sh",
    "git config diff.pdf.command ./t.sh",
    "git config diff.external ./t.sh",
    "git config filter.lfs.clean ./c.sh",
    "git config filter.lfs.smudge ./c.sh",
    "git config filter.lfs.process ./c.sh",
    "git config merge.ours.driver ./m.sh",
    "git config alias.x '!sh -c id'",
    "git config Alias.y '-c core.pager=./p.sh log'",
    "git config pager.log ./p.sh",
    "git config include.path ./other.cfg",
    "git config includeIf.gitdir:~/x/.path ./other.cfg",
    "git config core.hooksPath h",
]

ONE_COMMAND = [
    "git -c core.pager=./p.sh log",
    "git -c core.sshCommand='ssh -i k' fetch",
    "git -C . -c diff.external=./t.sh diff",
    "git --config-env=core.editor=ED commit",
    "git -c alias.x='!id' x",
    "GIT_CONFIG_PARAMETERS=\"'core.pager'='./p.sh'\" git log",
    "GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=core.pager "
    "GIT_CONFIG_VALUE_0=./p.sh git log",
    "export GIT_CONFIG_KEY_0=Credential.Helper; git fetch",
    "sh -c 'git -c core.fsmonitor=./m.sh status'",
]


@pytest.mark.parametrize("cmd", PERSISTENT + ONE_COMMAND)
@pytest.mark.parametrize("mode", ["default", "acceptEdits", "bypassPermissions"])
def test_a_setting_that_runs_a_program_is_asked_in_every_mode(ws, cmd, mode):
    perms = _perms(ws, mode, answer=False)
    err = _gate(perms, cmd)
    assert perms.asked and "GIT CONFIGURATION" in str(perms.asked[0]), (
        cmd, mode, perms.asked)
    assert err, (cmd, mode)
    headless = _perms(ws, mode)
    err = _gate(headless, cmd)
    assert err and err.startswith("blocked:") and \
        "git configuration" in err.lower(), (cmd, mode, err)


@pytest.mark.parametrize("cmd", PERSISTENT[:3] + ONE_COMMAND[:1])
def test_an_approved_setting_is_asked_once(ws, cmd):
    perms = _perms(ws, "default", answer=True)
    _gate(perms, cmd)
    assert len(perms.asked) == 1, (cmd, perms.asked)


@pytest.mark.parametrize("cmd", [
    "git config --get core.pager",
    "git config get core.sshCommand",
    "git config --get-regexp alias",
    "git config --unset credential.helper",
    "git config core.pager",
    "git config -l",
    "git config user.name t",
    "git config --global user.email t@t",
    "git config alias.st status",
    "git config pager.branch false",
    "git config color.ui auto",
    "git -c color.ui=always log -1",
    "git -c user.name=t log -1",
    "git log --format='core.pager x'",
    "GIT_CONFIG_KEY_0=user.name git log -1",
    "echo 'git config core.pager less'",
])
def test_reading_and_harmless_keys_are_not_asked(ws, cmd):
    perms = _perms(ws, "bypassPermissions", answer=True)
    assert _gate(perms, cmd) is None, cmd
    assert not perms.asked, (cmd, perms.asked)
