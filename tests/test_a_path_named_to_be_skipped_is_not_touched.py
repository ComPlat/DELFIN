"""A path a command is told to SKIP is not a path it touches.

The secret scan reads every path-shaped token on a command line and
refuses the command if one matches a protected glob. It reads the path,
not its role — so on 2026-09-18 a session got:

    find agent_workspace/base-run -type f -not -path "*/.git/*" -delete
    → bash command references a secret-deny path ('/.git/*')

while the same intent, spelled the way grep spells it, ran:

    grep -rn foo --exclude-dir=.git .        → allowed

Same meaning, opposite answers, told apart by nothing a caller can see.
`.git` is on the protected list on purpose — its config can carry
credentials — and that is unchanged here: only the VALUE of an exclusion
is dropped before the scan. Every other token on the line is still read,
so a command that excludes one secret and reads another is caught on the
half that reads.
"""

from __future__ import annotations

import pathlib
import tempfile

import pytest

from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


@pytest.fixture()
def gate(tmp_path):
    (tmp_path / ".git").mkdir()
    perms = KitToolPermissions(workspace=tmp_path)
    perms.mode = "bypassPermissions"
    eng = _DocToolExecutor.__new__(_DocToolExecutor)
    eng._permissions = perms

    def _run(cmd: str):
        return eng._run_permission_gate("bash", {"command": cmd}, perms)
    return _run


# -- the reported shape, and its siblings ------------------------------------

@pytest.mark.parametrize("cmd", [
    'find . -type f -not -path "*/.git/*" -delete',
    'find . -type f ! -path "*/.git/*" -print',
    'find . ! -name "*.key" -print',
    "grep -rn foo --exclude-dir=.git .",
    "grep -rn foo --exclude=*.key .",
    "rsync -a src/ dst/ --exclude .env",
    "tar -cf out.tar . --exclude=.ssh",
    "git log -- . ':(exclude).git'",
])
def test_an_excluded_path_does_not_refuse_the_command(gate, cmd):
    assert gate(cmd) is None, cmd


# -- what must still be refused ---------------------------------------------

@pytest.mark.parametrize("cmd", [
    "ls .git/objects",
    "cat ~/.ssh/id_rsa",
    "cat .env",
    "head config/server.key",
])
def test_a_path_the_command_reads_is_still_refused(gate, cmd):
    assert gate(cmd) is not None, cmd


def test_excluding_one_secret_does_not_hide_another(gate):
    """The exclusion drops its own value and nothing else. A line that
    skips one protected path and reads a second is caught on the second."""
    assert gate("find . --exclude=.git -exec cat ~/.ssh/id_rsa ;") is not None


def test_the_exclusion_value_alone_is_dropped(gate):
    """`--exclude=.git` goes; the directory being searched does not."""
    assert gate("grep -rn x --exclude-dir=.git ~/.ssh") is not None
