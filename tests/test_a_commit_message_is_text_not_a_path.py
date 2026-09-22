"""A commit message is text. The gate must not read it as a path.

On 2026-09-17 a session wrote a module that skips secret-looking files
and said so in its commit message. The message named what it skips --
env files, keys, credentials -- and the bash gate refused the commit
twice, once for the word "Secrets" and once for "credentials", as if
the command were reaching for them. The session reworded its own commit
to get past its own gate.

A message is never run. What IS run inside one is a substitution, and
that still reads as a command. So:

  the words of a message        not paths, not commands
  a substitution inside it      still read, still refused
  the file flag is a file       -F names a path, and a path is scanned
  everything after the message  scanned as before
  only git's message flag       another command's -m is untouched
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor, _prose_blanked


@pytest.fixture
def perms(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    return KitToolPermissions(workspace=ws, mode="bypassPermissions",
                              confirm_callback=None)


def _gate(cmd: str, perms) -> str | None:
    return _doc_executor._run_permission_gate("bash", {"command": cmd}, perms)


#: The command as the field sent it, message and all.
FIELD_COMMAND = (
    'git add delfin/agent/report_junk.py tests/test_report_junk.py && '
    'git commit -q -m "Add read-only workspace junk report\n\n'
    'Walks a workspace without git and lists what git status hides:\n'
    'large files, __pycache__/*.pyc, ORCA/xtb leftovers, *.lock, and\n'
    'per-top-folder sizes. Secrets (.env*, *.key, *.pem, credentials*)\n'
    'are skipped entirely; nothing is read or deleted." && git log --oneline -1'
)


def test_the_commit_that_was_refused_twice_goes_through(perms):
    assert _gate(FIELD_COMMAND, perms) is None


@pytest.mark.parametrize("message", [
    "Secrets (.env*, *.key, *.pem, credentials*) are skipped",
    "read ~/.ssh/id_rsa is exactly what this must never do",
    "document .git-credentials handling",
    "rm -rf / was the footgun this removes",
    "run this on shutdown",
])
def test_words_in_a_message_are_words(perms, message):
    assert _gate(f'git commit -m "{message}"', perms) is None


@pytest.mark.parametrize("cmd", [
    'git commit -m "$(cat ~/.ssh/id_rsa)"',
    'git commit -m "see `cat ~/.ssh/id_rsa`"',
    'git commit -m "${x:-$(cat ~/.ssh/id_rsa)}"',
])
def test_a_substitution_inside_a_message_still_runs_and_is_refused(perms, cmd):
    """These read the file and put what they read in the commit."""
    assert _gate(cmd, perms) is not None


def test_a_plain_expansion_writes_a_name_and_reads_nothing(perms):
    """``${HOME}/.ssh/id_rsa`` in a message is the NAME of a file: the
    shell expands a variable, nothing opens anything. Refusing it would
    be the same mistake one layer down."""
    assert _gate('git commit -m "never read ${HOME}/.ssh/id_rsa"', perms) is None


@pytest.mark.parametrize("cmd", [
    "cat ~/.ssh/id_rsa",
    "git commit -F ~/.ssh/id_rsa",
    "git commit --file ~/.ssh/id_rsa",
    'git commit -m "fine" && cat ~/.ssh/id_rsa',
    'git commit -m "fine"; cat .env',
    'echo -m ~/.ssh/id_rsa',
])
def test_everything_else_is_read_as_before(perms, cmd):
    assert _gate(cmd, perms) is not None


def test_a_command_after_a_message_is_still_a_command(perms):
    assert _gate('git commit -m "tidy up" && rm -rf /', perms) is not None


@pytest.mark.parametrize("cmd", [
    'git tag -m "credentials handling" v1',
    'git notes add -m "keys are skipped"',
    "git commit -m wip-on-credentials",
    'GIT_EDITOR=true git commit -m "credentials"',
])
def test_the_other_places_git_carries_a_message(perms, cmd):
    assert _gate(cmd, perms) is None


def test_a_stash_is_refused_for_being_a_stash_not_for_its_message(perms):
    """`git stash push -m "wip on .env parsing"` stood in the list above
    until the stash went on the deny list: the stack belongs to the
    repository, and sessions in sibling worktrees share it.

    It is still worth a case, for the half this file is about — the
    refusal must be about the stash, never about the .env in the
    message."""
    err = _gate('git stash push -m "wip on .env parsing"', perms)
    assert err is not None
    assert "stash" in err
    assert ".env" not in err and "secret" not in err.lower()


# -- the blanking itself ----------------------------------------------------

def test_the_blanking_keeps_every_other_offset_where_it_was():
    cmd = 'git commit -m "secret words here" && ls'
    blanked = _prose_blanked(cmd)
    assert len(blanked) == len(cmd)
    assert blanked.endswith('" && ls')
    assert "secret" not in blanked


def test_a_command_without_a_message_is_untouched():
    for cmd in ("cat .env", "git status", "git log --oneline -5"):
        assert _prose_blanked(cmd) == cmd
