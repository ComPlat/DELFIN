"""A deny-listed path is refused with the reason that actually applies.

The path deny list holds two different kinds of file: things that hold a
secret (``.ssh``, ``.env``, ``*.key``, ``credentials*``) and git's own
state under ``.git/``, which holds no secret at all and is blocked
because writing it changes what git does next -- a hook runs on the next
git command, and ``.git/info/exclude`` hides files from ``git status``,
so work can leave review unseen.

One message for both told a caller reaching for ``.git/info/exclude``
that it had touched a secret. That reads as a mistake in the rule rather
than a rule, so it is retried in another spelling: measured 8 retries of
that one path in a single agent run, against a guard working exactly as
intended.

Input: a bash command naming a deny-listed path.
Output: the refusal text.

Semantics: the refusal still refuses -- the test asserts the command is
denied in both cases. Only which rule it names changes.
"""
from __future__ import annotations

import pytest

from delfin.agent.api_client import _why_path_is_denied


@pytest.mark.parametrize("path", [
    ".git/info/exclude",
    ".git/config",
    ".git/hooks/pre-commit",
    "sub/project/.git/info/exclude",
    ".git",
])
def test_git_state_is_not_called_a_secret(path):
    why = _why_path_is_denied(path)
    assert "not a secret" in why
    # The two mechanisms, so the caller knows why a harmless-looking
    # file is refused and does not look for a spelling that passes.
    assert "hook" in why
    assert "git status" in why


@pytest.mark.parametrize("path", [
    ".ssh/id_ed25519",
    ".env",
    "deploy.key",
    "credentials.json",
    "secrets.yaml",
])
def test_a_real_secret_still_says_secret(path):
    why = _why_path_is_denied(path)
    assert ".ssh/.env/*.key/credentials" in why
    assert "not a secret" not in why


def test_a_path_named_gitignore_is_not_git_state():
    """A tracked ``.gitignore`` is the route the message points at; it
    must not be described as git's own state."""
    assert "not a secret" not in _why_path_is_denied(".gitignore")


@pytest.mark.parametrize("path", ["", None, "   "])
def test_nothing_named_still_answers(path):
    assert _why_path_is_denied(path)


def test_the_refusal_still_refuses(tmp_path):
    """The message is the only thing this changed: the command is denied."""
    from delfin.agent.api_client import KitToolPermissions, _doc_executor
    perms = KitToolPermissions(workspace=tmp_path, mode="bypassPermissions",
                               confirm_callback=None)
    refusal = _doc_executor._run_permission_gate(
        "bash", {"command": "echo hide-me >> .git/info/exclude"}, perms)
    assert refusal, "writing git's own state must stay refused"
    assert "not a secret" in refusal
