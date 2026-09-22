"""`-F /dev/null` is hygiene on a laptop and a dead end on a cluster.

Measured on 2026-09-18. A session pushed its branch successfully, then
ran

    GIT_SSH_COMMAND="ssh -F /dev/null" git fetch origin main
    → ssh: Could not resolve hostname github.com

and spent the next eight calls on it: reading /etc/ssh/ssh_config.d,
copying it somewhere writable, hitting "Bad owner or permissions", and
finally "connect to host 140.82.121.3 port 22: Permission denied". All
of it solving a problem the flag had created — the login node reaches
the outside through a ProxyCommand that lives in the very file the flag
discarded.

Nothing is refused here. The result simply says what happened, at the
moment the model is deciding what to try next, which is the only moment
it helps.
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import _discarded_ssh_config_note


RESOLVE = "ssh: Could not resolve hostname github.com: Name or service not known"
REMOTE = "fatal: Could not read from remote repository."
DENIED = "ssh: connect to host 140.82.121.3 port 22: Permission denied"


@pytest.mark.parametrize("cmd,err", [
    ('GIT_SSH_COMMAND="ssh -F /dev/null" git fetch origin main', RESOLVE),
    ("git -c core.sshCommand='ssh -F /dev/null' pull", REMOTE),
    ('GIT_SSH_COMMAND="ssh -o StrictHostKeyChecking=no" git push', RESOLVE),
    ('ssh -F /dev/null git@github.com', DENIED),
])
def test_the_discarded_config_is_named(cmd, err):
    note = _discarded_ssh_config_note(cmd, 128, err)
    assert "ProxyCommand" in note
    assert "WITHOUT overriding" in note
    assert "Do not build a replacement config" in note, (
        "the session built one, and that is the eight calls")


def test_a_command_that_kept_the_config_says_nothing():
    """The same failure without the override is a different problem —
    a real outage, a wrong remote — and must not be explained away."""
    assert _discarded_ssh_config_note("git fetch origin main", 128,
                                      RESOLVE) == ""


def test_a_command_that_succeeded_says_nothing():
    assert _discarded_ssh_config_note(
        'GIT_SSH_COMMAND="ssh -F /dev/null" git fetch', 0, "") == ""


def test_an_unrelated_failure_says_nothing():
    assert _discarded_ssh_config_note(
        'GIT_SSH_COMMAND="ssh -F /dev/null" git fetch', 1,
        "fatal: couldn't find remote ref nope") == ""


def test_the_hint_reaches_the_result(tmp_path):
    """Through the tool, not only the helper: a note nobody is handed is
    a note nobody reads."""
    import json
    from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

    perms = KitToolPermissions(workspace=tmp_path)
    perms.mode = "bypassPermissions"
    eng = _DocToolExecutor.__new__(_DocToolExecutor)
    eng._permissions = perms

    out = json.loads(eng._execute_bash({
        "command": ('GIT_SSH_COMMAND="ssh -F /dev/null" sh -c '
                    '"echo ssh: Could not resolve hostname github.com >&2; '
                    'exit 128"'),
        "description": "fetch"}, perms))

    assert out.get("exit_code") == 128
    assert "ProxyCommand" in str(out.get("hint") or ""), out
