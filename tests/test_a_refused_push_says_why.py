"""A push that cannot work says so once, instead of thirteen times.

Input: the text git wrote when a push failed. Output: the cause, whose it
is, and the one thing that ends it.

Measured across one wave of six agent sessions, 2026-09-28: two of them
spent eleven and thirteen attempts on a push that could not work from
that host. Every message was legible and none was read -- the sessions
varied the command instead, trying other transports, other proxies, other
temp directories, because nothing told them the first answer was final.
The user pushed by hand in the end.

Every pattern in the table comes from a message one of those sessions
actually received. None of this pushes, reaches a network or relaxes
anything: an agent that learns in one attempt what it learned in thirteen
asks the user sooner, which is what the push rule wants anyway.

Deliberately NOT here: a workaround for the host's broken ssh config.
Ignoring it means ignoring its host-key policy too, and a session that
tried exactly that had to add StrictHostKeyChecking=accept-new to get
past. Trading host verification for convenience is not a fix.
"""

from __future__ import annotations

import pytest

from delfin.agent.push_diagnosis import diagnose, explain

#: The real messages, as the sessions received them.
_SSH_CONFIG = ("Bad owner or permissions on /etc/ssh/ssh_config.d/"
               "50-redhat.conf\nfatal: Could not read from remote "
               "repository.")
_NO_DNS = "ssh: Could not resolve hostname github.com: Name or service not known"
_CREDS = ("Missing or invalid credentials.\nError: connect EACCES "
          "/run/user/507211/vscode-git-93588b7db8.sock")
_TMPDIR = ("mktemp: failed to create file via template "
           "'/scratch/tmp.XXXXXXXXXX': Read-only file system")
_GATE = ("blocked: `git push` publishes to a shared remote, and the user "
         "has not asked for a push since their last message")


@pytest.mark.parametrize("text, owner", [
    (_SSH_CONFIG, "host"),
    (_NO_DNS, "host"),
    (_CREDS, "host"),
    (_TMPDIR, "repo"),
    (_GATE, "policy"),
])
def test_a_real_message_is_placed_with_its_owner(text, owner):
    found = diagnose(text)
    assert found is not None, "a message a session really got is unrecognised"
    assert found.owner == owner
    assert found.remedy, "no remedy is named"


def test_only_the_repo_kind_is_the_agents_to_fix():
    """The other two end with somebody else, and saying so is the point:
    an agent that knows it cannot fix this asks rather than retries."""
    assert diagnose(_TMPDIR).agent_can_fix
    assert not diagnose(_SSH_CONFIG).agent_can_fix
    assert not diagnose(_NO_DNS).agent_can_fix
    assert not diagnose(_GATE).agent_can_fix


def test_the_specific_message_wins_over_the_generic_one():
    """git prints "Could not read from remote repository" on top of most
    of these; matching that first would bury every real cause."""
    assert diagnose(_SSH_CONFIG).cause.startswith("the host's SSH config")


def test_an_unknown_message_is_not_given_an_invented_reason():
    assert diagnose("some failure nobody has seen") is None
    assert explain("some failure nobody has seen") == ""


def test_empty_input_is_not_a_diagnosis():
    assert diagnose("") is None and diagnose(None) is None


def test_the_line_a_person_reads_names_all_three_parts():
    line = explain(_CREDS)
    assert "push refused:" in line
    assert "credentials" in line
    assert "machine or the account" in line, "it does not say whose this is"


def test_nothing_here_runs_or_reaches_anything():
    """A classifier that shells out would be a second way to push."""
    import inspect

    from delfin.agent import push_diagnosis

    src = inspect.getsource(push_diagnosis)
    for forbidden in ("subprocess", "os.system", "socket", "urllib",
                      "requests", "Popen"):
        assert forbidden not in src, forbidden


# -- the doctor asks the question before the work, not after --------------

def test_the_doctor_checks_the_remote():
    import inspect

    from delfin.agent import doctor

    src = inspect.getsource(doctor)
    assert "_check_push" in src
    assert '("git remote", "_check_push")' in src, (
        "the check exists but is not in the table, so it never runs")
    assert "ls-remote" in src, "the probe pushes instead of asking"


def test_the_doctor_probe_publishes_nothing():
    """`git ls-remote` asks what refs exist. Anything that writes would
    make a diagnostic into a second push path."""
    import inspect

    from delfin.agent.doctor import _check_push

    src = inspect.getsource(_check_push)
    # The property is which COMMAND runs, not which word appears: the body
    # talks about pushing in its comments, and asserting on the word made
    # this test fail while the probe was correct.
    import re as _re
    commands = _re.findall(r"\[([^\]]*?)\]", src)
    argv = " ".join(c for c in commands if '"git"' in c)
    assert "ls-remote" in argv, "the probe does not ask the remote"
    assert '"push"' not in argv, "the probe pushes instead of asking"
    assert '"fetch"' not in argv and '"clone"' not in argv, (
        "the probe writes to the checkout")


def test_a_foreign_credential_helper_is_named_not_rewritten():
    """It can never work -- the socket belongs to another uid -- and it is
    still the user's configuration."""
    import inspect

    from delfin.agent.doctor import _check_push

    src = inspect.getsource(_check_push)
    assert "credential.helper" in src
    assert "--unset" in src, "no remedy is offered"
    assert "config --unset" not in src.replace('"git config --unset '
                                               'credential.helper in this '
                                               'checkout, "', ""), (
        "the check rewrites the user's git configuration")
