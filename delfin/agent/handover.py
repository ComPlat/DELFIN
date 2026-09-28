"""The command the user can run, when the agent is not the one to run it.

Input: a refused command and WHY it was refused. Output: a line offering
that command for the user to run with ``!``, or "" when offering it would
be wrong.

Why this exists. On 2026-09-28 two agent sessions spent eleven and
thirteen attempts each on a push that could not work from that host, and
the user pushed by hand in the end. Nothing in a refusal said "this one
is yours", so the sessions kept looking for a way it could be theirs --
which is the behaviour a refusal must never provoke.

The line that decides everything here: a refusal about CONSENT may be
handed over, a refusal about CONTENT may not.

  consent   the user is entitled to do this and the agent is not without
            being asked -- publishing to a shared remote, sending data
            out, acting under an installation role it does not hold.
            Handing the command over asks the person who may decide.

  content   the command itself is the problem -- a deny-pattern match, a
            path that escapes the workspace, a read of somebody's keys.
            Offering it to the user is inviting the same harm by another
            hand, politely. Nothing here does that, and a test asserts it
            for every content kind by name.

Nothing in this module runs, spawns or reaches anything: it formats a
string. The decision that produced the refusal is untouched and is made
before this is called.
"""

from __future__ import annotations

#: Refusals about who may decide. The command is unchanged and the user
#: is the one entitled to run it.
_CONSENT_KINDS: frozenset = frozenset({
    "push_unrequested",       # publishing to a shared remote, unasked
    "git_role",               # an installation role the agent does not hold
})

#: Refusals about the command itself. Never handed over, and named here
#: so that adding a kind is a decision somebody makes on purpose.
_CONTENT_KINDS: frozenset = frozenset({
    "deny_pattern",           # rm -rf and its family
    "secret_path",            # .ssh, .env, credentials
    "path_escape",            # out of the workspace
    "workspace_stray",        # a path belonging to somebody else
    "filesystem_walk",        # a walk over a shared file system
    "data_egress",            # the user's content leaving the machine
    # These two read as consent -- both refusals say a human would make
    # them fine -- and are deliberately NOT handed over. A push names a
    # branch on the project's own remote and is bounded by that. A
    # generic outbound transfer or a write to another machine is bounded
    # by nothing: handing one to a person who pastes it is the same
    # transfer, made by somebody who did not see what was in it. A rule
    # that can be satisfied by asking somebody else is not a rule.
    "outbound_unattended",
    "remote_write",
})


def kinds() -> dict:
    """Every kind this module knows, and whether it may be handed over."""
    return {**{k: True for k in _CONSENT_KINDS},
            **{k: False for k in _CONTENT_KINDS}}


def for_user(command: str, *, kind: str) -> str:
    """The handover line, or "".

    An unknown kind is treated as content: a refusal nobody has classified
    is not one to pass to the user, and the safe default for an unknown
    rule is the strict one.
    """
    text = " ".join(str(command or "").split())
    if not text or str(kind) not in _CONSENT_KINDS:
        return ""
    # One line, because it is read inside a refusal. The "!" prefix is how
    # this session runs a command the user types.
    return (f"\n\nIf you want this done, you can run it yourself:\n"
            f"  ! {text}")
