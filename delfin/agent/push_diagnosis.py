"""What a refused `git push` actually means, in one line.

Input: the text git wrote when a push failed. Output: the cause, and the
one thing that resolves it -- or "" when the text says nothing this knows.

Why. Across one wave of six agent sessions on 2026-09-28, two of them
spent eleven and thirteen attempts on a push that could not work from
that host. Every failure was legible, and none of them was read: the
sessions varied the command instead -- other transports, other proxies,
other temp directories -- because nothing told them the first message was
final. The user pushed by hand in the end.

This reads the message. It does not push, does not reach a network, does
not relax anything: an agent that learns in one attempt what it learned
in thirteen asks the user sooner, which is the outcome the push gate
wants anyway.

Every pattern here comes from a message a session actually received.
"""

from __future__ import annotations

import re
from typing import NamedTuple


class Diagnosis(NamedTuple):
    """What went wrong, whose it is, and what ends it."""
    cause: str
    #: "host" — the machine or the account, not something the agent can
    #: change from inside a repository. "policy" — DELFIN refused on
    #: purpose. "repo" — ordinary git state the agent can fix itself.
    owner: str
    remedy: str

    @property
    def agent_can_fix(self) -> bool:
        return self.owner == "repo"


#: (pattern, cause, owner, remedy). Ordered: the first match wins, so the
#: specific messages come before the generic "could not read from remote",
#: which git prints on top of most of them.
_PATTERNS: tuple[tuple[str, str, str, str], ...] = (
    # The tool itself is absent or logged out. Both arrived raw -- bash
    # exit 127, or gh's own stderr -- with nothing here mapping them, while
    # the write gate routes changes to the default branch through gh. The
    # first two patterns before the transport ones: a missing binary is not
    # a network problem, and the generic matches below would not fire on it
    # at all.
    (r"\bgh: (?:command )?not found|command not found: gh|"
     r"\bgh\b.*No such file or directory",
     "the GitHub CLI is not installed on this host",
     "host",
     "install gh (https://cli.github.com), or hand the user a compare URL "
     "to open the pull request themselves"),
    (r"gh auth login|not logged in(?:to| to) |"
     r"To get started with GitHub CLI, please run",
     "the GitHub CLI has no active login",
     "host",
     "the user runs 'gh auth login'; an agent cannot hold the credential "
     "for them"),
    (r"\bgit: (?:command )?not found|command not found: git",
     "git is not installed on this host",
     "host",
     "install git; nothing about branching, committing or pushing works "
     "without it"),
    (r"Please tell me who you are|empty ident name|"
     r"unable to auto-detect email address",
     "git has no commit identity on this host",
     "repo",
     "git config user.name and user.email, then commit again"),
    (r"Bad owner or permissions on .*ssh",
     "the host's SSH config is rejected by ssh itself",
     "host",
     "an administrator fixes the mode or owner of that file; no transport "
     "or proxy works around it"),
    (r"Could not resolve hostname",
     "this node has no DNS for the remote",
     "host",
     "push from a node with outbound access, or ask the user to push"),
    (r"Missing or invalid credentials|could not read Username|"
     r"Authentication failed|terminal prompts disabled",
     "no usable credentials for the remote",
     "host",
     "the user logs the credential helper in; an agent cannot supply a "
     "password it must not hold"),
    (r"EACCES .*vscode-git.*\.sock|vscode-git-[0-9a-f]+\.sock",
     "git is asking another user's VS Code for the credential",
     "host",
     "unset the credential helper for this checkout, or ask the user to "
     "push from their own session"),
    (r"Read-only file system.*tmp|mktemp: failed to create file",
     "TMPDIR points at a read-only filesystem",
     "repo",
     "set TMPDIR to a writable directory for the push"),
    (r"publishes to a shared remote",
     "DELFIN refused the push: the user has not asked for one",
     "policy",
     "say what would be pushed, to which branch, and whether its tests "
     "are green -- then push when they answer"),
    (r"\bnon-fast-forward\b|Updates were rejected",
     "the remote has commits this branch does not",
     "repo",
     "fetch and rebase onto the remote branch, then push again"),
    (r"protected branch|refusing to allow|GH00[0-9]",
     "the remote's branch protection refused it",
     "host",
     "push to a work branch instead, or ask the user to land it"),
    (r"Permission denied \(publickey\)",
     "the remote does not accept this key",
     "host",
     "the user adds the key to the account, or pushes themselves"),
    (r"Could not read from remote repository",
     "the remote could not be reached, for a reason git printed above",
     "host",
     "read the line before this one; varying the transport will not "
     "change it"),
)


def diagnose(text: str) -> Diagnosis | None:
    """The first known cause in *text*, or None.

    None means "this message is not one of the known ones" and not "the
    push was fine": a caller that gets None should show the text itself
    rather than inventing a reason for it.
    """
    body = str(text or "")
    if not body.strip():
        return None
    for pattern, cause, owner, remedy in _PATTERNS:
        if re.search(pattern, body, re.I):
            return Diagnosis(cause=cause, owner=owner, remedy=remedy)
    return None


def explain(text: str) -> str:
    """One line for a human, or "" when nothing is recognised."""
    found = diagnose(text)
    if found is None:
        return ""
    whose = {"host": "This is the machine or the account, not the "
                     "repository",
             "policy": "This is DELFIN's own rule",
             "repo": "This one is fixable from here"}[found.owner]
    return f"push refused: {found.cause}. {whose} — {found.remedy}."
