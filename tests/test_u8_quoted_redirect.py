"""U8: a quoted '>' in a commit message is not a redirection target.

Operator finding (wave 14a, 2026-10-09): ``_bash_write_targets`` read words
of a QUOTED commit message as redirection targets -- ``git commit -m
"fix -> apply_proposal"`` was refused as "would write to 'apply_proposal'".
The ``_redirect`` regex (api_client._bash_write_targets) ran on the raw
segment before shlex, so a ``>`` inside quotes was misread as an operator.

This is the phase-3 red test for the quoted-'>' fix (protected patch to
api_client.py). Both directions are pinned:

  * a ``>`` / ``->`` / ``=>`` INSIDE quotes is text, never a write target;
  * a real ``>`` / ``>>`` / ``2>`` OUTSIDE quotes is still a target
    (fail-closed must survive the fix).
"""

from __future__ import annotations

from delfin.agent.api_client import _bash_write_targets


def test_quoted_arrow_in_commit_message_is_not_a_target():
    cmd = 'git commit -m "fix: use -> apply_proposal here"'
    assert _bash_write_targets(cmd) == []


def test_quoted_equals_arrow_is_not_a_target():
    cmd = 'git commit -m "calc => new_value"'
    assert _bash_write_targets(cmd) == []


def test_quoted_greater_than_is_not_a_target():
    cmd = 'git commit -m "a > b stays literal"'
    assert _bash_write_targets(cmd) == []


def test_quoted_arrow_in_single_quoted_message_is_not_a_target():
    cmd = "git commit -m 'fix -> apply_proposal'"
    assert _bash_write_targets(cmd) == []


def test_unquoted_real_redirection_stays_a_target():
    assert _bash_write_targets("echo hi > out.txt") == ["out.txt"]


def test_unquoted_append_stays_a_target():
    assert _bash_write_targets("ls >> log.txt") == ["log.txt"]


def test_unquoted_fd_redirect_stays_a_target():
    assert _bash_write_targets("cmd 2> err.txt") == ["err.txt"]


def test_heredoc_body_arrow_is_not_a_target():
    # git commit -F - reads the message from stdin; a '>' in the body is text.
    cmd = "git commit -F -"
    assert _bash_write_targets(cmd) == []


def test_unquoted_arrow_after_command_still_refused():
    # Without quotes `>` IS a real shell redirection; the name after it is a
    # target. Fail-closed must not be lost.
    targets = _bash_write_targets("git commit -m x -> apply_proposal")
    assert "apply_proposal" in targets
