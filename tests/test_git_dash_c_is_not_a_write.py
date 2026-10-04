"""``git -C <dir>`` names where git runs, not a file it writes.

The write-target parser read git's -C like tar's, so a read such as
``git -C <worktree> branch --show-current`` reported a write to the
worktree root. Inside the workspace that was silent; with a repository
write scope it became a refusal, and three of them ended a session.
"""
from __future__ import annotations

import pytest

from delfin.agent.api_client import _bash_write_targets


@pytest.mark.parametrize("cmd", [
    "git -C /x/y branch --show-current",
    "git -C /x symbolic-ref --short HEAD",
    "git -C ../other log --oneline -3",
])
def test_git_dash_c_is_not_a_target(cmd):
    assert _bash_write_targets(cmd) == []


@pytest.mark.parametrize("cmd, target", [
    ("git diff --output=o.txt", "o.txt"),
    ("git -C /x diff --output=o.txt", "o.txt"),
    ("tar -C /d -xf a.tar", "/d"),
])
def test_real_destinations_stay_targets(cmd, target):
    assert target in _bash_write_targets(cmd)
