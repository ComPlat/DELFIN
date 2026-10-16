"""A preview server does not overwrite the note that says where the dashboard runs.

Input: the launcher's environment. Output of ``cli_voila._announces_where``:
True for an ordinary dashboard, False for a branch preview (the
``DELFIN_PREVIEW_*`` markers ``branch_preview.preview_env`` sets). The note
(``~/.delfin/dashboard_here.json``) is the way back to *the* dashboard;
a preview pointing it at itself would send ``delfin-agent where`` and a
later login to a server that closes with the preview.
"""

from __future__ import annotations

from delfin import cli_voila
from delfin.dashboard import branch_preview as BP


def test_an_ordinary_launch_announces_itself():
    assert cli_voila._announces_where({}) is True


def test_a_preview_launch_does_not():
    env = {BP.PREVIEW_ENV_BRANCH: "feat", BP.PREVIEW_ENV_WORKTREE: "/x"}
    assert cli_voila._announces_where(env) is False
