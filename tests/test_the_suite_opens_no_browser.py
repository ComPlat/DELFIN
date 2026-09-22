"""The test suite never opens the developer's browser.

delfin-voila opens its dashboard by itself when it believes it runs in a VS
Code terminal. A suite started in one inherits the variables that say so,
and on 2026-09-15 a test that launched a real server opened dashboard tabs
on the developer's machine.
"""

from __future__ import annotations

import os

from delfin import cli_voila


def test_no_test_sees_a_vscode_browser_hook():
    assert "VSCODE_IPC_HOOK_CLI" not in os.environ
    assert "BROWSER" not in os.environ
    assert os.environ.get("TERM_PROGRAM") != "vscode"
    assert cli_voila._is_vscode_session() is False
