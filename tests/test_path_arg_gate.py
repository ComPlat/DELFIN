"""The permission gate must accept the file_path alias, like the executors do.

Bug 2026-06-25 (ka_ew7404, kit.qwen3.5-397b): the model called write_file with
{"file_path": ...} (the Claude-Code convention). The permission gate read only
"path" and rejected it with "path is required" BEFORE the executor (which DOES
accept file_path via _get_path_arg) could run — so the agent fell back to bash
heredoc writes on every file.
"""

import tempfile
from pathlib import Path

from delfin.agent.api_client import _DocToolExecutor, KitToolPermissions


def _perms():
    return KitToolPermissions(workspace=Path(tempfile.mkdtemp()), mode="default",
                              confirm_callback=lambda *a, **k: True)


def test_write_file_accepts_file_path_through_gate():
    out = _DocToolExecutor().execute(
        "write_file", {"file_path": "a.py", "content": "hi"}, permissions=_perms())
    assert "path is required" not in out
    assert "a.py" in out


def test_write_file_still_accepts_path():
    out = _DocToolExecutor().execute(
        "write_file", {"path": "b.py", "content": "hi"}, permissions=_perms())
    assert "path is required" not in out


def test_edit_file_accepts_file_path_through_gate():
    ex, perms = _DocToolExecutor(), _perms()
    ex.execute("write_file", {"path": "e.py", "content": "x=1"}, permissions=perms)
    out = ex.execute("edit_file",
                     {"file_path": "e.py", "old_string": "x=1", "new_string": "x=2"},
                     permissions=perms)
    assert "path is required" not in out
