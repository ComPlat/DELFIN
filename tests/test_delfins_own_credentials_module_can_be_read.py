"""A committed, unmodified Python module named like a secret can be read.

Driven 2026-09-16: grep_file on delfin/agent/credentials.py -- the code that
manages ~/.delfin/credentials.json -- was refused by the "credentials*" name
glob. Only that case opens: the file must be Python source, tracked by git
and unmodified, and match no other deny glob.
"""
import json
import subprocess

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor


def _git(ws, *args):
    subprocess.run(["git", "-C", str(ws), *args], check=True, capture_output=True)


@pytest.fixture
def repo(tmp_path):
    ws = tmp_path / "ws"
    (ws / "pkg").mkdir(parents=True)
    _git(ws, "init", "-q")
    _git(ws, "config", "user.email", "t@t")
    _git(ws, "config", "user.name", "t")
    (ws / "pkg" / "credentials.py").write_text("def load_credential(name):\n    return ''\n")
    (ws / "pkg" / "credentials.json").write_text(json.dumps({"k": "SECRET"}))
    (ws / "pkg" / ".env").write_text("TOKEN=SECRET\n")
    _git(ws, "add", "-f", "pkg/credentials.py", "pkg/credentials.json", "pkg/.env")
    _git(ws, "commit", "-q", "-m", "init")
    return ws


def _read(ws, path):
    perms = KitToolPermissions(workspace=ws, mode="bypassPermissions", confirm_callback=None)
    return _doc_executor._execute_read_file({"path": path}, perms)


def test_a_tracked_clean_module_named_credentials_is_readable(repo):
    assert "load_credential" in _read(repo, "pkg/credentials.py")


def test_a_modified_or_untracked_one_is_still_denied(repo):
    (repo / "pkg" / "credentials.py").write_text("TOKEN = 'SECRET'\n")
    assert "read denied" in _read(repo, "pkg/credentials.py")
    (repo / "pkg" / "secrets.py").write_text("x = 1\n")
    assert "read denied" in _read(repo, "pkg/secrets.py")


def test_a_tracked_json_store_and_env_stay_denied(repo):
    assert "read denied" in _read(repo, "pkg/credentials.json")
    assert "read denied" in _read(repo, "pkg/.env")


def test_a_directory_search_reads_the_module_but_not_the_store(repo):
    perms = KitToolPermissions(workspace=repo, mode="bypassPermissions", confirm_callback=None)
    out = _doc_executor._execute_grep_file({"pattern": "load_credential|SECRET", "path": "pkg"}, perms)
    assert "credentials.py:" in out
    assert "credentials.json:" not in out and ".env:" not in out
